#include "vj_query_processing.hpp"
#include "vj_hits_filter.hpp"
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    std::shared_ptr<BaseFillFixCropProcessor> VJQueryProcessor::GetFillFixCropProcessor() {
        return std::shared_ptr<BaseFillFixCropProcessor>(
                new AggressiveFillFixCropProcessor(params_.fix_crop_fill_params,
                                                   read_archive_));
    }

    ProcessedVJHits VJQueryProcessor::ComputeFilteringResults(core::Read &read, VJHits vj_hits) {
        ProcessedVJHits processed_hits(read);
        processed_hits.vj_hits = vj_hits;
        if(params_.filtering_params.enable_filtering) {
            VersatileVjFilter vj_filter(params_.filtering_params);
            processed_hits.filtering_info = vj_filter.Filter(vj_hits);
        }
        return processed_hits;
    }

    ProcessedVJHits VJQueryProcessor::Process(core::Read &read) {
        VJQueryAligner vj_query_aligner(params_, v_db_, j_db_);
        VJHits vj_hits = vj_query_aligner.Align(read);
//        if (vj_hits.NumVHits() > 1) {
//            std::cout << vj_hits.NumVHits() << " V hits were computed" << std::endl;
//        }
        ProcessedVJHits hits_after_fitering = ComputeFilteringResults(read, vj_hits);
        if(hits_after_fitering.ReadToBeFiltered()) {
            return hits_after_fitering;
        }
//        if (hits_after_fitering.vj_hits.NumVHits() > 0) {
//            std::cout << hits_after_fitering.vj_hits.NumVHits() << " processed V hits were computed" << std::endl;
//        }


        // refinement

        hits_after_fitering.vj_hits = RefineAlignment(hits_after_fitering.vj_hits);

        // end refinement



        auto fix_fill_crop_processor = GetFillFixCropProcessor();
        hits_after_fitering.vj_hits = fix_fill_crop_processor->Process(hits_after_fitering.vj_hits);
        return hits_after_fitering;
    }

    VJHits VJQueryProcessor::RefineAlignment(const VJHits& vj_hits) {
        const auto& alignment_params = params_.scoring_params.v_scoring;
        size_t best_index = 0;
        size_t num_vhits = vj_hits.NumVHits();
        std::vector<int> exact_alignment_scores(num_vhits);

        using namespace seqan;

        for (size_t i = 0; i < num_vhits; ++i) {
            auto v_hit = vj_hits.GetVHitByIndex(i);
            const auto& gene_seq = v_hit.ImmuneGene().seq();
            const auto& read_seq = v_hit.Read().seq;

            int score = globalAlignmentScore(gene_seq,
                                             read_seq,
                                             Score<int, Simple>(alignment_params.match_reward,
                                                                -alignment_params.mismatch_opening_cost,
                                                                -alignment_params.gap_extention_cost,
                                                                -alignment_params.gap_opening_cost),
                                             AlignConfig<false, true, true, false>(),
                                             AffineGaps());
            exact_alignment_scores[i] = score;
        }
        if (std::max_element(exact_alignment_scores.begin(),
                             exact_alignment_scores.end()) != exact_alignment_scores.begin()) {
            ++num_refined_alignments_;
        }
//        std::cout << "vjf score: " << v_hit1.IntScore()
//                  << ", exact score: " << score1 << "\n";
//
//        auto alignment = alignment_converter_.ConvertToAlignment(v_hit1.ImmuneGene(),
//                                                                 vj_hits.Read(),
//                                                                 v_hit1.BlockAlignment());
//        std::cout << "\nVJF gene\n" << row(alignment.Alignment(), 0);
//        std::cout << "\nVJF read\n" << row(alignment.Alignment(), 1);
//        std::cout << "\nSeqAn gene\n" << row(align, 0);
//        std::cout << "\nSeqAn read\n" << row(align, 1) << "\n \n";
//        std::cout << std::endl << std::endl;

        return vj_hits;

    }
}