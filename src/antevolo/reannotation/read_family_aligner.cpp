#include "read_family_aligner.hpp"
#include <algorithm>

namespace antevolo {
    void ReadFamilyAligner::Run() {

        using namespace seqan;

        size_t best_v_index;
        size_t best_v_index_rc;
        size_t best_j_index;
        int v_score;
        int v_score_rc;
        int j_score;

        core::Read read_rc = read_.ReverseComplement();

        std::tie(best_v_index, v_score) = ComputeBestGeneIndex<true>(read_);
        if (config_.cdr_labeler_config.vj_finder_config.algorithm_params.aligner_params.fix_strand) {
            std::tie(best_v_index_rc, v_score_rc) = ComputeBestGeneIndex<true>(read_rc);
            if (v_score_rc > v_score) {
                read_ = read_rc;
                best_v_index = best_v_index_rc;
                v_score = v_score_rc;
            }
        }

        //TODO: remove code duplication
        auto align_config = GetAlignConfig();

        // V
        auto v_scoring_scheme = GetScoringScheme(GetScoringParams<true>());
        assignSource(row(v_alignment_, 0), v_db_[best_v_index].seq());
        assignSource(row(v_alignment_, 1), read_.seq);

        int new_v_score = globalAlignment(v_alignment_,
                                    v_scoring_scheme,
                                    align_config,
                                    AffineGaps());
        VERIFY(new_v_score == v_score);


        std::string gene_string = AnythingToStdString(row(v_alignment_, 0));

//        auto gene_string = Dna5StringToStdString(row(v_alignment_, 0));
        auto read_string = AnythingToStdString(row(v_alignment_, 1));
        size_t fix_left = config_.cdr_labeler_config.vj_finder_config.algorithm_params.fix_crop_fill_params.fix_left;
        std::tie(gene_string, read_string) = FixSides(gene_string, read_string, fix_left, true);
        read_.seq = GetReadFromAlignment(read_string);

        // J
        std::tie(best_j_index, j_score) = ComputeBestGeneIndex<false>(read_);

        auto j_scoring_scheme = GetScoringScheme(GetScoringParams<false>());
        assignSource(row(j_alignment_, 0), j_db_[best_j_index].seq());
        assignSource(row(j_alignment_, 1), read_.seq);
        int new_j_score = globalAlignment(j_alignment_,
                                          j_scoring_scheme,
                                          align_config,
                                          AffineGaps());
        VERIFY(new_j_score == j_score);
        gene_string = AnythingToStdString(row(j_alignment_, 0));
        read_string = AnythingToStdString(row(j_alignment_, 1));
        size_t fix_right = config_.cdr_labeler_config.vj_finder_config.algorithm_params.fix_crop_fill_params.fix_right;
        std::tie(gene_string, read_string) = FixSides(gene_string, read_string, fix_right, false);
        read_.seq = GetReadFromAlignment(read_string);
        assignSource(row(j_alignment_, 0), j_db_[best_j_index].seq());
        assignSource(row(j_alignment_, 1), read_.seq);
        int final_j_score = globalAlignment(j_alignment_,
                                            j_scoring_scheme,
                                            align_config,
                                            AffineGaps());
        size_t j_begin_view_position = toViewPosition(row(j_alignment_, 0), 0);
        setClippedBeginPosition(row(j_alignment_, 0), j_begin_view_position);
        setClippedBeginPosition(row(j_alignment_, 1), j_begin_view_position);


        assignSource(row(v_alignment_, 0), v_db_[best_v_index].seq());
        assignSource(row(v_alignment_, 1), read_.seq);
        int final_v_score = globalAlignment(v_alignment_,
                                            v_scoring_scheme,
                                            align_config,
                                            AffineGaps());
        size_t v_end_view_position = toViewPosition(row(v_alignment_, 0),
                                                    length(v_db_[best_v_index].seq()));
        setClippedEndPosition(row(v_alignment_, 0), v_end_view_position);
        setClippedEndPosition(row(v_alignment_, 1), v_end_view_position);

        v_id_ = best_v_index;
        j_id_ = best_j_index;

//        std::cout << row(v_alignment_, 1) << std::endl
//                  << row(v_alignment_, 0) << std::endl << "+++++++++++++++++" << std::endl;
//        std::cout << row(j_alignment_, 1) << std::endl
//                  << row(j_alignment_, 0) << std::endl << "=================" << std::endl;

    }

    template<bool is_v>
    std::pair<size_t, int> ReadFamilyAligner::ComputeBestGeneIndex(const core::Read& read) {
        using namespace seqan;

        std::vector<int> scores;
        auto scoring_scheme = GetScoringScheme(GetScoringParams<is_v>());
        auto align_config = GetAlignConfig();
        const auto& db = is_v ? v_db_ : j_db_;

        for (size_t i = 0; i < db.size(); ++i) {
            const auto& gene = db[i];
            int score = globalAlignmentScore(gene.seq(),
                                             read.seq,
                                             scoring_scheme,
                                             align_config,
                                             AffineGaps());
            scores.push_back(score);
        }
        auto best_score_it = std::max_element(scores.begin(), scores.end());
        size_t best_index = static_cast<size_t>(std::distance(scores.begin(), best_score_it));
        return std::make_pair(best_index, *best_score_it);
    }

    std::pair<std::string, std::string> ReadFamilyAligner::FixSides(std::string gene_alignment,
                                                                    std::string read_alignment,
                                                                    size_t length_to_fix,
                                                                    bool is_left) {
        if (!is_left) {
            std::reverse(gene_alignment.begin(), gene_alignment.end());
            std::reverse(read_alignment.begin(), read_alignment.end());
        }
        VERIFY(gene_alignment.length() == read_alignment.length());
        size_t length = gene_alignment.length();
        //find gene start
        size_t gene_start = 0;
        while (gene_alignment[gene_start] == '-' && gene_start < length) { ++gene_start; }
        size_t cropped_length = gene_alignment.length() - gene_start;
        std::string gene_res = gene_alignment.substr(gene_start, cropped_length);
        std::string read_res = read_alignment.substr(gene_start, cropped_length);
        size_t current_pos = 0;
        while (current_pos < length_to_fix && current_pos < cropped_length) {
            auto c = gene_res[current_pos];
            if (c != '-') {
                read_res[current_pos] = c;
            }
            ++current_pos;
        }
        if (!is_left) {
            std::reverse(gene_res.begin(), gene_res.end());
            std::reverse(read_res.begin(), read_res.end());
        }
        return std::make_pair(gene_res, read_res);
    }

    seqan::Dna5String ReadFamilyAligner::GetReadFromAlignment(const std::string& alignment) {
        std::string res;
        for (char c : alignment) {
            if (c != '-') {
                res.push_back(c);
            }
        }
        return AnythingToStdString(res);
    }
}

