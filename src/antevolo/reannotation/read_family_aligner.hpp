#pragma once

#include <tuple>
#include "antevolo_config.hpp"
#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include <vj_finder_config.hpp>


namespace antevolo {

    class ReadFamilyAligner {
        const AntEvoloConfig& config_;
        const germline_utils::CustomGeneDatabase& v_db_;
        const germline_utils::CustomGeneDatabase& j_db_;
        core::Read read_;
        seqan::Align<seqan::Dna5String> v_alignment_;
        seqan::Align<seqan::Dna5String> j_alignment_;
        size_t v_id_;
        size_t j_id_;
        bool to_be_filtered_;

        seqan::AlignConfig<false, true, true, false> GetAlignConfig() {
            return seqan::AlignConfig<false, true, true, false>();
        };

        typedef vj_finder::VJFinderConfig::AlgorithmParams::ScoringParams::VScoringParams VScoringParams;
        typedef vj_finder::VJFinderConfig::AlgorithmParams::ScoringParams::JScoringParams JScoringParams;

        template<bool is_v>
        typename std::enable_if<is_v, VScoringParams>::type GetScoringParams() {
            return config_.cdr_labeler_config.vj_finder_config.algorithm_params.scoring_params.v_scoring;
        }
        template<bool is_v>
        typename std::enable_if<!is_v, JScoringParams>::type GetScoringParams() {
            return config_.cdr_labeler_config.vj_finder_config.algorithm_params.scoring_params.j_scoring;
        }
        template<class T>
        seqan::Score<int, seqan::Simple> GetScoringScheme(T params) {
            return seqan::Score<int, seqan::Simple>(params.match_reward,
                                                    -params.mismatch_opening_cost,
                                                    -params.gap_extention_cost,
                                                    -params.gap_opening_cost);
        }


        std::pair<std::string, std::string> FixSides(std::string gene_alignment,
                                                     std::string read_alignment,
                                                     size_t length_to_fix,
                                                     bool is_left);

//        std::pair<std::string, std::string> FixSides(seqan::Dna5String gene_alignment,
//                                                     seqan::Dna5String  read_alignment,
//                                                     size_t length_to_fix,
//                                                     bool is_left) {
//            return FixSides(AnythingToStdString(gene_alignment),
//                            AnythingToStdString(read_alignment),
//                            length_to_fix,
//                            is_left);
//        };

        seqan::Dna5String GetReadFromAlignment(const std::string& alignment);

        template<class T>
        std::string AnythingToStdString(const T& seq) {
            std::stringstream ss;
            ss << seq;
            return ss.str();
        }

    public:
        ReadFamilyAligner(const AntEvoloConfig& config,
                          const germline_utils::CustomGeneDatabase& v_db,
                          const germline_utils::CustomGeneDatabase& j_db,
                          core::Read& read) :
                config_(config),
                v_db_(v_db),
                j_db_(j_db),
                read_(read),
                to_be_filtered_(false) {
            using namespace seqan;
            resize(rows(v_alignment_), 2);
            resize(rows(j_alignment_), 2);
        }

        seqan::Align<seqan::Dna5String> GetVAlignment() const { return v_alignment_; }
        seqan::Align<seqan::Dna5String> GetJAlignment() const { return j_alignment_; }
        size_t GetVId() const { return v_id_; }
        size_t GetJId() const { return j_id_; }
        core::Read GetRead() const { return read_; }

        template<bool is_v>
        std::pair<size_t, int> ComputeBestGeneIndex(const core::Read& read);

        //TODO: filter

        void Run();
    };
}




