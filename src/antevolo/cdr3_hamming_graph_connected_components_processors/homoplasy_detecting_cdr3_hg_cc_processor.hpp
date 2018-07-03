#pragma once

#include "edmonds_cdr3_hg_cc_processor.hpp"
#include <boost/unordered_map.hpp>
#include <algorithm>


namespace antevolo {
    class Homoplasy_Detecting_CDR3_HG_CC_Processor : public Edmonds_CDR3_HG_CC_Processor {
        const size_t EDGE_LENGTH_THRESHOLD = 100;
        const double EDGE_COEF = 2;

    public:
        Homoplasy_Detecting_CDR3_HG_CC_Processor(CloneSetWithFakesPtr clone_set_ptr,
                                                 const AntEvoloConfig::AlgorithmParams &config,
                                                 const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                                 CDR3HammingGraphComponentInfo& hamming_graph_info,
                                                 size_t current_fake_clone_index,
                                                 const ShmModelEdgeWeightCalculator& edge_weight_calculator)
                : Edmonds_CDR3_HG_CC_Processor(clone_set_ptr,
                                               config,
                                               clone_by_read_constructor,
                                               hamming_graph_info,
                                               current_fake_clone_index,
                                               edge_weight_calculator) {}

        EvolutionaryTree Process() override;

        void TopLevelDfs(size_t v,
                         const boost::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>>& graph_out,
                         const boost::unordered_set<std::pair<size_t, size_t>>& all_edges,
                         boost::unordered_set<std::pair<size_t, size_t>>& edges_to_be_removed);
        void InnerLevelDfs(size_t main_v,
                           size_t current_v,
                           const boost::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>>& graph_out,
                           const boost::unordered_set<std::pair<size_t, size_t>>& all_edges,
                           boost::unordered_set<std::pair<size_t, size_t>>& edges_to_be_removed,
                           boost::unordered_set<size_t>& passed);

    };
}




