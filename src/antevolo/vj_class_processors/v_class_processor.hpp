#pragma once

#include "base_gene_class_processor.hpp"
#include <cdr3_hamming_graph_component_info.hpp>
#include <clone_set_decomposers/clone_set_decomposer.hpp>

namespace antevolo {
    class VClassProcessor : public BaseGeneClassProcessor {
        std::vector<size_t> jdifference_positions_;

        EvolutionaryTree ProcessComponentWithEdmonds(SparseGraphPtr hg_component, size_t component_id,
                                                     const ShmModelEdgeWeightCalculator &edge_weight_calculator);
    public:

        VClassProcessor(CloneSetWithFakesPtr clone_set_ptr,
                        const core::DecompositionClass& decomposition_class,
                        const AntEvoloConfig& config,
                        const GeneDbInfo& gene_db_info,
                        size_t current_fake_clone_index) :
                BaseGeneClassProcessor(clone_set_ptr,
                                       decomposition_class,
                                       config,
                                       gene_db_info,
                                       current_fake_clone_index) {
            auto chain = clone_set_ptr->operator[](*decomposition_class_.cbegin()).ChainType().Chain();
            if (chain == germline_utils::ImmuneChainType::HeavyIgChain) {
                jdifference_positions_ = {17, 18, 19, 22, 25, 26, 27}; // FIXME: move to config!!
            } else if (chain == germline_utils::ImmuneChainType::KappaIgChain) {
                jdifference_positions_ = {12, 13, 14, 15, 16, 22, 23, 24};
            } else if (chain == germline_utils::ImmuneChainType::LambdaIgChain) {
                jdifference_positions_ = {12, 15, 19, 22, 23, 24, 25};
            } else {
                jdifference_positions_ = {};
            }
            VERIFY(jdifference_positions_.size() != 0);

//            auto current_representative_name = CloneSetDecomposer::GetGeneBaseName(clone_set_ptr->operator[](*decomposition_class.cbegin()).VGene().name());
//            if (current_representative_name == "IGHV1-18") {
//                for (size_t i = 0;
//                     i < gene_db_info.representative_to_db_map_.find("IGHV1-18")->second.first.size(); ++i) {
//                    INFO(gene_db_info.representative_to_db_map_.find("IGHV1-18")->second.first[i].name());
//                }
//            }


        }

        void CreateUniqueCDR3Map() override;
        std::vector<SparseGraphPtr> ComputeConnectedComponents() override;

        void ChangeJgene(const germline_utils::ImmuneGene &v_gene,
                         const germline_utils::ImmuneGene &j_gene);
        void ChangeJgeneToMax(CDR3HammingGraphComponentInfo hamming_graph_info);
        void ChangeVJgenesToMax(CDR3HammingGraphComponentInfo hamming_graph_info);

        template<class T> T GetMode(std::vector<T> v) {
            std::map<T, size_t> freq_map;
            for (auto i : v) {
                if (freq_map.find(i) == freq_map.end()) {
                    freq_map.insert({i, 0});
                }
                ++freq_map[i];
            }
            T most_frequent = std::max_element(
                    std::begin(freq_map),
                    std::end(freq_map),
                    [](const std::pair<T, size_t>& p1,
                       const std::pair<T, size_t>& p2) {
                        return p1.second < p2.second;
                    })->first;
            return most_frequent;
        }
    };
}