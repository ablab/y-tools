#pragma once

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <shm_model_utils/shm_model_edge_weight_calculator.hpp>
#include "antevolo_config.hpp"
#include "annotation_utils/annotated_clone_set.hpp"
#include "annotated_clone_by_read_constructor.hpp"
#include "evolutionary_tree_storage.hpp"
#include "gene_db_info.hpp"

namespace antevolo {
    class AntEvoloLaunch {
        const AntEvoloConfig& config_;

        ShmModelEdgeWeightCalculator ShmModelPosteriorCalculation(
                const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>&);
        std::vector<boost::unordered_set<size_t>> ReadClusters(
                const boost::unordered_map<std::string, size_t>& read_name_to_index);

        void LaunchDefault(const GeneDbInfo& gene_db_info,
                           annotation_utils::CDRAnnotatedCloneSet& annotated_clone_set,
                           size_t total_number_of_reads);

        void LaunchEvoQuast(const annotation_utils::CDRAnnotatedCloneSet& clone_set);

        void AnalyzeParallelEvolution(const EvolutionaryTreeStorage& trees);

        annotation_utils::CDRAnnotatedCloneSet ComputeCompressedCloneSet(
                const annotation_utils::CDRAnnotatedCloneSet& uncompressed_annotated_clone_set);

        std::string GetGeneBaseName(seqan::CharString name) const;

        std::map<std::string, std::pair<germline_utils::CustomGeneDatabase,
                                        cdr_labeler::DbCDRLabeling>> ComputeRepresentativeVToDBMap(
                const germline_utils::CustomGeneDatabase& representatives_v_db);

        GeneDbInfo ComputeGeneDbInfo();

    public:
        AntEvoloLaunch(const AntEvoloConfig& config) : config_(config) { }

        void Launch();
    };
}