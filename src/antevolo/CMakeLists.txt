############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(antevolo CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${GRAPH_UTILS})
include_directories(${CDR_LABELER_DIR})

add_executable(antevolo
        ../cdr_labeler/cdr_config.cpp
        ../vj_finder/vj_finder_config.cpp
        ../vj_finder/germline_db_generator.cpp
        ../vj_finder/vj_alignment_structs.cpp
        ../vj_finder/vj_query_aligner.cpp
        ../vj_finder/vj_hits_filter.cpp
        ../vj_finder/vj_alignment_info.cpp
        ../vj_finder/vj_query_fix_fill_crop.cpp
        ../vj_finder/vj_query_processing.cpp
        ../vj_finder/vj_parallel_processor.cpp
        ../cdr_labeler/germline_db_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/single_loop_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr1_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr2_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr3_v_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr3_j_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/single_loop_helper.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/immune_gene_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/annotated_gene_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/immune_gene_labeling_helper.cpp
        ../cdr_labeler/germline_db_labeling.cpp
        ../cdr_labeler/read_labeler.cpp
#        ../cdr_labeler/cdr_annotator.cpp
#        compressed_cdr_set.cpp
#        cdr_output.cpp
#        diversity_analyser.cpp
#        cdr_launch.cpp
        clone_set_decomposers/vj_clone_set_decomposer.cpp
        clonally_related_candidates_calculators/clonally_related_candidates.cpp
        clonally_related_candidates_calculators/undirectred_first_tree_calculator.cpp
        clonally_related_candidates_calculators/edmonds_tarjan_DMST_calculator.cpp
        model_utils/shm_model.cpp
        evolutionary_graph_utils/evolutionary_edge.cpp
        evolutionary_graph_utils/evolutionary_edge_constructor.cpp
        evolutionary_graph_utils/evolutionary_tree.cpp
        evolutionary_graph_utils/evolutionary_graph.cpp
        evolutionary_graph_utils/evolutionary_graph_constructor.cpp
        antevolo_config.cpp
        naive_antevolo_processing.cpp
        antevolo_launch.cpp
        main.cpp)

target_link_libraries(antevolo core algorithms vdj_utils graph_utils input boost_program_options ${COMMON_LIBRARIES})
