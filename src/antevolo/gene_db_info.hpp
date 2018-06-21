#pragma once

//#include "annotated_clone_by_read_constructor.hpp"

#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include <germline_db_labeling.hpp>

namespace antevolo {
    struct GeneDbInfo {
        germline_utils::CustomGeneDatabase representative_v_db_;
        germline_utils::CustomGeneDatabase j_db_;
        std::map<std::string, std::pair<germline_utils::CustomGeneDatabase,
                cdr_labeler::DbCDRLabeling>> representative_to_db_map_;
        cdr_labeler::DbCDRLabeling v_labeling_;
        cdr_labeler::DbCDRLabeling j_labeling_;

    public:
        GeneDbInfo(const germline_utils::CustomGeneDatabase& representative_v_db,
                   const germline_utils::CustomGeneDatabase& j_db,
                   const std::map<std::string, std::pair<germline_utils::CustomGeneDatabase,
                           cdr_labeler::DbCDRLabeling>>& representative_to_db_map,
                   const cdr_labeler::DbCDRLabeling& v_labeling,
                   const cdr_labeler::DbCDRLabeling& j_labeling) :
                representative_v_db_(representative_v_db),
                j_db_(j_db),
                representative_to_db_map_(representative_to_db_map),
                v_labeling_(v_labeling),
                j_labeling_(j_labeling){}
    };
}




