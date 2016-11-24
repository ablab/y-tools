#pragma once

#include "cdr3_hamming_graph_info.hpp"
#include <annotation_utils/annotated_clone.hpp>

namespace antevolo {

    // the constructors should not be used! The interface is factory function getRelatedClonesIterator
    class ClonesSharingCDR3sIterator {
        enum State {
            sameCDR3, similarCDR3, isDone
        };

        size_t old_index_;
        CDR3HammingGraphInfo& hamming_graph_info_;
        SparseGraph::EdgesIterator similar_cdr3s_it_;
        SparseGraph::EdgesIterator similar_cdr3s_it_end_;
        State state_;

    public:
        ClonesSharingCDR3sIterator(const ClonesSharingCDR3sIterator& oth) = default;
        ClonesSharingCDR3sIterator(CDR3HammingGraphInfo& hamming_graph_info,
                                   size_t old_index) :
                old_index_(old_index),
                hamming_graph_info_(hamming_graph_info),
                similar_cdr3s_it_(hamming_graph_info_.GetSimilarCDR3sBeginByOldIndex(old_index_)),
                similar_cdr3s_it_end_(hamming_graph_info_.GetSimilarCDR3sEndByOldIndex(old_index_)),
                state_(State::sameCDR3) {}

        bool HasNext();
        const std::vector<size_t>& operator*() const;
        ClonesSharingCDR3sIterator& operator++();
        ClonesSharingCDR3sIterator operator++(int);

    };


    class RelatedClonesIterator {
        ClonesSharingCDR3sIterator vectors_iterator_;
        const std::vector<size_t>& clones_sharing_current_cdr3_;
        std::vector<size_t>::const_iterator current_vector_it_;

    public:
        RelatedClonesIterator(const ClonesSharingCDR3sIterator& vectors_iterator) :
            vectors_iterator_(vectors_iterator),
            clones_sharing_current_cdr3_(*vectors_iterator_),
            current_vector_it_(clones_sharing_current_cdr3_.cbegin()) {}

        bool HasNext();
        size_t Next();

    };

    RelatedClonesIterator getRelatedClonesIterator(CDR3HammingGraphInfo& hamming_graph_info, const std::string& cdr3);

    RelatedClonesIterator getRelatedClonesIterator(CDR3HammingGraphInfo& hamming_graph_info,
                                                   const annotation_utils::AnnotatedClone& clone);


}