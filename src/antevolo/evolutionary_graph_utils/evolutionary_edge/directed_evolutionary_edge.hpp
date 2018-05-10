#pragma once

#include "base_evolutionary_edge.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace  antevolo {

    class DirectedEvolutionaryEdge : public BaseEvolutionaryEdge {
        size_t num_added_v_shms;
        size_t num_intersected_v_shms;
        size_t num_added_j_shms;
        size_t num_intersected_j_shms;
        size_t num_added_shms;
        size_t num_intersected_shms;


    public:
        DirectedEvolutionaryEdge(const annotation_utils::AnnotatedClone &src_clone,
                                 const annotation_utils::AnnotatedClone &dst_clone,
                                 size_t src_num, size_t dst_num) : BaseEvolutionaryEdge(src_clone, dst_clone,
                                                                                        src_num, dst_num) {
            edge_type = EvolutionaryEdgeType ::DirectedEdgeType;
            cdr3_distance = HammingDistance(this->src_clone->CDR3(), this->dst_clone->CDR3());
            VERIFY_MSG(this->dst_clone->VSHMs().size() + this->dst_clone->JSHMs().size() >
                       this->src_clone->VSHMs().size() + this->src_clone->JSHMs().size(),
                       "# SHMs in destination clone ("
                               << this->dst_clone->VSHMs().size() + this->dst_clone->JSHMs().size()
                               << ") does not exceed # SHMs in source clone ("
                               << this->src_clone->VSHMs().size() + this->src_clone->JSHMs().size() << ")");
            num_added_v_shms = this->dst_clone->VSHMs().size() - this->src_clone->VSHMs().size();
            num_intersected_v_shms = this->src_clone->VSHMs().size();
            num_added_j_shms = this->dst_clone->JSHMs().size() - this->src_clone->JSHMs().size();
            num_intersected_j_shms = this->src_clone->JSHMs().size();
            num_added_shms = num_added_v_shms + num_added_j_shms;
            num_intersected_shms = num_intersected_v_shms + num_intersected_j_shms;
        }
        size_t Length() const override {
            return num_added_shms + cdr3_distance;
        }

        bool IsDirected() const override { return true; }

        std::string TypeString() const override { return "directed"; }

        size_t NumAddedShms() const override { return num_added_shms; }

        size_t NumSharedShms() const override { return num_intersected_shms; }

        virtual void appendAddedVSHMsInModernFormat(std::ostream& out) const {
            auto v_blocks = annotation_utils::SHMComparator::SHMs1BlocksNotPresentInSHMs2(
                    dst_clone->VSHMs(),
                    src_clone->VSHMs());
            for (auto p : v_blocks) {
                for (size_t i = 0; i < p.second; ++i) {
                    out << "+";
                    auto shm_it = dst_clone->VSHMs().cbegin() + p.first + i;
                    shm_it->AppendInModernFormat(out);
                    out << ";";
                }
            }
        }
        virtual void appendAddedJSHMsInModernFormat(std::ostream& out) const {
            auto j_blocks = annotation_utils::SHMComparator::SHMs1BlocksNotPresentInSHMs2(
                    dst_clone->JSHMs(),
                    src_clone->JSHMs());
            for (auto p : j_blocks) {
                for (size_t i = 0; i < p.second; ++i) {
                    out << "+";
                    auto shm_it = dst_clone->JSHMs().cbegin() + p.first + i;
                    shm_it->AppendInModernFormat(out);
                    out << ";";
                }
            }
        }

    };

}
