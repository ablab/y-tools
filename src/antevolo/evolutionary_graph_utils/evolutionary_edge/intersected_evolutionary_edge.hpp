#pragma once

#include "base_evolutionary_edge.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace  antevolo {

    class IntersectedEvolutionaryEdge : public BaseEvolutionaryEdge {
        size_t num_individual_v_shms_;
        size_t num_intersected_v_shms_;
        size_t num_individual_j_shms_;
        size_t num_intersected_j_shms_;
        size_t num_individual_shms_;
        size_t num_intersected_shms_;
        bool is_double_mytated_;

    public:
        IntersectedEvolutionaryEdge(const annotation_utils::AnnotatedClone &src_clone,
                                    const annotation_utils::AnnotatedClone &dst_clone,
                                    size_t src_num, size_t dst_num,
                                    size_t num_individual_v_shms, size_t num_intersected_v_shms,
                                    size_t num_individual_j_shms, size_t num_intersected_j_shms)
                : BaseEvolutionaryEdge(src_clone,
                                       dst_clone,
                                       src_num,
                                       dst_num),
                  num_individual_v_shms_(num_individual_v_shms),
                  num_intersected_v_shms_(num_intersected_v_shms),
                  num_individual_j_shms_(num_individual_j_shms),
                  num_intersected_j_shms_(num_intersected_j_shms) {
            edge_type = EvolutionaryEdgeType ::IntersectedEdgeType;
            //sum
            num_individual_shms_ = num_individual_v_shms_ + num_individual_v_shms_;
            num_intersected_shms_ = num_intersected_v_shms_ + num_intersected_j_shms_;
            is_double_mytated_ = annotation_utils::SHMComparator::IndividualSHMsAreIdenticallyPositioned(
                                         src_clone.VSHMs(),
                                         dst_clone.VSHMs())   &&
                                 annotation_utils::SHMComparator::IndividualSHMsAreIdenticallyPositioned(
                                         src_clone.JSHMs(),
                                         dst_clone.JSHMs());
        }
        size_t Length() const override {
            return cdr3_distance+num_individual_shms_;
        }

        bool IsIntersected() const override { return true; }

        bool IsDoubleMutated() const {
            return is_double_mytated_;
        }

        std::string TypeString() const override { return IsDoubleMutated() ? "double_mutated" : "intersected"; }

        size_t NumAddedShms() const override { return num_individual_shms_; }

        size_t NumSharedShms() const override { return num_intersected_shms_; }

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
            auto reverse_v_blocks = annotation_utils::SHMComparator::SHMs1BlocksNotPresentInSHMs2(
                    src_clone->VSHMs(),
                    dst_clone->VSHMs());
            for (auto p : reverse_v_blocks) {
                for (size_t i = 0; i < p.second; ++i) {
                    out << "-";
                    auto shm_it = src_clone->VSHMs().cbegin() + p.first + i;
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
                VERIFY(p.first + p.second <= dst_clone->JSHMs().size());
                for (size_t i = 0; i < p.second; ++i) {
                    out << "+";
                    auto shm_it = dst_clone->JSHMs().cbegin() + p.first + i;
                    shm_it->AppendInModernFormat(out);
                    out << ";";
                }
            }
            auto reverse_j_blocks = annotation_utils::SHMComparator::SHMs1BlocksNotPresentInSHMs2(
                    src_clone->JSHMs(),
                    dst_clone->JSHMs());
            for (auto p : reverse_j_blocks) {
                VERIFY(p.first + p.second <= src_clone->JSHMs().size());
                for (size_t i = 0; i < p.second; ++i) {
                    out << "-";
                    auto shm_it = src_clone->JSHMs().cbegin() + p.first + i;
                    shm_it->AppendInModernFormat(out);
                    out << ";";
                }
            }
        }

    };

}