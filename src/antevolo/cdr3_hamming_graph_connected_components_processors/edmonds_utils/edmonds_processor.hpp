#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "../../../algorithms/edmonds_dmst/edmonds_optimum_branching.hpp"

namespace antevolo {

    class EdmondsProcessor {
    public:
        typedef boost::property<boost::edge_weight_t, double>       EdgeProperty;
        typedef boost::adjacency_list<boost::listS,
                boost::vecS,
                boost::directedS,
                boost::no_property,
                EdgeProperty>                 Graph;
        typedef boost::graph_traits<Graph>::vertex_descriptor       Vertex;
        typedef boost::graph_traits<Graph>::edge_descriptor         Edge;

        struct WeightedEdge {
            size_t src_;
            size_t dst_;
            double weight_;

            WeightedEdge() :
                    src_(size_t(-1)),
                    dst_(size_t(-1)),
                    weight_(400) { }

            WeightedEdge(size_t src, size_t dst, double weight) :
                    src_(src),
                    dst_(dst),
                    weight_(weight) { }

            WeightedEdge(const WeightedEdge &oth) :
                    src_(oth.src_),
                    dst_(oth.dst_),
                    weight_(oth.weight_) { }

            WeightedEdge &operator=(const WeightedEdge& oth) {
                src_ = oth.src_;
                dst_ = oth.dst_;
                weight_ = oth.weight_;
                return *this;
            }
            bool operator== (const WeightedEdge& oth) const {
                return (src_ == oth.src_ && dst_ == oth.dst_ && weight_ == oth.weight_);
            }
            bool operator!= (const WeightedEdge& oth) const {
                return !(*this == oth);
            }
            bool operator> (const WeightedEdge& oth) const {
                return weight_ > oth.weight_;
            }
        };


        static std::vector<WeightedEdge> process_edge_list(const std::vector<WeightedEdge>& input_edges);
    };

} //antevolo



