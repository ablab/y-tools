#include "homoplasy_detecting_cdr3_hg_cc_processor.hpp"


namespace antevolo {
    EvolutionaryTree Homoplasy_Detecting_CDR3_HG_CC_Processor::Process() {
        auto& clone_set  = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();
        boost::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> graph_in;
        boost::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> graph_out;
        boost::unordered_set<std::pair<size_t, size_t>> all_edges;
        boost::unordered_set<std::pair<size_t, size_t>> edges_to_be_removed;

        EvolutionaryTree tree(clone_set_ptr_);

        vertices_nums_ = boost::unordered_set<size_t>(hamming_graph_info_.GetAllClones());
        for (auto it = vertices_nums_.begin(); it != vertices_nums_.end(); ) {
            if (clone_set[*it].AnnotationIsNotValid()) {
                vertices_nums_.erase(it++);
            }
            else {
                ++it;
            }
        }

        SetShortestDirectedParentEdges();
        auto input_edges = PrepareEdgeVector();
        auto branching_edges = EdmondsProcessor().process_edge_list(input_edges);
        SetEdges(tree, branching_edges);
        for (auto it = tree.begin(); it != tree.end(); ++it) {
            auto& edge = *it;
            if (edge->IsUndirected()) {
                graph_in[edge->DstNum()].push_back({edge->SrcNum(), edge->Length()});
                graph_out[edge->SrcNum()].push_back({edge->DstNum(), edge->Length()});
                all_edges.insert({edge->SrcNum(), edge->DstNum()});
            }
        }

        for (size_t dst_num : vertices_nums_) {
            const auto& dst_clone = clone_set[dst_num];
            graph_in[dst_num] = std::vector<std::pair<size_t, size_t>>();
            auto& edges_in = graph_in[dst_num];
            std::vector<std::pair<size_t, size_t>> edges;
            for (size_t src_num : vertices_nums_) {
                const auto& src_clone = clone_set[src_num];
                auto edge = edge_constructor->ConstructEdge(src_clone, dst_clone,
                                                            src_num, dst_num);
                if (edge->IsDirected()) {
                    edges.push_back({src_num, edge->Length()});
                }
            }
            std::sort(edges.begin(), edges.end(), [](std::pair<size_t, size_t> p1,
                                                     std::pair<size_t, size_t> p2) -> bool {
                return p1.second < p2.second;
            });
            size_t min_len = std::min(edges.empty() ? EDGE_LENGTH_THRESHOLD
                                                    : edges.front().second,
                                      edges_in.empty() ? EDGE_LENGTH_THRESHOLD
                                                       : edges_in.front().second);
            for (auto p : edges) {
                if (static_cast<double>(p.second) >= static_cast<double>(min_len) * EDGE_COEF) {
                    break;
                }
                edges_in.push_back(p);
                if (graph_out.find(p.first) == graph_out.end()) {
                    graph_out[p.first] = std::vector<std::pair<size_t, size_t>>();
                }
                graph_out[p.first].push_back({dst_num, p.second});
                all_edges.insert({p.first, dst_num});
            }
        }

//        INFO("in edges");
//        for (auto p : graph_in) {
//            std::cout << p.first << " <-\n";
//            for (auto p2 : p.second) {
//                std::cout << "\t" << p2.first;
//            }
//            std::cout << std::endl;
//        }
//        INFO("out edges");
//        for (auto p : graph_out) {
//            std::cout << p.first << " ->\n";
//            for (auto p2 : p.second) {
//                std::cout << "\t" << p2.first;
//            }
//            std::cout << std::endl;
//        }
//        INFO("all edges");
//        for (auto p : all_edges) {
//            std::cout << p.first << " -> " << p.second;
//            std::cout << std::endl;
//        }


        // removing transitive edges: an edge (u -> v) is transitive
        // if v is reachable from w and the edge (u -> w) exists
        for (size_t v : vertices_nums_) {
            TopLevelDfs(v,
                        graph_out,
                        all_edges,
                        edges_to_be_removed);
        }
        for (auto p : edges_to_be_removed) {
            all_edges.erase(p);
        }


        EvolutionaryTree res_tree(clone_set_ptr_);
        for (auto p : all_edges) {
            auto edge_ptr = edge_constructor->ConstructEdge(clone_set[p.first],
                                                            clone_set[p.second],
                                                            p.first,
                                                            p.second);
//            std::cout << p.first << " " << p.second << std::endl;
//            std::cout << edge_ptr->SrcNum() << " " << edge_ptr->DstNum() << "\nh\n" << std::endl;
            res_tree.AddEdge(edge_ptr->DstNum(), edge_ptr);
        }
//        for (auto it = res_tree.cbegin(); it != res_tree.cend(); ++it) {
//            auto edge_ptr = *it;
//            std::cout << edge_ptr->SrcNum() << " " << edge_ptr->DstNum() << "\nh\n" << std::endl;
//        }
        return res_tree;
    }

    void Homoplasy_Detecting_CDR3_HG_CC_Processor::TopLevelDfs(size_t v,
                     const boost::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>>& graph_out,
                     const boost::unordered_set<std::pair<size_t, size_t>>& all_edges,
                     boost::unordered_set<std::pair<size_t, size_t>>& edges_to_be_removed) {
        auto out_edges_it = graph_out.find(v);
        if (out_edges_it == graph_out.end()) {
            return;
        }
        for (auto p : out_edges_it->second) {
            boost::unordered_set<size_t> passed;
            size_t dst_num = p.first;
            InnerLevelDfs(v,
                          dst_num,
                          graph_out,
                          all_edges,
                          edges_to_be_removed,
                          passed);

        }
    }
    void Homoplasy_Detecting_CDR3_HG_CC_Processor::InnerLevelDfs(size_t main_v,
                       size_t current_v,
                       const boost::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>>& graph_out,
                       const boost::unordered_set<std::pair<size_t, size_t>>& all_edges,
                       boost::unordered_set<std::pair<size_t, size_t>>& edges_to_be_removed,
                       boost::unordered_set<size_t>& passed) {
        passed.insert(current_v);
        auto out_edges_it = graph_out.find(current_v);
        if (out_edges_it == graph_out.end()) {
            return;
        }
        for (auto p : out_edges_it->second) {
            size_t dst_num = p.first;
            if (passed.find(dst_num) != passed.end()) {
                continue;
            }
            if (all_edges.find({main_v, dst_num}) != all_edges.end()) {
                edges_to_be_removed.insert({main_v, dst_num});
            }
            InnerLevelDfs(main_v,
                          dst_num,
                          graph_out,
                          all_edges,
                          edges_to_be_removed,
                          passed);
        }
    }
}