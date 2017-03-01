#include "evolutionary_tree_splitter.hpp"

namespace antevolo {
    EvolutionaryTree ConnectedTreeSplitter::GetTreeByRoot(const EvolutionaryTree& tree,
                                                                           size_t root_id) {
        std::queue<size_t> vertex_queue;
        vertex_queue.push(root_id);
        EvolutionaryTree connected_tree(tree.GetCloneSetPtr());
        while(!vertex_queue.empty()) {
            size_t cur_vertex = vertex_queue.front();
            vertex_queue.pop();
            if(!tree.IsRoot(cur_vertex)) {
                EvolutionaryEdgePtr edge = tree.GetParentEdge(cur_vertex);
                VERIFY_MSG(edge->IsDirected() || edge->IsUndirected(), "Edge from " << edge->DstNum() << " -> " << edge->SrcNum() <<
                                               " is not directed or undirected");
                connected_tree.AddEdge(edge->DstNum(), edge);
            }
            if(!tree.IsLeaf(cur_vertex)) {
                auto outgoing_edges = tree.OutgoingEdges(cur_vertex);
                for(auto it = outgoing_edges.begin(); it != outgoing_edges.end(); it++) {
                    EvolutionaryEdgePtr edge = *it;
                    vertex_queue.push(edge->DstNum());
                }
            }
        }
        return connected_tree;
    }

    // todo: add tree_index (3rd) assigning
    std::vector<EvolutionaryTree> ConnectedTreeSplitter::Split(const EvolutionaryTree& tree) {
        std::vector<EvolutionaryTree> connected_trees;
        if(tree.GetRootNumber() == 1) {
            connected_trees.push_back(tree);
            return connected_trees;
        }
        auto roots = tree.GetRoots();
        for(auto it = roots.begin(); it != roots.end(); it++) {
            EvolutionaryTree connected_tree = GetTreeByRoot(tree, *it);
            connected_trees.push_back(connected_tree);
        }
        VERIFY_MSG(connected_trees.size() == tree.GetRootNumber(), "ERROR: number of connected trees (" <<
                connected_trees.size() << " does not match with number of roots (" << tree.GetRootNumber() <<
                ") in the original forest");
        return connected_trees;
    }
}