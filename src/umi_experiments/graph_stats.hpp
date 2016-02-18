#pragma once

#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>

struct GraphStats {
public:
    static const size_t LARGE_COMPONENT_THRESHOLD = 100;
    static const string ACGT;

    const SparseGraphPtr graph;
    const size_t size;
    const size_t single;
    const size_t doublets;
    const size_t stars;
    const size_t vert_in_stars;
    const size_t unclassified_vert;
    const vector<vector<size_t> > unclassified_comp;
    const vector<size_t> acgt;
    const vector<size_t> acgt_in_large_components;

    GraphStats(const GraphStats& other) = default;

    GraphStats(const SparseGraphPtr graph, size_t size, size_t single, size_t doublets, size_t stars,
               size_t vert_in_stars, size_t unclassified_vert, const vector<vector<size_t>> unclassified_comp,
               vector<size_t> acgt, vector<size_t> acgt_in_large_components)
            : graph(graph), size(size), single(single), doublets(doublets), stars(stars), vert_in_stars(vert_in_stars),
              unclassified_vert(unclassified_vert), unclassified_comp(unclassified_comp),
              acgt(acgt), acgt_in_large_components(acgt_in_large_components) {
        assert(single + doublets * 2 + vert_in_stars + unclassified_vert == graph->N());
    }

    string ToString(size_t min_unclassified_size) const;

    string ToString() const { return ToString(size + 1); }

    static GraphStats GetStats(const vector<seqan::Dna5String>& reads, const SparseGraphPtr graph);

private:

    static void mark_component(const SparseGraphPtr graph, const size_t v, vector<bool> &visited, vector<size_t> &cc);
    static bool is_star(const SparseGraphPtr  graph, const size_t v);
    static bool is_doublet(const SparseGraphPtr  graph, const size_t v);
};

const string GraphStats::ACGT = "ACGT";

GraphStats GraphStats::GetStats(const vector<seqan::Dna5String>& reads, const SparseGraphPtr  graph) {
    vector<bool> visited(graph->N(), false);
    size_t single = 0;
    size_t stars = 0;
    size_t stars_sum = 0;
    size_t doublets = 0;
    size_t left = 0;
    vector<vector<size_t> > unclassified_comp;
    vector<size_t> acgt;
    vector<size_t> acgt_in_large_components;
    for (size_t i = 0; i < graph->N(); i++) {
        if (visited[i]) continue;

        vector<size_t> component;
        mark_component(graph, i, visited, component);
        size_t size = component.size();

        for (size_t v : component) {
            for (auto nt : reads[v]) {
                auto nt_idx = ACGT.find(nt);
                acgt[nt_idx] ++;
                if (size >= LARGE_COMPONENT_THRESHOLD) {
                    acgt_in_large_components[nt_idx] ++;
                }
            }
        }

        if (graph->Degree(i) == 0) {
            single ++;
            continue;
        }
        if (graph->Degree(i) > 1) {
            if (is_star(graph, i)) {
                stars++;
                stars_sum += size;
                if (size < 2) {
                    ERROR("star of " << size << " vertex");
                }
                continue;
            }
        }
        {
            size_t center = *graph->VertexEdges(i).begin();
            if (graph->Weight()[i] > graph->Weight()[center]) {
                center = i;
            }
            if (is_star(graph, center)) {
                stars ++;
                stars_sum += size;
                continue;
            }
            if (is_doublet(graph, i)) {
                doublets ++;
                continue;
            }
        }

        left += size;
        unclassified_comp.push_back(component);
    }
    return GraphStats(graph, graph->N(), single, doublets, stars, stars_sum, left, unclassified_comp, acgt, acgt_in_large_components);
}

void GraphStats::mark_component(const SparseGraphPtr graph, const size_t v, vector<bool> &visited, vector<size_t> &cc) {
    if (visited[v]) return;
    visited[v] = true;
    cc.push_back(v);
    for (auto u : graph->VertexEdges(v)) {
        mark_component(graph, u, visited, cc);
    }
}

bool GraphStats::is_star(const SparseGraphPtr  graph, const size_t v) {
    bool result = true;
    for (auto u : graph->VertexEdges(v)) {
        result &= graph->Weight()[u] < graph->Weight()[v] && graph->Degree(u) == 1;
    }
    return result;
}

bool GraphStats::is_doublet(const SparseGraphPtr  graph, const size_t v) {
    return graph->Degree(v) == 1 && graph->Degree(*graph->VertexEdges(v).begin()) == 1;
}

string GraphStats::ToString(size_t min_unclassified_size) const {
    stringstream ss;
    ss << boost::format("Found:\ntotal unique UMIs: %d\nsingletons: %d\nstars: %d\ndoublets: %d\nunclassified vertices: %d\naverage star size: %f\naverage unclassified component size: %f\n\n")
            % graph->N() % single % stars % doublets % unclassified_vert
            % (static_cast<double>(vert_in_stars) / static_cast<double>(stars))
            % (static_cast<double>(unclassified_vert) / static_cast<double>(unclassified_comp.size())) << endl;

    ss << "ACGT counts in reads: ";
    for (auto cnt : acgt) {
        ss << cnt << " ";
    }
    ss << endl;
    ss << "ACGT counts in reads forming large components: ";
    for (auto cnt : acgt_in_large_components) {
        ss << cnt << " ";
    }
    ss << endl;

    for (auto component : unclassified_comp) {
        if (component.size() < min_unclassified_size) continue;
        ss << "component vertex weights: ";
        for (auto u : component) {
            ss << graph->Weight()[u] << " ";
        }
        ss << endl << "edges: ";
        for (size_t u = 0; u < component.size(); u ++) {
            for (size_t v = 0; v < u; v ++) {
                if (graph->HasEdge(component[v], component[u])) {
                    ss << "(" << v << ", " << u << ")  ";
                }
            }
        }
        ss << endl;
    }
    return ss.str();
}
