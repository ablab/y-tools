#pragma once

#include <logger/logger.hpp>
#include "sparse_graph.hpp"

class GraphReader {
    std::string graph_filename;

public:
    GraphReader(std::string graph_filename) {
        this->graph_filename = graph_filename;
    }

    SparseGraphPtr CreateGraph();

private:
    DECL_LOGGER("GraphReader");
};

class GraphWriter {
    std::string graph_filename;

public:
    GraphWriter(std::string graph_filename) {
        this->graph_filename = graph_filename;
    }

    void PrintGraph(SparseGraphPtr graph_ptr);

private:
    DECL_LOGGER("GraphWriter");
};