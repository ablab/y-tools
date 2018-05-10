//
// Created by Andrew Bzikadze on 4/14/17.
//

#include "exporters.hpp"

namespace ig_simulator {

void ExportTree(const Tree& tree, size_t tree_ind, std::ostream& full, std::ostream& included)
{
    const auto& sequences = tree.Sequences();
    for (size_t i = 0; i < sequences.size(); ++i) {
        std::stringstream id_ss;
        id_ss << ">tree_" << tree_ind << "_antibody_" << i;
        std::string id { id_ss.str() };
        full << id << '\n' << sequences[i] << '\n';
        if (tree.IsNodeIncluded(i)) {
            included << id << '\n' << sequences[i] << '\n';
        }
    }
}

void ExportDerivedRepertoire(const DerivedRepertoire& derived_repertoire, std::ostream& full, std::ostream& included) {
    for (size_t i = 0; i < derived_repertoire.size(); ++i) {
        ExportTree(derived_repertoire[i], i, full, included);
    }
}

void ExportEdgeLists(const DerivedRepertoire& derived_repertoire, const IgSimulatorConfig::IOParams::OutputParams& config) {
    std::string path = path::append_path(config.output_dir, config.trees_dir);
    path::make_dir(path);
    for (size_t i = 0; i < derived_repertoire.size(); ++i) {
        std::stringstream filename;
        filename << "tree_" << i << ".dot";
        std::string full_filename = path::append_path(path, filename.str());
        std::ofstream out;
        out.open(full_filename);
        out << derived_repertoire[i];
        out.close();
    }
}

} // End namespace ig_simulator
