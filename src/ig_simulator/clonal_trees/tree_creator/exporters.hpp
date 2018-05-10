//
// Created by Andrew Bzikadze on 4/14/17.
//

#pragma once

#include <fstream>
#include <clonal_trees/tree/tree.hpp>
#include "ig_simulator_config.hpp"
#include "derived_repertoire_creator.hpp"

namespace ig_simulator {

void ExportTree(const Tree& tree, size_t tree_ind, std::ostream& full, std::ostream& included);
void ExportDerivedRepertoire(const DerivedRepertoire& derived_repertoire, std::ostream& full, std::ostream& included);

void ExportEdgeLists(const DerivedRepertoire& derived_repertoire, const IgSimulatorConfig::IOParams::OutputParams& config);

} // End namespace ig_simulator
