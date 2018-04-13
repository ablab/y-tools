#pragma once

#include <base_repertoire/base_repertoire.hpp>
#include "tree_creator.hpp"

namespace ig_simulator {

    using DerivedRepertoire = std::vector<Tree>;

    class DerivedRepertoireCreator {
    private:
        const TreeCreator tree_creator;

    public:
        DerivedRepertoireCreator(const vj_finder::VJFinderConfig& vjf_config,
                                 const ClonalTreeSimulatorParams& config):
                tree_creator(vjf_config, config)
        { }

        DerivedRepertoireCreator(const DerivedRepertoireCreator&) = delete;
        DerivedRepertoireCreator(DerivedRepertoireCreator&&) = delete;
        DerivedRepertoireCreator& operator=(const DerivedRepertoireCreator&) = delete;
        DerivedRepertoireCreator& operator=(DerivedRepertoireCreator&&) = delete;

        DerivedRepertoire GenerateDerivedRepertoire(const BaseRepertoire &repertoire,
                                                    const PoolManagerStrategy pool_manager_strategy) const {
            DerivedRepertoire storage;
            for(const auto& cluster : repertoire) {
                storage.emplace_back(tree_creator.GenerateTree(cluster.MetarootPtr().get(), pool_manager_strategy));
            }
            return storage;
        }
    };

} // End namespace ig_simulator
