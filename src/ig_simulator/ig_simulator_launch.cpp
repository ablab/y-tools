//
// Created by Andrew Bzikadze on 3/15/17.
//

#include <chrono>

#include <clonal_trees/tree_creator/tree_creator.hpp>
#include <clonal_trees/tree_creator/exporters.hpp>
#include "ig_simulator_launch.hpp"
#include "base_repertoire/base_repertoire_simulator.hpp"

using namespace germline_utils;

namespace ig_simulator {

germline_utils::ChainType IgSimulatorLaunch::GetLaunchChainType() const {
    auto v_chain_type = germline_utils::LociParam::ConvertIntoChainTypes(config_.germline_params.loci);
    VERIFY_MSG(v_chain_type.size() == 1, "Only specific chain type is allowed");
    return v_chain_type[0];
}

std::vector<germline_utils::CustomGeneDatabase>
IgSimulatorLaunch::GetDB(const germline_utils::ChainType chain_type) const
{
    GermlineDbGenerator db_generator(config_.io_params.input_params.germline_input,
                                     config_.germline_params);
    INFO("Generation of DB for variable segments...");
    germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
    INFO("Generation of DB for diversity segments...");
    germline_utils::CustomGeneDatabase d_db = db_generator.GenerateDiversityDb();
    INFO("Generation of DB for join segments...");
    germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();

    std::vector<germline_utils::CustomGeneDatabase> db;
    db.emplace_back(std::move(v_db));
    if (chain_type.IsVDJ())
        db.emplace_back(std::move(d_db));
    db.emplace_back(std::move(j_db));
    return db;
}

BaseRepertoire
IgSimulatorLaunch::GetBaseRepertoire(const germline_utils::ChainType chain_type,
                                     std::vector<germline_utils::CustomGeneDatabase>& db) const
{
    INFO("== Base Repertoire starts ==");
    BaseRepertoireSimulator base_repertoire_simulator{config_.simulation_params.base_repertoire_params,
                                                      chain_type,
                                                      db};
    auto base_repertoire =
        base_repertoire_simulator.Simulate(config_.simulation_params.base_repertoire_params.number_of_metaroots);
    std::ofstream base_repertoire_fasta;
    std::ofstream base_repertoire_info;
    base_repertoire_fasta.open(path::append_path(config_.io_params.output_params.output_dir,
                                                 config_.io_params.output_params.base_repertoire_filename));
    base_repertoire_info.open(path::append_path(config_.io_params.output_params.output_dir,
                                                config_.io_params.output_params.base_repertoire_info));
    print_base_repertoire(base_repertoire, base_repertoire_fasta, base_repertoire_info);
    base_repertoire_fasta.close();
    base_repertoire_info.close();
    INFO("== Base Repertoire ends ==");
    return base_repertoire;
}

DerivedRepertoire IgSimulatorLaunch::__GetDerivedRepertoire(const BaseRepertoire &base_repertoire,
                                                            const PoolManagerStrategy pool_manager_strategy) const {
    INFO("== Derived repertoire generation starts ==");
    const auto& vjf_config = config_.simulation_params.base_repertoire_params.metaroot_simulation_params.
                             cdr_labeler_config.vj_finder_config;

    DerivedRepertoireCreator derived_repertoire_creator(vjf_config,
                                                        config_.simulation_params.clonal_tree_simulator_params);
    auto derived_repertoire = derived_repertoire_creator.GenerateDerivedRepertoire(base_repertoire,
                                                                                   pool_manager_strategy);
    INFO("== Derived repertoire generation ends ==");

    INFO("== Derived repertoire export starts ==");
    INFO("== Full and filtered pool export start");
    std::ofstream full, included;
    full.open(path::append_path(config_.io_params.output_params.output_dir,
                                config_.io_params.output_params.full_pool));
    included.open(path::append_path(config_.io_params.output_params.output_dir,
                                    config_.io_params.output_params.filtered_pool));
    DerivedRepertoireExporter(derived_repertoire, full, included);
    full.close();
    included.close();
    INFO("== Full and filtered pool export ends");

    if (!config_.io_params.output_params.trees_dir.empty()) {
        INFO("'" << config_.io_params.output_params.trees_dir << "'")
        INFO("== Edge lists export starts");
        EdgeListsExporters(derived_repertoire, config_.io_params.output_params);
        INFO("== Edge lists export ends");
    }
    INFO("== Derived repertoire export ends ==");
    return derived_repertoire;
}

DerivedRepertoire IgSimulatorLaunch::GetDerivedRepertoire(const BaseRepertoire &base_repertoire) const
{
    const auto pool_manager_strategy = config_.simulation_params.clonal_tree_simulator_params.pool_manager_strategy;
    return __GetDerivedRepertoire(base_repertoire, pool_manager_strategy);
}

void IgSimulatorLaunch::Setup() const {
    const auto& policy = config_.simulation_params.random_seed_policy;
    if (policy.reset) {
        MTSingleton::SetSeed(policy.seed);
    }
}

void IgSimulatorLaunch::Run() {
    Setup();

    INFO("== IgSimulator starts ==");

    germline_utils::ChainType chain_type = GetLaunchChainType();
    std::vector<germline_utils::CustomGeneDatabase> db { GetDB(chain_type) };

    const BaseRepertoire base_repertoire = GetBaseRepertoire(chain_type, db);
    const DerivedRepertoire derived_repertoire = GetDerivedRepertoire(base_repertoire);

    INFO("== IgSimulator ends ==");
}

} // End namespace ig_simulator
