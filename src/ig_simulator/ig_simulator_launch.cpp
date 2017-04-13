//
// Created by Andrew Bzikadze on 3/15/17.
//

#include <chrono>

#include <clonal_trees/tree_creator/tree_creator.hpp>
#include "ig_simulator_launch.hpp"
#include "germline_utils/germline_db_generator.hpp"
// #include "random_generator.hpp"
// #include "gene_chooser/uniform_gene_chooser.hpp"
// #include "nucleotides_remover/uniform_nucleotides_remover.hpp"
// #include "p_nucleotides_creator/uniform_nucleotides_creator.hpp"
// #include "n_nucleotides_inserter/uniform_n_nucleotides_inserter.hpp"
// #include "metaroot_creator/metaroot_creator.hpp"
#include "base_repertoire/gene_chooser/config_based_getter.hpp"
#include "base_repertoire/n_nucleotides_inserter/config_based_getter.hpp"
#include "base_repertoire/nucleotides_remover/config_based_getter.hpp"
#include "base_repertoire/p_nucleotides_creator/config_based_getter.hpp"
#include "base_repertoire/base_repertoire_simulator.hpp"
#include "clonal_trees/tree_creator/pool_manager.hpp"
#include "clonal_trees/tree_creator/shm_creator.hpp"
#include "clonal_trees/tree_creator/tree_size_generator.hpp"
#include "clonal_trees/tree_creator/tree_creator.hpp"
#include "clonal_trees/tree_creator/forest_creator.hpp"
#include "clonal_trees/tree_creator/pool_manager.hpp"

using namespace germline_utils;

namespace ig_simulator {

germline_utils::ChainType IgSimulatorLaunch::GetLaunchChainType() const {
    auto v_chain_type = germline_utils::LociParam::ConvertIntoChainTypes(config_.algorithm_params.germline_params.loci);
    VERIFY_MSG(v_chain_type.size() == 1, "Only specific chain type is allowed");
    return v_chain_type[0];
}

void IgSimulatorLaunch::Run() {
    std::cout << config_.simulation_params.base_repertoire_params.
        metaroot_simulation_params.nucleotides_remover_params.
        uniform_remover_params.max_remove_v_gene << std::endl;

    MTSingleton::SetSeed(1);
    INFO("== IgSimulator starts ==");

    germline_utils::ChainType chain_type = GetLaunchChainType();

    GermlineDbGenerator db_generator(config_.io_params.input_params.germline_input,
                                     config_.algorithm_params.germline_params);
    INFO("Generation of DB for variable segments...");
    germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
    INFO("Generation of DB for diversity segments...");
    germline_utils::CustomGeneDatabase d_db = db_generator.GenerateDiversityDb();
    INFO("Generation of DB for join segments...");
    germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();

    std::vector<germline_utils::CustomGeneDatabase*> db;
    db.push_back(&v_db);
    if (chain_type.IsVDJ())
        db.push_back(&d_db);
    db.push_back(&j_db);

    BaseRepertoireSimulator base_repertoire_simulator{config_.simulation_params.base_repertoire_params,
                                                      chain_type,
                                                      db};
    auto base_repertoire = base_repertoire_simulator.Simulate(100);
    std::ofstream base_repertoire_out;
    base_repertoire_out.open(config_.io_params.output_params.output_dir + "/test.fa");
    base_repertoire_out << base_repertoire;
    base_repertoire_out.close();
    std::cout << config_.io_params.output_params.output_dir + "/test.fa\n";

    UniformPoolManager manager(0.5);
    // manager.GetIndex();
    // for (size_t i = 0; i < 10000000; i++) {
    //     std::cout << i << " " << manager.GetIndex().second << std::endl;
    // }

    AbstractShmCreatorCPtr shm_creator(new PoissonShmCreator(1));
    AbstractTreeSizeGeneratorCPtr tree_size_generator(new GeometricTreeSizeGenerator(0.001));
    // TreeCreator tree_creator(std::move(shm_creator), std::move(tree_size_generator), 0.5);
    ForestCreator forest_creator(std::move(shm_creator), std::move(tree_size_generator), 0.5);
    //std::cout << tree;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    // auto tree = tree_creator.GenerateTree<DeepTreePoolManager>(base_repertoire[0].MetarootPtr().get());
    auto forest = forest_creator.GenerateForest<DeepTreePoolManager>(base_repertoire[1]);
    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    // std::cout << forest << "\n";

    // auto loci = germline_utils::LociParam::ConvertIntoChainTypes(config_.algorithm_params.germline_params.loci);
    // VERIFY_MSG(loci.size() == 1, "Simulation only one locus");
    // auto locus = loci[0];

    // AbstractVDJGeneChooserCPtr gene_chooser(new UniformVDJGeneChooser(v_db, d_db, j_db));
    // AbstractNucleotidesRemoverCPtr nucl_remover(new UniformNucleotidesRemover());
    // AbstractPNucleotidesCreatorCPtr nucl_creator(new UniformPNucleotidesCreator());
    // AbstractNNucleotidesInserterCPtr nucl_inserter(new UniformNNucleotidesInserter());

    // VDJMetarootCreator metaroot_creator(v_db, d_db, j_db,
    //                                     std::move(gene_chooser),
    //                                     std::move(nucl_remover),
    //                                     std::move(nucl_creator),
    //                                     std::move(nucl_inserter),
    //                                     locus);
    // auto root = metaroot_creator.CreateRoot();
    // std::cout << root << std::endl;
    // std::cout << root.Sequence() << std::endl;

    // MTSingleton::SetSeed(1);
    // auto& g = MTSingleton::GetInstance();
    // auto distr = std::uniform_int_distribution<size_t>();
    // std::cout << distr(g) << std::endl;
    // MTSingleton::SetSeed(1);
    // std::cout << distr(g) << std::endl;
    // MTSingleton::SetSeed(1);
    // auto& g1 = MTSingleton::GetInstance();
    // std::cout << distr(g1) << std::endl;
}

} // End namespace ig_simulator
