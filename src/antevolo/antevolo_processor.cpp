#include <logger/logger.hpp>

#include "antevolo_processor.hpp"
#include "clone_set_decomposers/v_clone_set_decomposer.hpp"
#include "vj_class_processors/vj_class_processor.hpp"


namespace antevolo {

    EvolutionaryTreeStorage AntEvoloProcessor::JoinEvolutionaryStoragesFromThreads() {
        EvolutionaryTreeStorage resulting_tree_storage ;
        for(auto it = thread_tree_storages_.begin(); it != thread_tree_storages_.end(); it++) {
            resulting_tree_storage.AppendArchive(*it);
        }
        auto const& tree = *resulting_tree_storage.cbegin();
        return resulting_tree_storage;
    }

    EvolutionaryTreeStorage AntEvoloProcessor::Process() {

        VCloneSetDecomposer clone_set_decomposer(clone_set_); // storage for reconstructed fake vertices
        auto v_decomposition = clone_set_decomposer.CreateDecomposition();
        INFO("V decomposition containing " << v_decomposition.Size() << " classes was created.");
        INFO("Largest class contains " << v_decomposition.MaxClassSize() << " clone(s)");
        omp_set_num_threads(config_.run_params.num_threads);
        INFO("Construction of clonal trees starts");

        /*
	std::ifstream in_genes("/home/aslabodkin/Ð");
        std::map<std::string, std::pair<std::string, std::string>> gene_names;
        for (size_t i = 0; i < 31055; ++i) {
            std::string read_name;
            std::string v_gene;
            std::string j_gene;
            in_genes >> read_name >> v_gene >> j_gene;
            gene_names[read_name] = std::make_pair(v_gene, j_gene);
        }
        in_genes.close();
        size_t check = 0;
	for (auto it = clone_set_.begin(); it != clone_set_.end(); it++) {
            std::string v_gene;
            for (size_t i = 0; i < seqan::length(it->VGene().name()); ++i) {
                v_gene.push_back(it->VGene().name()[i]);
            }
            std::string j_gene;
            for (size_t i = 0; i < seqan::length(it->JGene().name()); ++i) {
                j_gene.push_back(it->JGene().name()[i]);
            }
            if (gene_names.find(it->Read().name) == gene_names.end()) {
                ++check;
                gene_names[it->Read().name] = std::make_pair(v_gene, j_gene);
            }
        }
        INFO(check);*/
        std::vector<std::map<size_t, size_t>> cdr3_dist_map(config_.run_params.num_threads);
        std::vector<std::map<size_t, size_t>> cdr3_dist_percent_map(config_.run_params.num_threads);
        std::vector<std::map<std::pair<size_t, size_t>, size_t>> cdr3_dist_cdr3_len_map(config_.run_params.num_threads);
        std::vector<std::map<std::pair<size_t, size_t>, size_t>> all_pairs_cdr3_dist_cdr3_len_map(config_.run_params.num_threads);
        std::vector<std::map<size_t, size_t>> all_pairs_cdr3_dist_map(config_.run_params.num_threads);


#pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < v_decomposition.Size(); i++) {
            auto edge_constructor = std::shared_ptr<EvolutionaryEdgeConstructor>(
                    new VJEvolutionaryEdgeConstructor(config_.algorithm_params.edge_construction_params));
            size_t thread_id = omp_get_thread_num();
            auto v_class = v_decomposition.GetClass(i);
            for (auto it1 = v_class.begin(); it1 != v_class.end(); it1++) {
                for (auto it2 = it1; it2 != v_class.end(); it2++) {
                    if ((it1 == it2 || clone_set_[*it1].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3)
                         || clone_set_[*it2].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3)) ||
                        std::min(clone_set_[*it1].VSHMs().size() + clone_set_[*it1].JSHMs().size(),
                                 clone_set_[*it2].VSHMs().size() + clone_set_[*it2].JSHMs().size()) < 0) {
                        continue;
                    }
		    // INFO("here")

                    auto edge = edge_constructor->ConstructEdge(clone_set_[*it1],
                                                                clone_set_[*it2],
                                                                *it1,
                                                                *it2);
                    auto edge_r = edge_constructor->ConstructEdge(clone_set_[*it2],
                                                                  clone_set_[*it1],
                                                                  *it2,
                                                                  *it1);
                    // if (clone_set_[*it1].CDR3Range().length() == clone_set_[*it2].CDR3Range().length()) {
                    //     size_t cdr3_dist = HammingDistance(clone_set_[*it1].CDR3(), clone_set_[*it2].CDR3());
                    //     size_t cdr3_len = clone_set_[*it1].CDR3Range().length();
                    //     auto p = std::make_pair(cdr3_dist, cdr3_len);
                    //     if (all_pairs_cdr3_dist_cdr3_len_map[thread_id].find(p) == all_pairs_cdr3_dist_cdr3_len_map[thread_id].end()) {
                    //         all_pairs_cdr3_dist_cdr3_len_map[thread_id][p] = 0;
                    //     }
                    //     all_pairs_cdr3_dist_cdr3_len_map[thread_id][p] += 1;

                    //     if (all_pairs_cdr3_dist_map[thread_id].find(cdr3_dist) == all_pairs_cdr3_dist_map[thread_id].end()) {
                    //         all_pairs_cdr3_dist_map[thread_id][cdr3_dist] = 0;
                    //     }
                    //     all_pairs_cdr3_dist_map[thread_id][cdr3_dist] += 1;

                    // }
                    if ((edge->IsDirected() || edge_r->IsDirected() || edge->IsUndirected()) &&
                        clone_set_[*it1].CDR3Range().length() == clone_set_[*it2].CDR3Range().length()) {
                        size_t cdr3_dist = edge->CDR3Distance();
                        VERIFY_MSG(cdr3_dist == HammingDistance(clone_set_[*it1].CDR3(), clone_set_[*it2].CDR3()),
                                   cdr3_dist << " " << HammingDistance(clone_set_[*it1].CDR3(), clone_set_[*it2].CDR3()));
                        size_t cdr3_len = clone_set_[*it1].CDR3Range().length();
                        size_t percent = static_cast<double>(cdr3_dist) * 100 / cdr3_len;
                        if (cdr3_dist_map[thread_id].find(cdr3_dist) == cdr3_dist_map[thread_id].end()) {
                            cdr3_dist_map[thread_id][cdr3_dist] = 0;
                        }
                        cdr3_dist_map[thread_id][cdr3_dist] += 1;

                        if (cdr3_dist_percent_map[thread_id].find(percent) == cdr3_dist_percent_map[thread_id].end()) {
                            cdr3_dist_percent_map[thread_id][percent] = 0;
                        }
                        cdr3_dist_percent_map[thread_id][percent] += 1;

                        auto p = std::make_pair(cdr3_dist, cdr3_len);
                        if (cdr3_dist_cdr3_len_map[thread_id].find(p) == cdr3_dist_cdr3_len_map[thread_id].end()) {
                            cdr3_dist_cdr3_len_map[thread_id][p] = 0;
                        }
                        ++cdr3_dist_cdr3_len_map[thread_id][p];
//                        if (HammingDistance(clone_set_[*it1].VGene().name(), clone_set_[*it2].VGene().name()) == 0 &&
//                            HammingDistance(clone_set_[*it1].JGene().name(), clone_set_[*it2].JGene().name()) == 0)
                    }
                }
            }
        }
        std::map<size_t, size_t> cdr3_dist_;
        std::map<size_t, size_t> cdr3_dist_percent_;
        std::map<std::pair<size_t, size_t>, size_t> cdr3_dist_cdr3_len_;
        std::map<std::pair<size_t, size_t>, size_t> all_pairs_cdr3_dist_cdr3_len_;
        std::map<size_t, size_t> all_pairs_cdr3_dist_;
        for (size_t i = 0; i < cdr3_dist_map.size(); ++i) {
            for (auto p : cdr3_dist_map[i]) {
                if (cdr3_dist_.find(p.first) == cdr3_dist_.end()) {
                    cdr3_dist_.insert(p);
                }
                else {
                    cdr3_dist_[p.first] += p.second;
                }
            }
            for (auto p : cdr3_dist_percent_map[i]) {
                if (cdr3_dist_percent_.find(p.first) == cdr3_dist_percent_.end()) {
                    cdr3_dist_percent_.insert(p);
                }
                else {
                    cdr3_dist_percent_[p.first] += p.second;
                }
            }
            for (auto p : cdr3_dist_cdr3_len_map[i]) {
                if (cdr3_dist_cdr3_len_.find(p.first) == cdr3_dist_cdr3_len_.end()) {
                    cdr3_dist_cdr3_len_.insert(p);
                }
                else {
                    cdr3_dist_cdr3_len_[p.first] += p.second;
                }
            }
            for (auto p : all_pairs_cdr3_dist_cdr3_len_map[i]) {
                if (all_pairs_cdr3_dist_cdr3_len_.find(p.first) == all_pairs_cdr3_dist_cdr3_len_.end()) {
                    all_pairs_cdr3_dist_cdr3_len_.insert(p);
                }
                else {
                    all_pairs_cdr3_dist_cdr3_len_[p.first] += p.second;
                }
            }
            for (auto p : all_pairs_cdr3_dist_map[i]) {
                if (all_pairs_cdr3_dist_.find(p.first) == all_pairs_cdr3_dist_.end()) {
                    all_pairs_cdr3_dist_.insert(p);
                }
                else {
                    all_pairs_cdr3_dist_[p.first] += p.second;
                }
            }
        }
        std::ofstream out(config_.output_params.output_dir + "/cdr_dist.txt");
        for (auto p : cdr3_dist_) {
            out << p.first << "\t" << p.second << "\n";
        }
        out.close();
        std::ofstream out2(config_.output_params.output_dir + "/cdr_percent.txt");
        for (auto p : cdr3_dist_percent_) {
            out2 << p.first << "\t" << p.second << "\n";
        }
        out2.close();
        std::ofstream out3(config_.output_params.output_dir + "/cdr_dist_len.txt");
        for (auto p : cdr3_dist_cdr3_len_) {
            out3 << p.first.first << "\t" << p.first.second << "\t" << p.second << "\n";
        }
        out3.close();
        std::ofstream out4(config_.output_params.output_dir + "/all_pairs_cdr_dist_len.txt");
        for (auto p : all_pairs_cdr3_dist_cdr3_len_) {
            out4 << p.first.first << "\t" << p.first.second << "\t" << p.second << "\n";
        }
        out4.close();
        std::ofstream out5(config_.output_params.output_dir + "/all_pairs_cdr_dist.txt");
        for (auto p : all_pairs_cdr3_dist_) {
            out5 << p.first << "\t" << p.second << "\n";
        }
        out5.close();
        return EvolutionaryTreeStorage();
    }
}
