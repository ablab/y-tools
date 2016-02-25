#pragma once

#include <cassert>
#include <vector>
#include <stdexcept>
#include <memory>

using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include <boost/program_options.hpp>
#include "fast_ig_tools.hpp"
using path::make_dirs;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

#include "ig_block_alignment.hpp"

namespace fast_ig_tools {
using std::string;

struct GermlineLociVJDB {
    std::vector<Dna5String> v_reads;
    std::vector<CharString> v_ids;
    std::vector<Dna5String> j_reads;
    std::vector<CharString> j_ids;
};

class GermlineDB {
public:
    std::vector<GermlineLociVJDB> locus_databases;
    GermlineLociVJDB all_loci;

    std::vector<std::string> locus_names;
    const std::string db_directory;
    const std::string organism;
    bool pseudogenes;
    std::vector<Dna5String> &all_v_reads = all_loci.v_reads;
    std::vector<Dna5String> &all_j_reads = all_loci.j_reads;
    std::vector<CharString> &all_v_ids = all_loci.v_ids;
    std::vector<CharString> &all_j_ids = all_loci.j_ids;
    std::vector<size_t> locus_index;

    std::unique_ptr<BlockAligner> valigner;
    std::unique_ptr<BlockAligner> jaligner;
    std::vector<std::unique_ptr<BlockAligner>> jaligners;


    static std::vector<std::string> expand_loci(const std::vector<std::string> &loci) {
        std::vector<std::string> IG = { "IGH", "IGK", "IGL" };
        std::vector<std::string> TR = { "TRA", "TRB", "TRD", "TRG" };
        std::vector<std::string> ALL(IG);
        ALL.insert(ALL.end(), TR.cbegin(), TR.cend());

        std::vector<std::string> result;
        for (const std::string &locus : loci) {
            if (locus == "IG") {
                result.insert(result.end(), IG.cbegin(), IG.cend());
            } else if (locus == "TR") {
                result.insert(result.end(), TR.cbegin(), TR.cend());
            } else if (locus == "all") {
                result.insert(result.end(), ALL.cbegin(), ALL.cend());
            } else if (std::find(ALL.cbegin(), ALL.cend(), locus) != ALL.cend()) {
                result.push_back(locus);
            } else {
                throw std::invalid_argument("Invalid locus name: " + locus);
            }
        }

        // Remove duplicates
        remove_duplicates(result);

        return result;
    }

    std::string gene_file_name(const std::string &locus,
                               const std::string &gene,
                               bool pseudo = false) const {
        return db_directory + "/" + organism + "/" + locus + gene + (pseudo ? "-allP" : "") + ".fa";
    }

    template<typename Tparam>
    GermlineDB(const Tparam &param) : db_directory{param.db_directory},
                                      organism{param.organism},
                                      pseudogenes{param.pseudogenes},
                                      locus_names{expand_loci({ param.loci })} {
        for (const std::string &locus : locus_names) {
            GermlineLociVJDB db;

            {
                std::string v_file = gene_file_name(locus, "V", pseudogenes);
                SeqFileIn seqFileIn(v_file.c_str());
                readRecords(db.v_ids, db.v_reads, seqFileIn);
            }

            {
                std::string j_file = gene_file_name(locus, "J", pseudogenes);
                SeqFileIn seqFileIn(j_file.c_str());
                readRecords(db.j_ids, db.j_reads, seqFileIn);
            }

            locus_databases.push_back(db);
        }

        // Join V genes
        for (size_t i = 0; i < locus_databases.size(); ++i) {
            const auto &db = locus_databases[i];

            all_v_reads.insert(all_v_reads.end(),
                               db.v_reads.cbegin(), db.v_reads.cend());
            all_v_ids.insert(all_v_ids.end(),
                               db.v_ids.cbegin(), db.v_ids.cend());

            all_j_reads.insert(all_j_reads.end(),
                               db.j_reads.cbegin(), db.j_reads.cend());
            all_j_ids.insert(all_j_ids.end(),
                               db.j_ids.cbegin(), db.j_ids.cend());

            for (const auto &_ : db.v_reads) {
                locus_index.push_back(i);
            }
        }

        valigner.reset(new BlockAligner(all_v_reads, param.K, param.max_global_gap, param.left_uncoverage_limit,
                                        100500,
                                        param.max_local_insertions, param.max_local_deletions, param.min_k_coverage));
        jaligner.reset(new BlockAligner(all_j_reads, param.word_size_j,
                                        param.max_global_gap, 100000, 10000,
                                        param.max_local_insertions, param.max_local_deletions, param.min_k_coverage_j));

        for (const auto db : locus_databases) {
            BlockAligner *p = new BlockAligner(all_j_reads, param.word_size_j,
                                               param.max_global_gap, 100000, 10000,
                                               param.max_local_insertions, param.max_local_deletions, param.min_k_coverage_j);
            jaligners.push_back(std::unique_ptr<BlockAligner>(p));
        }
    }

    struct VJAlignment {
        size_t locus_id;
        std::string locus;

        int strand;
        Dna5String read;
        GermlineLociVJDB *vbase, *jbase;

        std::vector<BlockAligner::Alignment> v_hits, j_hits;

        bool empty() const {
            return v_hits.empty() || j_hits.empty();
        }

        // add fix
        // add align
        // add various getters
    };


    template<typename T>
    static bool all_equal(const T &v) {
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] != v[0]) {
                return false;
            }
        }

        return true;
    }

    template<typename Titer, typename Tf>
    static auto max_map(Titer b, Titer e, const Tf &f) -> decltype(f(*b)) { // TODO Add decay
        assert(b != e);

        auto result = f(*b++);
        for(;b != e; ++b) {
            auto cur = f(*b);
            if (result <  cur) {
                result = cur;
            }
        }

        return result;
    }

    template<typename Tparam>
    void query(const Dna5String &read,
               bool fix_strand, // TODO Naming?
               const Tparam &param) const {
        // assert(limit_j > 0);
        // if (limit_j == 0) {
        //     limit_j = limit_v;
        // }

        VJAlignment result_alignment;
        // Fix strand if asked

        Dna5String read_rc = read;
        reverseComplement(read_rc);

        auto result_pstrand = valigner->query(read, param.max_candidates);
        auto result_nstrand = valigner->query(read_rc, param.max_candidates);

        int pscore = (result_pstrand.size() > 0) ? result_pstrand[0].kp_coverage : 0;
        int nscore = (result_nstrand.size() > 0) ? result_nstrand[0].kp_coverage : 0;

        int strand = result_alignment.strand = (pscore >= nscore) ? 1 : -1;
        const auto &stranded_read = (strand == 1) ? read : read_rc;
        const auto &result = (strand == 1) ? result_pstrand : result_nstrand;

        result_alignment.read = stranded_read;
        result_alignment.vbase = &all_loci;
        // Aling V --- already done


        // Identify locus, in case of misconsensus report empty alignment
        std::vector<size_t> locus_ids;
        locus_ids.reserve(result.size());
        for (const auto &align : result) {
            size_t loc_ind = locus_index[align.needle_index];

            locus_ids.push_back(loc_ind);
        }

        // Check equality
        if (!all_equal(locus_ids)) {
            return result_alignment; // Empty TODO Add marker
        }

        // Save V
        result_alignment.v_hits = result;

        size_t locus_id = result_alignment.locus_id = locus_ids[0];
        result_alignment.locus = locus_names[locus_id];
        result_alignment.vbase = locus_databases[locus_id];

        // Aling V --- already done
        // Find minimum suffix after V gene
        int end_of_v = max_map(result.cbegin(), result.cend(),
                               [](const decltype(result) &align) -> int { return align().last_match().read_pos + align().last_match().length; });

        // Align J
        auto jaligns = jaligners[locus_id]->query(seqan::suffix(read, end_of_v), param.max_candidates_j);

        result_alignment.j_hits = jaligns;

        // Report
        return result_alignment;
    }
};



} // namespace fast_ig_tools

// vim: ts=4:sw=4
