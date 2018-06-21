#pragma once

#include <utils/fastq_reader.hpp>
#include <utils/sequence_tools.hpp>
#include <utils/string_tools.hpp>

#include "logger/log_writers.hpp"

namespace SequenceMerger {

    struct merger_setting {
        size_t min_overlap;
        double max_mismatch_rate;
        bool simulated_mode;

        merger_setting() :
                min_overlap(50),
                max_mismatch_rate(.1),
                simulated_mode(false) {}

        void print();
    };

    class SequenceMerger {
        merger_setting setting_;

        std::string reverse_quality_string(std::string qual);

        std::pair<size_t, size_t> FindBestOverlap(const std::string& seq1, const std::string& seq2);

        std::pair<std::string, std::string> MergeQualifiedSeq(const PairedFastqRead& paired_read);

        std::string MergeNames(size_t index, const std::string& name1, const std::string&);

    public:
        explicit SequenceMerger(merger_setting setting) : setting_(setting) {}

        FastqRead Merge(size_t index, const PairedFastqRead& paired_read);
    };

    class PairedReadsMerger {
        SequenceMerger seq_merger_;
    public:
        explicit PairedReadsMerger(merger_setting settings) : seq_merger_(settings) {}

        std::vector<FastqRead> Merge(const std::vector<PairedFastqRead>& reads);
    };

}