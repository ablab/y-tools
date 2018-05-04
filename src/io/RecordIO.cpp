#include <path_helper.hpp>
#include "RecordIO.hpp"
#include "../ig_tools/utils/string_tools.hpp"
#include <seqan/seq_io.h>

namespace YTools {
    namespace IO {

        RecordSet IgrecClusterSizeIO::ReadRecords(std::istream& input_stream) const {
            std::vector<seqan::CharString> read_headers;
            std::vector<seqan::Dna5String> read_seqs;
            std::vector<seqan::CharString> read_quals;
            seqan::SeqFileIn file(input_stream);
            seqan::readRecords(read_headers, read_seqs, read_quals, file);

            RecordSet result(read_headers.size());
            for (size_t i = 0; i < read_seqs.size(); i ++) {
                auto record = result[i].Create().Set(SEQUENCE, read_seqs[i]).Set(QUALITY, read_quals[i]);
                ParseMeta(read_headers[i], record);
                result[i] = record.Build();
            }
            return result;
        }

        void IgrecClusterSizeIO::WriteRecords(std::ostream& output_stream, const RecordSet& records) const {
            const bool has_quality = records.size() == 0 || records[0].HasField(QUALITY);
            if (has_quality) {
                seqan::SeqFileOut file(output_stream, seqan::Fastq());
                for (const auto& record : records) {
                    VERIFY(record.HasField(QUALITY));
                    seqan::writeRecord(
                            file,
                            "cluster___" + record.GetField(CLUSTER_ID) +
                            "___size___" + number_to_string(record.GetField(CLUSTER_SIZE)),
                            record.GetField(SEQUENCE),
                            record.GetField(QUALITY)
                    );
                }
            } else {
                seqan::SeqFileOut file(output_stream, seqan::Fasta());
                for (const auto& record : records) {
                    VERIFY(!record.HasField(QUALITY));
                    seqan::writeRecord(
                            file,
                            "cluster___" + record.GetField(CLUSTER_ID) +
                            "___size___" + number_to_string(record.GetField(CLUSTER_SIZE)),
                            record.GetField(SEQUENCE)
                    );
                }
            }
        }
    }
}