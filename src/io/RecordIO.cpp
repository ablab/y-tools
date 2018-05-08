#include <path_helper.hpp>
#include "RecordIO.hpp"
#include "../ig_tools/utils/string_tools.hpp"
#include <seqan/seq_io.h>

namespace YTools {
    namespace IO {

        void IgrecClusterSizeRecordIO::WriteRecords(std::ostream& output_stream, const RecordSet& records) const {
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

        RecordSet RecordReaderHelper::ReadRecords(std::istream& input_stream) const {
            std::vector<seqan::CharString> read_headers;
            std::vector<seqan::Dna5String> read_seqs;
            std::vector<seqan::CharString> read_quals;
            seqan::SeqFileIn file(input_stream);
            seqan::readRecords(read_headers, read_seqs, read_quals, file);

            RecordSet result(read_headers.size());
            for (size_t i = 0; i < read_seqs.size(); i ++) {
                AddRecord(result, i, read_headers[i], read_seqs[i], read_quals[i]);
            }
            return result;
        }

        template <typename RecordHelper>
        bool IgrecClusterSizeRecordIO::ParseMeta(const seqan::CharString& meta, RecordHelper& record) {
            auto meta_s = seqan_string_to_string(meta);
            auto parts = split(meta_s, "___");
            if (parts.size() != 4) return false;
            if (parts[0] != "cluster") return false;
            if (parts[2] != "size") return false;
            record.Set(CLUSTER_ID, parts[1]);
            const auto size = try_string_to_number<size_t>(parts[3]);
            if (!size) return false;
            record.Set(CLUSTER_SIZE, size.get());
            return true;
        }

        bool IgrecClusterSizeRecordIO::AddRecord(RecordSet& record_set, const size_t idx, const seqan::CharString& header,
                                                 const seqan::Dna5String& sequence, const seqan::CharString& quality) const {
            auto record = record_set[idx].Create().Set(SEQUENCE, sequence);
            if (length(quality) > 0) {
                record.Set(QUALITY, quality);
            }
            if (!ParseMeta(header, record)) return false;
            record_set[idx] = record.Build();
            return true;
        };

        bool BarcodeInMetaRecordReader::AddRecord(RecordSet& record_set, const size_t idx, const seqan::CharString& header,
                                                 const seqan::Dna5String& sequence, const seqan::CharString& quality) const {
            auto record = record_set[idx].Create().Set(SEQUENCE, sequence);
            if (length(quality) > 0) {
                record.Set(QUALITY, quality);
            }
            std::vector<std::string> split_by_marker;
            for (const auto& marker : barcode_marker_) {
                const auto s = seqan_string_to_string(header);
                split_by_marker = split_ignore_case(s, marker);
                if (split_by_marker.size() == 2) break;
            }
            if (split_by_marker.size() != 2) return false;

            std::string meta = split_by_marker[0];
            boost::algorithm::trim(meta);
            std::string umi_info = split_by_marker[1].substr(1);
            if (umi_info.empty()) return false;
            size_t colon = umi_info.find(':');
            if (colon == std::string::npos) {
                record.Set(BARCODE, umi_info);
            } else {
                auto umi = umi_info.substr(0, colon);
                auto qual = umi_info.substr(colon + 1);
                VERIFY_MSG(umi.length() == qual.length(), "UMI and its quality are of different lengths: " << umi_info);
                record.Set(BARCODE, umi).Set(BARCODE_QUALITY, qual);
            }
            record_set[idx] = record.Build();
            return true;
        };

    }
}