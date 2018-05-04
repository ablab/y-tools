#pragma once

#include <vector>
#include "Record.hpp"
#include "RecordSet.hpp"
#include <seqan/seq_io.h>
#include "../ig_tools/utils/string_tools.hpp"

namespace YTools {
    namespace IO {

        class RecordReader {
        public:
            virtual RecordSet ReadRecords(std::istream& input_stream) const = 0;
        };

        class RecordWriter {
        public:
            virtual void WriteRecords(std::ostream& output_stream, const RecordSet& records) const = 0;
        };

        class RecordIO : public RecordReader, public RecordWriter {};

        // Use for compatibility reasons only. We try to get rid of this format.
        class IgrecClusterSizeIO : RecordIO {
        public:

            RecordSet ReadRecords(std::istream& input_stream) const override;

            void WriteRecords(std::ostream& output_stream, const RecordSet& records) const override;

        private:
            template <typename RecordHelper>
            void ParseMeta(const seqan::CharString& meta, RecordHelper& record) const {
                auto meta_s = seqan_string_to_string(meta);
                auto parts = split(meta_s, "___");
                VERIFY(parts.size() == 4);
                VERIFY(parts[0] == "cluster");
                VERIFY(parts[2] == "size");
                record.Set(CLUSTER_ID, parts[1]);
                const auto size = try_string_to_number<size_t>(parts[3]);
                VERIFY(size);
                record.Set(CLUSTER_SIZE, size.get());
            }
        };
    }
}
