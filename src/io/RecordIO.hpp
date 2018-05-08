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


        class RecordReaderHelper : public RecordReader {
        public:

            RecordSet ReadRecords(std::istream& input_stream) const override;

        protected:

            virtual bool AddRecord(RecordSet&, size_t,
                                   const seqan::CharString&, const seqan::Dna5String&, const seqan::CharString&) const = 0;
        };

        // Use for compatibility reasons only. We try to get rid of this format.
        class IgrecClusterSizeRecordIO : public RecordReaderHelper, public RecordWriter {
        public:

            void WriteRecords(std::ostream& output_stream, const RecordSet& records) const override;

        protected:
            bool AddRecord(RecordSet& record_set, size_t idx, const seqan::CharString& header,
                           const seqan::Dna5String& sequence, const seqan::CharString& quality) const override;

        private:
            template <typename RecordHelper>
            static bool ParseMeta(const seqan::CharString& meta, RecordHelper& record);
        };

        class BarcodeInMetaRecordReader : public RecordReaderHelper {
        public:

            explicit BarcodeInMetaRecordReader(std::vector<std::string> barcode_marker = { "umi", "barcode" })
                    : barcode_marker_(std::move(barcode_marker)) {}

        protected:
            bool AddRecord(RecordSet& record_set, size_t idx, const seqan::CharString& header,
                           const seqan::Dna5String& sequence, const seqan::CharString& quality) const override;

        private:

            const std::vector<std::string> barcode_marker_;
        };
    }
}
