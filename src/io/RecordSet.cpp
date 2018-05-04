#include "RecordSet.hpp"

namespace YTools {
    namespace IO {

        RecordSet::RecordSet(size_t size) : records_(size) {}

        internal::Record& RecordSet::operator[](size_t i) {
            return records_[i];
        }

        const internal::Record& RecordSet::operator[](size_t i) const {
            return records_[i];
        }

    }
}