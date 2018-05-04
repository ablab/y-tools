#include "Record.hpp"

namespace YTools {
    namespace IO {
        namespace internal {
            RecordHelper& RecordHelper::Set(const std::string &field_name, const FieldValueType& value) {
                fields_[field_name] = value;
                return *this;
            }

            Record RecordHelper::Build() {
                return Record(fields_);
            }

            std::vector<std::string> Record::GetFieldNames() const {
                std::vector<std::string> result;
                result.reserve(fields_.size());
                for (const auto& entry : fields_) {
                    result.push_back(entry.first);
                }
                return result;
            }
        }
    }
}
