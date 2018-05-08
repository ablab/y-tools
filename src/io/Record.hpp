#pragma once

#include <string>
#include <map>
#include <boost/variant.hpp>
#include <set>
#include <utility>
#include <seqan/seq_io.h>

namespace YTools{
    namespace IO {
        template <typename T>
        class RecordField;

        namespace internal {
            // Try to keep the list of types short. Add new types only if it provides better performance for critical code.
            typedef boost::variant<size_t, std::string, seqan::Dna5String, seqan::CharString> FieldValueType;

            class Record;

            class RecordHelper {
            public:

                RecordHelper& Set(const std::string& field_name, const FieldValueType& value);

                template <typename T, typename S>
                RecordHelper& Set(const RecordField<T>& field, const S& value) {
                    static_assert(std::is_assignable<T&, S>::value);

                    return Set(field.GetName(), static_cast<T>(value));
                }

                Record Build();

            private:
                std::map<std::string, internal::FieldValueType> fields_;

            };

            // This class should remain an implementation detail in internal namespace.
            // In case of performance problems it should be replaced with more lightweight implementation.
            class Record {
            public:

                Record() = default;
                explicit Record(std::map<std::string, internal::FieldValueType> fields) : fields_(std::move(fields)) {}

                static RecordHelper Create() {
                    return RecordHelper();
                }

                template <typename T>
                T GetField(const std::string& field_name) const {
                    return boost::get<T>(fields_.at(field_name));
                }

                template <typename T>
                T GetField(const RecordField<T>& field) const {
                    return boost::get<T>(fields_.at(field.GetName()));
                }

                template <typename T>
                bool HasField(const std::string& field_name) const {
                    return fields_.find(field_name) != fields_.end();
                }

                template <typename T>
                bool HasField(const RecordField<T>& field) const {
                    return fields_.find(field.GetName()) != fields_.end();
                }

                std::vector<std::string> GetFieldNames() const;

            private:

                std::map<std::string, FieldValueType> fields_;

            };

        }

        class RecordFieldBase {
        public:

            const std::string& GetName() const {
                return name_;
            }

        protected:

            explicit RecordFieldBase(std::string name) : name_(std::move(name)) {}

            const std::string name_;

        private:

            RecordFieldBase() = default;

        };

        template <typename T>
        class RecordField : public RecordFieldBase {
        public:

            explicit RecordField(std::string name) : RecordFieldBase(std::move(name)) {}

        };

        static RecordField<seqan::Dna5String> SEQUENCE("sequence");
        static RecordField<seqan::CharString> QUALITY("quality");
        static RecordField<std::string> CLUSTER_ID("cluster_id");
        static RecordField<size_t> CLUSTER_SIZE("cluster_size");

        static std::vector<RecordFieldBase> PREDEFINED_FIELDS { SEQUENCE, QUALITY, CLUSTER_ID, CLUSTER_SIZE };

    }
}
