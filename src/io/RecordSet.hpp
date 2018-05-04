#pragma once

#include "Record.hpp"

namespace YTools {
    namespace IO {

        class RecordSet {
        public:

            template <typename R, typename Itr>
            class Iterator_;

            using Iterator = Iterator_<internal::Record, std::vector<internal::Record>::iterator>;
            using ConstIterator = Iterator_<const internal::Record, std::vector<internal::Record>::const_iterator>;

            explicit RecordSet(size_t size);

            internal::Record& operator[](size_t i);

            const internal::Record& operator[](size_t i) const;

            size_t size() const {
                return records_.size();
            }

            Iterator begin() {
                return Iterator(records_.begin());
            }

            Iterator end() {
                return Iterator(records_.end());
            }

            ConstIterator begin() const {
                return ConstIterator(records_.cbegin());
            }

            ConstIterator end() const {
                return ConstIterator(records_.cend());
            }

            template <typename R, typename Itr>
            class Iterator_ {
            public:

                explicit Iterator_(const Itr& itr) : itr_(itr) {}

                R& operator*() {
                    return *itr_;
                }

                Iterator_& operator++() {
                    ++ itr_;
                    return *this;
                }

                const Iterator_ operator++(int) {
                    Iterator_ copy(itr_);
                    ++ itr_;
                    return copy;
                }

                bool operator==(const Iterator_& other) const {
                    return itr_ == other.itr_;
                }

                bool operator!=(const Iterator_& other) const {
                    return itr_ != other.itr_;
                }

            private:

                Itr itr_;
            };

        private:

            std::vector<internal::Record> records_;
        };

    }
}
