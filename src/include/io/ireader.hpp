//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
//todo rename to reader
#pragma once

#include <boost/noncopyable.hpp>

namespace io {

struct ReadStreamStat {
    size_t read_count_;
    size_t max_len_;
    uint64_t total_len_;


    ReadStreamStat(): read_count_(0), max_len_(0), total_len_(0) { }

    void write(std::ostream& stream) const {
        stream.write((const char *) &read_count_, sizeof(read_count_));
        stream.write((const char *) &max_len_, sizeof(max_len_));
        stream.write((const char *) &total_len_, sizeof(total_len_));
    }

    void read(std::istream& stream) {
        stream.read((char *) &read_count_, sizeof(read_count_));
        stream.read((char *) &max_len_, sizeof(max_len_));
        stream.read((char *) &total_len_, sizeof(total_len_));
    }

    template<class Read>
    void increase(const Read& read) {
        size_t len = read.size();

        ++read_count_;
        if (max_len_ < len) {
            max_len_ = len;
        }
        total_len_ += read.nucl_count();
    }

    void merge(const ReadStreamStat& stat) {
        read_count_ += stat.read_count_;
        if (max_len_ < stat.max_len_) {
            max_len_ = stat.max_len_;
        }
        total_len_ += stat.total_len_;
    }

    bool valid() const {
        return read_count_ != 0;
    }

};

/**
 * Reader is the interface for all other readers and reader wrappers.
 */
template<typename ReadType>
class ReadStream: boost::noncopyable {
 public:
  typedef ReadType ReadT;

  /*
   * Default destructor.
   */
  virtual ~ReadStream() {}

  /*
   * Check whether the stream is opened.
   *
   * @return true if the stream is opened and false otherwise.
   */
  virtual bool is_open() = 0;

  /*
   * Check whether we've reached the end of stream.
   *
   * @return true if the end of the stream is reached and false
   * otherwise.
   */
  virtual bool eof() = 0;

  /*
   * Read SingleRead or PairedRead from stream (according to ReadType).
   *
   * @param read The SingleRead or PairedRead that will store read data.
   *
   * @return Reference to this stream.
   */
  virtual ReadStream& operator>>(ReadType& read) = 0;

  /*
   * Close the stream.
   */
  virtual void close() = 0;

  /*
   * Close the stream and open it again.
   */
  virtual void reset() = 0;

  virtual ReadStreamStat get_stat() const = 0;

};

template<class Read>
class PredictableReadStream: public ReadStream<Read> {
public:
    virtual size_t size() const = 0;
};

}
