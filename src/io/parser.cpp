//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    parser.cpp
 * @author  Mariya Fomkina
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Parser is the parent class for all streams that read data from
 * different file types (fastq, fasta, sam etc).
 * This file contains functions that are used to select exact parser
 * according to extension.
 */

#include "standard_base.hpp"
#include "logger/logger.hpp"
#include "io/parser.hpp"
#include "io/fasta_fastq_gz_parser.hpp"


namespace io {

/*
 * Get extension from filename.
 *
 * @param filename The name of the file to read from.
 *
 * @return File extension (e.g. "fastq", "fastq.gz").
 */
std::string GetExtension(const std::string& filename) {
  std::string name = filename;
  size_t pos = name.find_last_of(".");
  std::string ext = "";
  if (pos != std::string::npos) {
    ext = name.substr(name.find_last_of(".") + 1);
    if (ext == "gz") {
      ext = name.substr(name.find_last_of
                        (".", name.find_last_of(".") - 1) + 1);
    }
  }
  return ext;
}

/*
 * Select parser type according to file extension.
 *
 * @param filename The name of the file to be opened.
 * @param offset The offset of the read quality.

 * @return Pointer to the new parser object with these filename and
 * offset.
 */
Parser* SelectParser(const std::string& filename,
                     OffsetType offset_type /*= PhredOffset*/) {
  std::string ext = GetExtension(filename);
  VERIFY_MSG(ext != "bam", "BAM format is not supported");

  return new FastaFastqGzParser(filename, offset_type);
  /*
  if ((ext == "fastq") || (ext == "fastq.gz") ||
      (ext == "fasta") || (ext == "fasta.gz") ||
      (ext == "fa") || (ext == "fq.gz") ||
      (ext == "fq") || (ext == "fa.gz") ||
      (ext == "seq") || (ext == "seq.gz")) {
    return new FastaFastqGzParser(filename, offset_type);
  }

  ERROR("Unknown file extention in input!");
  return NULL; */
}

void first_fun(int) {
}

}
