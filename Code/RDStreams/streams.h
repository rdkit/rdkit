//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_STREAMS_H
#define RD_STREAMS_H

#include <istream>
#include <memory>
#include <string>

namespace RDKit {

//! Input stream that transparently decompresses a gzip-compressed file.
/*!
  Backed by zlib. The decompressing stream buffer is held by \c d_buf so that
  zlib stays out of this header and callers only ever see a std::istream.
*/
class RDKIT_RDSTREAMS_EXPORT gzstream : public std::istream {
 public:
  explicit gzstream(const std::string &fname);
  ~gzstream() override;

  gzstream(const gzstream &) = delete;
  gzstream &operator=(const gzstream &) = delete;

 private:
  // owns the decompressing streambuf, which in turn owns the underlying file
  std::unique_ptr<std::streambuf> d_buf;
};

#ifdef RDK_USE_BZIP2
//! Input stream that transparently decompresses a bzip2-compressed file.
/*!
  Backed by libbz2, with the decompressing stream buffer held by \c d_buf so
  that bzlib stays out of this header (see \c gzstream).
*/
class RDKIT_RDSTREAMS_EXPORT bz2stream : public std::istream {
 public:
  explicit bz2stream(const std::string &fname);
  ~bz2stream() override;

  bz2stream(const bz2stream &) = delete;
  bz2stream &operator=(const bz2stream &) = delete;

 private:
  std::unique_ptr<std::streambuf> d_buf;
};
#endif

}  // namespace RDKit
#endif
