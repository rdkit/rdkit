//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "streams.h"

#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <zlib.h>
#ifdef RDK_USE_BZIP2
#include <bzlib.h>
#endif

namespace RDKit {
namespace {

constexpr std::size_t STREAM_BUF_SIZE = 16384;

//! a std::streambuf that inflates gzip/zlib data read from a file.
//!
//! The buffer owns the underlying std::ifstream so that a wrapping std::istream
//! only has to manage this single object.
class gzip_streambuf : public std::streambuf {
 public:
  explicit gzip_streambuf(const std::string &fname)
      : d_file(fname, std::ios_base::in | std::ios_base::binary),
        d_inBuf(STREAM_BUF_SIZE),
        d_outBuf(STREAM_BUF_SIZE) {
    std::memset(&d_zs, 0, sizeof(d_zs));
    // 15 window bits + 32 enables automatic zlib/gzip header detection
    if (inflateInit2(&d_zs, 15 + 32) != Z_OK) {
      throw std::runtime_error("could not initialize zlib inflate stream");
    }
    d_inited = true;
    setg(d_outBuf.data(), d_outBuf.data(), d_outBuf.data());
  }
  ~gzip_streambuf() override {
    if (d_inited) {
      inflateEnd(&d_zs);
    }
  }
  bool is_open() const { return d_file.is_open(); }

 protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    if (d_streamEnd) {
      return traits_type::eof();
    }

    d_zs.next_out = reinterpret_cast<Bytef *>(d_outBuf.data());
    d_zs.avail_out = static_cast<uInt>(d_outBuf.size());

    // keep inflating until we produce at least one byte or hit the end
    while (d_zs.avail_out == d_outBuf.size() && !d_streamEnd) {
      if (d_zs.avail_in == 0) {
        d_file.read(d_inBuf.data(), static_cast<std::streamsize>(d_inBuf.size()));
        d_zs.next_in = reinterpret_cast<Bytef *>(d_inBuf.data());
        d_zs.avail_in = static_cast<uInt>(d_file.gcount());
        if (d_zs.avail_in == 0) {
          // input exhausted before a clean end-of-stream marker
          d_streamEnd = true;
          break;
        }
      }
      int ret = inflate(&d_zs, Z_NO_FLUSH);
      if (ret == Z_STREAM_END) {
        if (d_zs.avail_in == 0 && d_file.eof()) {
          d_streamEnd = true;
        } else {
          // a concatenated gzip member follows; reset and keep going
          inflateReset(&d_zs);
        }
      } else if (ret == Z_BUF_ERROR) {
        // no progress possible without more input; loop to read more (or stop)
        if (d_zs.avail_in == 0 && d_file.eof()) {
          d_streamEnd = true;
        }
      } else if (ret != Z_OK) {
        throw std::runtime_error(std::string("zlib inflate error: ") +
                                 (d_zs.msg ? d_zs.msg : "unknown"));
      }
    }

    std::size_t produced = d_outBuf.size() - d_zs.avail_out;
    setg(d_outBuf.data(), d_outBuf.data(), d_outBuf.data() + produced);
    if (produced == 0) {
      return traits_type::eof();
    }
    return traits_type::to_int_type(*gptr());
  }

 private:
  std::ifstream d_file;
  z_stream d_zs;
  std::vector<char> d_inBuf;
  std::vector<char> d_outBuf;
  bool d_inited = false;
  bool d_streamEnd = false;
};

#ifdef RDK_USE_BZIP2
//! a std::streambuf that decompresses bzip2 data read from a file.
class bzip2_streambuf : public std::streambuf {
 public:
  explicit bzip2_streambuf(const std::string &fname)
      : d_file(fname, std::ios_base::in | std::ios_base::binary),
        d_inBuf(STREAM_BUF_SIZE),
        d_outBuf(STREAM_BUF_SIZE) {
    std::memset(&d_bs, 0, sizeof(d_bs));
    if (BZ2_bzDecompressInit(&d_bs, 0, 0) != BZ_OK) {
      throw std::runtime_error("could not initialize bzip2 decompress stream");
    }
    d_inited = true;
    setg(d_outBuf.data(), d_outBuf.data(), d_outBuf.data());
  }
  ~bzip2_streambuf() override {
    if (d_inited) {
      BZ2_bzDecompressEnd(&d_bs);
    }
  }
  bool is_open() const { return d_file.is_open(); }

 protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }
    if (d_streamEnd) {
      return traits_type::eof();
    }

    d_bs.next_out = d_outBuf.data();
    d_bs.avail_out = static_cast<unsigned int>(d_outBuf.size());

    while (d_bs.avail_out == d_outBuf.size() && !d_streamEnd) {
      if (d_bs.avail_in == 0) {
        d_file.read(d_inBuf.data(), static_cast<std::streamsize>(d_inBuf.size()));
        d_bs.next_in = d_inBuf.data();
        d_bs.avail_in = static_cast<unsigned int>(d_file.gcount());
        if (d_bs.avail_in == 0) {
          d_streamEnd = true;
          break;
        }
      }
      int ret = BZ2_bzDecompress(&d_bs);
      if (ret == BZ_STREAM_END) {
        if (d_bs.avail_in == 0 && d_file.eof()) {
          d_streamEnd = true;
        } else {
          // concatenated bzip2 stream: restart the decompressor
          BZ2_bzDecompressEnd(&d_bs);
          if (BZ2_bzDecompressInit(&d_bs, 0, 0) != BZ_OK) {
            throw std::runtime_error("could not reinitialize bzip2 stream");
          }
        }
      } else if (ret != BZ_OK) {
        throw std::runtime_error("bzip2 decompress error");
      }
    }

    std::size_t produced = d_outBuf.size() - d_bs.avail_out;
    setg(d_outBuf.data(), d_outBuf.data(), d_outBuf.data() + produced);
    if (produced == 0) {
      return traits_type::eof();
    }
    return traits_type::to_int_type(*gptr());
  }

 private:
  std::ifstream d_file;
  bz_stream d_bs;
  std::vector<char> d_inBuf;
  std::vector<char> d_outBuf;
  bool d_inited = false;
  bool d_streamEnd = false;
};
#endif

}  // namespace

gzstream::gzstream(const std::string &fname) : std::istream(nullptr) {
  auto buf = std::make_unique<gzip_streambuf>(fname);
  if (!buf->is_open()) {
    setstate(std::ios_base::failbit);
  }
  d_buf = std::move(buf);
  rdbuf(d_buf.get());
}
gzstream::~gzstream() = default;

#ifdef RDK_USE_BZIP2
bz2stream::bz2stream(const std::string &fname) : std::istream(nullptr) {
  auto buf = std::make_unique<bzip2_streambuf>(fname);
  if (!buf->is_open()) {
    setstate(std::ios_base::failbit);
  }
  d_buf = std::move(buf);
  rdbuf(d_buf.get());
}
bz2stream::~bz2stream() = default;
#endif

}  // namespace RDKit
