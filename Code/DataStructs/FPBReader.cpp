//
// Copyright (c) 2015 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Implementation details here are taken from the file fpb_io.py from chemfp
// (www.chemfp.org)
// Many thanks to Andrew Dalke for creating such great software and for
// helping explain the FPB implementation

#include <Datastructs/ExplicitBitVect.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/StreamOps.h>
#include "FPBReader.h"
#include <boost/cstdint.hpp>

namespace RDKit {

namespace detail {
const unsigned int magicSize = 8;
const std::string FPB_MAGIC("FPB1\r\n\0\0", 8);
const unsigned int tagNameSize = 4;

// the caller is responsible for calling delete[] on `data`
void readChunk(std::istream &istrm, std::string &nm, boost::uint64_t &sz,
               boost::uint8_t *&data) {
  streamRead(istrm, sz);
  char tag[tagNameSize + 1];
  tag[tagNameSize] = 0;
  istrm.read(tag, tagNameSize);
  nm = tag;
  if (sz) {
    data = new boost::uint8_t[sz];
    istrm.read((char *)data, sz);
  } else {
    data = NULL;
  }
  std::cerr << "  CHUNKSZ: " << sz << " name: " << nm << std::endl;
}

struct FPBReader_impl {
  unsigned int len;
  unsigned int nBits;
  // we're assuming that nothing practical has more than 65K bits set
  std::vector<boost::uint16_t> popCounts;
  boost::uint8_t *fpData;
};

}  // end of detail namespace

void FPBReader::init() {
  PRECONDITION(dp_istrm, "no stream");
  dp_impl = new detail::FPBReader_impl;
  // STUB
  dp_impl->len = 100;
  dp_impl->nBits = 2048;

  char magic[detail::magicSize];
  dp_istrm->read(magic, detail::magicSize);
  if (detail::FPB_MAGIC != std::string(magic, detail::magicSize)) {
    throw BadFileException("Invalid FPB magic");
  }
  while (1) {
    if (dp_istrm->eof()) throw BadFileException("EOF hit before FEND record");
    std::string chunkNm;
    boost::uint64_t chunkSz;
    boost::uint8_t *chunk = 0;
    detail::readChunk(*dp_istrm, chunkNm, chunkSz, chunk);
    if (chunkNm == "FEND") {
      break;
    } else if (chunkNm == "POPC") {
      // popcounts
    } else if (chunkNm == "AREN") {
      dp_impl->fpData = chunk;
      chunk = NULL;
    } else if (chunkNm == "FPID") {
    }
    delete[] chunk;
  }
  df_init = true;
};

void FPBReader::destroy() {
  if (dp_impl) delete[] dp_impl->fpData;
  delete dp_impl;
};

ExplicitBitVect *FPBReader::getFP(unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(dp_impl, "no impl");
  URANGE_CHECK(idx, dp_impl->len);

  // STUB
  ExplicitBitVect *res = new ExplicitBitVect(dp_impl->nBits);
  for (unsigned int i = 0; i < 17; ++i) res->setBit(i);
  return res;
};

std::string FPBReader::getId(unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(dp_impl, "no impl");
  URANGE_CHECK(idx, dp_impl->len);
  // STUB
  std::string res = "ZINC00902219";
  return res;
};
unsigned int FPBReader::length() const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(dp_impl, "no impl");
  return dp_impl->len;
};
}  // end of RDKit namespace
