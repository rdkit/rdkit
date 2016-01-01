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

#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/StreamOps.h>
#include "FPBReader.h"

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
  // std::cerr << "  CHUNKSZ: " << sz << " name: " << nm << std::endl;
}

struct FPBReader_impl {
  unsigned int len;
  unsigned int nBits;
  boost::uint32_t numBytesStoredPerFingerprint;
  std::vector<boost::uint32_t> popCountOffsets;
  const boost::uint8_t *dp_fpData;      // do not free this
  const boost::uint8_t *dp_arenaChunk;  // this is what should be freed
  boost::uint32_t num4ByteElements, num8ByteElements;  // for finding ids
  const boost::uint8_t *dp_idOffsets;                  // do not free this
  const boost::uint8_t *dp_idChunk;                    // free this
};

void extractPopCounts(FPBReader_impl *dp_impl, boost::uint64_t sz,
                      const boost::uint8_t *chunk) {
  PRECONDITION(dp_impl, "bad pointer");
  /* this section of the FPB format is under-documented in Andrew's code,
   * fortunately it looks pretty simple
  */
  if (sz % 4)
    throw ValueErrorException("POPC chunk size must be a multiple of 4 bytes");
  unsigned int nEntries = sz / 4;
  if (nEntries < 9)
    throw ValueErrorException("POPC must contain at least 9 offsets");

  dp_impl->popCountOffsets.reserve(nEntries);
  for (unsigned int i = 0; i < nEntries; ++i) {
    dp_impl->popCountOffsets.push_back(*(boost::uint32_t *)chunk);
    chunk += 4;
  }
  // FIX: Finish this;
};
void extractArena(FPBReader_impl *dp_impl, boost::uint64_t sz,
                  const boost::uint8_t *chunk) {
  PRECONDITION(dp_impl, "bad pointer");
  /* Documentation from Andrew's code on the structure of the arena:
The 'AREN'a starts with a header:
  <num_bytes: 4 bytes>  -- the number of bytes in a fingerprint
  <storage_size: 4 bytes>  -- number of bytes in fingerprint + extra bytes
  <spacer_size: 1 byte>   -- the number of spacer bytes used so the fingerprint
         chunk starts on an aligned file position.
  <spacer : $spacer_size> NUL bytes> -- up to 255 NUL bytes, used for alignment.
The fingerprints are N fingerprint fields, ordered sequentially.
  <fp0: $storage_size bytes> -- the first fingerprint
  <fp1: $storage_size bytes> -- the second fingerprint
     ...
The last fingerprint ends at the last byte of the arena chunk.

Each fingerprint contains:
  <fingerprint: $num_bytes bytes> -- the actual fingerprint data
  <extra: $storage_size-$num_bytes bytes> -- the 'extra' NULL padding bytes
      used so storage_size is a multiple of the alignment.

To get the number of fingerprints in the arena:
   (len(arena content) - 4 - 4 - 1 - $spacer_size) // $storage_size
   */
  boost::uint32_t numBytesPerFingerprint = *((boost::uint32_t *)chunk);
  dp_impl->nBits = numBytesPerFingerprint * 8;

  chunk += sizeof(boost::uint32_t);
  dp_impl->numBytesStoredPerFingerprint = *((boost::uint32_t *)chunk);
  chunk += sizeof(boost::uint32_t);
  boost::uint8_t spacer = *((boost::uint8_t *)chunk);
  chunk += 1;
  // now move forward the length of the spacer
  chunk += spacer;

  dp_impl->dp_fpData = chunk;
  dp_impl->len = (sz - 9 - spacer) / dp_impl->numBytesStoredPerFingerprint;
};

// the caller is responsible for delete'ing this
ExplicitBitVect *extractFP(const FPBReader_impl *dp_impl, unsigned int which) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(dp_impl->dp_fpData, "bad fpdata pointer");

  if (which >= dp_impl->len) {
    throw ValueErrorException("bad index");
  }
  const boost::uint8_t *fpData =
      dp_impl->dp_fpData + which * dp_impl->numBytesStoredPerFingerprint;
  boost::dynamic_bitset<boost::uint8_t> *fpbs =
      new boost::dynamic_bitset<boost::uint8_t>(fpData,
                                                fpData + dp_impl->nBits / 8);
  return new ExplicitBitVect((boost::dynamic_bitset<> *)fpbs);
};
// the caller is responsible for delete[]'ing this
boost::uint8_t *extractBytes(const FPBReader_impl *dp_impl,
                             unsigned int which) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(dp_impl->dp_fpData, "bad fpdata pointer");

  if (which >= dp_impl->len) {
    throw ValueErrorException("bad index");
  }
  const boost::uint8_t *fpData =
      dp_impl->dp_fpData + which * dp_impl->numBytesStoredPerFingerprint;
  boost::uint8_t *res =
      new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  memcpy((void *)res, (void *)fpData, dp_impl->numBytesStoredPerFingerprint);
  return res;
};

double tanimoto(const FPBReader_impl *dp_impl, unsigned int which,
                const ::boost::uint8_t *bv) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(bv, "bad bv pointer");
  if (which >= dp_impl->len) {
    throw ValueErrorException("bad index");
  }
  const boost::uint8_t *fpData =
      dp_impl->dp_fpData + which * dp_impl->numBytesStoredPerFingerprint;
  return CalcBitmapTanimoto(fpData, bv, dp_impl->numBytesStoredPerFingerprint);
};

void extractIds(FPBReader_impl *dp_impl, boost::uint64_t sz,
                const boost::uint8_t *chunk) {
  PRECONDITION(dp_impl, "bad pointer");
  /* Documentation from Andrew's code on the structure of the arena:

  The actual layout inside of the chunk is:
   <num_4byte_elements: 4 bytes> -- the number of 4 byte offsets.
   <num_8byte_elements: 4 bytes> -- the number of 8 byte offsets
   Note: the number of indicies is num_4byte_elements + num_8byte_elements + 1
   because even with no elements there will be the initial '\0\0\0\0'.

   <id 0> + NUL   -- the first string, with an added NUL terminator
   <id 1> + NUL   -- the second string, with an added NUL terminator
       ....
   <id N> + NUL   -- the last string, with an added NUL terminator

   <offset 0: 4 bytes>    -- the offset relative to the start of <text 0>.
                    (This always contains the 4 bytes "\0\0\0\0")
      ...
   <offset num_4byte_elements: 4 bytes>   -- the last offset stored in 4 bytes

        (Note: This next section exists only when <num 8 byte offsets> > 0)
   <offset num_4byte_elements+1: 8 bytes>    -- the first offset stored in 8
  bytes
      ...
   <offset num_4byte_elements+num_8byte_elements: 8 bytes>   -- the last offset
  stored in 8 bytes

  To get the identifier for record at position P >= 0:
    chunk_size = size of the chunk
    num_4byte_elements = decode bytes[0:4] as uint32
    num_8byte_elements = decode bytes[4:8] as uint32
    if P >= num_4byte_elements + num_8byte_elements:
        record does not exist
    offset_start = chunk_size - num_4byte_elements*4 - num_8byte_elements*8
    if P < num_4byte_elements:
      start, end = decode bytes[offset_start:offset_start+8] as (uint32, uint32)
    elif P == N4:
      start, end = decode bytes[offset_start:offset_start+12] as (uint32,
  uint64)
    else:
      start, end = decode bytes[offset_start:offset_start+16] as (uint64,
  uint64)
    id = bytes[start:end-1]
   */
  dp_impl->dp_idChunk = chunk;
  dp_impl->num4ByteElements = *((boost::uint32_t *)chunk);
  chunk += sizeof(boost::uint32_t);
  dp_impl->num8ByteElements = *((boost::uint32_t *)chunk);
  chunk += sizeof(boost::uint32_t);
  dp_impl->dp_idOffsets = dp_impl->dp_idChunk + sz -
                          (dp_impl->num4ByteElements + 1) * 4 -
                          dp_impl->num8ByteElements * 8;
};

std::string extractId(const FPBReader_impl *dp_impl, unsigned int which) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(dp_impl->dp_idOffsets, "bad idOffsets pointer");

  if (which >= dp_impl->num4ByteElements + dp_impl->num8ByteElements) {
    throw ValueErrorException("bad index");
  }
  std::string res;

  boost::uint64_t offset = 0, len = 0;
  if (which < dp_impl->num4ByteElements) {
    offset = *(boost::uint32_t *)(dp_impl->dp_idOffsets + which * 4);
    len = *(boost::uint32_t *)(dp_impl->dp_idOffsets + (which + 1) * 4);
  } else if (which == dp_impl->num4ByteElements) {
    // FIX: this code path is not yet tested
    offset = *(boost::uint32_t *)(dp_impl->dp_idOffsets + which * 4);
    len = *(boost::uint64_t *)(dp_impl->dp_idOffsets + (which + 1) * 4);
  } else {
    // FIX: this code path is not yet tested
    offset = *(boost::uint64_t *)(dp_impl->dp_idOffsets +
                                  dp_impl->num4ByteElements * 4 + which * 8);
    len = *(boost::uint64_t *)(dp_impl->dp_idOffsets +
                               dp_impl->num4ByteElements * 4 + (which + 1) * 8);
  }
  len -= offset;
  res = std::string((const char *)(dp_impl->dp_idChunk + offset), len);
  return res;
};
}  // end of detail namespace

void FPBReader::init() {
  PRECONDITION(dp_istrm, "no stream");
  dp_impl = new detail::FPBReader_impl;

  char magic[detail::magicSize];
  dp_istrm->read(magic, detail::magicSize);
  if (detail::FPB_MAGIC != std::string(magic, detail::magicSize)) {
    throw BadFileException("Invalid FPB magic");
  }
  while (1) {
    if (dp_istrm->eof()) throw BadFileException("EOF hit before FEND record");
    std::string chunkNm;
    boost::uint64_t chunkSz;
    boost::uint8_t *chunk = NULL;
    detail::readChunk(*dp_istrm, chunkNm, chunkSz, chunk);
    if (chunkNm == "FEND") {
      break;
    } else if (chunkNm == "POPC") {
      detail::extractPopCounts(dp_impl, chunkSz, chunk);
    } else if (chunkNm == "AREN") {
      dp_impl->dp_arenaChunk = chunk;
      detail::extractArena(dp_impl, chunkSz, chunk);
      chunk = NULL;
    } else if (chunkNm == "FPID") {
      detail::extractIds(dp_impl, chunkSz, chunk);
      chunk = NULL;
    }
    delete[] chunk;
  }
  if (!dp_impl->dp_arenaChunk) throw BadFileException("No AREN record found");
  if (!dp_impl->dp_idChunk) throw BadFileException("No FPID record found");

  df_init = true;
};

void FPBReader::destroy() {
  if (dp_impl) {
    delete[] dp_impl->dp_arenaChunk;
    dp_impl->dp_arenaChunk = NULL;

    delete[] dp_impl->dp_idChunk;
    dp_impl->dp_idChunk = NULL;

    dp_impl->dp_fpData = NULL;
    dp_impl->dp_idOffsets = NULL;
  }
  delete dp_impl;
};

ExplicitBitVect *FPBReader::getFP(unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");

  ExplicitBitVect *res = detail::extractFP(dp_impl, idx);
  return res;
};
boost::uint8_t *FPBReader::getBytes(unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");

  boost::uint8_t *res = detail::extractBytes(dp_impl, idx);
  return res;
};

std::string FPBReader::getId(unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");

  std::string res = detail::extractId(dp_impl, idx);
  return res;
};
unsigned int FPBReader::length() const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(dp_impl, "no impl");
  return dp_impl->len;
};
std::pair<unsigned int, unsigned int> FPBReader::getFPIdsInCountRange(
    unsigned int minCount, unsigned int maxCount) {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(dp_impl, "no impl");
  URANGE_CHECK(maxCount, dp_impl->nBits + 1);
  PRECONDITION(maxCount >= minCount, "max < min");
  if (dp_impl->popCountOffsets.size() == dp_impl->nBits + 2) {
    return std::make_pair(dp_impl->popCountOffsets[minCount],
                          dp_impl->popCountOffsets[maxCount + 1]);
  } else {
    // we don't have popcounts, so we have to work for it.
    // FIX: complete this
    return std::make_pair(0, 0);
  }
};
double FPBReader::getTanimoto(unsigned int idx,
                              const boost::uint8_t *bv) const {
  PRECONDITION(df_init, "not initialized");
  return detail::tanimoto(dp_impl, idx, bv);
}
}  // end of RDKit namespace
