//
// Copyright (c) 2016 Greg Landrum
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
#include <RDGeneral/Ranking.h>
#include "FPBReader.h"
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

namespace RDKit {

namespace detail {
const unsigned int magicSize = 8;
const std::string FPB_MAGIC("FPB1\r\n\0\0", 8);
const unsigned int tagNameSize = 4;

struct FPBReader_impl {
  unsigned int len;
  unsigned int nBits;
  boost::uint32_t numBytesStoredPerFingerprint;
  std::vector<boost::uint32_t> popCountOffsets;
  const boost::uint8_t *dp_fpData;  // do not free this
  boost::scoped_array<boost::uint8_t> dp_arenaChunk;
  boost::uint32_t num4ByteElements, num8ByteElements;  // for finding ids
  const boost::uint8_t *dp_idOffsets;                  // do not free this
  boost::scoped_array<boost::uint8_t> dp_idChunk;
  bool df_lazy;  // read the fp data lazily. In this case we use fpDataOffset
                 // and seek instead of using dp_fpData
  std::streampos fpDataOffset;   // file offset from tellg
  std::streampos idDataOffset;   // file offset from tellg
  std::streampos idChunkOffset;  // file offset from tellg
  std::istream *istrm;  // we don't own this, it's just used for the lazy reader
};

// the caller is responsible for calling delete[] on `data`
void readChunkDetails(std::istream &istrm, std::string &nm,
                      boost::uint64_t &sz) {
  streamRead(istrm, sz);
  char tag[tagNameSize + 1];
  tag[tagNameSize] = 0;
  istrm.read(tag, tagNameSize);
  nm = tag;
}
void readChunkData(std::istream &istrm, boost::uint64_t &sz,
                   boost::uint8_t *&data) {
  if (sz) {
    data = new boost::uint8_t[sz];
    istrm.read(reinterpret_cast<char *>(data), sz);
  } else {
    data = nullptr;
  }
  // std::cerr << "  CHUNKSZ: " << sz << " name: " << nm << std::endl;
}

void extractPopCounts(FPBReader_impl *dp_impl, boost::uint64_t sz,
                      const boost::uint8_t *chunk) {
  PRECONDITION(dp_impl, "bad pointer");
  /* this section of the FPB format is under-documented in Andrew's code,
   * fortunately it looks pretty simple
   */
  if (sz % 4) {
    throw ValueErrorException("POPC chunk size must be a multiple of 4 bytes");
  }
  unsigned int nEntries = sz / 4;
  if (nEntries < 9) {
    throw ValueErrorException("POPC must contain at least 9 offsets");
  }

  dp_impl->popCountOffsets.reserve(nEntries);
  for (unsigned int i = 0; i < nEntries; ++i) {
    dp_impl->popCountOffsets.push_back(
        *reinterpret_cast<const boost::uint32_t *>(chunk));
    chunk += 4;
  }
};

//-----------------------------------------------------
//  Arena processing

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
void extractArenaDetails(FPBReader_impl *dp_impl, boost::uint64_t sz) {
  PRECONDITION(dp_impl, "bad pointer");
  PRECONDITION(dp_impl->df_lazy, "should only be used in lazy mode");

  boost::uint32_t numBytesPerFingerprint;
  streamRead(*dp_impl->istrm, numBytesPerFingerprint);
  dp_impl->nBits = numBytesPerFingerprint * 8;

  boost::uint32_t numBytesStoredPerFingerprint;
  streamRead(*dp_impl->istrm, numBytesStoredPerFingerprint);
  dp_impl->numBytesStoredPerFingerprint = numBytesStoredPerFingerprint;
  boost::uint8_t spacer;
  streamRead(*dp_impl->istrm, spacer);
  dp_impl->len = (sz - 9 - spacer) / numBytesStoredPerFingerprint;

  // streamRead(*dp_impl->istrm, spacer);
  // now move forward the length of the spacer
  if (spacer) {
    dp_impl->istrm->seekg(static_cast<std::streamoff>(spacer),
                          std::ios_base::cur);
  }
  dp_impl->fpDataOffset = dp_impl->istrm->tellg();
  dp_impl->istrm->seekg(
      static_cast<std::streamoff>(numBytesStoredPerFingerprint * dp_impl->len),
      std::ios_base::cur);
}
void extractArena(FPBReader_impl *dp_impl, boost::uint64_t sz,
                  const boost::uint8_t *chunk) {
  PRECONDITION(dp_impl, "bad pointer");

  boost::uint32_t numBytesPerFingerprint =
      *reinterpret_cast<const boost::uint32_t *>(chunk);
  dp_impl->nBits = numBytesPerFingerprint * 8;

  chunk += sizeof(boost::uint32_t);
  dp_impl->numBytesStoredPerFingerprint =
      *reinterpret_cast<const boost::uint32_t *>(chunk);
  chunk += sizeof(boost::uint32_t);
  boost::uint8_t spacer = *reinterpret_cast<const boost::uint8_t *>(chunk);
  chunk += 1;
  // now move forward the length of the spacer
  chunk += spacer;

  dp_impl->dp_fpData = chunk;
  dp_impl->len = (sz - 9 - spacer) / dp_impl->numBytesStoredPerFingerprint;
};

// if dp_impl->df_lazy is true, we'll use the memory in fpData (should be large
// enough to hold the result!), otherwise
// we update it to a pointer to the memory dp_impl owns.
void extractBytes(const FPBReader_impl *dp_impl, unsigned int which,
                  boost::uint8_t *&fpData, unsigned int nToRead = 1) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION((dp_impl->df_lazy || dp_impl->dp_fpData), "bad fpdata pointer");
  PRECONDITION(!dp_impl->df_lazy || dp_impl->istrm, "no stream in lazy mode");
  PRECONDITION(!dp_impl->df_lazy || fpData, "no fpData in lazy mode");
  PRECONDITION(nToRead > 0, "bad nToRead");

  if (which + nToRead > dp_impl->len) {
    throw ValueErrorException("bad index");
  }
  boost::uint64_t offset = which * dp_impl->numBytesStoredPerFingerprint;
  if (!dp_impl->df_lazy) {
    fpData = const_cast<boost::uint8_t *>(dp_impl->dp_fpData) + offset;
  } else {
    dp_impl->istrm->seekg(dp_impl->fpDataOffset +
                          static_cast<std::streampos>(offset));
    dp_impl->istrm->read(reinterpret_cast<char *>(fpData),
                         nToRead * dp_impl->numBytesStoredPerFingerprint);
  }
};

// the caller is responsible for delete[]'ing this
boost::uint8_t *copyBytes(const FPBReader_impl *dp_impl, unsigned int which) {
  PRECONDITION(dp_impl, "bad reader pointer");
  boost::uint8_t *res;
  res = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  if (!dp_impl->df_lazy) {
    boost::uint8_t *fpData = nullptr;
    extractBytes(dp_impl, which, fpData);
    memcpy(static_cast<void *>(res), fpData,
           dp_impl->numBytesStoredPerFingerprint);
  } else {
    extractBytes(dp_impl, which, res);
  }
  return res;
};

// caller is responsible for delete'ing the result
RDKIT_DATASTRUCTS_EXPORT boost::dynamic_bitset<> *bytesToBitset(
    const boost::uint8_t *fpData, boost::uint32_t nBits) {
  unsigned int nBytes = nBits / 8;
  if (!(nBytes % sizeof(boost::dynamic_bitset<>::block_type))) {
    // I believe this could be faster (needs to be verified of course)
    unsigned int nBlocks = nBytes / sizeof(boost::dynamic_bitset<>::block_type);
    const auto *fpBlocks =
        reinterpret_cast<const boost::dynamic_bitset<>::block_type *>(fpData);
    return new boost::dynamic_bitset<>(fpBlocks, fpBlocks + nBlocks);
  } else {
    return reinterpret_cast<boost::dynamic_bitset<> *>(
        new boost::dynamic_bitset<boost::uint8_t>(fpData, fpData + nBytes));
  }
}

// caller is responsible for delete []'ing the result
RDKIT_DATASTRUCTS_EXPORT boost::uint8_t *bitsetToBytes(
    const boost::dynamic_bitset<> &bitset) {
  unsigned int nBits = bitset.size();
  unsigned int nBytes = nBits / 8;

  auto *res = new boost::uint8_t[nBytes];
  boost::to_block_range(
      bitset, reinterpret_cast<boost::dynamic_bitset<>::block_type *>(res));
  return res;
}

// the caller is responsible for delete'ing this
ExplicitBitVect *extractFP(const FPBReader_impl *dp_impl, unsigned int which) {
  PRECONDITION(dp_impl, "bad reader pointer");
  boost::uint8_t *fpData;
  if (dp_impl->df_lazy) {
    fpData = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  }
  extractBytes(dp_impl, which, fpData);
  boost::dynamic_bitset<> *resDBS = bytesToBitset(fpData, dp_impl->nBits);
  if (dp_impl->df_lazy) {
    delete[] fpData;
  }
  return new ExplicitBitVect(resDBS);
};

double tanimoto(const FPBReader_impl *dp_impl, unsigned int which,
                const ::boost::uint8_t *bv) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(bv, "bad bv pointer");
  if (which >= dp_impl->len) {
    throw ValueErrorException("bad index");
  }
  boost::uint8_t *fpData;
  if (dp_impl->df_lazy) {
    fpData = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  }
  extractBytes(dp_impl, which, fpData);
  double res =
      CalcBitmapTanimoto(fpData, bv, dp_impl->numBytesStoredPerFingerprint);
  if (dp_impl->df_lazy) {
    delete[] fpData;
  }
  return res;
};

double tversky(const FPBReader_impl *dp_impl, unsigned int which,
               const ::boost::uint8_t *bv, double ca, double cb) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(bv, "bad bv pointer");
  if (which >= dp_impl->len) {
    throw ValueErrorException("bad index");
  }
  boost::uint8_t *fpData;
  if (dp_impl->df_lazy) {
    fpData = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  }
  extractBytes(dp_impl, which, fpData);
  double res = CalcBitmapTversky(fpData, bv,
                                 dp_impl->numBytesStoredPerFingerprint, ca, cb);
  if (dp_impl->df_lazy) {
    delete[] fpData;
  }
  return res;
};

//-----------------------------------------------------
//  Id procesing
/* Documentation from Andrew's code on the structure of the arena:

The actual layout inside of the chunk is:
 <num_4byte_elements: 4 bytes> -- the number of 4 byte offsets.
 <num_8byte_elements: 4 bytes> -- the number of 8 byte offsets
 Note: the number of indices is num_4byte_elements + num_8byte_elements + 1
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
void extractIdsDetails(FPBReader_impl *dp_impl, boost::uint64_t sz) {
  PRECONDITION(dp_impl, "bad pointer");
  std::streampos start = dp_impl->istrm->tellg();
  dp_impl->idChunkOffset = start;
  streamRead(*dp_impl->istrm, dp_impl->num4ByteElements);
  streamRead(*dp_impl->istrm, dp_impl->num8ByteElements);

  dp_impl->idDataOffset = static_cast<boost::uint64_t>(start) + sz -
                          (dp_impl->num4ByteElements + 1) * 4 -
                          dp_impl->num8ByteElements * 8;
  dp_impl->istrm->seekg(start + static_cast<std::streampos>(sz),
                        std::ios_base::beg);
};

void extractIds(FPBReader_impl *dp_impl, boost::uint64_t sz,
                const boost::uint8_t *chunk) {
  PRECONDITION(dp_impl, "bad pointer");
  dp_impl->num4ByteElements = *reinterpret_cast<const boost::uint32_t *>(chunk);
  chunk += sizeof(boost::uint32_t);
  dp_impl->num8ByteElements = *reinterpret_cast<const boost::uint32_t *>(chunk);
  chunk += sizeof(boost::uint32_t);
  dp_impl->dp_idOffsets = dp_impl->dp_idChunk.get() + sz -
                          (dp_impl->num4ByteElements + 1) * 4 -
                          dp_impl->num8ByteElements * 8;
};

std::string extractId(const FPBReader_impl *dp_impl, unsigned int which) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION((dp_impl->df_lazy || dp_impl->dp_idOffsets),
               "bad idOffsets pointer");
  PRECONDITION(!dp_impl->df_lazy || dp_impl->istrm, "no stream in lazy mode");

  if (which >= dp_impl->num4ByteElements + dp_impl->num8ByteElements) {
    throw ValueErrorException("bad index");
  }
  std::string res;

  boost::uint64_t offset = 0, len = 0;
  if (which < dp_impl->num4ByteElements) {
    if (!dp_impl->df_lazy) {
      offset = *reinterpret_cast<const boost::uint32_t *>(
          dp_impl->dp_idOffsets + which * 4);
      len = *reinterpret_cast<const boost::uint32_t *>(dp_impl->dp_idOffsets +
                                                       (which + 1) * 4);
    } else {
      dp_impl->istrm->seekg(dp_impl->idDataOffset +
                            static_cast<std::streampos>(which * 4));
      dp_impl->istrm->read(reinterpret_cast<char *>(&offset), 4);
      dp_impl->istrm->read(reinterpret_cast<char *>(&len), 4);
    }
  } else if (which == dp_impl->num4ByteElements) {
    // FIX: this code path is not yet tested
    if (!dp_impl->df_lazy) {
      offset = *reinterpret_cast<const boost::uint32_t *>(
          dp_impl->dp_idOffsets + which * 4);
      len = *reinterpret_cast<const boost::uint64_t *>(dp_impl->dp_idOffsets +
                                                       (which + 1) * 4);
    } else {
      dp_impl->istrm->seekg(dp_impl->idDataOffset +
                            static_cast<std::streampos>(which * 4));
      dp_impl->istrm->read(reinterpret_cast<char *>(&offset), 4);
      dp_impl->istrm->read(reinterpret_cast<char *>(&len), 8);
    }
  } else {
    // FIX: this code path is not yet tested
    if (!dp_impl->df_lazy) {
      offset = *reinterpret_cast<const boost::uint64_t *>(
          dp_impl->dp_idOffsets + dp_impl->num4ByteElements * 4 + which * 8);
      len = *reinterpret_cast<const boost::uint64_t *>(
          dp_impl->dp_idOffsets + dp_impl->num4ByteElements * 4 +
          (which + 1) * 8);
    } else {
      dp_impl->istrm->seekg(dp_impl->idDataOffset +
                            static_cast<std::streampos>(
                                dp_impl->num4ByteElements * 4 + which * 8));
      dp_impl->istrm->read(reinterpret_cast<char *>(&offset), 8);
      dp_impl->istrm->read(reinterpret_cast<char *>(&len), 8);
    }
  }
  len -= offset;

  if (!dp_impl->df_lazy) {
    res = std::string(
        reinterpret_cast<const char *>(dp_impl->dp_idChunk.get() + offset),
        len);
  } else {
    boost::shared_array<char> buff(new char[len + 1]);
    buff[len] = 0;
    dp_impl->istrm->seekg(dp_impl->idChunkOffset +
                          static_cast<std::streampos>(offset));
    dp_impl->istrm->read(reinterpret_cast<char *>(buff.get()), len);
    res = std::string(reinterpret_cast<const char *>(buff.get()));
  }
  return res;
};

void tanimotoNeighbors(const FPBReader_impl *dp_impl, const boost::uint8_t *bv,
                       double threshold,
                       std::vector<std::pair<double, unsigned int>> &res,
                       bool usePopcountScreen, unsigned int readCache = 1000) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(bv, "bad bv");
  RANGE_CHECK(-1e-6, threshold, 1.0 + 1e-6);
  PRECONDITION(readCache > 0, "bad cache size");
  res.clear();
  boost::uint64_t probeCount =
      CalcBitmapPopcount(bv, dp_impl->numBytesStoredPerFingerprint);

  boost::uint64_t startScan = 0, endScan = dp_impl->len;
  if (usePopcountScreen &&
      dp_impl->popCountOffsets.size() == dp_impl->nBits + 2) {
    // figure out the bounds based on equation 24 from:
    // 1. Swamidass, S. J. & Baldi, P. Bounds and Algorithms for Fast Exact
    // Searches of Chemical Fingerprints in Linear and Sublinear Time. J. Chem.
    // Inf. Model. 47, 302–317 (2007).
    // http://pubs.acs.org/doi/abs/10.1021/ci600358f
    auto minDbCount =
        static_cast<boost::uint32_t>(floor(threshold * probeCount));
    boost::uint32_t maxDbCount =
        (threshold > 1e-6)
            ? static_cast<boost::uint32_t>(ceil(probeCount / threshold))
            : dp_impl->numBytesStoredPerFingerprint;
    // std::cerr << "probeCount: " << probeCount << " bounds: " << minDbCount
    //           << "-" << maxDbCount << std::endl;
    startScan = dp_impl->popCountOffsets[minDbCount];
    endScan = dp_impl->popCountOffsets[maxDbCount + 1];
    // std::cerr << " scan: " << startScan << "-" << endScan << std::endl;
  }
  boost::uint8_t *dbv;
  if (dp_impl->df_lazy) {
    dbv = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint * readCache];
  }
  for (boost::uint64_t i = startScan; i < endScan; i += readCache) {
    unsigned int toRead = readCache;
    if (i + toRead >= endScan) {
      toRead = endScan - i;
    }
    extractBytes(dp_impl, i, dbv, toRead);
    for (unsigned int j = 0; j < toRead; ++j) {
      double tani =
          CalcBitmapTanimoto(dbv + j * dp_impl->numBytesStoredPerFingerprint,
                             bv, dp_impl->numBytesStoredPerFingerprint);
      if (tani >= threshold) {
        res.emplace_back(tani, i + j);
      }
    }
  }
  if (dp_impl->df_lazy) {
    delete[] dbv;
  }
}

void tverskyNeighbors(const FPBReader_impl *dp_impl, const boost::uint8_t *bv,
                      double ca, double cb, double threshold,
                      std::vector<std::pair<double, unsigned int>> &res,
                      bool usePopcountScreen) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(bv, "bad bv");
  RANGE_CHECK(-1e-6, threshold, 1.0 + 1e-6);
  res.clear();
  boost::uint64_t probeCount =
      CalcBitmapPopcount(bv, dp_impl->numBytesStoredPerFingerprint);

  boost::uint64_t startScan = 0, endScan = dp_impl->len;
  if (usePopcountScreen &&
      dp_impl->popCountOffsets.size() == dp_impl->nBits + 2) {
    // figure out the bounds based on equation 25 from:
    // 1. Swamidass, S. J. & Baldi, P. Bounds and Algorithms for Fast Exact
    // Searches of Chemical Fingerprints in Linear and Sublinear Time. J. Chem.
    // Inf. Model. 47, 302–317 (2007).
    // http://pubs.acs.org/doi/abs/10.1021/ci600358f
    auto minDbCount = static_cast<boost::uint32_t>(floor(
        (threshold * probeCount * ca) / (1. - threshold + threshold * ca)));
    boost::uint32_t maxDbCount =
        ((threshold * cb) > 1e-6)
            ? static_cast<boost::uint32_t>(
                  ceil(probeCount * (1 - threshold + threshold * cb) /
                       (threshold * cb)))
            : dp_impl->numBytesStoredPerFingerprint;
    // std::cerr << "probeCount: " << probeCount << " bounds: " << minDbCount
    //          << "-" << maxDbCount << std::endl;
    startScan = dp_impl->popCountOffsets[minDbCount];
    endScan = dp_impl->popCountOffsets[maxDbCount + 1];
    // std::cerr << " scan: " << startScan << "-" << endScan << std::endl;
  }

  boost::uint8_t *dbv;
  if (dp_impl->df_lazy) {
    dbv = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  }
  for (boost::uint64_t i = startScan; i < endScan; ++i) {
    extractBytes(dp_impl, i, dbv);
    double sim = CalcBitmapTversky(
        dbv, bv, dp_impl->numBytesStoredPerFingerprint, ca, cb);
    // std::cerr << "  i:" << i << " " << tani << " ? " << threshold <<
    // std::endl;
    if (sim >= threshold) {
      res.emplace_back(sim, i);
    }
  }
  if (dp_impl->df_lazy) {
    delete[] dbv;
  }
}

void containingNeighbors(const FPBReader_impl *dp_impl,
                         const boost::uint8_t *bv,
                         std::vector<unsigned int> &res) {
  PRECONDITION(dp_impl, "bad reader pointer");
  PRECONDITION(bv, "bad bv");
  res.clear();
  boost::uint64_t probeCount =
      CalcBitmapPopcount(bv, dp_impl->numBytesStoredPerFingerprint);

  boost::uint64_t startScan = 0, endScan = dp_impl->len;
  if (dp_impl->popCountOffsets.size() == dp_impl->nBits + 2) {
    startScan = dp_impl->popCountOffsets[probeCount];
    // std::cerr << " scan: " << startScan << "-" << endScan << std::endl;
  }
  boost::uint8_t *dbv;
  if (dp_impl->df_lazy) {
    dbv = new boost::uint8_t[dp_impl->numBytesStoredPerFingerprint];
  }
  for (boost::uint64_t i = startScan; i < endScan; ++i) {
    extractBytes(dp_impl, i, dbv);
    if (CalcBitmapAllProbeBitsMatch(bv, dbv,
                                    dp_impl->numBytesStoredPerFingerprint)) {
      res.push_back(i);
    }
  }
  if (dp_impl->df_lazy) {
    delete[] dbv;
  }
}

}  // namespace detail

void FPBReader::init() {
  PRECONDITION(dp_istrm, "no stream");
  if (df_init) {
    return;
  }

  dp_impl = new detail::FPBReader_impl;
  dp_impl->istrm = dp_istrm;
  dp_impl->df_lazy = df_lazyRead;

  char magic[detail::magicSize];
  dp_istrm->read(magic, detail::magicSize);
  if (detail::FPB_MAGIC != std::string(magic, detail::magicSize)) {
    throw BadFileException("Invalid FPB magic");
  }
  while (1) {
    if (dp_istrm->eof()) {
      throw BadFileException("EOF hit before FEND record");
    }
    std::string chunkNm;
    boost::uint64_t chunkSz;
    boost::uint8_t *chunk = nullptr;
    detail::readChunkDetails(*dp_istrm, chunkNm, chunkSz);
    // std::cerr << " Chunk: " << chunkNm << " " << chunkSz << std::endl;
    if (!df_lazyRead || (chunkNm != "AREN" && chunkNm != "FPID")) {
      detail::readChunkData(*dp_istrm, chunkSz, chunk);
      if (chunkNm == "FEND") {
        break;
      } else if (chunkNm == "POPC") {
        detail::extractPopCounts(dp_impl, chunkSz, chunk);
      } else if (chunkNm == "AREN") {
        dp_impl->dp_arenaChunk.reset(chunk);
        detail::extractArena(dp_impl, chunkSz, chunk);
        chunk = nullptr;
      } else if (chunkNm == "FPID") {
        dp_impl->dp_idChunk.reset(chunk);
        detail::extractIds(dp_impl, chunkSz, chunk);
        chunk = nullptr;
      } else if (chunkNm == "META") {
        // currently ignored
      } else if (chunkNm == "HASH") {
        // currently ignored
      } else {
        BOOST_LOG(rdWarningLog)
            << "Unknown chunk: " << chunkNm << " ignored." << std::endl;
      }
      delete[] chunk;
    } else {
      // we are reading the AREN or FPID chunk in lazy mode, just get our
      // position in
      // the file.
      if (chunkNm == "AREN") {
        detail::extractArenaDetails(dp_impl, chunkSz);
      } else if (chunkNm == "FPID") {
        detail::extractIdsDetails(dp_impl, chunkSz);
      }
    }
  }
  if ((!df_lazyRead && !dp_impl->dp_arenaChunk) ||
      (df_lazyRead && !dp_impl->fpDataOffset)) {
    throw BadFileException("No AREN record found");
  }
  if ((!df_lazyRead && !dp_impl->dp_idChunk) ||
      (df_lazyRead && !dp_impl->idDataOffset)) {
    throw BadFileException("No FPID record found");
  }

  df_init = true;
};

void FPBReader::destroy() {
  if (dp_impl) {
    dp_impl->dp_arenaChunk.reset();
    dp_impl->dp_idChunk.reset();

    dp_impl->dp_fpData = nullptr;
    dp_impl->dp_idOffsets = nullptr;
  }
  delete dp_impl;
  dp_impl = nullptr;
};

boost::shared_ptr<ExplicitBitVect> FPBReader::getFP(unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");

  return boost::shared_ptr<ExplicitBitVect>(detail::extractFP(dp_impl, idx));
};
boost::shared_array<boost::uint8_t> FPBReader::getBytes(
    unsigned int idx) const {
  PRECONDITION(df_init, "not initialized");

  return boost::shared_array<boost::uint8_t>(detail::copyBytes(dp_impl, idx));
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
unsigned int FPBReader::nBits() const {
  PRECONDITION(df_init, "not initialized");
  PRECONDITION(dp_impl, "no impl");
  return dp_impl->nBits;
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
double FPBReader::getTanimoto(unsigned int idx,
                              const ExplicitBitVect &ebv) const {
  const boost::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  double res = getTanimoto(idx, bv);
  delete[] bv;
  return res;
}

std::vector<std::pair<double, unsigned int>> FPBReader::getTanimotoNeighbors(
    const boost::uint8_t *bv, double threshold, bool usePopcountScreen) const {
  PRECONDITION(df_init, "not initialized");
  std::vector<std::pair<double, unsigned int>> res;
  detail::tanimotoNeighbors(dp_impl, bv, threshold, res, usePopcountScreen);
  std::sort(res.begin(), res.end(), Rankers::pairGreater);
  return res;
}

std::vector<std::pair<double, unsigned int>> FPBReader::getTanimotoNeighbors(
    const ExplicitBitVect &ebv, double threshold,
    bool usePopcountScreen) const {
  const boost::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  std::vector<std::pair<double, unsigned int>> res =
      getTanimotoNeighbors(bv, threshold, usePopcountScreen);
  delete[] bv;
  return res;
}

double FPBReader::getTversky(unsigned int idx, const boost::uint8_t *bv,
                             double ca, double cb) const {
  PRECONDITION(df_init, "not initialized");
  return detail::tversky(dp_impl, idx, bv, ca, cb);
}
double FPBReader::getTversky(unsigned int idx, const ExplicitBitVect &ebv,
                             double ca, double cb) const {
  const boost::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  double res = getTversky(idx, bv, ca, cb);
  delete[] bv;
  return res;
}

std::vector<std::pair<double, unsigned int>> FPBReader::getTverskyNeighbors(
    const boost::uint8_t *bv, double ca, double cb, double threshold,
    bool usePopcountScreen) const {
  PRECONDITION(df_init, "not initialized");
  std::vector<std::pair<double, unsigned int>> res;
  detail::tverskyNeighbors(dp_impl, bv, ca, cb, threshold, res,
                           usePopcountScreen);
  std::sort(res.begin(), res.end(), Rankers::pairGreater);
  return res;
}

std::vector<std::pair<double, unsigned int>> FPBReader::getTverskyNeighbors(
    const ExplicitBitVect &ebv, double ca, double cb, double threshold,
    bool usePopcountScreen) const {
  const boost::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  std::vector<std::pair<double, unsigned int>> res =
      getTverskyNeighbors(bv, ca, cb, threshold, usePopcountScreen);
  delete[] bv;
  return res;
}

std::vector<unsigned int> FPBReader::getContainingNeighbors(
    const boost::uint8_t *bv) const {
  PRECONDITION(df_init, "not initialized");
  std::vector<unsigned int> res;
  detail::containingNeighbors(dp_impl, bv, res);
  std::sort(res.begin(), res.end());

  return res;
}

std::vector<unsigned int> FPBReader::getContainingNeighbors(
    const ExplicitBitVect &ebv) const {
  const boost::uint8_t *bv = detail::bitsetToBytes(*(ebv.dp_bits));
  std::vector<unsigned int> res = getContainingNeighbors(bv);
  delete[] bv;
  return res;
}

}  // namespace RDKit
