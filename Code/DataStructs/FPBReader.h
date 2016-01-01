//
// Copyright (c) 2015 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_FPBREADER_H_DEC2015
#define RD_FPBREADER_H_DEC2015

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <RDGeneral/BadFileException.h>
#include <boost/cstdint.hpp>

namespace RDKit {
namespace detail {
struct FPBReader_impl;
}
class FPBReader {
 public:
  FPBReader()
      : dp_istrm(NULL), dp_impl(NULL), df_owner(false), df_init(false){};
  FPBReader(const char *fname) { _initFromFilename(fname); };
  FPBReader(const std::string &fname) { _initFromFilename(fname.c_str()); };
  FPBReader(std::istream *inStream, bool takeOwnership = true)
      : dp_istrm(inStream), df_owner(takeOwnership), df_init(false){};
  ~FPBReader() {
    destroy();
    if (df_owner) delete dp_istrm;
    dp_istrm = NULL;
    df_init = false;
  };

  void init();
  // the caller is responsible for deleting the result
  ExplicitBitVect *getFP(unsigned int idx) const;

  std::string getId(unsigned int idx) const;
  // the caller is responsible for deleting the first element of the pair
  std::pair<ExplicitBitVect *, std::string> operator[](unsigned int idx) const {
    return std::make_pair(getFP(idx), getId(idx));
  };
  // returns the beginning and end index of fingerprints having on bit counts
  // within the range (including end points)
  std::pair<unsigned int, unsigned int> getFPIdsInCountRange(
      unsigned int minCount, unsigned int maxCount);

  unsigned int length() const;

  double getTanimoto(unsigned int idx, const boost::uint8_t *bv) const;

 private:
  std::istream *dp_istrm;
  detail::FPBReader_impl *dp_impl;  // implementation details
  bool df_owner;
  bool df_init;

  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  FPBReader(const FPBReader &);
  FPBReader &operator=(const FPBReader &);
  void destroy();
  void _initFromFilename(const char *fname) {
    std::istream *tmpStream = static_cast<std::istream *>(
        new std::ifstream(fname, std::ios_base::binary));
    if (!tmpStream || (!(*tmpStream)) || (tmpStream->bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fname;
      throw BadFileException(errout.str());
    }
    dp_istrm = tmpStream;
    df_owner = true;
    df_init = false;
  }
};
}
#endif
