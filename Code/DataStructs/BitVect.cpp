// $Id$
//
// Copyright (c) 2003-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BitVect.h"
#include <sstream>
#include <limits>
#include <RDGeneral/StreamOps.h>
#include "base64.h"
#ifdef WIN32
#include <ios>
#endif
#include <cstdint>

BitVect::~BitVect(){};  // must always implement virtual destructors

void BitVect::initFromText(const char *data, const unsigned int dataLen,
                           bool isBase64, bool allowOldFormat) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  if (isBase64) {
    unsigned int actualLen;
    char *decoded;
    decoded = Base64Decode((const char *)data, &actualLen);
    ss.write(decoded, actualLen);
    delete[] decoded;
  } else {
    ss.write(data, dataLen);
  }

  std::int32_t format = 0;
  std::int32_t version = 0;
  std::uint32_t nOn = 0;
  std::int32_t size;

  // earlier versions of the code did not have the version number encoded, so
  //  we'll use that to distinguish version 0
  RDKit::streamRead(ss, size);
  if (size < 0) {
    version = -1 * size;
    switch (version) {
      case 16:
        format = 1;
        break;
      case 32:
        format = 2;
        break;
      default:
        throw ValueErrorException("bad version in BitVect pickle");
        break;
    }
    RDKit::streamRead(ss, size);
  } else if (!allowOldFormat) {
    throw ValueErrorException("invalid BitVect pickle");
  }

  RDKit::streamRead(ss, nOn);
  _initForSize(static_cast<int>(size));

  // if the either have older version or or version 16 with ints for on bits
  if ((format == 0) ||
      ((format == 1) && (size >= std::numeric_limits<unsigned short>::max()))) {
    std::uint32_t tmp;
    for (unsigned int i = 0; i < nOn; i++) {
      RDKit::streamRead(ss, tmp);
      setBit(tmp);
    }
  } else if (format == 1) {  // version 16 and on bits stored as short ints
    std::uint16_t tmp;
    for (unsigned int i = 0; i < nOn; i++) {
      RDKit::streamRead(ss, tmp);
      setBit(tmp);
    }
  } else if (format == 2) {  // run length encoded format
    std::uint32_t curr = 0;
    for (unsigned int i = 0; i < nOn; i++) {
      curr += RDKit::readPackedIntFromStream(ss);
      setBit(curr);
      curr++;
    }
  }
}
