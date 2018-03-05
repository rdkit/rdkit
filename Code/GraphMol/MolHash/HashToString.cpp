// $Id$
//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <string>
#include <sstream>
#include <algorithm>

#include "MolHash.h"

using namespace boost::archive::iterators;

namespace RDKit {
namespace MolHash {
//=============================================================================
// MolHash Module API implementation:
//=============================================================================

std::string encode(const void* bin, size_t size) {
  typedef base64_from_binary<transform_width<const char*, 6, 8> >  // retrieve 6
                                                                   // bit
                                                                   // integers
                                                                   // from a
                                                                   // sequence
                                                                   // of 8 bit
                                                                   // bytes
      base64_enc;

  std::stringstream os;
  std::copy(base64_enc((const char*)bin), base64_enc(((const char*)bin) + size),
            ostream_iterator<char>(os));
  return os.str();
}
//=============================================================================
}
}
