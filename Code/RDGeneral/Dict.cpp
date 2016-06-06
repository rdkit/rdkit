// $Id$
//
// Copyright (C) 2015 Novartis Institute of BioMedical Research
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Dict.h"
#include "Exceptions.h"
#include "tags.h"

namespace RDKit {

RDTags RDKit::KeyIntPair::tagmap;


void Dict::getVal(int tag, std::string &res) const {
  //
  //  We're going to try and be somewhat crafty about this getVal stuff to make
  //  these
  //  containers a little bit more generic.  The normal behavior here is that
  //  the
  //  value being queried must be directly castable to type T.  We'll robustify
  //  a
  //  little bit by trying that and, if the cast fails, attempting a couple of
  //  other casts, which will then be lexically cast to type T.
  //
  for(size_t i=0; i< _data.size(); ++i) {
    if (_data[i].getKey() == tag) {
      rdvalue_tostring(_data[i].val, res);
      return;
    }

  }
  throw KeyErrorException(common_properties::getPropName(tag));    
}

bool Dict::getValIfPresent(int tag, std::string &res) const {
  //
  //  We're going to try and be somewhat crafty about this getVal stuff to make
  //  these
  //  containers a little bit more generic.  The normal behavior here is that
  //  the
  //  value being queried must be directly castable to type T.  We'll robustify
  //  a
  //  little bit by trying that and, if the cast fails, attempting a couple of
  //  other casts, which will then be lexically cast to type T.
  //
  for(size_t i=0; i< _data.size(); ++i) {
    if (_data[i].getKey() == tag) {
      rdvalue_tostring(_data[i].val, res);
      return true;
    }
  }
  return false;
}


namespace common_properties {

const char *getPropName(int v) {
  if (v>=0 && v<=MAX) {
    return common_properties::propnames[v];
  }
  if ((size_t)v < RDKit::KeyIntPair::tagmap.keys.size())
    return RDKit::KeyIntPair::tagmap.keys[v].c_str();

  throw KeyErrorException("Unknown tag");
}

}

}
