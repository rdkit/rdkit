// $Id$
//
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Dict.h"

namespace RDKit {

void Dict::getVal(const std::string &what, std::string &res) const {
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
  for (const auto &i : _data) {
    if (i.key == what) {
      rdvalue_tostring(i.val, res);
      return;
    }

  }
  throw KeyErrorException(what);    
}

bool Dict::getValIfPresent(const std::string &what, std::string &res) const {
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
  for (const auto &i : _data) {
    if (i.key == what) {
      rdvalue_tostring(i.val, res);
      return true;
    }
  }
  return false;
}

}
