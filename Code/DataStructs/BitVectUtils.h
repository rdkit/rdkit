//
// Copyright (c) 2002-2006 greg Landrum, Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_BITVECTS_UTILS_H__
#define __RD_BITVECTS_UTILS_H__

#include "BitVects.h"
#include <string>

//! \brief Construct a BitVect from the ASCII representation of a
//! Daylight fingerprint string
template <typename T>
void FromDaylightString(T &sbv,std::string s);

//! \brief Construct a BitVect from the ASCII representation of a
//! bit string (i.e. a bunch of zeros and ones)
template <typename T>
void FromBitString(T &sbv,const std::string &s);


//! Convert a SparseBitVector to an ExplicitBitVector
/*!
  \return a pointer to an ExplicitBitVector
  <b>Note:</b> the caller is responsible for <tt>delete</tt>ing this.

 */
ExplicitBitVect *convertToExplicit(const SparseBitVect *sbv);
#endif
