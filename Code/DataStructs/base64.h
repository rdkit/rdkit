//
// Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_BASE64_H__
#define __RD_BASE64_H__
/*! \file base64.h

  \brief Functionality for base64 encoding/decoding

*/

//! return the base64 encoding of an array of unsigned chars
/*!
   <b>Note:</b> The caller is responsible for calling \c delete[] on the
     char array returned by this function.
 */
char *Base64Encode(const unsigned char *,const unsigned int);

//! return the base64 encoding of an array of chars
/*!
   <b>Note:</b> The caller is responsible for calling \c delete[] on the
     char array returned by this function.
 */
char *Base64Encode(const char *,const unsigned int);

//! return the decoded version of a base64 encoded char array
/*!
   <b>Note:</b> The caller is responsible for calling \c delete[] on the
     char array returned by this function.
 */
char *Base64Decode(const char *,unsigned int *);

#endif
