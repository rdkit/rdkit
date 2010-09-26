//  $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDLib.h"
#include <string>

extern "C" {
  PG_FUNCTION_INFO_V1(rd_libversion);
   Datum rd_libversion(PG_FUNCTION_ARGS);
}

Datum rd_libversion(PG_FUNCTION_ARGS)
{
  std::string res="$Rev$";
  text *textRes;
  int totalSz=res.size()+VARHDRSZ;
  textRes = (text *)palloc(totalSz);
  VARATT_SIZEP(textRes)=totalSz;
  memcpy(VARDATA(textRes),res.c_str(),res.size());
  PG_RETURN_TEXT_P(textRes);
}

