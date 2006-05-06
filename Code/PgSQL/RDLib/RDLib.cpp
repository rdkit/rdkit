//  $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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

