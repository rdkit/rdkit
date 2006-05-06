//
//  Copyright (C) 2005 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

extern "C" {
#include "postgres.h"
#include "fmgr.h"				/* for argument/result macros */
#include "libpq/pqformat.h"		/* needed for send/recv functions */
#include "utils/varbit.h"		/* needed for varbit functions */
#include <stdlib.h>
}
#ifdef gettext
#undef gettext
#endif
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

template <typename T>
VarBit *BVToVarBits(const T *bv,int fpSize){
  VarBit *result;
  int len = VARBITTOTALLEN(fpSize);
  // initialize to zero:
  result = (VarBit *) palloc0(len);
  VARATT_SIZEP(result) = len;
  VARBITLEN(result) = fpSize;
  
  if(bv){
    int i=0;
    bits8 mask = BITHIGH;
    bits8 *r=VARBITS(result);
    while(i<fpSize){
      if(bv->GetBit(i)){
	*r |= mask;
      }
      mask >>=1;
      if(mask==0){
	mask=BITHIGH;
	++r;
      }
      ++i;
    }
  }
  return result;
}

