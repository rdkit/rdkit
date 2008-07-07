//  $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "RDLib.h"
#include <cmath>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>

extern "C" {
  PG_FUNCTION_INFO_V1(rd_fpsize);
   Datum rd_fpsize(PG_FUNCTION_ARGS);

  PG_FUNCTION_INFO_V1(rd_allprobebitsmatch_pkl);
   Datum rd_allprobebitsmatch_pkl(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_allprobebitsmatch_bits);
   Datum rd_allprobebitsmatch_bits(PG_FUNCTION_ARGS);

  PG_FUNCTION_INFO_V1(rd_tanimoto_pkl);
   Datum rd_tanimoto_pkl(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_dice_pkl);
   Datum rd_dice_pkl(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_cosine_pkl);
   Datum rd_cosine_pkl(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_tanimoto_bits);
   Datum rd_tanimoto_bits(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_dice_bits);
   Datum rd_dice_bits(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_cosine_bits);
   Datum rd_cosine_bits(PG_FUNCTION_ARGS);
}
extern char numByteOnBits[];
extern char byteByteDistMat[];


Datum rd_fpsize(PG_FUNCTION_ARGS)
{
  int res=0;
  if(!PG_ARGISNULL(0)) {
    bytea *buf=PG_GETARG_BYTEA_P(0);
    try{
      ExplicitBitVect bv(VARDATA(buf),VARSIZE(buf)-VARHDRSZ);
      res = bv.GetNumBits();
    } catch(...) {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("Could not construct bit vector from argument")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }
  PG_RETURN_INT32(res);

}

Datum rd_allprobebitsmatch_pkl(PG_FUNCTION_ARGS)
{
  bool res=false;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    bytea *probe=PG_GETARG_BYTEA_P(0);
    bytea *ref=PG_GETARG_BYTEA_P(1);

    try{
      res = AllProbeBitsMatch((const char *)VARDATA(probe),
    			      (const char *)VARDATA(ref));
    } catch(...) {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("Could not compare bit vectors from arguments")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }
  PG_RETURN_BOOL(res);

}

Datum rd_tanimoto_pkl(PG_FUNCTION_ARGS)
{
  ExplicitBitVect *bv1=0,*bv2=0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    bytea *buf1=PG_GETARG_BYTEA_P(0);
    bytea *buf2=PG_GETARG_BYTEA_P(1);
    try{
      bv1 = new ExplicitBitVect(VARDATA(buf1),VARSIZE(buf1)-VARHDRSZ);
      bv2 = new ExplicitBitVect(VARDATA(buf2),VARSIZE(buf2)-VARHDRSZ);
    } catch(...) {
      if(bv1) delete bv1;
      if(bv2) delete bv2;
      bv1=0;
      bv2=0;
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("Could not construct bit vectors from arguments")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }

  double res=0.0;
  if(bv1 && bv2) {
    if(bv1->GetNumBits() > bv2->GetNumBits()){
      unsigned int foldFact = bv1->GetNumBits() / bv2->GetNumBits();
      ExplicitBitVect *tmp = FoldFingerprint(*bv1,foldFact);
      delete bv1;
      bv1 = tmp;
    } else if (bv2->GetNumBits() > bv1->GetNumBits()) {
      unsigned int foldFact = bv2->GetNumBits() / bv1->GetNumBits();
      ExplicitBitVect *tmp = FoldFingerprint(*bv2,foldFact);
      delete bv2;
      bv2 = tmp;
    }
    res = TanimotoSimilarity(*bv1,*bv2);

    delete bv1;
    delete bv2;
  }

  PG_RETURN_FLOAT8(res);
}

Datum rd_dice_pkl(PG_FUNCTION_ARGS)
{
  ExplicitBitVect *bv1=0,*bv2=0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    bytea *buf1=PG_GETARG_BYTEA_P(0);
    bytea *buf2=PG_GETARG_BYTEA_P(1);
    try{
      bv1 = new ExplicitBitVect(VARDATA(buf1),VARSIZE(buf1)-VARHDRSZ);
      bv2 = new ExplicitBitVect(VARDATA(buf2),VARSIZE(buf2)-VARHDRSZ);
    } catch(...) {
      if(bv1) delete bv1;
      if(bv2) delete bv2;
      bv1=0;
      bv2=0;
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("Could not construct bit vectors from arguments")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }

  double res=0.0;
  if(bv1 && bv2) {
    if(bv1->GetNumBits() > bv2->GetNumBits()){
      int foldFact = bv1->GetNumBits() / bv2->GetNumBits();
      ExplicitBitVect *tmp = FoldFingerprint(*bv1,foldFact);
      delete bv1;
      bv1 = tmp;
    } else if (bv2->GetNumBits() > bv1->GetNumBits()) {
      int foldFact = bv2->GetNumBits() / bv1->GetNumBits();
      ExplicitBitVect *tmp = FoldFingerprint(*bv2,foldFact);
      delete bv2;
      bv2 = tmp;
    }
    res = DiceSimilarity(*bv1,*bv2);

    delete bv1;
    delete bv2;
  }

  PG_RETURN_FLOAT8(res);
}

Datum rd_cosine_pkl(PG_FUNCTION_ARGS)
{
  ExplicitBitVect *bv1=0,*bv2=0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    bytea *buf1=PG_GETARG_BYTEA_P(0);
    bytea *buf2=PG_GETARG_BYTEA_P(1);
    try{
      bv1 = new ExplicitBitVect(VARDATA(buf1),VARSIZE(buf1)-VARHDRSZ);
      bv2 = new ExplicitBitVect(VARDATA(buf2),VARSIZE(buf2)-VARHDRSZ);
    } catch(...) {
      if(bv1) delete bv1;
      if(bv2) delete bv2;
      bv1=0;
      bv2=0;
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("Could not construct bit vectors from arguments")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }

  double res=0.0;
  if(bv1 && bv2) {
    if(bv1->GetNumBits() > bv2->GetNumBits()){
      int foldFact = bv1->GetNumBits() / bv2->GetNumBits();
      ExplicitBitVect *tmp = FoldFingerprint(*bv1,foldFact);
      delete bv1;
      bv1 = tmp;
    } else if (bv2->GetNumBits() > bv1->GetNumBits()) {
      int foldFact = bv2->GetNumBits() / bv1->GetNumBits();
      ExplicitBitVect *tmp = FoldFingerprint(*bv2,foldFact);
      delete bv2;
      bv2 = tmp;
    }
    res = CosineSimilarity(*bv1,*bv2);

    delete bv1;
    delete bv2;
  }

  PG_RETURN_FLOAT8(res);
}

Datum rd_allprobebitsmatch_bits(PG_FUNCTION_ARGS)
{
  bool res=false;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    VarBit *probe=PG_GETARG_VARBIT_P(0);
    VarBit *ref=PG_GETARG_VARBIT_P(1);
    int bitlen1 = VARBITLEN(probe);
    if(bitlen1 != VARBITLEN(ref)){
      res = false;
    } else {
      bits8 *p1,*p2;
      p1 = VARBITS(probe);
      p2 = VARBITS(ref);
      res = true;
      for(int i=0;i<VARBITBYTES(probe);i++){
	if(*p1!=(*p1 & *p2)){
	  res=false;
	  break;
	}
	p1++;
	p2++;
      }
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }
  PG_RETURN_BOOL(res);

}

void getBitsInfo(VarBit *probe,VarBit *ref,int &numDifferent,
		 int &numA,int &numB){

  numDifferent=0;
  numA=0;
  numB=0;
  if(!probe || !ref) return;
  
  // make sure the probe isn't longer than the ref:
  if(VARBITLEN(probe)>VARBITLEN(ref)){
    VarBit *tmp=probe;
    probe=ref;
    ref=tmp;
  }

  int foldMultiplier=VARBITLEN(ref)/VARBITLEN(probe);
    
  for(int i=0;i<VARBITBYTES(probe);i++){
    bits8 v1 = VARBITS(probe)[i];
    bits8 v2 = VARBITS(ref)[i];
    for(int j=1;j<foldMultiplier;j++){
      v2 |= VARBITS(ref)[i+j*VARBITBYTES(probe)];
    }
    numDifferent += byteByteDistMat[v1*256+v2];
    numA += numByteOnBits[v1];
    numB += numByteOnBits[v2];
  }
}

Datum rd_tanimoto_bits(PG_FUNCTION_ARGS)
{
  double res=0.0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    int numDifferent=0,numA=0,numB=0;
    
    VarBit *probe=PG_GETARG_VARBIT_P(0);
    VarBit *ref=PG_GETARG_VARBIT_P(1);

    getBitsInfo(probe,ref,numDifferent,numA,numB);

    double denom = numA + numB + numDifferent;
    if(denom > 0){
      res = (numA + numB - numDifferent)/denom;
    } else {
      res = 0.0;
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }
  PG_RETURN_FLOAT8(res);

}

Datum rd_cosine_bits(PG_FUNCTION_ARGS)
{
  double res=0.0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    int numDifferent=0,numA=0,numB=0;
    
    VarBit *probe=PG_GETARG_VARBIT_P(0);
    VarBit *ref=PG_GETARG_VARBIT_P(1);

    getBitsInfo(probe,ref,numDifferent,numA,numB);

    double denom = sqrt(numA*numB);
    if(denom > 0){
      res = .5*(numA+numB-numDifferent)/denom;
    } else {
      res = 0.0;
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }
  PG_RETURN_FLOAT8(res);

}


Datum rd_dice_bits(PG_FUNCTION_ARGS)
{
  double res=0.0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    int numDifferent=0,numA=0,numB=0;
    
    VarBit *probe=PG_GETARG_VARBIT_P(0);
    VarBit *ref=PG_GETARG_VARBIT_P(1);

    getBitsInfo(probe,ref,numDifferent,numA,numB);

    double denom = numA+numB;
    if(denom > 0){
      res = 2.*0.5*(numA+numB-numDifferent)/denom;
    } else {
      res = 0.0;
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in function call")));
  }
  PG_RETURN_FLOAT8(res);

}
