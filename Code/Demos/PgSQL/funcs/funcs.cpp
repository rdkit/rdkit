//  $Id$
//
//  Copyright (C) 2005 Rational Discovery LLC
//   All Rights Reserved
//

extern "C" {
#include "postgres.h"
#include "fmgr.h"				/* for argument/result macros */
#include "libpq/pqformat.h"		/* needed for send/recv functions */
#include <stdlib.h>
}
#ifdef gettext
#undef gettext
#endif
#include <iostream>


#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>


/* These prototypes just prevent possible warnings from gcc. */

extern "C" {
  PG_FUNCTION_INFO_V1(substruct);
  Datum substruct(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(substructfp);
  Datum substructfp(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(fpsize);
  Datum fpsize(PG_FUNCTION_ARGS);
}


Datum substruct(PG_FUNCTION_ARGS)
{
  text *smi = PG_GETARG_TEXT_P(0);
  text *sma = PG_GETARG_TEXT_P(1);

  RDKit::RWMol *probe=0;

  std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);
  std::string smaS(VARDATA(sma),VARSIZE(sma)-VARHDRSZ);
  bool res = false;
  //elog(INFO,"smi(%d,%d,%d): %s\n",smi->vl_len,VARSIZE(smi),VARHDRSZ,smiS.c_str());
  RDKit::RWMol *mol=RDKit::SmilesToMol(smiS);
  if(mol) {
    probe=RDKit::SmartsToMol(smaS);
    if(probe){
      RDKit::MatchVectType matchVect;
      res = RDKit::SubstructMatch(mol,probe,matchVect);
    }
    delete mol;
  }
  //elog(INFO,"%d %s",res,smiS.c_str());
  PG_RETURN_INT32(res);
}

Datum substructfp(PG_FUNCTION_ARGS)
{
  text *smi = PG_GETARG_TEXT_P(0);
  int fpSize=2048;
  std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

  ExplicitBitVect *bv=0;
  RDKit::RWMol *mol=RDKit::SmilesToMol(smiS);
  if(mol) {
    try{
      bv = DaylightFingerprintMol(mol,1,7,fpSize);
    } catch (...) {
      bv = new ExplicitBitVect(fpSize);
    }
    delete mol;
  } else {
    bv = new ExplicitBitVect(fpSize);
  }
  std::string pkl=bv->ToString();
  delete bv;

  StringInfoData buf;
  pq_begintypsend(&buf);
  pq_sendtext(&buf,pkl.c_str(),pkl.size());
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

Datum fpsize(PG_FUNCTION_ARGS)
{
  int res=0;
  if(!PG_ARGISNULL(0)) {
    bytea *buf=PG_GETARG_BYTEA_P(0);
    std::string pkl(VARDATA(buf),VARSIZE(buf)-VARHDRSZ);
    ExplicitBitVect bv(pkl);
    res = bv.GetNumBits();
  }
  PG_RETURN_INT32(res);

}
