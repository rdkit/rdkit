//  $Id: MolOps.cpp 5061 2006-03-08 00:36:29Z glandrum $
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "RDLib.h"
#include <string>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>

extern "C" {
  PG_FUNCTION_INFO_V1(rd_canonsmiles);
   Datum rd_canonsmiles(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_molpickle);
   Datum rd_molpickle(PG_FUNCTION_ARGS);

  PG_FUNCTION_INFO_V1(rd_hassubstruct_smi);
   Datum rd_hassubstruct_smi(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_hassubstruct_pkl);
   Datum rd_hassubstruct_pkl(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_substructcount_smi);
   Datum rd_substructcount_smi(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_substructcount_pkl);
   Datum rd_substructcount_pkl(PG_FUNCTION_ARGS);

  PG_FUNCTION_INFO_V1(rd_substructfp);
   Datum rd_substructfp(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_substructfp_bits);
   Datum rd_substructfp_bits(PG_FUNCTION_ARGS);

  PG_FUNCTION_INFO_V1(rd_similarityfp);
   Datum rd_similarityfp(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_similarityfp_bits);
   Datum rd_similarityfp_bits(PG_FUNCTION_ARGS);

  PG_FUNCTION_INFO_V1(rd_mollogp_smi);
   Datum rd_mollogp_smi(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_mollogp_pkl);
   Datum rd_mollogp_pkl(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_amw_smi);
   Datum rd_amw_smi(PG_FUNCTION_ARGS);
  PG_FUNCTION_INFO_V1(rd_amw_pkl);
   Datum rd_amw_pkl(PG_FUNCTION_ARGS);


}

Datum rd_canonsmiles(PG_FUNCTION_ARGS)
{
  std::string res="";
  if(!PG_ARGISNULL(0)){
    text *smi = PG_GETARG_TEXT_P(0);
    RDKit::RWMol *mol=0;
    std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);
    try {
      mol=RDKit::SmilesToMol(smiS);
    } catch (...) {
      mol=0;
    }
    if(mol) {
      try {
	res=RDKit::MolToSmiles(*mol,true);
      } catch (...) {
	ereport(ERROR,(errcode(ERRCODE_INTERNAL_ERROR),
		       errmsg("molecule constructed from SMILES '%s' could not be converted back to SMILES",
			      smiS.c_str())));
      }
      delete mol;
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("SMILES '%s' could not be converted into a molecule",
			    smiS.c_str())));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty SMILES in rd_canonsmiles function call")));
  }
  text *textRes;
  int totalSz=res.size()+VARHDRSZ;
  textRes = (text *)palloc(totalSz);
  VARATT_SIZEP(textRes)=totalSz;
  memcpy(VARDATA(textRes),res.c_str(),res.size());
  PG_RETURN_TEXT_P(textRes);
}


Datum rd_molpickle(PG_FUNCTION_ARGS)
{
  StringInfoData buf;
  pq_begintypsend(&buf);
  std::string pkl="";
  if(!PG_ARGISNULL(0)){
    text *smi = PG_GETARG_TEXT_P(0);
    RDKit::RWMol *mol=0;
    std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);
    try {
      mol=RDKit::SmilesToMol(smiS);
    } catch (...) {
      mol=0;
    }
    if(mol) {
      RDKit::MolPickler::pickleMol(*mol,pkl);
      delete mol;
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("SMILES '%s' could not be converted into a molecule",
			    smiS.c_str())));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty SMILES in rd_molpickle function call")));
  }
  pq_sendtext(&buf,pkl.c_str(),pkl.size());
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}


Datum rd_hassubstruct_smi(PG_FUNCTION_ARGS)
{
  if(PG_ARGISNULL(0) || PG_ARGISNULL(1)) {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_hassubstruct_smi function call")));
  }
  text *sma = PG_GETARG_TEXT_P(0);
  text *smi = PG_GETARG_TEXT_P(1);

  RDKit::RWMol *probe=0;

  std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);
  std::string smaS(VARDATA(sma),VARSIZE(sma)-VARHDRSZ);
  bool res = false;
  RDKit::RWMol *mol;
  try {
    mol=RDKit::SmilesToMol(smiS);
  } catch (...) {
    mol=0;
  }
  if(mol) {
    try {
      probe=RDKit::SmartsToMol(smaS);
    } catch (...) {
      probe=0;
    }
    if(probe){
      RDKit::MatchVectType matchVect;
      res = RDKit::SubstructMatch(*mol,*probe,matchVect);
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("SMILES/SMARTS '%s' could not be converted into a query",
			    smaS.c_str())));
    }
    delete mol;
  } else {
    ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		   errmsg("SMILES '%s' could not be converted into a molecule",
			  smiS.c_str())));
  }
  PG_RETURN_BOOL(res);
}



Datum rd_hassubstruct_pkl(PG_FUNCTION_ARGS)
{
  bool res = false;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    bytea *probePkl=PG_GETARG_BYTEA_P(0);
    bytea *refPkl=PG_GETARG_BYTEA_P(1);
    RDKit::ROMol *probe=0,*ref=0;
    try {
      probe = new RDKit::ROMol(std::string(VARDATA(probePkl),VARSIZE(probePkl)-VARHDRSZ));
    } catch (...) {
      probe=0;
    }
    try {
      ref = new RDKit::ROMol(std::string(VARDATA(refPkl),VARSIZE(refPkl)-VARHDRSZ));
    } catch (...) {
      ref=0;
    }
    if(probe && ref){
      RDKit::MatchVectType matchVect;
      res = RDKit::SubstructMatch(*ref,*probe,matchVect);
    } else {
      if(!probe){
	ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		       errmsg("problems converting argument 1 to a molecule")));
      }
      if(!ref){
	ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		       errmsg("problems converting argument 2 to a molecule")));
      }
    }
    if(probe) delete probe;
    if(ref) delete ref;
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_hassubstruct_pkl function call")));
  }
  //elog(INFO,"%d %s",res,smiS.c_str());
  PG_RETURN_BOOL(res);
}


Datum rd_substructcount_smi(PG_FUNCTION_ARGS)
{
  if(PG_ARGISNULL(0) || PG_ARGISNULL(1)) {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_substructcount_smi function call")));
  }
  text *sma = PG_GETARG_TEXT_P(0);
  text *smi = PG_GETARG_TEXT_P(1);

  RDKit::RWMol *probe=0;

  std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);
  std::string smaS(VARDATA(sma),VARSIZE(sma)-VARHDRSZ);
  int res = 0;
  RDKit::RWMol *mol;
  try {
    mol=RDKit::SmilesToMol(smiS);
  } catch (...) {
    mol=0;
  }
  if(mol) {
    try {
      probe=RDKit::SmartsToMol(smaS);
    } catch (...) {
      probe=0;
    }
    if(probe){
      std::vector<RDKit::MatchVectType> matchVect;
      res = RDKit::SubstructMatch(*mol,*probe,matchVect,true);
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("SMILES/SMARTS '%s' could not be converted into a query",
			    smaS.c_str())));
    }
    delete mol;
  } else {
    ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		   errmsg("SMILES '%s' could not be converted into a molecule",
			  smiS.c_str())));
  }
  PG_RETURN_INT32(res);
}



Datum rd_substructcount_pkl(PG_FUNCTION_ARGS)
{
  int res = 0;
  if(!PG_ARGISNULL(0) && !PG_ARGISNULL(1)) {
    bytea *probePkl=PG_GETARG_BYTEA_P(0);
    bytea *refPkl=PG_GETARG_BYTEA_P(1);
    RDKit::ROMol *probe=0,*ref=0;
    try {
      probe = new RDKit::ROMol(std::string(VARDATA(probePkl),VARSIZE(probePkl)-VARHDRSZ));
    } catch (...) {
      probe=0;
    }
    try {
      ref = new RDKit::ROMol(std::string(VARDATA(refPkl),VARSIZE(refPkl)-VARHDRSZ));
    } catch (...) {
      ref=0;
    }
    if(probe && ref){
      std::vector<RDKit::MatchVectType> matchVect;
      res = RDKit::SubstructMatch(*ref,*probe,matchVect,true);
    } else {
      if(!probe){
	ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		       errmsg("problems converting argument 1 to a molecule")));
      }
      if(!ref){
	ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		       errmsg("problems converting argument 2 to a molecule")));
      }
    }
    if(probe) delete probe;
    if(ref) delete ref;
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_substructcount_pkl function call")));
  }
  //elog(INFO,"%d %s",res,smiS.c_str());
  PG_RETURN_INT32(res);
}


Datum rd_substructfp(PG_FUNCTION_ARGS)
{
  StringInfoData buf;
  pq_begintypsend(&buf);
  ExplicitBitVect *bv=0;
  int fpSize=2048;
  if(!PG_ARGISNULL(0)){
    text *smi = PG_GETARG_TEXT_P(0);
    std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

    RDKit::RWMol *mol;
    try {
      mol=RDKit::SmilesToMol(smiS);
    } catch (...) {
      mol=0;
    }
    if(mol) {
      try{
	bv = DaylightFingerprintMol(*mol,1,7,fpSize);
      } catch (...) {
	bv = new ExplicitBitVect(fpSize);
      }
      delete mol;
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("SMILES '%s' could not be converted into a molecule",
			    smiS.c_str())));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty SMILES in rd_substructfp function call")));
  }
  std::string pkl=bv->ToString();
  delete bv;
  pq_sendtext(&buf,pkl.c_str(),pkl.size());
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}


Datum rd_substructfp_bits(PG_FUNCTION_ARGS)
{
  int fpSize=2048;

  ExplicitBitVect *bv=0;
  if(!PG_ARGISNULL(0)){
    text *smi = PG_GETARG_TEXT_P(0);
    std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

    RDKit::RWMol *mol;
    try {
      mol=RDKit::SmilesToMol(smiS);
    } catch (...) {
      mol=0;
    }
    if(mol) {
      try{
	bv = DaylightFingerprintMol(*mol,1,7,fpSize);
      } catch (...) {
	bv = 0;
      }
      delete mol;
    } else {
      bv = 0;
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty SMILES in rd_substructfp_bits function call")));
  }
  // BVToVarBits handles null pointers w/o problems so we can go ahead and
  // call it directly:
  VarBit *result = BVToVarBits(bv,fpSize);
  if(bv){
    delete bv;
  }
  PG_RETURN_VARBIT_P(result);
}

Datum rd_similarityfp(PG_FUNCTION_ARGS)
{
  StringInfoData buf;
  pq_begintypsend(&buf);
  ExplicitBitVect *bv=0;
  int fpSize=2048;
  if(!PG_ARGISNULL(0)){
    text *smi = PG_GETARG_TEXT_P(0);
    std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

    double tgtDensity=0.3;
    RDKit::RWMol *mol;
    try {
      mol=RDKit::SmilesToMol(smiS);
    } catch (...) {
      mol=0;
    }
    if(mol) {
      try{
	bv = DaylightFingerprintMol(*mol,1,7,fpSize);
	double density=static_cast<double>(bv->GetNumOnBits())/bv->GetNumBits();
	while(density<tgtDensity && bv->GetNumBits()>8){
	  ExplicitBitVect *tmp=FoldFingerprint(*bv,2);
	  delete bv;
	  bv = tmp;
	  density=static_cast<double>(bv->GetNumOnBits())/bv->GetNumBits();
	}
      } catch (...) {
	bv = new ExplicitBitVect(fpSize);
      }
      delete mol;
    } else {
      bv = new ExplicitBitVect(fpSize);
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty SMILES in rd_similarityfp function call")));
  }
  std::string pkl=bv->ToString();
  delete bv;
  pq_sendtext(&buf,pkl.c_str(),pkl.size());
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}


Datum rd_similarityfp_bits(PG_FUNCTION_ARGS)
{
  int fpSize=2048;
  ExplicitBitVect *bv=0;

  if(!PG_ARGISNULL(0)){
    text *smi = PG_GETARG_TEXT_P(0);
    std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

    double tgtDensity=0.3;
    
    RDKit::RWMol *mol;
    try {
      mol=RDKit::SmilesToMol(smiS);
    } catch (...) {
      mol=0;
    }
    if(mol) {
      try{
	bv = DaylightFingerprintMol(*mol,1,7,fpSize);
	double density=static_cast<double>(bv->GetNumOnBits())/bv->GetNumBits();
	while(density<tgtDensity && bv->GetNumBits()>8){
	  ExplicitBitVect *tmp=FoldFingerprint(*bv,2);
	  delete bv;
	  bv = tmp;
	  density=static_cast<double>(bv->GetNumOnBits())/bv->GetNumBits();
	}
      } catch (...) {
	bv = new ExplicitBitVect(fpSize);
      }
      delete mol;
    } else {
      bv = new ExplicitBitVect(fpSize);
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty SMILES in rd_similarityfp_bits function call")));
  }
  // BVToVarBits handles null pointers w/o problems so we can go ahead and
  // call it directly:
  if(bv) fpSize=bv->GetNumBits();
  VarBit *result = BVToVarBits(bv,fpSize);
  if(bv){
    delete bv;
  }
  PG_RETURN_VARBIT_P(result);
}


Datum rd_mollogp_smi(PG_FUNCTION_ARGS)
{
  //elog(NOTICE,"rd_mollogp_smi");
  if(PG_ARGISNULL(0)){
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_mollogp_smi function call")));
  }

  text *smi = PG_GETARG_TEXT_P(0);
  std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

  
  double res = 0.0;
  RDKit::RWMol *mol;
  try {
    mol=RDKit::SmilesToMol(smiS);
    
  } catch (...) {
    mol=0;
  }
  if(mol) {
    double tmp;
    RDKit::Descriptors::CalcCrippenDescriptors(*mol,res,tmp);
    delete mol;
  } else {
    ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		   errmsg("SMILES '%s' could not be converted into a molecule",
			  smiS.c_str())));
  }
  PG_RETURN_FLOAT8(res);
}



Datum rd_mollogp_pkl(PG_FUNCTION_ARGS)
{
  //elog(NOTICE,"rd_mollogp_pkl");
  double res=0;
  if(!PG_ARGISNULL(0)){
    bytea *molPkl=PG_GETARG_BYTEA_P(0);
    RDKit::ROMol *mol=0;
    try {
      mol = new RDKit::ROMol(std::string(VARDATA(molPkl),VARSIZE(molPkl)-VARHDRSZ));
    } catch (...) {
      mol=0;
    }
    if(mol){
      //elog(NOTICE,"   numatoms: %d",mol->getNumAtoms());
      double tmp;
      RDKit::Descriptors::CalcCrippenDescriptors(*mol,res,tmp);
      delete mol;
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("problems converting argument to a molecule in rd_mollogp_pkl")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_mollogp_pkl function call")));
  }
  PG_RETURN_FLOAT8(res);
}

Datum rd_amw_smi(PG_FUNCTION_ARGS)
{
  if(PG_ARGISNULL(0)){
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_amw_smi function call")));
  }

  text *smi = PG_GETARG_TEXT_P(0);
  std::string smiS(VARDATA(smi),VARSIZE(smi)-VARHDRSZ);

  
  double res = 0.0;
  RDKit::RWMol *mol;
  try {
    mol=RDKit::SmilesToMol(smiS);
    
  } catch (...) {
    mol=0;
  }
  if(mol) {
    res = RDKit::Descriptors::CalcAMW(*mol);
    delete mol;
  } else {
    ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		   errmsg("SMILES '%s' could not be converted into a molecule",
			  smiS.c_str())));
  }
  PG_RETURN_FLOAT8(res);
}



Datum rd_amw_pkl(PG_FUNCTION_ARGS)
{
  double res=0;
  if(!PG_ARGISNULL(0)){
    bytea *molPkl=PG_GETARG_BYTEA_P(0);
    RDKit::ROMol *mol=0;
    try {
      mol = new RDKit::ROMol(std::string(VARDATA(molPkl),VARSIZE(molPkl)-VARHDRSZ));
    } catch (...) {
      mol=0;
    }
    if(mol){
      res = RDKit::Descriptors::CalcAMW(*mol);
      delete mol;
    } else {
      ereport(ERROR,(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
		     errmsg("problems converting argument to a molecule in rd_amw_pkl")));
    }
  } else {
    ereport(ERROR,(errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
		   errmsg("missing or empty argument in rd_amw_pkl function call")));
  }
  PG_RETURN_FLOAT8(res);
}

