// $Id$
// 
// Copyright (C) 2007,2008 Greg Landrum
//
// @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <sqlite3ext.h>
SQLITE_EXTENSION_INIT1
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <boost/cstdint.hpp>
#include <string>
#include <map>

std::string stringFromTextArg(sqlite3_value *arg){
  const unsigned char *text=sqlite3_value_text(arg);
  int nBytes=sqlite3_value_bytes(arg);
  std::string res((const char *)text,nBytes);
  return res;
}

std::string stringFromBlobArg(sqlite3_value *arg){
  const void *blob=sqlite3_value_blob(arg);
  int nBytes=sqlite3_value_bytes(arg);
  std::string res((const char *)blob,nBytes);
  return res;
}

RDKit::ROMol *molFromBlobArg(sqlite3_value *arg){
  std::string pkl=stringFromBlobArg(arg);
  RDKit::ROMol *m;
  try{
    m = new RDKit::ROMol(pkl);
  } catch (RDKit::MolPicklerException &){
    m=0;
  }
  return m;
}


ExplicitBitVect *ebvFromBlobArg(sqlite3_value *arg){
  std::string pkl=stringFromBlobArg(arg);
  ExplicitBitVect *ebv;
  try{
    ebv = new ExplicitBitVect(pkl);
  } catch (ValueErrorException &){
    ebv=0;
  }
  return ebv;
}

template <typename T>
RDKit::SparseIntVect<T> *sivFromBlobArg(sqlite3_value *arg){
  std::string pkl=stringFromBlobArg(arg);
  RDKit::SparseIntVect<T> *siv;
  try{
    siv = new RDKit::SparseIntVect<T>(pkl);
  } catch (ValueErrorException &){
    siv=0;
  }
  return siv;
}


/* ---------------------------------

  Benchmarking results.

    Database: 65385 pubchem compounds

    Simple access: select count(*) from molecules where length(molpkl)>40;
                   0.3s
    depickle     : select count(*) from molecules where rdk_molNumAtoms(molpkl)>40;
                   11.3s
 
    substruct1   : select count(*) from molecules where 
                   rdk_molHasSubstruct(molpkl,'c1ncncn1');
                   18.0s

    substruct2   : select count(*) from molecules where 
                   rdk_molHasSubstruct(molpkl,'[#6;r10]');
                   15.8

    3 Oct 2007:
    depickle     : select count(*) from molecules where rdk_molNumAtoms(molpkl)>40;
                   9.4s
    mw           : select count(*) from molecules where rdk_molAMW(molpkl)<200;
                   9.7s


		   
 --------------------------------- */


static void numAtomsFunc(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(m){
    int res=m->getNumAtoms();
    delete m;
    sqlite3_result_int(context, res);
  } else {
    std::string errorMsg="BLOB could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
  }
}

static void molWtFunc(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(m){
    double res=RDKit::Descriptors::CalcAMW(*m);
    delete m;
    sqlite3_result_double(context, res);
  } else {
    std::string errorMsg="BLOB could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
  }
}

static void molLogPFunc(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(m){
    double res,tmp;
    RDKit::Descriptors::CalcCrippenDescriptors(*m,res,tmp);
    delete m;
    sqlite3_result_double(context, res);
  } else {
    std::string errorMsg="BLOB could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
  }
}

static void smilesToBlob(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  std::string smiles=stringFromTextArg(argv[0]);
  RDKit::ROMol *m=0;
  try{
    m=RDKit::SmilesToMol(smiles);
  } catch(RDKit::MolSanitizeException &){
    m=0;
  }
  if(m){
    std::string text;
    RDKit::MolPickler::pickleMol(*m,text);
    delete m;
    sqlite3_result_blob(context, text.c_str(), text.length(), SQLITE_TRANSIENT );
  } else {
    std::string errorMsg="SMILES could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
  }
}

static void molHasSubstruct(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(!m){
    std::string errorMsg="BLOB (argument 1) could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }

  std::string smarts=stringFromTextArg(argv[1]);

  std::map<std::string,boost::any> &molMap=
    *static_cast<std::map<std::string,boost::any> *>(sqlite3_user_data(context));
  RDKit::ROMol *patt=0;
  if(molMap.find(smarts)!=molMap.end()){
    patt=boost::any_cast<RDKit::ROMOL_SPTR>(molMap[smarts]).get();
  } else {
    patt=static_cast<RDKit::ROMol *>(RDKit::SmartsToMol(smarts));
    molMap[smarts]=boost::any(RDKit::ROMOL_SPTR(patt));
  }
  if(!patt){
    std::string errorMsg="SMARTS (argument 2) could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  RDKit::MatchVectType match;
  int res=RDKit::SubstructMatch(*m,*patt,match,true,false,true);
  delete m;
  sqlite3_result_int(context, res);
}

static void molSubstructCount(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(!m){
    std::string errorMsg="BLOB (argument 1) could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }

  std::string smarts=stringFromTextArg(argv[1]);

  std::map<std::string,boost::any> &molMap=
    *static_cast<std::map<std::string,boost::any> *>(sqlite3_user_data(context));
  RDKit::ROMol *patt=0;
  if(molMap.find(smarts)!=molMap.end()){
    patt=boost::any_cast<RDKit::ROMOL_SPTR>(molMap[smarts]).get();
  } else {
    patt=static_cast<RDKit::ROMol *>(RDKit::SmartsToMol(smarts));
    molMap[smarts]=boost::any(RDKit::ROMOL_SPTR(patt));
  }
  if(!patt){
    std::string errorMsg="SMARTS (argument 2) could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  std::vector<RDKit::MatchVectType> matches;
  int res=RDKit::SubstructMatch(*m,*patt,matches,true,true,false);
  delete m;
  sqlite3_result_int(context, res);
}

static void blobToRDKitFingerprint(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(!m){
    std::string errorMsg="BLOB (argument 1) could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  ExplicitBitVect *fp=RDKit::DaylightFingerprintMol(*m,1,7,2048,4,true,0.3,128);
  std::string text=fp->toString();
  delete fp;
  delete m;
  sqlite3_result_text(context, text.c_str(), text.length(), SQLITE_TRANSIENT );
}

static void blobToAtomPairFingerprint(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::ROMol *m=molFromBlobArg(argv[0]);
  if(!m){
    std::string errorMsg="BLOB (argument 1) could not be converted into a molecule";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  RDKit::SparseIntVect<int> *fp=RDKit::Descriptors::AtomPairs::getAtomPairFingerprint(*m);
  std::string text=fp->toString();
  delete fp;
  delete m;
  sqlite3_result_text(context, text.c_str(), text.length(), SQLITE_TRANSIENT );
}

static void bvTanimotoSim(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  ExplicitBitVect *bv1=ebvFromBlobArg(argv[0]);
  if(!bv1){
    std::string errorMsg="BLOB (argument 1) could not be converted into a bit vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  ExplicitBitVect *bv2=ebvFromBlobArg(argv[1]);
  if(!bv2){
    delete bv1;
    std::string errorMsg="BLOB (argument 2) could not be converted into a bit vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  double res=SimilarityWrapper(*bv1,*bv2,TanimotoSimilarity);
  delete bv1;
  delete bv2;
  sqlite3_result_double(context, res);
}
static void ucvTanimotoSim(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  // table from Andrew Dalke:
  static const unsigned int popCounts[] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
  };
  const unsigned char *t1=(const unsigned char *)sqlite3_value_blob(argv[0]);
  int nB1=sqlite3_value_bytes(argv[0]);
  const unsigned char *t2=(const unsigned char *)sqlite3_value_blob(argv[1]);
  int nB2=sqlite3_value_bytes(argv[1]);
  if(nB1!=nB2){
    std::string errorMsg="bit vectors not ths same length";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  unsigned int x=0,y=0,z=0;
  for(unsigned int i=0;i<(unsigned int)nB1;++i){
    y+= popCounts[*t1];
    z+= popCounts[*t2];
    x+= popCounts[(*t1)&(*t2)];
    ++t1;
    ++t2;
  }
  double res=0;
  if(y+z-x>0){
    res=double(x) / (y+z-x);
  }
  sqlite3_result_double(context, res);
}

#if 0
// Naive approach: actually construct two sparse int vects:
static void sivDiceSim(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  RDKit::SparseIntVect<int> *v1=sivFromBlobArg<int>(argv[0]);
  if(!v1){
    std::string errorMsg="BLOB (argument 1) could not be converted into an int vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  RDKit::SparseIntVect<int> *v2=sivFromBlobArg<int>(argv[1]);
  if(!v2){
    delete v1;
    std::string errorMsg="BLOB (argument 2) could not be converted into a bit vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  double res= RDKit::DiceSimilarity(*v1,*v2);
  delete v1;
  delete v2;
  sqlite3_result_double(context, res);
}
#else
// faster, just parse the format directly
static void sivDiceSim(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  const unsigned char *t1=(const unsigned char *)sqlite3_value_blob(argv[0]);
  int nB1=sqlite3_value_bytes(argv[0]);
  const unsigned char *t2=(const unsigned char *)sqlite3_value_blob(argv[1]);
  int nB2=sqlite3_value_bytes(argv[1]);

  // check the version flags:
  boost::uint32_t tmp;
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  if(tmp!=ci_SPARSEINTVECT_VERSION){
    std::string errorMsg="BLOB (argument 1) could not be converted into an int vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  if(tmp!=ci_SPARSEINTVECT_VERSION){
    std::string errorMsg="BLOB (argument 2) could not be converted into an int vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }

  // check the element size:
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  if(tmp!=sizeof(boost::uint32_t)){
    std::string errorMsg="BLOB (argument 1) could not be converted into an uint32_t vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  if(tmp!=sizeof(boost::uint32_t)){
    std::string errorMsg="BLOB (argument 2) could not be converted into an uint32_t vector";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }
 
  double res=0.;
  // start reading:
  boost::uint32_t len1,len2;
  len1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  len2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  if(len1!=len2){
    std::string errorMsg="attempt to compare fingerprints of different length";
    sqlite3_result_error(context,errorMsg.c_str(),errorMsg.length());
    return;
  }

  boost::uint32_t nElem1,nElem2;
  nElem1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  nElem2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);

  if(!nElem1 || !nElem2){
    res=0.0;
    sqlite3_result_double(context, res);
  }

  double v1Sum=0,v2Sum=0,numer=0;
  boost::uint32_t idx1=0;
  boost::int32_t v1;
  boost::uint32_t idx2=0;
  boost::int32_t v2;
  idx1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1+=sizeof(boost::uint32_t);
  v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
  t1+=sizeof(boost::int32_t);
  nElem1--;
  v1Sum += v1;

  idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2+=sizeof(boost::uint32_t);
  v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
  t2+=sizeof(boost::int32_t);
  nElem2--;
  v2Sum += v2;

  while(1){
    while(nElem2 && idx2<idx1){
      idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
      t2+=sizeof(boost::uint32_t);
      v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
      t2+=sizeof(boost::int32_t);
      nElem2--;
      v2Sum += v2;
    }
    if(idx2==idx1 ){
      //std::cerr<<"   --- "<<idx1<<" "<<v1<<" - "<<idx2<<" "<<v2<<std::endl;
      numer += std::min(v1,v2);
    }
    if(nElem1){
      idx1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
      t1+=sizeof(boost::uint32_t);
      v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
      t1+=sizeof(boost::int32_t);
      nElem1--;
      v1Sum += v1;
    } else {
      break;
    }
  }
  while(nElem2){
    idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
    t2+=sizeof(boost::uint32_t);
    v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
    t2+=sizeof(boost::int32_t);
    nElem2--;
    v2Sum += v2;
  }
  double denom=v1Sum+v2Sum;
  if(fabs(denom)<1e-6){
    res=0.0;
  } else {
    res = 2.*numer/denom;
  }
  //std::cerr<<" "<<v1Sum<<" "<<v2Sum<<" "<<numer<<" "<<res<<std::endl;
  sqlite3_result_double(context, res);
}
#endif

/* SQLite invokes this routine once when it loads the extension.
** Create new functions, collating sequences, and virtual table
** modules here.  This is usually the only exported symbol in
** the shared library.
*/
extern "C" int sqlite3_extension_init(
  sqlite3 *db,
  char **pzErrMsg,
  const sqlite3_api_routines *pApi
){
  SQLITE_EXTENSION_INIT2(pApi);
  std::map<std::string,boost::any> *molMap=new std::map<std::string,boost::any>();
  sqlite3_create_function(db, "rdk_molNumAtoms", 1, SQLITE_ANY, 0, numAtomsFunc, 0, 0);
  sqlite3_create_function(db, "rdk_molAMW", 1, SQLITE_ANY, 0, molWtFunc, 0, 0);
  sqlite3_create_function(db, "rdk_smilesToBlob", 1, SQLITE_ANY, 0, smilesToBlob, 0, 0);
  sqlite3_create_function(db, "rdk_molToRDKitFP", 1, SQLITE_ANY, 0,
			  blobToRDKitFingerprint, 0, 0);
  sqlite3_create_function(db, "rdk_bvTanimotoSim", 2, SQLITE_ANY, 0,
			  bvTanimotoSim, 0, 0);
  sqlite3_create_function(db, "rdk_ucvTanimotoSim", 2, SQLITE_ANY, 0,
			  ucvTanimotoSim, 0, 0);
  sqlite3_create_function(db, "rdk_molToAtomPairFP", 1, SQLITE_ANY, 0,
			  blobToAtomPairFingerprint, 0, 0);
  sqlite3_create_function(db, "rdk_sivDiceSim", 2, SQLITE_ANY, 0,
			  sivDiceSim, 0, 0);
  sqlite3_create_function(db, "rdk_sivDiceSim2", 2, SQLITE_ANY, 0,
			  sivDiceSim2, 0, 0);
  sqlite3_create_function(db, "rdk_molHasSubstruct", 2, SQLITE_ANY,
			  static_cast<void *>(molMap),
                          molHasSubstruct, 0, 0);
  sqlite3_create_function(db, "rdk_molSubstructCount", 2, SQLITE_ANY,
			  static_cast<void *>(molMap),
                          molSubstructCount, 0, 0);
  sqlite3_create_function(db, "rdk_molLogP", 1, SQLITE_ANY, 0, molLogPFunc, 0, 0);
  return 0;
}
