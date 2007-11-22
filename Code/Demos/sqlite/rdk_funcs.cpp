#include <sqlite3ext.h>
SQLITE_EXTENSION_INIT1
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <string>
#include <map>

/*
** The half() SQL function returns half of its input value.
*/
static void halfFunc(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  sqlite3_result_double(context, 0.5*sqlite3_value_double(argv[0]));
}

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
  int res=RDKit::SubstructMatch(*m,*patt,matches,true,false);
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
  std::string text=fp->ToString();
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
  sqlite3_create_function(db, "half", 1, SQLITE_ANY, 0, halfFunc, 0, 0);
  std::map<std::string,boost::any> *molMap=new std::map<std::string,boost::any>();
  sqlite3_create_function(db, "rdk_molNumAtoms", 1, SQLITE_ANY, 0, numAtomsFunc, 0, 0);
  sqlite3_create_function(db, "rdk_molAMW", 1, SQLITE_ANY, 0, molWtFunc, 0, 0);
  sqlite3_create_function(db, "rdk_smilesToBlob", 1, SQLITE_ANY, 0, smilesToBlob, 0, 0);
  sqlite3_create_function(db, "rdk_molToRDKitFP", 1, SQLITE_ANY, 0,
			  blobToRDKitFingerprint, 0, 0);
  sqlite3_create_function(db, "rdk_bvTanimotoSim", 2, SQLITE_ANY, 0,
			  bvTanimotoSim, 0, 0);
  sqlite3_create_function(db, "rdk_molHasSubstruct", 2, SQLITE_ANY,
			  static_cast<void *>(molMap),
                          molHasSubstruct, 0, 0);
  sqlite3_create_function(db, "rdk_molSubstructCount", 2, SQLITE_ANY,
			  static_cast<void *>(molMap),
                          molSubstructCount, 0, 0);
  sqlite3_create_function(db, "rdk_molLogP", 1, SQLITE_ANY, 0, molLogPFunc, 0, 0);
  return 0;
}
