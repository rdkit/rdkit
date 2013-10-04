#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>


#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include "inchi.h"

using namespace RDKit;


#ifdef RDK_TEST_MULTITHREADED
namespace {
  void runblock(const std::vector<ROMol *> &mols,unsigned int count,
                unsigned int idx,
                std::vector<std::string> &inchis,
                std::vector<std::string> &keys
                ){
    for(unsigned int j=0;j<50;j++){
      for(unsigned int i=0;i<mols.size();++i){
        if(i%count != idx) continue;
        ROMol *mol = mols[i];
        ExtraInchiReturnValues tmp;
        std::string inchi=MolToInchi(*mol,tmp);
        TEST_ASSERT(inchi==inchis[i]);
        std::string key=InchiToInchiKey(inchi);
        TEST_ASSERT(key==keys[i]);
        ROMol *mol2 = InchiToMol(inchi,tmp);
        TEST_ASSERT(mol2);
        ExtraInchiReturnValues tmp2;
        std::string inchi2=MolToInchi(*mol2,tmp2);
        TEST_ASSERT(inchi==inchi2);
        delete mol2;
      }
    }
  };
}

#include <boost/thread.hpp>  
void testMultiThread(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  std::cerr<<"reading molecules"<<std::endl;
  std::vector<ROMol *> mols;
  while(!suppl.atEnd()&&mols.size()<100){
    ROMol *mol=0;
    try{
      mol=suppl.next();
    } catch(...){
      continue;
    }
    if(!mol) continue;
    mols.push_back(mol);
  }
  std::cerr<<"generating reference data"<<std::endl;
  std::vector<std::string> inchis;
  std::vector<std::string> keys;
  BOOST_FOREACH(const ROMol *mol,mols){
    ExtraInchiReturnValues tmp;
    std::string inchi=MolToInchi(*mol,tmp);
    std::string key=InchiToInchiKey(inchi);
    inchis.push_back(inchi);
    keys.push_back(key);
  }
  
  boost::thread_group tg;
  std::cerr<<"processing"<<std::endl;
  unsigned int count=4;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock,mols,count,i,inchis,keys));
  }
  tg.join_all();

  for(unsigned int i=0;i<mols.size();++i) delete mols[i];


  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread(){
  {
  }
}
#endif

void testGithubIssue3(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) <<"testing github issue 3: bad inchis when mol has no stereoinfo" << std::endl;
  {
    std::string fName= getenv("RDBASE");
    fName += "/External/INCHI-API/test_data/github3.mol";
    ROMol *m = static_cast<ROMol *>(MolFileToMol(fName));
    TEST_ASSERT(m);
    std::string smi=MolToSmiles(*m,true);
    TEST_ASSERT(smi=="CNC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO");

    ExtraInchiReturnValues tmp;
    std::string inchi=MolToInchi(*m,tmp);
    TEST_ASSERT(inchi=="InChI=1S/C7H17NO5/c1-8-2-4(10)6(12)7(13)5(11)3-9/h4-13H,2-3H2,1H3/t4-,5+,6+,7+/m0/s1");

    // blow out the stereo information with a copy:
    RWMol *m2=new RWMol(*m);
    m2->clearComputedProps();
    MolOps::sanitizeMol(*m2);

    inchi=MolToInchi(*m2,tmp);
    TEST_ASSERT(inchi=="InChI=1S/C7H17NO5/c1-8-2-4(10)6(12)7(13)5(11)3-9/h4-13H,2-3H2,1H3/t4-,5+,6+,7+/m0/s1");

    delete m;
    delete m2;
    
  }
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void testGithubIssue8(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) <<"testing a consequence of the fix for github issue 8: bad mols from some inchis" << std::endl;
  {
    std::string fName= getenv("RDBASE");
    fName += "/External/INCHI-API/test_data/github8_extra.mol";
    ROMol *m = static_cast<ROMol *>(MolFileToMol(fName,true,false));
    TEST_ASSERT(m);

    ExtraInchiReturnValues tmp;
    std::string inchi=MolToInchi(*m,tmp);
    TEST_ASSERT(inchi=="InChI=1S/C20H13IN2O3/c21-15-8-14-18(9-16(15)23)26-17-7-10(22)5-6-13(17)19(14)11-3-1-2-4-12(11)20(24)25/h1-9,23H,22H2,(H,24,25)/b23-16+/i21-2");

    ExtraInchiReturnValues tmp2;
    ROMol *m2 = InchiToMol(inchi,tmp2);
    TEST_ASSERT(m2);
    std::string smi=MolToSmiles(*m2,true);
    TEST_ASSERT(smi=="N=c1cc2oc3cc(N)ccc3c(-c3ccccc3C(=O)O)c-2cc1[125I]");

    inchi=MolToInchi(*m2,tmp2);
    TEST_ASSERT(inchi=="InChI=1S/C20H13IN2O3/c21-15-8-14-18(9-16(15)23)26-17-7-10(22)5-6-13(17)19(14)11-3-1-2-4-12(11)20(24)25/h1-9,23H,22H2,(H,24,25)/i21-2");

    delete m;
    delete m2;
    
  }
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void testGithubIssue40(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) <<"testing github issue 40: bad MWs from inchis" << std::endl;
  {
    ExtraInchiReturnValues tmp;
    std::string inchi="InChI=1S/C10H9N3O/c1-7-11-10(14)9(13-12-7)8-5-3-2-4-6-8/h2-6H,1H3,(H,11,12,14)";
    ROMol *m = InchiToMol(inchi,tmp);
    TEST_ASSERT(m);

    double mw=Descriptors::calcAMW(*m);
    TEST_ASSERT(feq(mw,187.202));
    
    delete m;
  }
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void testGithubIssue67(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) <<"testing github issue 67: seg fault from inchi" << std::endl;
  {
    ExtraInchiReturnValues tmp;
    std::string inchi="InChI=1S/C18H17N3/c19-18(20)21-14-17-12-10-16(11-13-17)9-5-4-8-15-6-2-1-3-7-15/h1-3,6-13H,14H2,(H4,19,20,21)/b9-8+";
    ROMol *m = InchiToMol(inchi,tmp);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==21);
    
    delete m;
  }
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}

void testGithubIssue68(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) <<"testing github issue 68: hang while reading InChI" << std::endl;
  {
    ExtraInchiReturnValues tmp;
    std::string inchi="InChI=1S/C11H20NSSi2.Li/c1-14(2)8-9-15(3,4)12(14)10-11-6-5-7-13-11;/h6-7H,8-10H2,1-4H3;/q-1;+1";
    BOOST_LOG(rdInfoLog) <<"  parse 1:" << std::endl;
    ROMol *m = InchiToMol(inchi,tmp);
    BOOST_LOG(rdInfoLog) <<"  done" << std::endl;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==16);
    
    delete m;
  }
  {
    ExtraInchiReturnValues tmp;
    std::string inchi="InChI=1S/C12H22NSSi2.Li/c1-15(2)10-11-16(3,4)13(15)8-7-12-6-5-9-14-12;/h6,9H,7-8,10-11H2,1-4H3;/q-1;+1";
    BOOST_LOG(rdInfoLog) <<"  parse 2:" << std::endl;
    ROMol *m = InchiToMol(inchi,tmp);
    BOOST_LOG(rdInfoLog) <<"  done" << std::endl;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms()==17);
    
    delete m;
  }
  BOOST_LOG(rdInfoLog) <<"done" << std::endl;
}


//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(){
  RDLog::InitLogs();
  testMultiThread();
  testGithubIssue3();
  testGithubIssue8();
  testGithubIssue40();
  testGithubIssue67();
  testGithubIssue68();
}
