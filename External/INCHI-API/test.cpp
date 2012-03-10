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
}
#endif

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(){
  RDLog::InitLogs();
  testMultiThread();
}
