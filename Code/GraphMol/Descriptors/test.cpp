// $Id$
//
//  Copyright (C) 2004-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

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
using namespace RDKit;
using namespace RDKit::Descriptors;

void test1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Crippen parameter acquisition." << std::endl;

  CrippenParamCollection *params=CrippenParamCollection::getParams();
  TEST_ASSERT(params);
  
  CrippenParams p=*(params->begin());
  TEST_ASSERT(p.label=="C1");
  TEST_ASSERT(p.smarts=="[CH4]");

  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Crippen calculation." << std::endl;

  ROMol *mol;
  double logp,mr;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,0.6361));
  TEST_ASSERT(feq(mr,6.7310));
  // check that caching works:
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,0.6361));
  TEST_ASSERT(feq(mr,6.7310));
  calcCrippenDescriptors(*mol,logp,mr,true,true);
  TEST_ASSERT(feq(logp,0.6361));
  TEST_ASSERT(feq(mr,6.7310));

  // check that things work when we don't add Hs:
  calcCrippenDescriptors(*mol,logp,mr,false,true);
  TEST_ASSERT(feq(logp,0.1441));
  TEST_ASSERT(feq(mr,2.503));
  delete mol;

  mol = SmilesToMol("C=C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr,true);
  TEST_ASSERT(feq(logp,0.8022));
  TEST_ASSERT(feq(mr,11.2540));
  delete mol;

  mol = SmilesToMol("C#C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,0.2494));
  TEST_ASSERT(feq(mr,9.8900));
  delete mol;

  mol = SmilesToMol("CO");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,-0.3915));
  TEST_ASSERT(feq(mr,8.1428));
  delete mol;

  mol = SmilesToMol("C=O");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,-0.1849));
  TEST_ASSERT(feq(mr,7.121));
  delete mol;

  mol = SmilesToMol("C#[O+]");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,0.0059));
  TEST_ASSERT(feq(mr,5.6315));
  delete mol;

  mol = SmilesToMol("C(C)(C)C");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,1.6623));
  TEST_ASSERT(feq(mr,20.512));
  delete mol;

  mol = SmilesToMol("C(C)(C)(C)O");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,0.7772));
  TEST_ASSERT(feq(mr,21.9718));
  delete mol;



  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue262(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Issue262: problems with Crippen calculation from pickles." << std::endl;

  ROMol *mol,*mol2;
  RWMol *mol3;
  std::string pkl;
  double rlogp,rmr,logp,mr;

  mol = SmilesToMol("c1ncccc1");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,rlogp,rmr);
  
  MolPickler::pickleMol(*mol,pkl);

  mol2=new ROMol(pkl);
  TEST_ASSERT(mol2);
  calcCrippenDescriptors(*mol2,logp,mr);
  TEST_ASSERT(feq(logp,rlogp));
  TEST_ASSERT(feq(mr,rmr));

  mol3=new RWMol();
  TEST_ASSERT(mol3);
  MolPickler::molFromPickle(pkl,mol3);
  
  calcCrippenDescriptors(*mol3,logp,mr);
  TEST_ASSERT(feq(logp,rlogp));
  TEST_ASSERT(feq(mr,rmr));
  
  delete mol;
  delete mol2;
  delete mol3;

  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test3(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test AMW calculation." << std::endl;

  ROMol *mol,*mol2;
  double amw;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  amw = calcAMW(*mol);
  TEST_ASSERT(feq(amw,16.043,.001));
  amw = calcAMW(*mol,true);
  TEST_ASSERT(feq(amw,12.011,.001));
  mol2 = MolOps::addHs(*mol);
  amw = calcAMW(*mol2);
  TEST_ASSERT(feq(amw,16.043,.001));
  amw = calcAMW(*mol2,true);
  TEST_ASSERT(feq(amw,12.011,.001));
  delete mol;
  delete mol2;

  mol = SmilesToMol("[CH4]");
  TEST_ASSERT(mol);
  amw = calcAMW(*mol);
  TEST_ASSERT(feq(amw,16.043,.001));
  amw = calcAMW(*mol,true);
  TEST_ASSERT(feq(amw,12.011,.001));
  delete mol;

  mol = SmilesToMol("C[2H]");
  TEST_ASSERT(mol);
  amw = calcAMW(*mol);
  TEST_ASSERT(feq(amw,17.0,.1));


  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testLabute(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Labute ASA descriptors." << std::endl;
  ROMol *mol;
  double asa;

  mol = SmilesToMol("CO");
  asa=calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa,13.5335,.0001));
  asa=calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa,13.5335,.0001));
  asa=calcLabuteASA(*mol,true,true);
  TEST_ASSERT(feq(asa,13.5335,.0001));

  delete mol;
  mol = SmilesToMol("OC(=O)c1ccncc1C(=O)O");
  asa=calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa,67.2924,.0001));
  
  delete mol;
  mol = SmilesToMol("C1CCC(c2cccnc2)NC1");
  asa=calcLabuteASA(*mol);
  TEST_ASSERT(feq(asa,73.0198,.0001));
  
  delete mol;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testTPSA(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test TPSA descriptors." << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.tpsa.csv";
  std::ifstream inf(fName.c_str());
  TEST_ASSERT(inf && !inf.bad());

  while(!inf.eof()){
    std::string inl=getLine(inf);
    boost::trim(inl);
    if(inl.size()==0 || inl[0]=='#') continue;
    std::vector<std::string> tokens;
    boost::split(tokens,inl,boost::is_any_of(","));
    if(tokens.size()!=2) continue;
    std::string smiles=tokens[0];
    double oTPSA=boost::lexical_cast<double>(tokens[1]);
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    double nTPSA = calcTPSA(*mol);
    if(!feq(nTPSA,oTPSA,.0001)){
      std::cerr<<" TPSA ERR: "<<smiles<<" "<<oTPSA<<" "<<nTPSA<<std::endl;
      std::vector<double> contribs(mol->getNumAtoms());
      getTPSAAtomContribs(*mol,contribs);
      for(unsigned int i=0;i<mol->getNumAtoms();++i){
        std::cerr<<"\t"<<i<<"\t"<<contribs[i]<<std::endl;
      }
    }
    TEST_ASSERT(feq(nTPSA,oTPSA,.0001));

    // make sure that adding Hs doesn't affect the value
    // (this was issue 1969745)
    ROMol *mol2 = MolOps::addHs(*mol);
    double hTPSA=calcTPSA(*mol2);
    TEST_ASSERT(feq(nTPSA,hTPSA,.0001));
    
    delete mol2;
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testLipinski1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Lipinski parameters." << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  int idx=-1;
  while(!suppl.atEnd()){
    ROMol *mol=0;
    ++idx;
    try{
      mol=suppl.next();
    } catch(...){
      continue;
    }
    if(!mol) continue;

    unsigned int oVal,nVal;
    std::string foo;

    mol->getProp("NUM_HACCEPTORS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumHBA(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);

    mol->getProp("NUM_HDONORS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumHBD(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);

    mol->getProp("NUM_LIPINSKIHDONORS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcLipinskiHBD(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);

    mol->getProp("NUM_LIPINSKIHACCEPTORS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcLipinskiHBA(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);

    mol->getProp("NUM_RINGS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumRings(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);


    mol->getProp("NUM_HETEROATOMS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumHeteroatoms(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);

    mol->getProp("NUM_ROTATABLEBONDS",foo);
    oVal=boost::lexical_cast<unsigned int>(foo);
    nVal = calcNumRotatableBonds(*mol);
    if(oVal!=nVal){
      std::cerr<<"  failed: "<<idx<<" "<<oVal<<" "<<nVal<<std::endl;
    }
    TEST_ASSERT(oVal==nVal);

    
    
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(){
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  testIssue262();
  test3();
  testLabute();
  testTPSA();
  testLipinski1();
#endif
}
