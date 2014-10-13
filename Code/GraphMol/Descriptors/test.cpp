// $Id$
//
//  Copyright (C) 2004-2012 Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/SmilesParse/SmilesWrite.h>
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

  const CrippenParamCollection *params=CrippenParamCollection::getParams();
  TEST_ASSERT(params);
  
  CrippenParams p=*(params->begin());
  TEST_ASSERT(p.idx==0);
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

void test3a(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Exact MW calculation." << std::endl;

  ROMol *mol,*mol2;
  double mw;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw,16.031,.001));
  mw = calcExactMW(*mol,true);
  TEST_ASSERT(feq(mw,12.000,.001));
  mol2 = MolOps::addHs(*mol);
  mw = calcExactMW(*mol2);
  TEST_ASSERT(feq(mw,16.031,.001));
  mw = calcExactMW(*mol2,true);
  TEST_ASSERT(feq(mw,12.000,.001));
  delete mol;
  delete mol2;

  mol = SmilesToMol("[CH4]");
  TEST_ASSERT(mol);
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw,16.031,.001));
  mw = calcExactMW(*mol,true);
  TEST_ASSERT(feq(mw,12.000,.001));
  delete mol;

  mol = SmilesToMol("C[2H]");
  TEST_ASSERT(mol);
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw,17.037,.001));
  mw = calcExactMW(*mol,true);
  TEST_ASSERT(feq(mw,12.000,.001));
  delete mol;

  mol = SmilesToMol("Cl");
  TEST_ASSERT(mol);
  mw = calcAMW(*mol);
  TEST_ASSERT(feq(mw,35.453+1.008,.001));
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw,34.9688+1.0078,.001));
  delete mol;

  mol = SmilesToMol("[35ClH]");
  TEST_ASSERT(mol);
  mw = calcAMW(*mol);
  TEST_ASSERT(feq(mw,34.9688+1.008,.001));
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw,34.9688+1.0078,.001));
  delete mol;

  mol = SmilesToMol("[36ClH]");
  TEST_ASSERT(mol);
  mw = calcAMW(*mol);
  TEST_ASSERT(feq(mw,35.9683+1.008,.001));
  mw = calcExactMW(*mol);
  TEST_ASSERT(feq(mw,35.9683+1.0078,.001));
  delete mol;

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
    ROMol *mol = SmilesToMol(smiles);
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

void testVSADescriptors(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test VSA descriptors." << std::endl;

  {
    ROMol *mol;
    std::vector<double> vals;
  
    mol = SmilesToMol("CO");
    vals = calcSlogP_VSA(*mol);
    TEST_ASSERT(vals.size()==12);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 1:
        TEST_ASSERT(feq(vals[i],12.216,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;

    mol = SmilesToMol("CCO");
    vals = calcSlogP_VSA(*mol);
    TEST_ASSERT(vals.size()==12);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 1:
        TEST_ASSERT(feq(vals[i],11.713,.001));
        break;
      case 4:
        TEST_ASSERT(feq(vals[i],6.924,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;

    mol = SmilesToMol("Fc1ccccc1");
    vals = calcSlogP_VSA(*mol);
    TEST_ASSERT(vals.size()==12);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 3:
        TEST_ASSERT(feq(vals[i],5.817,.001));
        break;
      case 5:
        TEST_ASSERT(feq(vals[i],30.332,.001));
        break;
      case 9:
        TEST_ASSERT(feq(vals[i],4.390,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;
  }


  {
    ROMol *mol;
    std::vector<double> vals;
  
    mol = SmilesToMol("CO");
    vals = calcSMR_VSA(*mol);
    TEST_ASSERT(vals.size()==10);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 0:
        TEST_ASSERT(feq(vals[i],5.106,.001));
        break;
      case 5:
        TEST_ASSERT(feq(vals[i],7.110,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;

    mol = SmilesToMol("CCO");
    vals = calcSMR_VSA(*mol);
    TEST_ASSERT(vals.size()==10);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 0:
        TEST_ASSERT(feq(vals[i],5.106,.001));
        break;
      case 4:
        TEST_ASSERT(feq(vals[i],6.924,.001));
        break;
      case 5:
        TEST_ASSERT(feq(vals[i],6.607,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;
  }

  {
    ROMol *mol;
    std::vector<double> vals;
  
    mol = SmilesToMol("CO");
    vals = calcPEOE_VSA(*mol);
    TEST_ASSERT(vals.size()==14);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 0:
        TEST_ASSERT(feq(vals[i],5.106,.001));
        break;
      case 7:
        TEST_ASSERT(feq(vals[i],7.110,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;

    mol = SmilesToMol("CCO");
    vals = calcPEOE_VSA(*mol);
    TEST_ASSERT(vals.size()==14);
    for(unsigned int i=0;i<vals.size();++i){
      switch(i){
      case 0:
        TEST_ASSERT(feq(vals[i],5.106,.001));
        break;
      case 6:
        TEST_ASSERT(feq(vals[i],6.924,.001));
        break;
      case 7:
        TEST_ASSERT(feq(vals[i],6.607,.001));
        break;
      default:
        TEST_ASSERT(feq(vals[i],0,.001));
      }
    }
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMolFormula(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test molecular formula calculation." << std::endl;

  ROMol *mol,*mol2;
  std::string formula;

  mol = SmilesToMol("C");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH4");
  mol2 = MolOps::addHs(*mol);
  formula = calcMolFormula(*mol2);
  TEST_ASSERT(formula=="CH4");
  delete mol;
  delete mol2;

  mol = SmilesToMol("[CH4]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH4");
  delete mol;

  mol = SmilesToMol("CO");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH4O");
  mol2 = MolOps::addHs(*mol);
  formula = calcMolFormula(*mol2);
  TEST_ASSERT(formula=="CH4O");
  delete mol;
  delete mol2;

  mol = SmilesToMol("C(=O)N");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH3NO");
  mol2 = MolOps::addHs(*mol);
  formula = calcMolFormula(*mol2);
  TEST_ASSERT(formula=="CH3NO");
  delete mol;
  delete mol2;

  mol = SmilesToMol("C(=O)=O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CO2");
  delete mol;

  mol = SmilesToMol("C(=O)[O-]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CHO2-");
  delete mol;

  mol = SmilesToMol("C([O-])[O-]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH2O2-2");
  delete mol;

  mol = SmilesToMol("C([NH3+])[O-]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH5NO");
  delete mol;

  mol = SmilesToMol("C([NH3+])O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH6NO+");
  delete mol;

  mol = SmilesToMol("C([NH3+])[NH3+]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH8N2+2");
  delete mol;

  // H isotope tests
  mol = SmilesToMol("[2H]C([3H])O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="CH4O");
  formula = calcMolFormula(*mol,true);
  TEST_ASSERT(formula=="CH2DTO");
  formula = calcMolFormula(*mol,true,false);
  TEST_ASSERT(formula=="CH2[2H][3H]O");
  delete mol;

  // isotope test
  mol = SmilesToMol("[13CH3]C([2H])O");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="C2H6O");
  formula = calcMolFormula(*mol,true);
  TEST_ASSERT(formula=="C[13C]H5DO");
  formula = calcMolFormula(*mol,true,false);
  TEST_ASSERT(formula=="C[13C]H5[2H]O");
  delete mol;

  // isotope test
  mol = SmilesToMol("[13CH3]C[13CH2]C");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="C4H10");
  formula = calcMolFormula(*mol,true);
  TEST_ASSERT(formula=="C2[13C]2H10");
  formula = calcMolFormula(*mol,true,false);
  TEST_ASSERT(formula=="C2[13C]2H10");
  delete mol;

  // order test
  mol = SmilesToMol("[13CH3]C[13CH2]CB(O)O[2H]");
  TEST_ASSERT(mol);
  formula = calcMolFormula(*mol);
  TEST_ASSERT(formula=="C4H11BO2");
  formula = calcMolFormula(*mol,true);
  TEST_ASSERT(formula=="C2[13C]2H10DBO2");
  formula = calcMolFormula(*mol,true,false);
  TEST_ASSERT(formula=="C2[13C]2H10[2H]BO2");
  delete mol;


  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testIssue3415534(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Issue 3415534." << std::endl;

  {
    ROMol *mol= SmilesToMol("CN");
    TEST_ASSERT(mol);
    int nHBD=calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD==2);
    delete mol;
  }
  {
    ROMol *mol= SmilesToMol("CNC");
    TEST_ASSERT(mol);
    int nHBD=calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD==1);
    delete mol;
  }
  {
    ROMol *mol= SmilesToMol("C[NH3+]");
    TEST_ASSERT(mol);
    int nHBD=calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD==3);
    delete mol;
  }
  {
    ROMol *mol= SmilesToMol("CO");
    TEST_ASSERT(mol);
    int nHBD=calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD==1);
    delete mol;
  }
  {
    ROMol *mol= SmilesToMol("C[OH2+]");
    TEST_ASSERT(mol);
    int nHBD=calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD==2);
    delete mol;
  }
  {
    ROMol *mol= SmilesToMol("COC");
    TEST_ASSERT(mol);
    int nHBD=calcLipinskiHBD(*mol);
    TEST_ASSERT(nHBD==0);
    delete mol;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue3433771(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Issue3433771: Bad definition for Crippen atom type O11." << std::endl;

  ROMol *mol;
  double logp,mr;

  mol = SmilesToMol("O=C(NC)n1cccc1");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,0.6756,.001));
  
  delete mol;
  mol = SmilesToMol("O=C(n1cccc1)n1cccc1");
  TEST_ASSERT(mol);
  calcCrippenDescriptors(*mol,logp,mr);
  TEST_ASSERT(feq(logp,1.806,.001));
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
  void runblock(const std::vector<ROMol *> &mols,unsigned int count,unsigned int idx){
    for(unsigned int j=0;j<1000;j++){
      for(unsigned int i=0;i<mols.size();++i){
        if(i%count != idx) continue;
        ROMol *mol = mols[i];
        int nHBD=calcNumHBD(*mol);
        int nHBA=calcNumHBA(*mol);

        unsigned int oVal;
        std::string foo;
        mol->getProp("NUM_HACCEPTORS",foo);
        oVal=boost::lexical_cast<unsigned int>(foo);
        TEST_ASSERT(oVal==nHBA);
        mol->getProp("NUM_HDONORS",foo);
        oVal=boost::lexical_cast<unsigned int>(foo);
        TEST_ASSERT(oVal==nHBD);

        int nAmide=calcNumAmideBonds(*mol);
        double logp,mr;
        calcCrippenDescriptors(*mol,logp,mr);
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
  boost::thread_group tg;

  std::cerr<<"processing"<<std::endl;
  unsigned int count=4;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock,mols,count,i));
  }
  tg.join_all();

  for(unsigned int i=0;i<mols.size();++i) delete mols[i];


  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread(){
}
#endif
void testCrippenContribs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Crippen atom type calculations." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("n1ccccc1CO");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());
    std::vector<unsigned int> ts(mol->getNumAtoms());
    std::vector<std::string> ls(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true,&ts,&ls);
    TEST_ASSERT(ts[0]==59);
    TEST_ASSERT(ts[1]==25);
    TEST_ASSERT(ts[2]==25);
    TEST_ASSERT(ts[3]==25);
    TEST_ASSERT(ts[4]==25);
    TEST_ASSERT(ts[5]==28);
    TEST_ASSERT(ts[6]==17);
    TEST_ASSERT(ts[7]==69);

    TEST_ASSERT(ls[0]=="N11");
    TEST_ASSERT(ls[7]=="O2");
    delete mol;
  }  

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testIssue252(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Issue252: Bad definitions for Crippen atom types." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("O=[N+]([O-])C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],0.0335,.001));
    TEST_ASSERT(feq(logp[1],-0.3396,.001));
    TEST_ASSERT(feq(logp[2],0.0335,.001));
    TEST_ASSERT(feq(logp[3],-0.2035,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("CP");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2035,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(C)P");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2035,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(C)(C)P");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2051,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(C)(C)(C)P");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2051,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(=C)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],0.264,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(=C)(C)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],0.264,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C=C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],0.1551,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("O=S");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.3339,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("S=O");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.0024,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2035,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2051,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)(Cl)C");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.2051,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.0516,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],0.1193,.001));
    delete mol;
  }  
  {
    ROMol *mol;
    mol = SmilesToMol("C(Cl)(Cl)(Cl)c1ccccc1");
    TEST_ASSERT(mol);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());

    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[0],-0.0967,.001));
    delete mol;
  }  
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testChiVs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of ChiVs." << std::endl;

  {
    std::string sdata[] = {"CCCCCC",
                           "CCC(C)CC",
                           "CC(C)CCC",
                           "CC(C)C(C)C",
                           "CC(C)(C)CC",
                           "CCCCCO",
                           "CCC(O)CC",
                           "CC(O)(C)CC",
                           "c1ccccc1O",
                           "CCCl",
                           "CCBr",
                           "CCI",
                           "EOS"};
    double ddata[] = {4.828,
                      4.992,
                      4.992,
                      5.155,
                      5.207,
                      4.276,
                      4.439,
                      4.654,
                      3.834,
                      2.841,
                      3.671,
                      4.242};
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi0v(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {2.914,2.808,2.770,
                      2.643,2.561,
                      2.523,2.489,2.284,2.134};
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi1v(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {1.707,1.922,2.183,
                      2.488,2.914,
                      1.431,1.470,2.166,1.336};
    
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi2v(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {0.957,1.394,0.866,1.333,1.061,
                      0.762,0.943,0.865,0.756};
    
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi3v(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {0.500,0.289,0.577,
                      0.000,0.000,0.362,0.289,0.000,0.428};

    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi4v(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testChiNs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of ChiNs." << std::endl;

  {
    std::string sdata[] = {"CCCCCC",
                           "CCC(C)CC",
                           "CC(C)CCC",
                           "CC(C)C(C)C",
                           "CC(C)(C)CC",
                           "CCCCCO",
                           "CCC(O)CC",
                           "CC(O)(C)CC",
                           "c1ccccc1O",
                           "CCCl",
                           "CCBr",
                           "CCI",
                           "EOS"};
    double ddata[] = {4.828,
                      4.992,
                      4.992,
                      5.155,
                      5.207,
                      4.276,
                      4.439,
                      4.654,
                      3.834,
                      2.085,
                      2.085,
                      2.085};
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi0n(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "C=S",
                          "EOS"};

    double ddata[] = {2.914,2.808,2.770,
                      2.643,2.561,
                      2.523,2.489,2.284,2.134,
                      0.289
    };
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi1n(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {1.707,1.922,2.183,
                      2.488,2.914,
                      1.431,1.470,2.166,1.336};
    
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi2n(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {0.957,1.394,0.866,1.333,1.061,
                      0.762,0.943,0.865,0.756};
    
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi3n(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  
  {
    std::string sdata[] ={"CCCCCC","CCC(C)CC","CC(C)CCC",
                          "CC(C)C(C)C","CC(C)(C)CC",
                          "CCCCCO","CCC(O)CC","CC(O)(C)CC","c1ccccc1O",
                          "EOS"};

    double ddata[] = {0.500,0.289,0.577,
                      0.000,0.000,0.362,0.289,0.000,0.428};

    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcChi4n(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testHallKierAlpha(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of HallKierAlpha." << std::endl;

  {
    std::string sdata[] = {
      "C=O",
      "CCC1(CC)C(=O)NC(=O)N(C)C1=O",
      "OCC(O)C(O)C(O)C(O)CO",
      "OCC1OC(O)C(O)C(O)C1O",
      "Fc1c[nH]c(=O)[nH]c1=O",
      "OC1CNC(C(=O)O)C1",
      "CCCc1[nH]c(=S)[nH]c(=O)c1",
      "CN(CCCl)CCCl",
      "CBr",
      "CI",
      "EOS"
    };
    double ddata[] = {
      -0.3300,
      -1.3900,
      -0.2400,
      -0.2400,
      -1.3900,
      -0.6100,
      -0.9000,
      0.5400,
      0.480,
      0.730,
    };
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcHallKierAlpha(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testKappa1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of Kappa1." << std::endl;

  {
    std::string sdata[] = {
      "C12CC2C3CC13",
      "C1CCC12CC2",
      "C1CCCCC1",
      "CCCCCC",
      "CCC(C)C1CCC(C)CC1",
      "CC(C)CC1CCC(C)CC1",
      "CC(C)C1CCC(C)CCC1",
      "EOS"
    };
    double ddata[] = {
      2.344,
      3.061,
      4.167,
      6.000,
      9.091,
      9.091,
      9.091
    };
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcKappa1(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testKappa2(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of Kappa2." << std::endl;

  {
    std::string sdata[] = {
      "[C+2](C)(C)(C)(C)(C)C",
      "[C+](C)(C)(C)(C)(CC)",
      "C(C)(C)(C)(CCC)",
      "CC(C)CCCC",
      "CCCCCCC",
      "CCCCCC",
      "CCCCCCC",
      "C1CCCC1",
      "C1CCCC1C",
      "C1CCCCC1",
      "C1CCCCCC1",
      "CCCCC",
      "CC=CCCC",
      "C1=CN=CN1",
      "c1ccccc1",
      "c1cnccc1",
      "n1ccncc1",
      "CCCCF",
      "CCCCCl",
      "CCCCBr",
      "CCC(C)C1CCC(C)CC1",
      "CC(C)CC1CCC(C)CC1",
      "CC(C)C1CCC(C)CCC1",
      "EOS"
    };
    double ddata[] = {
      0.667000,
      1.240000,
      2.344400,
      4.167000,
      6.000000,
      5.000000,
      6.000000,
      1.440000,
      1.633000,
      2.222000,
      3.061000,
      4.000000,
      4.740000,
      0.884000,
      1.606000,
      1.552000,
      1.500000,
      3.930000,
      4.290000,
      4.480000,
      4.133000,
      4.133000,
      4.133000
    };
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcKappa2(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
void testKappa3(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test calculation of Kappa3." << std::endl;

  {
    std::string sdata[] = {
      "C[C+](C)(C)(C)C(C)(C)C",
      "CCC(C)C(C)(C)(CC)",
      "CCC(C)CC(C)CC",
      "CC(C)CCC(C)CC",
      "CC(C)CCCC(C)C",
      "CCC(C)C1CCC(C)CC1",
      "CC(C)CC1CCC(C)CC1",
      "CC(C)C1CCC(C)CCC1",
      "EOS"
    };
    double ddata[] = {
      2.000000,
      2.380000,
      4.500000,
      5.878000,
      8.000000,
      2.500000,
      3.265000,
      2.844000
    };
    unsigned int idx=0;
    while(sdata[idx]!="EOS"){
      ROMol *mol;
      mol = SmilesToMol(sdata[idx]);
      TEST_ASSERT(mol);
      double v=calcKappa3(*mol);
      TEST_ASSERT(feq(v,ddata[idx],0.002));
      ++idx;
      delete mol;
    }
  }  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRingDescriptors(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test ring descriptors" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/Descriptors/test_data/aid466.trunc.sdf";
  SDMolSupplier suppl(fName);
  while(!suppl.atEnd()){
    ROMol *mol=suppl.next();
    TEST_ASSERT(mol);
    unsigned int iv;
    mol->getProp("NumRings",iv);
    TEST_ASSERT(iv==calcNumRings(*mol));
    mol->getProp("NumAromaticRings",iv);
    TEST_ASSERT(iv==calcNumAromaticRings(*mol));
    mol->getProp("NumSaturatedRings",iv);
    TEST_ASSERT(iv==calcNumSaturatedRings(*mol));
    mol->getProp("NumAromaticHeterocycles",iv);
    TEST_ASSERT(iv==calcNumAromaticHeterocycles(*mol));
    mol->getProp("NumAromaticCarbocycles",iv);
    TEST_ASSERT(iv==calcNumAromaticCarbocycles(*mol));
    mol->getProp("NumSaturatedHeterocycles",iv);
    TEST_ASSERT(iv==calcNumSaturatedHeterocycles(*mol));
    mol->getProp("NumSaturatedCarbocycles",iv);
    TEST_ASSERT(iv==calcNumSaturatedCarbocycles(*mol));
    mol->getProp("NumAliphaticRings",iv);
    TEST_ASSERT(iv==calcNumAliphaticRings(*mol));
    mol->getProp("NumAliphaticHeterocycles",iv);
    TEST_ASSERT(iv==calcNumAliphaticHeterocycles(*mol));
    mol->getProp("NumAliphaticCarbocycles",iv);
    TEST_ASSERT(iv==calcNumAliphaticCarbocycles(*mol));
    mol->getProp("NumHeterocycles",iv);
    TEST_ASSERT(iv==calcNumHeterocycles(*mol));
    
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMiscCountDescriptors(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test other count descriptors." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("OCCO");
    TEST_ASSERT(feq(calcFractionCSP3(*mol),1.0,0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("OO");
    TEST_ASSERT(feq(calcFractionCSP3(*mol),0.0,0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("OC=CO");
    TEST_ASSERT(feq(calcFractionCSP3(*mol),0.0,0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("CCC=C");
    TEST_ASSERT(feq(calcFractionCSP3(*mol),0.5,0.001));
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMQNs(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MQN" << std::endl;

  unsigned int tgt[42]={98,  0,  4,  0,  0,  1,  0,  3,  9,  5,  4,  0, 29,  3,  0, 66, 35,
                      0, 25, 30, 21,  2,  2,  0,  0,  6, 12,  6,  0, 70, 26,  0,  0,  0,
                      2, 16,  0,  0,  0,  0, 10,  5};

  std::vector<unsigned int> accum(42,0);
  
  std::string fName = getenv("RDBASE");
  fName += "/Code/GraphMol/Descriptors/test_data/aid466.trunc.sdf";
  SDMolSupplier suppl(fName);
  while(!suppl.atEnd()){
    ROMol *mol=suppl.next();
    TEST_ASSERT(mol);
    std::vector<unsigned int> v = calcMQNs(*mol);
    TEST_ASSERT(v.size()==42);
    for(unsigned int i=0;i<42;++i) accum[i]+=v[i];
    delete mol;
  }
  for(unsigned int i=0;i<42;++i){
    if(accum[i] != tgt[i]){
      std::cerr<<" !! "<<i<<accum[i]<<"!="<<tgt[i]<<std::endl;
    }
    TEST_ASSERT(accum[i]==tgt[i]);
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue56(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test GitHub Issue 56." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("[H+]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw,1.008,0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,1.0078,0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("[2H+]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw,2.014,0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,2.014,0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("[H]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw,1.008,0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,1.0078,0.001));
    delete mol;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("[2H]");
    double mw;
    mw = calcAMW(*mol);
    TEST_ASSERT(feq(mw,2.014,0.001));
    mw = calcExactMW(*mol);
    TEST_ASSERT(feq(mw,2.014,0.001));
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue92(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Github92: Bad Crippen atom type for pyrrole H." << std::endl;

  {
    RWMol *mol;
    mol = SmilesToMol("c1cccn1[H]",0,0);
    TEST_ASSERT(mol);
    MolOps::sanitizeMol(*mol);
    TEST_ASSERT(mol->getNumAtoms()==6);
    std::vector<double> logp(mol->getNumAtoms());
    std::vector<double> mr(mol->getNumAtoms());
    getCrippenAtomContribs(*mol,logp,mr,true);
    TEST_ASSERT(feq(logp[5],0.2142,.001));
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
  test3a();
  testLabute();
  testTPSA();
  testLipinski1();
  testVSADescriptors();
  testMolFormula();
  testIssue3415534();
  testIssue3433771();
  testMultiThread();
  testIssue252();
  testChiVs();
  testChiNs();
  testHallKierAlpha();
  testKappa1();
  testKappa2();
  testKappa3();
  testCrippenContribs();
  testRingDescriptors();
  testMiscCountDescriptors();
  testMQNs();
#endif
  testGitHubIssue56();
  testGitHubIssue92();

}
