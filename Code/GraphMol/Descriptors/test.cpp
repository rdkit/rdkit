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

  const CrippenParamCollection *params=CrippenParamCollection::getParams();
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
#endif
  testIssue3433771();
  testMultiThread();
  testIssue252();
}
