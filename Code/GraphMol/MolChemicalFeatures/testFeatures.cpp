// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <fstream>
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <sstream>
#include <boost/shared_ptr.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureDef.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <GraphMol/MolChemicalFeatures/FeatureParser.h>

using namespace RDKit;

void test1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "AtomType parser testing." << std::endl;

  std::string inLine;
  std::map<std::string,std::string> atomTypeDefs;
  atomTypeDefs.clear();
  bool ok;

  inLine="AtomType donor [N,O]";
  Local::parseAtomType(inLine,atomTypeDefs,0);
  TEST_ASSERT(atomTypeDefs.count("{donor}"));
  TEST_ASSERT(!atomTypeDefs.count("{unsaturatedDonor}"));
  TEST_ASSERT(atomTypeDefs["{donor}"]=="$([N,O])");
  
  // check robustness to whitespace as we expand:
  inLine=" AtomType donor\t[S]";
  Local::parseAtomType(inLine,atomTypeDefs,0);
  TEST_ASSERT(atomTypeDefs.count("{donor}"));
  TEST_ASSERT(!atomTypeDefs.count("{unsaturatedDonor}"));
  TEST_ASSERT(atomTypeDefs["{donor}"]=="$([N,O,$([S])])");

  inLine="AtomType !donor [NH0]";
  Local::parseAtomType(inLine,atomTypeDefs,0);
  TEST_ASSERT(atomTypeDefs.count("{donor}"));
  TEST_ASSERT(!atomTypeDefs.count("{unsaturatedDonor}"));
  TEST_ASSERT(atomTypeDefs["{donor}"]=="$([!$([NH0]);N,O,$([S])])");

  inLine="AtomType unsaturatedDonor [$([{donor}]!-[*])]";
  Local::parseAtomType(inLine,atomTypeDefs,0);
  TEST_ASSERT(atomTypeDefs.count("{donor}"));
  TEST_ASSERT(atomTypeDefs.count("{unsaturatedDonor}"));
  TEST_ASSERT(atomTypeDefs["{donor}"]=="$([!$([NH0]);N,O,$([S])])");
  TEST_ASSERT(atomTypeDefs["{unsaturatedDonor}"]=="$([$([$([!$([NH0]);N,O,$([S])])]!-[*])])");

  ok=false;
  try {
    inLine="";
    Local::parseAtomType(inLine,atomTypeDefs,0);
  } catch (FeatureFileParseException &fpe) {
    ok=true;
    TEST_ASSERT(fpe.lineNo()==0);
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");
  
  ok=false;
  try {
    inLine="AtomTyped donor [N,O]";
    Local::parseAtomType(inLine,atomTypeDefs,0);
  } catch (FeatureFileParseException &fpe) {
    ok=true;
    TEST_ASSERT(fpe.lineNo()==0);
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");
  

  ok=false;
  try {
    inLine="AtomType";
    Local::parseAtomType(inLine,atomTypeDefs,0);
  } catch (FeatureFileParseException &fpe) {
    ok=true;
    TEST_ASSERT(fpe.lineNo()==0);
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");
  

  ok=false;
  try {
    inLine="AtomType donor";
    Local::parseAtomType(inLine,atomTypeDefs,0);
  } catch (FeatureFileParseException &fpe) {
    ok=true;
    TEST_ASSERT(fpe.lineNo()==0);
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");
  

  ok=false;
  try {
    inLine="AtomType donor [N,O";
    Local::parseAtomType(inLine,atomTypeDefs,0);
  } catch (FeatureFileParseException &fpe) {
    ok=true;
    TEST_ASSERT(fpe.lineNo()==0);
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");
  
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "FeatureDefinition parser testing." << std::endl;

  std::string inText,inLine;
  std::istringstream ss;
  MolChemicalFeatureDef *featDef;
  bool ok;

  std::map<std::string,std::string> atomTypeDefs;
  atomTypeDefs["{donor}"]="$([O,N])";  

  inText=
    "DefineFeature HDonor1 [O,N]\n"
    "  Family HDONOR\n"
    "  Weights 1.0\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  unsigned int tmpLine=2;
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonor1");
  TEST_ASSERT(featDef->getSmarts()=="[O,N]");
  TEST_ASSERT(featDef->getNumWeights()==1);
  TEST_ASSERT(feq(*featDef->beginWeights(),1.0));
  delete featDef;
  
  inText=
    "DefineFeature HDonor1 [{donor}]\n"
    "  Family HDONOR\n"
    "  Weights 2.0\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonor1");
  TEST_ASSERT(featDef->getSmarts()=="[$([O,N])]");
  TEST_ASSERT(featDef->getNumWeights()==1);
  TEST_ASSERT(feq(*featDef->beginWeights(),1.0));
  delete featDef;


  inText=
    "DefineFeature HDonorPair [{donor}][{donor}]\n"
    "  Family HDONOR\n"
    "  Weights 1.0,1\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonorPair");
  TEST_ASSERT(featDef->getSmarts()=="[$([O,N])][$([O,N])]");
  TEST_ASSERT(featDef->getNumWeights()==2);
  TEST_ASSERT(feq(*featDef->beginWeights(),0.5));
  TEST_ASSERT(feq(*++(featDef->beginWeights()),0.5));
  delete featDef;


  ok=false;
  inText=
    "DefineFeature HDonor1 [{donor}]\n"
    "  Family HDONOR\n"
    "  Weights 1.0\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  try{
    featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  } catch (FeatureFileParseException &){
    ok=true;
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");
  
  ok=false;
  inText=
    "DefineFeature HDonor1 [{donor}]\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  try{
    featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  } catch (FeatureFileParseException &){
    ok=true;
  }
  CHECK_INVARIANT(ok,"expected parse failure did not happen");

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void test3(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test pulling feat defs from a string." << std::endl;

  //std::string pathName=getenv("RDBASE");
  //pathName += "/Code/GraphMol/MolChemicalFeatures/test_data/smallsample.feats";
  std::string inText,inLine;
  std::istringstream ss;
  int res;
  MolChemicalFeatureDef::CollectionType featureDefs;


  inText=
    "AtomType donor [N,O]\n"
    "AtomType !donor [H0]\n"
    "DefineFeature HDonor1 [{donor}]\n"
    "  Family HBondDonor\n"
    "  Weights 1.0\n"
    "EndFeature\n";
  res=parseFeatureData(inText,featureDefs);
  TEST_ASSERT(!res);
  TEST_ASSERT(featureDefs.size()==1);
  TEST_ASSERT((*featureDefs.begin())->getSmarts()=="[$([!$([H0]);N,O])]");

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test4(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "FeatureDef functionality testing." << std::endl;

  ROMol *testMol;
  MatchVectType mv;
  std::string inText;
  int res;
  MolChemicalFeatureDef::CollectionType featureDefs;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;


  inText=
    "DefineFeature HDonor1 [N,O;!H0]\n"
    "  Family HBondDonor\n"
    "  Weights 1.0\n"
    "EndFeature\n"
    "DefineFeature HAcceptor1 [N,O]\n"
    "  Family HAcceptor\n"
    "  Weights 1.0\n"
    "EndFeature\n";
  res=parseFeatureData(inText,featureDefs);
  TEST_ASSERT(!res);
  TEST_ASSERT(featureDefs.size()==2);
  featDefIt=featureDefs.begin();
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getSmarts()=="[N,O;!H0]");
  TEST_ASSERT(featDef->getPattern());
  TEST_ASSERT(featDef->getPattern()->getNumAtoms()==1);
  TEST_ASSERT(featDef->getNumWeights()==1);
  
  featDef=*(++featureDefs.begin());
  TEST_ASSERT(featDef->getSmarts()=="[N,O]");
  TEST_ASSERT(featDef->getPattern());
  TEST_ASSERT(featDef->getPattern()->getNumAtoms()==1);
  TEST_ASSERT(featDef->getNumWeights()==1);

  testMol=SmilesToMol("COCN");
  TEST_ASSERT(testMol);
  featDef=*featureDefs.begin();
  TEST_ASSERT(SubstructMatch(*testMol,*featDef->getPattern(),mv));
  TEST_ASSERT(mv.size()==1);
  TEST_ASSERT(mv[0].first==0);
  TEST_ASSERT(mv[0].second==3);
  featDef=*(++featureDefs.begin());
  TEST_ASSERT(SubstructMatch(*testMol,*featDef->getPattern(),mv));
  TEST_ASSERT(mv.size()==1);
  TEST_ASSERT(mv[0].first==0);
  TEST_ASSERT(mv[0].second==1||mv[0].second==3);
  delete testMol;
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void test5(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "FeatureFactory testing." << std::endl;

  ROMol *testMol;
  MatchVectType mv;
  std::string inText;
  FeatSPtrList featSPtrs;
  boost::shared_ptr<MolChemicalFeature> featSPtr;


  MolChemicalFeatureFactory *factory;

  MolChemicalFeatureDef::CollectionType::value_type featDef;

  inText=
    "DefineFeature HDonor1 [N,O;!H0]\n"
    "  Family HBondDonor\n"
    "  Weights 1.0\n"
    "EndFeature\n"
    "DefineFeature HAcceptor1 [N,O]\n"
    "  Family HBondAcceptor\n"
    "  Weights 1.0\n"
    "EndFeature\n";

  factory = buildFeatureFactory(inText);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);

  testMol=SmilesToMol("COCN");
  TEST_ASSERT(testMol);
  featDef=*factory->beginFeatureDefs();
  TEST_ASSERT(SubstructMatch(*testMol,*featDef->getPattern(),mv));
  BOOST_LOG(rdErrorLog) << "1" << std::endl;
  TEST_ASSERT(mv.size()==1);
  TEST_ASSERT(mv[0].first==0);
  TEST_ASSERT(mv[0].second==3);
  featDef=*(++factory->beginFeatureDefs());
  TEST_ASSERT(SubstructMatch(*testMol,*featDef->getPattern(),mv));
BOOST_LOG(rdErrorLog) << "2" << std::endl;
  TEST_ASSERT(mv.size()==1);
  TEST_ASSERT(mv[0].first==0);
  TEST_ASSERT(mv[0].second==1||mv[0].second==3);

  // Test using the factory to find features:
  featSPtrs = factory->getFeaturesForMol(*testMol);
BOOST_LOG(rdErrorLog) << "3" << std::endl;
  TEST_ASSERT(featSPtrs.size()==3);
  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="HBondDonor");
  TEST_ASSERT(featSPtr->getType()=="HDonor1");
  featSPtr=*(++featSPtrs.begin());
  TEST_ASSERT(featSPtr->getFamily()=="HBondAcceptor");
  TEST_ASSERT(featSPtr->getType()=="HAcceptor1");
  featSPtr=*(++featSPtrs.begin());
  TEST_ASSERT(featSPtr->getFamily()=="HBondAcceptor");
  TEST_ASSERT(featSPtr->getType()=="HAcceptor1");

  // Test limiting stuff with includeOnly
  featSPtrs = factory->getFeaturesForMol(*testMol,"HBondAcceptor");
BOOST_LOG(rdErrorLog) << "4" << std::endl;
  TEST_ASSERT(featSPtrs.size()==2);
  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="HBondAcceptor");
  TEST_ASSERT(featSPtr->getType()=="HAcceptor1");
  featSPtr=*(++featSPtrs.begin());
  TEST_ASSERT(featSPtr->getFamily()=="HBondAcceptor");
  TEST_ASSERT(featSPtr->getType()=="HAcceptor1");

  featSPtrs = factory->getFeaturesForMol(*testMol,"HBondDonor");
BOOST_LOG(rdErrorLog) << "5" << std::endl;
  TEST_ASSERT(featSPtrs.size()==1);
  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="HBondDonor");
  TEST_ASSERT(featSPtr->getType()=="HDonor1");

  featSPtrs = factory->getFeaturesForMol(*testMol,"NotPresent");
BOOST_LOG(rdErrorLog) << "6" << std::endl;
  TEST_ASSERT(featSPtrs.size()==0);

  delete testMol;
  delete factory;
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void test6(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Feature Location testing." << std::endl;

  ROMol *testMol;
  Conformer *conf;
  std::string inText;
  FeatSPtrList featSPtrs;
  boost::shared_ptr<MolChemicalFeature> featSPtr;


  MolChemicalFeatureFactory *factory;

  MolChemicalFeatureDef::CollectionType::value_type featDef;

  inText=
    "DefineFeature HDonor1 [N,O;!H0]\n"
    "  Family HBondDonor\n"
    "  Weights 1.0\n"
    "EndFeature\n"
    "DefineFeature Carboxyl1 C(=O)[O;H1,-]\n"
    "  Family ZnBinder\n"
    "  Weights 1.0,1.0,1.0\n"
    "EndFeature\n";

  factory = buildFeatureFactory(inText);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);

  testMol=SmilesToMol("C(=O)O");
  TEST_ASSERT(testMol);
  conf = new Conformer(3);
  testMol->addConformer(conf);
  conf->setAtomPos(0,RDGeom::Point3D(0,0,0.0));
  conf->setAtomPos(1,RDGeom::Point3D(1.2,0,0.0));
  conf->setAtomPos(2,RDGeom::Point3D(0,1.5,0.0));


  featSPtrs = factory->getFeaturesForMol(*testMol);
  TEST_ASSERT(featSPtrs.size()==2);
  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="HBondDonor");
  TEST_ASSERT(featSPtr->getType()=="HDonor1");
  TEST_ASSERT(feq(featSPtr->getPos().x,0.0));
  TEST_ASSERT(feq(featSPtr->getPos().y,1.5));
  TEST_ASSERT(feq(featSPtr->getPos().z,0.0));
  

  featSPtr=*(++featSPtrs.begin());
  TEST_ASSERT(featSPtr->getFamily()=="ZnBinder");
  TEST_ASSERT(featSPtr->getType()=="Carboxyl1");
  TEST_ASSERT(feq(featSPtr->getPos().x,0.4));
  TEST_ASSERT(feq(featSPtr->getPos().y,0.5));
  TEST_ASSERT(feq(featSPtr->getPos().z,0.0));
  delete testMol;
  
  delete factory;
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test7(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test building a FeatureFactory from a file" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/MolChemicalFeatures/test_data/featDef.txt";

  MolChemicalFeatureFactory *factory;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;

  std::ifstream inStream(fName.c_str());
  TEST_ASSERT(inStream.is_open());
  
  std::istream &instrm = static_cast<std::istream &>(inStream);
  factory = buildFeatureFactory(instrm);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);
  featDefIt=factory->beginFeatureDefs();
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getFamily()=="HBondDonor");

  featDefIt++;
  featDef=*featDefIt;
  TEST_ASSERT(featDefIt!=factory->endFeatureDefs());
  TEST_ASSERT(featDef->getFamily()=="HBondAcceptor");

  featDefIt++;
  TEST_ASSERT(featDefIt==factory->endFeatureDefs());

  delete(factory);
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testIssue224(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue 224." << std::endl;

  std::string inText,inLine;
  std::istringstream ss;
  MolChemicalFeatureDef *featDef;


  std::map<std::string,std::string> atomTypeDefs;
  atomTypeDefs["{donor}"]="$([O,N])";  

  inText=
    "\n"
    "DefineFeature HDonor1 [O,N]\n"
    "  Family HDONOR\n"
    "  Weights 1.0\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  unsigned int tmpLine=2;
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonor1");
  TEST_ASSERT(featDef->getSmarts()=="[O,N]");
  TEST_ASSERT(featDef->getNumWeights()==1);
  TEST_ASSERT(feq(*featDef->beginWeights(),1.0));
  delete featDef;

  inText=
    "\n"
    "# comment, ignore me\n"
    "DefineFeature HDonor1 [O,N]\n"
    " \n"
    "  Family HDONOR\n"
    " \n"
    "# comment, ignore me\n"
    "  Weights 1.0\n"
    "# comment, ignore me\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonor1");
  TEST_ASSERT(featDef->getSmarts()=="[O,N]");
  TEST_ASSERT(featDef->getNumWeights()==1);
  TEST_ASSERT(feq(*featDef->beginWeights(),1.0));
  delete featDef;


  // repeat those tests using some continuation lines
  // (Issue 346):
  inText=
    "\n"
    "# comment, ignore me\n"
    "DefineFeature HDonor1 \\\n"
    "  [O,N]\n"
    "  Family HDONOR\n"
    " \n"
    "# comment, ignore me\n"
    "  Weights 1.0\n"
    "# comment, ignore me\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonor1");
  TEST_ASSERT(featDef->getSmarts()=="[O,N]");
  TEST_ASSERT(featDef->getNumWeights()==1);
  TEST_ASSERT(feq(*featDef->beginWeights(),1.0));
  delete featDef;

  inText=
    "\n"
    "# comment, ignore me\n"
    "DefineFeature HDonor1 \\\n"
    "# comment, ignore me\n"
    "  [O,N]\n"
    "  Family \\\n"
    "\tHDONOR\n"
    "# comment, ignore me\n"
    "  Weights \\\n"
    "1.0\n"
    "# comment, ignore me\n"
    "EndFeature\n";
  ss.clear();
  ss.str(inText);
  inLine = RDKit::getLine(ss);
  featDef = Local::parseFeatureDef(ss,inLine,tmpLine,atomTypeDefs);
  TEST_ASSERT(featDef);
  TEST_ASSERT(featDef->getFamily()=="HDONOR");
  TEST_ASSERT(featDef->getType()=="HDonor1");
  TEST_ASSERT(featDef->getSmarts()=="[O,N]");
  TEST_ASSERT(featDef->getNumWeights()==1);
  TEST_ASSERT(feq(*featDef->beginWeights(),1.0));
  delete featDef;

  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue225(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue 225." << std::endl;


  MatchVectType mv;
  std::string inText;
  int res;
  MolChemicalFeatureDef::CollectionType featureDefs;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;


  inText=
    "DefineFeature HDonor1 [N,O;!H0]\n"
    "  Family HBondDonor\n"
    "  Weights 1.0\n"
    "EndFeature\n"
    "\n"
    "DefineFeature HAcceptor1 [N,O]\n"
    "  Family HAcceptor\n"
    "  Weights 1.0\n"
    "EndFeature\n";
  res=parseFeatureData(inText,featureDefs);
  TEST_ASSERT(!res);
  TEST_ASSERT(featureDefs.size()==2);
  featDefIt=featureDefs.begin();
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getSmarts()=="[N,O;!H0]");
  TEST_ASSERT(featDef->getPattern());
  TEST_ASSERT(featDef->getPattern()->getNumAtoms()==1);
  TEST_ASSERT(featDef->getNumWeights()==1);
  
  featDef=*(++featureDefs.begin());
  TEST_ASSERT(featDef->getSmarts()=="[N,O]");
  TEST_ASSERT(featDef->getPattern());
  TEST_ASSERT(featDef->getPattern()->getNumAtoms()==1);
  TEST_ASSERT(featDef->getNumWeights()==1);

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testIssue346(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test Issue346" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/MolChemicalFeatures/test_data/featDef2.txt";

  MolChemicalFeatureFactory *factory;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;

  std::ifstream inStream(fName.c_str());
  TEST_ASSERT(inStream.is_open());
  
  std::istream &instrm = static_cast<std::istream &>(inStream);
  factory = buildFeatureFactory(instrm);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);
  featDefIt=factory->beginFeatureDefs();
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getFamily()=="HBondDonor");

  featDefIt++;
  featDef=*featDefIt;
  TEST_ASSERT(featDefIt!=factory->endFeatureDefs());
  TEST_ASSERT(featDef->getFamily()=="HBondAcceptor");

  featDefIt++;
  TEST_ASSERT(featDefIt==factory->endFeatureDefs());

  delete(factory);
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue347(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test Issue347" << std::endl;

  ROMol *testMol;
  MolChemicalFeatureFactory *factory;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;
  FeatSPtrList featSPtrs;
  boost::shared_ptr<MolChemicalFeature> featSPtr;


  std::string fdef;

  fdef=
    "DefineFeature CTriplet [C][C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0,1.0\n"
    "EndFeature\n"
    "DefineFeature CPair [C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0\n"
    "EndFeature\n";
  factory = buildFeatureFactory(fdef);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);
  featDefIt=factory->beginFeatureDefs();
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getFamily()=="LumpedHydrophobe");

  featDefIt++;
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getFamily()=="LumpedHydrophobe");

  featDefIt++;
  TEST_ASSERT(featDefIt==factory->endFeatureDefs());

  testMol=SmilesToMol("CCC");
  TEST_ASSERT(testMol);
  featSPtrs = factory->getFeaturesForMol(*testMol);
  TEST_ASSERT(featSPtrs.size()==1);
  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="LumpedHydrophobe");
  TEST_ASSERT(featSPtr->getType()=="CTriplet");

  // now reverse the order and we should get two matches:
  fdef=
    "DefineFeature CPair [C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0\n"
    "EndFeature\n"
    "DefineFeature CTriplet [C][C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0,1.0\n"
    "EndFeature\n";
  factory = buildFeatureFactory(fdef);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);
  featDefIt=factory->beginFeatureDefs();
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getFamily()=="LumpedHydrophobe");

  featDefIt++;
  featDef=*featDefIt;
  TEST_ASSERT(featDef->getFamily()=="LumpedHydrophobe");

  featDefIt++;
  TEST_ASSERT(featDefIt==factory->endFeatureDefs());
  featSPtrs = factory->getFeaturesForMol(*testMol);
  TEST_ASSERT(featSPtrs.size()==3);

  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="LumpedHydrophobe");
  TEST_ASSERT(featSPtr->getType()=="CPair");

  featSPtr=*(++featSPtrs.begin());
  TEST_ASSERT(featSPtr->getFamily()=="LumpedHydrophobe");
  TEST_ASSERT(featSPtr->getType()=="CPair");
  
  featSPtr=*(++++featSPtrs.begin());
  TEST_ASSERT(featSPtr->getFamily()=="LumpedHydrophobe");
  TEST_ASSERT(featSPtr->getType()=="CTriplet");


  
  delete(factory);
  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIssue348(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Issue 348." << std::endl;


  MatchVectType mv;
  std::string inText;
  int res;
  MolChemicalFeatureDef::CollectionType featureDefs;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;


  inText=
    "DefineFeature CTriplet \\\n"
    "[C][C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0,1.0\n"
    "EndFeature\n"
    "DefineFeature CTriplet [C[C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0\n"
    "EndFeature\n";

    
  bool ok;
  try{
    res=parseFeatureData(inText,featureDefs);
    ok=false;
  } catch (FeatureFileParseException &fpe){
    ok=true;
    TEST_ASSERT(fpe.lineNo()==6);
  }
  TEST_ASSERT(ok);

}

void testNestedAtomTypes(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Nested AtomType definitions." << std::endl;


  MatchVectType mv;
  std::string inText;
  int res;
  MolChemicalFeatureDef::CollectionType featureDefs;
  MolChemicalFeatureDef::CollectionType::const_iterator featDefIt;
  MolChemicalFeatureDef::CollectionType::value_type featDef;


  inText=
    "AtomType Carbon_NonPolar [C;!$(C=[O,N,P,S]);!$(C#N)]\n"
    "AtomType Hphobe [{Carbon_NonPolar},c]\n"
    "AtomType RingHphobe [{Hphobe};R]\n"
    "DefineFeature HP1 [{RingHphobe}]\n"
    "  Family Hydrophobe\n"
    "  Weights 1.0\n"
    "EndFeature\n";
    
  res=parseFeatureData(inText,featureDefs);
  TEST_ASSERT(!res);
  TEST_ASSERT(featureDefs.size()==1);
  std::cerr << (*featureDefs.begin())->getSmarts() << std::endl;
  

}

void testGithub252(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Testing Github Issue #252: crash when calling getPos() with no conformer." << std::endl;
  BOOST_LOG(rdErrorLog) << "     expect a precondition failure message below" << std::endl;

  ROMol *testMol;
  MolChemicalFeatureFactory *factory;
  FeatSPtrList featSPtrs;
  boost::shared_ptr<MolChemicalFeature> featSPtr;

  std::string fdef=
    "DefineFeature CTriplet [C][C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0,1.0\n"
    "EndFeature\n"
    "DefineFeature CPair [C][C]\n"
    "  Family LumpedHydrophobe\n"
    "  Weights 1.0,1.0\n"
    "EndFeature\n";
  factory = buildFeatureFactory(fdef);
  TEST_ASSERT(factory);
  TEST_ASSERT(factory->getNumFeatureDefs()==2);
  testMol=SmilesToMol("CCC");
  TEST_ASSERT(testMol);
  featSPtrs = factory->getFeaturesForMol(*testMol);
  TEST_ASSERT(featSPtrs.size()==1);
  featSPtr=*featSPtrs.begin();
  TEST_ASSERT(featSPtr->getFamily()=="LumpedHydrophobe");
  bool ok=false;
  try{
    featSPtr->getPos();
  } catch (const Invar::Invariant &i){
    ok=true;
  }
  TEST_ASSERT(ok);
  BOOST_LOG(rdErrorLog) << "   Done" << std::endl;
}


//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main(){
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  testIssue224();
  testIssue225();
  testIssue346();
  testIssue347();
  testIssue348();
#endif
  testNestedAtomTypes();
  testGithub252();
}
