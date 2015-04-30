// $Id$
//
//  Copyright (C) 2004-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/FileParseException.h>
#include "RDDepictor.h"
#include "DepictUtils.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/Conformer.h>
#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <RDGeneral/utils.h>
#include <stdlib.h>

#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;


using namespace RDKit;

void test1() {
  //std::string smiString = "OCC(F)(F)C(F)(F)C(F)(F)C(F)F";
  
  std::string smiString = "CN([CH](Cc3ccc(OS(c2cccc1c2ccnc1)(=O)=O)cc3)C(=O)N5CCN(c4ccccc4)CC5)S(c7cccc6c7ccnc6)(=O)=O \
                           C[n](n1)ccc1NC(=O)C(/C#N)=C/c2[s]cc(c2)Br O=C2C1CC3CC(C1)CC2C3 \
                           c8cc9cc3c5c(c(=O)cc4n1nc6c2c(ccc(c12)c(cc3)c45)cc7c6cc(=O)cc7)c9cc8 \
                           c1cccc(c1)CCCCNC(=O)C(/C#N)=C/c2ccc(O)c(c2)O \
                           COc(cc1)ccc1C(=C2)C=C(NC2=C)c3ccc(OC)cc3O \
                           CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O \
                           Cl.C[N+](C)(C)CCO O=C2C1CC3CC(C1)CC2C3 \
                           CN3CCC25[CH]4CCC(=O)C5(C)Oc1c(O)ccc(c12)C[CH]34 \
                           c1ccccc1\\C=C/C=C\\Cl.c1ccccc1C=C(Cl)C1SCCC1.Cl\\C=C1/SCCC1 Cl\\C=C/Br \
                           c1ccccc1\\C=C/C=C\\C=C/Cl c1ccccc1C=C(Cl)C1SCCC1 Cl\\C=C1/SCCC1 Cl\\C=C/Br \
                           CN2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3 \
                           N2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3 \
                           C2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3 \
                           ClC=C1SCCC1 C/C=C/C(C1CCCCCC1)=O c1ccccc1\\C=C/C=C\\C=C/Cl \
                           C[n](n1)ccc1NC(=O)C(/C#N)=C/c2[s]cc(c2)Br \
                           C1CC(C(Cl)(Br)F)CC1 \
                           C1(C2)CCC2C=C1 \
                           OC2C1(C)CCCCC1CCC2 \
                           CC12CCCC3OC1CCCC23 \
                           CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O \
                           ON=C(CC1=CC=CC=C1)[CH](C#N)C2=CC=CC=C2 \
                           COc(cc1)ccc1C(=C2/C#N)\\C=C(NC2=C(C#N)C#N)\\c3ccc(OC)cc3O \
                           COc(cc1)ccc1C(=C2)C=C(NC2=C)c3ccc(OC)cc3O";
                          
                          /*C1=CC=CC=C1CCCC2=CC(C=CC=C3)=C3C=C2 C/C=C/C(C1CCCCCC1)=O C/C=C/C(C1CCCCCC1)=O.c1ccccc1CC(Cl)=O \
                          [I-].CCC[N+](C)(CCC)CCC C1COCCN2CCOCCN(CCO1)CCOCCOCC2 C1CN2CCN1CC2 \
                          C(CCCCCCC)(=O)CCCCCC \
                           ClCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCF \
                           C(CCCCCCC)(=O)CCCCCC \
                           C/C=C/C(C1CCCCCC1)=O C1CC(CC12)C=C2 c1ccccc1\\C=C(Cl)/C#N N#C\\C=C(Cl)/Cc1ccccc1 \
                           CN(C)c(cc1)ccc1\\C=C(/C#N)c(n2)c(C#N)c([n]3cccc3)[n]2c(c4)cccc4 \
                           CCOC(=O)CN1C(=O)/C=C(/C)c(c12)c(C)n[n]2C C/12=C\\C(=O)c3cc(OC(F)(F)F)ccc3N2C(=O)/C(CC(=O)OC)=C\\1C(=O)OCC \
                           C/12=C(\\NC(N2)=O)NC(=O)NC1=O F\\C=C/Cl F/C=C/Cl c1ccccc1\\C=C/C c1ccccc1\\C=C\\C=C\\C=C\\Cl \
                           c1ccccc1\\C=C/C=C\\C=C/Cl c1ccccc1\\C=C\\C=C(O)\\C=C(Br)\\Cl \
                           c1ccccc1\\C=C/C=C(O)\\C=C(Br)/Cl \
                           CC#CCC O=C=O C1=CC=CC=C1 C1CCC1 C1CC1(C#CC) C1CCCCCCC1 C1=CC=CC=C1(CCCCCCC) \
                           C1=CC=CC=C1(CC(CCC)CCC(C)(C)C) \
                           C1=CC=CC=C1CCCC2=CC(C=CC=C3)=C3C=C2 \
                           C1=CC=CC=C1(CC(CCC)CCC(CCC)(CCC)C) \
                           C1=CC(C=CC=C2)=C2C=C1 C1=CC(C=C2)=C2C=C1 \
                           C1=CC(CCC2)=C2C=C1 C1(CCC3)CCCC2C1C3CCC2 \
                           C12=CC=C3C(C4=C(C=CC=C5)C5=C3)C1C(C=C4)=CC=C2 \
                           C12CCCC3(CCCCC3)C1CCCC2 \
                           C1CCCC2(CC2)C1 C12C3C4C1C5C2C3C45 \
                           C2(C=C(C=C5)NC5=C4)=CC(C=C2)=CC1=CC=C(C=C3C=CC4=N3)N1";
    */
  std::string rdbase = getenv("RDBASE");
  std::string ofile = rdbase + "/Code/GraphMol/Depictor/test_data/test1.out.sdf";
  SDWriter writer(ofile);

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();
      ++token){
    std::string smi=*token;
    RWMol *m = SmilesToMol(smi, 0, 1); 
    TEST_ASSERT(m)
    RDDepict::compute2DCoords(*m);
    writer.write(*m);
    delete m;
  }
}

void testCollisions() {
  std::string smiString = "CN([CH](Cc3ccc(OS(c2cccc1c2ccnc1)(=O)=O)cc3)C(=O)N5CCN(c4ccccc4)CC5)S(c7cccc6c7ccnc6)(=O)=O \
                           CC(C)C12CCC(C)(CC1)CC2 \
                           CC(C)C1CCC2(CC)CCCC1CCC2 \
                           CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O \
                           OC(=O)C1=C(C=CC=C1)C2=C3C=CC(=O)C(=C3OC4=C2C=CC(=C4Br)O)Br \
                           CC(CC(=O)OCC(F)(F)C(F)(F)C(F)(F)C(F)F)CC(=O)OCC(F)(F)C(F)(F)C(F)(F)C(F)F \
                           CC(=O)Nc1ccccc1C(=O)C1CCCC1C1CCCCC1 \
                           Cc1ccc(C(c2ccc(C)o2)c2ccccc2[N+](=O)[O-])o1 \
    CN2C3CC(OC(=O)C(CO)c1ccccc1)CC2CC3";
  // this used to be a test, but it's currently failing catastrophically:
  //      CC(C)(C)O[Si](OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F)(OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F)OC(C)(C)C

  
  std::string rdbase = getenv("RDBASE");
  std::string ofile = rdbase + "/Code/GraphMol/Depictor/test_data/collisions.out.sdf";
  SDWriter writer(ofile);

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();
      ++token){
    std::string smi=*token;
    RWMol *m = SmilesToMol(smi, 0, 1);
    TEST_ASSERT(m);
    unsigned int confId = RDDepict::compute2DCoords(*m);
    // check that there are no collisions in the molecules
    const Conformer &conf = m->getConformer(confId);

    int natms = m->getNumAtoms();
    for (int i = 0; i < natms; i++) {
      RDGeom::Point3D loci = conf.getAtomPos(i);
      for (int j = i+1; j < natms; j++) {
        RDGeom::Point3D locj = conf.getAtomPos(j);
        locj -= loci;
        CHECK_INVARIANT(locj.length() > 0.35, "");
      }
    }
    writer.write(*m);
    delete m;
  }
}

void testAddHs() {
  // test for issue 193
  std::string smi = "F[C@H](Cl)Br";
  RWMol *m = SmilesToMol(smi, 0, 1);
  MolOps::addHs(*m);
  RDDepict::compute2DCoords(*m);
  delete m;
}

void testIssue198() {
  std::string smi = "Cl.C[N+](C)(C)CCO";
  RWMol *m = SmilesToMol(smi, 0, 1);
  unsigned int confId = RDDepict::compute2DCoords(*m);
  RDKit::ROMol::AtomIterator ai;
  const Conformer &conf = m->getConformer(confId);
  for (ai = m->beginAtoms(); ai != m->endAtoms(); ai++) {
    RDGeom::Point3D loc = conf.getAtomPos((*ai)->getIdx());
    CHECK_INVARIANT(loc.x < 100.0, "");
    CHECK_INVARIANT(loc.y < 100.0, "");
    CHECK_INVARIANT(loc.z < 100.0, "");
  }
  delete m;
}

void test2() {
  std::string rdbase = getenv("RDBASE");
  std::string smifile = rdbase + "/Code/GraphMol/Depictor/test_data/first_200.tpsa.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, -1);
  
  std::string ofile = rdbase + "/Code/GraphMol/Depictor/test_data/first_200.out.sdf";
  SDWriter writer(ofile);
  ROMol *mol;
  while (1) {
    try {
      mol = smiSup.next();
      std::string mname;
      //wmol = static_cast<RWMol *>(mol);
      RDDepict::compute2DCoords(*mol);
      writer.write(*mol);
    } catch (FileParseException &) {
      break;
    }
  }
}

void test3() {
  std::string rdbase = getenv("RDBASE");
  std::string smifile = rdbase + "/Code/GraphMol/Depictor/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);
  
  std::string ofile = rdbase + "/Code/GraphMol/Depictor/test_data/cis_trans_cpp.out.sdf";
  SDWriter writer(ofile);
  ROMol *mol;
  while (1) {
    try {
      mol = smiSup.next();
      std::string mname;
      //wmol = static_cast<RWMol *>(mol);
      RDDepict::compute2DCoords(*mol);
      writer.write(*mol);
    } catch (FileParseException &) {
      break;
    }
  }
}

void _compareCoords(const ROMol *mol1, unsigned int cid1, 
                    const ROMol *mol2, unsigned int cid2,
                    double tol = 0.01) {
  unsigned int nat = mol1->getNumAtoms();
  CHECK_INVARIANT(nat == mol2->getNumAtoms(), "");
 
  const RDKit::Conformer &conf1 = mol1->getConformer(cid1);
  const RDKit::Conformer &conf2 = mol1->getConformer(cid2);

  for (unsigned int i = 0; i < nat; i++) {
    RDGeom::Point3D pt1 = conf1.getAtomPos(i);
    RDGeom::Point3D pt2 = conf2.getAtomPos(i);
    pt2 -= pt1;
    
    CHECK_INVARIANT(fabs(pt2.x) < tol, "");
    CHECK_INVARIANT(fabs(pt2.y) < tol, "");
    //if ((fabs(pt2.x) >= tol) || (fabs(pt2.y) >= tol) ) {
    //  BOOST_LOG(rdInfoLog)<< pt1 << " " << pt2 << "\n";
    //}
  }
}

void test4() {
  // test prespecified coordinates for various smiles
  RDGeom::INT_POINT2D_MAP crdMap;
  crdMap[0] = RDGeom::Point2D(3.52, 1.30); crdMap[1] = RDGeom::Point2D(2.77, 0.0);
  crdMap[2] = RDGeom::Point2D(1.27, 0.0); crdMap[3] = RDGeom::Point2D(0.39, 1.21);
  crdMap[4] = RDGeom::Point2D(-1.03, 0.75); crdMap[5] = RDGeom::Point2D(-1.03, -0.75);
  crdMap[6] = RDGeom::Point2D(0.39, -1.21);

  // self check
  std::string smi = "Cl\\C=C1/SCCC1";
  RWMol *m1 = SmilesToMol(smi, 0, 1); 
  unsigned int cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  RWMol *mref = SmilesToMol(smi, 0, 1);
  unsigned int cid2 = RDDepict::compute2DCoords(*mref);
  MolToMolFile(*mref, "junk.mol");
  //_compareCoords(m1, mref);

  delete m1;
  
  // now lets remove some of the coordinates
  crdMap.erase(crdMap.find(3));
  crdMap.erase(crdMap.find(4));
  m1 = SmilesToMol(smi, 0, 1); 
  cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  //MolToMolFile(m1, "junk.mol");
  //_compareCoords(m1, mref);
  delete m1;

  //little bit more complicated
  smi = "Cl\\C=C1/SCCC1(CC1CC1)";
  m1 = SmilesToMol(smi, 0, 1); 
  cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  //MolToMolFile(m1, "junk1.mol");
  delete mref;
  mref = SmilesToMol(smi, 0, 1);
  RDDepict::compute2DCoords(*mref, 0, false);
  _compareCoords(m1, cid1, mref, cid2);
  delete m1;

  // little more complicate we will specify coordinates from half of both the rings
  crdMap[7] = RDGeom::Point2D(0.85, -2.64);
  crdMap[8] = RDGeom::Point2D(-0.15, -3.75);
  crdMap[9] = RDGeom::Point2D(-0.46, -5.22);
  m1 = SmilesToMol(smi, 0, 1);
  cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  //MolToMolFile(m1, "junk1.mol");
  //MolToMolFile(mref, "junk2.mol");
  _compareCoords(m1, cid1, mref, cid2);
  delete m1;

  // one final test for a case we know is a pain
  smi = "C1CCCC2(Cl)(CCCCCC12)";
  crdMap.clear();
  crdMap[0]=RDGeom::Point2D(-0.83,-3.12); crdMap[1]=RDGeom::Point2D(0.19,-4.22);
  crdMap[2]=RDGeom::Point2D(1.66,-3.88); crdMap[3]=RDGeom::Point2D(2.10, -2.45); 
  crdMap[4]=RDGeom::Point2D(1.08,-1.35); crdMap[5]=RDGeom::Point2D(2.56, -1.12);
  crdMap[6]=RDGeom::Point2D(1.73,0.00); crdMap[7]=RDGeom::Point2D(1.08, 1.35);
  crdMap[8]=RDGeom::Point2D(-0.38,1.69); crdMap[9]=RDGeom::Point2D(-1.56, 0.75);
  crdMap[10]=RDGeom::Point2D(-1.56,-0.75); crdMap[11]=RDGeom::Point2D(-0.38, -1.69);
  
  delete mref;
  mref = SmilesToMol(smi, 0, 1);
  cid2 = RDDepict::compute2DCoords(*mref, 0, false);
  crdMap.erase(crdMap.find(5));
  //MolToMolFile(mref, "junk1.mol");
  m1 = SmilesToMol(smi, 0, 1);
  cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  _compareCoords(m1, cid1, mref, cid2);
  delete m1;
  delete mref;
}

void tempTest() {
  int i;
  RDGeom::Point3D pt1, pt2, pt3;
  double cosT, sinT;
  RDGeom::Transform3D trans;
  RDGeom::Point3D rotnAxis;
  double dt;
  for (i = 0; i < 100; i++) {
    pt1.x = (double)(rand()%1000);
    pt1.y = (double)(rand()%1000);
    pt1.z = (double)(rand()%1000);
    pt1.normalize();
    pt2.x = (double)(rand()%1000);
    pt2.y = (double)(rand()%1000);
    pt2.z = (double)(rand()%1000);
    pt2.normalize();
    
    cosT = -pt1.dotProduct(pt2);
    sinT = sqrt(1.0-cosT*cosT);
    rotnAxis = pt1.crossProduct(pt2);
    rotnAxis.normalize();
    trans.setToIdentity();
    trans.SetRotation(cosT,sinT, rotnAxis);
    pt3 = pt2;
    trans.TransformPoint(pt3);
    dt = pt1.dotProduct(pt3);
    if (fabs(dt + 1.0) > 1.0e-3) {
      BOOST_LOG(rdInfoLog)<< i << " (" << pt1 << ") (" << pt2 << ") (" 
                << pt3 << ") (" << rotnAxis << ") " << dt << "\n";
    }
  }
}

void testIssue248() {
  std::string smiString = "OCC(=O)C1(O)[CH](O)C[CH]2[CH]3CCC4=CC(=O)C=CC4(C)C3(F)[CH](O)CC21C \
   CC(C)[CH](OC(C)=O)c1nc2ccccc2c(=O)n1[CH]1CC2(OC1=O)c1c(cccc1)N1C(=O)C(C)(C)N(O)[CH]12 \
   CC(CC(=O)OCC(F)(F)C(F)(F)C(F)(F)C(F)F)CC(=O)OCC(F)(F)C(F)(F)C(F)(F)C(F)F \
   CC(O[CH]1C(=O)C2(C)[CH](O)C[CH]3OCC3(OC(C)=O)[CH]2[CH](OC(=O)c2ccccc2)C2(O)C[CH](OC([CH](O)[CH](NC(=O)c3ccccc3)c3ccccc3)=O)C(C)=C1C2(C)C)=O \
  COc1c2O[CH]3[CH](O)C=C[CH]4C33CCN(C)[CH]4Cc(c32)cc1 \
  CN1CCC23[CH]4Oc5c(O)ccc(c25)C[CH]1[CH]3C=C[CH]4O \
  CN1C2CCC1C[CH](OC(=O)C(CO)c1ccccc1)C2 \
  COC([CH]1C2CCC(C[CH]1OC(c1ccccc1)=O)N2C)=O \
  OCC(=O)C1(O)[CH](O)C[CH]2[CH]3CCC4=CC(=O)C=CC4(C)C3(F)[CH](O)CC21C";
  //COC12C3Oc4c(O)ccc5c4C33CCN(CC4CC4)[CH](C5)C3(C[CH]1C(O)(C)C(C)(C)C)CC2";

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();
      ++token){
    std::string smi=*token;
    RWMol *m = SmilesToMol(smi, 0, 1); 
    unsigned int confId = RDDepict::compute2DCoords(*m);
    // check that there are no collisions in the molecules
    int natms = m->getNumAtoms();
    int i, j;
    for (i = 0; i < natms; i++) {
      const Conformer &conf = m->getConformer(confId);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      for (j = i+1; j < natms; j++) {
        RDGeom::Point3D locj = conf.getAtomPos(j);
        locj -= loci;
        if(locj.length()<=0.30){
          std::cout << "mismatch: " << i << " " << j << " " << locj.length() << std::endl;
          std::cout << "\t" << smi << std::endl;
          std::cout << MolToMolBlock(*m,true,confId)<<std::endl;
        }
        CHECK_INVARIANT(locj.length() > 0.30, "");
      }
    }
    delete m;
  }
}

void testQueries() {
  std::string smaString = "C[C,N] c1cc[c,n]cc1C";

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smaString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();
      ++token){
    std::string sma=*token;
    RWMol *m = SmartsToMol(sma);
    unsigned int confId = RDDepict::compute2DCoords(*m);
    // check that there are no collisions in the molecules
    int natms = m->getNumAtoms();
    int i, j;
    for (i = 0; i < natms; i++) {
      const Conformer &conf = m->getConformer(confId);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      for (j = i+1; j < natms; j++) {
        RDGeom::Point3D locj = conf.getAtomPos(j);
        locj -= loci;
        if(locj.length()<=0.30){
          std::cout << "mismatch: " << i << " " << j << " " << locj.length() << std::endl;
          std::cout << "\t" << sma << std::endl;
          std::cout << MolToMolBlock(*m,true,confId)<<std::endl;
        }
        CHECK_INVARIANT(locj.length() > 0.30, "");
      }
    }
    delete m;
  }
}

void testRemoveHsCrash() {
  std::string rdbase = getenv("RDBASE");
  std::string molfile = rdbase + "/Code/GraphMol/Depictor/test_data/hs_crash.mol";
  RWMol *m=MolFileToMol(molfile,true,false);
  TEST_ASSERT(m);
  ROMol *newM=MolOps::removeHs(*static_cast<ROMol *>(m));
  delete m;
  unsigned int confId = RDDepict::compute2DCoords(*newM);
  TEST_ASSERT(confId>=0);
  delete newM;
}
    

void testIssue2091304() {
  // the problem here was a crash, so just finishing is success.
  RDGeom::INT_POINT2D_MAP crdMap;
  crdMap[0] = RDGeom::Point2D(0., 1.50);

  std::string smi = "COC";
  RWMol *m1 = SmilesToMol(smi);
  RDDepict::compute2DCoords(*m1, &crdMap, false);
  delete m1;
}

void testIssue2303566() {
  RDGeom::INT_POINT2D_MAP crdMap;
  crdMap[0]=RDGeom::Point2D(0.0781827, -0.112109);
  crdMap[1]=RDGeom::Point2D(1.13205, -0.415657);
  crdMap[2]=RDGeom::Point2D(0.676323, -1.51974);
  crdMap[3]=RDGeom::Point2D(0.735363, -2.62024);
  crdMap[4]=RDGeom::Point2D(0.974216, -1.5726);
  crdMap[5]=RDGeom::Point2D(1.02041, -1.56945);
  std::string mb = "\n\
     RDKit          2D\n\
\n\
  6  7  0  0  0  0  0  0  0  0999 V2000\n\
    0.0782   -0.1121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.1321   -0.4157    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.6763   -1.5197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    0.7354   -2.6202    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.0098   -1.5757    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n\
    1.0218   -1.5749    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n\
  1  2  1  0\n\
  1  3  1  0\n\
  3  2  1  0\n\
  3  4  1  0\n\
  5  6  2  0\n\
  2  4  1  0\n\
  1  4  1  0\n\
M  END\n";

  RWMol *m1 = MolBlockToMol(mb);
  // in the original bug, this led to a crash, so just having it complete
  // means the test passed.
  unsigned int cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  TEST_ASSERT(cid1>=0);
  delete m1;
}

void testIssue2821647() {
  {
    std::string smi = "CCCCC";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    double xx=0,yy=0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx+=loci.x*loci.x;
      yy+=loci.y*loci.y;
    }
    TEST_ASSERT(xx>yy);
    delete m1;
  }
  {
    std::string smi = "c1ccccc1CCCCCC1CC1";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    double xx=0,yy=0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx+=loci.x*loci.x;
      yy+=loci.y*loci.y;
    }
    TEST_ASSERT(xx>yy);
    delete m1;
  }
  {
    std::string smi = "c1ccc2c(c1)oc1c3ccccc3oc21";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    double xx=0,yy=0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx+=loci.x*loci.x;
      yy+=loci.y*loci.y;
    }
    TEST_ASSERT(xx>yy);
    delete m1;
  }
  {
    std::string smi = "[H]n1c2ccccc2c2n([H])c3ccccc3c12";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    double xx=0,yy=0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx+=loci.x*loci.x;
      yy+=loci.y*loci.y;
    }
    TEST_ASSERT(xx>yy);
    delete m1;
  }
}

void testIssue2948402() {
  {
    std::string smi = "C1C2CC3=CC=CC(C2)CC(O3)C1";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    TEST_ASSERT(cid1>=0);
    delete m1;
  }
}

void testIssue2995724() {
  {
    // the original problem from Thomas Heller:
    std::string smi = "OC(=O)[C@@H]1CCCN1";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    
    const Conformer &conf = m1->getConformer(cid1);
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      RDGeom::Point3D loci = conf.getAtomPos(i);
      // we're testing for NaNs here:
      TEST_ASSERT(loci.x>-5.0);
      TEST_ASSERT(loci.x<5.0);
      TEST_ASSERT(loci.y>-5.0);
      TEST_ASSERT(loci.y<5.0);
    }
    delete m1;
  }
  {
    // additional test cases from Kirk DeLisle:
    std::string smis[]={"CCC(N1CCN(CC1)C(CC)O)O.Cl[Pt](Cl)(Cl)(Cl)(Cl)Cl",
                        "CN(C)S(=O)(=O)N1CCN(CC1)S(=O)(=O)N(C)C",
                        "Cc1ccc(cc1C)Nc2nc(nc(n2)N)CN3CCN(CC3)C",
                        "CC1(OC(=C(C(=O)O1)C2=NCCC2)O)C"};
    for(unsigned int j=0;j<4;j++){
      RWMol *m1 = SmilesToMol(smis[j]);
      unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    
      const Conformer &conf = m1->getConformer(cid1);
      for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
        RDGeom::Point3D loci = conf.getAtomPos(i);
        TEST_ASSERT(loci.x>-7.0);
        TEST_ASSERT(loci.x<7.0);
        TEST_ASSERT(loci.y>-6.0);
        TEST_ASSERT(loci.y<6.0);
      }
      delete m1;
    }
  }
}

void testBondLengthChange() {
  {
    std::string smi = "CC";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    
    const Conformer &conf = m1->getConformer(cid1);
    TEST_ASSERT(feq(conf.getAtomPos(0).x,-0.75));
    TEST_ASSERT(feq(conf.getAtomPos(1).x,0.75));

    delete m1;
  }
  {
    std::string smi = "CC";
    RWMol *m1 = SmilesToMol(smi);
    RDDepict::BOND_LEN=1.0;
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    
    const Conformer &conf = m1->getConformer(cid1);
    TEST_ASSERT(feq(conf.getAtomPos(0).x,-0.5));
    TEST_ASSERT(feq(conf.getAtomPos(1).x,0.5));

    delete m1;
  }
}

void testIssue3122141() {
  {
    std::string smi = "ClC1=C(Cl)C2(Cl)C3C(Cl)C=CC3C1(Cl)C2(Cl)Cl";
    RWMol *m1 = SmilesToMol(smi);
    RDGeom::INT_POINT2D_MAP crdMap;
    crdMap[1]=RDGeom::Point2D(1.0,2.0);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
    TEST_ASSERT(cid1>=0);
    delete m1;
  }
}


void testIssue3135833() {
  {
    std::string smi = "N#CC(/C=[N+](/c1cccnc1)c1c2cc(C(=O)c3ncsc3-c3cccnc3)ccc2ccc1)=C(/[NH3+])[S-]";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    
    RDGeom::INT_POINT2D_MAP crdMap;
    crdMap[32]=RDGeom::Point2D(-2.8518,1.0270);
    crdMap[33]=RDGeom::Point2D(-2.7999,-0.4721);
    crdMap[4]=RDGeom::Point2D(-1.4239,-2.6758);
    crdMap[11]=RDGeom::Point2D(-1.4757,-1.1767);
    crdMap[12]=RDGeom::Point2D(-0.2034,-0.3823);
    crdMap[13]=RDGeom::Point2D(1.1208,-1.0869);
    crdMap[14]=RDGeom::Point2D(2.3931,-0.2924);
    crdMap[15]=RDGeom::Point2D(3.7173,-0.9971);
    crdMap[28]=RDGeom::Point2D(2.3412,1.2067);
    crdMap[29]=RDGeom::Point2D(1.0170,1.9113);
    crdMap[30]=RDGeom::Point2D(-0.2553,1.1168);
    crdMap[31]=RDGeom::Point2D(-1.5795,1.8215);

    unsigned int cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
    TEST_ASSERT(cid1>=0);
    delete m1;
  }
}

void testIssue3487469() {
  {
    std::string smi = "C*C";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(1)->getHybridization()==Atom::UNSPECIFIED);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1,0,true);
    const Conformer &conf = m1->getConformer(cid1);
    RDGeom::Point3D p0 = conf.getAtomPos(0);
    RDGeom::Point3D p1 = conf.getAtomPos(1);
    RDGeom::Point3D p2 = conf.getAtomPos(2);

    RDGeom::Point3D v1=p0-p1,v2=p2-p1;
    v1.normalize();
    v2.normalize();
    TEST_ASSERT(feq(v1.dotProduct(v2),-0.5,.01))
    delete m1;
  }
}

void testGitHubIssue8() {
  {
    std::string smi = "[I-].C[n+]1c(\\C=C\\2/C=CC=CN2CC=C)sc3ccccc13";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    RWMol *p = SmilesToMol("[I-]");
    TEST_ASSERT(p);
    ROMol *m2 = deleteSubstructs(*m1,*p);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getNumAtoms()==m1->getNumAtoms()-1);
    unsigned int cid1 = RDDepict::compute2DCoords(*m2);
    TEST_ASSERT(cid1==0);
    delete m1;
    delete m2;
    delete p;
  }
}

void testGitHubIssue78() {
  {  // the basic test: the smallest reproducible:
    std::string smi = "C3CCCC1C3C(O2)C2CC1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1);
    const Conformer &conf = m1->getConformer(cid1);
    RDGeom::Point3D p5 = conf.getAtomPos(5);
    RDGeom::Point3D p6 = conf.getAtomPos(6);
    RDGeom::Point3D p7 = conf.getAtomPos(7);
    RDGeom::Point3D p8 = conf.getAtomPos(8);
    RDGeom::Point3D p9 = conf.getAtomPos(9);

    RDGeom::Point3D p87 = p8-p7;
    RDGeom::Point3D p86 = p8-p6;
    RDGeom::Point3D p89 = p8-p9;
    TEST_ASSERT(p87.dotProduct(p86)*p87.dotProduct(p89)<0);
    
    RDGeom::Point3D p67 = p6-p7;
    RDGeom::Point3D p68 = p6-p8;
    RDGeom::Point3D p65 = p6-p5;
    TEST_ASSERT(p67.dotProduct(p68)*p67.dotProduct(p65)<0);

    delete m1;
  }
  {  // a collection of previous failures
    std::string smis[4] = {"C1=CC=C2C(=C1)C3=CC=CC=C3C4=C2C(C(C5C4O5)O)O",
                          "CC1=C2C=C(C=CC2=NC3=C1C4=C(O4)C5=CC=CC=C53)CO",
                          "CC1=C2C=CC3=C(C2=CC4=CC=CC=C14)C5C(O5)C(C3O)O",
                          "C1=CC=C2C(=C1)C=C3C=CC4=C5C3=C2C6C(C5=CC=C4)O6"};
    RWMol *p = SmartsToMol("[#6]~[#6]~1-[#8]-[#6]~1~[#6]");
    TEST_ASSERT(p);
    for(unsigned int i=0;i<4;++i){
      RWMol *m = SmilesToMol(smis[i]);
      TEST_ASSERT(m);
      MatchVectType mv;
      TEST_ASSERT(SubstructMatch(*m,*p,mv));
      TEST_ASSERT(mv.size()==5);
                  
      unsigned int cid1 = RDDepict::compute2DCoords(*m);
      const Conformer &conf = m->getConformer(cid1);
      RDGeom::Point3D v10 = conf.getAtomPos(mv[1].second)-conf.getAtomPos(mv[0].second);
      RDGeom::Point3D v12 = conf.getAtomPos(mv[1].second)-conf.getAtomPos(mv[2].second);
      RDGeom::Point3D v13 = conf.getAtomPos(mv[1].second)-conf.getAtomPos(mv[3].second);
      TEST_ASSERT(v12.dotProduct(v10)*v12.dotProduct(v13)<0);
      
      RDGeom::Point3D v31 = conf.getAtomPos(mv[3].second)-conf.getAtomPos(mv[1].second);
      RDGeom::Point3D v32 = conf.getAtomPos(mv[3].second)-conf.getAtomPos(mv[2].second);
      RDGeom::Point3D v34 = conf.getAtomPos(mv[3].second)-conf.getAtomPos(mv[4].second);
      TEST_ASSERT(v32.dotProduct(v32)*v32.dotProduct(v34)<0);
      
      
      delete m;
    }
    delete p;
  }
}

int main() { 
  RDLog::InitLogs();
#if 1
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   test1 \n";
  test1();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n\n";
  
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   testCollisions \n";
  testCollisions();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   testAddHs \n";
  testAddHs();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   test198 \n";
  testIssue198();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   test2 \n";
  test2();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   test3 \n";
  test3();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   test4 \n";
  test4();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   tempTest \n";
  tempTest();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Issue248 \n";
  testIssue248();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Queries \n";
  testQueries();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test crashes associated with RemoveHs \n";
  testRemoveHsCrash();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 2091304\n";
  testIssue2091304();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 2303566\n";
  testIssue2303566();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 2821647\n";
  testIssue2821647();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 2948402\n";
  testIssue2948402();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 2995724\n";
  testIssue2995724();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Testing a change of the depictor bond length\n";
  testBondLengthChange();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 3122141\n";
  testIssue3122141();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 3135833\n";
  testIssue3135833();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test Issue 3487469\n";
  testIssue3487469();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";

  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test GitHub Issue 8\n";
  testGitHubIssue8();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
#endif
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";
  BOOST_LOG(rdInfoLog)<< "   Test GitHub Issue 78\n";
  testGitHubIssue78();
  BOOST_LOG(rdInfoLog)<< "***********************************************************\n";


  return(0);
}

