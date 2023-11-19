//
//  Copyright (C) 2004-2018 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
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
#include "Templates.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <RDGeneral/utils.h>
#include <cstdlib>
#include <cmath>

#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

using namespace RDKit;

auto defaultRDKitBondLen = RDDepict::BOND_LEN;

void _compareCoords(const ROMol *mol1, unsigned int cid1, const ROMol *mol2,
                    unsigned int cid2, double tol = 0.01) {
  unsigned int nat = mol1->getNumAtoms();
  CHECK_INVARIANT(nat == mol2->getNumAtoms(), "");

  const RDKit::Conformer &conf1 = mol1->getConformer(cid1);
  const RDKit::Conformer &conf2 = mol2->getConformer(cid2);

  for (unsigned int i = 0; i < nat; i++) {
    RDGeom::Point3D pt1 = conf1.getAtomPos(i);
    RDGeom::Point3D pt2 = conf2.getAtomPos(i);
    pt2 -= pt1;

    if ((fabs(pt2.x) >= tol) || (fabs(pt2.y) >= tol)) {
      std::cerr << MolToMolBlock(*mol1, cid1) << std::endl;
      std::cerr << MolToMolBlock(*mol2, cid2) << std::endl;
      break;
    }
    CHECK_INVARIANT(fabs(pt2.x) < tol, "");
    CHECK_INVARIANT(fabs(pt2.y) < tol, "");
  }
}

void test1() {
  // std::string smiString = "OCC(F)(F)C(F)(F)C(F)(F)C(F)F";

  std::string smiString =
      "CN([CH](Cc3ccc(OS(c2cccc1c2ccnc1)(=O)=O)cc3)C(=O)N5CCN(c4ccccc4)CC5)S(c7cccc6c7ccnc6)(=O)=O \
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
                           COc(cc1)ccc1C(=C2)C=C(NC2=C)c3ccc(OC)cc3O \
                          C1=CC=CC=C1CCCC2=CC(C=CC=C3)=C3C=C2 C/C=C/C(C1CCCCCC1)=O \
                          C/C=C/C(C1CCCCCC1)=O.c1ccccc1CC(Cl)=O \
                          [I-].CCC[N+](C)(CCC)CCC C1COCCN2CCOCCN(CCO1)CCOCCOCC2 C1CN2CCN1CC2 \
                          C(CCCCCCC)(=O)CCCCCC \
                          ClCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCF \
                          C(CCCCCCC)(=O)CCCCCC \
                          C/C=C/C(C1CCCCCC1)=O C1CC(CC12)C=C2 c1ccccc1\\C=C(Cl)/C#N \
                          N#C\\C=C(Cl)/Cc1ccccc1 \
                          CN(C)c(cc1)ccc1\\C=C(/C#N)c(n2)c(C#N)c([n]3cccc3)[n]2c(c4)cccc4 \
                          CCOC(=O)CN1C(=O)/C=C(/C)c(c12)c(C)n[n]2C \
                          C/12=C\\C(=O)c3cc(OC(F)(F)F)ccc3N2C(=O)/C(CC(=O)OC)=C\\1C(=O)OCC \
                          C/12=C(\\NC(N2)=O)NC(=O)NC1=O F\\C=C/Cl F/C=C/Cl c1ccccc1\\C=C/C \
                          c1ccccc1\\C=C\\C=C\\C=C\\Cl \
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

  std::string rdbase = getenv("RDBASE");
  std::string ofile =
      rdbase + "/Code/GraphMol/Depictor/test_data/test1.out.sdf";
  SDWriter writer(ofile);
  std::string ifile = rdbase + "/Code/GraphMol/Depictor/test_data/test1.sdf";
  SDMolSupplier suppl(ifile);
  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString, spaceSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string smi = *token;
    std::unique_ptr<RWMol> m{SmilesToMol(smi)};
    TEST_ASSERT(m)
    RDDepict::compute2DCoords(*m);
    writer.write(*m);
    std::unique_ptr<ROMol> ref{suppl.next()};
    _compareCoords(m.get(), 0, ref.get(), 0);
  }
}

void testCollisions() {
  std::string smiString =
      "CN([CH](Cc3ccc(OS(c2cccc1c2ccnc1)(=O)=O)cc3)C(=O)N5CCN(c4ccccc4)CC5)S(c7cccc6c7ccnc6)(=O)=O \
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
  std::string ofile =
      rdbase + "/Code/GraphMol/Depictor/test_data/collisions.out.sdf";
  SDWriter writer(ofile);

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString, spaceSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string smi = *token;
    RWMol *m = SmilesToMol(smi, 0, 1);
    TEST_ASSERT(m);
    unsigned int confId = RDDepict::compute2DCoords(*m);
    // check that there are no collisions in the molecules
    const Conformer &conf = m->getConformer(confId);

    int natms = m->getNumAtoms();
    for (int i = 0; i < natms; i++) {
      RDGeom::Point3D loci = conf.getAtomPos(i);
      for (int j = i + 1; j < natms; j++) {
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
  std::string smifile =
      rdbase + "/Code/GraphMol/Depictor/test_data/first_200.tpsa.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, -1);

  std::string ofile =
      rdbase + "/Code/GraphMol/Depictor/test_data/first_200.out.sdf";
  SDWriter writer(ofile);
  ROMol *mol;
  while (1) {
    try {
      mol = smiSup.next();
      std::string mname;
      // wmol = static_cast<RWMol *>(mol);
      RDDepict::compute2DCoords(*mol);
      writer.write(*mol);
      delete mol;
    } catch (FileParseException &) {
      break;
    }
  }
}

void test3() {
  std::string rdbase = getenv("RDBASE");
  std::string smifile =
      rdbase + "/Code/GraphMol/Depictor/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);

  std::string ofile =
      rdbase + "/Code/GraphMol/Depictor/test_data/cis_trans_cpp.out.sdf";
  SDWriter writer(ofile);
  ROMol *mol;
  while (1) {
    try {
      mol = smiSup.next();
      std::string mname;
      // wmol = static_cast<RWMol *>(mol);
      RDDepict::compute2DCoords(*mol);
      writer.write(*mol);
      delete mol;
    } catch (FileParseException &) {
      break;
    }
  }
}

void test4() {
  // test prespecified coordinates for various smiles
  RDGeom::INT_POINT2D_MAP crdMap;
  crdMap[0] = RDGeom::Point2D(3.52, 1.30);
  crdMap[1] = RDGeom::Point2D(2.77, 0.0);
  crdMap[2] = RDGeom::Point2D(1.27, 0.0);
  crdMap[3] = RDGeom::Point2D(0.39, 1.21);
  crdMap[4] = RDGeom::Point2D(-1.03, 0.75);
  crdMap[5] = RDGeom::Point2D(-1.03, -0.75);
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
  // MolToMolFile(m1, "junk.mol");
  //_compareCoords(m1, mref);
  delete m1;

  // little bit more complicated
  smi = "Cl\\C=C1/SCCC1(CC1CC1)";
  m1 = SmilesToMol(smi, 0, 1);
  cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  // MolToMolFile(m1, "junk1.mol");
  delete mref;
  mref = SmilesToMol(smi, 0, 1);
  RDDepict::compute2DCoords(*mref, nullptr, false);
  _compareCoords(m1, cid1, mref, cid2);
  delete m1;

  // little more complicate we will specify coordinates from half of both the
  // rings
  crdMap[7] = RDGeom::Point2D(0.85, -2.64);
  crdMap[8] = RDGeom::Point2D(-0.15, -3.75);
  crdMap[9] = RDGeom::Point2D(-0.46, -5.22);
  m1 = SmilesToMol(smi, 0, 1);
  cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
  // MolToMolFile(m1, "junk1.mol");
  // MolToMolFile(mref, "junk2.mol");
  _compareCoords(m1, cid1, mref, cid2);
  delete m1;

  // one final test for a case we know is a pain
  smi = "C1CCCC2(Cl)(CCCCCC12)";
  crdMap.clear();
  crdMap[0] = RDGeom::Point2D(-0.83, 3.12);
  crdMap[1] = RDGeom::Point2D(0.19, 4.22);
  crdMap[2] = RDGeom::Point2D(1.66, 3.88);
  crdMap[3] = RDGeom::Point2D(2.10, 2.45);
  crdMap[4] = RDGeom::Point2D(1.08, 1.35);
  crdMap[5] = RDGeom::Point2D(2.56, 1.12);
  crdMap[6] = RDGeom::Point2D(1.73, 0.00);
  crdMap[7] = RDGeom::Point2D(1.08, -1.35);
  crdMap[8] = RDGeom::Point2D(-0.38, -1.69);
  crdMap[9] = RDGeom::Point2D(-1.56, -0.75);
  crdMap[10] = RDGeom::Point2D(-1.56, 0.75);
  crdMap[11] = RDGeom::Point2D(-0.38, 1.69);

  delete mref;
  mref = SmilesToMol(smi, 0, 1);
  cid2 = RDDepict::compute2DCoords(*mref, nullptr, false);
  crdMap.erase(crdMap.find(5));
  // MolToMolFile(mref, "junk1.mol");
  // std::cerr << MolToXYZBlock(*mref) << std::endl;
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
    pt1.x = (double)(rand() % 1000);
    pt1.y = (double)(rand() % 1000);
    pt1.z = (double)(rand() % 1000);
    pt1.normalize();
    pt2.x = (double)(rand() % 1000);
    pt2.y = (double)(rand() % 1000);
    pt2.z = (double)(rand() % 1000);
    pt2.normalize();

    cosT = -pt1.dotProduct(pt2);
    sinT = sqrt(1.0 - cosT * cosT);
    rotnAxis = pt1.crossProduct(pt2);
    rotnAxis.normalize();
    trans.setToIdentity();
    trans.SetRotation(cosT, sinT, rotnAxis);
    pt3 = pt2;
    trans.TransformPoint(pt3);
    dt = pt1.dotProduct(pt3);
    if (fabs(dt + 1.0) > 1.0e-3) {
      BOOST_LOG(rdInfoLog) << i << " (" << pt1 << ") (" << pt2 << ") (" << pt3
                           << ") (" << rotnAxis << ") " << dt << "\n";
    }
  }
}

void testIssue248() {
  std::string smiString =
      "OCC(=O)C1(O)[CH](O)C[CH]2[CH]3CCC4=CC(=O)C=CC4(C)C3(F)[CH](O)CC21C \
   CC(C)[CH](OC(C)=O)c1nc2ccccc2c(=O)n1[CH]1CC2(OC1=O)c1c(cccc1)N1C(=O)C(C)(C)N(O)[CH]12 \
   CC(CC(=O)OCC(F)(F)C(F)(F)C(F)(F)C(F)F)CC(=O)OCC(F)(F)C(F)(F)C(F)(F)C(F)F \
   CC(O[CH]1C(=O)C2(C)[CH](O)C[CH]3OCC3(OC(C)=O)[CH]2[CH](OC(=O)c2ccccc2)C2(O)C[CH](OC([CH](O)[CH](NC(=O)c3ccccc3)c3ccccc3)=O)C(C)=C1C2(C)C)=O \
  COc1c2O[CH]3[CH](O)C=C[CH]4C33CCN(C)[CH]4Cc(c32)cc1 \
  CN1CCC23[CH]4Oc5c(O)ccc(c25)C[CH]1[CH]3C=C[CH]4O \
  CN1C2CCC1C[CH](OC(=O)C(CO)c1ccccc1)C2 \
  COC([CH]1C2CCC(C[CH]1OC(c1ccccc1)=O)N2C)=O \
  OCC(=O)C1(O)[CH](O)C[CH]2[CH]3CCC4=CC(=O)C=CC4(C)C3(F)[CH](O)CC21C";
  // COC12C3Oc4c(O)ccc5c4C33CCN(CC4CC4)[CH](C5)C3(C[CH]1C(O)(C)C(C)(C)C)CC2";

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString, spaceSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string smi = *token;
    RWMol *m = SmilesToMol(smi, 0, 1);
    unsigned int confId = RDDepict::compute2DCoords(*m);
    // check that there are no collisions in the molecules
    int natms = m->getNumAtoms();
    int i, j;
    for (i = 0; i < natms; i++) {
      const Conformer &conf = m->getConformer(confId);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      for (j = i + 1; j < natms; j++) {
        RDGeom::Point3D locj = conf.getAtomPos(j);
        locj -= loci;
        if (locj.length() <= 0.30) {
          std::cout << "mismatch: " << i << " " << j << " " << locj.length()
                    << std::endl;
          std::cout << "\t" << smi << std::endl;
          std::cout << MolToMolBlock(*m, true, confId) << std::endl;
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
  tokenizer tokens(smaString, spaceSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string sma = *token;
    RWMol *m = SmartsToMol(sma);
    unsigned int confId = RDDepict::compute2DCoords(*m);
    // check that there are no collisions in the molecules
    int natms = m->getNumAtoms();
    int i, j;
    for (i = 0; i < natms; i++) {
      const Conformer &conf = m->getConformer(confId);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      for (j = i + 1; j < natms; j++) {
        RDGeom::Point3D locj = conf.getAtomPos(j);
        locj -= loci;
        if (locj.length() <= 0.30) {
          std::cout << "mismatch: " << i << " " << j << " " << locj.length()
                    << std::endl;
          std::cout << "\t" << sma << std::endl;
          std::cout << MolToMolBlock(*m, true, confId) << std::endl;
        }
        CHECK_INVARIANT(locj.length() > 0.30, "");
      }
    }
    delete m;
  }
}

void testRemoveHsCrash() {
  std::string rdbase = getenv("RDBASE");
  std::string molfile =
      rdbase + "/Code/GraphMol/Depictor/test_data/hs_crash.mol";
  RWMol *m = MolFileToMol(molfile, true, false);
  TEST_ASSERT(m);
  ROMol *newM = MolOps::removeHs(*static_cast<ROMol *>(m));
  delete m;
  RDDepict::compute2DCoords(*newM);
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
  crdMap[0] = RDGeom::Point2D(0.0781827, -0.112109);
  crdMap[1] = RDGeom::Point2D(1.13205, -0.415657);
  crdMap[2] = RDGeom::Point2D(0.676323, -1.51974);
  crdMap[3] = RDGeom::Point2D(0.735363, -2.62024);
  crdMap[4] = RDGeom::Point2D(0.974216, -1.5726);
  crdMap[5] = RDGeom::Point2D(1.02041, -1.56945);
  std::string mb =
      "\n\
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
  (void)cid1;
  // TEST_ASSERT(cid1>=0);
  delete m1;
}

void testIssue2821647() {
  {
    std::string smi = "CCCCC";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);
    double xx = 0, yy = 0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx += loci.x * loci.x;
      yy += loci.y * loci.y;
    }
    TEST_ASSERT(xx > yy);
    delete m1;
  }
  {
    std::string smi = "c1ccccc1CCCCCC1CC1";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);
    double xx = 0, yy = 0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx += loci.x * loci.x;
      yy += loci.y * loci.y;
    }
    TEST_ASSERT(xx > yy);
    delete m1;
  }
  {
    std::string smi = "c1ccc2c(c1)oc1c3ccccc3oc21";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);
    double xx = 0, yy = 0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx += loci.x * loci.x;
      yy += loci.y * loci.y;
    }
    TEST_ASSERT(xx > yy);
    delete m1;
  }
  {
    std::string smi = "[H]n1c2ccccc2c2n([H])c3ccccc3c12";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);
    double xx = 0, yy = 0;
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      const Conformer &conf = m1->getConformer(cid1);
      RDGeom::Point3D loci = conf.getAtomPos(i);
      xx += loci.x * loci.x;
      yy += loci.y * loci.y;
    }
    TEST_ASSERT(xx > yy);
    delete m1;
  }
}

void testIssue2948402() {
  {
    std::string smi = "C1C2CC3=CC=CC(C2)CC(O3)C1";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);
    (void)cid1;
    // TEST_ASSERT(cid1>=0);
    delete m1;
  }
}

void testIssue2995724() {
  {
    // the original problem from Thomas Heller:
    std::string smi = "OC(=O)[C@@H]1CCCN1";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);

    const Conformer &conf = m1->getConformer(cid1);
    for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
      RDGeom::Point3D loci = conf.getAtomPos(i);
      // we're testing for NaNs here:
      TEST_ASSERT(loci.x > -5.0);
      TEST_ASSERT(loci.x < 5.0);
      TEST_ASSERT(loci.y > -5.0);
      TEST_ASSERT(loci.y < 5.0);
    }
    delete m1;
  }
  {
    // additional test cases from Kirk DeLisle:
    std::string smis[] = {"CCC(N1CCN(CC1)C(CC)O)O.Cl[Pt](Cl)(Cl)(Cl)(Cl)Cl",
                          "CN(C)S(=O)(=O)N1CCN(CC1)S(=O)(=O)N(C)C",
                          "Cc1ccc(cc1C)Nc2nc(nc(n2)N)CN3CCN(CC3)C",
                          "CC1(OC(=C(C(=O)O1)C2=NCCC2)O)C"};
    for (const auto &smi : smis) {
      RWMol *m1 = SmilesToMol(smi);
      unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);

      const Conformer &conf = m1->getConformer(cid1);
      for (unsigned int i = 0; i < m1->getNumAtoms(); i++) {
        RDGeom::Point3D loci = conf.getAtomPos(i);
        TEST_ASSERT(loci.x > -7.0);
        TEST_ASSERT(loci.x < 7.0);
        TEST_ASSERT(loci.y > -6.0);
        TEST_ASSERT(loci.y < 6.0);
      }
      delete m1;
    }
  }
}

void testBondLengthChange() {
  {
    std::string smi = "CC";
    RWMol *m1 = SmilesToMol(smi);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);

    const Conformer &conf = m1->getConformer(cid1);
    TEST_ASSERT(feq(conf.getAtomPos(0).x, -0.75));
    TEST_ASSERT(feq(conf.getAtomPos(1).x, 0.75));

    delete m1;
  }
  {
    std::string smi = "CC";
    RWMol *m1 = SmilesToMol(smi);
    RDDepict::BOND_LEN = 1.0;
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);

    const Conformer &conf = m1->getConformer(cid1);
    TEST_ASSERT(feq(conf.getAtomPos(0).x, -0.5));
    TEST_ASSERT(feq(conf.getAtomPos(1).x, 0.5));

    delete m1;
  }
}

void testIssue3122141() {
  {
    std::string smi = "ClC1=C(Cl)C2(Cl)C3C(Cl)C=CC3C1(Cl)C2(Cl)Cl";
    RWMol *m1 = SmilesToMol(smi);
    RDGeom::INT_POINT2D_MAP crdMap;
    crdMap[1] = RDGeom::Point2D(1.0, 2.0);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
    (void)cid1;
    // TEST_ASSERT(cid1>=0);
    delete m1;
  }
}

void testIssue3135833() {
  {
    std::string smi =
        "N#CC(/C=[N+](/c1cccnc1)c1c2cc(C(=O)c3ncsc3-c3cccnc3)ccc2ccc1)=C(/"
        "[NH3+])[S-]";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);

    RDGeom::INT_POINT2D_MAP crdMap;
    crdMap[32] = RDGeom::Point2D(-2.8518, 1.0270);
    crdMap[33] = RDGeom::Point2D(-2.7999, -0.4721);
    crdMap[4] = RDGeom::Point2D(-1.4239, -2.6758);
    crdMap[11] = RDGeom::Point2D(-1.4757, -1.1767);
    crdMap[12] = RDGeom::Point2D(-0.2034, -0.3823);
    crdMap[13] = RDGeom::Point2D(1.1208, -1.0869);
    crdMap[14] = RDGeom::Point2D(2.3931, -0.2924);
    crdMap[15] = RDGeom::Point2D(3.7173, -0.9971);
    crdMap[28] = RDGeom::Point2D(2.3412, 1.2067);
    crdMap[29] = RDGeom::Point2D(1.0170, 1.9113);
    crdMap[30] = RDGeom::Point2D(-0.2553, 1.1168);
    crdMap[31] = RDGeom::Point2D(-1.5795, 1.8215);

    unsigned int cid1 = RDDepict::compute2DCoords(*m1, &crdMap, false);
    (void)cid1;
    // TEST_ASSERT(cid1>=0);
    delete m1;
  }
}

void testIssue3487469() {
  {
    std::string smi = "C*C";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getAtomWithIdx(1)->getHybridization() == Atom::UNSPECIFIED);
    unsigned int cid1 = RDDepict::compute2DCoords(*m1, nullptr, true);
    const Conformer &conf = m1->getConformer(cid1);
    RDGeom::Point3D p0 = conf.getAtomPos(0);
    RDGeom::Point3D p1 = conf.getAtomPos(1);
    RDGeom::Point3D p2 = conf.getAtomPos(2);

    RDGeom::Point3D v1 = p0 - p1, v2 = p2 - p1;
    v1.normalize();
    v2.normalize();
    TEST_ASSERT(feq(v1.dotProduct(v2), -0.5, .01))
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
    ROMol *m2 = deleteSubstructs(*m1, *p);
    TEST_ASSERT(m2);
    TEST_ASSERT(m2->getNumAtoms() == m1->getNumAtoms() - 1);
    unsigned int cid1 = RDDepict::compute2DCoords(*m2);
    TEST_ASSERT(cid1 == 0);
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

    RDGeom::Point3D p87 = p8 - p7;
    RDGeom::Point3D p86 = p8 - p6;
    RDGeom::Point3D p89 = p8 - p9;
    TEST_ASSERT(p87.dotProduct(p86) * p87.dotProduct(p89) < 0);

    RDGeom::Point3D p67 = p6 - p7;
    RDGeom::Point3D p68 = p6 - p8;
    RDGeom::Point3D p65 = p6 - p5;
    TEST_ASSERT(p67.dotProduct(p68) * p67.dotProduct(p65) < 0);

    delete m1;
  }
  {  // a collection of previous failures
    std::string smis[4] = {"C1=CC=C2C(=C1)C3=CC=CC=C3C4=C2C(C(C5C4O5)O)O",
                           "CC1=C2C=C(C=CC2=NC3=C1C4=C(O4)C5=CC=CC=C53)CO",
                           "CC1=C2C=CC3=C(C2=CC4=CC=CC=C14)C5C(O5)C(C3O)O",
                           "C1=CC=C2C(=C1)C=C3C=CC4=C5C3=C2C6C(C5=CC=C4)O6"};
    RWMol *p = SmartsToMol("[#6]~[#6]~1-[#8]-[#6]~1~[#6]");
    TEST_ASSERT(p);
    for (const auto &smi : smis) {
      RWMol *m = SmilesToMol(smi);
      TEST_ASSERT(m);
      MatchVectType mv;
      TEST_ASSERT(SubstructMatch(*m, *p, mv));
      TEST_ASSERT(mv.size() == 5);

      unsigned int cid1 = RDDepict::compute2DCoords(*m);
      const Conformer &conf = m->getConformer(cid1);
      RDGeom::Point3D v10 =
          conf.getAtomPos(mv[1].second) - conf.getAtomPos(mv[0].second);
      RDGeom::Point3D v12 =
          conf.getAtomPos(mv[1].second) - conf.getAtomPos(mv[2].second);
      RDGeom::Point3D v13 =
          conf.getAtomPos(mv[1].second) - conf.getAtomPos(mv[3].second);
      TEST_ASSERT(v12.dotProduct(v10) * v12.dotProduct(v13) < 0);

      RDGeom::Point3D v31 =
          conf.getAtomPos(mv[3].second) - conf.getAtomPos(mv[1].second);
      RDGeom::Point3D v32 =
          conf.getAtomPos(mv[3].second) - conf.getAtomPos(mv[2].second);
      RDGeom::Point3D v34 =
          conf.getAtomPos(mv[3].second) - conf.getAtomPos(mv[4].second);
      TEST_ASSERT(v32.dotProduct(v32) * v32.dotProduct(v34) < 0);

      delete m;
    }
    delete p;
  }
}

void testGitHubIssue910() {
  {
    // this is a ChEMBL molecule
    std::string smiles =
        "CSCC[C@H](NC(=O)[C@@H](CCC(N)=O)NC(=O)[C@@H](N)Cc1c[nH]c2ccccc12)C(=O)"
        "NCC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](CO)C(=O)O";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);

    // add chiral Hs, they were part of the problem
    std::vector<unsigned int> chiralAts;
    for (RWMol::AtomIterator atIt = m->beginAtoms(); atIt != m->endAtoms();
         ++atIt) {
      if ((*atIt)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
          (*atIt)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
        chiralAts.push_back((*atIt)->getIdx());
      }
    }
    MolOps::addHs(*m, false, false, &chiralAts);
    RDDepict::compute2DCoords(*m, nullptr, true);
#if 0
    m->setProp("_Name", "github910");
    std::cerr << MolToMolBlock(*m);
#endif
    // now look for close contacts.
    const Conformer &conf = m->getConformer();
    for (unsigned int i = 0; i < conf.getNumAtoms(); ++i) {
      for (unsigned int j = i + 1; j < conf.getNumAtoms(); ++j) {
        double l = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        TEST_ASSERT(l > 0.75);
      }
    }

    delete m;
  }
}

void testGitHubIssue1073() {
  // computeInitialCoords() should call the SSSR code before it calls
  // assignStereochemistry()
  {
    std::string smarts = "[a]12[a][a][a][a][a]1[a][a][a]2";
    RWMol *m = SmartsToMol(smarts);
    TEST_ASSERT(m);

    // compute2DCoords does ring finding internally, so there's no error
    RDDepict::compute2DCoords(*m);

    // but the molecule itself is not modified
    RingInfo *ri = m->getRingInfo();
    TEST_ASSERT(!ri->isInitialized());

    delete m;
  }
}

void testConstrainedCoords() {
  std::string rdbase = getenv("RDBASE");
  std::string ofile =
      rdbase + "/Code/GraphMol/Depictor/test_data/constrainedCoords.out.sdf";
  SDWriter writer(ofile);

  std::string templ_smiles = "c1nccc2n1ccc2";
  ROMol *templ = SmilesToMol(templ_smiles);
  TEST_ASSERT(templ);
  RDDepict::compute2DCoords(*templ);
  std::string smiles = "c1cccc2ncn3cccc3c21";
  ROMol *m = SmilesToMol(smiles);
  TEST_ASSERT(m);
  RDDepict::generateDepictionMatching2DStructure(*m, *templ);
  writer.write(*m);

  std::string smarts = "*1****2*1***2";
  ROMol *refPatt = SmartsToMol(smarts);
  RDDepict::generateDepictionMatching2DStructure(*m, *templ, -1, refPatt);
  writer.write(*m);

  delete templ;
  delete m;
  delete refPatt;

  std::string xp0_file =
      rdbase + "/Code/GraphMol/Depictor/test_data/1XP0_ligand.sdf";
  RDKit::ROMol *xp0_lig = RDKit::MolFileToMol(xp0_file);
  auto *xp0_lig_2d = new RDKit::ROMol(*xp0_lig);
  RDDepict::compute2DCoords(*xp0_lig_2d);
  writer.write(*xp0_lig_2d);
  RDDepict::generateDepictionMatching3DStructure(*xp0_lig_2d, *xp0_lig);
  writer.write(*xp0_lig_2d);

  delete xp0_lig;
  delete xp0_lig_2d;
}
void testGitHubIssue1112() {
  // Bad coordinate generation for H2
  {
    std::string smiles = "[H][H]";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);

    RDDepict::compute2DCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).x, 0));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).y, 0));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(0).z, 0));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(1).x, 0));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(1).y, -1));
    TEST_ASSERT(feq(m->getConformer().getAtomPos(1).z, 0));

    delete m;
  }
}

void testGitHubIssue1286() {
  // GenerateDepictionMatching2DStructure isn't matching 2D structure
  {  // the original report
    std::string smiles = "C(=O)C(C)NC=O";
    RWMol *templ = SmilesToMol(smiles);
    TEST_ASSERT(templ);
    TEST_ASSERT(templ->getNumAtoms() == 7);

    smiles = "C(=O)C(C)NC(=O)C1CC1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 10);

    RDDepict::compute2DCoords(*templ);
    TEST_ASSERT(templ->getNumConformers() == 1);
    RDDepict::generateDepictionMatching2DStructure(*mol, *templ);
    TEST_ASSERT(mol->getNumConformers() == 1);

    // std::cout << MolToMolBlock(*templ) << std::endl;
    // std::cout << MolToMolBlock(*mol) << std::endl;

    const Conformer &tconf = templ->getConformer();
    const Conformer &mconf = mol->getConformer();
    for (unsigned int i = 0; i < templ->getNumAtoms(); ++i) {
      const RDGeom::Point3D &tp = tconf.getAtomPos(i);
      const RDGeom::Point3D &mp = mconf.getAtomPos(i);
      // std::cerr << i << ": " << tp << " | " << mp << std::endl;
      TEST_ASSERT(feq(tp.x, mp.x));
      TEST_ASSERT(feq(tp.y, mp.y));
    }

    delete templ;
    delete mol;
  }
  {  // extremely crowded. This one tests bond shortening and angle opening
    std::string smiles = "CC(=O)C1=CC=CC2=C1C=CC=C2";
    RWMol *templ = SmilesToMol(smiles);
    TEST_ASSERT(templ);
    TEST_ASSERT(templ->getNumAtoms() == 13);

    smiles = "O=C(N)C1=C(C=CC2=C1C(=CC=C2)C(C)=O)C(C)(C)C";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 20);

    RDDepict::compute2DCoords(*templ);
    TEST_ASSERT(templ->getNumConformers() == 1);
    RDDepict::generateDepictionMatching2DStructure(*mol, *templ);
    TEST_ASSERT(mol->getNumConformers() == 1);

    // std::cout << MolToMolBlock(*templ) << std::endl;
    // std::cout << MolToMolBlock(*mol) << std::endl;

    MatchVectType mv;
    TEST_ASSERT(SubstructMatch(*mol, *templ, mv));

    const Conformer &tconf = templ->getConformer();
    const Conformer &mconf = mol->getConformer();
    for (unsigned int i = 0; i < templ->getNumAtoms(); ++i) {
      const RDGeom::Point3D &tp = tconf.getAtomPos(mv[i].first);
      const RDGeom::Point3D &mp = mconf.getAtomPos(mv[i].second);
      // std::cerr << i << ": " << tp << " | " << mp << std::endl;
      TEST_ASSERT(feq(tp.x, mp.x));
      TEST_ASSERT(feq(tp.y, mp.y));
    }

    delete templ;
    delete mol;
  }
}

void testGithub1691() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Github issue "
         "1691: Acetylenic hydrogens not given appropriate 2D coordinates"
      << std::endl;
#if 1
  {
    SmilesParserParams ps;
    ps.removeHs = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("C1#C2.[F]1.[F]2", ps));

    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(mol->getBondBetweenAtoms(1, 3));
    RDDepict::compute2DCoords(*mol);

    // std::cerr << MolToMolBlock(*mol) << std::endl;
    const Conformer &conf = mol->getConformer();
    RDGeom::Point3D v20 = conf.getAtomPos(2) - conf.getAtomPos(0);
    RDGeom::Point3D v10 = conf.getAtomPos(1) - conf.getAtomPos(0);
    RDGeom::Point3D v31 = conf.getAtomPos(3) - conf.getAtomPos(1);
    RDGeom::Point3D v01 = conf.getAtomPos(0) - conf.getAtomPos(1);
    // std::cerr << v20.dotProduct(v10) << std::endl;
    // std::cerr << v31.dotProduct(v01) << std::endl;
    TEST_ASSERT(v20.dotProduct(v10) <= -1.0);
    TEST_ASSERT(v31.dotProduct(v01) <= -1.0);
  }
#endif
  {
    SmilesParserParams ps;
    ps.removeHs = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("C1#C2.[H]1.[H]2", ps));

    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(mol->getBondBetweenAtoms(1, 3));
    RDDepict::compute2DCoords(*mol);

    // std::cerr << MolToMolBlock(*mol) << std::endl;
    const Conformer &conf = mol->getConformer();
    RDGeom::Point3D v20 = conf.getAtomPos(2) - conf.getAtomPos(0);
    RDGeom::Point3D v10 = conf.getAtomPos(1) - conf.getAtomPos(0);
    RDGeom::Point3D v31 = conf.getAtomPos(3) - conf.getAtomPos(1);
    RDGeom::Point3D v01 = conf.getAtomPos(0) - conf.getAtomPos(1);
    // std::cerr << v20.dotProduct(v10) << std::endl;
    // std::cerr << v31.dotProduct(v01) << std::endl;
    TEST_ASSERT(v20.dotProduct(v10) <= -1.0);
    TEST_ASSERT(v31.dotProduct(v01) <= -1.0);
  }
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("C#C"));

    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getNumAtoms() == 2);
    MolOps::addHs(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(mol->getBondBetweenAtoms(1, 3));
    RDDepict::compute2DCoords(*mol);

    // std::cerr << MolToMolBlock(*mol) << std::endl;
    const Conformer &conf = mol->getConformer();
    RDGeom::Point3D v20 = conf.getAtomPos(2) - conf.getAtomPos(0);
    RDGeom::Point3D v10 = conf.getAtomPos(1) - conf.getAtomPos(0);
    RDGeom::Point3D v31 = conf.getAtomPos(3) - conf.getAtomPos(1);
    RDGeom::Point3D v01 = conf.getAtomPos(0) - conf.getAtomPos(1);
    // std::cerr << v20.dotProduct(v10) << std::endl;
    // std::cerr << v31.dotProduct(v01) << std::endl;
    TEST_ASSERT(v20.dotProduct(v10) <= -1.0);
    TEST_ASSERT(v31.dotProduct(v01) <= -1.0);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub2027() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github issue "
                          "2027: \"linear\" fragments not canonically oriented"
                       << std::endl;
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("C#CC#CC#CC#CC#CC#CC#C"));
    TEST_ASSERT(mol);
    RDDepict::compute2DCoords(*mol, nullptr, true);
    const Conformer &conf = mol->getConformer();
    TEST_ASSERT(feq(conf.getAtomPos(0).y, 0.0));
    TEST_ASSERT(feq(conf.getAtomPos(1).y, 0.0));
    TEST_ASSERT(feq(conf.getAtomPos(2).y, 0.0));
  }
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=CC2=CC3=CC=CC=C3C=C12"));
    TEST_ASSERT(mol);
    RDDepict::compute2DCoords(*mol, nullptr, true);

    // a stupidly simple test to ensure that we're oriented along the x axis:
    const Conformer &conf = mol->getConformer();
    RDGeom::Point2D paccum(0, 0);
    for (const auto &pt : conf.getPositions()) {
      paccum.x += fabs(pt.x);
      paccum.y += fabs(pt.y);
    }
    TEST_ASSERT(paccum.x > paccum.y);
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGenerate2DDepictionRefPatternMatchVect() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Test "
         "generateDepictionMatching2DStructure with refPattern and matchVect"
      << std::endl;
  auto indazoleRef = R"RES(
     RDKit          2D

  9 10  0  0  0  0  0  0  0  0999 V2000
   -6.0878    2.4335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3867    1.6835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3867    0.1833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0878   -0.5666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7887    0.1833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7887    1.6835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4897   -0.5664    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1906    1.6833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1906    0.1835    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  8  9  2  0
  6  8  1  0
  7  9  1  0
  7  5  1  0
M  END)RES"_ctab;
  auto cycloheptylPyrazole = "c1cc(C2CCCCCC2)[nH]n1"_smiles;
  double msd;
  bool raised;

  // test using refPattern
  auto refPatt = "a1aan[nH]1"_smarts;
  RDDepict::generateDepictionMatching2DStructure(
      *cycloheptylPyrazole, *indazoleRef, -1, refPatt.get());
  TEST_ASSERT(cycloheptylPyrazole->getNumConformers() == 1);
  MatchVectType molMatchVect;
  TEST_ASSERT(SubstructMatch(*cycloheptylPyrazole, *refPatt, molMatchVect));
  MatchVectType refMatchVect;
  TEST_ASSERT(SubstructMatch(*indazoleRef, *refPatt, refMatchVect));
  TEST_ASSERT(molMatchVect.size() == refMatchVect.size());
  msd = 0.0;
  for (size_t i = 0; i < molMatchVect.size(); ++i) {
    msd += (indazoleRef->getConformer().getAtomPos(refMatchVect.at(i).second) -
            cycloheptylPyrazole->getConformer().getAtomPos(
                molMatchVect.at(i).second))
               .lengthSq();
  }
  msd /= static_cast<double>(molMatchVect.size());
  TEST_ASSERT(msd < 1.0e-4);
  // try with a pattern larger than the reference molecule
  auto hugePatt = "CCCCCCCCCCCCCCCCCCCCCCCCCCC"_smarts;
  raised = false;
  try {
    RDDepict::generateDepictionMatching2DStructure(
        *cycloheptylPyrazole, *indazoleRef, -1, hugePatt.get());
  } catch (const RDDepict::DepictException &) {
    raised = true;
  }
  TEST_ASSERT(raised);
  // try with an out of range confId
  raised = false;
  try {
    RDDepict::generateDepictionMatching2DStructure(
        *cycloheptylPyrazole, *indazoleRef, 1, refPatt.get());
  } catch (const RDKit::ConformerException &) {
    raised = true;
  }

  // test using matchVect directly
  cycloheptylPyrazole->removeConformer(0);
  MatchVectType matchVect;
  for (size_t i = 0; i < molMatchVect.size(); ++i) {
    matchVect.emplace_back(
        std::make_pair(refMatchVect.at(i).second, molMatchVect.at(i).second));
  }
  RDDepict::generateDepictionMatching2DStructure(*cycloheptylPyrazole,
                                                 *indazoleRef, matchVect);
  TEST_ASSERT(cycloheptylPyrazole->getNumConformers() == 1);
  msd = 0.0;
  for (const auto &pair : matchVect) {
    msd += (indazoleRef->getConformer().getAtomPos(pair.first) -
            cycloheptylPyrazole->getConformer().getAtomPos(pair.second))
               .lengthSq();
  }
  msd /= static_cast<double>(matchVect.size());
  TEST_ASSERT(msd < 1.0e-4);
  // try with a matchVect larger than the reference molecule
  MatchVectType matchVectHuge(matchVect);
  for (size_t i = 0; i < indazoleRef->getNumAtoms(); ++i) {
    matchVectHuge.emplace_back(std::make_pair(0, 0));
  }
  raised = false;
  try {
    RDDepict::generateDepictionMatching2DStructure(*cycloheptylPyrazole,
                                                   *indazoleRef, matchVectHuge);
  } catch (const RDDepict::DepictException &) {
    raised = true;
  }
  // try with a matchVect with out of range indices
  MatchVectType matchVectOutOfRange(matchVect);
  matchVectOutOfRange.emplace_back(std::make_pair(100, 100));
  raised = false;
  try {
    RDDepict::generateDepictionMatching2DStructure(
        *cycloheptylPyrazole, *indazoleRef, matchVectOutOfRange);
  } catch (const RDDepict::DepictException &) {
    raised = true;
  }
  // try with an out of range confId
  raised = false;
  try {
    RDDepict::generateDepictionMatching2DStructure(*cycloheptylPyrazole,
                                                   *indazoleRef, matchVect, 1);
  } catch (const RDKit::ConformerException &) {
    raised = true;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGenerate2DDepictionAllowRGroupsOrig() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Test "
         "generateDepictionMatching2DStructure with allowRGroups"
      << std::endl;
  auto templateRef = R"RES(
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END)RES"_ctab;
  auto orthoMeta = "c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1"_smiles;
  auto ortho = "c1ccc(-c2ccccc2-c2ccccc2)cc1"_smiles;
  auto meta = "c1ccc(-c2cccc(-c3ccccc3)c2)cc1"_smiles;
  auto biphenyl = "c1ccccc1-c1ccccc1"_smiles;
  auto phenyl = "c1ccccc1"_smiles;

  RDDepict::generateDepictionMatching2DStructure(*orthoMeta, *templateRef);
  TEST_ASSERT(orthoMeta->getNumConformers() == 1);

  for (auto mol : {ortho.get(), meta.get(), biphenyl.get(), phenyl.get()}) {
    // fails as does not match template
    bool raised = false;
    try {
      RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef);
    } catch (const RDDepict::DepictException &) {
      raised = true;
    }
    TEST_ASSERT(raised);

    // succeeds with allowRGroups = true
    auto matchVect = RDDepict::generateDepictionMatching2DStructure(
        *mol, *templateRef, -1, nullptr, false, false, true);
    TEST_ASSERT(mol->getNumConformers() == 1);
    double msd = 0.0;
    for (const auto &pair : matchVect) {
      msd += (templateRef->getConformer().getAtomPos(pair.first) -
              mol->getConformer().getAtomPos(pair.second))
                 .lengthSq();
    }
    msd /= static_cast<double>(matchVect.size());
    TEST_ASSERT(msd < 1.0e-4);
  }

  // test that using a refPattern with R groups and a reference without works
  auto pyridineRef = R"RES(
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END)RES"_ctab;
  auto genericRefPatternWithRGroups = "[*:3]a1a([*:1])aa([*:2])aa1"_smarts;
  std::unique_ptr<ROMol> pyridineRefHs(
      MolOps::addHs(static_cast<const ROMol &>(*pyridineRef)));

  for (auto mol : {ortho.get(), meta.get(), biphenyl.get(), phenyl.get()}) {
    auto matchVect = RDDepict::generateDepictionMatching2DStructure(
        *mol, *pyridineRef, -1, genericRefPatternWithRGroups.get(), false,
        false, true);
    TEST_ASSERT(mol->getNumConformers() == 1);
    double msd = 0.0;
    for (const auto &pair : matchVect) {
      msd += (pyridineRef->getConformer().getAtomPos(pair.first) -
              mol->getConformer().getAtomPos(pair.second))
                 .lengthSq();
    }
    msd /= static_cast<double>(matchVect.size());
    TEST_ASSERT(msd < 1.0e-4);
  }

  // test that using a reference with query atoms including H works
  auto scaffold = R"CTAB(
  MJ201100                      

 12 13  0  0  0  0  0  0  0  0999 V2000
   -0.5398    0.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3648    0.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7773   -0.6745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3649   -1.3889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5399   -1.3889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1273   -0.6744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6976   -0.6744    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9167    0.6531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6704    0.3176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5842   -0.5028    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3849    0.7302    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7451    1.4600    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
  8  9  2  0  0  0  0
  2  8  1  0  0  0  0
  9 10  1  0  0  0  0
  3 10  1  0  0  0  0
  9 11  1  0  0  0  0
  8 12  1  0  0  0  0
M  ALS   7 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  ALS  11 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  ALS  12 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  END
)CTAB"_ctab;
  auto mol = R"CTAB(
  MJ201100                      

 13 14  0  0  0  0  0  0  0  0999 V2000
   -0.6112    0.3665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3648    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4510   -0.7895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7836   -1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0299   -0.9389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0562   -0.1183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8099    0.2172    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1184    0.3666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6705   -0.2464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2580   -0.9608    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.6374   -1.4238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8961    1.0377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5512   -2.2443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  8  9  2  0  0  0  0
  2  8  1  0  0  0  0
  9 10  1  0  0  0  0
  3 10  1  0  0  0  0
  6  7  1  0  0  0  0
  5 11  1  0  0  0  0
  7 12  1  0  0  0  0
 11 13  1  0  0  0  0
M  END
)CTAB"_ctab;
  auto matchVect = RDDepict::generateDepictionMatching2DStructure(
      *mol, *scaffold, -1, nullptr, false, false, true);
  TEST_ASSERT(mol->getNumConformers() == 1);
  TEST_ASSERT(matchVect.size() == 10);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGenerate2DDepictionAllowRGroups() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Test "
         "generateDepictionMatching2DStructure with allowRGroups"
      << std::endl;
  auto templateRef = R"RES(
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END)RES"_ctab;
  TEST_ASSERT(templateRef);
  auto orthoMeta = "c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1"_smiles;
  auto ortho = "c1ccc(-c2ccccc2-c2ccccc2)cc1"_smiles;
  auto meta = "c1ccc(-c2cccc(-c3ccccc3)c2)cc1"_smiles;
  auto para = "c1ccc(-c2ccc(-c3ccccc3)cc2)cc1"_smiles;
  auto biphenyl = "c1ccccc1-c1ccccc1"_smiles;
  auto phenyl = "c1ccccc1"_smiles;

  auto prevBondLen = RDDepict::BOND_LEN;
  RDDepict::BOND_LEN = defaultRDKitBondLen;
  RDDepict::generateDepictionMatching2DStructure(*orthoMeta, *templateRef);
  TEST_ASSERT(orthoMeta->getNumConformers() == 1);
  for (bool alignOnly : {true, false}) {
    for (auto mol :
         {ortho.get(), meta.get(), para.get(), biphenyl.get(), phenyl.get()}) {
      TEST_ASSERT(mol);
      RDDepict::ConstrainedDepictionParams p;
      p.allowRGroups = true;
      p.alignOnly = alignOnly;
      // fails as does not match template
      bool raised = false;
      try {
        RDDepict::generateDepictionMatching2DStructure(*mol, *templateRef);
      } catch (const RDDepict::DepictException &) {
        raised = true;
      }
      TEST_ASSERT(raised);

      // succeeds with allowRGroups = true
      auto matchVect = RDDepict::generateDepictionMatching2DStructure(
          *mol, *templateRef, -1, nullptr, p);
      TEST_ASSERT(!matchVect.empty());
      TEST_ASSERT(mol->getNumConformers() == 1);
      double msd = 0.0;
      for (const auto &pair : matchVect) {
        msd += (templateRef->getConformer().getAtomPos(pair.first) -
                mol->getConformer().getAtomPos(pair.second))
                   .lengthSq();
      }
      msd /= static_cast<double>(matchVect.size());
      TEST_ASSERT(msd < 1.0e-4);
    }

    // test that using a refPattern with R groups and a reference missing one
    // works
    auto pyridineRef = R"RES(
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.5469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3395    0.7734    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3395   -0.7732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5469    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3395   -0.7732    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3395    0.7734    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -3.0938    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  4  8  1  0
M  END)RES"_ctab;
    TEST_ASSERT(pyridineRef);
    auto genericRefPatternWithRGroups = "[*:3]a1a([*:1])aa([*:2])aa1"_smarts;
    TEST_ASSERT(genericRefPatternWithRGroups);
    for (auto [numExpectedMatches, mol] :
         std::vector<std::pair<unsigned int, ROMol *>>{{8, orthoMeta.get()},
                                                       {7, ortho.get()},
                                                       {7, meta.get()},
                                                       {8, para.get()},
                                                       {7, biphenyl.get()},
                                                       {6, phenyl.get()}}) {
      RDDepict::ConstrainedDepictionParams p;
      p.allowRGroups = true;
      p.alignOnly = alignOnly;
      auto matchVect = RDDepict::generateDepictionMatching2DStructure(
          *mol, *pyridineRef, -1, genericRefPatternWithRGroups.get(), p);
      TEST_ASSERT(matchVect.size() == numExpectedMatches);
      TEST_ASSERT(mol->getNumConformers() == 1);
      double msd = 0.0;
      for (const auto &pair : matchVect) {
        msd += (pyridineRef->getConformer().getAtomPos(pair.first) -
                mol->getConformer().getAtomPos(pair.second))
                   .lengthSq();
      }
      msd /= static_cast<double>(matchVect.size());
      TEST_ASSERT(msd < (alignOnly ? 5.e-3 : 1.0e-4));
    }

    // test that using a reference with query atoms including H works
    auto scaffold = R"CTAB(
  MJ201100                      

 12 13  0  0  0  0  0  0  0  0999 V2000
   -0.5398    0.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3648    0.0400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7773   -0.6745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3649   -1.3889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5399   -1.3889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1273   -0.6744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6976   -0.6744    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9167    0.6531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6704    0.3176    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5842   -0.5028    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3849    0.7302    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7451    1.4600    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
  8  9  2  0  0  0  0
  2  8  1  0  0  0  0
  9 10  1  0  0  0  0
  3 10  1  0  0  0  0
  9 11  1  0  0  0  0
  8 12  1  0  0  0  0
M  ALS   7 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  ALS  11 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  ALS  12 10 F H   C   N   O   F   P   S   Cl  Br  I   
M  END
)CTAB"_ctab;
    TEST_ASSERT(scaffold);
    auto mol = R"CTAB(
  MJ201100                      

 13 14  0  0  0  0  0  0  0  0999 V2000
   -0.6112    0.3665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3648    0.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4510   -0.7895    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7836   -1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0299   -0.9389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0562   -0.1183    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8099    0.2172    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1184    0.3666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6705   -0.2464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2580   -0.9608    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.6374   -1.4238    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8961    1.0377    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5512   -2.2443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  8  9  2  0  0  0  0
  2  8  1  0  0  0  0
  9 10  1  0  0  0  0
  3 10  1  0  0  0  0
  6  7  1  0  0  0  0
  5 11  1  0  0  0  0
  7 12  1  0  0  0  0
 11 13  1  0  0  0  0
M  END
)CTAB"_ctab;
    TEST_ASSERT(mol);
    auto matchVect = RDDepict::generateDepictionMatching2DStructure(
        *mol, *scaffold, -1, nullptr, false, false, true);
    TEST_ASSERT(mol->getNumConformers() == 1);
    TEST_ASSERT(matchVect.size() == 10);
  }
  RDDepict::BOND_LEN = prevBondLen;
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testNormalizeStraighten() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Test normalize and straighten depiction"
      << std::endl;

  auto noradrenalineMJ = R"RES(
  MJ201100                      

 12 12  0  0  1  0  0  0  0  0999 V2000
    2.2687    1.0716    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0312    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437   -0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2062    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437   -1.0716    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2687    0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  2  1  0  0  0  0
  3  4  1  6  0  0  0
  3  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  2  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
 11 12  2  0  0  0  0
  5 12  1  0  0  0  0
M  END)RES"_ctab;
  {
    auto noradrenalineMJCopy =
        std::unique_ptr<RWMol>(new RWMol(*noradrenalineMJ));
    const auto &conformer0 = noradrenalineMJCopy->getConformer(0);
    auto conformer1 = new Conformer(conformer0);
    noradrenalineMJCopy->addConformer(conformer1, true);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  0) < 1.e-5);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  1) < 1.e-5);
    auto scalingFactor = RDDepict::normalizeDepiction(*noradrenalineMJCopy, 1);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  0) < 1.e-5);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  1) > 1.e-5);
    TEST_ASSERT(RDKit::feq(scalingFactor, 1.875, 1.e-3));
    auto conformer2 = new Conformer(*conformer1);
    noradrenalineMJCopy->addConformer(conformer2, true);
    auto bond10_11Conf0 = conformer0.getAtomPos(11) - conformer0.getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf0.x, 0.825, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf0.y, 0.0, 1.e-3));
    auto bond10_11Conf1 =
        conformer1->getAtomPos(11) - conformer1->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.x, 1.513, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.y, -0.321, 1.e-3));
    RDDepict::straightenDepiction(*noradrenalineMJCopy, 1);
    bond10_11Conf1 = conformer1->getAtomPos(11) - conformer1->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.x, 1.340, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.y, -0.773, 1.e-3));
    auto bond4_11Conf1 = conformer1->getAtomPos(11) - conformer1->getAtomPos(4);
    TEST_ASSERT(RDKit::feq(bond4_11Conf1.x, 0.0, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond4_11Conf1.y, 1.547, 1.e-3));
    RDDepict::straightenDepiction(*noradrenalineMJCopy, 2, true);
    auto bond10_11Conf2 =
        conformer2->getAtomPos(11) - conformer2->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf2.x, 1.547, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf2.y, 0.0, 1.e-3));
    auto bond4_11Conf2 = conformer2->getAtomPos(11) - conformer2->getAtomPos(4);
    TEST_ASSERT(RDKit::feq(bond4_11Conf2.x, -0.773, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond4_11Conf2.y, 1.339, 1.e-3));
  }
  {
    auto noradrenalineMJCopy =
        std::unique_ptr<RWMol>(new RWMol(*noradrenalineMJ));
    const auto &conformer0 = noradrenalineMJCopy->getConformer(0);
    auto conformer1 = new Conformer(conformer0);
    noradrenalineMJCopy->addConformer(conformer1, true);
    auto scalingFactor =
        RDDepict::normalizeDepiction(*noradrenalineMJCopy, 1, -1);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  0) < 1.e-5);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  1) > 1.e-5);
    TEST_ASSERT(RDKit::feq(scalingFactor, 1.875, 1.e-3));
    auto conformer2 = new Conformer(*conformer1);
    noradrenalineMJCopy->addConformer(conformer2, true);
    auto bond10_11Conf0 = conformer0.getAtomPos(11) - conformer0.getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf0.x, 0.825, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf0.y, 0.0, 1.e-3));
    auto bond10_11Conf1 =
        conformer1->getAtomPos(11) - conformer1->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.x, 0.321, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.y, 1.513, 1.e-3));
    RDDepict::straightenDepiction(*noradrenalineMJCopy, 1);
    bond10_11Conf1 = conformer1->getAtomPos(11) - conformer1->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.x, 0.0, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.y, 1.547, 1.e-3));
    RDDepict::straightenDepiction(*noradrenalineMJCopy, 2, true);
    auto bond10_11Conf2 =
        conformer2->getAtomPos(11) - conformer2->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf2.x, bond10_11Conf1.x, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf2.y, bond10_11Conf1.y, 1.e-3));
  }
  {
    auto noradrenalineMJCopy =
        std::unique_ptr<RWMol>(new RWMol(*noradrenalineMJ));
    const auto &conformer0 = noradrenalineMJCopy->getConformer(0);
    auto conformer1 = new Conformer(conformer0);
    noradrenalineMJCopy->addConformer(conformer1, true);
    auto scalingFactor =
        RDDepict::normalizeDepiction(*noradrenalineMJCopy, 1, 0, 3.0);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  0) < 1.e-5);
    TEST_ASSERT(MolAlign::CalcRMS(*noradrenalineMJ, *noradrenalineMJCopy, 0,
                                  1) > 1.e-5);
    TEST_ASSERT(RDKit::feq(scalingFactor, 3.0, 1.e-3));
    auto conformer2 = new Conformer(*conformer1);
    noradrenalineMJCopy->addConformer(conformer2, true);
    auto conformer3 = new Conformer(*conformer1);
    noradrenalineMJCopy->addConformer(conformer3, true);
    auto bond10_11Conf0 = conformer0.getAtomPos(11) - conformer0.getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf0.x, 0.825, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf0.y, 0.0, 1.e-3));
    auto bond10_11Conf1 =
        conformer1->getAtomPos(11) - conformer1->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.x, 2.475, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.y, 0.0, 1.e-3));
    RDDepict::straightenDepiction(*noradrenalineMJCopy, 1);
    bond10_11Conf1 = conformer1->getAtomPos(11) - conformer1->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.x, 2.143, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf1.y, -1.237, 1.e-3));
    auto bond4_11Conf1 = conformer1->getAtomPos(11) - conformer1->getAtomPos(4);
    TEST_ASSERT(RDKit::feq(bond4_11Conf1.x, 0.0, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond4_11Conf1.y, 2.475, 1.e-3));
    RDDepict::straightenDepiction(*noradrenalineMJCopy, 2, true);
    auto bond10_11Conf2 =
        conformer2->getAtomPos(11) - conformer2->getAtomPos(10);
    auto bond10_11Conf3 =
        conformer3->getAtomPos(11) - conformer3->getAtomPos(10);
    TEST_ASSERT(RDKit::feq(bond10_11Conf2.x, bond10_11Conf3.x, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond10_11Conf2.y, bond10_11Conf3.y, 1.e-3));
    auto bond4_11Conf2 = conformer2->getAtomPos(11) - conformer2->getAtomPos(4);
    auto bond4_11Conf3 = conformer3->getAtomPos(11) - conformer3->getAtomPos(4);
    TEST_ASSERT(RDKit::feq(bond4_11Conf2.x, bond4_11Conf3.x, 1.e-3));
    TEST_ASSERT(RDKit::feq(bond4_11Conf2.y, bond4_11Conf3.y, 1.e-3));
  }
  {
    auto zeroCoordCTab = R"RES(
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
)RES";
    std::unique_ptr<RWMol> zeroCoordBenzene(MolBlockToMol(zeroCoordCTab));
    auto res = RDDepict::normalizeDepiction(*zeroCoordBenzene);
    TEST_ASSERT(res < 0.);
    TEST_ASSERT(MolToMolBlock(*zeroCoordBenzene) == zeroCoordCTab);
  }
  {
    // cyclopentadiene which is already straight should not be biased
    // towards a 30-degree angle rotate since it has no bonds
    // whose angle with the X axis is multiple of 60 degrees
    auto cpSittingOnHorizontalBondCTab = R"RES(
  MJ201100                      

  5  5  0  0  0  0  0  0  0  0999 V2000
   -2.3660    0.3892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0334   -0.0957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7785   -0.8803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9535   -0.8803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6986   -0.0957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  5  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
M  END
)RES";
    std::unique_ptr<RWMol> cpSittingOnHorizontalBond(
        MolBlockToMol(cpSittingOnHorizontalBondCTab));
    std::unique_ptr<RWMol> cpSittingOnHorizontalBondCopy(
        new RWMol(*cpSittingOnHorizontalBond));
    RDDepict::straightenDepiction(*cpSittingOnHorizontalBond);
    TEST_ASSERT(MolAlign::CalcRMS(*cpSittingOnHorizontalBond,
                                  *cpSittingOnHorizontalBondCopy) < 1.e-3);
    RDGeom::Transform3D trans;
    // rotate by 90 degrees
    trans.SetRotation(0.5 * M_PI, RDGeom::Z_Axis);
    MolTransforms::transformConformer(cpSittingOnHorizontalBond->getConformer(),
                                      trans);
    cpSittingOnHorizontalBondCopy.reset(new RWMol(*cpSittingOnHorizontalBond));
    RDDepict::straightenDepiction(*cpSittingOnHorizontalBond);
    TEST_ASSERT(MolAlign::CalcRMS(*cpSittingOnHorizontalBond,
                                  *cpSittingOnHorizontalBondCopy) < 1.e-3);
  }
}

void testValidRingSystemTemplates() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Test that ring system templates are valid "
      << std::endl;
  constexpr double RDKIT_BOND_LEN = 1.5;
  for (auto &smiles : TEMPLATE_SMILES) {
    std::unique_ptr<ROMol> mol{SmilesToMol(smiles)};
    RDDepict::CoordinateTemplates::assertValidTemplate(*mol, smiles);

    // also check whether the bonds in the template are the correct length
    double avg_length = 0.0;
    const Conformer &conf = mol->getConformer();
    for (auto &bond : mol->bonds()) {
      auto bond_length = (conf.getAtomPos(bond->getBeginAtomIdx()) -
                          conf.getAtomPos(bond->getEndAtomIdx()))
                             .length();
      avg_length += bond_length;
    }
    avg_length /= mol->getNumBonds();

    // this is a loose tolerance, since some complicated ring systems may have
    // odd bond lengths
    bool valid_length = RDKit::feq(avg_length, RDKIT_BOND_LEN, 0.1);
    if (!valid_length) {
      BOOST_LOG(rdWarningLog)
          << "Template has invalid average bond "
          << "length of " << avg_length << ": " << smiles << std::endl;
    }
    TEST_ASSERT(valid_length);
  }
}

int main() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = false;
#endif

  RDLog::InitLogs();
#if 1
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   test1 \n";
  test1();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   testCollisions \n";
  testCollisions();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   testAddHs \n";
  testAddHs();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   test198 \n";
  testIssue198();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   test2 \n";
  test2();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   test3 \n";
  test3();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   test4 \n";
  test4();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   tempTest \n";
  tempTest();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Issue248 \n";
  testIssue248();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Queries \n";
  testQueries();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test crashes associated with RemoveHs \n";
  testRemoveHsCrash();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 2091304\n";
  testIssue2091304();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 2821647\n";
  testIssue2821647();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 2948402\n";
  testIssue2948402();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 2995724\n";
  testIssue2995724();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Testing a change of the depictor bond length\n";
  testBondLengthChange();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 3122141\n";
  testIssue3122141();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 3135833\n";
  testIssue3135833();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 3487469\n";
  testIssue3487469();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test GitHub Issue 8\n";
  testGitHubIssue8();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test GitHub Issue 78\n";
  testGitHubIssue78();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test GitHub Issue 910\n";
  testGitHubIssue910();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test GitHub Issue 1073\n";
  testGitHubIssue1073();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   testConstrainedCoords\n";
  testConstrainedCoords();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "   Test GitHub Issue 1112: Bad coordinate generation for H2\n";
  testGitHubIssue1112();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";

  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test Issue 2303566\n";
  testIssue2303566();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  BOOST_LOG(rdInfoLog) << "   Test GitHub Issue 1286: "
                          "GenerateDepictionMatching2DStructure isn't matching "
                          "2D structure\n";
  testGitHubIssue1286();
  BOOST_LOG(rdInfoLog)
      << "***********************************************************\n";
  testGithub1691();
#endif
  testGithub2027();
  testGenerate2DDepictionRefPatternMatchVect();
  testGenerate2DDepictionAllowRGroupsOrig();
  testGenerate2DDepictionAllowRGroups();
  testNormalizeStraighten();
  testValidRingSystemTemplates();

  return (0);
}
