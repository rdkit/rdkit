//
//  Copyright (C) 2004-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <vector>
#include <iostream>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <iostream>
#include "BoundsMatrixBuilder.h"
#include "Embedder.h"
#include <cstdlib>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ROMol.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <RDGeneral/FileParseException.h>
#include <ForceField/ForceField.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include <cmath>
#include <RDGeneral/Exceptions.h>

#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

using namespace RDKit;

void test1() {
  boost::logging::disable_logs("rdApp.warning");

  std::string smiString =
      "CC1=C(C(C)=CC=C2)C2=CC=C1 c1ccccc1C C/C=C/CC \
                           C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4 C1CCCCS1(=O)(=O) c1ccccc1 \
                           C1CCCC1 C1CCCCC1 \
                           C1CC1(C)C C12(C)CC1CC2";
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/initCoords.sdf";
  SDMolSupplier sdsup(fname);
  // SDWriter writer("foo.sdf");

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString, spaceSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string smi = *token;
    RWMol *m = SmilesToMol(smi, 0, 1);
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, true, false, 2, true, 1,
                                          nullptr, 1e-2);
    CHECK_INVARIANT(cid >= 0, "");
    ROMol *m2 = sdsup.next();
    // BOOST_LOG(rdInfoLog) << ">>> " << smi << std::endl;
    // writer.write(*m);
    // writer.flush();

    // ROMol *m2 = NULL;
    if (m2) {
      unsigned int nat = m->getNumAtoms();

      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
#if 0
      BOOST_LOG(rdInfoLog) << "-----------------------" << std::endl;
      BOOST_LOG(rdInfoLog) << MolToMolBlock(*m2) << std::endl;
      BOOST_LOG(rdInfoLog) << "---" << std::endl;
      BOOST_LOG(rdInfoLog) << MolToMolBlock(*m) << std::endl;
      BOOST_LOG(rdInfoLog) << "-----------------------" << std::endl;
#endif
      for (unsigned int i = 0; i < nat; i++) {
        RDGeom::Point3D pt1i = conf1.getAtomPos(i);
        RDGeom::Point3D pt2i = conf2.getAtomPos(i);
        for (unsigned int j = i + 1; j < nat; j++) {
          RDGeom::Point3D pt1j = conf1.getAtomPos(j);
          RDGeom::Point3D pt2j = conf2.getAtomPos(j);
          double d1 = (pt1j - pt1i).length();
          double d2 = (pt2j - pt2i).length();
          if (m->getBondBetweenAtoms(i, j)) {
            // BOOST_LOG(rdInfoLog) << ">1> " <<i<<","<<j<<":"<< d1 << " " << d2
            // << std::endl;
            TEST_ASSERT(fabs(d1 - d2) / d1 < 0.06);
          } else {
            // BOOST_LOG(rdInfoLog) << ">2> " <<i<<","<<j<<":"<< d1 << " " << d2
            // << " "<<fabs(d1-d2)/d1<<std::endl;
            TEST_ASSERT(fabs(d1 - d2) / d1 < 0.12);
          }
        }
      }
    }
    delete m;
    delete m2;
  }
  boost::logging::enable_logs("rdApp.warning");
}

void computeDistMat(const RDGeom::PointPtrVect &origCoords,
                    RDNumeric::DoubleSymmMatrix &distMat) {
  unsigned int N = origCoords.size();
  CHECK_INVARIANT(N == distMat.numRows(), "");
  unsigned int i, j;
  RDGeom::Point3D pti, ptj;
  double d;
  for (i = 1; i < N; i++) {
    pti = *(RDGeom::Point3D *)origCoords[i];
    for (j = 0; j < i; j++) {
      ptj = *(RDGeom::Point3D *)origCoords[j];
      ptj -= pti;
      d = ptj.length();
      distMat.setVal(i, j, d);
    }
  }
}

void computeMolDmat(ROMol &mol, RDNumeric::DoubleSymmMatrix &distMat) {
  RDGeom::PointPtrVect origCoords;
  unsigned int i, nat = mol.getNumAtoms();
  Conformer &conf = mol.getConformer(0);
  for (i = 0; i < nat; i++) {
    origCoords.push_back(&conf.getAtomPos(i));
  }
  computeDistMat(origCoords, distMat);
}

void test2() {
  boost::logging::disable_logs("rdApp.warning");
  // check for in ring distances, and distances containing two bonds in a ring
  std::string smi = "Cc1c(C=CC(C)N2)c2[nH]n1";
  ROMol *mol = SmilesToMol(smi, 0, 1);
  DistGeom::BoundsMatPtr bm;
  RDNumeric::DoubleSymmMatrix *dmat;
  int cid;
  unsigned int nat = mol->getNumAtoms();

#if 1
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);

  std::cerr << "go" << std::endl;
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);

  TEST_ASSERT((bm->getUpperBound(0, 9) - bm->getLowerBound(0, 9)) < 0.13);
  TEST_ASSERT((bm->getUpperBound(0, 9) - dmat->getVal(0, 9) > -0.1) &&
              (bm->getLowerBound(0, 9) - dmat->getVal(0, 9) < 0.10));

  TEST_ASSERT((bm->getUpperBound(10, 7) - bm->getLowerBound(10, 7)) < 0.13);
  TEST_ASSERT((bm->getUpperBound(10, 7) - dmat->getVal(10, 7) > -0.1) &&
              (bm->getLowerBound(10, 7) - dmat->getVal(10, 7) < 0.10));

  TEST_ASSERT((bm->getUpperBound(2, 5) - bm->getLowerBound(2, 5)) < 0.20);
  TEST_ASSERT((bm->getUpperBound(2, 5) - dmat->getVal(2, 5) > -0.1) &&
              (bm->getLowerBound(2, 5) - dmat->getVal(2, 5) < 0.10));

  TEST_ASSERT((bm->getUpperBound(8, 4) - bm->getLowerBound(8, 4)) > 1.);
  TEST_ASSERT((bm->getUpperBound(8, 4) - bm->getLowerBound(8, 4)) < 1.2);
  TEST_ASSERT((bm->getUpperBound(8, 4) - dmat->getVal(8, 4) > -0.1) &&
              (bm->getLowerBound(8, 4) - dmat->getVal(8, 4) < 0.10));

  TEST_ASSERT((bm->getUpperBound(8, 6) - bm->getLowerBound(8, 6)) > 1.0);
  TEST_ASSERT((bm->getUpperBound(8, 6) - bm->getLowerBound(8, 6)) < 1.2);
  TEST_ASSERT((bm->getUpperBound(8, 6) - dmat->getVal(8, 6) > -0.1) &&
              (bm->getLowerBound(8, 6) - dmat->getVal(8, 6) < 0.10));

  delete mol;
  delete dmat;

  // chain of singles
  smi = "CCCC";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) > 1.0);
  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < 1.3);
  TEST_ASSERT((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1) &&
              (bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));

  delete mol;
  delete dmat;

  smi = "C=C=C=C";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);

  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
  TEST_ASSERT(bm->getUpperBound(0, 3) > dmat->getVal(0, 3));
  // this is kinda goofy but this linear molecule doesn't satisfy the bounds
  // completely
  TEST_ASSERT(fabs(bm->getLowerBound(0, 3) - dmat->getVal(0, 3)) < 0.2);

  delete mol;
  delete dmat;
#endif
  std::cerr << "-------------------------------------\n\n";
  smi = "C/C=C/C";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  // std::cerr << "\n-----\n" << MolToMolBlock(mol,false,cid) << std::endl;;
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
  // std::cerr << bm->getUpperBound(0,3) << " " << dmat->getVal(0,3) << " " <<
  // bm->getLowerBound(0,3) << std::endl;
  TEST_ASSERT((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1) &&
              (bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.1));

  delete mol;
  delete dmat;

  smi = "C/C=C\\C";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);

  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
  TEST_ASSERT((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1) &&
              (bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));

  delete mol;
  delete dmat;

  smi = "CC=CC";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < 1.13);
  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) > 1.);
  TEST_ASSERT((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1) &&
              (bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));

  delete mol;
  delete dmat;

  smi = "O=S-S=O";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid > -1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
  TEST_ASSERT((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1) &&
              (bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));

  delete mol;
  delete dmat;

#if 0
  // this next block of tests all handle the special case that led to Issue284
  smi = "COC=O";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  //TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  double x,y,z;
  x = bm->getUpperBound(0,3);
  y = bm->getLowerBound(0,3);
  z = dmat->getVal(0,3);
  //TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
  //             && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10 ));

  delete mol;
  delete dmat;

  smi = "C[NH]C=O";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  //TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
  //             && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10));

  delete mol;
  delete dmat;
#endif
  boost::logging::enable_logs("rdApp.warning");
}

void test3() {
  boost::logging::disable_logs("rdApp.warning");

  // check embedding based based on distances calculated from previous created
  // (good) coordinates
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/combi_coords.sdf";
  SDMolSupplier sdsup(fname);

  unsigned int i, j, nat;
  bool gotCoords;
  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    RDGeom::PointPtrVect origCoords, newCoords;
    nat = mol->getNumAtoms();
    Conformer &conf = mol->getConformer(0);
    for (i = 0; i < nat; i++) {
      origCoords.push_back(&conf.getAtomPos(i));
    }
    RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
    computeDistMat(origCoords, distMat);

    gotCoords = DistGeom::computeInitialCoords(distMat, origCoords);

    CHECK_INVARIANT(gotCoords, "");
    RDNumeric::DoubleSymmMatrix distMatNew(nat, 0.0);
    computeDistMat(origCoords, distMatNew);

    for (i = 1; i < nat; i++) {
      for (j = 0; j < i; j++) {
        CHECK_INVARIANT(
            RDKit::feq(distMat.getVal(i, j), distMatNew.getVal(i, j), 0.01),
            "");
      }
    }
    delete mol;
  }
  boost::logging::enable_logs("rdApp.warning");
}

void test4() {
  boost::logging::disable_logs("rdApp.warning");

  std::string smi =
      "c1cc(C(F)(F)F)ccc1/C=N/NC(=O)c(n2)c[n]3cc(C(F)(F)F)cc(c23)Cl";
  ROMol *m = SmilesToMol(smi, 0, 1);
  DGeomHelpers::EmbedMolecule(*m, 10, 1);  // etCoords(*m, iter);
  std::string fname = "test.mol";
  MolToMolFile(*m, fname);
  delete m;
  boost::logging::enable_logs("rdApp.warning");
}

void test5() {
  // some real CDK2 molecules lets see how many fail
  std::string rdbase = getenv("RDBASE");
  std::string smifile =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);

  int cid;
  while (1) {
    try {
      std::unique_ptr<RWMol> mol{static_cast<RWMol *>(smiSup.next())};
      MolOps::addHs(*mol);
      cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);  // getCoords(*mol, iter);
      TEST_ASSERT(cid > -1);
    } catch (FileParseException &) {
      break;
    }
  }
}

void test6() {
  std::string smi;
  ROMol *m;
  DistGeom::BoundsMatPtr bm;

  m = SmilesToMol("CC");
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  TEST_ASSERT(bm);
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  TEST_ASSERT(feq(bm->getLowerBound(0, 1), 0.0));
  TEST_ASSERT(feq(bm->getLowerBound(1, 0), 0.0));
  TEST_ASSERT(feq(bm->getUpperBound(0, 1), 1000.0));
  TEST_ASSERT(feq(bm->getUpperBound(1, 0), 1000.0));

  DGeomHelpers::setTopolBounds(*m, bm);
  TEST_ASSERT(bm->getLowerBound(0, 1) > 0.0);
  TEST_ASSERT(bm->getUpperBound(0, 1) < 1000.0);
  TEST_ASSERT(bm->getLowerBound(0, 1) < bm->getUpperBound(0, 1));

  delete m;
}

void testIssue215() {
  std::string smi;
  ROMol *m;
  DistGeom::BoundsMatPtr bm;
  bool ok;

  m = SmilesToMol("C=C1C2CC1C2");
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  TEST_ASSERT(bm);
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*m, bm);

  // this was the specific problem:
  TEST_ASSERT(bm->getUpperBound(0, 4) < 100.0);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  delete m;
}

void _computeStats(const std::vector<double> &vals, double &mean,
                   double &stdDev) {
  unsigned int N = vals.size();
  mean = 0.0;
  double sqSum = 0.0;
  stdDev = 0.0;
  std::vector<double>::const_iterator vci;
  for (vci = vals.begin(); vci != vals.end(); vci++) {
    mean += (*vci);
    sqSum += (*vci) * (*vci);
  }
  mean /= N;
  sqSum /= N;
  stdDev = sqrt(sqSum - (mean * mean));
}

void testTemp() {
  // ROMol *m =
  // SmilesToMol("CN(C)c(cc1)ccc1\\C=C(/C#N)c(n2)c(C#N)c([n]3cccc3)[n]2c(c4)cccc4");
  // ROMol *m =
  // SmilesToMol("C12C(=O)N(Cc3ccccc3)C(=O)C1C(C(=O)OC)(NC2c4ccc(Cl)c(c4)[N+]([O-])=O)Cc5c6ccccc6[nH]c5");
  // c1(n2cccc2)n(c3ccccc3)ncc1
  // ROMol *m = SmilesToMol("N#CCc1c(C#N)cn(C)n1");
  // ROMol *m = SmilesToMol("Cc1c(C#N)cn(C)n1");
  std::string rdbase = getenv("RDBASE");
  std::string smifile =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);
  ROMol *m;
  unsigned int cnt = 0;
  while (!smiSup.atEnd()) {
    m = smiSup.next();
    unsigned int i;
    std::vector<double> energies;
    double energy;
    for (i = 0; i < 20; i++) {
      int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, false);
      if (cid >= 0) {
        ForceFields::ForceField *ff = UFF::constructForceField(*m, 10, cid);
        ff->initialize();
        energy = ff->calcEnergy();
        energies.push_back(energy);
        delete ff;
      }
    }
    double mean, stdDev;
    _computeStats(energies, mean, stdDev);
    std::string mname;
    cnt++;
    m->getProp(common_properties::_Name, mname);
    BOOST_LOG(rdDebugLog) << cnt << "," << mname << "," << mean << "," << stdDev
                          << "\n";
    delete m;
  }
}

void test15Dists() {
  ROMol *m = SmilesToMol("c1ccccc1C");
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DGeomHelpers::initBoundsMat(mat);
  DistGeom::BoundsMatPtr mmat(mat);
  DGeomHelpers::setTopolBounds(*m, mmat);
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(2, 6), 4.32, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(2, 6), 4.16, 0.01), "");
  delete m;

  m = SmilesToMol("CC1=C(C(C)=CC=C2)C2=CC=C1");
  nat = m->getNumAtoms();
  mmat.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(mmat);
  DGeomHelpers::setTopolBounds(*m, mmat);

  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0, 4), 2.31, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0, 4), 2.47, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(4, 11), 4.11, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(4, 11), 4.27, 0.01), "");

  delete m;

  m = SmilesToMol("C/C=C/C=C/C", 0, 1);
  nat = m->getNumAtoms();

  mmat.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(mmat);
  DGeomHelpers::setTopolBounds(*m, mmat);

  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0, 4), 4.1874), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0, 4), 4.924), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(1, 5), 4.1874), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(1, 5), 4.924), "");

  delete m;
  m = SmilesToMol("NCc(c1)cccc1");
  delete m;
}

void testMultipleConfs() {
  boost::logging::disable_logs("rdApp.warning");

  std::string smi = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC";
  ROMol *m = SmilesToMol(smi, 0, 1);
  INT_VECT cids =
      DGeomHelpers::EmbedMultipleConfs(*m, 10, 30, 100, true, false, -1);
  INT_VECT_CI ci;
  // SDWriter writer("junk.sdf");
  double energy;

  for (ci = cids.begin(); ci != cids.end(); ci++) {
    // writer.write(*m, *ci);
    ForceFields::ForceField *ff = UFF::constructForceField(*m, 10, *ci);
    ff->initialize();
    energy = ff->calcEnergy();
    // BOOST_LOG(rdInfoLog) << energy << std::endl;
    TEST_ASSERT(energy > 100.0);
    TEST_ASSERT(energy < 300.0);
    delete ff;
  }
  delete m;
  boost::logging::enable_logs("rdApp.warning");
}

void testMultipleConfsExpTors() {
  boost::logging::disable_logs("rdApp.warning");
  std::string smi = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC";
  ROMol *m = SmilesToMol(smi, 0, 1);
  INT_VECT cids = DGeomHelpers::EmbedMultipleConfs(
      *m, 10, 30, 100, true, false, -1, true, 1, -1.0, nullptr, 1e-3, false,
      true, false, false, false, 5.0, false, 1, false, false);

  INT_VECT_CI ci;
  // SDWriter writer("junk.sdf");
  double energy;

  for (ci = cids.begin(); ci != cids.end(); ci++) {
    // writer.write(*m, *ci);
    ForceFields::ForceField *ff = UFF::constructForceField(*m, 10, *ci);
    ff->initialize();
    energy = ff->calcEnergy();
    // BOOST_LOG(rdInfoLog) << energy << std::endl;
    TEST_ASSERT(energy > 50.0);
    TEST_ASSERT(energy < 300.0);
    delete ff;
  }
  delete m;
  boost::logging::enable_logs("rdApp.warning");
}

void testOrdering() {
  std::string smi = "CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DGeomHelpers::initBoundsMat(mat);
  DistGeom::BoundsMatPtr mmat(mat);
  DGeomHelpers::setTopolBounds(*m, mmat);
  delete m;

  std::string smi2 = "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2";
  ROMol *m2 = SmilesToMol(smi2, 0, 1);
  auto *mat2 = new DistGeom::BoundsMatrix(nat);
  DGeomHelpers::initBoundsMat(mat2);
  DistGeom::BoundsMatPtr mmat2(mat2);
  DGeomHelpers::setTopolBounds(*m2, mmat2);
  delete m2;
}

#if 1
void testIssue227() {
  std::string smi =
      "CCOP1(OCC)=CC(c2ccccc2)=C(c2ccc([N+]([O-])=O)cc2)N=C1c1ccccc1";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  bool ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;

  smi = "OC(=O)c1cc2cc(c1)-c1c(O)c(ccc1)-c1cc(C(O)=O)cc(c1)OCCOCCO2";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  auto *nmat = new DistGeom::BoundsMatrix(nat);
  bm.reset(nmat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;
}
#endif

void testIssue236() {
  std::string smi =
      "Cn1c2n(-c3ccccc3)c(=O)c3c(nc4ccc([N+]([O-])=O)cc4c3)c2c(=O)n(C)c1=O";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  bool ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;

  smi = "Cc1cccc2c1c(C3=CCC3)c(C)cc2";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  auto *nmat = new DistGeom::BoundsMatrix(nat);
  bm.reset(nmat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;
}

void testIssue244() {
  // the bug here was a crash during the setting of the topological
  // bounds, so just completing the calls to setTopolBounds() indicates
  // success
  std::string smi = "NC1C(=O)N2C1SCC(Cl)=C2C(O)=O";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  delete m;

  smi = "NC1C(=O)N2C1SCC(Cl)=C2C(O)=O";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  delete m;

  smi = "CC1(C)SC2N(C1C(O)=O)C(=O)C2N";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  delete m;
}

void testIssue251() {
  std::string smi = "COC=O";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  TEST_ASSERT(RDKit::feq(bm->getLowerBound(0, 3), 2.67, 0.01));
  TEST_ASSERT(RDKit::feq(bm->getUpperBound(0, 3), 2.79, 0.01));
  delete m;
}

void testIssue276() {
  bool ok;
  std::string smi = "CP1(C)=CC=CN=C1C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
#if 0
  for(unsigned int i=0;i<nat;i++){
    for(unsigned int j=0;j<nat;j++){
      if(i<j) std::cout << bm->getUpperBound(i,j) << " ";
      else if(i>j) std::cout << bm->getLowerBound(i,j) << " ";
      else std::cout << "0.00000" << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  delete m;
}

void testIssue284() {
  bool ok;
  std::string smi = "CNC(=O)C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  // amide bonds are cis-oid:
  TEST_ASSERT(bm->getLowerBound(0, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(0, 3) < 3.0);

  delete m;

  smi = "CN(C)C(=O)C";
  m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  nat = m->getNumAtoms();

  auto *mat2 = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm2(mat2);
  DGeomHelpers::initBoundsMat(bm2);
  DGeomHelpers::setTopolBounds(*m, bm2);

  ok = DistGeom::triangleSmoothBounds(bm2);
  TEST_ASSERT(ok);

  // we've got no information to tell us cis-oid vs trans-oid here, so
  // the windows are huge:
  TEST_ASSERT(bm2->getLowerBound(0, 4) < 3.0);
  TEST_ASSERT(bm2->getUpperBound(0, 4) > 3.5);
  TEST_ASSERT(bm2->getLowerBound(2, 4) < 3.0);
  TEST_ASSERT(bm2->getUpperBound(2, 4) > 3.5);
  TEST_ASSERT(bm->getLowerBound(0, 3) < bm2->getLowerBound(0, 4));
  TEST_ASSERT(bm->getUpperBound(0, 3) < bm2->getUpperBound(0, 4));
  TEST_ASSERT(bm->getLowerBound(0, 3) < bm2->getLowerBound(2, 4));
  TEST_ASSERT(bm->getUpperBound(0, 3) < bm2->getUpperBound(2, 4));

  delete m;
}

void testIssue285() {
  bool ok;
  std::string smi = "CNC(=O)C";
  RWMol *m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  MolOps::addHs(*m);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  unsigned int tgtNumber = 10;
  INT_VECT cids = DGeomHelpers::EmbedMultipleConfs(*m, tgtNumber);
  TEST_ASSERT(cids.size() == tgtNumber);

  std::vector<std::string> molBlocks;
  for (INT_VECT_CI cid = cids.begin(); cid != cids.end(); ++cid) {
    molBlocks.push_back(MolToMolBlock(*m, true, *cid));
  }
  for (std::vector<std::string>::const_iterator mbI = molBlocks.begin();
       mbI != molBlocks.end(); ++mbI) {
    for (auto mbJ = mbI + 1; mbJ != molBlocks.end(); ++mbJ) {
      TEST_ASSERT((*mbI) != (*mbJ));
    }
    // std::cerr << (*mbI) << "\n$$$$\n";
  }
  delete m;
}

void testIssue355() {
  bool ok;
  std::string smi = "CNC(=O)C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  unsigned int nat = m->getNumAtoms();
  auto *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  TEST_ASSERT(bm->getLowerBound(0, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(0, 3) < 3.0);

  TEST_ASSERT(bm->getUpperBound(0, 4) > 3.2);
  TEST_ASSERT(bm->getLowerBound(0, 4) > 3.2);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  delete m;

  smi = "CNC(=O)NC";
  m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  TEST_ASSERT(bm->getLowerBound(0, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(0, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(0, 4) > 3.2);
  TEST_ASSERT(bm->getLowerBound(0, 4) > 3.2);
  TEST_ASSERT(bm->getLowerBound(5, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(5, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(5, 1) > 3.2);
  TEST_ASSERT(bm->getLowerBound(5, 1) > 3.2);

  delete m;

  smi = "CNC(=O)Nc1ccccc1";
  m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  TEST_ASSERT(bm->getLowerBound(0, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(0, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(0, 4) > 3.2);
  TEST_ASSERT(bm->getLowerBound(0, 4) > 3.2);
  TEST_ASSERT(bm->getLowerBound(5, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(5, 3) < 3.0);
  TEST_ASSERT(bm->getUpperBound(5, 1) > 3.2);
  TEST_ASSERT(bm->getLowerBound(5, 1) > 3.2);

  delete m;
}

void testRandomCoords() {
  std::string smiString =
      "CC1=C(C(C)=CC=C2)C2=CC=C1 c1ccccc1C C/C=C/CC \
                           C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4 C1CCCCS1(=O)(=O) c1ccccc1 \
                           C1CCCC1 C1CCCCC1 \
                           C1CC1(C)C C12(C)CC1CC2";
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/initCoords.random.sdf";
  SDMolSupplier sdsup(fname, true, false);
  // SDWriter writer("foo.sdf");
  // SDWriter writer(fname);

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString, spaceSep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string smi = *token;
    // std::cerr << "SMI: " << smi << std::endl;
    ROMol *m = SmilesToMol(smi, 0, 1);
    auto *m2 = (RWMol *)MolOps::addHs(*m);
    delete m;
    m = m2;
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, true, true, 2, true, 1,
                                          nullptr, 1e-2);
    CHECK_INVARIANT(cid >= 0, "");
    // writer.write(*m);
    // writer.flush();
#if 1
    m2 = static_cast<RWMol *>(sdsup.next());
    // ROMol *m2 = NULL;
    if (m2) {
      TEST_ASSERT(m->getNumAtoms() == m2->getNumAtoms());
      unsigned int nat = m->getNumAtoms();

      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
#if 0
      BOOST_LOG(rdWarningLog) << "-----------------------" << std::endl;
      BOOST_LOG(rdWarningLog) << MolToMolBlock(*m2) << std::endl;
      BOOST_LOG(rdWarningLog) << "---" << std::endl;
      BOOST_LOG(rdWarningLog) << MolToMolBlock(*m) << std::endl;
      BOOST_LOG(rdWarningLog) << "-----------------------" << std::endl;
#endif
      for (unsigned int i = 0; i < nat; i++) {
        RDGeom::Point3D pt1i = conf1.getAtomPos(i);
        RDGeom::Point3D pt2i = conf2.getAtomPos(i);
        for (unsigned int j = i + 1; j < nat; j++) {
          RDGeom::Point3D pt1j = conf1.getAtomPos(j);
          RDGeom::Point3D pt2j = conf2.getAtomPos(j);
          double d1 = (pt1j - pt1i).length();
          double d2 = (pt2j - pt2i).length();
          if (m->getBondBetweenAtoms(i, j)) {
            TEST_ASSERT(fabs(d1 - d2) / d1 < 0.05);
          } else {
            TEST_ASSERT(fabs(d1 - d2) / d1 < 0.1);
          }
        }
      }
    }
    delete m2;
#endif
    delete m;
  }
}

void testIssue1989539() {
  {
    std::string smi = "c1ccccc1.Cl";
    ROMol *m = SmilesToMol(smi, 0, 1);
    auto *m2 = (RWMol *)MolOps::addHs(*m);
    delete m;
    m = m2;
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
  {
    std::string smi = "[Cl-].c1ccccc1C[NH3+]";
    RWMol *m = SmilesToMol(smi, 0, 1);
    MolOps::addHs(*m);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
}

void testConstrainedEmbedding() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf";
  SDMolSupplier sdsup(fname);

  ROMol *ref = sdsup.next();
  {
    auto *test = new RWMol(*ref);
    MolOps::addHs(*test);
    std::map<int, RDGeom::Point3D> coords;
    coords[0] = ref->getConformer().getAtomPos(0);
    coords[1] = ref->getConformer().getAtomPos(1);
    coords[2] = ref->getConformer().getAtomPos(2);
    coords[3] = ref->getConformer().getAtomPos(3);
    coords[4] = ref->getConformer().getAtomPos(4);

#if 1
    int cid = DGeomHelpers::EmbedMolecule(*test, 30, 22, true, false, 2., true,
                                          1, &coords);
    TEST_ASSERT(cid > -1);

    MatchVectType alignMap;
    alignMap.push_back(std::make_pair(0, 0));
    alignMap.push_back(std::make_pair(1, 1));
    alignMap.push_back(std::make_pair(2, 2));
    alignMap.push_back(std::make_pair(3, 3));
    alignMap.push_back(std::make_pair(4, 4));
    double ssd = MolAlign::alignMol(*test, *ref, -1, -1, &alignMap);
    BOOST_LOG(rdInfoLog) << "ssd: " << ssd << std::endl;
    TEST_ASSERT(ssd < 0.1);
#endif
    delete test;
  }

  {
    RWMol *test = static_cast<RWMol *>(sdsup.next());
    MolOps::addHs(*test);
    std::map<int, RDGeom::Point3D> coords;
    coords[4] = ref->getConformer().getAtomPos(0);
    coords[5] = ref->getConformer().getAtomPos(1);
    coords[6] = ref->getConformer().getAtomPos(2);
    coords[7] = ref->getConformer().getAtomPos(3);
    coords[8] = ref->getConformer().getAtomPos(4);
    int cid = DGeomHelpers::EmbedMolecule(*test, 30, 22, true, false, 2., true,
                                          1, &coords);
    TEST_ASSERT(cid > -1);

    MatchVectType alignMap;
    alignMap.push_back(std::make_pair(4, 0));
    alignMap.push_back(std::make_pair(5, 1));
    alignMap.push_back(std::make_pair(6, 2));
    alignMap.push_back(std::make_pair(7, 3));
    alignMap.push_back(std::make_pair(8, 4));
    double ssd = MolAlign::alignMol(*test, *ref, -1, -1, &alignMap);
    BOOST_LOG(rdInfoLog) << "ssd: " << ssd << std::endl;
    TEST_ASSERT(ssd < 0.1);
    delete test;
  }
  delete ref;
}

void testIssue2091864() {
  boost::logging::disable_logs("rdApp.warning");

  {
    std::string smi = "C1C2CC12";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    delete m;
  }
  {
    std::string smi = "C1CC2C3C1C23";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
  {
    std::string smi = "c1ccc2c(c1)C1C3C2C13";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
  boost::logging::enable_logs("rdApp.warning");
}

void testIssue2091974() {
  {
    std::string smi = "CCOC(OCC)(OCC)OCC";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 = MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    delete m2;
  }
  {
    std::string smi = "O=N(=O)OCC(CON(=O)=O)(CON(=O)=O)CON(=O)=O";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 = MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    delete m2;
  }
}

void testIssue2835784() {
  boost::logging::disable_logs("rdApp.warning");

#if 1
  {
    std::string smi = "C1C=C1";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
  {
    std::string smi = "C1C=C1";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 = MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m2, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m2;
  }
  {
    std::string smi = "C12=CCC1C2";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
#endif
  {
    std::string smi = "C12=CCC1C2";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 = MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m2, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m2;
  }
  boost::logging::enable_logs("rdApp.warning");
}

void testIssue3019283() {
  {
    std::string smi = "C1=C2C1C1CC21";
    RWMol *m = SmilesToMol(smi);
    MolOps::addHs(*m);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
}

void testIssue3238580() {
  {
    std::string smi = "C1CCC2=CC12";
    RWMol *m = SmilesToMol(smi);

    auto *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());
    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);

    delete m;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.1.mol";
    RWMol *m = MolFileToMol(molfile);

    auto *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());

    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);

    delete m;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.2.mol";
    RWMol *m = MolFileToMol(molfile);

    auto *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());
    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);

    delete m;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.3.mol";
    RWMol *m = MolFileToMol(molfile);

    auto *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());
    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);

    delete m;
  }
}

void testIssue3483968() {
  boost::logging::disable_logs("rdApp.warning");

  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3483968.mol";
    RWMol *m = MolFileToMol(molfile);
    TEST_ASSERT(m);

    int cid = DGeomHelpers::EmbedMolecule(*m, 0, -1, true, false, 2.0, true, 1,
                                          nullptr, 1e-3, true);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(
        *m, 10, 30, 1, true, false, 2.0, true, 1, -1.0, nullptr, 1e-3, true);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(), cids.end(), -1) == cids.end());
    delete m;
  }
  boost::logging::enable_logs("rdApp.warning");
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock(const std::vector<ROMol *> &mols,
              const std::vector<double> &energies, unsigned int count,
              unsigned int idx) {
  for (unsigned int j = 0; j < 100; j++) {
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      ROMol mol(*mols[i]);
      std::vector<int> cids =
          DGeomHelpers::EmbedMultipleConfs(mol, 10, 30, 0xFEED);
      TEST_ASSERT(cids.size() == 10);
      ForceFields::ForceField *field = nullptr;
      try {
        field = UFF::constructForceField(mol, 100.0, cids[0]);
      } catch (...) {
        field = nullptr;
      }
      TEST_ASSERT(field);
      field->initialize();
      double eng = field->calcEnergy();
      if (!feq(eng, energies[i])) {
        std::cerr << i << " iter " << j << " " << energies[i] << " != " << eng
                  << std::endl;
      }

      TEST_ASSERT(feq(eng, energies[i]));
      delete field;
    }
  }
};
}  // namespace

#include <thread>
#include <future>
void testMultiThread() {
  std::cerr << "building molecules" << std::endl;
  // std::string smi="C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4";
  std::string smi = "c1ccc2c(c1)C1C3C2C13";
  std::vector<ROMol *> mols;
  for (unsigned int i = 0; i < 100; ++i) {
    RWMol *m = SmilesToMol(smi);
    MolOps::addHs(*m);
    mols.push_back(m);
  }

  std::cerr << "generating reference data" << std::endl;
  std::vector<double> energies(mols.size(), 0.0);
  for (unsigned int i = 0; i < mols.size(); ++i) {
    ROMol mol(*mols[i]);
    std::vector<int> cids =
        DGeomHelpers::EmbedMultipleConfs(mol, 10, 30, 0xFEED);
    TEST_ASSERT(cids.size() == 10);
    ForceFields::ForceField *field = nullptr;
    try {
      field = UFF::constructForceField(mol, 100.0, cids[0]);
    } catch (...) {
      field = nullptr;
    }
    TEST_ASSERT(field);
    field->initialize();
    double eng = field->calcEnergy();
    TEST_ASSERT(eng != 0.0);
    energies[i] = eng;
    delete field;
  }

  std::cerr << "validating reference data" << std::endl;
  for (unsigned int i = 0; i < mols.size(); ++i) {
    ROMol mol(*mols[i]);
    std::vector<int> cids =
        DGeomHelpers::EmbedMultipleConfs(mol, 10, 30, 0xFEED);
    TEST_ASSERT(cids.size() == 10);
    ForceFields::ForceField *field = nullptr;
    try {
      field = UFF::constructForceField(mol, 100.0, cids[0]);
    } catch (...) {
      field = nullptr;
    }
    TEST_ASSERT(field);
    field->initialize();
    double eng = field->calcEnergy();
    TEST_ASSERT(feq(eng, energies[i]));
    delete field;
  }

  std::vector<std::future<void>> tg;
  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(
        std::async(std::launch::async, runblock, mols, energies, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }

  for (auto &mol : mols) {
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread() {}
#endif

void testGithub55() {
  boost::logging::disable_logs("rdApp.warning");

  {
    std::string smiles = "c1cnco1";
    RWMol *core = SmilesToMol(smiles);
    TEST_ASSERT(core);

    int cid = DGeomHelpers::EmbedMolecule(*core);
    TEST_ASSERT(cid >= 0);

    smiles = "o1cncc1C";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    std::map<int, RDGeom::Point3D> coords;
    coords[0] = core->getConformer().getAtomPos(4);
    coords[1] = core->getConformer().getAtomPos(3);
    coords[2] = core->getConformer().getAtomPos(2);
    coords[3] = core->getConformer().getAtomPos(1);
    coords[4] = core->getConformer().getAtomPos(0);
    cid = DGeomHelpers::EmbedMolecule(*mol, 50, 22, true, false, 2., true, 1,
                                      &coords);
    TEST_ASSERT(cid > -1);

    delete core;
    delete mol;
  }
  {
    std::string smiles = "c1cncs1";
    RWMol *core = SmilesToMol(smiles);
    TEST_ASSERT(core);

    int cid = DGeomHelpers::EmbedMolecule(*core);
    TEST_ASSERT(cid >= 0);

    smiles = "s1cncc1C";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    std::map<int, RDGeom::Point3D> coords;
    coords[0] = core->getConformer().getAtomPos(4);
    coords[1] = core->getConformer().getAtomPos(3);
    coords[2] = core->getConformer().getAtomPos(2);
    coords[3] = core->getConformer().getAtomPos(1);
    coords[4] = core->getConformer().getAtomPos(0);
    cid = DGeomHelpers::EmbedMolecule(*mol, 50, 22, true, false, 2., true, 1,
                                      &coords);
    TEST_ASSERT(cid > -1);

    delete core;
    delete mol;
  }
  boost::logging::enable_logs("rdApp.warning");
}

void testGithub256() {
  {
    auto *mol = new RWMol();
    TEST_ASSERT(mol);

    bool ok = false;
    try {
      DGeomHelpers::EmbedMolecule(*mol);
      ok = false;
    } catch (const ValueErrorException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
    delete mol;
  }
}

#ifdef RDK_TEST_MULTITHREADED
void testMultiThreadMultiConf() {
  boost::char_separator<char> sep("|");
  auto bldString = std::string(RDKit::rdkitBuild);
  tokenizer tokens(bldString, sep);
  std::vector<std::string> tokenVect(tokens.begin(), tokens.end());
  const double ENERGY_TOLERANCE = ((tokenVect[2] != "MINGW") ? 1.0e-6 : 1.0);
  const double MSD_TOLERANCE = ((tokenVect[2] != "MINGW") ? 1.0e-6 : 1.0e-5);
  std::string smi = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC";
  std::unique_ptr<RWMol> m{SmilesToMol(smi, 0, 1)};
  TEST_ASSERT(m);
  MolOps::addHs(*m);
  INT_VECT cids;
  ROMol m2(*m);
  DGeomHelpers::EmbedMultipleConfs(*m, cids, 200, 1, 30, 100, true, false, -1);
  DGeomHelpers::EmbedMultipleConfs(m2, cids, 200, 0, 30, 100, true, false, -1);
  INT_VECT_CI ci;

  for (ci = cids.begin(); ci != cids.end(); ci++) {
    ForceFields::ForceField *ff = UFF::constructForceField(*m, 100, *ci);
    ff->initialize();
    double e1 = ff->calcEnergy();
    const RDGeom::PointPtrVect &pVect = ff->positions();
    TEST_ASSERT(e1 > 100.0);
    TEST_ASSERT(e1 < 300.0);
    ForceFields::ForceField *ff2 = UFF::constructForceField(m2, 100, *ci);
    ff2->initialize();
    double e2 = ff2->calcEnergy();
    const RDGeom::PointPtrVect &p2Vect = ff2->positions();
    TEST_ASSERT(RDKit::feq(e1, e2, ENERGY_TOLERANCE));
    TEST_ASSERT(pVect.size() == p2Vect.size());
    double msd = 0.0;
    for (unsigned int i = 0; i < pVect.size(); ++i) {
      const auto *p = dynamic_cast<const RDGeom::Point3D *>(pVect[i]);
      const auto *p2 = dynamic_cast<const RDGeom::Point3D *>(p2Vect[i]);
      TEST_ASSERT(p && p2);
      msd += (*p - *p2).lengthSq();
    }
    msd /= static_cast<double>(pVect.size());
    TEST_ASSERT(msd < MSD_TOLERANCE);
    delete ff;
    delete ff2;
  }
}
#endif

void testGithub563() {
  {
    std::string smi = "[H][C@]1(C[NH3+])CC[C@]([H])(CNC)CC1";
    ROMol *m = SmilesToMol(smi);
    std::string csmi = MolToSmiles(*m, true);
    std::cerr << csmi << std::endl;
    for (unsigned int i = 1; i < 100; ++i) {
      RWMol m2 = *m;
      MolOps::addHs(m2);
      DGeomHelpers::EmbedMolecule(m2, 50, i);
      MolOps::assignChiralTypesFrom3D(m2);
      MolOps::removeHs(m2);
      std::string smi = MolToSmiles(m2, true);
      TEST_ASSERT(smi == csmi);
    }
    delete m;
  }
}

void testGithub568() {
  {
    // sample molecules (either from ChEMBL or derived from ChEMBL) that were
    // problematic
    std::string smis[] = {
        "C1CN2C[C@@H]1[C@@]1(CN=CO1)C2",
        "OC(=O)[C@@H]1[C@H]2CC[C@H](O2)[C@@H]1C(=O)O",
        "O=C(CCc1ccccc1)OC[C@H]2C[C@@H]3O[C@H]2[C@@H]4[C@H]3C(=O)OC4=O",
        "Cn1cc(C2=NC[C@]3(CN4CC[C@@H]3C4)O2)c5ccccc15",
        "Nc1ncnc2c1ncn2[C@@H]3C=C(CO)[C@@H](O)[C@H]3O",
        "CO[C@H]1[C@@H](N)[C@@H](O)[C@@H](O)[C@H]1O",
        "CO[C@@H]1[C@@H](N)[C@@H](O)[C@@H](O)[C@H]1O",
        "CCS[C@@H]1[C@@H](N)[C@@H](O)[C@@H](O)[C@H]1O",
        "CCCCC[C@H](O)\\C=C\\[C@@]1(F)[C@H](O)C[C@H](O)[C@@H]1C\\C=C/"
        "CCCC(=O)OC",
        "CCCCC[C@@H](O)\\C=C\\[C@]1(F)[C@@H](O)C[C@@H](O)[C@H]1C\\C=C/"
        "CCCC(=O)OC",
        "CNC(=O)Oc1cccc2[C@@H]3[C@@H](CCN3C)COc12",
        "CNC(=O)Oc1ccc2OC[C@@H]3CCN(C)[C@@H]3c2c1",
        "CN1CC[C@H]2COc3c(O)cccc3[C@@H]12",
        "CC(=O)C[C@@]12CCC[C@H]1[C@@H]3CCC4=CC(=O)CC[C@@H]4[C@H]3CC2",
        "OC(=O)[C@@H]1C[C@H]2C[C@@H](CCc3nn[nH]n3)CC[C@H]2CN1",
        "COc1ccc2[C@H](O)[C@@H](COc2c1)N3CCC(O)(CC3)c4ccc(F)cc4",
        "O[C@@H]1[C@@H](COc2cc(O)ccc12)N3CCC(Cc4ccccc4)CC3",
        "O[C@H]1[C@H](COc2cc(O)ccc12)N3CCC(O)(CC3)c4ccccc4",
        "CC(=O)c1ccc2C[C@@H]3[C@@H]4CCCC[C@]4(CCN3CC5CC5)c2c1",
        "CC(=O)O[C@@H]1CC[C@@H]2CN3CCc4cc5OCOc5cc4[C@@H]3C[C@@H]2C1",
        "COc1cc2CCN3C[C@H]4CC[C@@H](O)C[C@H]4C[C@H]3c2cc1OC",
        "O[C@H]1C[C@H]2C[C@@H]3N(CCc4cc5OCOc5cc34)C[C@H]2C[C@H]1C#N",
        "COC(=O)[C@@H]1C[C@@](C)(NC(=O)N1)C(=O)O",
        "CC(=O)OC1=C(C[C@H]2Cc3cc4cccc(O)c4c(O)c3C(=O)[C@H]2C1)Sc5ccccc5",
        "Nc1ncnc2c1ncn2[C@H]3C[C@H](O)[C@@H](CO)C3",
        "C[C@@H](N1CC[C@@]23CCCC[C@@H]2[C@@H]1Cc4ccc(OCc5cccc(F)c5)cc34)C(=O)N",
        "CN1C(=O)CC[C@@]2(C)C1=CCc3cc(Cl)ccc23",
        "Cc1nc(COc2ccc3OC[C@H](Cc4cccnc4)[C@H](O)c3c2)ccc1[N+](=O)[O-]",
        "EOS"};
    for (unsigned int idx = 0; smis[idx] != "EOS"; ++idx) {
      ROMol *m = SmilesToMol(smis[idx]);
      std::string csmi = MolToSmiles(*m, true);
      std::cerr << csmi << std::endl;
      for (unsigned int i = 1; i < 20;
           ++i) {  // increase the limit here to make this a real torture test
        RWMol m2 = ROMol(*m);
        MolOps::addHs(m2);
        int cid = DGeomHelpers::EmbedMolecule(m2, 50, i);
        TEST_ASSERT(cid >= 0);
        MolOps::assignChiralTypesFrom3D(m2);

        // m2.setProp("_Name",smis[idx]);
        // std::cerr<<MolToMolBlock(m2)<<std::endl;
        // TEST_ASSERT(0);
        MolOps::removeHs(m2);
        std::string smi = MolToSmiles(m2, true);
        if (smi != csmi) {
          std::cerr << "-------------" << std::endl;
          std::cerr << smis[idx] << " " << i << std::endl;
          std::cerr << smi << "\n" << csmi << std::endl;
          m2.setProp("_Name", smis[idx]);
          std::cerr << MolToMolBlock(m2) << std::endl;
        }
        TEST_ASSERT(smi == csmi);
      }
      delete m;
    }
  }
}

void testGithub696() {
  {
    std::string smi = "COc1ccc2CCCCCCCCCCc3ccc(OC)c(c3)-c1c2";
    ROMol *m = SmilesToMol(smi);
    // m->debugMol(std::cerr);
    DistGeom::BoundsMatPtr bm;

    bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm, false, false);

    TEST_ASSERT(bm->getUpperBound(2, 19) > bm->getLowerBound(2, 19));
    TEST_ASSERT(bm->getLowerBound(2, 19) > 2.0);
    TEST_ASSERT(bm->getUpperBound(2, 19) > 2.5);

    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);

    delete m;
  }
}

void testGithub697() {
  {  // a group of chembl molecules (and things derived from them), all of which
    // contain a c1cscn1 heterocycle
    std::string smis[] = {
        "C1SC2=NC1CCCCCC2",
        "C1CCCc2nc(CC1)cs2",
        "C1Cc2coc(n2)-c2coc(C1)n2",
        "C1Cc2coc(n2)-c2csc(C1)n2",
        "C1CCc2nc(cs2)-c2nc(C1)co2",
        "C1Cc2nc(co2)-c2nc(cs2)-c2nc1co2",
        "C1Cc2nc(co2)-c2nc(co2)-c2nc(cs2)-c2nc(co2)-c2nc1co2",
        "C1CNCc2coc(n2)-c2coc(n2)-c2csc(n2)-c2coc(n2)-c2coc(CNCCN1)n2",
        "C=C1NC(=O)C(C(C)C)NC(=O)C(C(C)CC)NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-"
        "c2csc(n2)-c2coc(n2)-c2coc1n2",
        "CCC(C)C1NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-c2csc(n2)-c2coc(n2)-c2coc("
        "n2)"
        "CNC(=O)C(C(C)C)NC1=O",
        "CCC(C)C1NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-c2csc(n2)-c2coc(n2)-c2coc("
        "n2)"
        "C(COC(C)=O)NC(=O)C(C(C)C)NC1=O",
        "CCC(C)C1NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-c2csc(n2)-c2coc(n2)-c2coc("
        "n2)"
        "C(COC(=O)COCCOCCOC)NC(=O)C(C(C)C)NC1=O",
        "C=C1NC(=O)C(C(C)C)NC(=O)C(C(C)CC)NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-"
        "c2csc(n2)-c2coc(n2)-c2nc1oc2C",
        "C=C1NC(=O)C(C(C)C)NC(=O)C(C(C)CC)NC(=O)c2nc(oc2-c2ccccc2)-c2nc(oc2C)-"
        "c2csc(n2)-c2coc(n2)-c2nc1oc2C",
        "CCC(C)C1NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-c2csc(n2)-c2coc(n2)-c2coc("
        "n2)"
        "C(COC(=O)COCCOC)NC(=O)C(C(C)C)NC1=O",
        "C=C1NC(=O)C(C(C)C)NC(=O)C(C(C)CC)NC(=O)c2nc(oc2-c2ccccc2)-c2nc(oc2C)-"
        "c2csc(n2)-c2coc(n2)-c2coc1n2",
        "CCC(C)C1NC(=O)c2nc(oc2-c2ccccc2)-c2coc(n2)-c2csc(n2)-c2coc(n2)-c2coc("
        "n2)"
        "C(COC(=O)C(N)CCCNC(=N)N)NC(=O)C(C(C)C)NC1=O",
        "EOS"};

    for (unsigned int idx = 0; smis[idx] != "EOS"; ++idx) {
      ROMol *m = SmilesToMol(smis[idx]);
      TEST_ASSERT(m);
      DistGeom::BoundsMatPtr bm;

      bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
      DGeomHelpers::initBoundsMat(bm);
      DGeomHelpers::setTopolBounds(*m, bm, false, false);
      bool ok = DistGeom::triangleSmoothBounds(bm);
      if (!ok) {
        m->debugMol(std::cerr);
        std::cerr << " FAILED: " << smis[idx] << std::endl;
      }
      TEST_ASSERT(ok);
      delete m;
    }
  }
}

namespace {
void compareConfs(const ROMol *m, const ROMol *expected, int molConfId = -1,
                  int expectedConfId = -1) {
  PRECONDITION(m, "bad pointer");
  PRECONDITION(expected, "bad pointer");
  TEST_ASSERT(m->getNumAtoms() == expected->getNumAtoms());
  const Conformer &conf1 = m->getConformer(molConfId);
  const Conformer &conf2 = expected->getConformer(expectedConfId);
  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    TEST_ASSERT(m->getAtomWithIdx(i)->getAtomicNum() ==
                expected->getAtomWithIdx(i)->getAtomicNum());

    RDGeom::Point3D pt1i = conf1.getAtomPos(i);
    RDGeom::Point3D pt2i = conf2.getAtomPos(i);
    TEST_ASSERT((pt1i - pt2i).length() < 0.05);
  }
}
}  // namespace

void testGithub971() {
  {
    // sample molecule found by Sereina
    std::string smi = "C/C(=C\\c1ccccc1)CN1C2CC[NH2+]CC1CC2";

    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    MolOps::addHs(*m);
    int cid = DGeomHelpers::EmbedMolecule(*m, 0, 0xf00d);
    TEST_ASSERT(cid >= 0);
    MolOps::removeHs(*m);
    std::string expectedMb = R"CTAB(
     RDKit          3D

 19 21  0  0  0  0  0  0  0  0999 V2000
    1.1886   -1.4168    0.8579 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0673   -0.0768    0.1995 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2169    0.4750   -0.1935 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5072    0.0051    0.3290 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3384    0.9895    0.8425 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6714    0.7376    1.1358 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1353   -0.5269    0.8881 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3159   -1.5094    0.3457 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1414   -1.0880   -0.2477 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0961    0.0762   -0.7307 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3281    0.3264   -0.0331 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1244    1.2097   -0.8827 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2111    1.8916   -0.0980 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9969    1.6208    1.3672 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3806    0.2980    1.7743 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4530   -0.6563    0.7070 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1661   -0.8675   -0.0283 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5534   -1.0960   -1.4823 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6831    0.3159   -1.9837 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  2 10  1  0
 10 11  1  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 16 17  1  0
 17 18  1  0
 18 19  1  0
  9  4  1  0
 17 11  1  0
 19 12  1  0
M  CHG  1  15   1
M  END)CTAB";
    RWMol *expected = MolBlockToMol(expectedMb);
    unsigned int nat = expected->getNumAtoms();
    TEST_ASSERT(nat == m->getNumAtoms());

    compareConfs(m, expected, 0, 0);
    delete m;
    delete expected;
  }
}

void testEmbedParameters() {
  std::string rdbase = getenv("RDBASE");
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    params.useExpTorsionAnglePrefs = true;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    params.useExpTorsionAnglePrefs = true;
    params.useBasicKnowledge = true;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("n1cccc(C)c1ON");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    params.useExpTorsionAnglePrefs = true;
    params.useBasicKnowledge = true;
    params.ETversion = 2;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    params.useBasicKnowledge = true;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  //------------
  // using the pre-defined parameter sets
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETDG);
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("OCCC");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params(DGeomHelpers::KDG);
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  // small ring torsions improvement test
  {
    std::string fname = rdbase +
                        "/Code/GraphMol/DistGeomHelpers/test_data/"
                        "simple_torsion.smallring.etkdgv3.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("C1CCCCC1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params(DGeomHelpers::srETKDGv3);
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  // macrocycles torsions backward compatibility test
  {
    std::string fname = rdbase +
                        "/Code/GraphMol/DistGeomHelpers/test_data/"
                        "simple_torsion.macrocycle.etkdg.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("O=C1NCCCCCCCCC1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
  // macrocycles torsions improvement test
  {
    std::string fname = rdbase +
                        "/Code/GraphMol/DistGeomHelpers/test_data/"
                        "simple_torsion.macrocycle.etkdgv3.mol";
    RWMol *ref = MolFileToMol(fname, true, false);
    TEST_ASSERT(ref);
    RWMol *mol = SmilesToMol("C1NCCCCCCCCC1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDGv3);
    params.randomSeed = 42;
    TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
    compareConfs(ref, mol);

    delete ref;
    delete mol;
  }
}

void testGithub1227() {
  {
    RWMol *m = SmilesToMol("CC");
    TEST_ASSERT(m);
    MolOps::addHs(*m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    INT_VECT cids;
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    params.randomSeed = 0xf00d;

    cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, params);
    TEST_ASSERT(cids.size() == 10);

    params.pruneRmsThresh = 0.5;
    cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, params);
    TEST_ASSERT(cids.size() == 1);

    params.onlyHeavyAtomsForRMS = false;  // the old default behavior
    params.useSymmetryForPruning = false;
    cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, params);
    TEST_ASSERT(cids.size() == 6);

    delete m;
  }
}

void testGithub1240() {
  {
    RWMol *mol = SmilesToMol("C1CCCCCCC1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    params.maxIterations = 1;
    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    delete mol;
  }
  {
    RWMol *mol = SmilesToMol("C1C3CC2CC(CC1C2)C3");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    params.maxIterations = 1;
    boost::logging::disable_logs("rdApp.warning");
    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    boost::logging::enable_logs("rdApp.warning");
    TEST_ASSERT(cid >= 0);
    delete mol;
  }
  {
    RWMol *mol = SmilesToMol("c1ccccccccc1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 20);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    params.maxIterations = 1;  // we should get this in one iteration

    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    delete mol;
  }
  {
    RWMol *mol = SmilesToMol("c1ccccccc1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 16);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    params.maxIterations = 1;  // we should get this in one iteration

    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    delete mol;
  }
  {
    RWMol *mol = SmilesToMol("c1ccccccccccc1");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 24);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    params.maxIterations = 1;  // we should get this in one iteration

    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    delete mol;
  }
  {
    // CHEMBL307150
    RWMol *mol =
        SmilesToMol("Cc1cc2ccn(C)c2c3c4C(=O)NC(=O)c4c5c6ccccc6[nH]c5c13");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    params = DGeomHelpers::ETKDG;
    params.randomSeed = 0xf00d;
    cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);

    // std::cerr << MolToMolBlock(*mol) << std::endl;
    delete mol;
  }

  {
    // CHEMBL43398
    RWMol *mol =
        SmilesToMol("C[C@@H]1[C@@H]2Cc3ccc(O)cc3[C@]1(C)CCN2CCN4CCCC4");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    params = DGeomHelpers::ETKDG;
    params.randomSeed = 0xf00d;
    cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);

    delete mol;
  }

  {
    // CHEMBL290986
    // std::cerr << "-----------------------------------" << std::endl;
    RWMol *mol = SmilesToMol(
        "COc1c(O)ccc2O[C@@H]([C@@H]3CCCC(=C3)C)c4c(ccc5NC(C)(C)C=C(C)c45)c12");
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);
    // std::cerr << "-----------------------------------" << std::endl;
    params = DGeomHelpers::ETKDG;
    params.randomSeed = 0xf00d;
    cid = DGeomHelpers::EmbedMolecule(*mol, params);
    TEST_ASSERT(cid >= 0);

    // std::cerr << MolToMolBlock(*mol) << std::endl;
    delete mol;
  }
}

void testGithubPullRequest1635() {
  {
    RWMol *m = SmilesToMol("C1(F)(F)CCC(CC1)COCC(C23CC4CC(C2)CC(C4)C3)N");
    TEST_ASSERT(m);
    MolOps::addHs(*m);
    const int expected_num_atoms = 54;
    TEST_ASSERT(m->getNumAtoms() == expected_num_atoms);

    RWMol firstMol(*m);
    RWMol secondMol(*m);
    delete m;

    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    params.randomSeed = MAX_INT;  // the largest possible random seed

    INT_VECT firstCids = DGeomHelpers::EmbedMultipleConfs(firstMol, 10, params);
    INT_VECT secondCids =
        DGeomHelpers::EmbedMultipleConfs(secondMol, 10, params);
    TEST_ASSERT(firstCids.size() == 10);
    TEST_ASSERT(secondCids.size() == 10);

    for (size_t i = 0; i < 10; i++) {
      TEST_ASSERT(firstCids[i] == secondCids[i]);

      int confIdx = firstCids[i];
      const Conformer &firstConf = firstMol.getConformer(confIdx);
      const Conformer &secondConf = secondMol.getConformer(confIdx);

      for (int atomIdx = 0; atomIdx < expected_num_atoms; ++atomIdx) {
        const RDGeom::Point3D &firstPoint = firstConf.getAtomPos(atomIdx);
        const RDGeom::Point3D &secondPoint = secondConf.getAtomPos(atomIdx);
        TEST_ASSERT(firstPoint.x == secondPoint.x);
        TEST_ASSERT(firstPoint.y == secondPoint.y);
        TEST_ASSERT(firstPoint.z == secondPoint.z);
      }
    }
  }
}

void testGithub1990() {
  boost::logging::disable_logs("rdApp.warning");
  {  // we saw the problem here (though it came from something in MolOps)
    std::unique_ptr<RWMol> mol(SmilesToMol("F/C=C/F"));
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    MolOps::removeHs(*mol);
    TEST_ASSERT(mol->getNumAtoms() == 4);
    int cid = DGeomHelpers::EmbedMolecule(*mol);
    TEST_ASSERT(cid >= 0);
  }
  {  // The original problem report
    std::unique_ptr<RWMol> mol(SmilesToMol(
        "CCCCCCCCCCCCCCCC(=O)O[C@@H]1CC(C)=C(/C=C/C(C)=C/C=C/C(C)=C/"
        "C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C2=C(C)C[C@@H](OC(=O)CCCCCCCCCCCCCCC)"
        "CC2(C)C)C(C)(C)C1"));
    TEST_ASSERT(mol);
    MolOps::addHs(*mol);
    MolOps::removeHs(*mol);
    int cid = DGeomHelpers::EmbedMolecule(*mol);
    TEST_ASSERT(cid >= 0);
  }
  boost::logging::enable_logs("rdApp.warning");
}

void testGithub2246() {
  {  // make sure the mechanics work
    std::vector<RDGeom::Point3D> pts = {{0, 0, 0}, {1.5, 0, 0}};
    auto m = "C1CC1C"_smiles;
    TEST_ASSERT(m);
    MolOps::addHs(*m);

    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    std::map<int, RDGeom::Point3D> coordMap;
    params.useRandomCoords = true;
    params.coordMap = &coordMap;
    params.maxIterations = 1;
    for (unsigned int i = 0; i < pts.size(); ++i) {
      coordMap[i] = pts[i];
    }
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*m, params);
    TEST_ASSERT(cid >= 0);
    for (unsigned int i = 0; i < pts.size(); ++i) {
      auto d = (m->getConformer().getAtomPos(i) - pts[i]).length();
      TEST_ASSERT(d < 1e-3);
    }
  }
  {  // a more complex example
    std::vector<RDGeom::Point3D> pts = {
        {0, 0, 0}, {1.5, 0, 0}, {1.5, 1.5, 0}, {0, 1.5, 0}};
    auto m = "C12C3CC1.O2C.C3CC"_smiles;
    TEST_ASSERT(m);
    MolOps::addHs(*m);
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    std::map<int, RDGeom::Point3D> coordMap;
    params.useRandomCoords = true;
    params.coordMap = &coordMap;
    params.maxIterations = 1;
    for (unsigned int i = 0; i < pts.size(); ++i) {
      coordMap[i] = pts[i];
    }
    for (unsigned int i = 0; i < 100; ++i) {
      params.randomSeed = i + 1;
      int cid = DGeomHelpers::EmbedMolecule(*m, params);
      TEST_ASSERT(cid >= 0);
      for (unsigned int i = 0; i < pts.size(); ++i) {
        auto d = (m->getConformer().getAtomPos(i) - pts[i]).length();
        TEST_ASSERT(d < 1e-3);
      }
    }
    MolOps::removeHs(*m);
  }
}

void testProvideBoundsMatrix() {
  boost::logging::disable_logs("rdApp.warning");
  {  // make sure the mechanics work
    auto m = "C1CCC1C"_smiles;
    TEST_ASSERT(m);
    auto nats = m->getNumAtoms();
    DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
    DGeomHelpers::initBoundsMat(mat);
    DGeomHelpers::setTopolBounds(*m, mat);

    // pick some silly bounds, just to make sure this works:
    mat->setUpperBound(3, 0, 1.21);
    mat->setLowerBound(3, 0, 1.2);
    mat->setUpperBound(3, 2, 1.21);
    mat->setLowerBound(3, 2, 1.2);
    mat->setUpperBound(3, 4, 1.21);
    mat->setLowerBound(3, 4, 1.2);
    DistGeom::triangleSmoothBounds(mat);

    DGeomHelpers::EmbedParameters params;
    params.useRandomCoords = true;
    params.boundsMat = mat;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*m, params);
    TEST_ASSERT(cid >= 0);

    const auto conf = m->getConformer(cid);
    TEST_ASSERT(
        feq((conf.getAtomPos(3) - conf.getAtomPos(0)).length(), 1.2, 0.05));
    TEST_ASSERT(
        feq((conf.getAtomPos(3) - conf.getAtomPos(2)).length(), 1.2, 0.05));
    TEST_ASSERT(
        feq((conf.getAtomPos(3) - conf.getAtomPos(4)).length(), 1.2, 0.05));
  }
  boost::logging::enable_logs("rdApp.warning");
}

void testDisableFragmentation() {
  {  // make sure the mechanics work
    auto m = "OO.OO"_smiles;
    TEST_ASSERT(m);
    MolOps::addHs(*m);
    DGeomHelpers::EmbedParameters params;
    params.embedFragmentsSeparately = false;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*m, params);
    TEST_ASSERT(cid >= 0);

    const auto conf = m->getConformer(cid);

    TEST_ASSERT((conf.getAtomPos(0) - conf.getAtomPos(2)).length() > 2.0);
    TEST_ASSERT((conf.getAtomPos(0) - conf.getAtomPos(3)).length() > 2.0);
    TEST_ASSERT((conf.getAtomPos(1) - conf.getAtomPos(2)).length() > 2.0);
    TEST_ASSERT((conf.getAtomPos(1) - conf.getAtomPos(3)).length() > 2.0);
  }
}

void testGithub3019() {
  {  // make sure the mechanics work
    std::unique_ptr<RWMol> m(SmilesToMol(std::string(2000, 'C')));
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2000);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*m, params);
    TEST_ASSERT(cid >= 0);
  }
}

namespace {
void throwerror(unsigned int) {
  throw ValueErrorException("embedder is abortable");
}
}  // namespace

void testGithub3667() {
  auto *mol = SmilesToMol(
      "c12c3c4c5c6c1c1c7c8c9c%10c%11c(c28)c3c2c3c4c4c5c5c8c6c1c1c6c7c9c7c9c%"
      "10c%10c%11c2c2c3c3c4c4c5c5c%11c%12c(c1c85)c6c7c1c%12c5c%11c4c3c3c5c(c91)"
      "c%10c23");
  TEST_ASSERT(mol);

  bool ok = false;
  try {
    DGeomHelpers::EmbedParameters params;
    params.callback = throwerror;
    DGeomHelpers::EmbedMolecule(*mol, params);
    ok = false;
  } catch (const ValueErrorException &) {
    ok = true;
  }
  TEST_ASSERT(ok);
  delete mol;
}

void testForceTransAmides() {
  auto mol = "CC(=O)NC"_smiles;
  TEST_ASSERT(mol);
  bool updateLabel = true;
  bool takeOwnership = true;
  mol->addAtom(new Atom(1), updateLabel, takeOwnership);
  mol->addBond(3, 5, Bond::BondType::SINGLE);
  MolOps::sanitizeMol(*mol);
  MolOps::addHs(*mol);
#if 0
  // worth leaving this here just to allow looking at the conformers if anything goes wrong later
  {
    DGeomHelpers::EmbedParameters params;
    params.forceTransAmides = true;
    params.randomSeed = 0xf00d;
    params.useExpTorsionAnglePrefs = false;
    params.useBasicKnowledge = true;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, params);
    SDWriter w("amide.sdf");
    for (auto cid : cids) {
      TEST_ASSERT(cid >= 0);
      auto conf = mol->getConformer(cid);
      auto tors = MolTransforms::getDihedralDeg(conf, 0, 1, 3, 4);
      if (fabs(fabs(tors) - 180) > 5) {
        w.write(*mol, cid);
        std::cerr << cid << " TORS: " << tors << std::endl;
        if (fabs(fabs(tors) - 180) > 40) {
          std::cerr << "---------- DM " << std::endl;
          double *dm = MolOps::get3DDistanceMat(*mol, cid);
          auto nAtoms = mol->getNumAtoms();
          for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = 0; j < nAtoms; ++j) {
              std::cerr << " " << std::setprecision(3) << std::setw(5)
                        << dm[i * nAtoms + j];
            }
            std::cerr << std::endl;
          }
        }
      }
      // TEST_ASSERT(fabs(fabs(tors) - 180) < 5);
    }
    w.flush();
  }
#endif
  {
    DGeomHelpers::EmbedParameters params;
    params.forceTransAmides = true;
    params.randomSeed = 0xf00d;
    params.useExpTorsionAnglePrefs = false;
    params.useBasicKnowledge = true;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, params);
    for (auto cid : cids) {
      TEST_ASSERT(cid >= 0);
      auto conf = mol->getConformer(cid);
      auto tors = MolTransforms::getDihedralDeg(conf, 0, 1, 3, 4);
      TEST_ASSERT(fabs(fabs(tors) - 180) < 37);
      tors = MolTransforms::getDihedralDeg(conf, 2, 1, 3, 5);
      TEST_ASSERT(fabs(fabs(tors) - 180) < 37);
    }
  }
  {  // make sure we can find at least one non-trans
    DGeomHelpers::EmbedParameters params;
    params.forceTransAmides = false;
    params.randomSeed = 0xf00d;
    params.useExpTorsionAnglePrefs = false;
    params.useBasicKnowledge = true;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, params);
    bool foundOne = false;
    for (auto cid : cids) {
      TEST_ASSERT(cid >= 0);
      auto conf = mol->getConformer(cid);
      auto tors = MolTransforms::getDihedralDeg(conf, 0, 1, 3, 4);
      if (fabs(fabs(tors) - 180) > 50) {
        foundOne = true;
        break;
      }
    }
    TEST_ASSERT(foundOne);
  }
}

void testSymmetryPruning() {
  auto mol = "CCOC(C)(C)C"_smiles;
  TEST_ASSERT(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters params;
  params.useSymmetryForPruning = true;
  params.onlyHeavyAtomsForRMS = true;
  params.pruneRmsThresh = 0.5;
  params.randomSeed = 0xf00d;
  auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 50, params);
  TEST_ASSERT(cids.size() == 1);

  params.useSymmetryForPruning = false;
  cids = DGeomHelpers::EmbedMultipleConfs(*mol, 50, params);
  TEST_ASSERT(cids.size() == 3);
}

void testMissingHsWarning() {
  auto mol = "CC"_smiles;
  TEST_ASSERT(mol);

  std::stringstream ss;
  rdWarningLog->SetTee(ss);
  DGeomHelpers::EmbedParameters params;
  DGeomHelpers::EmbedMolecule(*mol, params);
  rdWarningLog->ClearTee();
  TEST_ASSERT(ss.str().find("Molecule does not have explicit Hs") !=
              std::string::npos);
}

void testHydrogenBondBasics() {
  auto mol = "CC1O[H]O=C(C)C1 |H:4.3|"_smiles;
  TEST_ASSERT(mol);
  MolOps::addHs(*mol);

  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(mol->getNumAtoms()));
  DGeomHelpers::initBoundsMat(mat.get());
  DGeomHelpers::setTopolBounds(*mol, mat);
  DistGeom::triangleSmoothBounds(mat.get());
  TEST_ASSERT(mat->getVal(3, 4) > 1.8);
  TEST_ASSERT(mat->getVal(3, 4) < 2.2);
  TEST_ASSERT(mat->getVal(4, 3) > 1.0);
  TEST_ASSERT(mat->getVal(4, 3) < 1.5);

  DGeomHelpers::EmbedParameters params = DGeomHelpers::ETKDGv3;
  params.randomSeed = 0xf00d;
  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
  auto dist = MolTransforms::getBondLength(mol->getConformer(), 3, 4);
  TEST_ASSERT(dist < 1.5);
}

int main() {
  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing DistGeomHelpers\n";

#if 1
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test2 \n\n";
  test2();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test3 \n\n";
  test3();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test4 \n\n";
  test4();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test5 \n\n";
  test5();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test6 \n\n";
  test6();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test15Dists \n\n";
  test15Dists();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test Issue 215 \n\n";
  testIssue215();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testMultipleConfs \n\n";
  testMultipleConfs();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testMultipleConfsExpTors \n\n";
  testMultipleConfsExpTors();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue227 \n\n";
  testIssue227();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue236 \n\n";
  testIssue236();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testOrdering \n\n";
  testOrdering();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue244 \n\n";
  testIssue244();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue251 \n\n";
  testIssue251();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue276 \n";
  testIssue276();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue284 \n\n";
  testIssue284();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue285 \n\n";
  testIssue285();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue355 \n\n";
  testIssue355();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testRandomCoords \n\n";
  testRandomCoords();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 1989539 \n\n";
  testIssue1989539();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 2091864 \n\n";
  testIssue2091864();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test constrained embedding \n\n";
  testConstrainedEmbedding();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 2091974 \n\n";
  testIssue2091974();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 2835784 \n\n";
  testIssue2835784();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 3019283 \n\n";
  testIssue3019283();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 3238580 \n\n";
  testIssue3238580();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 3483968 \n\n";
  testIssue3483968();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test github issue 55 \n\n";
  testGithub55();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test1 \n\n";
  test1();
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test multi-threading \n\n";
  testMultiThread();
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog)
      << "\t test github issue 256: handling of zero-atom molecules\n\n";
  testGithub256();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test embedder callback function \n";

  testGithub3667();

#ifdef RDK_TEST_MULTITHREADED
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test multi-threaded multi-conf embedding \n\n";
  testMultiThreadMultiConf();
#endif

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test github issue 563: Incorrect ring "
                          "stereochemistry after embedding\n\n";
  testGithub563();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test github issue 568: Incorrect stereochemistry "
                          "after embedding\n\n";
  testGithub568();
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test github issue 696: Bad 1-4 bounds matrix "
                          "elements in highly constrained system\n";
  testGithub696();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog)
      << "\t More ChEMBL molecules failing bounds smoothing.\n";
  testGithub697();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t ugly conformations can be generated for highly "
                          "constrained ring systems.\n";
  testGithub971();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test embed parameters structure.\n";
  testEmbedParameters();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog)
      << "\t test github 1227: Hs being used in RMSD filtering.\n";
  testGithub1227();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Failure to embed larger aromatic rings.\n";
  testGithub1240();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Deterministic with large random seeds\n";
  testGithubPullRequest1635();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Github #1990: seg fault after RemoveHs\n";
  testGithub1990();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Github #2246: Use coordMap when starting "
                          "embedding from random coords\n";
  testGithub2246();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Providing a distance bounds matrix.\n";
  testProvideBoundsMatrix();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Disabling fragmentation.\n";
  testDisableFragmentation();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Using symmetry in conformation pruning.\n";
  testSymmetryPruning();

#endif
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Force trans amides.\n";
  testForceTransAmides();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test missing Hs warning.\n";
  testMissingHsWarning();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t basic H bond support.\n";
  testHydrogenBondBasics();

#ifdef EXECUTE_LONG_TESTS
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog)
      << "\t Github #3019: Seg fault for very large molecules.\n";
  testGithub3019();
#endif

  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";

  return (0);
}
