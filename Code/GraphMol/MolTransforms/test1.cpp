//
//   Copyright (C) 2003-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <Geometry/point.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <iostream>
#include <algorithm>

using namespace RDKit;
using namespace MolTransforms;
bool comparePts(const RDGeom::Point3D &pt1, const RDGeom::Point3D &pt2,
                double tol = 1.0e-4) {
  RDGeom::Point3D tpt = pt1;
  tpt -= pt2;
  return (tpt.length() < tol);
}

void test1Canonicalization() {
  ROMol *mol = SmilesToMol("C", 0, 1);
  auto *conf = new Conformer(1);
  conf->setAtomPos(0, RDGeom::Point3D(4.0, 5.0, 6.0));
  int cid = mol->addConformer(conf, true);
  CHECK_INVARIANT(cid >= 0, "")
  RDGeom::Point3D pt = computeCentroid(*conf);
  CHECK_INVARIANT(comparePts(pt, RDGeom::Point3D(4.0, 5.0, 6.0)), "");

  RDGeom::Transform3D *trans = computeCanonicalTransform(*conf);
  transformConformer(*conf, *trans);
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(0), RDGeom::Point3D(0.0, 0.0, 0.0)), "");

  conf->setAtomPos(0, RDGeom::Point3D(4.0, 5.0, 6.0));
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(0), RDGeom::Point3D(0.0, 0.0, 0.0)), "");

  delete mol;
  // delete conf;
  delete trans;
  // lets try two points now
  mol = SmilesToMol("CC", 0, 1);
  conf = new Conformer(2);
  conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
  cid = mol->addConformer(conf, true);
  trans = computeCanonicalTransform(*conf);
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");

  conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(0.0, 1.5, 0.0));
  delete trans;
  trans = computeCanonicalTransform(*conf);
  canonicalizeConformer(*conf);

  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");
  delete mol;
  delete trans;

  mol = SmilesToMol("CC", 0, 1);
  conf = new Conformer(2);
  conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
  cid = mol->addConformer(conf, true);
  trans = computeCanonicalTransform(*conf);
  transformConformer(*conf, *trans);
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  CHECK_INVARIANT(
      comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");
  delete mol;
  delete trans;

  mol = SmilesToMol("C1CC1", 0, 1);
  conf = new Conformer(3);
  conf->setAtomPos(0, RDGeom::Point3D(0.58, -0.66, -0.08));
  conf->setAtomPos(1, RDGeom::Point3D(-0.88, -0.18, -0.04));
  conf->setAtomPos(2, RDGeom::Point3D(.26, 0.82, 0.14));
  cid = mol->addConformer(conf, true);
  // trans = computeCanonicalTransform(*conf);
  // transformConformer(*conf, *trans);
  canonicalizeConformer(*conf);

// computeCanonicalTransform returns more approximate eigenvalues/eigencvectors
// when built against the native RDKit PowerEigenSolver, so unit test results
// differ slightly
#ifdef RDK_HAS_EIGEN3
  std::vector<RDGeom::Point3D> expected = {
      RDGeom::Point3D(0.8244, -0.3268, 0.0),
      RDGeom::Point3D(-0.6975, -0.5449, 0.0),
      RDGeom::Point3D(-0.1269, 0.8716, 0.0)};
#else
  std::vector<RDGeom::Point3D> expected = {
      RDGeom::Point3D(-0.6418, 0.6158, 0.0),
      RDGeom::Point3D(-0.2029, -0.8602, 0.0),
      RDGeom::Point3D(0.8447, 0.2445, 0.0)};
#endif
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), expected.at(0)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(1), expected.at(1)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(2), expected.at(2)), "");
  MolToMolFile(*mol, "junk.mol", true, 0);
  delete mol;

  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/MolTransforms/test_data/1oir.mol";
  mol = MolFileToMol(fname1);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolTransforms/test_data/1oir_canon.mol";
  ROMol *mol2 = MolFileToMol(fname2);

  Conformer &conf1 = mol->getConformer(0);
  canonicalizeConformer(conf1);

  Conformer &conf2 = mol2->getConformer();
  unsigned int i, nats = mol->getNumAtoms();
  for (i = 0; i < nats; ++i) {
    CHECK_INVARIANT(comparePts(conf1.getAtomPos(i), conf2.getAtomPos(i)), "");
  }

  delete mol;
  delete mol2;
}

void test1() {
  std::cout << " ----------> Test1 " << std::endl;

  std::cout << " Finished <---------- " << std::endl;
}

void testGetSetBondLength() {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase +
      "/Code/GraphMol/MolTransforms/test_data/3-cyclohexylpyridine.mol";
  RWMol *m = MolFileToMol(fName, true, false);
  TEST_ASSERT(m);
  Conformer &conf = m->getConformer();
  double dist = getBondLength(conf, 0, 19);
  TEST_ASSERT(RDKit::feq(dist, 1.36));
  setBondLength(conf, 0, 19, 2.5);
  dist = getBondLength(conf, 0, 19);
  TEST_ASSERT(RDKit::feq(dist, 2.5));
  setBondLength(conf, 19, 0, 3.0);
  dist = getBondLength(conf, 0, 19);
  TEST_ASSERT(RDKit::feq(dist, 3.0));
  delete m;
}

void testGetSetAngle() {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase +
      "/Code/GraphMol/MolTransforms/test_data/3-cyclohexylpyridine.mol";
  RWMol *m = MolFileToMol(fName, true, false);
  TEST_ASSERT(m);
  Conformer &conf = m->getConformer();
  double angle = getAngleDeg(conf, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(angle, 109.7, 0.05));
  setAngleDeg(conf, 0, 19, 21, 125.0);
  angle = getAngleDeg(conf, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(angle, 125.0));
  setAngleRad(conf, 21, 19, 0, M_PI / 2.);
  angle = getAngleRad(conf, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(angle, M_PI / 2.));
  angle = getAngleDeg(conf, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(angle, 90.0));
  delete m;
}

void testGetSetDihedral() {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase +
      "/Code/GraphMol/MolTransforms/test_data/3-cyclohexylpyridine.mol";
  RWMol *m = MolFileToMol(fName, true, false);
  TEST_ASSERT(m);
  Conformer &conf = m->getConformer();
  double dihedral = getDihedralDeg(conf, 0, 19, 21, 24);
  TEST_ASSERT(RDKit::feq(dihedral, 176.05, 0.05));
  setDihedralDeg(conf, 8, 0, 19, 21, 65.0);
  dihedral = getDihedralDeg(conf, 8, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(dihedral, 65.0));
  setDihedralDeg(conf, 8, 0, 19, 21, -130.0);
  dihedral = getDihedralDeg(conf, 8, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(dihedral, -130.0));
  setDihedralRad(conf, 21, 19, 0, 8, -2. / 3. * M_PI);
  dihedral = getDihedralRad(conf, 8, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(dihedral, -2. / 3. * M_PI));
  dihedral = getDihedralDeg(conf, 8, 0, 19, 21);
  TEST_ASSERT(RDKit::feq(dihedral, -120.0));
  delete m;
}

void testGetSetDihedralThroughTripleBond() {
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MolTransforms/test_data/github1262_2.mol";
  RWMol *m = MolFileToMol(fName, true, false);
  TEST_ASSERT(m);
  Conformer &conf = m->getConformer();
  setDihedralDeg(conf, 6, 1, 2, 9, 0.0);
  double dihedral = getDihedralDeg(conf, 6, 1, 2, 9);
  TEST_ASSERT(RDKit::feq(dihedral, 0.0));
  double dist = getBondLength(conf, 6, 9);
  setDihedralDeg(conf, 6, 1, 2, 9, 120.0);
  dihedral = getDihedralDeg(conf, 6, 1, 2, 9);
  TEST_ASSERT(RDKit::feq(dihedral, 120.0));
  double dist2 = getBondLength(conf, 6, 7);
  TEST_ASSERT(RDKit::feq(dist, dist2, 0.05));
  setDihedralDeg(conf, 6, 1, 2, 9, 180.0);
  dihedral = getDihedralDeg(conf, 6, 1, 2, 9);
  TEST_ASSERT(RDKit::feq(dihedral, 180.0));
  double dist3 = getBondLength(conf, 6, 9);
  TEST_ASSERT(!RDKit::feq(dist, dist3, 0.3));
  bool exceptionRaised = false;
  try {
    setDihedralDeg(conf, 6, 0, 3, 9, 0.0);
  } catch (ValueErrorException &) {
    exceptionRaised = true;
  }
  TEST_ASSERT(exceptionRaised);
  delete m;
}

#ifndef RDK_HAS_EIGEN3
void testGithub1262() {}
#else
void _calcAxesAndMoments(RWMol *m, Eigen::Matrix3d &axes,
                         Eigen::Vector3d &moments) {
  TEST_ASSERT(m);
  Conformer &conf = m->getConformer();
  std::vector<double> weights;
  weights.resize(m->getNumAtoms());
  for (ROMol::AtomIterator cai = m->beginAtoms(); cai != m->endAtoms(); ++cai) {
    weights[(*cai)->getIdx()] = (*cai)->getMass();
  }

  bool ignoreHs = false, force = true;
  computePrincipalAxesAndMoments(conf, axes, moments, ignoreHs, force,
                                 &weights);
}

void testGithub1262() {
  std::string rdbase = getenv("RDBASE");
  {  // a disc (benzene)
    std::string fName =
        rdbase + "/Code/GraphMol/MolTransforms/test_data/github1262_1.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    Eigen::Matrix3d axes;
    Eigen::Vector3d moments;
    _calcAxesAndMoments(m, axes, moments);
    TEST_ASSERT((moments(2) - moments(0)) > 10.);
    TEST_ASSERT((moments(2) - moments(1)) > 10.);
    TEST_ASSERT((moments(1) - moments(0)) < 1e-2);

    delete m;
  }
  {  // a rod
    std::string fName =
        rdbase + "/Code/GraphMol/MolTransforms/test_data/github1262_2.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    Eigen::Matrix3d axes;
    Eigen::Vector3d moments;
    _calcAxesAndMoments(m, axes, moments);
    TEST_ASSERT((moments(2) - moments(0)) > 10.);
    TEST_ASSERT((moments(2) - moments(1)) < 1e-2);
    TEST_ASSERT((moments(1) - moments(0)) > 10);

    delete m;
  }
  {  // adamantane
    std::string fName =
        rdbase + "/Code/GraphMol/MolTransforms/test_data/github1262_3.mol";
    RWMol *m = MolFileToMol(fName, true, false);
    TEST_ASSERT(m);
    Eigen::Matrix3d axes;
    Eigen::Vector3d moments;
    _calcAxesAndMoments(m, axes, moments);
    TEST_ASSERT((moments(2) - moments(0)) < 1e-2);
    TEST_ASSERT((moments(2) - moments(1)) < 1e-2);
    TEST_ASSERT((moments(1) - moments(0)) < 1e-2);

    delete m;
  }
}
#endif

void testGithub1908() {
  std::string rdbase = getenv("RDBASE");
  {  // a disc (benzene)
    std::string fName =
        rdbase + "/Code/GraphMol/MolTransforms/test_data/github1908_2.mol";
    std::unique_ptr<RWMol> m(MolFileToMol(fName));
    TEST_ASSERT(m);

    Conformer &conf = m->getConformer();
    double dist = getBondLength(conf, 0, 1);
    // std::cerr << " 1: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));
    dist = getBondLength(conf, 1, 2);
    // std::cerr << " 2: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));

    canonicalizeConformer(conf);

    dist = getBondLength(conf, 0, 1);
    // std::cerr << " 3: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));
    dist = getBondLength(conf, 1, 2);
    // std::cerr << " 4: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));
  }
  {  // a disc (benzene)
    std::string fName =
        rdbase + "/Code/GraphMol/MolTransforms/test_data/github1908_1.mol";
    std::unique_ptr<RWMol> m(MolFileToMol(fName));
    TEST_ASSERT(m);

    Conformer &conf = m->getConformer();
    double dist = getBondLength(conf, 0, 1);
    // std::cerr << " 1: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));
    dist = getBondLength(conf, 1, 2);
    // std::cerr << " 2: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));

    canonicalizeConformer(conf);

    dist = getBondLength(conf, 0, 1);
    // std::cerr << " 3: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));
    dist = getBondLength(conf, 1, 2);
    // std::cerr << " 4: " << dist << std::endl;
    TEST_ASSERT(feq(dist, 1.38, .02));
  }
}

void testGithub4302() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/MolTransforms/test_data/github4302.sdf";
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolTransforms/test_data/github4302_canon.sdf";
  RDKit::SDMolSupplier reader(fname1);
  RDKit::SDWriter writer(fname2);
  while (!reader.atEnd()) {
    std::unique_ptr<RDKit::ROMol> mol(reader.next());
    const RDKit::Conformer &conf = mol->getConformer();
    auto canonConf = new RDKit::Conformer(conf);
    auto cid = mol->addConformer(canonConf, true);
    canonicalizeConformer(*canonConf);
    // the native RDKit eigensolver comes up with non-canonical
    // or distorted coordinates with these conformations
#ifdef RDK_HAS_EIGEN3
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
      TEST_ASSERT(comparePts(canonConf->getAtomPos(i), conf.getAtomPos(i)));
    }
#endif
    writer.write(*mol, cid);
  }
  writer.close();
}

void testWeightedCentroid() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/MolTransforms/test_data/github4302.sdf";
  RDKit::SDMolSupplier reader(fname1);
  while (!reader.atEnd()) {
    std::unique_ptr<RDKit::ROMol> mol(reader.next());
    const RDKit::Conformer &conf = mol->getConformer();
    std::vector<double> weights;
    weights.reserve(mol->getNumAtoms());
    for (const auto a : mol->atoms()) {
      weights.push_back(
          PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum()));
    }
    RDGeom::Point3D ctd;
    double totalMass = 0.0;
    for (const auto a : mol->atoms()) {
      auto atomicMass =
          PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum());
      totalMass += atomicMass;
      ctd += conf.getAtomPos(a->getIdx()) * atomicMass;
    }
    ctd /= totalMass;
    TEST_ASSERT((ctd - computeCentroid(conf, true, &weights)).length() < 1.e-4);
  }
}

int main() {
  // test1();
  std::cout << "***********************************************************\n";
  std::cout << "Testing MolTransforms\n";

#if 1
  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1Canonicalization \n\n";
  test1Canonicalization();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGetSetBondLength \n\n";
  testGetSetBondLength();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGetSetAngle \n\n";
  testGetSetAngle();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGetSetDihedral \n\n";
  testGetSetDihedral();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGetSetDihedralThroughTripleBond \n\n";
  testGetSetDihedralThroughTripleBond();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGithub1262: PMI descriptors incorrect  \n\n";
  testGithub1262();
#endif
  std::cout << "\t---------------------------------\n";
  std::cout
      << "\t testGithub1908: CanonicalizeMol() distorting bond lengths\n\n";
  testGithub1908();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testGithub4302: native computeCanonicalTransform() "
               "generates non-canonical coords\n\n";
  testGithub4302();
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testWeightedCentroid \n\n";
  testWeightedCentroid();
  std::cout << "***********************************************************\n";
  return (0);
}
