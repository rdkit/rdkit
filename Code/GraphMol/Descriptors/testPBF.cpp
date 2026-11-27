//
//  Copyright (C) 2012-2016 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <memory>

#include <GraphMol/Descriptors/PBF.h>

void test1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic PBF tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::string fName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.out";
  std::ifstream instrm(fName.c_str());
  while (!reader.atEnd()) {
    RDKit::ROMol *m = reader.next();
    TEST_ASSERT(m);
    std::string nm;
    m->getProp("_Name", nm);
    double dpbf = RDKit::Descriptors::PBF(*m);

    std::string inm;
    double ref;
    instrm >> inm;
    instrm >> ref;
    TEST_ASSERT(inm == nm);
    if (fabs(ref - dpbf) > .001) {
      std::cerr << "value mismatch: " << inm << " " << ref << " " << dpbf
                << std::endl;
    }
    TEST_ASSERT(fabs(ref - dpbf) < 0.001);
    delete m;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testPBFEdges() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    PBF edge cases." << std::endl;

  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/linear.mol";

    RDKit::ROMol *m = RDKit::MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double dpbf = RDKit::Descriptors::PBF(*m);
    TEST_ASSERT(dpbf <= 1e-4);
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/linear_2atom.mol";

    RDKit::ROMol *m = RDKit::MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double dpbf = RDKit::Descriptors::PBF(*m);
    TEST_ASSERT(dpbf <= 1e-4);
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/planar.mol";

    RDKit::ROMol *m = RDKit::MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double dpbf = RDKit::Descriptors::PBF(*m);
    TEST_ASSERT(dpbf <= 1e-4);
    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName =
        pathName + "/Code/GraphMol/Descriptors/test_data/planar_3atom.mol";

    RDKit::ROMol *m = RDKit::MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double dpbf = RDKit::Descriptors::PBF(*m);
    TEST_ASSERT(dpbf <= 1e-4);
    delete m;
  }
  {
    RDKit::RWMol m;
    bool updateLabel = true;
    bool takeOwnership = true;
    m.addAtom(new RDKit::Atom(6), updateLabel, takeOwnership);
    m.addAtom(new RDKit::Atom(6), updateLabel, takeOwnership);
    m.addAtom(new RDKit::Atom(6), updateLabel, takeOwnership);
    m.addAtom(new RDKit::Atom(6), updateLabel, takeOwnership);
    m.addAtom(new RDKit::Atom(6), updateLabel, takeOwnership);
    m.addAtom(new RDKit::Atom(6), updateLabel, takeOwnership);
    m.addConformer(new RDKit::Conformer(m.getNumAtoms()));
    double dpbf = RDKit::Descriptors::PBF(m);
    TEST_ASSERT(dpbf <= 1e-4);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testPBFMultipleConformers() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    PBF with multiple conformers." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName =
      pathName + "/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName, true, false);
  std::unique_ptr<RDKit::ROMol> m(reader.next());
  TEST_ASSERT(m);
  TEST_ASSERT(m->getNumConformers() == 1);
  double refPBF = RDKit::Descriptors::PBF(*m);
  TEST_ASSERT(refPBF > 0.1);

  // add rotated copies of the conformer: PBF is rotation-invariant, so each
  // must give the same value as the original. The principal axes from the
  // first conformer must not be reused for the others.
  const double c = std::cos(0.5);
  const double s = std::sin(0.5);
  auto *conf1 = new RDKit::Conformer(m->getConformer());
  auto *conf2 = new RDKit::Conformer(m->getConformer());
  for (unsigned int i = 0; i < m->getNumAtoms(); ++i) {
    const auto &p = m->getConformer().getAtomPos(i);
    conf1->setAtomPos(i, RDGeom::Point3D(p.y, p.z, p.x));
    conf2->setAtomPos(
        i, RDGeom::Point3D(c * p.x - s * p.y, s * p.x + c * p.y, p.z));
  }
  auto cid1 = m->addConformer(conf1, true);
  auto cid2 = m->addConformer(conf2, true);

  m->clearComputedProps();
  for (auto cid : {0, static_cast<int>(cid1), static_cast<int>(cid2)}) {
    double dpbf = RDKit::Descriptors::PBF(*m, cid);
    TEST_ASSERT(fabs(dpbf - refPBF) < 1e-4);
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main() {
  RDLog::InitLogs();
  test1();
  testPBFEdges();
  testPBFMultipleConformers();
}
