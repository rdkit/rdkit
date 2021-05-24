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
  int nDone = 0;
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
    ++nDone;
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

int main() {
  RDLog::InitLogs();
  test1();
  testPBFEdges();
}
