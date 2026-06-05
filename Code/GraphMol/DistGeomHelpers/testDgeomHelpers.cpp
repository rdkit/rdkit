//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#include <RDGeneral/RDLog.h>
#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <algorithm>
#include <map>
#include <sstream>
#include <memory>
#include <vector>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
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
using tokenizer = boost::tokenizer<boost::char_separator<char>>;
using namespace RDKit;

namespace {
/**
 * @brief Computes pairwise distance matrix for 3D coordinates.
 * @param origCoords Input 3D coordinates.
 * @param distMat Output distance matrix.
 */
void computeDistMat(const RDGeom::PointPtrVect &origCoords,
                    RDNumeric::DoubleSymmMatrix &distMat) {
  unsigned int N = origCoords.size();
  REQUIRE(N == distMat.numRows());
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

/**
 * @brief Computes distance matrix for a molecule's first conformer.
 * @param mol Input molecule.
 * @param distMat Output distance matrix.
 */
void computeMolDmat(ROMol &mol, RDNumeric::DoubleSymmMatrix &distMat) {
  RDGeom::PointPtrVect origCoords;
  std::size_t nat = mol.getNumAtoms();
  Conformer &conf = mol.getConformer(0);
  for (std::size_t i = 0; i < nat; i++) {
    origCoords.push_back(&conf.getAtomPos(i));
  }
  computeDistMat(origCoords, distMat);
}

/**
 * @brief Generates distance geometry bounds matrix for a molecule.
 * @param m Input molecule pointer.
 * @param setTopolBounds If true, sets topological bounds.
 * @return Bounds matrix pointer.
 */
DistGeom::BoundsMatPtr _getBoundsMatrix(const std::unique_ptr<RWMol> &m,
                                        const bool setTopolBounds = true) {
  DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
  DGeomHelpers::initBoundsMat(bm);
  REQUIRE(bm);
  if (setTopolBounds) {
    DGeomHelpers::setTopolBounds(*m, bm);
    REQUIRE(bm);
  }
  return bm;
}

/**
 * @brief Generates distance geometry bounds matrix from a SMILES string.
 * @param smiles Input SMILES string.
 * @param setTopolBounds If true, sets topological bounds.
 * @return Bounds matrix pointer.
 */
DistGeom::BoundsMatPtr _getBoundsMatrix(const std::string &smiles,
                                        const bool setTopolBounds = true) {
  const auto m = v2::SmilesParse::MolFromSmiles(smiles);
  REQUIRE(m);
  return _getBoundsMatrix(m, setTopolBounds);
}

/**
 * @brief Asserts two molecule conformations are geometrically equivalent.
 * @param m Molecule to test.
 * @param expected Reference molecule.
 * @param molConfId Test conformer ID.
 * @param expectedConfId Reference conformer ID.
 */
void compareConfs(const RWMol *m, const RWMol *expected, int molConfId = -1,
                  int expectedConfId = -1) {
  PRECONDITION(m, "bad pointer");
  PRECONDITION(expected, "bad pointer");
  REQUIRE(m->getNumAtoms() == expected->getNumAtoms());
  const Conformer &conf1 = m->getConformer(molConfId);
  const Conformer &conf2 = expected->getConformer(expectedConfId);
  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    REQUIRE(m->getAtomWithIdx(i)->getAtomicNum() ==
            expected->getAtomWithIdx(i)->getAtomicNum());

    RDGeom::Point3D pt1i = conf1.getAtomPos(i);
    RDGeom::Point3D pt2i = conf2.getAtomPos(i);
    CHECK((pt1i - pt2i).length() < 0.05);
  }
}
}  // namespace

const std::string rdbase = getenv("RDBASE");

TEST_CASE("test1") {
  boost::logging::disable_logs("rdApp.warning");
  std::vector<std::string> smiString = {
      "CC1=C(C(C)=CC=C2)C2=CC=C1",
      "c1ccccc1C",
      "C/C=C/CC",
      "C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4",
      "C1CCCCS1(=O)(=O)",
      "c1ccccc1",
      "C1CCCC1",
      "C1CCCCC1",
      "C1CC1(C)C",
      "C12(C)CC1CC2"};
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/initCoords.sdf";
  v2::FileParsers::SDMolSupplier sdsup(fname);
  // SDWriter writer("foo.sdf");
  for (std::size_t i = 0; i < smiString.size(); ++i) {
    auto m = v2::SmilesParse::MolFromSmiles(smiString[i]);
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, true, false, 2, true, 1,
                                          nullptr, 1e-2);
    REQUIRE(cid >= 0);
    auto m2 = sdsup.next();
    // BOOST_LOG(rdInfoLog) << ">>> " << smi << std::endl;
    // writer.write(*m);
    // writer.flush();

    // ROMol *m2 = NULL;
    if (m2) {
      unsigned int nat = m->getNumAtoms();

      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
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
            CHECK(fabs(d1 - d2) / d1 < 0.06);
          } else {
            // BOOST_LOG(rdInfoLog) << ">2> " <<i<<","<<j<<":"<< d1 << " " << d2
            // << " "<<fabs(d1-d2)/d1<<std::endl;
            CHECK(fabs(d1 - d2) / d1 < 0.12);
          }
        }
      }
    }
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("test2") {
  auto getResults = [](const std::string &smiles) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    auto bm = _getBoundsMatrix(mol, true);
    auto cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
    REQUIRE(cid >= 0);
    auto dmat =
        std::make_unique<RDNumeric::DoubleSymmMatrix>(mol->getNumAtoms(), 0.0);
    computeMolDmat(*mol, *dmat);
    return std::make_tuple(std::move(mol), std::move(dmat), std::move(bm));
  };

  boost::logging::disable_logs("rdApp.warning");
  SECTION("Ring Distances") {
    const std::string smiles = "Cc1c(C=CC(C)N2)c2[nH]n1";
    auto [mol, dmat, bm] = getResults(smiles);

    CHECK((bm->getUpperBound(0, 9) - bm->getLowerBound(0, 9)) < 0.13);
    CHECK(bm->getUpperBound(0, 9) - dmat->getVal(0, 9) > -0.1);
    CHECK(bm->getLowerBound(0, 9) - dmat->getVal(0, 9) < 0.10);

    CHECK((bm->getUpperBound(10, 7) - bm->getLowerBound(10, 7)) < 0.13);
    CHECK(bm->getUpperBound(10, 7) - dmat->getVal(10, 7) > -0.1);
    CHECK(bm->getLowerBound(10, 7) - dmat->getVal(10, 7) < 0.10);

    CHECK((bm->getUpperBound(2, 5) - bm->getLowerBound(2, 5)) < 0.20);
    CHECK(bm->getUpperBound(2, 5) - dmat->getVal(2, 5) > -0.1);
    CHECK(bm->getLowerBound(2, 5) - dmat->getVal(2, 5) < 0.10);

    CHECK((bm->getUpperBound(8, 4) - bm->getLowerBound(8, 4)) > 1.);
    CHECK((bm->getUpperBound(8, 4) - bm->getLowerBound(8, 4)) < 1.2);
    CHECK((bm->getUpperBound(8, 4) - dmat->getVal(8, 4) > -0.1));
    CHECK((bm->getLowerBound(8, 4) - dmat->getVal(8, 4) < 0.10));

    CHECK((bm->getUpperBound(8, 6) - bm->getLowerBound(8, 6)) > 1.0);
    CHECK((bm->getUpperBound(8, 6) - bm->getLowerBound(8, 6)) < 1.2);
    CHECK((bm->getUpperBound(8, 6) - dmat->getVal(8, 6) > -0.1));
    CHECK((bm->getLowerBound(8, 6) - dmat->getVal(8, 6) < 0.10));
  }

  SECTION("chain of singles") {
    const std::string smiles = "CCCC";
    auto [mol, dmat, bm] = getResults(smiles);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) > 1.0);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < 1.3);
    CHECK((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1));
    CHECK(bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10);
  }

  SECTION("chain of doubles") {
    const std::string smiles = "C=C=C=C";
    auto [mol, dmat, bm] = getResults(smiles);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
    CHECK(bm->getUpperBound(0, 3) > dmat->getVal(0, 3));
    // this is kinda goofy but this linear molecule doesn't satisfy the bounds
    // completely
    CHECK(fabs(bm->getLowerBound(0, 3) - dmat->getVal(0, 3)) < 0.2);
  }
  std::cerr << "-------------------------------------\n\n";
  SECTION("Ethene trans") {
    const std::string smiles = "C/C=C/C";
    auto [mol, dmat, bm] = getResults(smiles);
    // std::cerr << "\n-----\n" << MolToMolBlock(mol,false,cid) << std::endl;;
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
    // std::cerr << bm->getUpperBound(0,3) << " " << dmat->getVal(0,3) << " " <<
    // bm->getLowerBound(0,3) << std::endl;
    CHECK((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1));
    CHECK(bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.1);
  }
  SECTION("Ethene cis") {
    const std::string smiles = "C/C=C\\C";
    auto [mol, dmat, bm] = getResults(smiles);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
    CHECK((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1));
    CHECK((bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));
  }

  SECTION("Ethene No Stereo") {
    const std::string smiles = "CC=CC";
    auto [mol, dmat, bm] = getResults(smiles);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < 1.13);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) > 1.);
    CHECK((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1));
    CHECK((bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));
  }
  SECTION("Disulfur Dioxide") {
    const std::string smiles = "O=S-S=O";
    auto [mol, dmat, bm] = getResults(smiles);
    CHECK((bm->getUpperBound(0, 3) - bm->getLowerBound(0, 3)) < .13);
    CHECK((bm->getUpperBound(0, 3) - dmat->getVal(0, 3) > -0.1));
    CHECK((bm->getLowerBound(0, 3) - dmat->getVal(0, 3) < 0.10));
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("test3") {
  boost::logging::disable_logs("rdApp.warning");
  // check embedding based based on distances calculated from previous created
  // (good) coordinates
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/combi_coords.sdf";
  v2::FileParsers::SDMolSupplier sdsup(fname);

  bool gotCoords;
  while (!sdsup.atEnd()) {
    auto mol = sdsup.next();
    std::string mname;
    mol->getProp(common_properties::_Name, mname);
    RDGeom::PointPtrVect origCoords, newCoords;
    const std::size_t nat = mol->getNumAtoms();
    Conformer &conf = mol->getConformer(0);
    for (std::size_t i = 0; i < nat; i++) {
      origCoords.push_back(&conf.getAtomPos(i));
    }
    RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
    computeDistMat(origCoords, distMat);

    gotCoords = DistGeom::computeInitialCoords(distMat, origCoords);

    REQUIRE(gotCoords);
    RDNumeric::DoubleSymmMatrix distMatNew(nat, 0.0);
    computeDistMat(origCoords, distMatNew);

    for (std::size_t i = 1; i < nat; i++) {
      for (std::size_t j = 0; j < i; j++) {
        CHECK(RDKit::feq(distMat.getVal(i, j), distMatNew.getVal(i, j), 0.01));
      }
    }
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("test4") {
  boost::logging::disable_logs("rdApp.warning");
  auto m =
      "c1cc(C(F)(F)F)ccc1/C=N/NC(=O)c(n2)c[n]3cc(C(F)(F)F)cc(c23)Cl"_smiles;
  DGeomHelpers::EmbedMolecule(*m, 10, 1);  // etCoords(*m, iter);
  // std::string fname = "test.mol";
  // MolToMolFile(*m, fname);
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("test5") {
  // some real CDK2 molecules lets see how many fail
  std::string smifile =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv";
  v2::FileParsers::SmilesMolSupplierParams params{.delimiter = ","};
  v2::FileParsers::SmilesMolSupplier smiSup(smifile, params);

  for (auto mol : smiSup) {
    MolOps::addHs(*mol);
    const auto cid =
        DGeomHelpers::EmbedMolecule(*mol, 10, 1);  // getCoords(*mol,
    CHECK(cid > -1);
  }
}

TEST_CASE("test6") {
  const auto m = "CC"_smiles;
  REQUIRE(m);
  DistGeom::BoundsMatPtr bm;
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  REQUIRE(bm);
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  CHECK(feq(bm->getLowerBound(0, 1), 0.0));
  CHECK(feq(bm->getLowerBound(1, 0), 0.0));
  CHECK(feq(bm->getUpperBound(0, 1), 1000.0));
  CHECK(feq(bm->getUpperBound(1, 0), 1000.0));

  DGeomHelpers::setTopolBounds(*m, bm);
  CHECK(bm->getLowerBound(0, 1) > 0.0);
  CHECK(bm->getUpperBound(0, 1) < 1000.0);
  CHECK(bm->getLowerBound(0, 1) < bm->getUpperBound(0, 1));
}

TEST_CASE("testIssue215") {
  const std::string smiles = "C=C1C2CC1C2";
  auto bm = _getBoundsMatrix(smiles);
  // this was the specific problem:
  CHECK(bm->getUpperBound(0, 4) < 100.0);

  const bool ok = DistGeom::triangleSmoothBounds(bm);
  CHECK(ok);
}

TEST_CASE("test15Dists") {
  SECTION("1") {
    std::string smiles = "c1ccccc1C";
    const auto mmat = _getBoundsMatrix(smiles);
    CHECK(RDKit::feq(mmat->getUpperBound(2, 6), 4.32, 0.01));
    CHECK(RDKit::feq(mmat->getLowerBound(2, 6), 4.16, 0.01));
  }
  SECTION("2") {
    const std::string smiles = "CC1=C(C(C)=CC=C2)C2=CC=C1";
    const auto mmat = _getBoundsMatrix(smiles);
    CHECK(RDKit::feq(mmat->getLowerBound(0, 4), 2.31, 0.01));
    CHECK(RDKit::feq(mmat->getUpperBound(0, 4), 2.47, 0.01));
    CHECK(RDKit::feq(mmat->getLowerBound(4, 11), 4.11, 0.01));
    CHECK(RDKit::feq(mmat->getUpperBound(4, 11), 4.27, 0.01));
  }
  SECTION("3") {
    const std::string smiles = "C/C=C/C=C/C";
    const auto mmat = _getBoundsMatrix(smiles);
    CHECK(RDKit::feq(mmat->getLowerBound(0, 4), 4.1874));
    CHECK(RDKit::feq(mmat->getUpperBound(0, 4), 4.924));
    CHECK(RDKit::feq(mmat->getLowerBound(1, 5), 4.1874));
    CHECK(RDKit::feq(mmat->getUpperBound(1, 5), 4.924));
  }
}

TEST_CASE("testMultipleConfs") {
  boost::logging::disable_logs("rdApp.warning");
  const auto m = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC"_smiles;
  SECTION("DG") {
    INT_VECT cids =
        DGeomHelpers::EmbedMultipleConfs(*m, 10, 30, 100, true, false, -1);
    // SDWriter writer("junk.sdf");
    for (const auto ci : cids) {
      std::unique_ptr<ForceFields::ForceField> ff{
          UFF::constructForceField(*m, 10, ci)};
      ff->initialize();
      double energy = ff->calcEnergy();
      // BOOST_LOG(rdInfoLog) << energy << std::endl;
      CHECK(energy > 100.0);
      CHECK(energy < 300.0);
    }
  }
  SECTION("ExpTors") {
    INT_VECT cids = DGeomHelpers::EmbedMultipleConfs(
        *m, 10, 30, 100, true, false, -1, true, 1, -1.0, nullptr, 1e-3, false,
        true, false, false, false, 5.0, false, 1, false, false);

    // SDWriter writer("junk.sdf");
    for (const auto ci : cids) {
      std::unique_ptr<ForceFields::ForceField> ff{
          UFF::constructForceField(*m, 10, ci)};
      ff->initialize();
      double energy = ff->calcEnergy();
      // BOOST_LOG(rdInfoLog) << energy << std::endl;
      CHECK(energy > 100.0);
      CHECK(energy < 300.0);
    }
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("test Set Topol Bounds") {
  SECTION("1") {
    const std::string smiles = "CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12";
    _getBoundsMatrix(smiles);
  }
  SECTION("2") {
    const std::string smiles = "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2";
    _getBoundsMatrix(smiles);
  }
  SECTION("Issue244 1") {
    const std::string smi = "NC1C(=O)N2C1SCC(Cl)=C2C(O)=O";
    _getBoundsMatrix(smi);
  }
  SECTION("Issue244 2") {
    const std::string smi = "NC1C(=O)N2C1SCC(Cl)=C2C(O)=O";
    _getBoundsMatrix(smi);
  }
  SECTION("Issue244 3") {
    const std::string smi = "CC1(C)SC2N(C1C(O)=O)C(=O)C2N";
    _getBoundsMatrix(smi);
  }
}

TEST_CASE("testTriangleSmoothing") {
  auto runTest = [](const std::string &smiles) {
    auto bm = _getBoundsMatrix(smiles);
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);
  };
  SECTION("testIssue227 1") {
    const std::string smiles =
        "CCOP1(OCC)=CC(c2ccccc2)=C(c2ccc([N+]([O-])=O)cc2)N=C1c1ccccc1";
    runTest(smiles);
  }

  SECTION("testIssue227 1") {
    const std::string smiles =
        "OC(=O)c1cc2cc(c1)-c1c(O)c(ccc1)-c1cc(C(O)=O)cc(c1)OCCOCCO2";
    runTest(smiles);
  }
  SECTION("testissue236 1") {
    const std::string smiles =
        "Cn1c2n(-c3ccccc3)c(=O)c3c(nc4ccc([N+]([O-])=O)cc4c3)c2c(=O)n(C)c1=O";
    runTest(smiles);
  }

  SECTION("testissue236 2") {
    const std::string smiles = "Cc1cccc2c1c(C3=CCC3)c(C)cc2";
    runTest(smiles);
  }
  SECTION("testIssue276") {
    const std::string smiles = "CP1(C)=CC=CN=C1C";
    runTest(smiles);
  }
}

TEST_CASE("testIssue251") {
  const auto m = "COC=O"_smiles;
  auto bm = _getBoundsMatrix(m);
  CHECK(RDKit::feq(bm->getLowerBound(0, 3), 2.67, 0.01));
  CHECK(RDKit::feq(bm->getUpperBound(0, 3), 2.79, 0.01));
}

TEST_CASE("testIssue284") {
  std::string smi = "CNC(=O)C";
  auto bm = _getBoundsMatrix(smi);
  bool ok = DistGeom::triangleSmoothBounds(bm);
  CHECK(ok);

  // amide bonds are cis-oid:
  CHECK(bm->getLowerBound(0, 3) < 3.0);
  CHECK(bm->getUpperBound(0, 3) < 3.0);

  const std::string smi2 = "CN(C)C(=O)C";
  auto bm2 = _getBoundsMatrix(smi2);
  ok = DistGeom::triangleSmoothBounds(bm2);
  CHECK(ok);

  // we've got no information to tell us cis-oid vs trans-oid here, so
  // the windows are huge:
  CHECK(bm2->getLowerBound(0, 4) < 3.0);
  CHECK(bm2->getUpperBound(0, 4) > 3.5);
  CHECK(bm2->getLowerBound(2, 4) < 3.0);
  CHECK(bm2->getUpperBound(2, 4) > 3.5);
  CHECK(bm->getLowerBound(0, 3) < bm2->getLowerBound(0, 4));
  CHECK(bm->getUpperBound(0, 3) < bm2->getUpperBound(0, 4));
  CHECK(bm->getLowerBound(0, 3) < bm2->getLowerBound(2, 4));
  CHECK(bm->getUpperBound(0, 3) < bm2->getUpperBound(2, 4));
}

TEST_CASE("testIssue285") {
  auto m = "CNC(=O)C"_smiles;
  REQUIRE(m);
  MolOps::addHs(*m);
  REQUIRE(m);
  auto bm = _getBoundsMatrix(m);
  bool ok = DistGeom::triangleSmoothBounds(bm);
  CHECK(ok);

  std::size_t tgtNumber = 10;
  INT_VECT cids = DGeomHelpers::EmbedMultipleConfs(*m, tgtNumber);
  REQUIRE(cids.size() == tgtNumber);

  std::vector<std::string> molBlocks;
  for (auto cid : cids) {
    molBlocks.push_back(MolToMolBlock(*m, true, cid));
  }
  for (auto mbI = molBlocks.begin(); mbI != molBlocks.end(); ++mbI) {
    for (auto mbJ = mbI + 1; mbJ != molBlocks.end(); ++mbJ) {
      CHECK((*mbI) != (*mbJ));
    }
    // std::cerr << (*mbI) << "\n$$$$\n";
  }
}

TEST_CASE("testIssue355") {
  SECTION("1") {
    const std::string smi = "CNC(=O)C";
    const auto bm = _getBoundsMatrix(smi);

    CHECK(bm->getLowerBound(0, 3) < 3.0);
    CHECK(bm->getUpperBound(0, 3) < 3.0);
    CHECK(bm->getUpperBound(0, 4) > 3.2);
    CHECK(bm->getLowerBound(0, 4) > 3.2);

    const bool ok = DistGeom::triangleSmoothBounds(bm);
    CHECK(ok);
  }

  SECTION("2") {
    const std::string smi = "CNC(=O)NC";
    const auto bm = _getBoundsMatrix(smi);
    CHECK(bm->getLowerBound(0, 3) < 3.0);
    CHECK(bm->getUpperBound(0, 3) < 3.0);
    CHECK(bm->getUpperBound(0, 4) > 3.2);
    CHECK(bm->getLowerBound(0, 4) > 3.2);
    CHECK(bm->getLowerBound(5, 3) < 3.0);
    CHECK(bm->getUpperBound(5, 3) < 3.0);
    CHECK(bm->getUpperBound(5, 1) > 3.2);
    CHECK(bm->getLowerBound(5, 1) > 3.2);
  }

  SECTION("3") {
    const std::string smi = "CNC(=O)Nc1ccccc1";
    const auto bm = _getBoundsMatrix(smi);
    CHECK(bm->getLowerBound(0, 3) < 3.0);
    CHECK(bm->getUpperBound(0, 3) < 3.0);
    CHECK(bm->getUpperBound(0, 4) > 3.2);
    CHECK(bm->getLowerBound(0, 4) > 3.2);
    CHECK(bm->getLowerBound(5, 3) < 3.0);
    CHECK(bm->getUpperBound(5, 3) < 3.0);
    CHECK(bm->getUpperBound(5, 1) > 3.2);
    CHECK(bm->getLowerBound(5, 1) > 3.2);
  }
}

TEST_CASE("testRandomCoords") {
  std::vector<std::string> smiString = {
      "CC1=C(C(C)=CC=C2)C2=CC=C1",
      "c1ccccc1C",
      "C/C=C/CC",
      "C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4",
      "C1CCCCS1(=O)(=O)",
      "c1ccccc1",
      "C1CCCC1",
      "C1CCCCC1",
      "C1CC1(C)C",
      "C12(C)CC1CC2"};
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/initCoords.random.sdf";
  v2::FileParsers::MolFileParserParams params{.removeHs = false};
  v2::FileParsers::SDMolSupplier sdsup(fname, params);
  // SDWriter writer("foo.sdf");
  // SDWriter writer(fname);
  for (const auto &smi : smiString) {
    // std::cerr << "SMI: " << smi << std::endl;
    auto m = v2::SmilesParse::MolFromSmiles(smi);
    MolOps::addHs(*m);
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, true, true, 2, true, 1,
                                          nullptr, 1e-2);
    CHECK_INVARIANT(cid >= 0, "");
    // writer.write(*m);
    // writer.flush();
    auto m2 = sdsup.next();
    // ROMol *m2 = NULL;
    if (m2) {
      CHECK(m->getNumAtoms() == m2->getNumAtoms());
      unsigned int nat = m->getNumAtoms();

      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
      for (unsigned int i = 0; i < nat; i++) {
        RDGeom::Point3D pt1i = conf1.getAtomPos(i);
        RDGeom::Point3D pt2i = conf2.getAtomPos(i);
        for (unsigned int j = i + 1; j < nat; j++) {
          RDGeom::Point3D pt1j = conf1.getAtomPos(j);
          RDGeom::Point3D pt2j = conf2.getAtomPos(j);
          double d1 = (pt1j - pt1i).length();
          double d2 = (pt2j - pt2i).length();
          if (m->getBondBetweenAtoms(i, j)) {
            CHECK(fabs(d1 - d2) / d1 < 0.05);
          } else {
            CHECK(fabs(d1 - d2) / d1 < 0.1);
          }
        }
      }
    }
  }
}

TEST_CASE("testIssue1989539") {
  {
    auto m = "c1ccccc1.Cl"_smiles;
    const int cid = DGeomHelpers::EmbedMolecule(*m);
    CHECK(cid >= 0);
    const std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    CHECK(cids.size() == 10);
    CHECK(std::find(cids.begin(), cids.end(), -1) == cids.end());
  }
  {
    auto m = "[Cl-].c1ccccc1C[NH3+]"_smiles;
    const int cid = DGeomHelpers::EmbedMolecule(*m);
    CHECK(cid >= 0);
    const std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    CHECK(cids.size() == 10);
    CHECK(std::find(cids.begin(), cids.end(), -1) == cids.end());
  }
}

TEST_CASE("testConstrainedEmbedding") {
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf";
  v2::FileParsers::SDMolSupplier sdsup(fname);
  auto ref = sdsup.next();
  SECTION("1") {
    auto test = std::make_unique<RWMol>(*ref);
    MolOps::addHs(*test);
    std::map<int, RDGeom::Point3D> coords;
    coords[0] = ref->getConformer().getAtomPos(0);
    coords[1] = ref->getConformer().getAtomPos(1);
    coords[2] = ref->getConformer().getAtomPos(2);
    coords[3] = ref->getConformer().getAtomPos(3);
    coords[4] = ref->getConformer().getAtomPos(4);

    const int cid = DGeomHelpers::EmbedMolecule(*test, 30, 22, true, false, 2.,
                                                true, 1, &coords);
    REQUIRE(cid > -1);

    MatchVectType alignMap;
    alignMap.push_back(std::make_pair(0, 0));
    alignMap.push_back(std::make_pair(1, 1));
    alignMap.push_back(std::make_pair(2, 2));
    alignMap.push_back(std::make_pair(3, 3));
    alignMap.push_back(std::make_pair(4, 4));
    const double ssd = MolAlign::alignMol(*test, *ref, -1, -1, &alignMap);
    BOOST_LOG(rdInfoLog) << "ssd: " << ssd << std::endl;
    CHECK(ssd < 0.1);
  }
  SECTION("2") {
    auto test = sdsup.next();
    MolOps::addHs(*test);
    std::map<int, RDGeom::Point3D> coords;
    coords[4] = ref->getConformer().getAtomPos(0);
    coords[5] = ref->getConformer().getAtomPos(1);
    coords[6] = ref->getConformer().getAtomPos(2);
    coords[7] = ref->getConformer().getAtomPos(3);
    coords[8] = ref->getConformer().getAtomPos(4);
    const int cid = DGeomHelpers::EmbedMolecule(*test, 30, 22, true, false, 2.,
                                                true, 1, &coords);
    REQUIRE(cid > -1);

    MatchVectType alignMap;
    alignMap.push_back(std::make_pair(4, 0));
    alignMap.push_back(std::make_pair(5, 1));
    alignMap.push_back(std::make_pair(6, 2));
    alignMap.push_back(std::make_pair(7, 3));
    alignMap.push_back(std::make_pair(8, 4));
    const double ssd = MolAlign::alignMol(*test, *ref, -1, -1, &alignMap);
    BOOST_LOG(rdInfoLog) << "ssd: " << ssd << std::endl;
    CHECK(ssd < 0.1);
  }
}

TEST_CASE("Check Mols embedddable") {
  auto check = [](RWMol &m) {
    int cid = DGeomHelpers::EmbedMolecule(m);
    CHECK(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(m, 10);
    CHECK(cids.size() == 10);
    CHECK(std::find(cids.begin(), cids.end(), -1) == cids.end());
  };
  boost::logging::disable_logs("rdApp.warning");
  SECTION("testIssue2091864 1") {
    auto m = "C1C2CC12"_smiles;
    int cid = DGeomHelpers::EmbedMolecule(*m);
    CHECK(cid >= 0);
  }
  SECTION("testIssue2091864 2") {
    auto m = "C1CC2C3C1C23"_smiles;
    check(*m);
  }
  SECTION("testIssue2091864 3") {
    auto m = "c1ccc2c(c1)C1C3C2C13"_smiles;
    check(*m);
  }
  SECTION("testIssue3019283") {
    auto m = "C1=C2C1C1CC21"_smiles;
    MolOps::addHs(*m);
    check(*m);
  }
  SECTION("testIssue2091974 1") {
    auto m = "CCOC(OCC)(OCC)OCC"_smiles;
    MolOps::addHs(*m);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    CHECK(cid >= 0);
  }
  SECTION("testIssue2091974 2") {
    auto m = "O=N(=O)OCC(CON(=O)=O)(CON(=O)=O)CON(=O)=O"_smiles;
    MolOps::addHs(*m);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    CHECK(cid >= 0);
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("testIssue2835784") {
  boost::logging::disable_logs("rdApp.warning");
  auto runTest = [](const std::string &smiles, const bool addHs) {
    auto m = v2::SmilesParse::MolFromSmiles(smiles);
    if (addHs) {
      MolOps::addHs(*m);
    }
    int cid = DGeomHelpers::EmbedMolecule(*m);
    CHECK(cid >= 0);
    std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(*m, 10);
    CHECK(cids.size() == 10);
    CHECK(std::find(cids.begin(), cids.end(), -1) == cids.end());
  };

  SECTION("1") {
    std::string smi = "C1C=C1";
    auto addHs = GENERATE(false, true);
    runTest(smi, addHs);
  }
  SECTION("2") {
    std::string smi = "C12=CCC1C2";
    auto addHs = GENERATE(false, true);
    runTest(smi, addHs);
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("testIssue3238580") {
  SECTION("Basic") {
    std::string smi = "C1CCC2=CC12";
    auto bm = _getBoundsMatrix(smi);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    CHECK(ok);
  }
  SECTION("Issue3238580 1") {
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.1.mol";
    auto m = v2::FileParsers::MolFromMolFile(molfile);
    auto bm = _getBoundsMatrix(m);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    CHECK(ok);
  }
  SECTION("Issue3238580 2") {
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.2.mol";
    auto m = v2::FileParsers::MolFromMolFile(molfile);
    auto bm = _getBoundsMatrix(m);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    CHECK(ok);
  }
  SECTION("Issue3238580 3") {
    std::string rdbase = getenv("RDBASE");
    std::string molfile =
        rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.3.mol";
    auto m = v2::FileParsers::MolFromMolFile(molfile);
    auto bm = _getBoundsMatrix(m);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    CHECK(ok);
  }
}

TEST_CASE("testIssue3483968") {
  boost::logging::disable_logs("rdApp.warning");
  std::string molfile =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3483968.mol";
  auto m = v2::FileParsers::MolFromMolFile(molfile);
  REQUIRE(m);
  int cid = DGeomHelpers::EmbedMolecule(*m, 0, -1, true, false, 2.0, true, 1,
                                        nullptr, 1e-3, true);
  CHECK(cid >= 0);
  std::vector<int> cids = DGeomHelpers::EmbedMultipleConfs(
      *m, 10, 30, 1, true, false, 2.0, true, 1, -1.0, nullptr, 1e-3, true);
  CHECK(cids.size() == 10);
  CHECK(std::find(cids.begin(), cids.end(), -1) == cids.end());
  boost::logging::enable_logs("rdApp.warning");
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock(const std::vector<std::shared_ptr<RWMol>> &mols,
              const std::vector<double> &energies, unsigned int count,
              unsigned int idx) {
  for (std::size_t j = 0; j < 100; j++) {
    for (std::size_t i = 0; i < mols.size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      auto &mol = mols[i];
      std::vector<int> cids =
          DGeomHelpers::EmbedMultipleConfs(*mol, 10, 30, 0xFEED);
      REQUIRE(cids.size() == 10);
      std::unique_ptr<ForceFields::ForceField> field(
          UFF::constructForceField(*mol, 100, cids[0]));
      REQUIRE(field);
      field->initialize();
      double eng = field->calcEnergy();
      if (!feq(eng, energies[i])) {
        std::cerr << i << " iter " << j << " " << energies[i] << " != " << eng
                  << std::endl;
      }
      CHECK(feq(eng, energies[i]));
    }
  }
}
}  // namespace

#include <thread>
#include <future>
TEST_CASE("testMultiThread") {
  std::cerr << "building molecules" << std::endl;
  // std::string smi="C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4";
  std::string smi = "c1ccc2c(c1)C1C3C2C13";
  std::vector<std::shared_ptr<RWMol>> mols;
  mols.reserve(100);
  auto mol = v2::SmilesParse::MolFromSmiles(smi);
  MolOps::addHs(*mol);
  for (unsigned int i = 0; i < 100; ++i) {
    mols.push_back(std::make_shared<RWMol>(*mol));
  }

  std::cerr << "generating reference data" << std::endl;
  std::vector<double> energies(mols.size(), 0.0);
  for (unsigned int i = 0; i < mols.size(); ++i) {
    auto &mol = mols[i];
    std::vector<int> cids =
        DGeomHelpers::EmbedMultipleConfs(*mol, 10, 30, 0xFEED);
    REQUIRE(cids.size() == 10);
    std::unique_ptr<ForceFields::ForceField> field(
        UFF::constructForceField(*mol, 100, cids[0]));
    REQUIRE(field);
    field->initialize();
    double eng = field->calcEnergy();
    CHECK(eng != 0.0);
    energies[i] = eng;
  }

  std::cerr << "validating reference data" << std::endl;
  for (unsigned int i = 0; i < mols.size(); ++i) {
    auto &mol = mols[i];
    std::vector<int> cids =
        DGeomHelpers::EmbedMultipleConfs(*mol, 10, 30, 0xFEED);
    REQUIRE(cids.size() == 10);
    std::unique_ptr<ForceFields::ForceField> field(
        UFF::constructForceField(*mol, 100, cids[0]));
    REQUIRE(field);
    field->initialize();
    double eng = field->calcEnergy();
    CHECK(feq(eng, energies[i]));
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

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread() {}
#endif

TEST_CASE("testGithub55") {
  boost::logging::disable_logs("rdApp.warning");
  SECTION("basic") {
    auto core = "c1cnco1"_smiles;
    REQUIRE(core);

    int cid = DGeomHelpers::EmbedMolecule(*core);
    CHECK(cid >= 0);

    auto mol = "o1cncc1C"_smiles;
    REQUIRE(mol);

    std::map<int, RDGeom::Point3D> coords;
    coords[0] = core->getConformer().getAtomPos(4);
    coords[1] = core->getConformer().getAtomPos(3);
    coords[2] = core->getConformer().getAtomPos(2);
    coords[3] = core->getConformer().getAtomPos(1);
    coords[4] = core->getConformer().getAtomPos(0);
    cid = DGeomHelpers::EmbedMolecule(*mol, 50, 22, true, false, 2., true, 1,
                                      &coords);
    CHECK(cid > -1);
  }
  SECTION("With Sulfur") {
    auto core = "c1cncs1"_smiles;
    REQUIRE(core);

    int cid = DGeomHelpers::EmbedMolecule(*core);
    CHECK(cid >= 0);

    auto mol = "s1cncc1C"_smiles;
    REQUIRE(mol);

    std::map<int, RDGeom::Point3D> coords;
    coords[0] = core->getConformer().getAtomPos(4);
    coords[1] = core->getConformer().getAtomPos(3);
    coords[2] = core->getConformer().getAtomPos(2);
    coords[3] = core->getConformer().getAtomPos(1);
    coords[4] = core->getConformer().getAtomPos(0);
    cid = DGeomHelpers::EmbedMolecule(*mol, 50, 22, true, false, 2., true, 1,
                                      &coords);
    CHECK(cid > -1);
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("testGithub256") {
  auto mol = std::make_unique<RWMol>();
  REQUIRE(mol);
  CHECK_THROWS_AS(DGeomHelpers::EmbedMolecule(*mol), ValueErrorException);
}

#ifdef RDK_TEST_MULTITHREADED
TEST_CASE("testMultiThreadMultiConf") {
  boost::char_separator<char> sep("|");
  auto bldString = std::string(RDKit::rdkitBuild);
  tokenizer tokens(bldString, sep);
  std::vector<std::string> tokenVect(tokens.begin(), tokens.end());
  const double ENERGY_TOLERANCE = ((tokenVect[2] != "MINGW") ? 1.0e-6 : 1.0);
  const double MSD_TOLERANCE = ((tokenVect[2] != "MINGW") ? 1.0e-6 : 1.0e-5);
  auto m = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC"_smiles;
  REQUIRE(m);
  MolOps::addHs(*m);
  RWMol m2(*m);
  INT_VECT cids;
  DGeomHelpers::EmbedMultipleConfs(*m, cids, 200, 1, 30, 100, true, false, -1);
  DGeomHelpers::EmbedMultipleConfs(m2, cids, 200, 0, 30, 100, true, false, -1);
  for (auto ci : cids) {
    std::unique_ptr<ForceFields::ForceField> ff(
        UFF::constructForceField(*m, 100, ci));
    ff->initialize();
    double e1 = ff->calcEnergy();
    const RDGeom::PointPtrVect &pVect = ff->positions();
    CHECK(e1 > 100.0);
    CHECK(e1 < 300.0);
    std::unique_ptr<ForceFields::ForceField> ff2(
        UFF::constructForceField(m2, 100, ci));
    ff2->initialize();
    double e2 = ff2->calcEnergy();
    const RDGeom::PointPtrVect &p2Vect = ff2->positions();
    CHECK(RDKit::feq(e1, e2, ENERGY_TOLERANCE));
    CHECK(pVect.size() == p2Vect.size());
    double msd = 0.0;
    for (unsigned int i = 0; i < pVect.size(); ++i) {
      const auto *p = dynamic_cast<const RDGeom::Point3D *>(pVect[i]);
      const auto *p2 = dynamic_cast<const RDGeom::Point3D *>(p2Vect[i]);
      REQUIRE(p);
      REQUIRE(p2);
      msd += (*p - *p2).lengthSq();
    }
    msd /= static_cast<double>(pVect.size());
    CHECK(msd < MSD_TOLERANCE);
  }
}
#endif

TEST_CASE("testGithub563") {
  {
    std::string smi = "[H][C@]1(C[NH3+])CC[C@]([H])(CNC)CC1";
    auto m = v2::SmilesParse::MolFromSmiles(smi);
    std::string csmi = MolToSmiles(*m, true);
    std::cerr << csmi << std::endl;
    for (unsigned int i = 1; i < 100; ++i) {
      RWMol m2 = *m;
      MolOps::addHs(m2);
      DGeomHelpers::EmbedMolecule(m2, 50, i);
      MolOps::assignChiralTypesFrom3D(m2);
      MolOps::removeHs(m2);
      std::string smi = MolToSmiles(m2, true);
      CHECK(smi == csmi);
    }
  }
}

TEST_CASE("testGithub568") {
  {
    // sample molecules (either from ChEMBL or derived from ChEMBL) that were
    // problematic
    std::vector<std::string> smis = {
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
        "Cc1nc(COc2ccc3OC[C@H](Cc4cccnc4)[C@H](O)c3c2)ccc1[N+](=O)[O-]"};
    for (const auto &smi : smis) {
      auto m = v2::SmilesParse::MolFromSmiles(smi);
      std::string csmi = MolToSmiles(*m, true);
      std::cerr << csmi << std::endl;
      // increase the limit here to make this a real torture test
      for (unsigned int i = 1; i < 20; ++i) {
        RWMol m2 = *m;
        MolOps::addHs(m2);
        int cid = DGeomHelpers::EmbedMolecule(m2, 50, i);
        CHECK(cid >= 0);
        MolOps::assignChiralTypesFrom3D(m2);

        // m2.setProp("_Name",smis[idx]);
        // std::cerr<<MolToMolBlock(m2)<<std::endl;
        // TEST_ASSERT(0);
        MolOps::removeHs(m2);
        std::string smi2 = MolToSmiles(m2, true);
        if (smi2 != csmi) {
          std::cerr << "-------------" << std::endl;
          std::cerr << smi << " " << i << std::endl;
          std::cerr << smi2 << "\n" << csmi << std::endl;
          m2.setProp("_Name", smi);
          std::cerr << MolToMolBlock(m2) << std::endl;
        }
        CHECK(smi2 == csmi);
      }
    }
  }
}

TEST_CASE("testGithub696") {
  {
    auto m = "COc1ccc2CCCCCCCCCCc3ccc(OC)c(c3)-c1c2"_smiles;
    // m->debugMol(std::cerr);
    DistGeom::BoundsMatPtr bm;
    bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm, false, false);

    CHECK(bm->getUpperBound(2, 19) > bm->getLowerBound(2, 19));
    CHECK(bm->getLowerBound(2, 19) > 2.0);
    CHECK(bm->getUpperBound(2, 19) > 2.5);

    bool ok = DistGeom::triangleSmoothBounds(bm);
    CHECK(ok);
  }
}

TEST_CASE("testGithub697") {
  // a group of chembl molecules (and things derived from them), all of which
  // contain a c1cscn1 heterocycle
  std::vector<std::string> smis = {
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
      "C(COC(=O)C(N)CCCNC(=N)N)NC(=O)C(C(C)C)NC1=O"};

  for (const auto &smi : smis) {
    ROMol *m = SmilesToMol(smi);
    REQUIRE(m);
    DistGeom::BoundsMatPtr bm;
    bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm, false, false);
    const bool ok = DistGeom::triangleSmoothBounds(bm);
    if (!ok) {
      m->debugMol(std::cerr);
      std::cerr << " FAILED: " << smi << std::endl;
    }
    CHECK(ok);
  }
}

TEST_CASE("testGithub971") {
  // sample molecule found by Sereina
  auto m = "C/C(=C\\c1ccccc1)CN1C2CC[NH2+]CC1CC2"_smiles;
  REQUIRE(m);
  MolOps::addHs(*m);
  int cid = DGeomHelpers::EmbedMolecule(*m, 0, 0xf00d);
  CHECK(cid >= 0);
  MolOps::removeHs(*m);
  std::string expectedMb = R"CTAB(
     RDKit          3D

 19 21  0  0  0  0  0  0  0  0999 V2000
    1.7021   -0.5378    1.1920 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1689   -0.7539   -0.1895 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0808   -1.0678   -1.1107 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4918   -0.7501   -0.9384 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3219   -1.4658   -0.0847 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2429   -0.7644    0.6840 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7365    0.4471    0.2159 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1309    0.9731   -0.9166 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0771    0.3483   -1.5564 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2344   -1.1991   -0.2964 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1636   -0.1335   -0.5561 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3130    0.7928    0.5372 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1519    1.9502    0.1036 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6156    1.8065    0.3639 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3065    0.9392   -0.5468 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4579    0.1053   -1.3418 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4944   -0.7376   -0.5875 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8236   -0.9671    0.8498 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9211   -0.0480    1.6167 C   0  0  0  0  0  0  0  0  0  0  0  0
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
  auto expected = v2::FileParsers::MolFromMolBlock(expectedMb);
  unsigned int nat = expected->getNumAtoms();
  CHECK(nat == m->getNumAtoms());

  compareConfs(m.get(), expected.get(), 0, 0);
}

TEST_CASE("testEmbedParameters") {
  auto runTest = [](const std::string &smiles, const std::string &fname,
                    DGeomHelpers::EmbedParameters &params) {
    auto ps = v2::FileParsers::MolFileParserParams{.removeHs = false};
    auto ref = v2::FileParsers::MolFromMolFile(fname, ps);
    REQUIRE(ref);
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(ref->getNumAtoms() == mol->getNumAtoms());
    params.randomSeed = 42;
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
    // std::cerr << MolToMolBlock(*ref) << std::endl;
    // std::cerr << MolToMolBlock(*mol) << std::endl;
    // std::cerr << fname << std::endl;
  };
  SECTION("default params") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 42;
    runTest(smiles, fname, params);
  }
  SECTION("default etdg") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params;
    params.useExpTorsionAnglePrefs = true;
    runTest(smiles, fname, params);
  }
  SECTION("ETKDGv1") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params;
    params.useExpTorsionAnglePrefs = true;
    params.useBasicKnowledge = true;
    runTest(smiles, fname, params);
  }
  SECTION("ETKDGv2") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol";
    std::string smiles = "n1cccc(C)c1ON";
    DGeomHelpers::EmbedParameters params;
    params.useExpTorsionAnglePrefs = true;
    params.useBasicKnowledge = true;
    params.ETversion = 2;
    runTest(smiles, fname, params);
  }
  SECTION("KDG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params;
    params.useBasicKnowledge = true;
    runTest(smiles, fname, params);
  }
  //------------
  // using the pre-defined parameter sets
  SECTION("predefined ETDG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETDG);
    runTest(smiles, fname, params);
  }
  SECTION("predefined ETKDG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    runTest(smiles, fname, params);
  }
  SECTION("predefined KDG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol";
    std::string smiles = "OCCC";
    DGeomHelpers::EmbedParameters params(DGeomHelpers::KDG);
    runTest(smiles, fname, params);
  }
  SECTION("predefined srETKDGv3") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.smallring.etkdgv3.mol";
    std::string smiles = "C1CCCCC1";
    DGeomHelpers::EmbedParameters params(DGeomHelpers::srETKDGv3);
    runTest(smiles, fname, params);
  }
  SECTION("predefined ETKDG - macrocycle") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdg.mol";
    std::string smiles = "O=C1NCCCCCCCCC1";
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
    runTest(smiles, fname, params);
  }
  SECTION("predefined ETKDGv3 - macrocycle") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdgv3.mol";
    std::string smiles = "C1NCCCCCCCCC1";
    DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDGv3);
    runTest(smiles, fname, params);
  }
}

TEST_CASE("testGithub1227") {
  auto m = "CC"_smiles;
  REQUIRE(m);
  MolOps::addHs(*m);
  CHECK(m->getNumAtoms() == 8);
  INT_VECT cids;
  DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
  params.randomSeed = 0xf00d;
  cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, params);
  CHECK(cids.size() == 10);

  params.pruneRmsThresh = 0.5;
  cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, params);
  CHECK(cids.size() == 1);

  params.onlyHeavyAtomsForRMS = false;  // the old default behavior
  params.useSymmetryForPruning = false;
  cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, params);
  CHECK(cids.size() == 6);
}

TEST_CASE("testGithub1240") {
  auto runTest = [](const std::string &smiles,
                    DGeomHelpers::EmbedParameters &params,
                    const int randomSeed) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(mol);
    params.randomSeed = randomSeed;
    int cid = DGeomHelpers::EmbedMolecule(*mol, params);
    CHECK(cid >= 0);
  };
  SECTION("1") {
    const std::string smiles = "C1CCCCCCC1";
    DGeomHelpers::EmbedParameters params;
    params.maxIterations = 1;
    runTest(smiles, params, 42);
  }
  SECTION("2") {
    const std::string smiles = "C1C3CC2CC(CC1C2)C3";
    DGeomHelpers::EmbedParameters params;
    params.maxIterations = 1;
    boost::logging::disable_logs("rdApp.warning");
    runTest(smiles, params, 42);
    boost::logging::enable_logs("rdApp.warning");
  }
  SECTION("3") {
    const std::string smiles =
        GENERATE("c1ccccccccc1", "c1ccccccc1", "c1ccccccccccc1");
    DGeomHelpers::EmbedParameters params;
    params.maxIterations = 1;
    runTest(smiles, params, 0xf00d);
  }
  SECTION("CHEMBL307150") {
    const std::string smiles =
        "Cc1cc2ccn(C)c2c3c4C(=O)NC(=O)c4c5c6ccccc6[nH]c5c13";
    DGeomHelpers::EmbedParameters params;
    runTest(smiles, params, 0xf00d);
    params = DGeomHelpers::ETKDG;
    runTest(smiles, params, 0xf00d);
  }

  SECTION("CHEMBL43398") {
    const std::string smiles =
        "C[C@@H]1[C@@H]2Cc3ccc(O)cc3[C@]1(C)CCN2CCN4CCCC4";
    DGeomHelpers::EmbedParameters params;
    runTest(smiles, params, 0xf00d);
    params = DGeomHelpers::ETKDG;
    runTest(smiles, params, 0xf00d);
  }
  SECTION("CHEMBL290986") {
    const std::string smiles =
        "COc1c(O)ccc2O[C@@H]([C@@H]3CCCC(=C3)C)c4c(ccc5NC(C)(C)C=C(C)c45)c12";
    DGeomHelpers::EmbedParameters params;
    runTest(smiles, params, 0xf00d);
    params = DGeomHelpers::ETKDG;
    runTest(smiles, params, 0xf00d);
  }
}

TEST_CASE("testGithubPullRequest1635") {
  auto m = "C1(F)(F)CCC(CC1)COCC(C23CC4CC(C2)CC(C4)C3)N"_smiles;
  REQUIRE(m);
  MolOps::addHs(*m);
  const int expected_num_atoms = 54;
  CHECK(m->getNumAtoms() == expected_num_atoms);
  RWMol firstMol(*m);
  RWMol secondMol(*m);

  DGeomHelpers::EmbedParameters params(DGeomHelpers::ETKDG);
  params.randomSeed = MAX_INT;  // the largest possible random seed

  INT_VECT firstCids = DGeomHelpers::EmbedMultipleConfs(firstMol, 10, params);
  INT_VECT secondCids = DGeomHelpers::EmbedMultipleConfs(secondMol, 10, params);
  CHECK(firstCids.size() == 10);
  CHECK(secondCids.size() == 10);

  for (size_t i = 0; i < 10; i++) {
    CHECK(firstCids[i] == secondCids[i]);

    int confIdx = firstCids[i];
    const Conformer &firstConf = firstMol.getConformer(confIdx);
    const Conformer &secondConf = secondMol.getConformer(confIdx);

    for (int atomIdx = 0; atomIdx < expected_num_atoms; ++atomIdx) {
      const RDGeom::Point3D &firstPoint = firstConf.getAtomPos(atomIdx);
      const RDGeom::Point3D &secondPoint = secondConf.getAtomPos(atomIdx);
      CHECK(firstPoint.x == secondPoint.x);
      CHECK(firstPoint.y == secondPoint.y);
      CHECK(firstPoint.z == secondPoint.z);
    }
  }
}

TEST_CASE("testGithub1990") {
  boost::logging::disable_logs("rdApp.warning");
  SECTION("Mol ops check") {  // we saw the problem here (though it came from
                              // something in MolOps)
    auto mol = "F/C=C/F"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    MolOps::removeHs(*mol);
    CHECK(mol->getNumAtoms() == 4);
    int cid = DGeomHelpers::EmbedMolecule(*mol);
    CHECK(cid >= 0);
  }
  SECTION("Reported Problem") {  // The original problem report
    auto mol =
        "CCCCCCCCCCCCCCCC(=O)O[C@@H]1CC(C)=C(/C=C/C(C)=C/C=C/C(C)=C/"
        "C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\C2=C(C)C[C@@H](OC(=O)CCCCCCCCCCCCCCC)"
        "CC2(C)C)C(C)(C)C1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    MolOps::removeHs(*mol);
    int cid = DGeomHelpers::EmbedMolecule(*mol);
    CHECK(cid >= 0);
  }
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("testGithub2246") {
  using PointVect = std::vector<RDGeom::Point3D>;
  SECTION("simple") {  // make sure the mechanics work
    PointVect pts = {{0, 0, 0}, {1.5, 0, 0}};
    auto m = "C1CC1C"_smiles;
    REQUIRE(m);
    MolOps::addHs(*m);
    REQUIRE(m);
    auto params = DGeomHelpers::ETKDG;
    std::map<int, RDGeom::Point3D> coordMap;
    params.useRandomCoords = true;
    params.coordMap = &coordMap;
    params.maxIterations = 1;
    for (unsigned int i = 0; i < pts.size(); ++i) {
      coordMap[i] = pts[i];
    }
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*m, params);
    CHECK(cid >= 0);
    for (unsigned int i = 0; i < pts.size(); ++i) {
      auto d = (m->getConformer().getAtomPos(i) - pts[i]).length();
      CHECK(d < 1e-3);
    }
  }
  SECTION("complex") {  // a more complex example
    std::vector<RDGeom::Point3D> pts = {
        {0, 0, 0}, {1.5, 0, 0}, {1.5, 1.5, 0}, {0, 1.5, 0}};
    auto m = "C12C3CC1.O2C.C3CC"_smiles;
    REQUIRE(m);
    MolOps::addHs(*m);
    REQUIRE(m);
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

TEST_CASE("testProvideBoundsMatrix") {
  boost::logging::disable_logs("rdApp.warning");
  auto m = "C1CCC1C"_smiles;
  REQUIRE(m);
  auto mat = _getBoundsMatrix(m);

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
  REQUIRE(cid >= 0);

  const auto conf = m->getConformer(cid);
  CHECK(feq((conf.getAtomPos(3) - conf.getAtomPos(0)).length(), 1.2, 0.05));
  CHECK(feq((conf.getAtomPos(3) - conf.getAtomPos(2)).length(), 1.2, 0.05));
  CHECK(feq((conf.getAtomPos(3) - conf.getAtomPos(4)).length(), 1.2, 0.05));
  boost::logging::enable_logs("rdApp.warning");
}

TEST_CASE("testDisableFragmentation") {
  auto m = "OO.OO"_smiles;
  REQUIRE(m);
  MolOps::addHs(*m);
  REQUIRE(m);
  DGeomHelpers::EmbedParameters params;
  params.embedFragmentsSeparately = false;
  params.randomSeed = 0xf00d;
  int cid = DGeomHelpers::EmbedMolecule(*m, params);
  REQUIRE(cid >= 0);

  const auto conf = m->getConformer(cid);

  CHECK((conf.getAtomPos(0) - conf.getAtomPos(2)).length() > 2.0);
  CHECK((conf.getAtomPos(0) - conf.getAtomPos(3)).length() > 2.0);
  CHECK((conf.getAtomPos(1) - conf.getAtomPos(2)).length() > 2.0);
  CHECK((conf.getAtomPos(1) - conf.getAtomPos(3)).length() > 2.0);
}

TEST_CASE("testGithub3019") {
  {  // make sure the mechanics work
    auto m = v2::SmilesParse::MolFromSmiles(std::string(2000, 'C'));
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2000);
    DGeomHelpers::EmbedParameters params;
    params.randomSeed = 0xf00d;
    int cid = DGeomHelpers::EmbedMolecule(*m, params);
    CHECK(cid >= 0);
  }
}

TEST_CASE("testGithub3667") {
  auto throwError = [](unsigned int) {
    throw ValueErrorException("embedder is abortable");
  };
  auto mol =
      "c12c3c4c5c6c1c1c7c8c9c%10c%11c(c28)c3c2c3c4c4c5c5c8c6c1c1c6c7c9c7c9c%"
      "10c%10c%11c2c2c3c3c4c4c5c5c%11c%12c(c1c85)c6c7c1c%12c5c%11c4c3c3c5c(c91)"
      "c%10c23"_smiles;
  REQUIRE(mol);

  DGeomHelpers::EmbedParameters params;
  params.callback = throwError;
  CHECK_THROWS_AS(DGeomHelpers::EmbedMolecule(*mol, params),
                  ValueErrorException);
}

TEST_CASE("testForceTransAmides") {
  auto mol = "CC(=O)NC"_smiles;
  REQUIRE(mol);
  bool updateLabel = true;
  bool takeOwnership = true;
  mol->addAtom(new Atom(1), updateLabel, takeOwnership);
  mol->addBond(3, 5, Bond::BondType::SINGLE);
  MolOps::sanitizeMol(*mol);
  MolOps::addHs(*mol);
  SECTION("Get Trans") {
    DGeomHelpers::EmbedParameters params;
    params.forceTransAmides = true;
    params.randomSeed = 0xf00d;
    params.useExpTorsionAnglePrefs = false;
    params.useBasicKnowledge = true;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, params);
    for (auto cid : cids) {
      REQUIRE(cid >= 0);
      auto conf = mol->getConformer(cid);
      auto tors = MolTransforms::getDihedralDeg(conf, 0, 1, 3, 4);
      CHECK(fabs(fabs(tors) - 180) < 35);
      tors = MolTransforms::getDihedralDeg(conf, 2, 1, 3, 5);
      CHECK(fabs(fabs(tors) - 180) < 35);
    }
  }
  SECTION("Get Cis") {  // make sure we can find at least one non-trans
    DGeomHelpers::EmbedParameters params;
    params.forceTransAmides = false;
    params.randomSeed = 0xf00d;
    params.useExpTorsionAnglePrefs = false;
    params.useBasicKnowledge = true;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, params);
    bool foundOne = false;
    for (auto cid : cids) {
      REQUIRE(cid >= 0);
      auto conf = mol->getConformer(cid);
      auto tors = MolTransforms::getDihedralDeg(conf, 0, 1, 3, 4);
      if (fabs(fabs(tors) - 180) > 50) {
        foundOne = true;
        break;
      }
    }
    CHECK(foundOne);
  }
}

TEST_CASE("testSymmetryPruning") {
  auto mol = "CCOC(C)(C)C"_smiles;
  REQUIRE(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters params;
  params.useSymmetryForPruning = true;
  params.onlyHeavyAtomsForRMS = true;
  params.pruneRmsThresh = 0.5;
  params.randomSeed = 0xf00d;
  auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 50, params);
  CHECK(cids.size() == 2);

  params.useSymmetryForPruning = false;
  cids = DGeomHelpers::EmbedMultipleConfs(*mol, 50, params);
  CHECK(cids.size() == 8);
}

TEST_CASE("testMissingHsWarning") {
  auto mol = "CC"_smiles;
  REQUIRE(mol);

  std::stringstream ss;
  rdWarningLog->SetTee(ss);
  DGeomHelpers::EmbedParameters params;
  DGeomHelpers::EmbedMolecule(*mol, params);
  rdWarningLog->ClearTee();
  CHECK(ss.str().find("Molecule does not have explicit Hs") !=
        std::string::npos);
}

TEST_CASE("testHydrogenBondBasics") {
  auto mol = "CC1O[H]O=C(C)C1 |H:4.3|"_smiles;
  REQUIRE(mol);
  MolOps::addHs(*mol);
  REQUIRE(mol);
  auto mat = _getBoundsMatrix(mol);
  DistGeom::triangleSmoothBounds(mat.get());
  CHECK(mat->getVal(3, 4) > 1.8);
  CHECK(mat->getVal(3, 4) < 2.2);
  CHECK(mat->getVal(4, 3) > 1.0);
  CHECK(mat->getVal(4, 3) < 1.5);

  DGeomHelpers::EmbedParameters params = DGeomHelpers::ETKDGv3;
  params.randomSeed = 0xf00d;
  REQUIRE(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
  auto dist = MolTransforms::getBondLength(mol->getConformer(), 3, 4);
  CHECK(dist < 1.5);
}