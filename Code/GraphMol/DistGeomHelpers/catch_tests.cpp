//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "catch.hpp"

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include "Embedder.h"
#include "BoundsMatrixBuilder.h"
#include <tuple>

using namespace RDKit;

TEST_CASE("Torsions not found in fused macrocycles", "[macrocycles]") {
  RDLog::InitLogs();
  SECTION("reported") {
    // this is 6VY8 from the PDB
    auto mol1 =
        "CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@@H]2CCCN2C(=O)[C@H]([C@@H](C)CC)NC(=O)[C@H](CO)NC(=O)[C@H](CCCC[NH3+])NC(=O)[C@H]([C@@H](C)O)NC(O)[C@@H]2CN3NNC[C@H]3C[C@H](NC1=O)C(O)N[C@@H](Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@@H](CCCNC(N)=[NH2+])C(=O)N2"_smiles;
    REQUIRE(mol1);
    MolOps::addHs(*mol1);
    ForceFields::CrystalFF::CrystalFFDetails details;
    bool useExpTorsions = true;
    bool useSmallRingTorsions = false;
    bool useMacrocycleTorsions = true;
    bool useBasicKnowledge = true;
    unsigned int version = 2;
    bool verbose = true;
    std::stringstream sstrm;
    rdInfoLog->SetTee(sstrm);
    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol1, details, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version, verbose);
    rdInfoLog->ClearTee();
    auto txt = sstrm.str();
    CHECK(txt.find("{9-}") != std::string::npos);
  }
  SECTION("edges") {
    std::vector<std::tuple<std::string, bool, unsigned int>> tests{
        {"O=C1CNC(=O)C2CCC(N1)NC(=O)CNC2=O", true, 15},  // 9-9
        {"O=C1NC2CCC(C(=O)N1)C(=O)NCC(=O)N2", true, 4},  // 9-8
        {"O=C1NC2CCC(C(=O)N1)C(=O)NC(=O)N2", false, 0},  // 8-8
        {"O=C1CC(=O)NC2NC(=O)CC(=O)NC(N1)NC(=O)CC(=O)N2", true,
         18}};  // 12-12-12
    for (const auto &tpl : tests) {
      std::unique_ptr<RWMol> m{SmilesToMol(std::get<0>(tpl))};
      REQUIRE(m);
      MolOps::addHs(*m);
      ForceFields::CrystalFF::CrystalFFDetails details;
      bool useExpTorsions = true;
      bool useSmallRingTorsions = false;
      bool useMacrocycleTorsions = true;
      bool useBasicKnowledge = true;
      unsigned int version = 2;
      bool verbose = true;
      std::stringstream sstrm;
      rdInfoLog->SetTee(sstrm);
      std::cerr << "-----------" << std::endl;
      ForceFields::CrystalFF::getExperimentalTorsions(
          *m, details, useExpTorsions, useSmallRingTorsions,
          useMacrocycleTorsions, useBasicKnowledge, version, verbose);
      rdInfoLog->ClearTee();
      auto txt = sstrm.str();
      if (std::get<1>(tpl)) {
        CHECK(txt.find("{9-}") != std::string::npos);
      } else {
        CHECK(txt.find("{9-}") == std::string::npos);
      }
      CHECK(details.expTorsionAngles.size() == std::get<2>(tpl));
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

TEST_CASE("update parameters from JSON") {
  std::string rdbase = getenv("RDBASE");
  SECTION("DG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol";
    std::unique_ptr<RWMol> ref{MolFileToMol(fname, true, false)};
    REQUIRE(ref);
    std::unique_ptr<RWMol> mol{SmilesToMol("OCCC")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    CHECK(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    std::string json = R"JSON({"randomSeed":42})JSON";
    DGeomHelpers::updateEmbedParametersFromJSON(params, json);
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
  }
  SECTION("ETKDG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    std::unique_ptr<RWMol> ref{MolFileToMol(fname, true, false)};
    REQUIRE(ref);
    std::unique_ptr<RWMol> mol{SmilesToMol("OCCC")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    CHECK(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    std::string json = R"JSON({"randomSeed":42,
    "useExpTorsionAnglePrefs":true,
    "useBasicKnowledge":true})JSON";
    DGeomHelpers::updateEmbedParametersFromJSON(params, json);
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
  }
  SECTION("ETKDGv2") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol";
    std::unique_ptr<RWMol> ref{MolFileToMol(fname, true, false)};
    REQUIRE(ref);
    std::unique_ptr<RWMol> mol{SmilesToMol("n1cccc(C)c1ON")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    CHECK(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    std::string json = R"JSON({"randomSeed":42,
    "useExpTorsionAnglePrefs":true,
    "useBasicKnowledge":true,
    "ETversion":2})JSON";
    DGeomHelpers::updateEmbedParametersFromJSON(params, json);
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
  }

  SECTION("setting atommap") {
    std::unique_ptr<RWMol> mol{SmilesToMol("OCCC")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    {
      DGeomHelpers::EmbedParameters params;
      std::string json = R"JSON({"randomSeed":42,
    "coordMap":{"0":[0,0,0],"1":[0,0,1.5],"2":[0,1.5,1.5]}})JSON";
      DGeomHelpers::updateEmbedParametersFromJSON(params, json);
      CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
      delete params.coordMap;
      auto conf = mol->getConformer();
      auto v1 = conf.getAtomPos(0) - conf.getAtomPos(1);
      auto v2 = conf.getAtomPos(2) - conf.getAtomPos(1);
      CHECK(v1.angleTo(v2) == Approx(M_PI / 2).margin(0.15));
    }
  }
}

TEST_CASE(
    "github #4346: Specified cis/trans stereo being ignored during "
    "conformation generation in macrocycles") {
  SECTION("basics 1") {
    auto m1 = "C1C/C=C/CCCCCCCC1"_smiles;
    REQUIRE(m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::addHs(*m1);
    DGeomHelpers::EmbedParameters params = DGeomHelpers::KDG;
    params.randomSeed = 0xf00d;
    CHECK(DGeomHelpers::EmbedMolecule(*m1, params) != -1);
    MolOps::assignStereochemistryFrom3D(*m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOE);
  }
  SECTION("basics 2") {
    auto m1 = "C1C/C=C\\CCCCCCCC1"_smiles;
    REQUIRE(m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOZ);
    MolOps::addHs(*m1);
    DGeomHelpers::EmbedParameters params = DGeomHelpers::KDG;
    params.randomSeed = 0xf00d;
    CHECK(DGeomHelpers::EmbedMolecule(*m1, params) != -1);
    MolOps::assignStereochemistryFrom3D(*m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOZ);
  }
}
TEST_CASE("nontetrahedral stereo", "[nontetrahedral]") {
  SECTION("bounds matrix basics") {
    {
      auto m = "Cl[Pt@SP1]([35Cl])([36Cl])[37Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 3);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(2))
                ->getIdx() == 4);
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      // std::cerr << *bm << std::endl;
      CHECK(bm->getLowerBound(0, 3) - bm->getLowerBound(0, 2) > 1.0);
      CHECK(bm->getUpperBound(0, 3) - bm->getUpperBound(0, 2) > 1.0);
    }

    {
      auto m = "Cl[Pt@SP1]([35Cl])[36Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 3);
      CHECK(!Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                            m->getAtomWithIdx(2)));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      // std::cerr << *bm << std::endl;
      CHECK(bm->getLowerBound(0, 3) - bm->getLowerBound(0, 2) > 1.0);
      CHECK(bm->getUpperBound(0, 3) - bm->getUpperBound(0, 2) > 1.0);
    }

    {
      // note that things aren't quite as nice here since we don't actually have
      // TBP UFF parameters
      auto m = "Cl[Pt@TB1]([35Cl])([36Cl])([37Cl])[38Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 5);
      CHECK(!Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                            m->getAtomWithIdx(2)));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(5)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(3), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(120, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      CHECK(bm->getLowerBound(0, 5) - bm->getLowerBound(0, 2) > 0.5);
      CHECK(bm->getUpperBound(0, 5) - bm->getUpperBound(0, 2) > 0.5);
      CHECK(bm->getLowerBound(0, 5) - bm->getLowerBound(2, 3) > 0.5);
      CHECK(bm->getUpperBound(0, 5) - bm->getUpperBound(2, 3) > 0.5);
      CHECK(bm->getLowerBound(2, 3) - bm->getLowerBound(0, 2) > 0.5);
      CHECK(bm->getUpperBound(2, 3) - bm->getUpperBound(0, 2) > 0.5);
    }
    {
      auto m = "Cl[Th@OH1]([35Cl])([36Cl])([37Cl])([38Cl])[39Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 6);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(2))
                ->getIdx() == 4);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(3))
                ->getIdx() == 5);

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(6)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(4), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(180, 0.001));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(3), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      CHECK(bm->getLowerBound(0, 6) - bm->getLowerBound(0, 2) > 0.5);
      CHECK(bm->getUpperBound(0, 6) - bm->getUpperBound(0, 3) > 0.5);
      CHECK(bm->getLowerBound(0, 6) - bm->getLowerBound(2, 3) > 0.5);
      CHECK(bm->getUpperBound(0, 6) - bm->getUpperBound(2, 4) < 0.01);
      CHECK(bm->getLowerBound(2, 4) - bm->getLowerBound(2, 3) > 0.5);
    }
  }
#if 1
  SECTION("Embedding") {
    {
      auto m = "Cl[Pt@SP1](<-N)(<-N)[Cl]"_smiles;
      REQUIRE(m);
      m->setProp("_Name", "cis platin");
      MolOps::addHs(*m);
      CHECK(DGeomHelpers::EmbedMolecule(*m) == 0);
      auto mb = MolToV3KMolBlock(*m);
      // std::cerr << mb << std::endl;
      std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
      MolOps::assignStereochemistryFrom3D(*m2);
      CHECK(m2->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_SQUAREPLANAR);
      unsigned int perm = 100;
      CHECK(m2->getAtomWithIdx(1)->getPropIfPresent(
          common_properties::_chiralPermutation, perm));
      CHECK(perm == 1);
    }
    {
      auto m = "Cl[Pt@SP3](<-N)(<-N)[Cl]"_smiles;
      REQUIRE(m);
      m->setProp("_Name", "trans platin");
      MolOps::addHs(*m);
      CHECK(DGeomHelpers::EmbedMolecule(*m) == 0);
      auto mb = MolToV3KMolBlock(*m);
      // std::cerr << mb << std::endl;
      std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
      MolOps::assignStereochemistryFrom3D(*m2);
      CHECK(m2->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_SQUAREPLANAR);
      unsigned int perm = 100;
      CHECK(m2->getAtomWithIdx(1)->getPropIfPresent(
          common_properties::_chiralPermutation, perm));
      CHECK(perm == 3);
    }
  }
#endif
}