//
//  Copyright (C) 2021-2024 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#include <RDGeneral/RDLog.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include "Embedder.h"
#include "BoundsMatrixBuilder.h"
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

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
      CHECK(v1.angleTo(v2) == Catch::Approx(M_PI / 2).margin(0.15));
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

TEST_CASE("problems with bounds matrix smoothing and aromatic sulfur") {
  SECTION("basics") {
    auto core = R"CTAB(test structure - renumbered
     RDKit          3D

  7  7  0  0  0  0  0  0  0  0999 V2000
   48.6842  -14.8137    0.1450 C   0  0  0  0  0  0  0  0  0  0  0  0
   48.0829  -13.5569    0.6868 C   0  0  0  0  0  0  0  0  0  0  0  0
   48.0162  -12.0909   -0.1327 S   0  0  0  0  0  0  0  0  0  0  0  0
   47.1565  -11.3203    1.0899 C   0  0  0  0  0  0  0  0  0  0  0  0
   46.9350  -12.2470    2.1088 C   0  0  0  0  0  0  0  0  0  0  0  0
   46.1942  -11.9293    3.3432 C   0  0  0  0  0  0  0  0  0  0  0  0
   47.4440  -13.4879    1.8745 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  7  2  0
  2  3  1  0
  7  5  1  0
  5  4  2  0
  5  6  1  0
  4  3  1  0
M  END)CTAB"_ctab;
    REQUIRE(core);
    auto thiaz = "Cc1scc(C)n1"_smiles;
    REQUIRE(thiaz);
    MolOps::addHs(*thiaz);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    const auto conf = core->getConformer();
    std::map<int, RDGeom::Point3D> cmap;
    for (unsigned i = 0; i < core->getNumAtoms(); ++i) {
      cmap[i] = conf.getAtomPos(i);
    }
    ps.coordMap = &cmap;
    ps.randomSeed = 0xf00d;
    auto cid = DGeomHelpers::EmbedMolecule(*thiaz, ps);
    CHECK(cid >= 0);
  }
  SECTION("bulk") {
    // run a bunch of molecules with S-containing aromatic heterocycles
    std::vector<std::string> smileses = {
        "[O-][S+](c1ccccn1)c1cncs1",
        "Cn1cccc1C(=O)Nc1nccs1",
        "Cc1csc(=N)n1-c1ccc(Cl)cc1",
        "Nc1ncc([S+]([O-])c2ncccn2)s1",
        "CCCN1CCC=C(c2csc(N)n2)C1",
        "CNc1ncc([S+]([O-])c2ccccn2)s1",
        "Cn1nnnc1SCc1nc2ccccc2s1",
        "CCCC(C(=O)Nc1nccs1)c1ccccc1",
        "Cc1ccc(NC(=O)c2sc(Cl)nc2C)c(C)c1",
        "CCc1nc(-c2ccc(Cl)cc2)sc1C(=O)OC",
        "Cc1nc(CNS(=O)(=O)c2ccc(Cl)cc2)cs1",
        "Cc1ccc2sc(C)[n+](CCC(C)S(=O)(=O)[O-])c2c1",
        "Nc1nc2c(s1)-c1ccccc1Sc1ccccc1-2",
        "COc1ccccc1OCC(=O)Nc1nc(C)c(C)s1",
        "COc1ccc(NC(=O)Nc2sc(=S)n(C)c2C)cc1",
        "C=CCNc1nc(-c2c[nH]c3c(CC)cccc23)cs1",
    };
    auto patt = "[s]1*c[!#6]c1"_smarts;
    REQUIRE(patt);
    for (const auto &smi : smileses) {
      INFO(smi);
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      MolOps::addHs(*mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 0xf00d;
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      REQUIRE(cid >= 0);
      UFF::UFFOptimizeMolecule(*mol);

      auto match = SubstructMatch(*mol, *patt);
      REQUIRE(match.size() >= 1);

      const auto conf = mol->getConformer();
      std::map<int, RDGeom::Point3D> cmap;
      for (auto &mi : match[0]) {
        cmap[mi.second] = conf.getAtomPos(mi.second);
      }
      ps.coordMap = &cmap;
      auto cid2 = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid2 >= 0);
    }
  }
  SECTION("phosphorous") {
    std::vector<std::string> smileses = {
        "CCOC(=O)c1pc(P(Cl)Cl)c2n1[C@@H](C)C(=O)Nc1ccc(C)cc1-2",
        "N(c1c(O)ccc2c(P(Cl)Cl)pc(C(=O)O)n12)[N+](=O)[O-]",
    };
    auto patt = "[p]1*c[!#6]c1"_smarts;
    REQUIRE(patt);
    for (const auto &smi : smileses) {
      INFO(smi);
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      MolOps::addHs(*mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 0xf00d;
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      REQUIRE(cid >= 0);
      UFF::UFFOptimizeMolecule(*mol);

      auto match = SubstructMatch(*mol, *patt);
      REQUIRE(match.size() >= 1);

      const auto conf = mol->getConformer();
      std::map<int, RDGeom::Point3D> cmap;
      for (auto &mi : match[0]) {
        cmap[mi.second] = conf.getAtomPos(mi.second);
      }
      ps.coordMap = &cmap;
      auto cid2 = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid2 >= 0);
    }
  }
}

TEST_CASE("double bond stereo not honored in conformer generator") {
  SECTION("mol 1 basics") {
    // this test used to fail
    // from the platinum set
    auto m = "O=C1OCC/C=C/CC/C=C/C(=N/OCC(=O)N2CCCCC2)Cc2cc(O)cc(O)c21"_smiles;
    REQUIRE(m);
    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0xf00d + 81;
    auto cid = DGeomHelpers::EmbedMolecule(cp, ps);
    REQUIRE(cid >= 0);
    MolOps::assignStereochemistryFrom3D(cp);
    // std::cerr << MolToMolBlock(cp) << std::endl;
    for (const auto bnd : cp.bonds()) {
      if (bnd->getBondType() == Bond::BondType::DOUBLE) {
        INFO(bnd->getIdx());
        CHECK(bnd->getStereo() ==
              m->getBondWithIdx(bnd->getIdx())->getStereo());
      }
    }
  }
  SECTION("mol 1 multiple loops") {
    // from the platinum set
    auto m = "O=C1OCC/C=C/CC/C=C/C(=N/OCC(=O)N2CCCCC2)Cc2cc(O)cc(O)c21"_smiles;
    REQUIRE(m);
    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    for (unsigned int iter = 0; iter < 10; ++iter) {
      RWMol lcp(cp);
      ps.randomSeed = 0xf00d + iter;
      auto cid = DGeomHelpers::EmbedMolecule(lcp, ps);
      REQUIRE(cid >= 0);
      MolOps::assignStereochemistryFrom3D(lcp);
      // std::cerr << MolToMolBlock(cp) << std::endl;
      for (const auto bnd : lcp.bonds()) {
        if (bnd->getBondType() == Bond::BondType::DOUBLE) {
          INFO(iter);
          CHECK(bnd->getStereo() ==
                m->getBondWithIdx(bnd->getIdx())->getStereo());
        }
      }
    }
  }

  SECTION("github #5913") {
    auto m = "[H]/C(F)=C([H])\\C([H])=C(/[H])Br"_smiles;
    REQUIRE(m);

    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    for (unsigned int iter = 0; iter < 50; ++iter) {
      RWMol lcp(cp);
      ps.randomSeed = 0 + iter;
      auto cid = DGeomHelpers::EmbedMolecule(lcp, ps);
      REQUIRE(cid >= 0);
      MolOps::assignStereochemistryFrom3D(lcp);
      // std::cerr << MolToMolBlock(cp) << std::endl;
      for (const auto bnd : lcp.bonds()) {
        if (bnd->getBondType() == Bond::BondType::DOUBLE) {
          INFO(iter);
          CHECK(bnd->getStereo() ==
                m->getBondWithIdx(bnd->getIdx())->getStereo());
        }
      }
    }
  }

  SECTION("github #5283") {
    UseLegacyStereoPerceptionFixture useLegacy(false);
    auto m =
        "Cc3nn(CC(=O)N2CCN(c1ccccc1)CC2)c(C)c3/N=N\\c6ccc(CNC(=O)CCC(=O)Nc4cccc5C(=O)NCc45)cc6"_smiles;
    REQUIRE(m);
    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.enforceChirality = false;
    for (unsigned int iter = 0; iter < 10; ++iter) {
      INFO(iter);
      RWMol lcp(cp);
      ps.randomSeed = 140 + iter;
      auto cid = DGeomHelpers::EmbedMolecule(lcp, ps);
      REQUIRE(cid >= 0);
      MolOps::assignStereochemistryFrom3D(lcp, cid, true);
      auto bnd = lcp.getBondBetweenAtoms(22, 23);
      REQUIRE(bnd);
      REQUIRE(bnd->getBondType() == Bond::BondType::DOUBLE);
      CHECK(bnd->getStereo() == m->getBondWithIdx(bnd->getIdx())->getStereo());
    }
  }
}

TEST_CASE("tracking failure causes"){SECTION("basics"){
    auto mol =
        "C=CC1=C(N)Oc2cc1c(-c1cc(C(C)O)cc(=O)cc1C1NCC(=O)N1)c(OC)c2OC"_smiles;
REQUIRE(mol);
MolOps::addHs(*mol);
DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
ps.randomSeed = 0xf00d;
ps.trackFailures = true;
ps.maxIterations = 50;
ps.randomSeed = 42;
auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
CHECK(cid < 0);

CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::INITIAL_COORDS] > 5);
CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::ETK_MINIMIZATION] > 10);

auto fail_cp = ps.failures;
// make sure we reset the counts each time
cid = DGeomHelpers::EmbedMolecule(*mol, ps);
CHECK(ps.failures == fail_cp);
}
SECTION("chirality") {
  auto mol = R"CTAB(
  Ketcher  1102315302D 1   1.00000     0.00000     0

 10 11  0  0  1  0  0  0  0  0999 V2000
   10.1340  -11.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.1340  -12.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0000  -12.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.8660  -12.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.8660  -11.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0000  -10.5250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.0000  -11.5250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   11.2588  -12.4909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2680  -10.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.7629  -12.4673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  1  0     0  0
  1  2  1  0     0  0
  2  3  1  0     0  0
  3  4  1  0     0  0
  4  5  1  0     0  0
  5  6  1  0     0  0
  1  7  1  0     0  0
  7  8  1  0     0  0
  8  4  1  0     0  0
  1  9  1  1     0  0
  4 10  1  1     0  0
M  END
)CTAB"_ctab;
  REQUIRE(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = 0xf00d;
  ps.trackFailures = true;
  ps.maxIterations = 50;
  auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
  CHECK(cid < 0);
  CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::INITIAL_COORDS] > 5);
  CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::FINAL_CHIRAL_BOUNDS] > 5);
}

#ifdef RDK_TEST_MULTITHREADED
SECTION("multithreaded") {
  auto mol =
      "C=CC1=C(N)Oc2cc1c(-c1cc(C(C)O)cc(=O)cc1C1NCC(=O)N1)c(OC)c2OC"_smiles;
  REQUIRE(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = 0xf00d;
  ps.trackFailures = true;
  ps.maxIterations = 10;
  ps.randomSeed = 42;
  auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 20, ps);

  DGeomHelpers::EmbedParameters ps2 = ps;
  ps2.numThreads = 4;

  auto cids2 = DGeomHelpers::EmbedMultipleConfs(*mol, 20, ps2);
  CHECK(cids2 == cids);

  CHECK(ps.failures == ps2.failures);
}
#endif
}

TEST_CASE("Github #5883: confgen failing for chiral N in a three ring") {
  SECTION("basics1") {
    auto mol = "N1[C@H-]C1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 42;
    ps.maxIterations = 1;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
  }
  SECTION("basics2") {
    auto mol = "N1[N@H]C1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 42;
    ps.maxIterations = 1;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
  }
  SECTION("no ring") {
    auto mol = "N[C@H-]C"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 42;
    ps.maxIterations = 1;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
  }
}

TEST_CASE("Github #6365: cannot generate conformers for PF6- or SF6") {
  SECTION("basics") {
    std::vector<std::string> smileses = {"S(F)(F)(F)(F)(F)F",
                                         "[P-](F)(F)(F)(F)(F)F"};
    for (const auto &smi : smileses) {
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 42;
      ps.useRandomCoords = true;
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid >= 0);
    }
  }
}

TEST_CASE("Sequential random seeds") {
  SECTION("basics") {
    auto mol = "CCCCCCCCCCCC"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);

    RWMol mol2(*mol);

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.enableSequentialRandomSeeds = true;
    ps.useRandomCoords = true;
    ps.randomSeed = 0xf00d;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, ps);
    CHECK(cids.size() == 10);
    ps.randomSeed = 0xf00d + 5;
    auto cids2 = DGeomHelpers::EmbedMultipleConfs(mol2, 5, ps);
    CHECK(cids2.size() == 5);

    compareConfs(mol.get(), &mol2, 5, 0);
  }
}

TEST_CASE("Macrocycle bounds matrix") {
  SECTION("basics") {
    auto mol = "C1/C=C/C=C/CCCCCCCCC1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);

    DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(mol->getNumAtoms())};
    DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
    DGeomHelpers::setTopolBounds(*mol, bm, true, false, true);
    CHECK(bm->getLowerBound(1, 18) > 2.6);
    CHECK(bm->getLowerBound(1, 18) < 2.7);
    CHECK(bm->getLowerBound(4, 17) > 2.6);
    CHECK(bm->getLowerBound(4, 17) < 2.7);

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0;

    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
    const auto conf = mol->getConformer(cid);
    RDGeom::Point3D pos_1 = conf.getAtomPos(1);
    RDGeom::Point3D pos_4 = conf.getAtomPos(4);
    CHECK((pos_1 - pos_4).length() < 3.6);
    CHECK((pos_1 - pos_4).length() > 3.5);
  }
}

TEST_CASE("atropisomers and embedding") {
  SECTION("basics") {
    auto mol =
        "Cc1cccc(O)c1-c1c(N)cccc1Cl |(-8.88571,2.09707,;-8.17143,3.33425,;-6.74286,3.33425,;-6.02857,4.57143,;-6.74286,5.80861,;-8.17143,5.80861,;-8.88571,7.04579,;-8.88571,4.57143,;-10.3143,4.57143,;-11.0286,5.80861,;-10.3143,7.04579,;-12.4571,5.80861,;-13.1714,4.57143,;-12.4571,3.33425,;-11.0286,3.33425,;-10.3143,2.09707,),wU:8.15|"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondWithIdx(7)->getBondType() == Bond::BondType::SINGLE);
    REQUIRE(mol->getBondWithIdx(7)->getStereo() ==
            Bond::BondStereo::STEREOATROPCCW);
    MolOps::addHs(*mol);
    // mol->debugMol(std::cerr);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0xf00d;
    {
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      REQUIRE(cid >= 0);
      const auto conf = mol->getConformer(cid);

      Atropisomers::AtropAtomAndBondVec abvs[2];
      REQUIRE(Atropisomers::getAtropisomerAtomsAndBonds(mol->getBondWithIdx(7),
                                                        abvs, *mol));
      auto pos_1 = conf.getAtomPos(7);
      auto pos_2 = conf.getAtomPos(8);
      auto pos_3 = conf.getAtomPos(1);
      auto pos_4 = conf.getAtomPos(9);
      auto v2 = pos_2 - pos_1;
      auto v3 = pos_3 - pos_1;
      auto v4 = pos_4 - pos_1;
      auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
      CHECK(chiralVol < 0);
    }
    {
      RWMol mol2(*mol);
      mol2.getBondWithIdx(7)->setStereo(Bond::BondStereo::STEREOATROPCW);

      auto cid = DGeomHelpers::EmbedMolecule(mol2, ps);
      REQUIRE(cid >= 0);
      const auto conf = mol2.getConformer(cid);

      Atropisomers::AtropAtomAndBondVec abvs[2];
      REQUIRE(Atropisomers::getAtropisomerAtomsAndBonds(mol2.getBondWithIdx(7),
                                                        abvs, mol2));
      auto pos_1 = conf.getAtomPos(7);
      auto pos_2 = conf.getAtomPos(8);
      auto pos_3 = conf.getAtomPos(1);
      auto pos_4 = conf.getAtomPos(9);
      auto v2 = pos_2 - pos_1;
      auto v3 = pos_3 - pos_1;
      auto v4 = pos_4 - pos_1;
      auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
      CHECK(chiralVol > 0);
    }
  }
}

TEST_CASE("atropisomers bulk") {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/atropisomers.sdf";
  SDMolSupplier sdsup(fname);

  auto params = DGeomHelpers::ETKDGv3;
  params.randomSeed = 0xf00d + 1;

  for (auto i = 0u; i < sdsup.length(); ++i) {
    std::unique_ptr<RWMol> mol(static_cast<RWMol *>(sdsup[i]));
    REQUIRE(mol);
    auto bondIdx = mol->getProp<unsigned int>("atrop bond");
    REQUIRE((mol->getBondWithIdx(bondIdx)->getStereo() ==
                 Bond::BondStereo::STEREOATROPCCW ||
             mol->getBondWithIdx(bondIdx)->getStereo() ==
                 Bond::BondStereo::STEREOATROPCW));
    auto atropInfo = mol->getProp<std::string>("atrop volume");
    std::vector<std::string> tokens;
    boost::split(tokens, atropInfo, boost::is_any_of(" \t"));
    REQUIRE(tokens.size() == 5);
    std::vector<unsigned int> atropAtoms(4);
    for (auto j = 0u; j < 4u; ++j) {
      atropAtoms[j] = std::stol(tokens[j]);
    }
    int vol = std::stol(tokens[4]);

    MolOps::addHs(*mol);
    unsigned int nconfs = 20;
    {
      auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, nconfs, params);
      CHECK(cids.size() == nconfs);
      for (auto cid : cids) {
        const auto conf = mol->getConformer(cid);
        std::vector<RDGeom::Point3D> pts;
        for (auto idx : atropAtoms) {
          pts.push_back(conf.getAtomPos(idx));
        }
        auto v2 = pts[1] - pts[0];
        auto v3 = pts[2] - pts[0];
        auto v4 = pts[3] - pts[0];
        auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
        INFO(cid << MolToV3KMolBlock(*mol, true, cid));
        CHECK(chiralVol * vol > 0);
      }
    }  // now swap the stereo and see if it still works
    mol->getBondWithIdx(bondIdx)->setStereo(
        mol->getBondWithIdx(bondIdx)->getStereo() ==
                Bond::BondStereo::STEREOATROPCCW
            ? Bond::BondStereo::STEREOATROPCW
            : Bond::BondStereo::STEREOATROPCCW);
    {
      auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, nconfs, params);
      CHECK(cids.size() == nconfs);
      for (auto cid : cids) {
        const auto conf = mol->getConformer(cid);
        std::vector<RDGeom::Point3D> pts;
        for (auto idx : atropAtoms) {
          pts.push_back(conf.getAtomPos(idx));
        }
        auto v2 = pts[1] - pts[0];
        auto v3 = pts[2] - pts[0];
        auto v4 = pts[3] - pts[0];
        auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
        INFO(cid << MolToV3KMolBlock(*mol, true, cid));
        CHECK(chiralVol * vol < 0);
      }
    }
  }
}

TEST_CASE(
    "Github #7109: wrong stereochemistry in ring from stereospecific SMILES") {
  SECTION("basics") {
    auto m = "C1[C@H](C#CC#C)CC[C@H](C#CC#C)C1"_smiles;
    REQUIRE(m);
    MolOps::addHs(*m);
    REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    REQUIRE(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::KDG;
    {  // this always worked
      ps.randomSeed = 0xC0FFEE;
      auto cid = DGeomHelpers::EmbedMolecule(*m, ps);
      CHECK(cid >= 0);
      MolOps::assignStereochemistryFrom3D(*m, cid);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
      CHECK(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    }
    {  // this failed
      ps.randomSeed = 0xC0FFEE + 123;
      auto cid = DGeomHelpers::EmbedMolecule(*m, ps);
      CHECK(cid >= 0);
      MolOps::assignStereochemistryFrom3D(*m, cid);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
      CHECK(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    }
  }
}

TEST_CASE("Github #7181: ET terms applied to constrained atoms") {
  SECTION("basics") {
    auto templ =
        "CNc1ccc(OC)cc1 |(-3.3363,0.129414,1.28582;-2.44714,-0.687978,0.507453;-1.11383,-0.29452,0.197587;-0.622766,0.911164,0.645083;0.652332,1.29026,0.350281;1.45603,0.462513,-0.400278;2.7718,0.891528,-0.684446;3.83908,0.0736652,-0.224516;0.984393,-0.734112,-0.850528;-0.300532,-1.12218,-0.556656)|"_smiles;
    REQUIRE(templ);
    auto mol = "COc1ccc(NC(C)C)cc1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    auto matches = SubstructMatch(*mol, *templ);
    REQUIRE(matches.size() == 1);

    auto tconf = templ->getConformer();
    std::map<int, RDGeom::Point3D> cmap;
    for (auto [ti, mi] : matches[0]) {
      cmap[mi] = tconf.getAtomPos(ti);
    }

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0xC0FFEE;
    ps.coordMap = &cmap;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
    auto imatch = matches[0];
    for (auto &[ti, mi] : imatch) {
      std::swap(ti, mi);
    }
    auto rmsd = MolAlign::alignMol(*mol, *templ, cid, -1, &imatch);
    CHECK(rmsd < 0.2);
  }
}

TEST_CASE("terminal groups in pruning") {
  SECTION("basics") {
    std::vector<std::string> smiles = {"FCC(=O)O", "FCC(=O)[O-]",
                                       "FCC(=N)[NH-]", "FCS(=O)(=O)O",
                                       "FCP(=O)(O)O"};
    for (const auto &smi : smiles) {
      auto mol = v2::SmilesParse::MolFromSmiles(smi);
      REQUIRE(mol);
      MolOps::addHs(*mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 0xc0ffee;
      ps.pruneRmsThresh = 0.5;

      auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 50, ps);
      CHECK(cids.size() == 1);

      ps.symmetrizeConjugatedTerminalGroupsForPruning = false;
      cids = DGeomHelpers::EmbedMultipleConfs(*mol, 50, ps);
      CHECK(cids.size() >= 2);
    }
  }
}