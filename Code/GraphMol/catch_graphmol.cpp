//
//  Copyright (C) 2018-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <algorithm>
#include <limits>

#include <boost/format.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit;
#if 1
TEST_CASE("SMILES Parsing works", "[molops]") {
  std::unique_ptr<RWMol> mol(SmilesToMol("C1CC1"));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 3);
}

TEST_CASE("Sanitization tests", "[molops]") {
  std::unique_ptr<RWMol> mol(SmilesToMol("C1=CC=CC=C1Cc2ccccc2", false, false));
  REQUIRE(mol);
  REQUIRE(mol->getNumAtoms() == 13);

  SECTION("properties") {
    mol->updatePropertyCache();
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 1);
    CHECK(!mol->getAtomWithIdx(0)->getIsAromatic());
    CHECK(mol->getAtomWithIdx(10)->getIsAromatic());
    SECTION("aromaticity") {
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed, MolOps::SANITIZE_SETAROMATICITY);
      // mol->debugMol(std::cerr);
      CHECK(mol->getAtomWithIdx(10)->getIsAromatic());
      // blocked by #1730
      // CHECK(mol->getAtomWithIdx(0)->getIsAromatic());
    }
    SECTION("kekulize") {
      unsigned int opThatFailed;
      MolOps::sanitizeMol(*mol, opThatFailed, MolOps::SANITIZE_KEKULIZE);
      CHECK(!mol->getAtomWithIdx(0)->getIsAromatic());
      CHECK(!mol->getAtomWithIdx(10)->getIsAromatic());
    }
  }
}

TEST_CASE("Github #2062", "[bug][molops]") {
  SmilesParserParams ps;
  ps.removeHs = false;
  ps.sanitize = true;
  std::unique_ptr<RWMol> mol(SmilesToMol("[C:1][C:2]([H:3])([H])[O:4][H]", ps));
  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 6);
  mol->getAtomWithIdx(1)->setProp("intProp", 42);
  MolOps::mergeQueryHs(*mol);
  CHECK(mol->getNumAtoms() == 3);
  SECTION("basics") { CHECK(mol->getAtomWithIdx(1)->getAtomMapNum() == 2); }
  SECTION("other props") {
    REQUIRE(mol->getAtomWithIdx(1)->hasProp("intProp"));
    CHECK(mol->getAtomWithIdx(1)->getProp<int>("intProp") == 42);
  }
}

TEST_CASE("Github #2086", "[bug][molops]") {
  SECTION("reported version") {
    auto mol = "C1CCCC1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    REQUIRE(mol->getNumAtoms() == 15);
    mol->removeBond(4, 13);
    MolOps::removeHs(*mol);
    REQUIRE(mol->getNumAtoms() == 6);
  }
}

TEST_CASE("github #299", "[bug][molops][SSSR]") {
  SECTION("simplified") {
    auto mol =
        "C13%13%14.C124%18.C25%13%15.C368%17.C4679.C75%10%17.C8%11%14%16.C9%11%12%18.C%10%12%15%16"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 9);
  }

  SECTION("old example from molopstest") {
    auto mol = "C123C45C11C44C55C22C33C14C523"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 9);
  }

  SECTION("carborane") {
    std::unique_ptr<RWMol> mol(
        SmilesToMol("[B]1234[B]567[B]118[B]229[B]33%10[B]454[B]656[B]711[B]822["
                    "C]933[B]%1045[C]6123",
                    0, false));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 12);
    mol->updatePropertyCache(false);
    MolOps::findSSSR(*mol);
    REQUIRE(mol->getRingInfo()->isInitialized());
  }
  SECTION("original report from ChEbI") {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    std::unique_ptr<RWMol> mol(
        MolFileToMol(pathName + "ChEBI_50252.mol", false));
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 80);
    mol->updatePropertyCache(false);
    MolOps::findSSSR(*mol);
    REQUIRE(mol->getRingInfo()->isInitialized());
  }
}

TEST_CASE("github #2224", "[bug][molops][removeHs][query]") {
  SECTION("the original report") {
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    std::unique_ptr<RWMol> mol(MolFileToMol(pathName + "github2224_1.mol"));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 7);
  }
  SECTION("basics") {
    SmilesParserParams ps;
    ps.removeHs = false;
    ps.sanitize = true;
    std::unique_ptr<ROMol> mol(SmilesToMol("C[H]", ps));
    REQUIRE(mol);
    REQUIRE(mol->getNumAtoms() == 2);
    {  // The H without a query is removed
      std::unique_ptr<ROMol> m2(MolOps::removeHs(*mol));
      CHECK(m2->getNumAtoms() == 1);
    }
    {  // but if we add a query feature it's not removed
      RWMol m2(*mol);
      auto *qa = new QueryAtom(1);
      m2.replaceAtom(1, qa);
      m2.getAtomWithIdx(1)->setAtomicNum(1);
      MolOps::removeHs(m2);
      CHECK(m2.getNumAtoms() == 2);
      delete qa;
    }
  }
}

TEST_CASE(
    "github #2268: Recognize N in three-membered rings as potentially chiral",
    "[bug][stereo]") {
  SECTION("basics: N in a 3 ring") {
    const auto mol = "C[N@]1CC1C"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("basics: N in a 4 ring") {
    const auto mol = "C[N@]1CCC1C"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  }
  SECTION("the original molecule") {
    std::string mb = R"CTAB(
  Mrv1810 02131915062D

 18 20  0  0  1  0            999 V2000
   -0.7207   -1.3415    0.0000 N   0  0  1  0  0  0  0  0  0  0  0  0
   -0.0583   -0.8416    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.0083   -1.7540    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3956   -0.8666    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.3250   -0.0667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1955   -0.6499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1499   -0.0792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6541   -0.4292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7830   -1.2291    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6081   -1.6623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4080    0.1500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3665   -0.8374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6416    0.3958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1996    0.3708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4121    1.1624    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3498    0.8207    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0790   -0.4167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0665    0.4083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  1  3  1  1  0  0  0
  4  1  1  0  0  0  0
  5  2  1  0  0  0  0
  4  6  1  0  0  0  0
  7  4  1  0  0  0  0
  2  8  1  6  0  0  0
  9  6  2  0  0  0  0
  4 10  1  1  0  0  0
 11  6  1  0  0  0  0
 12  8  2  0  0  0  0
 13  8  1  0  0  0  0
 14 11  1  0  0  0  0
 15 14  1  0  0  0  0
 16 13  2  0  0  0  0
 17 12  1  0  0  0  0
 18 16  1  0  0  0  0
  2  3  1  0  0  0  0
  5  7  1  0  0  0  0
 17 18  2  0  0  0  0
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
}

TEST_CASE("github #2244", "[bug][molops][stereo]") {
  SECTION("the original report") {
    auto mol = "CC=CC=CC"_smiles;
    REQUIRE(mol);
    MolOps::findPotentialStereoBonds(*mol, true);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOANY);
    CHECK(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOANY);
    mol->getBondWithIdx(3)->setStereo(Bond::STEREONONE);
    MolOps::findPotentialStereoBonds(*mol, true);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::STEREOANY);
    CHECK(mol->getBondWithIdx(3)->getStereo() == Bond::STEREOANY);
  }
}

TEST_CASE(
    "github #2258: heterocycles with exocyclic bonds not failing valence check",
    "[bug][molops]") {
  SECTION("the original report") {
    std::vector<std::string> smiles = {"C=n1ccnc1", "C#n1ccnc1"};
    for (auto smi : smiles) {
      CHECK_THROWS_AS(SmilesToMol(smi), MolSanitizeException);
    }
  }
}

TEST_CASE("github #908: AddHs() using 3D coordinates with 2D conformations",
          "[bug][molops]") {
  SECTION("basics: single atom mols") {
    std::vector<std::string> smiles = {"Cl", "O", "N", "C"};
    for (auto smi : smiles) {
      // std::cerr << smi << std::endl;
      std::unique_ptr<RWMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      auto conf = new Conformer(1);
      conf->set3D(false);
      conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
      mol->addConformer(conf, true);
      bool explicitOnly = false;
      bool addCoords = true;
      MolOps::addHs(*mol, explicitOnly, addCoords);
      for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
        // std::cerr << "   " << i << " " << conf->getAtomPos(i) << std::endl;
        CHECK(conf->getAtomPos(i).z == 0.0);
      }
    }
  }
}

TEST_CASE(
    "github #2437: Canon::rankMolAtoms results in crossed double bonds in "
    "rings",
    "[bug][molops]") {
  SECTION("underlying problem") {
    std::string molb = R"CTAB(testmol
  Mrv1824 05081910082D

  4  4  0  0  0  0            999 V2000
    6.9312   -8.6277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9312   -9.4527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7562   -8.6277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7562   -9.4527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  3  4  1  0  0  0  0
  2  4  2  0  0  0  0
M  END
    )CTAB";
    bool sanitize = false;
    bool removeHs = false;
    std::unique_ptr<RWMol> mol(MolBlockToMol(molb, sanitize, removeHs));
    REQUIRE(mol);
    mol->updatePropertyCache();
    CHECK(mol->getBondWithIdx(3)->getBondType() == Bond::BondType::DOUBLE);
    CHECK(mol->getBondWithIdx(3)->getBondDir() == Bond::BondDir::NONE);
    std::vector<unsigned int> ranks;
    CHECK(!mol->getRingInfo()->isInitialized());
    Canon::rankMolAtoms(*mol, ranks);
  }

  SECTION("as discovered") {
    std::string molb = R"CTAB(testmol
  Mrv1824 05081910082D

  4  4  0  0  0  0            999 V2000
    6.9312   -9.4527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7562   -8.6277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.7562   -9.4527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9312   -8.6277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  3  4  1  0  0  0  0
  2  4  2  0  0  0  0
M  END
    )CTAB";
    bool sanitize = false;
    bool removeHs = false;
    std::unique_ptr<RWMol> mol(MolBlockToMol(molb, sanitize, removeHs));
    REQUIRE(mol);
    mol->updatePropertyCache();
    CHECK(mol->getBondWithIdx(3)->getBondType() == Bond::BondType::DOUBLE);
    CHECK(mol->getBondWithIdx(3)->getBondDir() == Bond::BondDir::NONE);
    auto nmb = MolToMolBlock(*mol);
    CHECK(nmb.find("2  4  2  3") == std::string::npos);
    CHECK(nmb.find("2  4  2  0") != std::string::npos);
    std::vector<unsigned int> ranks;
    Canon::rankMolAtoms(*mol, ranks);
    nmb = MolToMolBlock(*mol);
    CHECK(nmb.find("2  4  2  3") == std::string::npos);
    CHECK(nmb.find("2  4  2  0") != std::string::npos);
  }
}
TEST_CASE(
    "github #2423: Incorrect assignment of explicit Hs to Al+3 read from mol "
    "block",
    "[bug][molops]") {
  SECTION("basics: single atom mols") {
    std::string mb = R"CTAB(2300
  -OEChem-01301907122D

  1  0  0     0  0  0  0  0  0999 V2000
  -66.7000  999.0000    0.0000 Al  0  1  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   3
M  END)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getFormalCharge() == 3);
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 0);
  }
}

TEST_CASE("Specialized exceptions for sanitization errors", "[molops]") {
  SECTION("AtomValenceException") {
    std::vector<std::pair<std::string, unsigned int>> smiles = {
        {"C=n1ccnc1", 1}, {"CCO(C)C", 2}};
    for (auto pr : smiles) {
      CHECK_THROWS_AS(SmilesToMol(pr.first), AtomValenceException);
      try {
        auto m = SmilesToMol(pr.first);
        RDUNUSED_PARAM(m);
      } catch (const AtomValenceException &e) {
        CHECK(e.getType() == "AtomValenceException");
        CHECK(e.getAtomIdx() == pr.second);
      }
    }
  }
  SECTION("AtomKekulizeException") {
    std::vector<std::pair<std::string, unsigned int>> smiles = {
        {"CCcc", 2},
    };
    for (auto pr : smiles) {
      CHECK_THROWS_AS(SmilesToMol(pr.first), AtomKekulizeException);
      try {
        auto m = SmilesToMol(pr.first);
        RDUNUSED_PARAM(m);
      } catch (const AtomKekulizeException &e) {
        CHECK(e.getType() == "AtomKekulizeException");
        CHECK(e.getAtomIdx() == pr.second);
      }
    }
  }
  SECTION("KekulizeException") {
    std::vector<std::pair<std::string, std::vector<unsigned int>>> smiles = {
        {"c1cccc1", {0, 1, 2, 3, 4}},
        {"Cc1cc1", {1, 2, 3}},
        {"C1:c:CC1", {0, 1, 2}}};
    for (auto pr : smiles) {
      CHECK_THROWS_AS(SmilesToMol(pr.first), KekulizeException);
      try {
        auto m = SmilesToMol(pr.first);
        RDUNUSED_PARAM(m);
      } catch (const KekulizeException &e) {
        CHECK(e.getType() == "KekulizeException");
        CHECK(e.getAtomIndices() == pr.second);
      }
    }
  }
}

TEST_CASE("detectChemistryProblems", "[molops]") {
  SECTION("Basics") {
    SmilesParserParams ps;
    ps.sanitize = false;
    auto m = std::unique_ptr<ROMol>(SmilesToMol("CO(C)CFCc1cc1", ps));
    REQUIRE(m);
    auto res = MolOps::detectChemistryProblems(*m);
    REQUIRE(res.size() == 3);

    CHECK(res[0]->getType() == "AtomValenceException");
    REQUIRE(dynamic_cast<AtomValenceException *>(res[0].get()));
    CHECK(dynamic_cast<AtomSanitizeException *>(res[0].get())->getAtomIdx() ==
          1);

    CHECK(res[1]->getType() == "AtomValenceException");
    REQUIRE(dynamic_cast<AtomSanitizeException *>(res[1].get()));
    CHECK(dynamic_cast<AtomSanitizeException *>(res[1].get())->getAtomIdx() ==
          4);

    CHECK(res[2]->getType() == "KekulizeException");
    REQUIRE(dynamic_cast<KekulizeException *>(res[2].get()));
    CHECK(dynamic_cast<KekulizeException *>(res[2].get())->getAtomIndices() ==
          std::vector<unsigned int>({6, 7, 8}));
  }
  SECTION("No problems") {
    SmilesParserParams ps;
    ps.sanitize = false;
    auto m = std::unique_ptr<ROMol>(SmilesToMol("c1ccccc1", ps));
    REQUIRE(m);
    auto res = MolOps::detectChemistryProblems(*m);
    REQUIRE(res.size() == 0);
  }
}

TEST_CASE(
    "github #2606: Bad valence corrections on Pb, Sn"
    "[bug][molops]") {
  SECTION("basics-Pb") {
    std::string mb = R"CTAB(
  Mrv1810 08141905562D

  5  0  0  0  0  0            999 V2000
   -3.6316   -0.4737    0.0000 Pb  0  0  0  0  0  0  0  0  0  0  0  0
   -3.6541    0.3609    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.4586   -0.5188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -3.6992   -1.5338    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -4.5789   -0.4286    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  CHG  5   1   4   2  -1   3  -1   4  -1   5  -1
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getFormalCharge() == 4);
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 0);
  }
  SECTION("basics-Sn") {
    std::string mb = R"CTAB(
  Mrv1810 08141905562D

  5  0  0  0  0  0            999 V2000
   -3.6316   -0.4737    0.0000 Sn  0  0  0  0  0  0  0  0  0  0  0  0
   -3.6541    0.3609    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.4586   -0.5188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -3.6992   -1.5338    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -4.5789   -0.4286    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  CHG  5   1   4   2  -1   3  -1   4  -1   5  -1
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getFormalCharge() == 4);
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 0);
  }
  SECTION("basics-Ge") {
    std::string mb = R"CTAB(
  Mrv1810 08141905562D

  5  0  0  0  0  0            999 V2000
   -3.6316   -0.4737    0.0000 Ge  0  0  0  0  0  0  0  0  0  0  0  0
   -3.6541    0.3609    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.4586   -0.5188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -3.6992   -1.5338    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -4.5789   -0.4286    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  CHG  5   1   4   2  -1   3  -1   4  -1   5  -1
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getFormalCharge() == 4);
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 0);
  }
}
TEST_CASE(
    "github #2607: Pb, Sn should support valence 2"
    "[bug][molops]") {
  SECTION("basics-Pb") {
    std::string mb = R"CTAB(
  Mrv1810 08141905562D

  3  0  0  0  0  0            999 V2000
   -3.6316   -0.4737    0.0000 Pb  0  0  0  0  0  0  0  0  0  0  0  0
   -3.6541    0.3609    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.4586   -0.5188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  CHG  3   1   2   2  -1   3  -1
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getFormalCharge() == 2);
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 0);
  }
  SECTION("basics-Sn") {
    std::string mb = R"CTAB(
  Mrv1810 08141905562D

  3  0  0  0  0  0            999 V2000
   -3.6316   -0.4737    0.0000 Sn  0  0  0  0  0  0  0  0  0  0  0  0
   -3.6541    0.3609    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.4586   -0.5188    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  CHG  3   1   2   2  -1   3  -1
M  END
)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->getFormalCharge() == 2);
    CHECK(mol->getAtomWithIdx(0)->getTotalNumHs() == 0);
  }
}

TEST_CASE(
    "github #2649: Allenes read from mol blocks have crossed bonds assigned"
    "[bug][stereochemistry]") {
  SECTION("basics") {
    std::string mb = R"CTAB(mol
  Mrv1824 09191901002D

  6  5  0  0  0  0            999 V2000
   -1.6986   -7.4294    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2522   -6.8245    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1438   -8.0357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8095   -6.2156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3374   -7.8470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6162   -6.3886    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  2  0  0  0  0
  2  1  2  0  0  0  0
  3  5  1  0  0  0  0
  4  2  2  0  0  0  0
  6  4  1  0  0  0  0
M  END)CTAB";
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);
    CHECK(mol->getBondWithIdx(3)->getStereo() == Bond::STEREONONE);
    auto outmolb = MolToMolBlock(*mol);
    // std::cerr<<outmolb<<std::endl;
    CHECK(outmolb.find("1  3  2  0") != std::string::npos);
    CHECK(outmolb.find("2  1  2  0") != std::string::npos);
    CHECK(outmolb.find("4  2  2  0") != std::string::npos);
  }
}

TEST_CASE(
    "GitHub 2712: setBondStereoFromDirections() returning incorrect results"
    "[stereochemistry]") {
  SECTION("basics 1a") {
    std::string mb = R"CTAB(
  Mrv1810 10141909562D

  4  3  0  0  0  0            999 V2000
    3.3412   -2.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5162   -2.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1037   -3.7112    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7537   -2.2823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
)CTAB";
    bool sanitize = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb, sanitize));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    MolOps::setBondStereoFromDirections(*mol);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOTRANS);
  }
  SECTION("basics 1b") {
    std::string mb = R"CTAB(
  Mrv1810 10141909562D

  4  3  0  0  0  0            999 V2000
    3.3412   -2.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5162   -2.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1037   -3.7112    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7537   -2.2823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  4  1  1  0  0  0  0
M  END
)CTAB";
    bool sanitize = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb, sanitize));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    MolOps::setBondStereoFromDirections(*mol);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOTRANS);
  }
  SECTION("basics 2a") {
    std::string mb = R"CTAB(
  Mrv1810 10141909582D

  4  3  0  0  0  0            999 V2000
    3.4745   -5.2424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6495   -5.2424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2370   -5.9569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8870   -5.9569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
)CTAB";
    bool sanitize = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb, sanitize));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    MolOps::setBondStereoFromDirections(*mol);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOCIS);
  }
  SECTION("basics 2b") {
    std::string mb = R"CTAB(
  Mrv1810 10141909582D

  4  3  0  0  0  0            999 V2000
    3.4745   -5.2424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6495   -5.2424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2370   -5.9569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8870   -5.9569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  4  1  1  0  0  0  0
M  END
)CTAB";
    bool sanitize = false;
    std::unique_ptr<ROMol> mol(MolBlockToMol(mb, sanitize));
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(0)->getBondType() == Bond::DOUBLE);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREONONE);
    MolOps::setBondStereoFromDirections(*mol);
    CHECK(mol->getBondWithIdx(0)->getStereo() == Bond::STEREOCIS);
  }
}

TEST_CASE("removeHs screwing up double bond stereo", "[bug][removeHs]") {
  SECTION("example1") {
    std::string molblock = R"CTAB(molblock = """
  SciTegic12221702182D

 47 51  0  0  0  0            999 V2000
    0.2962    6.2611    0.0000 C   0  0
   -3.9004    4.4820    0.0000 C   0  0
    1.4195    5.2670    0.0000 C   0  0
   -3.8201   -7.4431    0.0000 C   0  0
   -4.9433   -6.4490    0.0000 C   0  0
   -2.3975   -6.9674    0.0000 C   0  0
    3.5921   -3.5947    0.0000 C   0  0
   -3.1475    2.3700    0.0000 C   0  0
    2.1695   -4.0705    0.0000 C   0  0
   -2.0242    1.3759    0.0000 C   0  0
   -4.6440   -4.9792    0.0000 C   0  0
    2.7681   -1.1308    0.0000 C   0  0
   -5.8626    1.1332    0.0000 C   0  0
    3.0674    0.3391    0.0000 C   0  0
    3.6660    3.2787    0.0000 C   0  0
    8.1591   -0.6978    0.0000 C   0  0
    7.3351    1.7662    0.0000 C   0  0
   -6.3876    3.5028    0.0000 C   0  0
   -0.6756   -5.0219    0.0000 C   0  0
    7.0358    0.2964    0.0000 C   0  0
    3.8914   -2.1249    0.0000 C   0  0
   -2.0982   -5.4976    0.0000 C   0  0
   -4.5701    1.8943    0.0000 C   0  0  1  0  0  0
   -6.9859    2.1273    0.0000 C   0  0  1  0  0  0
    4.4900    0.8148    0.0000 C   0  0
    1.3455   -1.6065    0.0000 C   0  0
    4.7893    2.2846    0.0000 C   0  0
    1.9442    1.3332    0.0000 C   0  0
    1.0462   -3.0763    0.0000 C   0  0
    2.2435    2.8030    0.0000 C   0  0
   -0.6017    1.8516    0.0000 C   0  0
    5.6132   -0.1794    0.0000 C   0  0
    0.2223   -0.6124    0.0000 Cl  0  0
    9.2823   -1.6919    0.0000 N   0  0
   -3.2215   -4.5035    0.0000 N   0  0
    6.2119    2.7603    0.0000 N   0  0
    5.3139   -1.6492    0.0000 N   0  0
    0.5216    0.8575    0.0000 N   0  0
   -4.8945    3.3588    0.0000 N   0  0
   -8.2913    2.8662    0.0000 O   0  0
   -0.3024    3.3214    0.0000 O   0  0
    1.1202    3.7971    0.0000 O   0  0
   -0.3763   -3.5520    0.0000 O   0  0
   -2.8482    3.8398    0.0000 H   0  0
   -2.3235   -0.0940    0.0000 H   0  0
   -3.9483    0.5292    0.0000 H   0  0
   -7.8572    0.9063    0.0000 H   0  0
  1  3  1  0
  2 39  1  0
  3 42  1  0
  4  5  2  0
  4  6  1  0
  5 11  1  0
  6 22  2  0
  7  9  2  0
  7 21  1  0
  8 44  1  0
  8 10  2  0
  8 23  1  0
  9 29  1  0
 10 45  1  0
 10 31  1  0
 11 35  2  0
 12 21  2  0
 12 26  1  0
 13 23  1  0
 13 24  1  0
 14 25  2  0
 14 28  1  0
 15 27  2  0
 15 30  1  0
 16 20  1  0
 16 34  3  0
 17 20  2  0
 17 36  1  0
 18 24  1  0
 18 39  1  0
 19 22  1  0
 19 43  1  0
 20 32  1  0
 21 37  1  0
 22 35  1  0
 23 46  1  6
 23 39  1  0
 24 47  1  1
 24 40  1  0
 25 27  1  0
 25 32  1  0
 26 29  2  0
 26 33  1  0
 27 36  1  0
 28 30  2  0
 28 38  1  0
 29 43  1  0
 30 42  1  0
 31 38  2  0
 31 41  1  0
 32 37  2  3
M  END
"""

)CTAB";
    bool sanitize = false;
    bool removeHs = false;
    std::unique_ptr<RWMol> m(MolBlockToMol(molblock, sanitize, removeHs));
    REQUIRE(m);
    m->updatePropertyCache();
    MolOps::setBondStereoFromDirections(*m);
    CHECK(m->getBondWithIdx(10)->getBondType() == Bond::DOUBLE);
    CHECK(m->getBondWithIdx(10)->getStereo() == Bond::STEREOTRANS);
    REQUIRE(m->getBondWithIdx(10)->getStereoAtoms().size() == 2);
    CHECK(m->getBondWithIdx(10)->getStereoAtoms()[0] == 43);
    CHECK(m->getBondWithIdx(10)->getStereoAtoms()[1] == 44);

    MolOps::removeHs(*m);  // implicitOnly,updateExplicitCount,sanitize);
    // m->debugMol(std::cerr);
    CHECK(m->getBondWithIdx(9)->getBondType() == Bond::DOUBLE);
    CHECK(m->getBondWithIdx(9)->getStereo() == Bond::STEREOTRANS);
    REQUIRE(m->getBondWithIdx(9)->getStereoAtoms().size() == 2);
    CHECK(m->getBondWithIdx(9)->getStereoAtoms()[0] == 22);
    CHECK(m->getBondWithIdx(9)->getStereoAtoms()[1] == 30);
  }
}

TEST_CASE("setDoubleBondNeighborDirections()", "[stereochemistry][bug]") {
  SECTION("basics cis") {
    auto m = "CC=CC"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->getStereoAtoms() = {0, 3};
    m->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    MolOps::setDoubleBondNeighborDirections(*m);
    CHECK(m->getBondWithIdx(0)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::ENDDOWNRIGHT);
    CHECK(MolToSmiles(*m) == "C/C=C\\C");
  }
  SECTION("basics trans") {
    auto m = "CC=CC"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->getStereoAtoms() = {0, 3};
    m->getBondWithIdx(1)->setStereo(Bond::STEREOTRANS);
    MolOps::setDoubleBondNeighborDirections(*m);
    CHECK(m->getBondWithIdx(0)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(MolToSmiles(*m) == "C/C=C/C");
  }
  SECTION("swap (Github #3322)") {
    auto m = "CC=CC"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->getStereoAtoms() = {0, 3};
    m->getBondWithIdx(1)->setStereo(Bond::STEREOTRANS);
    MolOps::setDoubleBondNeighborDirections(*m);
    CHECK(m->getBondWithIdx(0)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(MolToSmiles(*m) == "C/C=C/C");

    m->clearComputedProps();
    m->getBondWithIdx(1)->getStereoAtoms() = {0, 3};
    m->getBondWithIdx(1)->setStereo(Bond::STEREOCIS);
    MolOps::setDoubleBondNeighborDirections(*m);
    CHECK(m->getBondWithIdx(0)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::ENDDOWNRIGHT);
    CHECK(MolToSmiles(*m) == "C/C=C\\C");
  }
}

TEST_CASE("github #2782: addHs() fails on atoms with 'bad' valences", "[bug]") {
  SECTION("basics") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m(
        static_cast<RWMol *>(SmilesToMol("C=C1=CC=CC=C1", ps)));
    REQUIRE(m);
    bool strict = false;
    m->updatePropertyCache(strict);
    CHECK(m->getNumAtoms() == 7);
    MolOps::addHs(*m);
    CHECK(m->getNumAtoms() == 14);
    // this doesn't change the fact that there's still a bad valence present:
    CHECK_THROWS_AS(m->updatePropertyCache(), AtomValenceException);
  }
}

TEST_CASE(
    "Github #2784: Element symbol lookup for some transuranics returns "
    "incorrect results",
    "[transuranics][bug]") {
  auto pt = PeriodicTable::getTable();
  SECTION("number to symbol") {
    std::vector<std::pair<unsigned int, std::string>> data = {
        {113, "Nh"}, {114, "Fl"}, {115, "Mc"},
        {116, "Lv"}, {117, "Ts"}, {118, "Og"}};
    for (const auto &pr : data) {
      CHECK(pt->getElementSymbol(pr.first) == pr.second);
    }
  }
  SECTION("symbol to number") {
    std::vector<std::pair<int, std::string>> data = {
        {113, "Nh"}, {114, "Fl"}, {115, "Mc"},  {116, "Lv"},
        {117, "Ts"}, {118, "Og"}, {113, "Uut"}, {115, "Uup"}};
    for (const auto &pr : data) {
      CHECK(pt->getAtomicNumber(pr.second) == pr.first);
    }
  }
}
TEST_CASE("github #2775", "[valence][bug]") {
  SECTION("basics") {
    std::string molblock = R"CTAB(bismuth citrate
  Mrv1810 11111908592D

 14 12  0  0  0  0            999 V2000
    7.4050   -0.5957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6906   -1.0082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9761   -0.5957    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6906   -1.8332    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    7.4050    0.2293    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5800    0.2293    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1675    0.9438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.5800    1.6583    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.3425    0.9438    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    8.2300    0.2293    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6425   -0.4851    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.6425    0.9438    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    7.4050    1.0543    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.5175    0.9438    0.0000 Bi  0  1  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
  1  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
  5 10  1  0  0  0  0
 10 11  2  0  0  0  0
 10 12  1  0  0  0  0
  5 13  1  0  0  0  0
M  CHG  4   4  -1   9  -1  12  -1  14   3
M  END
)CTAB";
    std::unique_ptr<RWMol> m(MolBlockToMol(molblock));
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(13)->getSymbol() == "Bi");
    CHECK(m->getAtomWithIdx(13)->getNumImplicitHs() == 0);
  }
}

TEST_CASE("RemoveHsParameters", "[molops]") {
  SmilesParserParams smilesPs;
  smilesPs.removeHs = false;

  SECTION("H-H") {
    std::unique_ptr<RWMol> m{SmilesToMol("[H][H].[H]O[H]", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 5);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 3);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 3);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeOnlyHNeighbors = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 1);
    }
  }

  SECTION("dummies") {
    std::unique_ptr<RWMol> m{SmilesToMol("[H][*]O[H]", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 4);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 3);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 3);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeDummyNeighbors = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 2);
    }
  }

  SECTION("chiralHs") {
    std::unique_ptr<RWMol> m{SmilesToMol("[C@]12([H])CCC1CO2", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 7);
    // artificial wedging since we don't have a conformer
    m->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BEGINWEDGE);

    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 6);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 6);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeWithWedgedBond = false;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 7);
    }
  }

  SECTION("degree zero") {
    std::unique_ptr<RWMol> m{SmilesToMol("[F-].[H+]", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 2);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 2);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeDegreeZero = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 1);
    }
  }

  SECTION("isotopes") {
    std::unique_ptr<RWMol> m{SmilesToMol("F[2H]", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 2);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 2);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeIsotopes = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 1);
    }
  }

  SECTION("defining bond stereo") {
    std::unique_ptr<RWMol> m{SmilesToMol("F/C=N/[H]", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 4);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 4);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 4);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeDefiningBondStereo = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 3);
    }
  }
  SECTION("Query atoms") {
    std::unique_ptr<RWMol> m{SmartsToMol("O[#1]")};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 2);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 2);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeWithQuery = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 1);
    }
  }
  SECTION("higher degree") {
    // this is a silly example
    std::unique_ptr<RWMol> m{SmilesToMol("F[H-]F", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 3);
    }
    {
      MolOps::RemoveHsParameters ps;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 3);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeHigherDegrees = true;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 2);
    }
  }
  SECTION("mapped Hs") {
    std::unique_ptr<RWMol> m{SmilesToMol("[H:1]O[H]", smilesPs)};
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    {
      RWMol cp(*m);
      MolOps::removeHs(cp);
      CHECK(cp.getNumAtoms() == 1);
    }
    {
      MolOps::RemoveHsParameters ps;
      ps.removeMapped = false;
      RWMol cp(*m);
      MolOps::removeHs(cp, ps);
      CHECK(cp.getNumAtoms() == 2);
    }
  }
  SECTION("allHs") {
    std::unique_ptr<RWMol> m{SmilesToMol(
        "[C@]12([H])CCC1CO2.[H+].F[H-]F.[H][H].[H]*.F/C=C/[H]", smilesPs)};
    REQUIRE(m);
    // artificial wedging since we don't have a conformer
    m->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BEGINWEDGE);
    RWMol cp(*m);
    MolOps::removeAllHs(cp);
    for (auto atom : cp.atoms()) {
      CHECK(atom->getAtomicNum() != 1);
    }
  }
  SECTION("allHs2") {
    std::unique_ptr<ROMol> m{SmilesToMol(
        "[C@]12([H])CCC1CO2.[H+].F[H-]F.[H][H].[H]*.F/C=C/[H]", smilesPs)};
    REQUIRE(m);
    // artificial wedging since we don't have a conformer
    m->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BEGINWEDGE);
    std::unique_ptr<ROMol> cp{MolOps::removeAllHs(*m)};
    for (auto atom : cp->atoms()) {
      CHECK(atom->getAtomicNum() != 1);
    }
  }
}
#endif
TEST_CASE("github #2895: acepentalene aromaticity perception ",
          "[molops][bug][aromaticity]") {
  SECTION("acepentalene") {
    std::unique_ptr<RWMol> m{SmilesToMol("C1=CC2=CC=C3C2=C1C=C3")};
    REQUIRE(m);
    auto smi = MolToSmiles(*m);
    CHECK(smi == "C1=CC2=C3C1=CC=C3C=C2");
  }
}

TEST_CASE("github #3256: fused ring aromaticity perception",
          "[molops][bug][aromaticity]") {
  SECTION("nitrogen only central ring") {
    auto mol = "C1=CN2C3=CC=CN3C3=CC=CN3C2=C1"_smiles;
    REQUIRE(mol);
    for (const auto b : mol->bonds()) {
      CHECK(b->getBondType() == Bond::AROMATIC);
    }
    auto smi = MolToSmiles(*mol);
    CHECK(smi == "c1cc2n(c1)c1cccn1c1cccn21");
  }
}

TEST_CASE("phosphine and arsine chirality", "[Chirality]") {
  SECTION("chiral center recognized") {
    auto mol1 = "C[P@](C1CCCC1)C1=CC=CC=C1"_smiles;
    auto mol2 = "C[As@](C1CCCC1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(mol2->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center selective") {
    auto mol1 = "C[P@](C)C1CCCCC1"_smiles;
    auto mol2 = "C[As@](C)C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(mol2->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center specific: P") {
    auto mol1 = "C[P@](C1CCCC1)C1=CC=CC=C1"_smiles;
    auto mol2 = "C[P@@](C1CCCC1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
  SECTION("chiral center specific: As") {
    auto mol1 = "C[As@](C1CCCC1)C1=CC=CC=C1"_smiles;
    auto mol2 = "C[As@@](C1CCCC1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
  SECTION("chiral center, implicit H: P") {
    auto mol1 = "C[P@H]C1CCCCC1"_smiles;
    auto mol2 = "C[P@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center, implicit H: As") {
    auto mol1 = "C[As@H]C1CCCCC1"_smiles;
    auto mol2 = "C[As@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center specific, implicit H: P") {
    auto mol1 = "C[P@H]C1CCCCC1"_smiles;
    auto mol2 = "C[P@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
  SECTION("chiral center specific, implicit H: As") {
    auto mol1 = "C[As@H]C1CCCCC1"_smiles;
    auto mol2 = "C[As@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
}

TEST_CASE("github #2890", "[bug][molops][stereo]") {
  auto mol = "CC=CC"_smiles;
  REQUIRE(mol);

  auto bond = mol->getBondWithIdx(1);
  bond->setStereo(Bond::STEREOANY);
  REQUIRE(bond->getStereoAtoms().empty());

  MolOps::findPotentialStereoBonds(*mol);
  CHECK(bond->getStereo() == Bond::STEREOANY);
  CHECK(bond->getStereoAtoms().size() == 2);
}

TEST_CASE("github #3150 MolOps::removeHs removes hydrides", "[bug][molops]") {
  SmilesParserParams smilesPs;
  smilesPs.removeHs = false;

  SECTION("Hydride ion remove Hydrides false") {
    std::unique_ptr<RWMol> m{SmilesToMol("[H-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = false;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    // H atom not removed in this case because by default H atoms with degree 0
    // are not removed
    CHECK(cp.getNumAtoms() == 1);
    CHECK(MolOps::getFormalCharge(cp) == -1);
  }

  SECTION("Hydride ion remove Hydrides true") {
    std::unique_ptr<RWMol> m{SmilesToMol("[H-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = true;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    // H atom not removed in this case because by default H atoms with degree 0
    // are not removed
    CHECK(cp.getNumAtoms() == 1);
    CHECK(MolOps::getFormalCharge(cp) == -1);
  }

  SECTION("Water") {
    std::unique_ptr<RWMol> m{SmilesToMol("[OH+][H-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = false;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    CHECK(cp.getNumAtoms() == 2);
    CHECK(MolOps::getFormalCharge(cp) == 0);
  }

  SECTION("Water remove Hydrides true") {
    std::unique_ptr<RWMol> m{SmilesToMol("[OH+][H-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = true;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    CHECK(cp.getNumAtoms() == 1);
    CHECK(MolOps::getFormalCharge(cp) == 1);
  }

  SECTION("Iron Hydride") {
    std::unique_ptr<RWMol> m{SmilesToMol("[Fe+2]<-[H-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = false;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    CHECK(cp.getNumAtoms() == 2);
    CHECK(MolOps::getFormalCharge(cp) == 1);
  }

  SECTION("Iron Hydride remove Hydrides") {
    std::unique_ptr<RWMol> m{SmilesToMol("[Fe+2]<-[H-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = true;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    CHECK(cp.getNumAtoms() == 1);
    CHECK(MolOps::getFormalCharge(cp) == 2);
  }

  SECTION("Ferrous Hydroxide") {
    std::unique_ptr<RWMol> m{SmilesToMol("[Fe+2]<-[OH-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = false;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    CHECK(cp.getNumAtoms() == 2);
    CHECK(MolOps::getFormalCharge(cp) == 1);
  }

  SECTION("Ferrous Hydroxide remove Hydrides") {
    std::unique_ptr<RWMol> m{SmilesToMol("[Fe+2]<-[OH-]", smilesPs)};
    REQUIRE(m);
    MolOps::RemoveHsParameters ps;
    ps.removeHydrides = true;
    RWMol cp(*m);
    MolOps::removeHs(cp, ps);
    CHECK(cp.getNumAtoms() == 2);
    CHECK(MolOps::getFormalCharge(cp) == 1);
  }

  SECTION("Remove All Hs in Hydrides Ferrous Hydride") {
    std::unique_ptr<RWMol> m{SmilesToMol("[Fe+2]<-[H-]", smilesPs)};
    REQUIRE(m);

    RWMol cp(*m);
    MolOps::removeAllHs(cp);
    CHECK(cp.getNumAtoms() == 1);
    CHECK(MolOps::getFormalCharge(cp) == 2);
  }

  SECTION("Remove All Hs in Hydrides Water") {
    std::unique_ptr<RWMol> m{SmilesToMol("[OH+][H-]", smilesPs)};
    REQUIRE(m);

    RWMol cp(*m);
    MolOps::removeAllHs(cp);
    CHECK(cp.getNumAtoms() == 1);
    CHECK(MolOps::getFormalCharge(cp) == 1);
  }
}

TEST_CASE("hybridization of unknown atom types", "[bug][molops]") {
  SECTION("Basics") {
    auto m = "[U][U][U]"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getHybridization() ==
          Atom::HybridizationType::S);
    CHECK(m->getAtomWithIdx(1)->getHybridization() ==
          Atom::HybridizationType::SP);
    CHECK(m->getAtomWithIdx(2)->getHybridization() ==
          Atom::HybridizationType::S);
  }
  SECTION("comprehensive") {
    std::string smiles = "F";
    for (unsigned int i = 89; i <= 118; ++i) {
      smiles += (boost::format("[#%d]") % i).str();
    }
    smiles += "F";
    std::unique_ptr<ROMol> m(SmilesToMol(smiles));
    REQUIRE(m);
    for (auto i = 1u; i < m->getNumAtoms() - 1; ++i) {
      CHECK(m->getAtomWithIdx(i)->getHybridization() ==
            Atom::HybridizationType::SP);
    }
  }
}

TEST_CASE("Github #3470: Hydrogen is incorrectly identified as an early atom",
          "[bug][chemistry]") {
  SECTION("Basics") {
    RWMol m;
    m.addAtom(new Atom(1), true, true);
    m.getAtomWithIdx(0)->setFormalCharge(-1);
    m.updatePropertyCache();
    CHECK(m.getAtomWithIdx(0)->getNumImplicitHs() == 0);
    m.getAtomWithIdx(0)->setFormalCharge(1);
    m.updatePropertyCache();
    CHECK(m.getAtomWithIdx(0)->getNumImplicitHs() == 0);
    m.getAtomWithIdx(0)->setFormalCharge(0);
    m.updatePropertyCache();
    CHECK(m.getAtomWithIdx(0)->getNumImplicitHs() == 1);

    // make sure we still generate errors for stupid stuff
    m.getAtomWithIdx(0)->setFormalCharge(-2);
    CHECK_THROWS_AS(m.updatePropertyCache(), AtomValenceException);
    CHECK(m.getAtomWithIdx(0)->getNumImplicitHs() == 1);
  }
  SECTION("confirm with SMILES") {
    RWMol m;
    bool updateLabel = false;
    bool takeOwnership = true;
    m.addAtom(new Atom(1), updateLabel, takeOwnership);
    m.getAtomWithIdx(0)->setFormalCharge(-1);
    m.updatePropertyCache();
    CHECK(MolToSmiles(m) == "[H-]");
    m.getAtomWithIdx(0)->setFormalCharge(+1);
    m.updatePropertyCache();
    CHECK(MolToSmiles(m) == "[H+]");
    m.getAtomWithIdx(0)->setFormalCharge(0);
    m.updatePropertyCache();
    CHECK(MolToSmiles(m) == "[HH]");  // ugly, but I think [H] would be worse
  }
}

TEST_CASE("Additional oxidation states", "[chemistry]") {
  SECTION("Basics") {
    std::vector<std::string> smiles = {"F[Po](F)(F)(F)", "F[Po](F)(F)(F)(F)F",
                                       "F[Xe](F)(F)(F)", "F[Xe](F)(F)(F)(F)F",
                                       "F[I](F)F",       "F[I](F)(F)(F)F",
                                       "F[At](F)F",      "F[At](F)(F)(F)F"};
    for (const auto &smi : smiles) {
      std::unique_ptr<ROMol> m(SmilesToMol(smi));
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    }
  }
}

TEST_CASE("Github #3805: radicals on [He]", "[chemistry]") {
  SECTION("Basics") {
    {
      auto m = "[He]"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
      CHECK(m->getAtomWithIdx(0)->getTotalNumHs() == 0);
    }
    {
      auto m = "[Ne]"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
      CHECK(m->getAtomWithIdx(0)->getTotalNumHs() == 0);
    }
  }
  SECTION("Basics") {
    {
      auto m = "[He+]"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
      CHECK(m->getAtomWithIdx(0)->getTotalNumHs() == 0);
    }
    {
      auto m = "[Ne+]"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
      CHECK(m->getAtomWithIdx(0)->getTotalNumHs() == 0);
    }
  }
}

TEST_CASE("needsHs function", "[chemistry]") {
  SECTION("basics") {
    const auto m = "CC"_smiles;
    REQUIRE(m);
    CHECK(MolOps::needsHs(*m));

    // add a single H:
    bool updateLabel = false;
    bool takeOwnership = true;
    m->addAtom(new Atom(1), updateLabel, takeOwnership);
    m->addBond(0, 2, Bond::BondType::SINGLE);
    MolOps::sanitizeMol(*m);
    CHECK(MolOps::needsHs(*m));

    // now add all the Hs:
    MolOps::addHs(*m);
    CHECK(!MolOps::needsHs(*m));
  }
  SECTION("radical") {
    const auto m = "[O][O]"_smiles;
    REQUIRE(m);
    CHECK(!MolOps::needsHs(*m));
  }
  SECTION("none needed") {
    const auto m = "FF"_smiles;
    REQUIRE(m);
    CHECK(!MolOps::needsHs(*m));
  }
}

TEST_CASE(
    "github #3330: incorrect number of radicals electrons calculated for "
    "metals",
    "[chemistry][metals]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned int>> data = {
        {"[Mn+2]", 1}, {"[Mn+1]", 0}, {"[Mn]", 1}, {"[Mn-1]", 0},
        {"[C]", 4},    {"[C+1]", 3},  {"[C-1]", 3}};
    for (const auto &pr : data) {
      std::unique_ptr<ROMol> m(SmilesToMol(pr.first));
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == pr.second);
    }
  }
}

TEST_CASE("github #3879: bad H coordinates on fused rings", "[addhs]") {
  SECTION("reported") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.500000 2.598076 0.000000 0
M  V30 2 N 0.750000 1.299038 0.000000 0
M  V30 3 C 1.500000 -0.000000 0.000000 0
M  V30 4 C 0.750000 -1.299038 0.000000 0
M  V30 5 C 0.382772 -0.562069 0.000000 0
M  V30 6 C -0.295379 0.612525 0.000000 0
M  V30 7 C -0.750000 1.299038 0.000000 0
M  V30 8 C -1.500000 0.000000 0.000000 0
M  V30 9 O -0.750000 -1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 3 CFG=3
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 7 6 CFG=3
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 7 2
M  V30 10 1 9 4
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    bool explicitOnly = false;
    bool addCoords = true;
    UINT_VECT onlyOnAtoms = {3, 6};
    MolOps::addHs(*m, explicitOnly, addCoords, &onlyOnAtoms);
    const auto &conf = m->getConformer();
    // check that the H atoms bisect the angle correctly
    {
      REQUIRE(m->getAtomWithIdx(9)->getAtomicNum() == 1);
      REQUIRE(m->getBondBetweenAtoms(9, 3));
      REQUIRE(m->getBondBetweenAtoms(3, 4));
      REQUIRE(m->getBondBetweenAtoms(3, 2));
      REQUIRE(m->getBondBetweenAtoms(3, 8));
      auto v1 = conf.getAtomPos(9) - conf.getAtomPos(3);
      auto v2 = conf.getAtomPos(4) - conf.getAtomPos(3);
      auto v3 = conf.getAtomPos(2) - conf.getAtomPos(3);
      auto v4 = conf.getAtomPos(8) - conf.getAtomPos(3);
      CHECK(v1.angleTo(v3) < v1.angleTo(v2));
      CHECK(v1.angleTo(v4) < v1.angleTo(v2));
      CHECK(fabs(v1.angleTo(v3) - v1.angleTo(v4)) < 1e-4);
      CHECK(v1.dotProduct(v3) < -1e-4);
      CHECK(v1.dotProduct(v4) < -1e-4);
    }
    {
      REQUIRE(m->getAtomWithIdx(10)->getAtomicNum() == 1);
      REQUIRE(m->getBondBetweenAtoms(10, 6));
      REQUIRE(m->getBondBetweenAtoms(5, 6));
      REQUIRE(m->getBondBetweenAtoms(6, 1));
      REQUIRE(m->getBondBetweenAtoms(6, 7));
      auto v1 = conf.getAtomPos(10) - conf.getAtomPos(6);
      auto v2 = conf.getAtomPos(5) - conf.getAtomPos(6);
      auto v3 = conf.getAtomPos(1) - conf.getAtomPos(6);
      auto v4 = conf.getAtomPos(7) - conf.getAtomPos(6);
      CHECK(v1.angleTo(v3) < v1.angleTo(v2));
      CHECK(v1.angleTo(v4) < v1.angleTo(v2));
      CHECK(fabs(v1.angleTo(v3) - v1.angleTo(v4)) < 1e-4);
      CHECK(v1.dotProduct(v3) < -1e-4);
      CHECK(v1.dotProduct(v4) < -1e-4);
    }
  }
  SECTION("non-chiral version") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.500000 2.598076 0.000000 0
M  V30 2 N 0.750000 1.299038 0.000000 0
M  V30 3 C 1.500000 -0.000000 0.000000 0
M  V30 4 C 0.750000 -1.299038 0.000000 0
M  V30 5 C 0.382772 -0.562069 0.000000 0
M  V30 6 C -0.295379 0.612525 0.000000 0
M  V30 7 C -0.750000 1.299038 0.000000 0
M  V30 8 C -1.500000 0.000000 0.000000 0
M  V30 9 O -0.750000 -1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 3
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 7 6
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 7 2
M  V30 10 1 9 4
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    bool explicitOnly = false;
    bool addCoords = true;
    UINT_VECT onlyOnAtoms = {3, 6};
    MolOps::addHs(*m, explicitOnly, addCoords, &onlyOnAtoms);
    const auto &conf = m->getConformer();
    {
      REQUIRE(m->getAtomWithIdx(9)->getAtomicNum() == 1);
      REQUIRE(m->getBondBetweenAtoms(9, 3));
      REQUIRE(m->getBondBetweenAtoms(3, 4));
      REQUIRE(m->getBondBetweenAtoms(3, 2));
      REQUIRE(m->getBondBetweenAtoms(3, 8));
      auto v1 = conf.getAtomPos(9) - conf.getAtomPos(3);
      auto v2 = conf.getAtomPos(4) - conf.getAtomPos(3);
      auto v3 = conf.getAtomPos(2) - conf.getAtomPos(3);
      auto v4 = conf.getAtomPos(8) - conf.getAtomPos(3);
      CHECK(v1.angleTo(v3) < v1.angleTo(v2));
      CHECK(v1.angleTo(v4) < v1.angleTo(v2));
      CHECK(fabs(v1.angleTo(v3) - v1.angleTo(v4)) < 1e-4);
      CHECK(v1.dotProduct(v3) < -1e-4);
      CHECK(v1.dotProduct(v4) < -1e-4);
    }
    {
      REQUIRE(m->getAtomWithIdx(10)->getAtomicNum() == 1);
      REQUIRE(m->getBondBetweenAtoms(10, 6));
      REQUIRE(m->getBondBetweenAtoms(5, 6));
      REQUIRE(m->getBondBetweenAtoms(6, 1));
      REQUIRE(m->getBondBetweenAtoms(6, 7));
      auto v1 = conf.getAtomPos(10) - conf.getAtomPos(6);
      auto v2 = conf.getAtomPos(5) - conf.getAtomPos(6);
      auto v3 = conf.getAtomPos(1) - conf.getAtomPos(6);
      auto v4 = conf.getAtomPos(7) - conf.getAtomPos(6);
      CHECK(v1.angleTo(v3) < v1.angleTo(v2));
      CHECK(v1.angleTo(v4) < v1.angleTo(v2));
      CHECK(fabs(v1.angleTo(v3) - v1.angleTo(v4)) < 1e-4);
      CHECK(v1.dotProduct(v3) < -1e-4);
      CHECK(v1.dotProduct(v4) < -1e-4);
    }
  }
  SECTION("a simpler system") {
    auto m = R"CTAB(
  Mrv2014 03092106042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.3533 6.6867 0 0
M  V30 2 C -4.3533 5.1467 0 0 CFG=1
M  V30 3 O -2.8133 6.6867 0 0
M  V30 4 C -2.8133 5.1467 0 0 CFG=1
M  V30 5 C -3.5833 3.813 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 1 2 4
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 5 1 4 5 CFG=1
M  V30 6 1 2 1 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    bool explicitOnly = false;
    bool addCoords = true;
    UINT_VECT onlyOnAtoms = {3, 1};
    MolOps::addHs(*m, explicitOnly, addCoords, &onlyOnAtoms);
    const auto &conf = m->getConformer();
    {
      REQUIRE(m->getAtomWithIdx(5)->getAtomicNum() == 1);
      REQUIRE(m->getBondBetweenAtoms(5, 1));
      REQUIRE(m->getBondBetweenAtoms(1, 3));
      REQUIRE(m->getBondBetweenAtoms(1, 0));
      REQUIRE(m->getBondBetweenAtoms(1, 4));
      auto v1 = conf.getAtomPos(5) - conf.getAtomPos(1);
      auto v2 = conf.getAtomPos(3) - conf.getAtomPos(1);
      auto v3 = conf.getAtomPos(0) - conf.getAtomPos(1);
      auto v4 = conf.getAtomPos(4) - conf.getAtomPos(1);
      CHECK(v1.angleTo(v3) < v1.angleTo(v2));
      CHECK(v1.angleTo(v4) < v1.angleTo(v2));
      CHECK(fabs(v1.angleTo(v3) - v1.angleTo(v4)) < 1e-4);
      CHECK(v1.dotProduct(v3) < -1e-4);
      CHECK(v1.dotProduct(v4) < -1e-4);
    }
  }
  SECTION("#3932: followup from #3879") {
    auto m = R"CTAB(
     RDKit          2D

 21 22  0  0  0  0  0  0  0  0999 V2000
   -6.9959    0.0617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5212    0.3365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0219    1.7509    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5460   -0.8032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0713   -0.5284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0961   -1.6681    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6214   -1.3933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2270   -0.1562    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4640   -1.0046    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.6531   -0.0903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0918   -0.5150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1788    0.5186    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4434   -1.9732    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1510    1.3231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6092    1.6747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9538    2.8101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6516    1.2824    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.7678    2.7779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8236    1.5543    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6156   -2.2417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8904   -3.7163    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  1  0
  8  7  1  6
  8  9  1  0
  9 10  1  0
 10 11  1  6
 11 12  2  0
 11 13  1  0
 10 14  1  0
 14 15  1  0
 14 16  1  0
 14 17  1  0
 17 18  2  0
 17 19  2  0
  9 20  1  0
 20 21  2  0
 20  7  1  0
 17  8  1  0
M  END)CTAB"_ctab;
    REQUIRE(m);

    bool explicitOnly = false;
    bool addCoords = true;
    UINT_VECT onlyOnAtoms = {7};
    MolOps::addHs(*m, explicitOnly, addCoords, &onlyOnAtoms);
    const auto &conf = m->getConformer();
    {
      REQUIRE(m->getAtomWithIdx(21)->getAtomicNum() == 1);
      REQUIRE(m->getBondBetweenAtoms(21, 7));
      REQUIRE(m->getBondBetweenAtoms(7, 8));
      REQUIRE(m->getBondBetweenAtoms(7, 16));
      REQUIRE(m->getBondBetweenAtoms(7, 6));
      auto v1 = conf.getAtomPos(21) - conf.getAtomPos(7);
      auto v2 = conf.getAtomPos(6) - conf.getAtomPos(7);
      auto v3 = conf.getAtomPos(16) - conf.getAtomPos(7);
      auto v4 = conf.getAtomPos(8) - conf.getAtomPos(7);
      CHECK(v1.angleTo(v2) < v1.angleTo(v4));
      CHECK(v1.angleTo(v3) < v1.angleTo(v4));
      CHECK(fabs(v1.angleTo(v2) - v1.angleTo(v3)) < 1e-4);
    }
  }
}

TEST_CASE("batch edits", "[editing]") {
  SECTION("removeAtom") {
    auto m = "C1CCCO1"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->removeAtom(2);
    m->removeAtom(3);
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "CCO");
  }
  SECTION("removeAtom + removeBond") {
    auto m = "C1CCCO1"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->removeAtom(3);
    m->removeBond(4, 0);
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "CCC.O");
  }
  SECTION("rollback") {
    auto m = "C1CCCO1"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->removeAtom(2);
    m->removeAtom(3);
    m->rollbackBatchEdit();
    CHECK(MolToSmiles(*m) == "C1CCOC1");
  }
  SECTION("adding atoms while in a batch") {
    auto m = "CCCO"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->removeAtom(2);
    bool updateLabel = false;
    bool takeOwnership = true;
    m->addAtom(new Atom(7), updateLabel, takeOwnership);
    m->removeAtom(1);
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "C.N.O");
  }
  SECTION("removing added atoms while in a batch") {
    auto m = "CCCO"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->removeAtom(2);
    bool updateLabel = false;
    bool takeOwnership = true;
    m->addAtom(new Atom(7), updateLabel, takeOwnership);
    m->removeAtom(4);
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "CC.O");
  }
  SECTION("adding bonds while in a batch") {
    auto m = "CCCO"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->removeBond(2, 3);
    m->addBond(0, 3, Bond::BondType::SINGLE);
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "CCCO");
  }
  SECTION("removing added bonds while in a batch") {
    auto m = "CCCO"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    m->addBond(0, 3, Bond::BondType::SINGLE);
    m->removeBond(2, 3);
    m->removeBond(0, 3);
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "CCC.O");
  }
  SECTION("some details") {
    auto m = "CCCO"_smiles;
    REQUIRE(m);
    m->beginBatchEdit();
    CHECK_THROWS_AS(m->beginBatchEdit(), ValueErrorException);
    m->removeAtom(0U);
    // copying includes the edit status:
    RWMol m2(*m);
    CHECK_THROWS_AS(m2.beginBatchEdit(), ValueErrorException);

    // without a commit, the mols haven't changed
    CHECK(MolToSmiles(*m) == "CCCO");
    CHECK(MolToSmiles(m2) == "CCCO");
    m->commitBatchEdit();
    CHECK(MolToSmiles(*m) == "CCO");
    m2.commitBatchEdit();
    CHECK(MolToSmiles(m2) == "CCO");
  }
}

TEST_CASE("github #4122: segfaults in commitBatchEdit()", "[editing]][bug]") {
  SECTION("as reported, no atoms") {
    RWMol m;
    m.beginBatchEdit();
    m.addAtom();
    m.commitBatchEdit();
  }
  SECTION("no bonds") {
    auto m = "C.C"_smiles;
    m->beginBatchEdit();
    m->addBond(0, 1, Bond::BondType::SINGLE);
    m->commitBatchEdit();
  }
  SECTION("after add atom") {
    auto m = "CC"_smiles;
    m->beginBatchEdit();
    m->addAtom(6);
    m->removeAtom(0u);
    m->addAtom(6);
    m->commitBatchEdit();
  }
  SECTION("remove a just-added atom") {
    auto m = "CC"_smiles;
    m->beginBatchEdit();
    m->addAtom(6);
    m->removeAtom(2);
    m->commitBatchEdit();
  }
}

TEST_CASE("github #3912: cannot draw atom lists from SMARTS", "[query][bug]") {
  SECTION("original") {
    auto m = "C(-[N,O])-[#7,#8]"_smarts;
    REQUIRE(m);
    CHECK(isAtomListQuery(m->getAtomWithIdx(1)));
    CHECK(isAtomListQuery(m->getAtomWithIdx(2)));

    std::vector<int> vals;
    getAtomListQueryVals(m->getAtomWithIdx(2)->getQuery(), vals);
    CHECK(vals == std::vector<int>{7, 8});
    vals.clear();
    getAtomListQueryVals(m->getAtomWithIdx(1)->getQuery(), vals);
    CHECK(vals == std::vector<int>{7, 8});
  }
}

TEST_CASE("github #4496: cannot draw aromatic atom lists from SMARTS",
          "[query][bug]") {
  SECTION("original") {
    auto m = "[c,n]1[c,n][c,n][c,n][c,n][c,n]1"_smarts;
    REQUIRE(m);
    std::vector<int> expected({6, 7});
    for (const auto a : m->atoms()) {
      CHECK(isAtomListQuery(a));
      std::vector<int> vals;
      getAtomListQueryVals(a->getQuery(), vals);
      CHECK(vals == expected);
    }
  }
}

TEST_CASE("bridgehead queries", "[query]") {
  SECTION("basics") {
    {
      auto m = "CC12CCN(CC1)C2"_smiles;
      REQUIRE(m);
      for (const auto atom : m->atoms()) {
        auto test = queryIsAtomBridgehead(atom);
        if (atom->getIdx() == 1 || atom->getIdx() == 4) {
          CHECK(test == 1);
        } else {
          CHECK(test == 0);
        }
      }
    }
    {
      auto m = "CC12CCC(C)(CC1)CC2"_smiles;
      REQUIRE(m);
      for (const auto atom : m->atoms()) {
        auto test = queryIsAtomBridgehead(atom);
        if (atom->getIdx() == 1 || atom->getIdx() == 4) {
          CHECK(test == 1);
        } else {
          CHECK(test == 0);
        }
      }
    }
    {  // no bridgehead
      auto m = "C1CCC2CCCCC2C1"_smiles;
      REQUIRE(m);
      for (const auto atom : m->atoms()) {
        auto test = queryIsAtomBridgehead(atom);
        CHECK(test == 0);
      }
    }
  }
  SECTION("Github #6049") {
    auto m = "C1C=CC=C2CCCCC3CC(C3)N21"_smiles;
    REQUIRE(m);
    CHECK(!queryIsAtomBridgehead(m->getAtomWithIdx(13)));
    CHECK(!queryIsAtomBridgehead(m->getAtomWithIdx(4)));
    CHECK(queryIsAtomBridgehead(m->getAtomWithIdx(11)));
    CHECK(queryIsAtomBridgehead(m->getAtomWithIdx(9)));
  }
}

TEST_CASE("replaceAtom/Bond should not screw up bookmarks", "[RWMol]") {
  SECTION("atom basics") {
    auto m = "CCC"_smiles;
    REQUIRE(m);
    m->setAtomBookmark(m->getAtomWithIdx(2), 1);
    auto origAt2 = m->getAtomWithIdx(2);
    CHECK(m->getUniqueAtomWithBookmark(1) == origAt2);
    Atom O(8);
    m->replaceAtom(2, &O);
    auto at2 = m->getAtomWithIdx(2);
    CHECK(at2 != origAt2);
    CHECK(m->getUniqueAtomWithBookmark(1) == at2);
  }
  SECTION("bond basics") {
    auto m = "CCCC"_smiles;
    REQUIRE(m);
    m->setBondBookmark(m->getBondWithIdx(2), 1);
    auto origB2 = m->getBondWithIdx(2);
    CHECK(m->getUniqueBondWithBookmark(1) == origB2);
    Bond single(Bond::BondType::SINGLE);
    m->replaceBond(2, &single);
    auto b2 = m->getBondWithIdx(2);
    CHECK(b2 != origB2);
    CHECK(m->getUniqueBondWithBookmark(1) == b2);
  }
}

TEST_CASE("github #4071: StereoGroups not preserved by RenumberAtoms()",
          "[molops]") {
  SECTION("basics") {
    auto mol =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,o1:7,&1:1,&2:5,r|"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getStereoGroups().size() == 4);

    std::vector<unsigned int> aindices(mol->getNumAtoms());
    std::iota(aindices.begin(), aindices.end(), 0);
    std::reverse(aindices.begin(), aindices.end());
    std::unique_ptr<ROMol> nmol(MolOps::renumberAtoms(*mol, aindices));
    REQUIRE(nmol);
    CHECK(nmol->getStereoGroups().size() == 4);
    for (size_t i = 0; i < nmol->getStereoGroups().size(); ++i) {
      CHECK(nmol->getStereoGroups()[i].getGroupType() ==
            mol->getStereoGroups()[i].getGroupType());
    }
    CHECK(MolToCXSmiles(*nmol) ==
          "C[C@H]([C@@H](C)[C@@H](C)O)[C@@H](C)O |o1:7,&1:1,&2:2,&3:4|");
  }
}

TEST_CASE("github #4127: SEGV in ROMol::getAtomDegree if atom is not in graph",
          "[graphmol]") {
  // also includes tests for some related edge cases found as part of that bug
  // fix
  Atom atom(6);
  RWMol mol1;
  auto mol2 = "CCC"_smiles;
  SECTION("getAtomDegree") {
    CHECK_THROWS_AS(mol1.getAtomDegree(nullptr), Invar::Invariant);
    CHECK_THROWS_AS(mol1.getAtomDegree(&atom), Invar::Invariant);
    CHECK_THROWS_AS(mol1.getAtomDegree(mol2->getAtomWithIdx(0)),
                    Invar::Invariant);
  }
  SECTION("getAtomNeighbors") {
    CHECK_THROWS_AS(mol1.getAtomNeighbors(nullptr), Invar::Invariant);
    CHECK_THROWS_AS(mol1.getAtomNeighbors(&atom), Invar::Invariant);
    CHECK_THROWS_AS(mol1.getAtomNeighbors(mol2->getAtomWithIdx(0)),
                    Invar::Invariant);
  }
  SECTION("getAtomBonds") {
    CHECK_THROWS_AS(mol1.getAtomBonds(nullptr), Invar::Invariant);
    CHECK_THROWS_AS(mol1.getAtomBonds(&atom), Invar::Invariant);
    CHECK_THROWS_AS(mol1.getAtomBonds(mol2->getAtomWithIdx(0)),
                    Invar::Invariant);
  }
  SECTION("addAtom from another molecule") {
    RWMol mol1cp(mol1);
    CHECK_THROWS_AS(mol1cp.addAtom(nullptr), Invar::Invariant);
    bool updateLabel = false;
    bool takeOwnership = true;
    CHECK_THROWS_AS(
        mol1cp.addAtom(mol2->getAtomWithIdx(0), updateLabel, takeOwnership),
        Invar::Invariant);
    takeOwnership = false;
    CHECK(mol1cp.addAtom(mol2->getAtomWithIdx(0), updateLabel, takeOwnership) ==
          0);
  }
  SECTION("addBond from another molecule") {
    auto mol3 = "C.C.C"_smiles;
    bool takeOwnership = true;
    CHECK_THROWS_AS(mol3->addBond(mol2->getBondWithIdx(0), takeOwnership),
                    Invar::Invariant);
    takeOwnership = false;
    CHECK(mol3->addBond(mol2->getBondWithIdx(0), takeOwnership) == 1);
  }
}

TEST_CASE(
    "github #4128: SEGV from unsigned integer overflow in "
    "Conformer::setAtomPos",
    "[graphmol]") {
  Conformer conf;
  RDGeom::Point3D pt(0, 0, 0);
  CHECK_THROWS_AS(conf.setAtomPos(std::numeric_limits<unsigned>::max(), pt),
                  ValueErrorException);
}

TEST_CASE("KekulizeFragment", "[graphmol]") {
  SECTION("basics") {
    auto mol = "CCc1ccccc1"_smiles;
    REQUIRE(mol);
    boost::dynamic_bitset<> atomsInPlay(mol->getNumAtoms());
    for (auto aidx : std::vector<size_t>{0, 1, 2, 3}) {
      atomsInPlay.set(aidx);
    }
    boost::dynamic_bitset<> bondsInPlay(mol->getNumBonds());
    for (auto bidx : std::vector<size_t>{0, 1, 2}) {
      bondsInPlay.set(bidx);
    }
    MolOps::details::KekulizeFragment(*mol, atomsInPlay, bondsInPlay);
    CHECK(!mol->getAtomWithIdx(2)->getIsAromatic());
    CHECK(mol->getAtomWithIdx(4)->getIsAromatic());
    CHECK(!mol->getBondWithIdx(2)->getIsAromatic());
    // at the moment that bond still has an aromatic bond order, which isn't
    // optimal, but that will have to wait until we add a feature to allow
    // kekulization of conjugated chains.
    CHECK(mol->getBondWithIdx(2)->getBondType() == Bond::AROMATIC);
    CHECK(mol->getBondWithIdx(4)->getIsAromatic());
  }
}

TEST_CASE(
    "github #4266: fallback ring finding failing on molecules with multiple "
    "fragments",
    "[graphmol]") {
  SECTION("case1") {
    auto m = "C123C45C16C21C34C561.c1ccccc1"_smiles;
    REQUIRE(m);
    ROMol m2(*m);
    m2.getRingInfo()->reset();
    MolOps::fastFindRings(m2);
    CHECK(m->getRingInfo()->numRings() == m2.getRingInfo()->numRings());
  }
  SECTION("case2") {
    auto m = "c1ccccc1.C123C45C16C21C34C561"_smiles;
    REQUIRE(m);
    ROMol m2(*m);
    m2.getRingInfo()->reset();
    MolOps::fastFindRings(m2);
    CHECK(m->getRingInfo()->numRings() == m2.getRingInfo()->numRings());
  }
}

TEST_CASE("QueryBond valence contribs") {
  {
    auto m = "CO"_smarts;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(0)->getValenceContrib(m->getAtomWithIdx(0)) == 0.0);
    CHECK(m->getBondWithIdx(0)->getValenceContrib(m->getAtomWithIdx(1)) == 0.0);
  }
  {
    auto m = "C-O"_smarts;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(0)->getValenceContrib(m->getAtomWithIdx(0)) == 1.0);
    CHECK(m->getBondWithIdx(0)->getValenceContrib(m->getAtomWithIdx(1)) == 1.0);
  }
}

TEST_CASE(
    "github #4311: unreasonable calculation of implicit valence for atoms with "
    "query bonds",
    "[graphmol]") {
  SECTION("basics") {
    auto m = "C-,=O"_smarts;
    REQUIRE(m);
    m->updatePropertyCache();
    CHECK(m->getAtomWithIdx(0)->getTotalNumHs() == 0);
    CHECK(m->getAtomWithIdx(1)->getTotalNumHs() == 0);
    CHECK(MolToSmiles(*m) == "CO");
    CHECK(MolToSmarts(*m) == "C-,=O");
  }
}

TEST_CASE("atom copy ctor") {
  auto m = "CO"_smiles;
  REQUIRE(m);
  for (const auto atom : m->atoms()) {
    Atom cp(*atom);
    CHECK(cp.getAtomicNum() == atom->getAtomicNum());
    CHECK(!cp.hasOwningMol());
  }
}

TEST_CASE("bond copy ctor") {
  auto m = "COC"_smiles;
  REQUIRE(m);
  for (const auto bond : m->bonds()) {
    Bond cp(*bond);
    CHECK(cp.getBondType() == bond->getBondType());
    CHECK(!cp.hasOwningMol());
  }
}

TEST_CASE("valence edge") {
  {  // this is, of course, absurd:
    auto m = "[H-2]"_smiles;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setNoImplicit(false);
    m->updatePropertyCache(false);
    CHECK(m->getAtomWithIdx(0)->getFormalCharge() == -2);
    CHECK(m->getAtomWithIdx(0)->getImplicitValence() == 0);
  }
  {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m{SmilesToMol("CFC", ps)};
    REQUIRE(m);
    CHECK_THROWS_AS(m->getAtomWithIdx(1)->calcImplicitValence(true),
                    AtomValenceException);
  }
}

TEST_CASE("SetQuery on normal atoms") {
  auto m = "CC"_smiles;
  REQUIRE(m);
  auto qry = makeAtomAliphaticQuery();
  CHECK_THROWS_AS(m->getAtomWithIdx(0)->setQuery(qry), std::runtime_error);
  CHECK_THROWS_AS(m->getAtomWithIdx(0)->expandQuery(qry), std::runtime_error);
  delete qry;
}

TEST_CASE("SetQuery on normal bonds") {
  auto m = "CC"_smiles;
  REQUIRE(m);
  auto qry = makeBondOrderEqualsQuery(Bond::BondType::SINGLE);
  CHECK_THROWS_AS(m->getBondWithIdx(0)->setQuery(qry), std::runtime_error);
  CHECK_THROWS_AS(m->getBondWithIdx(0)->expandQuery(qry), std::runtime_error);
  delete qry;
}

TEST_CASE("additional atom props") {
  auto m = "CC"_smiles;
  REQUIRE(m);
  auto atom = m->getAtomWithIdx(0);
  {
    CHECK(!atom->hasProp(common_properties::_MolFileRLabel));
    setAtomRLabel(atom, 1);
    CHECK(atom->hasProp(common_properties::_MolFileRLabel));
    setAtomRLabel(atom, 0);
    CHECK(!atom->hasProp(common_properties::_MolFileRLabel));
  }
  {
    CHECK(!atom->hasProp(common_properties::molFileAlias));
    setAtomAlias(atom, "foo");
    CHECK(atom->hasProp(common_properties::molFileAlias));
    setAtomAlias(atom, "");
    CHECK(!atom->hasProp(common_properties::molFileAlias));
  }
  {
    CHECK(!atom->hasProp(common_properties::molFileValue));
    setAtomValue(atom, "foo");
    CHECK(atom->hasProp(common_properties::molFileValue));
    setAtomValue(atom, "");
    CHECK(!atom->hasProp(common_properties::molFileValue));
  }
  {
    CHECK(!atom->hasProp(common_properties::_supplementalSmilesLabel));
    setSupplementalSmilesLabel(atom, "foo");
    CHECK(atom->hasProp(common_properties::_supplementalSmilesLabel));
    setSupplementalSmilesLabel(atom, "");
    CHECK(!atom->hasProp(common_properties::_supplementalSmilesLabel));
  }
}

TEST_CASE("getBondTypeAsDouble()") {
  SECTION("plain") {
    std::vector<std::pair<Bond::BondType, double>> vals{
        {Bond::BondType::IONIC, 0},
        {Bond::BondType::ZERO, 0},
        {Bond::BondType::SINGLE, 1},
        {Bond::BondType::DOUBLE, 2},
        {Bond::BondType::TRIPLE, 3},
        {Bond::BondType::QUADRUPLE, 4},
        {Bond::BondType::QUINTUPLE, 5},
        {Bond::BondType::HEXTUPLE, 6},
        {Bond::BondType::ONEANDAHALF, 1.5},
        {Bond::BondType::TWOANDAHALF, 2.5},
        {Bond::BondType::THREEANDAHALF, 3.5},
        {Bond::BondType::FOURANDAHALF, 4.5},
        {Bond::BondType::FIVEANDAHALF, 5.5},
        {Bond::BondType::AROMATIC, 1.5},
        {Bond::BondType::DATIVEONE, 1.0},
        {Bond::BondType::DATIVE, 1.0},
        {Bond::BondType::HYDROGEN, 0}

    };
    for (const auto &pr : vals) {
      Bond bnd(pr.first);
      CHECK(bnd.getBondType() == pr.first);
      CHECK(bnd.getBondTypeAsDouble() == pr.second);
    }
  }
  SECTION("twice") {
    std::vector<std::pair<Bond::BondType, std::uint8_t>> vals{
        {Bond::BondType::IONIC, 0},         {Bond::BondType::ZERO, 0},
        {Bond::BondType::SINGLE, 2},        {Bond::BondType::DOUBLE, 4},
        {Bond::BondType::TRIPLE, 6},        {Bond::BondType::QUADRUPLE, 8},
        {Bond::BondType::QUINTUPLE, 10},    {Bond::BondType::HEXTUPLE, 12},
        {Bond::BondType::ONEANDAHALF, 3},   {Bond::BondType::TWOANDAHALF, 5},
        {Bond::BondType::THREEANDAHALF, 7}, {Bond::BondType::FOURANDAHALF, 9},
        {Bond::BondType::FIVEANDAHALF, 11}, {Bond::BondType::AROMATIC, 3},
        {Bond::BondType::DATIVEONE, 2},     {Bond::BondType::DATIVE, 2},
        {Bond::BondType::HYDROGEN, 0}

    };
    for (const auto &pr : vals) {
      Bond bnd(pr.first);
      CHECK(bnd.getBondType() == pr.first);
      CHECK(getTwiceBondType(bnd) == pr.second);
    }
  }
}

TEST_CASE("getValenceContrib()") {
  const auto m = "CO->[Fe]"_smiles;
  REQUIRE(m);
  CHECK(m->getBondWithIdx(1)->getValenceContrib(m->getAtomWithIdx(0)) == 0);
  CHECK(m->getBondWithIdx(1)->getValenceContrib(m->getAtomWithIdx(1)) == 0);
  CHECK(m->getBondWithIdx(1)->getValenceContrib(m->getAtomWithIdx(2)) == 1);
}

TEST_CASE("conformer details") {
  const auto m = "CC"_smiles;
  REQUIRE(m);
  Conformer *conf = new Conformer(m->getNumAtoms());
  CHECK(!conf->hasOwningMol());
  m->addConformer(conf);
  CHECK(conf->hasOwningMol());
  auto cid = conf->getId();
  *conf = *conf;
  CHECK(conf->hasOwningMol());
  CHECK(conf->getId() == cid);
}

#if !defined(_WIN32) || !defined(RDKIT_DYN_LINK)
namespace RDKit {
namespace Canon {
namespace details {
bool atomHasFourthValence(const Atom *atom);
bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q);
}  // namespace details
void switchBondDir(Bond *bond);
}  // namespace Canon
}  // namespace RDKit
TEST_CASE("canon details") {
  SECTION("h queries") {
    std::vector<std::pair<std::string, bool>> examples{
        {"C[CHD3](F)Cl", true},  {"C[CD3H](F)Cl", true},
        {"C[CH3D](F)Cl", false}, {"C[CDH3](F)Cl", false},
        {"C[CDR4H](F)Cl", true},
    };
    for (const auto &pr : examples) {
      std::unique_ptr<RWMol> m{SmartsToMol(pr.first)};
      REQUIRE(m);
      CHECK(RDKit::Canon::details::hasSingleHQuery(
                m->getAtomWithIdx(1)->getQuery()) == pr.second);
      CHECK(RDKit::Canon::details::atomHasFourthValence(m->getAtomWithIdx(1)) ==
            pr.second);
      // artificial, but causes atomHasFourthValence to always return true
      m->getAtomWithIdx(1)->setNumExplicitHs(1);
      CHECK(RDKit::Canon::details::atomHasFourthValence(m->getAtomWithIdx(1)));
    }
  }
}
TEST_CASE("switchBondDir") {
  auto m = "C/C=C/C"_smiles;
  REQUIRE(m);
  auto bond = m->getBondWithIdx(0);
  CHECK(bond->getBondDir() == Bond::BondDir::ENDUPRIGHT);
  Canon::switchBondDir(bond);
  CHECK(bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT);
  bond->setBondDir(Bond::BondDir::UNKNOWN);
  Canon::switchBondDir(bond);
  CHECK(bond->getBondDir() == Bond::BondDir::UNKNOWN);
}
#endif

TEST_CASE("allow 5 valent N/P/As to kekulize", "[kekulization]") {
  std::vector<std::pair<std::string, std::string>> tests = {
      {"O=n1ccccc1", "O=N1=CC=CC=C1"},
      {"O=p1ccccc1", "O=P1=CC=CC=C1"},
      {"O=[as]1ccccc1", "O=[As]1=CC=CC=C1"}};
  SmilesParserParams ps;
  ps.sanitize = false;
  SECTION("kekulization") {
    for (const auto &pr : tests) {
      INFO(pr.first);
      std::unique_ptr<RWMol> m{SmilesToMol(pr.first, ps)};
      REQUIRE(m);
      m->updatePropertyCache(false);
      MolOps::Kekulize(*m);
      CHECK(MolToSmiles(*m) == pr.second);
    }
  }
  SECTION("sanitization") {
    for (const auto &pr : tests) {
      INFO(pr.first);
      std::unique_ptr<RWMol> m{SmilesToMol(pr.first, ps)};
      REQUIRE(m);
      m->updatePropertyCache(false);
      unsigned int failed;
      unsigned int flags = MolOps::SanitizeFlags::SANITIZE_ALL ^
                           MolOps::SanitizeFlags::SANITIZE_CLEANUP ^
                           MolOps::SanitizeFlags::SANITIZE_PROPERTIES;
      MolOps::sanitizeMol(*m, failed, flags);
      CHECK(!failed);
      CHECK(MolToSmiles(*m) == pr.second);
    }
  }
}

TEST_CASE("KekulizeIfPossible") {
  SECTION("basics: molecules with failures") {
    std::vector<std::string> smis = {
        "c1cccn1",
        "c1ccccc1-c1cccn1",
    };
    for (const auto &smi : smis) {
      SmilesParserParams ps;
      ps.sanitize = false;
      std::unique_ptr<RWMol> m{SmilesToMol(smi, ps)};
      REQUIRE(m);
      m->updatePropertyCache(false);
      // confirm that we normally fail:
      {
        RWMol m2(*m);
        CHECK_THROWS_AS(MolOps::Kekulize(m2), MolSanitizeException);
      }
      {
        RWMol m2(*m);
        CHECK(!MolOps::KekulizeIfPossible(m2));
        for (unsigned i = 0; i < m2.getNumAtoms(); ++i) {
          CHECK(m2.getAtomWithIdx(i)->getIsAromatic() ==
                m->getAtomWithIdx(i)->getIsAromatic());
        }
        for (unsigned i = 0; i < m2.getNumBonds(); ++i) {
          CHECK(m2.getBondWithIdx(i)->getIsAromatic() ==
                m->getBondWithIdx(i)->getIsAromatic());
          CHECK(m2.getBondWithIdx(i)->getBondType() ==
                m->getBondWithIdx(i)->getBondType());
        }
      }
    }
  }
  SECTION("basics: molecules without failures") {
    std::vector<std::string> smis = {"c1ccc[nH]1", "c1ccccc1"};
    for (const auto &smi : smis) {
      SmilesParserParams ps;
      ps.sanitize = false;
      std::unique_ptr<RWMol> m{SmilesToMol(smi, ps)};
      REQUIRE(m);
      m->updatePropertyCache(false);
      {
        RWMol m2(*m);
        MolOps::Kekulize(m2);
        RWMol m3(*m);
        CHECK(MolOps::KekulizeIfPossible(m3));
        for (unsigned i = 0; i < m2.getNumAtoms(); ++i) {
          CHECK(m2.getAtomWithIdx(i)->getIsAromatic() ==
                m3.getAtomWithIdx(i)->getIsAromatic());
        }
        for (unsigned i = 0; i < m2.getNumBonds(); ++i) {
          CHECK(m2.getBondWithIdx(i)->getIsAromatic() ==
                m3.getBondWithIdx(i)->getIsAromatic());
          CHECK(m2.getBondWithIdx(i)->getBondType() ==
                m3.getBondWithIdx(i)->getBondType());
        }
      }
    }
  }
}

TEST_CASE("Github #4535: operator<< for AtomPDBResidue", "[PDB]") {
  SECTION("basics") {
    bool sanitize = true;
    int flavor = 0;
    std::unique_ptr<RWMol> mol(SequenceToMol("KY", sanitize, flavor));
    REQUIRE(mol);
    REQUIRE(mol->getAtomWithIdx(0)->getMonomerInfo());
    auto res = static_cast<AtomPDBResidueInfo *>(
        mol->getAtomWithIdx(0)->getMonomerInfo());
    REQUIRE(res);
    std::stringstream oss;
    oss << *res << std::endl;
    res = static_cast<AtomPDBResidueInfo *>(
        mol->getAtomWithIdx(mol->getNumAtoms() - 1)->getMonomerInfo());
    REQUIRE(res);
    oss << *res << std::endl;
    auto tgt = R"FOO(1  N   LYS A 1
22  OXT TYR A 2
)FOO";
    CHECK(oss.str() == tgt);
  }
}

TEST_CASE("isAromaticAtom") {
  SECTION("basics") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> mol(SmilesToMol("C:C:C", ps));
    REQUIRE(mol);
    CHECK(!mol->getAtomWithIdx(0)->getIsAromatic());
    CHECK(mol->getBondWithIdx(0)->getIsAromatic());
    CHECK(isAromaticAtom(*mol->getAtomWithIdx(0)));
  }
}

TEST_CASE(
    "Github #4785: aromatic bonds no longer set aromatic flags on atoms") {
  SECTION("basics1") {
    auto m = "C1:C:C:C:1"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "C1=CC=C1");
  }
  SECTION("basics2") {
    auto m = "C1:C:C:C:C:C:1"_smiles;
    REQUIRE(m);
    CHECK(MolToSmiles(*m) == "c1ccccc1");
  }
  SECTION("can still get kekulization errors") {
    CHECK_THROWS_AS(SmilesToMol("C1:C:C:C:C:1"), KekulizeException);
  }
}

namespace {
void check_dest(RWMol *m1, const ROMol &m2) {
  CHECK(m2.getNumAtoms() == 8);
  CHECK(m2.getNumBonds() == 7);
  for (const auto atom : m2.atoms()) {
    CHECK(&atom->getOwningMol() == &m2);
    CHECK(&atom->getOwningMol() != m1);
  }
  for (const auto bond : m2.bonds()) {
    CHECK(&bond->getOwningMol() == &m2);
    CHECK(&bond->getOwningMol() != m1);
  }
  CHECK(m2.getStereoGroups().size() == 2);
  CHECK(m2.getStereoGroups()[0].getAtoms().size() == 2);
  CHECK(m2.getStereoGroups()[0].getAtoms()[0]->getIdx() == 1);
  CHECK(m2.getStereoGroups()[0].getAtoms()[1]->getIdx() == 5);
  CHECK(m2.getStereoGroups()[1].getAtoms().size() == 1);
  CHECK(m2.getStereoGroups()[1].getAtoms()[0]->getIdx() == 3);

  const auto &sgs = getSubstanceGroups(m2);
  CHECK(sgs.size() == 1);
  CHECK(sgs[0].getAtoms().size() == 1);
  CHECK(sgs[0].getAtoms()[0] == 4);

  // check the state of m1:
  CHECK(m1->getNumAtoms() == 0);
  CHECK(m1->getNumBonds() == 0);
  CHECK(m1->getPropList().empty());
  CHECK(m1->getDict().getData().empty());
  CHECK(m1->getStereoGroups().empty());
  CHECK(getSubstanceGroups(*m1).empty());
  CHECK(m1->getRingInfo() == nullptr);

  // make sure we can still do something with m1:
  *m1 = m2;
  CHECK(!m1->getDict().getData().empty());
  CHECK(m1->getNumAtoms() == 8);
  CHECK(m1->getNumBonds() == 7);
  CHECK(m1->getRingInfo() != nullptr);
  CHECK(m1->getRingInfo()->isInitialized() ==
        m2.getRingInfo()->isInitialized());
}
}  // namespace
TEST_CASE("moves") {
  auto m1 =
      "C[C@H](O)[C@H](F)[C@@H](C)O |o2:1,5,&1:3,SgD:4:atom_data:foo::::|"_smiles;
  REQUIRE(m1);
  CHECK(m1->getStereoGroups().size() == 2);
  m1->setProp("foo", 1u);
  SECTION("molecule move") {
    ROMol m2 = std::move(*m1);
    check_dest(m1.get(), m2);
  }
  SECTION("molecule move-assign") {
    ROMol m2;
    m2 = std::move(*m1);
    check_dest(m1.get(), m2);
  }
}

TEST_CASE("query moves") {
  auto m1 =
      "C[C@H](O)[C@H](F)[C@@H](C)O |o2:1,5,&1:3,SgD:4:atom_data:foo::::|"_smarts;
  REQUIRE(m1);
  CHECK(m1->getStereoGroups().size() == 2);
  m1->setProp("foo", 1u);
  SECTION("molecule move") {
    ROMol m2 = std::move(*m1);
    check_dest(m1.get(), m2);
  }
  SECTION("molecule move-assign") {
    ROMol m2;
    m2 = std::move(*m1);
    check_dest(m1.get(), m2);
  }
}

TEST_CASE("moves with conformer") {
  auto m1 = R"CTAB(
  Mrv2108 01192209042D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.31 4.001 0 0
M  V30 2 C 1.54 2.6674 0 0 CFG=2
M  V30 3 O -0 2.6674 0 0
M  V30 4 C 2.31 1.3337 0 0 CFG=1
M  V30 5 F 3.85 1.3337 0 0
M  V30 6 C 1.54 0 0 0 CFG=2
M  V30 7 C 2.31 -1.3337 0 0
M  V30 8 O 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=3
M  V30 3 1 4 2
M  V30 4 1 4 5 CFG=1
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 7 1 6 8 CFG=3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 5) FIELDNAME=atom_data -
M  V30 FIELDDISP="    4.6200    0.5637    DA    ALL  0       0" FIELDDATA=foo
M  V30 END SGROUP
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEREL2 ATOMS=(2 2 6)
M  V30 MDLV30/STERAC1 ATOMS=(1 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
  REQUIRE(m1);
  CHECK(m1->getStereoGroups().size() == 2);
  m1->setProp("foo", 1u);
  SECTION("molecule move") {
    ROMol m2 = std::move(*m1);
    check_dest(m1.get(), m2);
  }
  SECTION("molecule move-assign") {
    ROMol m2;
    m2 = std::move(*m1);
    check_dest(m1.get(), m2);
  }
}

TEST_CASE("Github #5055") {
  SECTION("as reported") {
    auto m =
        "CC1(C)NC(=O)CN2C=C(C[C@H](C(=O)NC)NC(=O)CN3CCN(C(=O)[C@H]4Cc5c([nH]c6ccccc56)CN4C(=O)CN4CN(c5ccccc5)C5(CCN(CC5)C1=O)C4=O)[C@@H](Cc1ccccc1)C3=O)[N-][NH2+]2"_smiles;
    REQUIRE(m);
  }
}
TEST_CASE("Iterators") {
  auto m = "CCCCCCC=N"_smiles;
  REQUIRE(m);
  SECTION("Atom Iterator") {
    auto atoms = m->atoms();
    auto n_atom = std::find_if(atoms.begin(), atoms.end(), [](const auto &a) {
      return a->getAtomicNum() == 7;
    });
    REQUIRE(n_atom != atoms.end());
    CHECK((*n_atom)->getIdx() == 7);
  }
  SECTION("Atom Neighbor Iterator") {
    auto nbrs = m->atomNeighbors(m->getAtomWithIdx(6));
    auto n_atom = std::find_if(nbrs.begin(), nbrs.end(), [](const auto &a) {
      return a->getAtomicNum() == 7;
    });
    REQUIRE(n_atom != nbrs.end());
    CHECK((*n_atom)->getIdx() == 7);
  }
  SECTION("Bond Iterator") {
    auto bonds = m->bonds();
    auto dbl_bnd = std::find_if(bonds.begin(), bonds.end(), [](const auto &b) {
      return b->getBondType() == Bond::DOUBLE;
    });
    REQUIRE(dbl_bnd != bonds.end());
    CHECK((*dbl_bnd)->getIdx() == 6);
  }
  SECTION("Atom Bond Iterator") {
    auto nbr_bonds = m->atomBonds(m->getAtomWithIdx(6));
    auto dbl_bnd = std::find_if(
        nbr_bonds.begin(), nbr_bonds.end(),
        [](const auto &b) { return b->getBondType() == Bond::DOUBLE; });
    REQUIRE(dbl_bnd != nbr_bonds.end());
    CHECK((*dbl_bnd)->getIdx() == 6);
  }
}

TEST_CASE("general hybridization") {
  auto m = "CS(=O)(=O)C"_smiles;
  REQUIRE(m);
  CHECK(m->getAtomWithIdx(0)->getHybridization() ==
        Atom::HybridizationType::SP3);
  CHECK(m->getAtomWithIdx(1)->getHybridization() ==
        Atom::HybridizationType::SP3);
  CHECK(m->getAtomWithIdx(4)->getHybridization() ==
        Atom::HybridizationType::SP3);
}

TEST_CASE("metal hybridization") {
  SECTION("square planar") {
    {
      auto m = "F[U@SP](F)(F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP2D);
    }
    {
      auto m = "F[U@SP](F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP2D);
    }
    {
      auto m = "F[U@SP](F)(F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP2D);
    }
    {
      auto m = "F[U@TH](F)(F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP3);
    }
  }
  SECTION("octahedral") {
    {
      auto m = "F[S@OH](F)(F)(F)(F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP3D2);
    }
    {
      auto m = "F[P@OH](F)(F)(F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP3D2);
    }
    {
      auto m = "F[P](F)(F)(F)F"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getHybridization() ==
            Atom::HybridizationType::SP3D);
    }
  }
}

TEST_CASE(
    "github #5462: Invalid number of radical electrons calculated for [Pr+4]") {
  SECTION("as reported") {
    auto m = "[Pr+4]"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  }
  SECTION("also reported") {
    auto m = "[U+5]"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  }
  SECTION("extreme") {
    auto m = "[Fe+30]"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 0);
  }
}

TEST_CASE(
    "github #5505: Running kekulization on mols with query bonds will either fail or return incorrect results") {
  SECTION("as reported") {
    auto m = "[#6]-c1cccc(-[#6])c1"_smarts;
    REQUIRE(m);
    REQUIRE(m->getBondWithIdx(2)->hasQuery());
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::BondType::AROMATIC);
    MolOps::Kekulize(*m);
    REQUIRE(m->getBondWithIdx(2)->hasQuery());
    REQUIRE(m->getBondWithIdx(2)->getBondType() == Bond::BondType::AROMATIC);
  }
  SECTION("partial kekulization works") {
    auto m1 = "c1ccccc1"_smiles;
    REQUIRE(m1);
    auto m2 = "c1ccccc1"_smarts;
    REQUIRE(m2);
    m1->insertMol(*m2);
    MolOps::findSSSR(*m1);
    REQUIRE(!m1->getBondWithIdx(1)->hasQuery());
    REQUIRE(m1->getBondWithIdx(1)->getBondType() == Bond::BondType::AROMATIC);
    REQUIRE(m1->getBondWithIdx(6)->hasQuery());
    REQUIRE(m1->getBondWithIdx(6)->getBondType() == Bond::BondType::AROMATIC);
    MolOps::Kekulize(*m1);
    REQUIRE(!m1->getBondWithIdx(1)->hasQuery());
    REQUIRE(m1->getBondWithIdx(1)->getBondType() != Bond::BondType::AROMATIC);
    REQUIRE(m1->getBondWithIdx(6)->hasQuery());
    REQUIRE(m1->getBondWithIdx(6)->getBondType() == Bond::BondType::AROMATIC);
  }
  SECTION("kekulization with non-bond-type queries works") {
    auto m1 = "c1ccccc1"_smiles;
    REQUIRE(m1);
    QueryBond qbond;
    qbond.setBondType(Bond::BondType::AROMATIC);
    qbond.setIsAromatic(true);
    qbond.setQuery(makeBondIsInRingQuery());
    m1->replaceBond(0, &qbond);
    REQUIRE(m1->getBondWithIdx(0)->hasQuery());
    REQUIRE(m1->getBondWithIdx(0)->getBondType() == Bond::BondType::AROMATIC);
    MolOps::Kekulize(*m1);
    REQUIRE(m1->getBondWithIdx(0)->hasQuery());
    REQUIRE(m1->getBondWithIdx(0)->getBondType() != Bond::BondType::AROMATIC);
  }
  SECTION("make sure single-atom molecules and mols without rings still fail") {
    std::vector<std::string> expectedFailures = {"p", "c:c"};
    for (const auto &smi : expectedFailures) {
      INFO(smi);
      CHECK_THROWS_AS(SmilesToMol(smi), MolSanitizeException);
    }
  }
}

TEST_CASE("extended valences for alkali earths") {
  SECTION("valence of 2") {
    // make sure the valence of two works
    std::vector<std::pair<std::string, unsigned int>> cases{
        {"C[Be]", 1},  {"C[Mg]", 1},  {"C[Ca]", 1},  {"C[Sr]", 1},
        {"C[Ba]", 1},  {"C[Ra]", 1},  {"C[Be]C", 0}, {"C[Mg]C", 0},
        {"C[Ca]C", 0}, {"C[Sr]C", 0}, {"C[Ba]C", 0}, {"C[Ra]C", 0}};

    for (const auto &pr : cases) {
      INFO(pr.first);
      std::unique_ptr<RWMol> m{SmilesToMol(pr.first)};
      REQUIRE(m);
      m->getAtomWithIdx(1)->setNoImplicit(false);
      m->getAtomWithIdx(1)->setNumRadicalElectrons(0);
      MolOps::sanitizeMol(*m);
      CHECK(m->getAtomWithIdx(1)->getTotalNumHs() == pr.second);
    }
  }
  SECTION("higher valence") {
    // make sure the valence of two works
    std::vector<std::pair<std::string, unsigned int>> cases{{"C[Mg](C)C", 0},
                                                            {"C[Ca](C)C", 0},
                                                            {"C[Sr](C)C", 0},
                                                            {"C[Ba](C)C", 0},
                                                            {"C[Ra](C)C", 0}};

    for (const auto &pr : cases) {
      INFO(pr.first);
      std::unique_ptr<RWMol> m{SmilesToMol(pr.first)};
      REQUIRE(m);
      m->getAtomWithIdx(1)->setNoImplicit(false);
      m->getAtomWithIdx(1)->setNumRadicalElectrons(0);
      MolOps::sanitizeMol(*m);
      CHECK(m->getAtomWithIdx(1)->getTotalNumHs() == pr.second);
    }
  }
  SECTION("everybody loves grignards") {
    auto m = "CC(C)O[Mg](Cl)(<-O1CCCC1)<-O1CCCC1"_smiles;
    REQUIRE(m);
  }
}

TEST_CASE("Github #5849: aromatic tag allows bad valences to pass") {
  SECTION("basics") {
    std::vector<std::string> smis = {
        "Cc1(C)=NCCCC1",  // one bogus aromatic atom
    };
    for (const auto &smi : smis) {
      INFO(smi);
      CHECK_THROWS_AS(SmilesToMol(smi), AtomValenceException);
    }
  }
  SECTION("as reported") {
    std::vector<std::string> smis = {
        "c1c(ccc2NC(CN=c(c21)(C)C)=O)O",
        "c1c(c2n(C(NC(Oc2cc1)C)c1c[nH]cn1)=O)O",
        "c1c2C(c4c(cc(c(C(=O)O)c4)O)O)Oc(c2c(c(O)c1O)O)(=O)O",
        "c12C(CN(C)C)CCCCN=c(c2cc2c1cc[nH]c2)(C)C"};
    for (const auto &smi : smis) {
      INFO(smi);
      CHECK_THROWS_AS(SmilesToMol(smi), MolSanitizeException);
    }
  }
  SECTION("edge cases atoms") {
    std::vector<std::string> smis = {
        "c1c(ccc2NC(CN=c(c21)=C)=O)O",     // exocyclic double bond
        "c1c(ccc2NC(Cn=c(c21)(C)C)=O)O",   // even more bogus aromatic atoms
        "c1c(ccc2NC(CN=c(c21)(:c)C)=O)O",  // even more bogus aromatic atoms
    };
    for (const auto &smi : smis) {
      INFO(smi);
      CHECK_THROWS_AS(SmilesToMol(smi), AtomValenceException);
    }
  }
  SECTION("edge cases kekulization") {
    std::vector<std::string> smis = {
        "c12ccccc1=NCCC2",  // adjacent to an actual aromatic ring
        "Cc1(C)=NCCCc1",    // two bogus aromatic atoms
        "Cc1(C)nCCCC=1",    // two bogus aromatic atoms
        "Cc1(C)nCCCc1",     // three bogus aromatic atoms
        "Cc:1(C):nCCCc:1",  // three bogus aromatic atoms
    };
    for (const auto &smi : smis) {
      INFO(smi);
      CHECK_THROWS_AS(SmilesToMol(smi), KekulizeException);
    }
  }
}

TEST_CASE("molecules with more than 255 rings produce a bad pickle") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/";
  bool sanitize = false;
  std::unique_ptr<RWMol> mol(
      MolFileToMol(pathName + "mol_with_pickle_error.mol", sanitize));
  REQUIRE(mol);
  mol->updatePropertyCache(false);
  unsigned int opThatFailed = 0;
  MolOps::sanitizeMol(*mol, opThatFailed,
                      MolOps::SanitizeFlags::SANITIZE_ALL ^
                          MolOps::SanitizeFlags::SANITIZE_PROPERTIES ^
                          MolOps::SANITIZE_KEKULIZE);
  CHECK(mol->getRingInfo()->numRings() > 300);
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl);
  RWMol nMol(pkl);
  CHECK(nMol.getRingInfo()->numRings() == mol->getRingInfo()->numRings());
}