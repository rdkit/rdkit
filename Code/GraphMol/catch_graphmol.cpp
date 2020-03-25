//
//
//  Copyright (C) 2018-2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
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

TEST_CASE("Github #2062", "[bug, molops]") {
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

TEST_CASE("Github #2086", "[bug, molops]") {
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

TEST_CASE("github #299", "[bug, molops, SSSR]") {
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

TEST_CASE("github #2224", "[bug, molops, removeHs, query]") {
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
    "[bug,stereo]") {
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

TEST_CASE("github #2244", "[bug, molops, stereo]") {
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
    "[bug, molops]") {
  SECTION("the original report") {
    std::vector<std::string> smiles = {"C=n1ccnc1", "C#n1ccnc1"};
    for (auto smi : smiles) {
      CHECK_THROWS_AS(SmilesToMol(smi), MolSanitizeException);
    }
  }
}

TEST_CASE("github #908: AddHs() using 3D coordinates with 2D conformations",
          "[bug, molops]") {
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
    "[bug, molops]") {
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
    CHECK(!mol->getRingInfo()->isInitialized());
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
    "[bug, molops]") {
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
      } catch (const AtomValenceException &e) {
        CHECK(e.getType() == "AtomValenceException");
        CHECK(e.getAtomIdx() == pr.second);
      }
    }
  }
  SECTION("AtomKekulizeException") {
    std::vector<std::pair<std::string, unsigned int>> smiles = {
        {"CCcc", 2}, {"C1:c:CC1", 0}};
    for (auto pr : smiles) {
      CHECK_THROWS_AS(SmilesToMol(pr.first), AtomKekulizeException);
      try {
        auto m = SmilesToMol(pr.first);
      } catch (const AtomKekulizeException &e) {
        CHECK(e.getType() == "AtomKekulizeException");
        CHECK(e.getAtomIdx() == pr.second);
      }
    }
  }
  SECTION("KekulizeException") {
    std::vector<std::pair<std::string, std::vector<unsigned int>>> smiles = {
        {"c1cccc1", {0, 1, 2, 3, 4}}, {"Cc1cc1", {1, 2, 3}}};
    for (auto pr : smiles) {
      CHECK_THROWS_AS(SmilesToMol(pr.first), KekulizeException);
      try {
        auto m = SmilesToMol(pr.first);
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
    "[bug, molops]") {
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
    "[bug, molops]") {
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
    "[bug, stereochemistry]") {
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

TEST_CASE("removeHs screwing up double bond stereo", "[bug,removeHs]") {
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

TEST_CASE("setDoubleBondNeighborDirections()", "[stereochemistry,bug]") {
  SECTION("basics") {
    auto m = "CC=CC"_smiles;
    REQUIRE(m);
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
    "[transuranics,bug]") {
  auto pt = PeriodicTable::getTable();
  SECTION("number to symbol") {
    std::vector<std::pair<unsigned int, std::string>> data = {
        {113, "Nh"}, {114, "Fl"}, {115, "Mc"},
        {116, "Lv"}, {117, "Ts"}, {118, "Og"}};
    for (auto pr : data) {
      CHECK(pt->getElementSymbol(pr.first) == pr.second);
    }
  }
  SECTION("symbol to number") {
    std::vector<std::pair<unsigned int, std::string>> data = {
        {113, "Nh"}, {114, "Fl"}, {115, "Mc"},  {116, "Lv"},
        {117, "Ts"}, {118, "Og"}, {113, "Uut"}, {115, "Uup"}};
    for (auto pr : data) {
      CHECK(pt->getAtomicNumber(pr.second) == pr.first);
    }
  }
}
TEST_CASE("github #2775", "[valence,bug]") {
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
          "[molops,bug,aromaticity]") {
  SECTION("acepentalene") {
    std::unique_ptr<RWMol> m{SmilesToMol("C1=CC2=CC=C3C2=C1C=C3")};
    REQUIRE(m);
    auto smi = MolToSmiles(*m);
    CHECK(smi == "C1=CC2=C3C1=CC=C3C=C2");
  }
}

TEST_CASE("handling of bondStereoCare in updateQueryProperties") {
  SECTION("fully specified") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
  }
  SECTION("fully unspecified") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREONONE);
  }
  SECTION("partially unspecified") {
    std::vector<std::string> mbs = {R"CTAB(keep
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(keep
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(remove
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=0
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(remove
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(remove
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"};
    for (const auto &mb : mbs) {
      std::unique_ptr<RWMol> mol{MolBlockToMol(mb)};
      REQUIRE(mol);
      REQUIRE(mol->getBondBetweenAtoms(0, 1));
      CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
            Bond::BondStereo::STEREOE);
      MolOps::AdjustQueryParameters ps;
      ps.useStereoCareForBonds = true;
      MolOps::adjustQueryProperties(*mol, &ps);
      if (mol->getProp<std::string>(common_properties::_Name) == "keep") {
        CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
              Bond::BondStereo::STEREOE);
      } else {
        CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
              Bond::BondStereo::STEREONONE);
      }
    }
  }
  SECTION("V2000") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292015042D          

  4  3  0  0  0  0            999 V2000
   -3.7669    1.1053    0.0000 C   0  0  0  0  1  0  0  0  0  0  0  0
   -3.0524    1.5178    0.0000 C   0  0  0  0  1  0  0  0  0  0  0  0
   -2.3380    1.1053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4814    1.5178    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  2  2  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->hasProp(common_properties::molStereoCare));
    CHECK(mol->getAtomWithIdx(1)->hasProp(common_properties::molStereoCare));
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    // property added by the CTAB parser:
    CHECK(mol->getBondBetweenAtoms(0, 1)->hasProp(
        common_properties::molStereoCare));
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
  }
  SECTION("molecule from SMILES") {
    auto mol = "C/C=C/C"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(2, 1));
    CHECK(mol->getBondBetweenAtoms(2, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    // since stereoCare is not set on the bond from SMILES,
    // stereochem will be removed:
    {
      RWMol molcp(*mol);
      MolOps::adjustQueryProperties(molcp, &ps);
      CHECK(molcp.getBondBetweenAtoms(2, 1)->getStereo() ==
            Bond::BondStereo::STEREONONE);
    }
    // but we can preserve it by setting the property:
    {
      RWMol molcp(*mol);
      molcp.getBondBetweenAtoms(2, 1)->setProp(common_properties::molStereoCare,
                                               1);
      MolOps::adjustQueryProperties(molcp, &ps);
      CHECK(molcp.getBondBetweenAtoms(2, 1)->getStereo() ==
            Bond::BondStereo::STEREOE);
    }
  }
}

TEST_CASE("updateQueryParameters from JSON") {
  SECTION("basics") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.makeAtomsGeneric == false);
    CHECK(ps.makeBondsGeneric == false);
    CHECK(ps.makeBondsGenericFlags == MolOps::ADJUST_IGNORENONE);

    std::string json = R"JSON({"makeAtomsGeneric":true})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);

    CHECK(ps.makeAtomsGeneric == true);
    CHECK(ps.makeBondsGeneric == false);
    // the parsing updates the parameters, it doesn't replace them:

    json = R"JSON({"makeBondsGeneric":true,
      "makeBondsGenericFlags":"IGNOREDUMMIES|IGNORECHAINS"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);

    CHECK(ps.makeAtomsGeneric == true);
    CHECK(ps.makeBondsGeneric == true);
    CHECK(ps.makeBondsGenericFlags ==
          (MolOps::ADJUST_IGNOREDUMMIES | MolOps::ADJUST_IGNORECHAINS));
  }
  SECTION("useStereoCare") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.useStereoCareForBonds == false);

    std::string json = R"JSON({"useStereoCareForBonds":true})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.useStereoCareForBonds == true);
    json = R"JSON({"useStereoCareForBonds":false})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.useStereoCareForBonds == false);
  }
  SECTION("bogus contents") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.adjustDegree == true);
    CHECK(ps.adjustDegreeFlags ==
          (MolOps::ADJUST_IGNOREDUMMIES | MolOps::ADJUST_IGNORECHAINS));

    std::string json = R"JSON({"bogosity":true})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustDegree == true);

    json = R"JSON({"adjustDegree":"foo"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustDegree == true);

    json = R"JSON({"adjustDegreeFlags":"IGNORENONE|bogus"})JSON";
    // clang-format off
    CHECK_THROWS_AS(MolOps::parseAdjustQueryParametersFromJSON(ps, json),ValueErrorException);
  }
}

TEST_CASE("phosphine and arsine chirality", "[Chirality]") {
  SECTION("chiral center recognized"){
    auto mol1 = "C[P@](C1CCCC1)C1=CC=CC=C1"_smiles;
    auto mol2 = "C[As@](C1CCCC1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(mol2->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center selective"){
    auto mol1 = "C[P@](C)C1CCCCC1"_smiles;
    auto mol2 = "C[As@](C)C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(mol2->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center specific: P"){
    auto mol1 = "C[P@](C1CCCC1)C1=CC=CC=C1"_smiles;
    auto mol2 = "C[P@@](C1CCCC1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol1); 
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
  SECTION("chiral center specific: As"){
    auto mol1 = "C[As@](C1CCCC1)C1=CC=CC=C1"_smiles;
    auto mol2 = "C[As@@](C1CCCC1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol1); 
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
  SECTION("chiral center, implicit H: P"){
    auto mol1 = "C[P@H]C1CCCCC1"_smiles;
    auto mol2 = "C[P@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center, implicit H: As"){
    auto mol1 = "C[As@H]C1CCCCC1"_smiles;
    auto mol2 = "C[As@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("chiral center specific, implicit H: P"){
    auto mol1 = "C[P@H]C1CCCCC1"_smiles;
    auto mol2 = "C[P@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
  SECTION("chiral center specific, implicit H: As"){
    auto mol1 = "C[As@H]C1CCCCC1"_smiles;
    auto mol2 = "C[As@@H]C1CCCCC1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    CHECK(MolToSmiles(*mol1) != MolToSmiles(*mol2));
  }
}

TEST_CASE("github #2890", "[bug, molops, stereo]") {
    auto mol = "CC=CC"_smiles;
    REQUIRE(mol);

    auto bond = mol->getBondWithIdx(1);
    bond->setStereo(Bond::STEREOANY);
    REQUIRE(bond->getStereoAtoms().empty());

    MolOps::findPotentialStereoBonds(*mol);
    CHECK(bond->getStereo() == Bond::STEREOANY);
    CHECK(bond->getStereoAtoms().size() == 2);
}
