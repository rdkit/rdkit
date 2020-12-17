//
//  Copyright (c) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
///
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;
using std::unique_ptr;

TEST_CASE("Github #1039", "[]") {
  SECTION("double bond") {
    auto m1 = "C/C=C/C=C/C"_smiles;
    REQUIRE(m1);
    std::vector<unsigned int> bonds = {2};
    std::unique_ptr<ROMol> pieces(MolFragmenter::fragmentOnBonds(*m1, bonds));
    REQUIRE(pieces);
    CHECK(pieces->getNumAtoms() == 8);
    REQUIRE(pieces->getBondBetweenAtoms(3, 6));
    REQUIRE(pieces->getBondBetweenAtoms(2, 7));
    CHECK(pieces->getBondBetweenAtoms(3, 6)->getBondType() == Bond::SINGLE);
    CHECK(pieces->getBondBetweenAtoms(3, 6)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(pieces->getBondBetweenAtoms(2, 7)->getBondType() == Bond::SINGLE);
    CHECK(pieces->getBondBetweenAtoms(2, 7)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(MolToSmiles(*pieces) == "[2*]/C=C/C.[3*]/C=C/C");
  }
  SECTION("atomic stereo") {
    auto m1 = "C(C)(F)(Cl)O"_smiles;
    REQUIRE(m1);
    m1->getBondWithIdx(0)->setBondDir(Bond::BEGINWEDGE);
    std::vector<unsigned int> bonds = {0};
    std::unique_ptr<ROMol> pieces(MolFragmenter::fragmentOnBonds(*m1, bonds));
    REQUIRE(pieces);
    CHECK(pieces->getNumAtoms() == 7);
    REQUIRE(pieces->getBondBetweenAtoms(0, 6));
    REQUIRE(pieces->getBondBetweenAtoms(1, 5));
    CHECK(pieces->getBondBetweenAtoms(0, 6)->getBondDir() == Bond::BEGINWEDGE);
    CHECK(pieces->getBondBetweenAtoms(1, 5)->getBondDir() == Bond::NONE);
    // no actual stereo in the SMILES here since we haven't assigned it (need a
    // conformer to do that using wedging)
    CHECK(MolToSmiles(*pieces) == "*C.[1*]C(O)(F)Cl");
  }
  SECTION("bond stereo") {
      auto m =  "O/C=N/C=C"_smiles;
      std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1,1}};
      std::vector<unsigned int> bonds{0};
      auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
      CHECK(MolToSmiles(*resa) == "*/C=N/C=C.[1*]O");
      // make sure we still have stereo atoms
      std::vector<std::vector<int>> expected_stereo_atoms {
          {5,3}, // 5 is the new dummy atom, it was 0 before
          {},
          {},
          {},
          {},
      };
      std::vector<std::vector<int>> received_stereo;
      for(auto *bond: resa->bonds()) {
          received_stereo.push_back(bond->getStereoAtoms());
      }
      CHECK(received_stereo==expected_stereo_atoms);
  }
  { // break non stereo atom bond
    auto m =  "C/C(O)=N/C=C"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1,1}};
    std::vector<unsigned int> bonds{0};
    auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
    CHECK(MolToSmiles(*resa) == "*/C(O)=N/C=C.[1*]C");
    // make sure we still have stereo atoms
    std::vector<std::vector<int>> expected_stereo_atoms {
							 {},
							 {2,4},
							 {},
							 {},
							 {},
							 {}
    };
    std::vector<std::vector<int>> received_stereo;
    for(auto *bond: resa->bonds()) {
      received_stereo.push_back(bond->getStereoAtoms());
    }
    CHECK(received_stereo==expected_stereo_atoms);
  }  
  { // bond stereo should only be removed when deleting the double bond with E/Z
    auto m =  "O/C=N/C=C"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1,1}};
    std::vector<std::string> expected = {
					 "*/C=N/C=C.[1*]O",
					 "[1*]=NC=C.[2*]=CO", // bond stereo gone
					 "[2*]C=C.[3*]/N=C/O",
					 "[3*]=C.[4*]=C/N=C/O"
    };
    for(unsigned int i=0;i<m->getNumBonds();++i) {
      std::vector<unsigned int> bonds{i};
      auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
      auto smiles = MolToSmiles(*resa);
      CHECK(smiles == expected[i]);
    }
  }
  { // bond stereo should only be removed when deleting the double bond with E/Z
    // chiral stereo should stay
    auto m =  "O/C=N/[C@H](I)F"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1,1}};
    std::vector<std::string> expected = {
					 "*/C=N/[C@@H](F)I.[1*]O",
					 "[1*]=N[C@@H](F)I.[2*]=CO", // bond stereo gone
					 "[2*][C@@H](F)I.[3*]/N=C/O",
					 "[3*]I.[4*][C@H](F)/N=C/O",
					 "[3*]F.[5*][C@@H](I)/N=C/O"
    };
    for(unsigned int i=0;i<m->getNumBonds();++i) {
      std::vector<unsigned int> bonds{i};
      auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
      auto smiles = MolToSmiles(*resa);
      CHECK(smiles == expected[i]);
    }
  }
}

