//
//  Copyright (C) 2025 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/test_fixtures.h>

using namespace RDKit;

TEST_CASE("explicit tests") {
  // making some of the other atropisomer tests explicitly test what we expect
  std::string rdbase = getenv("RDBASE");
  REQUIRE(rdbase != "");
  SECTION("basics") {
    std::string file =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/atropisomers/AtropManyChirals.sdf";
    auto mol = v2::FileParsers::MolFromMolFile(file);
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(7)->getStereo() ==
          Bond::BondStereo::STEREOATROPCCW);
  }
  SECTION("macrocycle") {
    std::string file =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/atropisomers/macrocycle-5-meta-Cl-ortho-hash.sdf";
    auto mol = v2::FileParsers::MolFromMolFile(file);
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(15)->getStereo() ==
          Bond::BondStereo::STEREOATROPCW);
  }
  SECTION("with stereo groups and other chiral centers") {
    REQUIRE(rdbase != "");
    std::string file =
        rdbase +
        "/Code/GraphMol/FileParsers/test_data/atropisomers/AtropManyChiralsEnhanced.sdf";
    auto mol = v2::FileParsers::MolFromMolFile(file);
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(7)->getStereo() ==
          Bond::BondStereo::STEREOATROPCCW);
    CHECK(mol->getAtomWithIdx(6)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(9)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(13)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(16)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(10)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);

    REQUIRE(mol->getStereoGroups().size() == 2);
    {
      const auto &sg = mol->getStereoGroups()[0];
      CHECK(sg.getAtoms().size() == 2);
      CHECK(sg.getBonds().size() == 1);
    }
    {
      const auto &sg = mol->getStereoGroups()[1];
      CHECK(sg.getAtoms().size() == 2);
      CHECK(sg.getBonds().size() == 0);
    }
  }
}

TEST_CASE("atropisomerism and dependent chirality") {
  SECTION("chirality allows stereoisomers") {
    {
      auto m = "C=C(C)C1=CC([C@H](C)O)=CC([C@@H](C)O)=C1 |wD:3.14|"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(2)->getStereo() ==
            Bond::BondStereo::STEREOATROPCW);
    }
#if 0
    // FIX: this test does not work, this is Github #8340
    {
      // the two stereocenters are identical, we should not have atropisomerism
      auto m = "C=C(C)C1=CC([C@H](C)O)=CC([C@H](C)O)=C1 |wD:3.14|"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(2)->getStereo() == Bond::BondStereo::STEREONONE);
    }
#endif
  }

  SECTION("atropisomers allow chirality") {
#if 0
    // FIX: this test does not work
    {
      // here the atropisomeric bonds have opposite stereochemistry, so the
      // central atom should be chiral, this is Github #8341
      auto m =
          "C=C(C1C=C(C)C=CC=1)[C@H](C)C(=C)C1C=C(C)C=CC=1 |wD:2.2,wD:13.14|"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREOATROPCCW);
      CHECK(m->getBondWithIdx(12)->getStereo() ==
            Bond::BondStereo::STEREOATROPCCW);
      CHECK(m->getAtomWithIdx(9)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
#endif
    {
      // here the atropisomeric bonds have the same stereochemistry, so the
      // central atom should not be chiral
      auto m =
          "C=C(C1C=C(C)C=CC=1)[C@H](C)C(=C)C1C=C(C)C=CC=1 |wD:2.2,wU:13.14|"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREOATROPCCW);
      CHECK(m->getBondWithIdx(12)->getStereo() ==
            Bond::BondStereo::STEREOATROPCW);
      CHECK(m->getAtomWithIdx(9)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
}

TEST_CASE("Github #8602: atropisomer with slight z coordinates") {
  std::string rdbase = getenv("RDBASE");
  rdbase += "/Code/GraphMol/FileParsers/test_data/";
  SECTION("as reported") {
    auto m = v2::FileParsers::MolFromMolFile(rdbase + "atropWith3d_8602.mol");
    REQUIRE(m);
      for (auto bond : m->bonds()) {
        // we read in the molecule, but the bonds should be co-planar
        // in the current implementation, where 2D mol in 3D does not
        // fall back to a 2D bond explicit annotation check
        CHECK(bond->getStereo() != Bond::BondStereo::STEREOATROPCW);
        CHECK(bond->getStereo() != Bond::BondStereo::STEREOATROPCCW);
      }
  }
}
