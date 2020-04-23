//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("possible stereochemistry on bonds", "[molops]") {
  SECTION("simplest") {
    {
      auto mol = "CC=CC"_smiles;
      REQUIRE(mol);
      bool cleanIt = true;
      bool force = true;
      bool flagPossibleStereoCenters = true;
      MolOps::assignStereochemistry(*mol, cleanIt, force,
                                    flagPossibleStereoCenters);
      CHECK(mol->getBondWithIdx(1)->hasProp(
          common_properties::_ChiralityPossible));
    }
    {
      auto mol = "CC=C(C)C"_smiles;
      REQUIRE(mol);
      bool cleanIt = true;
      bool force = true;
      bool flagPossibleStereoCenters = true;
      MolOps::assignStereochemistry(*mol, cleanIt, force,
                                    flagPossibleStereoCenters);
      CHECK(!mol->getBondWithIdx(1)->hasProp(
          common_properties::_ChiralityPossible));
    }
    {
      auto mol = "CC=C"_smiles;
      REQUIRE(mol);
      bool cleanIt = true;
      bool force = true;
      bool flagPossibleStereoCenters = true;
      MolOps::assignStereochemistry(*mol, cleanIt, force,
                                    flagPossibleStereoCenters);
      CHECK(!mol->getBondWithIdx(1)->hasProp(
          common_properties::_ChiralityPossible));
    }
  }
}

TEST_CASE("para-stereocenters and assignStereochemistry", "[molops]") {
  SECTION("simplest") {
    auto mol = "CC(F)CC(C)F"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool force = true;
    bool flagPossibleStereoCenters = true;
    MolOps::assignStereochemistry(*mol, cleanIt, force,
                                  flagPossibleStereoCenters);
    CHECK(
        mol->getAtomWithIdx(1)->hasProp(common_properties::_ChiralityPossible));
    CHECK(
        mol->getAtomWithIdx(4)->hasProp(common_properties::_ChiralityPossible));
    CHECK(
        mol->getAtomWithIdx(3)->hasProp(common_properties::_ChiralityPossible));
  }

  SECTION("including bonds") {
    // thanks to Salome Rieder for this nasty example
    auto mol = "CC=CC(C=CC)C(C)C(C=CC)C=CC"_smiles;
    REQUIRE(mol);
    CHECK(
        mol->getAtomWithIdx(7)->hasProp(common_properties::_ChiralityPossible));
    CHECK(
        mol->getAtomWithIdx(3)->hasProp(common_properties::_ChiralityPossible));
    CHECK(
        mol->getAtomWithIdx(9)->hasProp(common_properties::_ChiralityPossible));
    CHECK(mol->getBondBetweenAtoms(13, 14)->hasProp(
        common_properties::_ChiralityPossible));
    CHECK(mol->getBondBetweenAtoms(4, 5)->hasProp(
        common_properties::_ChiralityPossible));
    CHECK(mol->getBondBetweenAtoms(1, 2)->hasProp(
        common_properties::_ChiralityPossible));
    CHECK(mol->getBondBetweenAtoms(10, 11)->hasProp(
        common_properties::_ChiralityPossible));
  }
}
