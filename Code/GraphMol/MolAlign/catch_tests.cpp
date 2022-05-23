//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "catch.hpp"
#include "AlignMolecules.h"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

using namespace RDKit;

TEST_CASE("symmetric functional groups") {
  SECTION("basics") {
    auto m1 =
        "CCC(=O)[O-] "
        "|(-1.11,0.08,-0.29;0.08,-0.18,0.58;1.34,0.03,-0.16;1.74,1.22,-0.32;2.06,-1.04,-0.66)|"_smiles;
    REQUIRE(m1);

    // swap the bond orders to the Os
    RWMol m2(*m1);
    m2.getAtomWithIdx(3)->setFormalCharge(-1);
    m2.getAtomWithIdx(4)->setFormalCharge(0);
    m2.getBondBetweenAtoms(2, 3)->setBondType(Bond::BondType::SINGLE);
    m2.getBondBetweenAtoms(2, 4)->setBondType(Bond::BondType::DOUBLE);

    {
      auto rmsd = MolAlign::getBestRMS(m2, *m1);
      CHECK(rmsd == Approx(0.0).margin(1e-3));
    }
    {
      // previous behavior
      int probeId = -1;
      int refId = -1;
      std::vector<MatchVectType> mp;
      int maxMatches = 1e6;
      bool symmetrize = false;
      auto rmsd = MolAlign::getBestRMS(m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize);
      CHECK(rmsd == Approx(0.747).margin(1e-3));
    }
  }
  SECTION("terminal sulfate1") {
    auto m1 =
        "CS(=O)(=O)[O-] |(-0.93,-0.06,-0.04;0.82,0.07,0.13;1.27,-0.04,1.54;1.21,1.40,-0.48;1.53,-1.11,-0.82)|"_smiles;
    REQUIRE(m1);

    // swap the bond orders to the Os
    RWMol m2(*m1);
    m2.getAtomWithIdx(2)->setFormalCharge(-1);
    m2.getAtomWithIdx(4)->setFormalCharge(0);
    m2.getBondBetweenAtoms(1, 2)->setBondType(Bond::BondType::SINGLE);
    m2.getBondBetweenAtoms(1, 4)->setBondType(Bond::BondType::DOUBLE);

    {
      auto rmsd = MolAlign::getBestRMS(m2, *m1);
      CHECK(rmsd == Approx(0.0).margin(1e-3));
    }
    {
      // previous behavior
      int probeId = -1;
      int refId = -1;
      std::vector<MatchVectType> mp;
      int maxMatches = 1e6;
      bool symmetrize = false;
      auto rmsd = MolAlign::getBestRMS(m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize);
      CHECK(rmsd == Approx(0.097).margin(1e-3));
    }
  }
  SECTION("terminal sulfate2") {
    auto m1 =
        "CS(=O)(=O)[O-] |(-0.93,-0.06,-0.04;0.82,0.07,0.13;1.27,-0.04,1.54;1.21,1.40,-0.48;1.53,-1.11,-0.82)|"_smiles;
    REQUIRE(m1);

    // swap the bond orders to the Os
    RWMol m2(*m1);
    m2.getAtomWithIdx(3)->setFormalCharge(-1);
    m2.getAtomWithIdx(4)->setFormalCharge(0);
    m2.getBondBetweenAtoms(1, 3)->setBondType(Bond::BondType::SINGLE);
    m2.getBondBetweenAtoms(1, 4)->setBondType(Bond::BondType::DOUBLE);

    {
      auto rmsd = MolAlign::getBestRMS(m2, *m1);
      CHECK(rmsd == Approx(0.0).margin(1e-3));
    }
    {
      // previous behavior
      int probeId = -1;
      int refId = -1;
      std::vector<MatchVectType> mp;
      int maxMatches = 1e6;
      bool symmetrize = false;
      auto rmsd = MolAlign::getBestRMS(m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize);
      CHECK(rmsd == Approx(0.097).margin(1e-3));
    }
  }
}