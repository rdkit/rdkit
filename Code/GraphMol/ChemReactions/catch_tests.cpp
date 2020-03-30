//
//  Copyright (c) 2018 Greg Landrum
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
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

using namespace RDKit;
using std::unique_ptr;

TEST_CASE("Github #1632", "[Reaction,PDB,bug]") {
  SECTION("basics") {
    bool sanitize = true;
    int flavor = 0;
    std::unique_ptr<RWMol> mol(SequenceToMol("K", sanitize, flavor));
    REQUIRE(mol);
    REQUIRE(mol->getAtomWithIdx(0)->getMonomerInfo());
    auto res = static_cast<AtomPDBResidueInfo*>(
        mol->getAtomWithIdx(0)->getMonomerInfo());
    CHECK(res->getResidueNumber() == 1);
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[O:1]=[CX3:2]-[CX4:3]-[NX3:4]>>[O:1]=[CX3:2]-[CX4:3]-[NX3:4]-[C]"));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    MOL_SPTR_VECT reacts;
    reacts.push_back(ROMOL_SPTR(new ROMol(*mol)));
    auto prods = rxn->runReactants(reacts);
    CHECK(prods.size() == 1);
    CHECK(prods[0].size() == 1);
    auto p = prods[0][0];
    CHECK(p->getNumAtoms() == mol->getNumAtoms() + 1);
    REQUIRE(p->getAtomWithIdx(0)->getMonomerInfo());
    auto pres = static_cast<AtomPDBResidueInfo*>(
        p->getAtomWithIdx(0)->getMonomerInfo());
    CHECK(pres->getResidueNumber() == 1);
    REQUIRE(!p->getAtomWithIdx(4)->getMonomerInfo());
  }
}

static void clearAtomMappingProps(ROMol& mol) {
  for (auto&& a : mol.atoms()) {
    a->clear();
  }
}

TEST_CASE("Github #2366 Enhanced Stereo", "[Reaction,StereoGroup,bug]") {
  SECTION("Reaction Preserves Stereo") {
    ROMOL_SPTR mol("F[C@H](Cl)Br |o1:1|"_smiles);
    REQUIRE(mol);
    unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("[C@:1]>>[C@:1]"));
    REQUIRE(rxn);

    MOL_SPTR_VECT reactants = {mol};

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reactants);
    REQUIRE(prods.size() == 1);
    REQUIRE(prods[0].size() == 1);
    auto p = prods[0][0];

    clearAtomMappingProps(*p);
    CHECK(MolToCXSmiles(*p) == "F[C@H](Cl)Br |o1:1|");
  }
  SECTION("Reaction destroys one center in StereoGroup") {
    ROMOL_SPTR mol("F[C@H](Cl)[C@@H](Cl)Br |&1:1,3|"_smiles);
    REQUIRE(mol);
    unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("[C@:1]F>>[C:1]F"));
    REQUIRE(rxn);

    MOL_SPTR_VECT reactants = {mol};

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reactants);
    REQUIRE(prods.size() == 1);
    REQUIRE(prods[0].size() == 1);
    auto p = prods[0][0];

    clearAtomMappingProps(*p);
    CHECK(MolToCXSmiles(*p) == "FC(Cl)[C@@H](Cl)Br |&1:3|");
  }
  SECTION("Reaction splits StereoGroup") {
    ROMOL_SPTR mol("F[C@H](Cl)[C@@H](Cl)Br |&1:1,3|"_smiles);
    REQUIRE(mol);
    unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[F:1][C@:2][C@:3][Cl:4]>>[F:1][C@:2]O.O[C@:3][Cl:4]"));
    REQUIRE(rxn);

    MOL_SPTR_VECT reactants = {mol};

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reactants);
    REQUIRE(prods.size() == 1);
    REQUIRE(prods[0].size() == 2);
    auto p0 = prods[0][0];
    auto p1 = prods[0][1];

    clearAtomMappingProps(*p0);
    clearAtomMappingProps(*p1);
    CHECK(MolToCXSmiles(*p0) == "O[C@@H](F)Cl |&1:1|");
    CHECK(MolToCXSmiles(*p1) == "O[C@@H](Cl)Br |&1:1|");
  }
  SECTION("Reaction combines StereoGroups") {
    ROMOL_SPTR mol1("F[C@H](Cl)O |&1:1|"_smiles);
    REQUIRE(mol1);
    ROMOL_SPTR mol2("Cl[C@H](Br)O |&1:1|"_smiles);
    REQUIRE(mol2);
    unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[F:1][C@:2]O.O[C@:3][Cl:4]>>[F:1][C@:2][C@:3][Cl:4]"));
    REQUIRE(rxn);

    MOL_SPTR_VECT reactants = {mol1, mol2};

    rxn->initReactantMatchers();
    auto prods = rxn->runReactants(reactants);
    REQUIRE(prods.size() == 1);
    REQUIRE(prods[0].size() == 1);
    auto p0 = prods[0][0];

    clearAtomMappingProps(*p0);
    CHECK(MolToCXSmiles(*p0) == "F[C@@H](Cl)[C@H](Cl)Br |&1:1,&2:3|");
  }
}

TEST_CASE("Github #2427 cannot set maxProducts>1000 in runReactants",
          "[Reaction,bug]") {
  SECTION("Basics") {
    std::string smi = "[C]";
    for (unsigned int i = 0; i < 49; ++i) {
      smi += ".[C]";
    }
    ROMOL_SPTR mol(SmilesToMol(smi));
    REQUIRE(mol);
    unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("([#6:1].[#6:2])>>[#6:1]-[#6:2]"));
    REQUIRE(rxn);

    MOL_SPTR_VECT reactants = {mol};

    rxn->initReactantMatchers();
    // by default we only get 1000 products:
    {
      auto prods = rxn->runReactants(reactants);
      CHECK(prods.size() == 1000);
      CHECK(prods[0].size() == 1);
    }
    {
      auto prods = rxn->runReactants(reactants, 2000);
      CHECK(prods.size() == 2000);
      CHECK(prods[0].size() == 1);
    }
  }
}

TEST_CASE("negative charge queries. Part of testing changes for github #2604",
          "[Reaction]") {
  SECTION("no redundancy") {
    unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("[N+{1-}:1]>>[#0-1:1]"));
    REQUIRE(rxn);

    // we don't have a way to directly create NegativeFormalCharge queries, so
    // make one by hand
    REQUIRE(rxn->getProducts()[0]->getAtomWithIdx(0)->hasQuery());
    static_cast<QueryAtom*>(rxn->getProducts()[0]->getAtomWithIdx(0))
        ->expandQuery(makeAtomNegativeFormalChargeQuery(1));
    unsigned nWarnings = 0;
    unsigned nErrors = 0;
    CHECK(rxn->validate(nWarnings, nErrors));
    CHECK(nWarnings == 0);
    CHECK(nErrors == 0);
  }
  SECTION("no redundancy2") {
    unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("[N+{1-}:1]>>[#0+1:1]"));
    REQUIRE(rxn);

    // we don't have a way to directly create NegativeFormalCharge queries, so
    // make one by hand
    REQUIRE(rxn->getProducts()[0]->getAtomWithIdx(0)->hasQuery());
    static_cast<QueryAtom*>(rxn->getProducts()[0]->getAtomWithIdx(0))
        ->expandQuery(makeAtomNegativeFormalChargeQuery(
            -1));  // a bit kludgy, but we need to check
    unsigned nWarnings = 0;
    unsigned nErrors = 0;
    CHECK(rxn->validate(nWarnings, nErrors));
    CHECK(nWarnings == 0);
    CHECK(nErrors == 0);
  }
  SECTION("redundancy") {
    unique_ptr<ChemicalReaction> rxn(
        // RxnSmartsToChemicalReaction("[N+{1-}:1]>>[#0+1+2:1]"));
        RxnSmartsToChemicalReaction("[N+{1-}:1]>>[#0-1:1]"));
    REQUIRE(rxn);
    // we don't have a way to directly create NegativeFormalCharge queries, so
    // make one by hand
    REQUIRE(rxn->getProducts()[0]->getAtomWithIdx(0)->hasQuery());
    static_cast<QueryAtom*>(rxn->getProducts()[0]->getAtomWithIdx(0))
        ->expandQuery(makeAtomNegativeFormalChargeQuery(2));
    unsigned nWarnings = 0;
    unsigned nErrors = 0;
    CHECK(rxn->validate(nWarnings, nErrors));
    CHECK(nWarnings == 1);
    CHECK(nErrors == 0);
  }
}

TEST_CASE("GithHub #2954: Reaction Smarts with Dative Bonds not parsed",
          "[Reaction, Bug]") {
  
  SECTION("Rxn Smart Processing with Dative Bond") {
    unique_ptr<ChemicalReaction> rxn(
      RxnSmartsToChemicalReaction("[O:1].[H+]>>[O:1]->[H+]")
    );
    REQUIRE(rxn);
    auto k = rxn->getProducts()[0]->getNumAtoms();
    CHECK(k == 2);
  }

}