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
    auto res = static_cast<AtomPDBResidueInfo *>(
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
    auto pres = static_cast<AtomPDBResidueInfo *>(
        p->getAtomWithIdx(0)->getMonomerInfo());
    CHECK(pres->getResidueNumber() == 1);
    REQUIRE(!p->getAtomWithIdx(4)->getMonomerInfo());
  }
}

// static void clearAtomProps()

TEST_CASE("Github #2366 Enhanced Stereo", "[Reaction,StereoGroup,bug]") {
    SECTION("Reaction Preserves Stereo") {
        ROMOL_SPTR mol("F[C@H](Cl)Br |o1:1|"_smiles);
        REQUIRE(mol);
        unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction("[C@:1]>>[C@:1]"));
        REQUIRE(rxn);

        MOL_SPTR_VECT reactants = {mol};

        rxn->initReactantMatchers();
        auto prods = rxn->runReactants(reactants);
        REQUIRE(prods.size() == 1);
        REQUIRE(prods[0].size() == 1);
        auto p = prods[0][0];

        // clear mapping properties
        for (auto&& a: p->atoms()) {
            a->clear();
        }
        CHECK(MolToCXSmiles(*p) == "F[C@H](Cl)Br |o1:1|");
    }
    SECTION("Reaction destroys one center in StereoGroup") {
        ROMOL_SPTR mol("F[C@H](Cl)[C@@H](Cl)Br |&1:1,3|"_smiles);
        REQUIRE(mol);
        unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction("[C@:1]F>>[C:1]F"));
        REQUIRE(rxn);

        MOL_SPTR_VECT reactants = {mol};

        rxn->initReactantMatchers();
        auto prods = rxn->runReactants(reactants);
        REQUIRE(prods.size() == 1);
        REQUIRE(prods[0].size() == 1);
        auto p = prods[0][0];
        for (auto&& a: p->atoms()) {
            a->clear();
        }
        CHECK(MolToCXSmiles(*p) == "FC(Cl)[C@@H](Cl)Br |&1:3|");

    }
}
