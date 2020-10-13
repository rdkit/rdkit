//
//  Copyright (c) 2018-2020 Greg Landrum
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
#include <GraphMol/FileParsers/PNGParser.h>

using namespace RDKit;
using std::unique_ptr;

TEST_CASE("Github #1632", "[Reaction][PDB][bug]") {
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

TEST_CASE("Github #2366 Enhanced Stereo", "[Reaction][StereoGroup][bug]") {
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
    CHECK(MolToCXSmiles(*p0) == "O[C@H](F)Cl |&1:1|");
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
    CHECK(MolToCXSmiles(*p0) == "F[C@H](Cl)[C@H](Cl)Br |&1:1,&2:3|");
  }
}

TEST_CASE("Github #2427 cannot set maxProducts>1000 in runReactants",
          "[Reaction][bug]") {
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
          "[Reaction][Bug]") {
  SECTION("Rxn Smart Processing with Dative Bond in Product") {
    unique_ptr<ChemicalReaction> rxn1(
        RxnSmartsToChemicalReaction("[O:1].[H+]>>[O:1]->[H+]"));

    REQUIRE(rxn1);
    auto k = rxn1->getProducts()[0]->getNumAtoms();
    CHECK(k == 2);
  }

  SECTION("Rxn Smart Processing with Dative Bond in Reactant") {
    unique_ptr<ChemicalReaction> rxn2(
        RxnSmartsToChemicalReaction("[O:1]->[H+]>>[O:1].[H+]"));

    REQUIRE(rxn2);

    auto k = rxn2->getReactants()[0]->getNumAtoms();
    CHECK(k == 2);
  }

  SECTION("Rxm Smart Processing with Dative Bond in Agent") {
    unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("[O:1][H]>N->[Cu]>[O:1].[H]"));
    REQUIRE(rxn);

    auto k = rxn->getAgents()[0]->getNumAtoms();
    CHECK(k == 2);
  }
}

TEST_CASE("GithHub #3119: partial reacting atom detection", "[Reaction][Bug]") {
  SECTION("as reported") {
    bool useSmiles = true;
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[CH:1]1=[CH:2][CH:3]=[C:4]([C:5](=[CH:6]1)[CH2:7]Cl)[O:8][CH2:9]"
        "[CH2:10]Cl>>[CH:1]1=[CH:2][CH:3]=[C:4]([C:5](=[CH:6]1)[CH2:7]I)["
        "O:8][CH2:9][CH2:10]I",
        nullptr, useSmiles));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    bool mappedAtomsOnly = true;
    auto rAtoms = getReactingAtoms(*rxn, mappedAtomsOnly);
    REQUIRE(rAtoms.size() == 1);
    CHECK(rAtoms[0].size() == 2);
    CHECK(rAtoms[0][0] == 6);
    CHECK(rAtoms[0][1] == 10);
  }
}

TEST_CASE("reaction data in PNGs 1", "[Reaction][PNG]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ChemReactions/testData/";
  std::string sma =
      "[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[Br].[#6:7]B(O)O>>[cH:1]1[cH:2]["
      "cH:3][cH:4][cH:5][c:6]1-[#6:7]";
  SECTION("add/read metadata") {
    std::vector<std::pair<std::string, std::string>> metadata;
    metadata.push_back(std::make_pair(PNGData::rxnSmartsTag, sma));
    std::string fname = pathName + "reaction1.no_metadata.png";
    bool compressed = false;
    auto pngData = addMetadataToPNGFile(fname, metadata, compressed);
    metadata = PNGStringToMetadata(pngData);
    auto iter =
        std::find_if(metadata.begin(), metadata.end(),
                     [](const std::pair<std::string, std::string>& val) {
                       return val.first == PNGData::rxnSmartsTag;
                     });
    REQUIRE(iter != metadata.end());
  }
  SECTION("read from file") {
    std::string fname = pathName + "reaction1.smarts.png";
    std::unique_ptr<ChemicalReaction> rxn(PNGFileToChemicalReaction(fname));
    REQUIRE(rxn);
    CHECK(rxn->getNumReactantTemplates() == 2);
    CHECK(rxn->getNumProductTemplates() == 1);
  }
  SECTION("add/read reaction to/from file") {
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(sma));
    REQUIRE(rxn);
    std::string fname = pathName + "reaction1.no_metadata.png";
    {
      auto pngData = addChemicalReactionToPNGFile(*rxn, fname);
      std::unique_ptr<ChemicalReaction> nrxn(
          PNGStringToChemicalReaction(pngData));
      REQUIRE(nrxn);
      CHECK(nrxn->getNumReactantTemplates() == 2);
      CHECK(nrxn->getNumProductTemplates() == 1);
    }
    {
      bool includePkl = true;
      bool includeSmiles = false;
      bool includeSmarts = false;
      bool includeRxn = false;
      auto pngData = addChemicalReactionToPNGFile(
          *rxn, fname, includePkl, includeSmiles, includeSmarts, includeRxn);
      std::unique_ptr<ChemicalReaction> nrxn(
          PNGStringToChemicalReaction(pngData));
      REQUIRE(nrxn);
      CHECK(nrxn->getNumReactantTemplates() == 2);
      CHECK(nrxn->getNumProductTemplates() == 1);
    }
    {
      bool includePkl = false;
      bool includeSmiles = true;
      bool includeSmarts = false;
      bool includeRxn = false;
      auto pngData = addChemicalReactionToPNGFile(
          *rxn, fname, includePkl, includeSmiles, includeSmarts, includeRxn);
      std::unique_ptr<ChemicalReaction> nrxn(
          PNGStringToChemicalReaction(pngData));
      REQUIRE(nrxn);
      CHECK(nrxn->getNumReactantTemplates() == 2);
      CHECK(nrxn->getNumProductTemplates() == 1);
    }
    {
      bool includePkl = false;
      bool includeSmiles = false;
      bool includeSmarts = true;
      bool includeRxn = false;
      auto pngData = addChemicalReactionToPNGFile(
          *rxn, fname, includePkl, includeSmiles, includeSmarts, includeRxn);
      std::unique_ptr<ChemicalReaction> nrxn(
          PNGStringToChemicalReaction(pngData));
      REQUIRE(nrxn);
      CHECK(nrxn->getNumReactantTemplates() == 2);
      CHECK(nrxn->getNumProductTemplates() == 1);
    }
    {
      bool includePkl = false;
      bool includeSmiles = false;
      bool includeSmarts = false;
      bool includeRxn = true;
      auto pngData = addChemicalReactionToPNGFile(
          *rxn, fname, includePkl, includeSmiles, includeSmarts, includeRxn);
      std::unique_ptr<ChemicalReaction> nrxn(
          PNGStringToChemicalReaction(pngData));
      REQUIRE(nrxn);
      CHECK(nrxn->getNumReactantTemplates() == 2);
      CHECK(nrxn->getNumProductTemplates() == 1);
    }
  }
}

TEST_CASE("Github #2891", "[Reaction][chirality][bug]") {
  SECTION("reaction parsing inversion logic") {
    std::vector<std::pair<std::string, int>> tests{
        {"[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])[S:3]", 2},
        {"[C:4][C@@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])[S:3]", 1},
        {"[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@:2]([S:3])[F:1]", 1},
        {"[C:4][C@@:2]([F:1])[Br:3]>>[C:4][C@:2]([S:3])[F:1]", 2},
        // add mapped substituents
        {"[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([S:3])[Cl:5]", 2},
        {"[C:4][C@@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([S:3])[Cl:5]", 1},
        {"[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([Cl:5])[S:3]", 1},
        {"[C:4][C@@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([Cl:5])[S:3]", 2},
        // remove mapped substituents
        {"[C:4][C@:2]([F:1])([Br:3])[Cl:5]>>[C:4][C@:2]([F:1])[S:3]", 2},
        {"[C:4][C@@:2]([F:1])([Br:3])[Cl:5]>>[C:4][C@:2]([F:1])[S:3]", 1},
        {"[C:4][C@:2]([F:1])([Cl:5])[Br:3]>>[C:4][C@:2]([F:1])[S:3]", 1},
        {"[C:4][C@@:2]([F:1])([Cl:5])[Br:3]>>[C:4][C@:2]([F:1])[S:3]", 2},
        // add unmapped substituents
        {"[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([S:3])[Cl]", 2},
        {"[C:4][C@@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([S:3])[Cl]", 1},
        {"[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([Cl])[S:3]", 1},
        {"[C:4][C@@:2]([F:1])[Br:3]>>[C:4][C@:2]([F:1])([Cl])[S:3]", 2},
        // remove unmapped substituents
        {"[C:4][C@:2]([F:1])([Br:3])[Cl]>>[C:4][C@:2]([F:1])[S:3]", 2},
        {"[C:4][C@@:2]([F:1])([Br:3])[Cl]>>[C:4][C@:2]([F:1])[S:3]", 1},
        {"[C:4][C@:2]([F:1])([Cl])[Br:3]>>[C:4][C@:2]([F:1])[S:3]", 1},
        {"[C:4][C@@:2]([F:1])([Cl])[Br:3]>>[C:4][C@:2]([F:1])[S:3]", 2},
    };
    for (const auto& pr : tests) {
      std::unique_ptr<ChemicalReaction> rxn(
          RxnSmartsToChemicalReaction(pr.first));
      REQUIRE(rxn);
      auto minv = rxn->getProducts()[0]->getAtomWithIdx(1)->getProp<int>(
          common_properties::molInversionFlag);
      if (minv != pr.second) {
        std::cerr << pr.first << std::endl;
      }
      CHECK(minv == pr.second);
    }
  }
  SECTION("simplified") {
    auto r1_1 = ROMOL_SPTR(SmilesToMol("C[C@H](F)Br"));
    auto r1_2 = ROMOL_SPTR(SmilesToMol("C[C@@H](F)Br"));
    auto r1_3 = ROMOL_SPTR(SmilesToMol("C[C@@H](Br)F"));
    auto r1_4 = ROMOL_SPTR(SmilesToMol("C[C@H](Br)F"));

    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@@:2]([F:1])[S:3]"));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    {
      MOL_SPTR_VECT reacts{r1_1};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_3};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_4};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
  }

  SECTION("simplified2") {
    auto r1_1 = ROMOL_SPTR(SmilesToMol("C[C@H](F)Br"));
    auto r1_2 = ROMOL_SPTR(SmilesToMol("C[C@@H](F)Br"));
    auto r1_3 = ROMOL_SPTR(SmilesToMol("C[C@@H](Br)F"));
    auto r1_4 = ROMOL_SPTR(SmilesToMol("C[C@H](Br)F"));

    // makes sure we also handle swaps in the atom ordering properly
    // this isn't an inversion, despite going @->@@, because the atom order
    // also changes
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@@:2]([S:3])[F:1]"));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    CHECK(rxn->getProducts()[0]->getAtomWithIdx(1)->getProp<int>(
              common_properties::molInversionFlag) == 2);
    {
      MOL_SPTR_VECT reacts{r1_1};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_3};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_4};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@H](F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
  }

  SECTION("simplified3") {
    auto r1_1 = ROMOL_SPTR(SmilesToMol("C[C@H](F)Br"));
    auto r1_2 = ROMOL_SPTR(SmilesToMol("C[C@@H](F)Br"));
    auto r1_3 = ROMOL_SPTR(SmilesToMol("C[C@@H](Br)F"));
    auto r1_4 = ROMOL_SPTR(SmilesToMol("C[C@H](Br)F"));

    // adding a bond in the products.
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[C:4][C@:2]([F:1])[Br:3]>>[C:4][C@@:2](O)([F:1])[S:3]"));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    CHECK(rxn->getProducts()[0]->getAtomWithIdx(1)->getProp<int>(
              common_properties::molInversionFlag) == 1);
    {
      MOL_SPTR_VECT reacts{r1_1};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@](O)(F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@](O)(F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_3};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@](O)(F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_4};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@](O)(F)S"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
  }

  SECTION("reported1") {
    auto r1_1 = ROMOL_SPTR(SmilesToMol("O[C@@H](Br)c1ccccc1"));
    auto r1_2 = ROMOL_SPTR(SmilesToMol("O[C@H](Br)c1ccccc1"));
    auto r1_3 = ROMOL_SPTR(SmilesToMol("O[C@@H](c1ccccc1)Br"));
    auto r1_4 = ROMOL_SPTR(SmilesToMol("O[C@H](c1ccccc1)Br"));
    auto r2 = ROMOL_SPTR(SmilesToMol("SCC"));

    std::unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction("[O:5][C@:1]([Br:2])[*:6].[S:3][C:4]>>[*:"
                                    "5][C@@:1]([S:3][*:4])[*:6]"));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    {
      MOL_SPTR_VECT reacts{r1_1, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("O[C@H](SCC)c1ccccc1"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_2, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("O[C@@H](SCC)c1ccccc1"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_3, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("O[C@H](c1ccccc1)SCC"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_4, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("O[C@@H](c1ccccc1)SCC"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
  }

  SECTION("reported2") {
    auto r1_1 = ROMOL_SPTR(SmilesToMol("C[C@@H](Br)c1ccccc1"));
    auto r1_2 = ROMOL_SPTR(SmilesToMol("C[C@H](Br)c1ccccc1"));
    auto r1_3 = ROMOL_SPTR(SmilesToMol("C[C@@H](c1ccccc1)Br"));
    auto r1_4 = ROMOL_SPTR(SmilesToMol("C[C@H](c1ccccc1)Br"));
    auto r2 = ROMOL_SPTR(SmilesToMol("SCC"));
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[C@:1][Br:2].[S:3][C:4]>>[C@@:1][*:3][*:4]"));
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    {
      MOL_SPTR_VECT reacts{r1_1, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@H](SCC)c1ccccc1"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_2, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@H](SCC)c1ccccc1"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_3, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@H](c1ccccc1)SCC"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
    {
      MOL_SPTR_VECT reacts{r1_4, r2};
      auto ps = rxn->runReactants(reacts);
      CHECK(ps.size() == 1);
      CHECK(ps[0].size() == 1);
      auto tsmi = MolToSmiles(*("C[C@@H](c1ccccc1)SCC"_smiles));
      CHECK(MolToSmiles(*ps[0][0]) == tsmi);
    }
  }
}
