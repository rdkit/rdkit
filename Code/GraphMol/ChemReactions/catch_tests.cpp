//
//  Copyright (c) 2018-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
///
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/FileParsers/FileParserUtils.h>

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

TEST_CASE("reaction literals") {
  {
    auto rxn = "[O:1]>>[N:1]"_rxnsmarts;
    CHECK(rxn != nullptr);
  }
  {
    auto rxn = "[O:1]>>[N:1]"_rxnsmiles;
    CHECK(rxn != nullptr);
  }
  {
    auto rxn = "CC1>>CC1"_rxnsmarts;
    CHECK(rxn == nullptr);
  }
  {
    auto rxn = "CC1>>CC1"_rxnsmiles;
    CHECK(rxn == nullptr);
  }
  {
    auto rxn = "Cc1cc1>>CCC"_rxnsmarts;
    CHECK(rxn != nullptr);
  }
  {
    auto rxn = "Cc1cc1>>CCC"_rxnsmiles;
    CHECK(rxn != nullptr);
  }
}
TEST_CASE("one-component reactions") {
  SECTION("removing atoms") {
    auto rxn = "CC[N:1]>>[N:1]"_rxnsmarts;
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    {  // molecule which does not match:
      auto mol = "CCO"_smiles;
      REQUIRE(mol);
      CHECK(!rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CCO");
    }
    {
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 1);
      CHECK(MolToSmiles(*mol) == "N");
    }
    {
      auto mol = "CCNC"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 2);
      CHECK(MolToSmiles(*mol) == "CN");
    }
    {
      auto mol = "NCC"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 1);
      CHECK(MolToSmiles(*mol) == "N");
    }
    {
      auto mol = "CNCC"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 2);
      CHECK(MolToSmiles(*mol) == "CN");
    }
    {
      auto mol = "CCCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 1);
      CHECK(MolToSmiles(*mol) == "N");
    }
    {
      auto mol = "CCCN(C)CO"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 4);
      CHECK(MolToSmiles(*mol) == "CNCO");
    }
    {  // multiple matches, we only modify one (and it's arbitrary which)
      auto mol = "CCNCNCC"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 5);
      CHECK(MolToSmiles(*mol) == "CCNCN");
    }
    {  // fragments don't pass through:
      auto mol = "CCN.Cl"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 1);
      CHECK(MolToSmiles(*mol) == "N");
    }
  }
  SECTION("removing atoms 2") {
    auto rxn = "[C:2][N:1]>>[N:1]"_rxnsmarts;
    REQUIRE(rxn);
    rxn->initReactantMatchers();
    {
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 1);
      CHECK(MolToSmiles(*mol) == "N");
    }
  }
  SECTION("unmapped atoms in the product is an error") {
    {
      auto rxn = "[N:1]>>[N:1]CC"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "N"_smiles;
      REQUIRE(mol);
      CHECK_THROWS_AS(rxn->runReactant(*mol), ChemicalReactionException);
    }
    {
      auto rxn = "[N:1]>>[N:1][C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "N"_smiles;
      REQUIRE(mol);
      CHECK_THROWS_AS(rxn->runReactant(*mol), ChemicalReactionException);
    }
  }
  SECTION("modifying atoms") {
    {
      auto rxn = "[C:2][N:1]>>[C:2][O:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      MolOps::sanitizeMol(*mol);
      CHECK(MolToSmiles(*mol) == "CCO");
    }
#if 0
// does not currently work properly either here or in the main reaction code
    {
      auto rxn = "[C:2][N+1:1]>>[C:2][N+0:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[NH3+]"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      MolOps::sanitizeMol(*mol);
      CHECK(MolToSmiles(*mol) == "CCN");
    }
#endif
    {
      auto rxn = "[C:2][N+0:1]>>[C:2][N+1:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      MolOps::sanitizeMol(*mol);
      CHECK(MolToSmiles(*mol) == "CC[NH3+]");
    }
    {
      auto rxn = "[C:2][N:1]>>[15N:1][C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      MolOps::sanitizeMol(*mol);
      CHECK(MolToSmiles(*mol) == "CC[15NH2]");
    }
    {
      auto rxn = "[C:2][15N:1]>>[0N:1][C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[15NH2]"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      MolOps::sanitizeMol(*mol);
      CHECK(MolToSmiles(*mol) == "CCN");
    }
  }
  SECTION("modifying bonds") {
    {
      auto rxn = "[C:2]-[N:1]>>([N:1].[C:2])"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CC.N");
    }
    {
      auto rxn = "([C:2].[N:1])>>[N:1]-[C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC.N"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CCN");
    }
    {
      auto rxn = "([CH3:2].[N:1])>>[N:1]-[C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "C1CN1");
    }
    {
      auto rxn = "([CH2:2].[N:1])>>[N:1]-[C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(!rxn->runReactant(*mol));
    }
    {
      auto rxn = "([CH2:2].[N:1])>>[N:1]=[C:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CC=N");
    }
    {
      auto rxn = "[C:2]-[N:1]>>[C:2]=[O:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CC=O");
    }
    {
      auto rxn = "[C:2]-[N:1]>>[C:2]~[O:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CCO");
    }
    {
      auto rxn = "[C:2]~[N:1]>>[C:2]=[O:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CCN"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(mol->getNumAtoms() == 3);
      CHECK(MolToSmiles(*mol) == "CC=O");
    }
  }
  SECTION("atom stereo") {
    {  // remove
      auto rxn =
          "[C:1][C@:2]([F:3])([Cl:4])[I:5]>>[C:1][C:2]([F:3])([Cl:4])[I:5]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[C@](F)(Cl)I"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(MolToSmiles(*mol) == "CCC(F)(Cl)I");
    }
    {  // create
      auto rxn =
          "[C:1][C:2]([F:3])([Cl:4])[I:5]>>[C:1][C@@:2]([F:3])([Cl:4])[I:5]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[C@](F)(Cl)I"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(MolToSmiles(*mol) == "CC[C@@](F)(Cl)I");
    }
    {  // create, swap order
      auto rxn =
          "[C:1][C:2]([F:3])([Cl:4])[I:5]>>[C:1][C@@:2]([Cl:4])([F:3])[I:5]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[C@](F)(Cl)I"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(MolToSmiles(*mol) == "CC[C@](F)(Cl)I");
    }

    {  // invert
      auto rxn =
          "[C:1][C@:2]([F:3])([Cl:4])[I:5]>>[C:1][C@@:2]([F:3])([Cl:4])[I:5]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[C@](F)(Cl)I"_smiles;
      REQUIRE(mol);
      CHECK(rxn->runReactant(*mol));
      CHECK(MolToSmiles(*mol) == "CC[C@@](F)(Cl)I");
    }
    {  // preserve, but order swap
      auto rxn =
          "[C:1][C@:2]([F:3])([Cl:4])[I:5]>>[C:1][C@@:2]([Cl:4])([F:3])[I:5]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "CC[C@](F)(Cl)I"_smiles;
      REQUIRE(mol);
      CHECK(!rxn->runReactant(*mol));
      CHECK(MolToSmiles(*mol) == "CC[C@](F)(Cl)I");
    }
  }
  SECTION("check reactant/product count") {
    {
      auto rxn = "[N:1].[O:2]>>[N:1]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "N"_smiles;
      REQUIRE(mol);
      CHECK_THROWS_AS(rxn->runReactant(*mol), ChemicalReactionException);
    }
    {
      auto rxn = "[N:1]>>[N:1].[O:2]"_rxnsmarts;
      REQUIRE(rxn);
      rxn->initReactantMatchers();
      auto mol = "N"_smiles;
      REQUIRE(mol);
      CHECK_THROWS_AS(rxn->runReactant(*mol), ChemicalReactionException);
    }
  }
}

TEST_CASE("Github #4759 Reaction parser fails when CX extensions are present") {
  std::string sma = "[C:1]Br.[C:2]O>>[C:2][C:1] |$Aryl;;;;;Aryl$|";
  SECTION("SMARTS") {
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(sma));
    REQUIRE(rxn);
    // make sure we have a product and that it didn't end up with a name:
    CHECK(rxn->getProducts().size() == 1);
    CHECK(!rxn->getProducts()[0]->hasProp(common_properties::_Name));
    CHECK(rxn->getProducts()[0]->getNumAtoms() == 2);
    CHECK(rxn->getProducts()[0]->getAtomWithIdx(1)->hasProp(
        common_properties::atomLabel));
    CHECK(rxn->getProducts()[0]->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "Aryl");
    CHECK(rxn->getReactants()[0]->getAtomWithIdx(0)->hasProp(
        common_properties::atomLabel));
    CHECK(rxn->getReactants()[0]->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "Aryl");
  }
  SECTION("SMILES") {
    bool useSmiles = true;
    std::unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction(sma, nullptr, useSmiles));
    REQUIRE(rxn);
    CHECK(rxn->getProducts().size() == 1);
    CHECK(!rxn->getProducts()[0]->hasProp(common_properties::_Name));
    CHECK(rxn->getProducts()[0]->getNumAtoms() == 2);
  }
  SECTION("disabling CXSMILES") {
    bool useSmiles = false;
    bool allowCXSMILES = false;
    std::unique_ptr<ChemicalReaction> rxn(
        RxnSmartsToChemicalReaction(sma, nullptr, useSmiles, allowCXSMILES));
    REQUIRE(rxn);
    CHECK(rxn->getProducts().size() == 1);
    CHECK(!rxn->getProducts()[0]->hasProp(common_properties::_Name));
    CHECK(rxn->getProducts()[0]->getNumAtoms() == 2);
    CHECK(!rxn->getProducts()[0]->getAtomWithIdx(1)->hasProp(
        common_properties::atomLabel));
    CHECK(!rxn->getReactants()[0]->getAtomWithIdx(0)->hasProp(
        common_properties::atomLabel));
  }
  SECTION("Ensure we still handle spaces before/after the >>") {
    auto rxn = " [C:1]Br.[C:2]O >> [C:2][C:1] |$Aryl;;;;;Aryl$|"_rxnsmarts;
    // make sure we have a product and that it didn't end up with a name:
    CHECK(rxn->getProducts().size() == 1);
    CHECK(!rxn->getProducts()[0]->hasProp(common_properties::_Name));
    CHECK(rxn->getProducts()[0]->getNumAtoms() == 2);
  }
  SECTION("advanced space removal") {
    // clang-format off
    auto rxn =
        " [C:1]Br  . [C:2]O    >  CCO  > [C:2][C:1] .   [Cl]    |$Aryl;;;;;Aryl$|"_rxnsmarts;
    // clang-format n
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    CHECK(rxn->getAgents().size() ==1);
    // make sure we have a product and that it didn't end up with a name:
    CHECK(rxn->getProducts().size() == 2);
    CHECK(!rxn->getProducts()[0]->hasProp(common_properties::_Name));
    CHECK(rxn->getProducts()[0]->getNumAtoms() == 2);
  }
  SECTION("not a cxsmiles") {
    // clang-format off
    auto rxn =
        "[C:1]Br.[C:2]O>CCO>[C:2][C:1].[Cl]  reaction_name"_rxnsmarts;
    // clang-format n
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    CHECK(rxn->getAgents().size() ==1);
    // make sure we have a product and that it didn't end up with a name:
    CHECK(rxn->getProducts().size() == 2);
    CHECK(!rxn->getProducts()[0]->hasProp(common_properties::_Name));
    CHECK(rxn->getProducts()[0]->getNumAtoms() == 2);
  }
}

TEST_CASE("CXSMILES for reactions", "[cxsmiles]") {
  SECTION("basics") {
    // clang-format off
    auto rxn = "[CH3:1][CH:2]([CH3:3])[*:4].[OH:5][CH2:6][*:7]>>[CH3:1][CH:2]([CH3:3])[CH2:6][OH:5] |$;;;_AP1;;;_AP1;;;;;$|"_rxnsmarts;
    // clang-format on
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    std::string alabel;
    CHECK(rxn->getReactants()[0]->getAtomWithIdx(3)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
    CHECK(rxn->getReactants()[1]->getAtomWithIdx(2)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
  }
  SECTION("basics with agents") {
    // clang-format off
    auto rxn = "[CH3:1][CH:2]([CH3:3])[*:4].[OH:5][CH2:6][*:7]>O=C=O>[CH3:1][CH:2]([CH3:3])[CH2:6][OH:5] |$;;;_AP1;;;_AP1;;;;;;;;$|"_rxnsmarts;
    // clang-format on
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    std::string alabel;
    CHECK(rxn->getReactants()[0]->getAtomWithIdx(3)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
    CHECK(rxn->getReactants()[1]->getAtomWithIdx(2)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
  }
  SECTION("missing products") {
    // clang-format off
    auto rxn="[CH3:1][CH:2]([CH3:3])[*:4].[OH:5][CH2:6][*:7]>> |$;;;_AP1;;;_AP1$|"_rxnsmarts;
    // clang-format on
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    std::string alabel;
    CHECK(rxn->getReactants()[0]->getAtomWithIdx(3)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
    CHECK(rxn->getReactants()[1]->getAtomWithIdx(2)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
  }
  SECTION("coordinate bonds and sgroups") {
    // when initially writing this, coordinate bonds were not properly parsed
    // from SMARTS, so we use SMILES
    // clang-format off
    auto rxn = "[CH3:1][CH:2]([CH3:3])[*:4].[Fe:8][OH:5][CH2:6][*:7]>>[Fe:8][OH:5][CH2:6][CH2:1][CH:2]([CH3:3])[*:4] "
    "|$;;;_AP1;;;;_AP1;;;;;;;_AP1$,C:5.3,9.6,SgD:6:foo:bar::::,SgD:10:bar:baz::::|"_rxnsmiles;
    // clang-format on
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    CHECK(rxn->getProducts().size() == 1);
    std::string alabel;
    CHECK(rxn->getReactants()[0]->getAtomWithIdx(3)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
    CHECK(rxn->getReactants()[1]->getAtomWithIdx(3)->getPropIfPresent(
        common_properties::atomLabel, alabel));
    CHECK(alabel == "_AP1");
    CHECK(getSubstanceGroups(*rxn->getReactants()[0]).empty());
    CHECK(getSubstanceGroups(*rxn->getReactants()[1]).size() == 1);

    const auto p0 = rxn->getProducts()[0];
    CHECK(p0->getAtomWithIdx(6)->getPropIfPresent(common_properties::atomLabel,
                                                  alabel));
    CHECK(alabel == "_AP1");
    CHECK(getSubstanceGroups(*p0).size() == 1);
  }
  SECTION("sgroup hierarchy") {
    // clang-format off
    auto rxn = "[CH3:6][O:5][CH:3](-*)[O:2]-*>>[CH3:6][NH:5][CH:3](-*)[O:2]-* "
    "|$;;;star_e;;star_e;;;;star_e;;star_e$,SgD:1,0:foo:bar::::,SgD:7,6:foo:baz::::,Sg:n:4,2,1,0::ht,Sg:n:10,8,7,6::ht,SgH:2:0,3:1|"_rxnsmiles;
    // clang-format on
    REQUIRE(rxn);
    CHECK(getSubstanceGroups(*rxn->getReactants()[0]).size() == 2);
    CHECK(getSubstanceGroups(*rxn->getProducts()[0]).size() == 2);
  }
  SECTION("link nodes") {
    // clang-format off
    auto rxn = "CO.OC1CCC(F)C1>>COC1CC(O)CC1F |LN:3:1.3.4.8,13:2.5.12.15|"_rxnsmarts;
    // clang-format on
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    CHECK(rxn->getProducts().size() == 1);
    CHECK(
        !rxn->getReactants()[0]->hasProp(common_properties::molFileLinkNodes));
    std::string lns;
    CHECK(rxn->getReactants()[1]->getPropIfPresent(
        common_properties::molFileLinkNodes, lns));
    CHECK(lns == "1 3 2 2 3 2 7");
    CHECK(rxn->getProducts()[0]->getPropIfPresent(
        common_properties::molFileLinkNodes, lns));
    CHECK(lns == "2 5 2 5 4 5 7");
  }
#if 1
  // note that these only work with the current parser if the
  // variable-attachment point part is grouped with the molecule it's attached
  // to. This probably isn't the end of the world
  SECTION("variable attachment points") {
    // clang-format off
    auto rxn = "CN.(CO*.CC1=CN=CC=C1)>>(CNC1=C(C)C=CC=N1.CO*) |m:4:11.9.10,23:17.19.18|"_rxnsmarts;
    // clang-format on
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 2);
    CHECK(rxn->getProducts().size() == 1);
    auto bnd = rxn->getReactants()[1]->getBondBetweenAtoms(1, 2);
    REQUIRE(bnd);
    CHECK(bnd->hasProp(common_properties::_MolFileBondAttach));
    CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondAttach) ==
          "ANY");
    CHECK(bnd->hasProp(common_properties::_MolFileBondEndPts));
    CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondEndPts) ==
          "(3 10 8 9)");

    bnd = rxn->getProducts()[0]->getBondBetweenAtoms(10, 11);
    REQUIRE(bnd);
    CHECK(bnd->hasProp(common_properties::_MolFileBondAttach));
    CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondAttach) ==
          "ANY");
    CHECK(bnd->hasProp(common_properties::_MolFileBondEndPts));
    CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondEndPts) ==
          "(3 6 8 7)");
  }
#endif
  SECTION("cis/trans markers") {
    auto rxn =
        "C1C=CC=CCC=CC=C1>>C1C=CC=CNC=CC=C1 |c:1,3,6,8,11,13,16,18|"_rxnsmiles;
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 1);
    CHECK(rxn->getProducts().size() == 1);
    CHECK(rxn->getReactants()[0]->getBondWithIdx(1)->getStereo() ==
          Bond::BondStereo::STEREOCIS);
    CHECK(rxn->getProducts()[0]->getBondWithIdx(1)->getStereo() ==
          Bond::BondStereo::STEREOCIS);
  }
  SECTION("wedged bonds") {
    auto rxn = "CC(O)(F)Cl>>CC(N)(F)Cl |w:1.0,6.5|"_rxnsmiles;
    REQUIRE(rxn);
    CHECK(rxn->getReactants().size() == 1);
    CHECK(rxn->getProducts().size() == 1);
    unsigned int bondcfg;
    CHECK(rxn->getReactants()[0]->getBondWithIdx(0)->getPropIfPresent(
        "_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 2);
    CHECK(rxn->getProducts()[0]->getBondWithIdx(1)->getPropIfPresent(
        "_MolFileBondCfg", bondcfg));
    CHECK(bondcfg == 2);
  }
}

TEST_CASE("V3K rxn blocks") {
  SECTION("writing basics") {
    // clang-format off
    auto rxn =
        "[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[Br].[#6:7]B(O)O>>[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[#6:7]"_rxnsmarts;
    // clang-format off
    REQUIRE(rxn);
    auto rxnb = ChemicalReactionToV3KRxnBlock(*rxn);
    bool separateAgents=false;
    bool forceV3000=true;
    auto rxnb2 = ChemicalReactionToRxnBlock(*rxn,separateAgents,forceV3000);
    CHECK(rxnb==rxnb2);

    std::unique_ptr<ChemicalReaction> rxn2{RxnBlockToChemicalReaction(rxnb)};
    REQUIRE(rxn2);
    CHECK(rxn->getNumAgentTemplates()==rxn2->getNumAgentTemplates());
    CHECK(rxn->getNumReactantTemplates()==rxn2->getNumReactantTemplates());
    CHECK(rxn->getNumProductTemplates()==rxn2->getNumProductTemplates());   
  }

  SECTION("github5324") {
      // Test sgroup in a ring - this example failed with improperr tail crossings
      auto mol = "C-1-C-C-C-C-O-1 |Sg:n:4:n:ht|"_smarts;
      MolOps::findSSSR(*mol);
      auto mbk = FileParserUtils::getV3000CTAB(*mol, -1);
      CHECK(mbk.find("ATOMS=(1 5) XBONDS=(2 4 5) XBHEAD=(1 4) XBCORR=(2 4 5)")!=std::string::npos);
      std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
                                  ">>C-1-C-C-C-C-O-1 |Sg:n:4:n:ht|"));
      auto rxnb = ChemicalReactionToV3KRxnBlock(*rxn);
      CHECK(rxnb.find("ATOMS=(1 5) XBONDS=(2 4 5) XBHEAD=(1 4) XBCORR=(2 4 5)")!=std::string::npos);
    }
}

TEST_CASE("CDXML Parser") {
  std::string cdxmlbase = std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("CDXML REACTION") {
      auto fname = cdxmlbase + "rxn2.cdxml";
      std::vector<std::string> expected = {"Cl[c:1]1[cH:4][cH:3][cH:2][cH:6][cH:5]1",
           "OC(O)B[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
           "[CH:1]1=[CH:5][C:6]([C:7]2=[CH:12][CH:11]=[CH:10][CH:9]=[CH:8]2)=[CH:2][CH:3]=[CH:4]1"};
      
       auto rxns = CDXMLFileToChemicalReactions(fname);
       CHECK(rxns.size() == 1);
       int i=0;
       int count = 0;
       for(auto &mol : rxns[0]->getReactants()) {
           CHECK(mol->getProp<unsigned int>("CDX_SCHEME_ID") == 397);
           CHECK(mol->getProp<unsigned int>("CDX_STEP_ID") == 398);
           CHECK(mol->getProp<unsigned int>("CDX_REAGENT_ID") == i++);
           CHECK(MolToSmiles(*mol) == expected[count++]);
       }
       i = 0;
       for(auto &mol : rxns[0]->getProducts()) {
           CHECK(mol->getProp<unsigned int>("CDX_SCHEME_ID") == 397);
           CHECK(mol->getProp<unsigned int>("CDX_STEP_ID") == 398);
           CHECK(mol->getProp<unsigned int>("CDX_PRODUCT_ID") == i++);
           CHECK(MolToSmiles(*mol) == expected[count++]);
       }
   
       auto smarts = ChemicalReactionToRxnSmarts(*rxns[0]);
       CHECK(smarts == "[#6&D2:2]1:[#6&D2:3]:[#6&D2:4]:[#6&D3:1](:[#6&D2:5]:[#6&D2:6]:1)-[#17&D1].[#6&D3](-[#5&D2]-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1)(-[#8&D1])-[#8&D1]>>[#6:1]1=[#6:5]-[#6:6](=[#6:2]-[#6:3]=[#6:4]-1)-[#6:7]1-[#6:8]=[#6:9]-[#6:10]=[#6:11]-[#6:12]=1");
  }
}

TEST_CASE("Github #5785: separateAgents ignored for V3000 RXN files") {
  SECTION("general separateAgents parse testing: V2000"){
    std::string rxnb = R"RXN($RXN

      RDKit

  1  1  1
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0
M  END
$MOL

     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
    1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
  1  2  2  0
M  END
$MOL

     RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
M  END
)RXN";
    std::unique_ptr<ChemicalReaction> rxn(RxnBlockToChemicalReaction(rxnb));
    REQUIRE(rxn);
    CHECK(rxn->getNumReactantTemplates()==1);
    CHECK(rxn->getNumProductTemplates()==1);
    CHECK(rxn->getNumAgentTemplates()==1);

    auto orxn = ChemicalReactionToRxnBlock(*rxn,true);
    CHECK(orxn.find("  1  1  1") != std::string::npos);
    
  }
  SECTION("general separateAgents parse testing: V3000"){
    std::string rxnb = R"RXN($RXN V3000

      Mrv2211  121520220816

M  V30 COUNTS 1 1 1
M  V30 BEGIN REACTANT
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.8977 -0.385 0 1
M  V30 2 O -3.564 0.385 0 2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  V30 END REACTANT
M  V30 BEGIN PRODUCT
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 4.444 -0.385 0 1
M  V30 2 O 5.7777 0.385 0 2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  V30 END PRODUCT
M  V30 BEGIN AGENT
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Pt 0 1.54 0 0
M  V30 END ATOM
M  V30 END CTAB
M  V30 END AGENT
M  END
)RXN";
std::unique_ptr<ChemicalReaction> rxn(RxnBlockToChemicalReaction(rxnb));
    REQUIRE(rxn);
    CHECK(rxn->getNumReactantTemplates()==1);
    CHECK(rxn->getNumProductTemplates()==1);
    CHECK(rxn->getNumAgentTemplates()==1);
   
    {    // with separate agents
      auto orxn = ChemicalReactionToV3KRxnBlock(*rxn,true);
      CHECK(orxn.find("COUNTS 1 1 1") != std::string::npos);
      CHECK(orxn.find("BEGIN AGENT") != std::string::npos);
      CHECK(orxn.find("END AGENT") != std::string::npos);
    }
    {    // without separate agents
      auto orxn = ChemicalReactionToV3KRxnBlock(*rxn,false);
      CHECK(orxn.find("COUNTS 2 1") != std::string::npos);
      CHECK(orxn.find("BEGIN AGENT") == std::string::npos);
      CHECK(orxn.find("END AGENT") == std::string::npos);
    }
  }
}