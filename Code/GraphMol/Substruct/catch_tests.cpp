//
//  Copyright (C) 2019-2023 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Tests of substructure searching
//

#include "catch.hpp"

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;
typedef std::tuple<std::string, std::string, size_t> matchCase;

class _IsSubstructOf : public Catch::MatcherBase<const ROMol &> {
  ROMol const *m_mol;
  SubstructMatchParameters m_ps;

 public:
  _IsSubstructOf(const ROMol &m) : m_mol(&m) {}

  _IsSubstructOf(const ROMol &m, SubstructMatchParameters ps)
      : m_mol(&m), m_ps(std::move(ps)) {}

  bool match(const ROMol &query) const override {
    return !SubstructMatch(*m_mol, query, m_ps).empty();
  }

  std::string describe() const override {
    std::ostringstream ss;
    ss << "is not a substructure of " << MolToCXSmiles(*m_mol);
    return ss.str();
  }
};

static _IsSubstructOf IsSubstructOf(const ROMol &m,
                                    const SubstructMatchParameters &ps) {
  return _IsSubstructOf(m, ps);
}

static _IsSubstructOf IsSubstructOf(const ROMol &m) {
  return _IsSubstructOf(m);
}

namespace Catch {
// ""_smiles returns an RWMol.
template <>
struct StringMaker<RDKit::RWMol> {
  static std::string convert(RDKit::RWMol const &m) { return MolToCXSmiles(m); }
};
}  // namespace Catch

TEST_CASE("substructure parameters", "[substruct]") {
  SECTION("chirality") {
    auto mol1 = "CCC[C@@H]1CN(CCC)CCN1"_smiles;
    auto mol2 = "CCC[C@H]1CN(CCC)CCN1"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);

    SubstructMatchParameters ps;
    // default is to ignore chirality:
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 1);
    CHECK(SubstructMatch(*mol1, *mol1, ps).size() == 1);

    ps.useChirality = true;
    CHECK_THAT(*mol2, !IsSubstructOf(*mol1, ps));
    CHECK(SubstructMatch(*mol1, *mol1, ps).size() == 1);
  }
  SECTION("conjugated matching aromaticity 1") {
    auto mol1 = "C1=COC=C1"_smiles;
    REQUIRE(mol1);
    RWMol mol2(*mol1);
    MolOps::Kekulize(mol2);
    SubstructMatchParameters ps;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 0);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 0);

    ps.aromaticMatchesConjugated = true;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 1);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 1);
  }
  SECTION("conjugated matching aromaticity 2") {
    auto mol1 = "c1ccccc1"_smiles;
    REQUIRE(mol1);
    RWMol mol2(*mol1);
    MolOps::Kekulize(mol2);
    SubstructMatchParameters ps;
    CHECK_THAT(mol2, !IsSubstructOf(*mol1));
    CHECK_THAT(*mol1, !IsSubstructOf(mol2));

    ps.aromaticMatchesConjugated = true;
    CHECK(SubstructMatch(*mol1, mol2, ps).size() == 1);
    CHECK(SubstructMatch(mol2, *mol1, ps).size() == 1);
  }

  SECTION("conjugated matching aromaticity bulk") {
    std::vector<matchCase> examples;
    examples.push_back(
        std::make_tuple(std::string("c1ccccc1"), std::string("C1CCCCC1"), 0));
    examples.push_back(
        std::make_tuple(std::string("C1CCCCC1"), std::string("c1ccccc1"), 0));
    examples.push_back(std::make_tuple(std::string("O=C1C=CC(=O)C=C1"),
                                       std::string("c1ccccc1"), 1));
    SubstructMatchParameters ps;
    ps.aromaticMatchesConjugated = true;
    for (const auto &example : examples) {
      // std::cerr << "   " << std::get<0>(example) << " - "
      //           << std::get<1>(example) << std::endl;
      std::unique_ptr<RWMol> m1(SmilesToMol(std::get<0>(example)));
      REQUIRE(m1);
      std::unique_ptr<RWMol> m2(SmilesToMol(std::get<1>(example)));
      CHECK(SubstructMatch(*m1, *m2, ps).size() == std::get<2>(example));
    }
  }
  SECTION("looping") {
    auto mol1 = "CC(=O)C(=O)C(=O)"_smiles;
    auto mol2 = "C=O"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    for (auto match : SubstructMatch(*mol1, *mol2)) {
      CHECK(match.size() == 2);
    }
  }
}

namespace {
bool no_match(const ROMol &mol, const std::vector<unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(ids);
  return false;
}
bool always_match(const ROMol &mol, const std::vector<unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(ids);
  return true;
}
bool bigger(const ROMol &mol, const std::vector<unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  return std::accumulate(ids.begin(), ids.end(), 0) > 5;
}
}  // namespace
TEST_CASE("providing a final match function", "[substruct]") {
  SECTION("basics") {
    auto mol1 = "CCOC"_smiles;
    auto mol2 = "CCO"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    SubstructMatchParameters ps;
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 1);
    ps.extraFinalCheck = &no_match;
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 0);
    ps.extraFinalCheck = &always_match;
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 1);
  }
  SECTION("test 2") {
    auto mol1 = "CCOCC"_smiles;
    auto mol2 = "CCO"_smiles;
    REQUIRE(mol1);
    REQUIRE(mol2);
    SubstructMatchParameters ps;
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 2);
    ps.extraFinalCheck = &bigger;
    CHECK(SubstructMatch(*mol1, *mol2, ps).size() == 1);
  }
}

TEST_CASE("Enhanced stereochemistry", "[substruct][StereoGroup]") {
  // Chirality specifications.
  // 1. An achiral molecule: CC(O)C(CC)F means unknown/all stereoisomers
  // 2. A chiral molecule: C[C@H](O)[C@H](CC)F means 1 stereoisomer
  // 3. A chiral molecule with an AND specifier: C[C@H](O)[C@H](CC)F |a1:1,3|
  // means both stereoisomers
  // 4. A chiral molecule with an OR specifier: C[C@H](O)[C@H](CC)F |o1:1,3|
  // means one of the two stereoisomers
  auto mol_achiral = "CC(O)C(CC)F"_smiles;
  auto mol_chiral = "C[C@H](O)[C@H](CC)F"_smiles;
  auto mol_and = "C[C@H](O)[C@H](CC)F |&1:1,3|"_smiles;
  auto mol_or = "C[C@H](O)[C@H](CC)F |o1:1,3|"_smiles;
  auto mol_absolute = "C[C@H](O)[C@H](CC)F |a:1,3|"_smiles;
  auto diastereomer = "C[C@H](O)[C@@H](CC)F"_smiles;

  SubstructMatchParameters ps;
  ps.useChirality = true;
  ps.useEnhancedStereo = true;

  SECTION("achiral search matches anything") {
    CHECK_THAT(*mol_achiral, IsSubstructOf(*mol_chiral, ps));
    CHECK_THAT(*mol_achiral, IsSubstructOf(*mol_and, ps));
    CHECK_THAT(*mol_achiral, IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_achiral, IsSubstructOf(*mol_absolute, ps));
    CHECK_THAT(*mol_achiral, IsSubstructOf(*diastereomer, ps));
  }
  SECTION("chiral molecule is a substructure of AND or OR") {
    CHECK_THAT(*mol_chiral, !IsSubstructOf(*mol_achiral, ps));
    CHECK_THAT(*mol_chiral, IsSubstructOf(*mol_and, ps));
    CHECK_THAT(*mol_chiral, IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_chiral, !IsSubstructOf(*diastereomer, ps));
    CHECK_THAT(*mol_absolute, !IsSubstructOf(*mol_achiral, ps));
    CHECK_THAT(*mol_absolute, IsSubstructOf(*mol_and, ps));
    CHECK_THAT(*mol_absolute, IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_absolute, !IsSubstructOf(*diastereomer, ps));
  }
  SECTION("AND query only matches AND") {
    // because it means BOTH, and only AND includes both.
    CHECK_THAT(*mol_and, !IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_and, IsSubstructOf(*mol_and, ps));
    CHECK_THAT(*mol_and, !IsSubstructOf(*mol_absolute, ps));
    CHECK_THAT(*mol_and, !IsSubstructOf(*mol_chiral, ps));
    CHECK_THAT(*mol_and, !IsSubstructOf(*mol_achiral, ps));
  }
  SECTION("An OR query matches AND and OR") {
    // because AND is both, so it's a superset of the molecules described in
    // the OR
    CHECK_THAT(*mol_or, !IsSubstructOf(*mol_chiral, ps));
    CHECK_THAT(*mol_or, !IsSubstructOf(*mol_absolute, ps));
    CHECK_THAT(*mol_or, !IsSubstructOf(*diastereomer, ps));
    CHECK_THAT(*mol_or, IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_or, IsSubstructOf(*mol_and, ps));
  }
  SECTION("AND and OR match their enantiomer") {
    // This is, like, the point of And/Or
    auto enantiomer = "C[C@@H](O)[C@@H](CC)F"_smiles;
    CHECK_THAT(*enantiomer, IsSubstructOf(*mol_and, ps));
    CHECK_THAT(*enantiomer, IsSubstructOf(*mol_or, ps));
  }
  SECTION("But not some arbitrary diastereomer") {
    CHECK_THAT(*diastereomer, !IsSubstructOf(*mol_and, ps));
    CHECK_THAT(*diastereomer, !IsSubstructOf(*mol_or, ps));
  }
  SECTION("Mixed stereo groups include single stereo groups") {
    auto mol_mixed_or = "C[C@H](O)[C@H](CC)F |o1:1,o2:3|"_smiles;
    CHECK_THAT(*mol_mixed_or, !IsSubstructOf(*mol_or, ps));
    // OR refers to two of the 4 molecules that mol_mixed_or
    CHECK_THAT(*mol_or, IsSubstructOf(*mol_mixed_or, ps));

    auto mol_mixed_or2 = "C[C@H](O)[C@@H](CC)F |o1:1,o2:3|"_smiles;
    CHECK_THAT(*mol_mixed_or2, !IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_or, IsSubstructOf(*mol_mixed_or2, ps));

    // I'm not sure about these ones, but they should be symmetric:
    auto mol_mixed_or_and_abs = "C[C@H](O)[C@H](CC)F |o1:1|"_smiles;
    CHECK_THAT(*mol_mixed_or_and_abs, !IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_or, !IsSubstructOf(*mol_mixed_or_and_abs, ps));

    auto mol_mixed_or_and_abs2 = "C[C@@H](O)[C@H](CC)F |o1:1|"_smiles;
    CHECK_THAT(*mol_mixed_or_and_abs2, !IsSubstructOf(*mol_or, ps));
    CHECK_THAT(*mol_or, !IsSubstructOf(*mol_mixed_or_and_abs, ps));
  }
  SECTION("It's OK to match part of a stereo group, though") {
    auto mol_and_long = "F[C@@H](O)C[C@@H](CC)F |&1:1,3|"_smiles;
    auto mol_and_partial = "F[C@@H](O)C |&1:1|"_smiles;
    auto mol_or_long = "F[C@@H](O)C[C@@H](CC)F |o1:1,3|"_smiles;
    auto mol_or_partial = "F[C@@H](O)C |o1:1|"_smiles;

    CHECK_THAT(*mol_and_partial, IsSubstructOf(*mol_and_long, ps));
    CHECK_THAT(*mol_or_partial, IsSubstructOf(*mol_or_long, ps));
    CHECK_THAT(*mol_or_partial, IsSubstructOf(*mol_and_long, ps));
    CHECK_THAT(*mol_and_partial, !IsSubstructOf(*mol_or_long, ps));
  }
}

TEST_CASE("Github #4138: empty query produces non-empty results",
          "[substruct][bug]") {
  auto mol = "C1CCCCO1"_smiles;
  auto emol = ""_smiles;
  auto qry = "C"_smarts;
  auto eqry = ""_smarts;
  REQUIRE(mol);
  REQUIRE(qry);
  SECTION("empty query") {
    {
      auto matches = SubstructMatch(*mol, *eqry);
      CHECK(matches.empty());
    }
    {
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*mol, *eqry, matches));
      CHECK(matches.empty());
    }
    {
      MatchVectType match;
      CHECK(!SubstructMatch(*mol, *eqry, match));
      CHECK(match.empty());
    }
  }
  SECTION("empty mol") {
    {
      auto matches = SubstructMatch(*emol, *qry);
      CHECK(matches.empty());
    }
    {
      std::vector<MatchVectType> matches;
      CHECK(!SubstructMatch(*emol, *qry, matches));
      CHECK(matches.empty());
    }
    {
      MatchVectType match;
      CHECK(!SubstructMatch(*emol, *qry, match));
      CHECK(match.empty());
    }
  }
}

TEST_CASE("Github #4558: GetSubstructMatches() loops at 43690 iterations",
          "[substruct][bug]") {
  // We need LOTS of water molecules here.
  auto num_mols = 22000u;
  std::stringstream smi;
  for (auto i = 1u; i < num_mols; ++i) {
    smi << "[H]O[H].";
  }
  smi << "[H]O[H]";  // last one (notice we started at 1)

  int debug = 0;
  bool sanitize = false;  // don't sanitize, it takes too long.
  std::unique_ptr<ROMol> mol(SmilesToMol(smi.str(), debug, sanitize));
  REQUIRE(mol);

  auto qry = "[H]O[H]"_smarts;

  SubstructMatchParameters ps;
  ps.uniquify = false;           // don't uniquify, it takes too long.
  ps.maxMatches = 3 * num_mols;  // exceed the numer of matches we expect

  auto matches = SubstructMatch(*mol, *qry, ps);
  CHECK(matches.size() == num_mols * 2);
}

TEST_CASE(
    "Github #888: GetSubstructMatches uniquify and maxMatches don't work well together ") {
  SECTION("Basics") {
    auto m = "CCCCCC"_smiles;
    auto q = "CC"_smarts;
    REQUIRE(m);
    REQUIRE(q);
    SubstructMatchParameters ps;
    ps.uniquify = false;
    ps.maxMatches = 4;
    auto matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 4);

    ps.uniquify = true;
    matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 4);

    ps.useChirality = true;
    matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 4);
  }
  SECTION("interaction with chirality") {
    auto m = "C/C=C/C=C/C=C\\C=C/C=C/C"_smiles;
    auto q = "C/C=C/C"_smarts;
    REQUIRE(m);
    REQUIRE(q);
    SubstructMatchParameters ps;
    ps.uniquify = false;
    ps.maxMatches = 2;
    auto matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 2);

    ps.uniquify = true;
    matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 2);

    ps.useChirality = true;
    matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 2);

    ps.maxMatches = 1000;
    matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 3);
  }
  SECTION("interactions with recursive SMARTS") {
    {
      auto m = "N#CC#N"_smiles;
      REQUIRE(m);
      auto q = "[!$(*#*)]"_smarts;
      REQUIRE(q);
      SubstructMatchParameters ps;
      ps.uniquify = true;
      ps.maxMatches = 1;
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.empty());
    }
    {
      auto m = "N#CC#N"_smiles;
      REQUIRE(m);
      auto q = "[$(*#*)&!D1]"_smarts;
      REQUIRE(q);
      SubstructMatchParameters ps;
      ps.uniquify = true;
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 2);
    }
    {
      auto m = "N#CCC#N"_smiles;
      REQUIRE(m);
      auto q = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"_smarts;
      REQUIRE(q);
      SubstructMatchParameters ps;
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.empty());
    }
  }
}
