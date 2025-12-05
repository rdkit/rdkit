//
//  Copyright (C) 2019-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Tests of substructure searching
//

#include <catch2/catch_all.hpp>

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/MolPickler.h>

using namespace RDKit;
typedef std::tuple<std::string, std::string, size_t> matchCase;

class _IsSubstructOf : public Catch::Matchers::MatcherBase<const ROMol &> {
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

  SECTION("atom properties") {
    std::vector<matchCase> examples;
    examples.push_back(
        std::make_tuple(std::string("CCCCCCCCC"), std::string("CCC"), 7));
    examples.push_back(
        std::make_tuple(std::string("CCCCCCCCC |atomProp:0.test_prop.1|"),
                        std::string("CCC |atomProp:0.test_prop.1|"), 1));
    examples.push_back(
        std::make_tuple(std::string("CCCCCCCCC |atomProp:0.test_prop.1|"),
                        std::string("CCC"), 6));
    examples.push_back(
        std::make_tuple(std::string("CCCCCCCCC"),
                        std::string("CCC |atomProp:0.test_prop.1|"), 0));
    examples.push_back(
        std::make_tuple(std::string("CCCCCCCCC |atomProp:0.test_prop.1|"),
                        std::string("CCC |atomProp:0.test_prop.2|"), 0));
    SubstructMatchParameters ps;
    ps.atomProperties = {"test_prop"};
    for (const auto &example : examples) {
      std::unique_ptr<RWMol> m1(SmilesToMol(std::get<0>(example)));
      REQUIRE(m1);
      std::unique_ptr<RWMol> m2(SmilesToMol(std::get<1>(example)));
      REQUIRE(m2);
      CHECK(SubstructMatch(*m1, *m2, ps).size() == std::get<2>(example));
    }
  }

  SECTION("bond properties") {
    std::unique_ptr<RWMol> m(SmilesToMol("CCCCCCCCC"));
    std::unique_ptr<RWMol> m_with_prop(SmilesToMol("CCCCCCCCC"));
    m_with_prop->getBondWithIdx(0)->setProp("test_prop", "1");

    std::unique_ptr<RWMol> q(SmilesToMol("CCC"));
    std::unique_ptr<RWMol> q_with_prop(SmilesToMol("CCC"));
    std::unique_ptr<RWMol> q_with_prop2(SmilesToMol("CCC"));
    q_with_prop->getBondWithIdx(0)->setProp("test_prop", "1");
    q_with_prop2->getBondWithIdx(0)->setProp("test_prop", "2");

    SubstructMatchParameters ps;
    ps.bondProperties = {"test_prop"};
    CHECK(SubstructMatch(*m, *q, ps).size() == 7);
    CHECK(SubstructMatch(*m_with_prop, *q_with_prop, ps).size() == 1);
    CHECK(SubstructMatch(*m_with_prop, *q_with_prop2, ps).size() == 0);
    CHECK(SubstructMatch(*m_with_prop, *q, ps).size() == 6);
    CHECK(SubstructMatch(*m, *q_with_prop, ps).size() == 0);

    // now check with bond and atom properties
    m_with_prop->getAtomWithIdx(0)->setProp("test_prop", "1");
    q_with_prop->getAtomWithIdx(0)->setProp("test_prop", "1");

    ps.atomProperties = {"test_prop"};
    CHECK(SubstructMatch(*m_with_prop, *q_with_prop, ps).size() == 1);
    CHECK(SubstructMatch(*m_with_prop, *q_with_prop2, ps).size() == 0);
    CHECK(SubstructMatch(*m_with_prop, *q, ps).size() == 6);
    CHECK(SubstructMatch(*m, *q_with_prop, ps).size() == 0);

    // Currently, a property set as an int will match to a property
    // set as a different type if they cast to the same value
    // TODO: Ensure property types are the same in substructure matching
    q_with_prop->getBondWithIdx(0)->clearProp("test_prop");
    q_with_prop->getBondWithIdx(0)->setProp<int>("test_prop", 1);
    CHECK(SubstructMatch(*m_with_prop, *q_with_prop, ps).size() == 1);
  }
}

namespace {
bool no_match(const ROMol &mol, const std::span<const unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(ids);
  return false;
}
bool always_match(const ROMol &mol, const std::span<const unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(ids);
  return true;
}
bool bigger(const ROMol &mol, const std::span<const unsigned int> &ids) {
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

TEST_CASE("Github #6017: add maxRecursiveMatches to SubstructMatchParameters") {
  SECTION("Basics") {
    auto m = "OCC(O)C(O)C(O)C(O)CO"_smiles;
    auto q = "[$(CO)][$(CO)]"_smarts;
    REQUIRE(m);
    REQUIRE(q);
    SubstructMatchParameters ps;
    ps.uniquify = true;
    {
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 5);
    }

    // if maxRecursiveMatches isn't larger than maxMatches this will fail
    ps.maxMatches = 3;
    {
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 3);
    }
  }
  SECTION("maxMatches larger than maxRecursiveMatches") {
    auto m = "OCC(O)C(O)C(O)C(O)CO"_smiles;
    auto q = "[$(CO)]C"_smarts;
    REQUIRE(m);
    REQUIRE(q);
    SubstructMatchParameters ps;
    ps.uniquify = true;
    ps.maxMatches = 3;
    ps.maxRecursiveMatches = 2;
    {
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 3);
    }
  }
}

TEST_CASE(
    "GitHub Issue #6983: SubstructMatch maxRecursiveMatches is not being honored",
    "[bug][substruct]") {
  constexpr unsigned num_atoms = 1005;
  std::string smiles;
  smiles.reserve(num_atoms * 2);

  // 'smiles' already contains 1 O, so start from 1
  // so we end up with 'num_atoms' water mols
  smiles += "O";
  for (unsigned i = 1; i < num_atoms; ++i) {
    smiles += ".O";
  }
  std::unique_ptr<RWMol> m(SmilesToMol(smiles));
  auto q = "[$(O)]"_smarts;
  REQUIRE(m);
  REQUIRE(q);

  SubstructMatchParameters ps;
  ps.maxMatches = num_atoms * 2;
  ps.maxRecursiveMatches = ps.maxMatches;
  {
    auto matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == num_atoms);
  }
}

TEST_CASE("pickling HasPropWithValue queries") {
  SubstructMatchParameters ps;
  SECTION("basics int") {
    auto mol = "CC"_smarts;
    auto target = "CC"_smiles;
    REQUIRE(mol);
    REQUIRE(target);

    RWMol mol2(*mol);
    RWMol mol3(*mol);

    mol->getAtomWithIdx(0)->expandQuery(makePropQuery<Atom, int>("foo", 1, 1));
    mol->getBondWithIdx(0)->expandQuery(makePropQuery<Bond, int>("bar", 1, 0));

    mol2.getAtomWithIdx(0)->expandQuery(makePropQuery<Atom, int>("foo", 1, 0));
    mol2.getBondWithIdx(0)->expandQuery(makePropQuery<Bond, int>("bar", 1, 0));
    mol3.getAtomWithIdx(0)->expandQuery(makePropQuery<Atom, int>("foo", 2, 0));
    mol3.getBondWithIdx(0)->expandQuery(makePropQuery<Bond, int>("bar", 2, 0));

    CHECK(SubstructMatch(*target, *mol, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol3, ps).size() == 0);
    target->getAtomWithIdx(0)->setProp<int>("foo", 2);
    target->getBondWithIdx(0)->setProp<int>("bar", 1);
    CHECK(SubstructMatch(*target, *mol, ps).size() == 1);
    CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol3, ps).size() == 0);

    {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 1);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, *mol, ps).size() == 1);
    }
    {
      std::string pkl;
      MolPickler::pickleMol(mol2, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 0);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    }
    {
      std::string pkl;
      MolPickler::pickleMol(mol3, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 0);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, mol3, ps).size() == 0);
    }
  }
  SECTION("basics string") {
    auto mol = "CC"_smarts;
    auto target = "CC"_smiles;
    REQUIRE(mol);
    REQUIRE(target);

    RWMol mol2(*mol);
    RWMol mol3(*mol);

    mol->getAtomWithIdx(0)->expandQuery(
        makePropQuery<Atom, std::string>("foo", "asdfs"));
    mol->getBondWithIdx(0)->expandQuery(
        makePropQuery<Bond, std::string>("bar", "dsafasdf"));

    mol2.getAtomWithIdx(0)->expandQuery(
        makePropQuery<Atom, std::string>("foo", "asdfs"));
    mol2.getBondWithIdx(0)->expandQuery(
        makePropQuery<Bond, std::string>("bar", "dsa"));
    mol3.getAtomWithIdx(0)->expandQuery(
        makePropQuery<Atom, std::string>("foo", "asd"));
    mol3.getBondWithIdx(0)->expandQuery(
        makePropQuery<Bond, std::string>("bar", "dsafasdf"));

    CHECK(SubstructMatch(*target, *mol, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol3, ps).size() == 0);
    target->getAtomWithIdx(0)->setProp<std::string>("foo", "asdfs");
    target->getBondWithIdx(0)->setProp<std::string>("bar", "dsafasdf");
    CHECK(SubstructMatch(*target, *mol, ps).size() == 1);
    CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol3, ps).size() == 0);

    {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getBondWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 1);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, *mol, ps).size() == 1);
    }
    {
      std::string pkl;
      MolPickler::pickleMol(mol2, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getBondWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 0);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    }
    {
      std::string pkl;
      MolPickler::pickleMol(mol3, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getBondWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 0);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, mol3, ps).size() == 0);
    }
  }
  SECTION("basics EBV") {
    auto mol = "CC"_smarts;
    auto target = "CC"_smiles;
    REQUIRE(mol);
    REQUIRE(target);

    ExplicitBitVect bv(10);
    bv.setBit(0);
    bv.setBit(2);
    bv.setBit(3);
    bv.setBit(5);
    bv.setBit(7);
    bv.setBit(9);

    RWMol mol2(*mol);
    RWMol mol3(*mol);

    mol->getAtomWithIdx(0)->expandQuery(
        makePropQuery<Atom, ExplicitBitVect>("foo", bv, 0.0));
    mol->getBondWithIdx(0)->expandQuery(
        makePropQuery<Bond, ExplicitBitVect>("bar", bv, 0.0));
    ExplicitBitVect bv2(bv);
    bv2.unsetBit(9);
    mol2.getAtomWithIdx(0)->expandQuery(
        makePropQuery<Atom, ExplicitBitVect>("foo", bv2, 0.0));
    mol3.getBondWithIdx(0)->expandQuery(
        makePropQuery<Bond, ExplicitBitVect>("bar", bv2, 0.0));

    CHECK(SubstructMatch(*target, *mol, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol3, ps).size() == 0);

    target->getAtomWithIdx(0)->setProp<ExplicitBitVect>("foo", bv);
    target->getBondWithIdx(0)->setProp<ExplicitBitVect>("bar", bv);
    CHECK(SubstructMatch(*target, *mol, ps).size() == 1);
    CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    CHECK(SubstructMatch(*target, mol3, ps).size() == 0);
    {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getBondWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 1);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, *mol, ps).size() == 1);
    }
    {
      std::string pkl;
      MolPickler::pickleMol(mol2, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getBondWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 0);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, mol2, ps).size() == 0);
    }
    {
      std::string pkl;
      MolPickler::pickleMol(mol3, pkl);
      RWMol pklmol(pkl);
      REQUIRE(pklmol.getAtomWithIdx(0)->hasQuery());
      REQUIRE(pklmol.getBondWithIdx(0)->hasQuery());
      CHECK(SubstructMatch(*target, pklmol, ps).size() == 0);
      // make sure we are idempotent in pickling
      CHECK(SubstructMatch(*target, mol3, ps).size() == 0);
    }
  }
}

TEST_CASE("specified query matches unspecified atom") {
  SECTION("atom basics") {
    auto q = "F[C@](Cl)(Br)C"_smarts;
    REQUIRE(q);

    auto m1 = "F[C@](Cl)(Br)C"_smiles;
    REQUIRE(m1);
    auto m2 = "FC(Cl)(Br)C"_smiles;
    REQUIRE(m2);
    auto m3 = "F[C@@](Cl)(Br)C"_smiles;
    REQUIRE(m3);

    SubstructMatchParameters ps;
    ps.useChirality = true;
    CHECK(SubstructMatch(*m1, *q, ps).size() == 1);
    CHECK(SubstructMatch(*m2, *q, ps).empty());
    CHECK(SubstructMatch(*m3, *q, ps).empty());

    ps.specifiedStereoQueryMatchesUnspecified = true;
    CHECK(SubstructMatch(*m1, *q, ps).size() == 1);
    CHECK(SubstructMatch(*m2, *q, ps).size() == 1);
    CHECK(SubstructMatch(*m3, *q, ps).empty());
  }
  SECTION("bond basics") {
    auto q = "F/C=C/Br"_smarts;
    REQUIRE(q);

    auto m1 = "F/C=C/Br"_smiles;
    REQUIRE(m1);
    auto m2 = "FC=CBr"_smiles;
    REQUIRE(m2);
    auto m3 = "F/C=C\\Br"_smiles;
    REQUIRE(m3);

    SubstructMatchParameters ps;
    ps.useChirality = true;
    CHECK(SubstructMatch(*m1, *q, ps).size() == 1);
    CHECK(SubstructMatch(*m2, *q, ps).empty());
    CHECK(SubstructMatch(*m3, *q, ps).empty());

    ps.specifiedStereoQueryMatchesUnspecified = true;
    CHECK(SubstructMatch(*m1, *q, ps).size() == 1);
    CHECK(SubstructMatch(*m2, *q, ps).size() == 1);
    CHECK(SubstructMatch(*m3, *q, ps).empty());
  }
}

TEST_CASE(
    "Github 8162: conjugated triple bonds match aromatic bonds with aromaticMatchesConjugated") {
  SECTION("as reported") {
    auto qry = R"CTAB(ACS Document 1996
  ChemDraw01092510212D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.357236 0.412500 0.000000 0
M  V30 2 C 1.071707 0.000000 0.000000 0
M  V30 3 C -0.357236 0.000000 0.000000 0
M  V30 4 N -1.071707 -0.412500 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 3 3 4
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(qry);
    auto mol = "CN1C=NC=C1"_smiles;
    REQUIRE(mol);
    SubstructMatchParameters ps;
    ps.aromaticMatchesConjugated = true;
    CHECK(SubstructMatch(*mol, *qry, ps).empty());
  }

  SECTION("details") {
    auto m_triple = "C#N"_smiles;
    REQUIRE(m_triple);
    m_triple->getBondWithIdx(0)->setIsConjugated(true);
    auto m_double = "C=N"_smiles;
    REQUIRE(m_double);
    m_double->getBondWithIdx(0)->setIsConjugated(true);
    auto m_single = "CN"_smiles;
    REQUIRE(m_single);
    m_single->getBondWithIdx(0)->setIsConjugated(true);
    auto m_aromatic = "CN"_smiles;
    REQUIRE(m_aromatic);
    m_aromatic->getBondWithIdx(0)->setBondType(Bond::AROMATIC);
    m_aromatic->getBondWithIdx(0)->setIsAromatic(true);
    m_aromatic->getBondWithIdx(0)->setIsConjugated(true);

    SubstructMatchParameters ps;
    ps.aromaticMatchesConjugated = true;
    CHECK(SubstructMatch(*m_triple, *m_aromatic, ps).empty());
    CHECK(SubstructMatch(*m_aromatic, *m_triple, ps).empty());

    CHECK(SubstructMatch(*m_double, *m_aromatic, ps).size() == 1);
    CHECK(SubstructMatch(*m_aromatic, *m_double, ps).size() == 1);
    CHECK(SubstructMatch(*m_single, *m_aromatic, ps).size() == 1);
    CHECK(SubstructMatch(*m_aromatic, *m_single, ps).size() == 1);
  }
}

TEST_CASE("Github #7295", "CIS/TRANS in aromatic ring") {
  SECTION("Smart CIS/TRANS bonds should match single or aromatic") {
    SubstructMatchParameters ps;

    auto query =
        "[O:1]=[c:2]1/[c:3](=[C:4]/[c:5]2:[c:10]:[c:9]:[c:8]:[c:7]:[c:6]:2):[s:11]:[c:12]2:[n:13]:1-[N:14]-[C:15]-[N:20]-[N:21]=2"_smarts;
    auto mol =
        "[O:1]=[c:2]1/[c:3](=[CH:4]/[c:5]2[cH:6][cH:7][cH:8][cH:9][cH:10]2)[s:11][c:12]2[n:13]1[NH:14][C:15]1([CH2:16][CH2:17][CH2:18][CH2:19]1)[NH:20][N:21]=2"_smiles;
    CHECK(SubstructMatch(*mol, *query, ps).size() == 1);
    ps.useChirality = true;
    CHECK(SubstructMatch(*mol, *query, ps).size() == 1);
    auto mol2 =
        "[O:1]=[c:2]1\\[c:3](=[CH:4]/[c:5]2[cH:6][cH:7][cH:8][cH:9][cH:10]2)[s:11][c:12]2[n:13]1[NH:14][C:15]1([CH2:16][CH2:17][CH2:18][CH2:19]1)[NH:20][N:21]=2"_smiles;
    CHECK(SubstructMatch(*mol2, *query, ps).size() == 0);
  }
}

TEST_CASE(
    "Github #8485: Allow single/double bonds to match aromatic in substructure search") {
  SECTION("basics") {
    auto m = "C1=CC=CC=C1"_smiles;
    REQUIRE(m);
    SubstructMatchParameters ps;
    ps.aromaticMatchesSingleOrDouble = true;

    {
      auto q = "[#6]:[#6]"_smarts;
      REQUIRE(q);
      CHECK(!SubstructMatch(*m, *q).empty());
      CHECK(!SubstructMatch(*m, *q, ps).empty());
    }
    {
      auto q = "C-C"_smiles;
      REQUIRE(q);
      CHECK(SubstructMatch(*m, *q).empty());
      CHECK(!SubstructMatch(*m, *q, ps).empty());
    }
    {
      auto q = "C=C"_smiles;
      REQUIRE(q);
      CHECK(SubstructMatch(*m, *q).empty());
      CHECK(!SubstructMatch(*m, *q, ps).empty());
    }
  }
  SECTION("as reported") {
    auto m = "C1=CC2=C(OC1=O)C(=CC=C2)C(=O)O"_smiles;
    REQUIRE(m);
    auto q = "CC(=O)OC1=CC=CC=C1C(=O)O"_smiles;
    REQUIRE(q);
    SubstructMatchParameters ps;
    ps.aromaticMatchesSingleOrDouble = true;
    CHECK(SubstructMatch(*m, *q).empty());
    CHECK(!SubstructMatch(*m, *q, ps).empty());
  }
  SECTION("symmetry") {
    auto m1 = "C1=CC=CC=C1"_smiles;
    REQUIRE(m1);
    auto m2 = "C1CCCCC1"_smiles;
    REQUIRE(m2);
    CHECK(SubstructMatch(*m1, *m2).empty());
    CHECK(SubstructMatch(*m2, *m1).empty());
    SubstructMatchParameters ps;
    ps.aromaticMatchesSingleOrDouble = true;
    CHECK(!SubstructMatch(*m1, *m2, ps).empty());
    CHECK(!SubstructMatch(*m2, *m1, ps).empty());
  }
}

TEST_CASE("extra atom and bond queries") {
  SECTION("basics") {
    auto m = "CCCC"_smiles;
    REQUIRE(m);
    m->getAtomWithIdx(1)->setFlags(0x3);
    m->getAtomWithIdx(2)->setFlags(0x5);
    m->getBondWithIdx(1)->setFlags(0x7);

    auto q = "CC"_smiles;
    REQUIRE(q);
    q->getAtomWithIdx(0)->setFlags(0x5);
    q->getAtomWithIdx(1)->setFlags(0x3);
    q->getBondWithIdx(0)->setFlags(0x7);

    SubstructMatchParameters ps;
    auto matches = SubstructMatch(*m, *q, ps);
    CHECK(matches.size() == 3);

    {
      SubstructMatchParameters ps;
      auto atomQuery = [](const Atom &queryAtom,
                          const Atom &targetAtom) -> bool {
        return queryAtom.getFlags() == targetAtom.getFlags();
      };
      ps.extraAtomCheck = atomQuery;
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 1);
      CHECK(matches[0][0].second == 2);
      CHECK(matches[0][1].second == 1);
    }
    {
      SubstructMatchParameters ps;
      auto bondQuery = [](const Bond &query, const Bond &target) -> bool {
        return query.getFlags() == target.getFlags();
      };
      ps.extraBondCheck = bondQuery;
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 1);
      CHECK(matches[0][0].second == 1);
      CHECK(matches[0][1].second == 2);
    }
  }
  SECTION("3D") {
    auto m = "CCCC |(0,0,0;1,0,0;2,0,0;3,0,0)|"_smiles;
    REQUIRE(m);
    auto q = "CC |(3,0,0;2,0,0)|"_smiles;
    REQUIRE(q);
    {
      SubstructMatchParameters ps;
      auto atomQuery = [](const Atom &queryAtom,
                          const Atom &targetAtom) -> bool {
        auto qconf = queryAtom.getOwningMol().getConformer();
        auto tconf = targetAtom.getOwningMol().getConformer();
        auto qpos = qconf.getAtomPos(queryAtom.getIdx());
        auto tpos = tconf.getAtomPos(targetAtom.getIdx());
        auto dist = (qpos - tpos).length();
        return dist < 0.1;
      };
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 3);

      ps.extraAtomCheck = atomQuery;
      matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 1);
      CHECK(matches[0][0].second == 3);
      CHECK(matches[0][1].second == 2);
    }
  }
  SECTION("AtomCoordsMatchFunctor") {
    auto m = "CCCC |(0,0,0;1,0,0;2,0,0;3,0,0)|"_smiles;
    REQUIRE(m);
    auto q = "CC |(3,0,0;2,0,0)|"_smiles;
    REQUIRE(q);
    SubstructMatchParameters ps;
    AtomCoordsMatchFunctor atomQuery;
    // I "<heart>" C++ syntax
    ps.extraAtomCheck =
        std::bind(&AtomCoordsMatchFunctor::operator(), &atomQuery,
                  std::placeholders::_1, std::placeholders::_2);
    {
      auto matches = SubstructMatch(*m, *q, ps);
      CHECK(matches.size() == 1);
      CHECK(matches[0][0].second == 3);
      CHECK(matches[0][1].second == 2);
    }
    {
      ROMol mcp(*m);
      mcp.clearConformers();
      auto matches = SubstructMatch(mcp, *q, ps);
      CHECK(matches.empty());
    }
    {
      ROMol qcp(*q);
      qcp.clearConformers();
      auto matches = SubstructMatch(*m, qcp, ps);
      CHECK(matches.empty());
    }
    {
      ROMol mcp(*m);
      mcp.clearConformers();
      ROMol qcp(*q);
      qcp.clearConformers();
      auto matches = SubstructMatch(mcp, qcp, ps);
      CHECK(matches.empty());
    }
    {
      // specifying conformer ID on the molecule
      ROMol mcp(*m);
      Conformer *conf = new Conformer(mcp.getConformer());
      mcp.getConformer().getAtomPos(3).z += 10;
      auto cid = mcp.addConformer(conf, true);
      auto matches = SubstructMatch(mcp, *q, ps);
      CHECK(matches.empty());

      SubstructMatchParameters ps2;
      AtomCoordsMatchFunctor atomQuery2(cid, -1);
      ps2.extraAtomCheck =
          std::bind(&AtomCoordsMatchFunctor::operator(), &atomQuery2,
                    std::placeholders::_1, std::placeholders::_2);
      matches = SubstructMatch(mcp, *q, ps2);
      CHECK(matches.size() == 1);
      CHECK(matches[0][0].second == 3);
      CHECK(matches[0][1].second == 2);
    }
    {
      // specifying conformer ID on the query
      ROMol qcp(*q);
      Conformer *conf = new Conformer(qcp.getConformer());
      qcp.getConformer().getAtomPos(0).z += 10;
      auto cid = qcp.addConformer(conf, true);
      auto matches = SubstructMatch(*m, qcp, ps);
      CHECK(matches.empty());

      SubstructMatchParameters ps2;
      AtomCoordsMatchFunctor atomQuery2(-1, cid);
      ps2.extraAtomCheck =
          std::bind(&AtomCoordsMatchFunctor::operator(), &atomQuery2,
                    std::placeholders::_1, std::placeholders::_2);
      matches = SubstructMatch(*m, qcp, ps2);
      CHECK(matches.size() == 1);
      CHECK(matches[0][0].second == 3);
      CHECK(matches[0][1].second == 2);
    }
  }
}