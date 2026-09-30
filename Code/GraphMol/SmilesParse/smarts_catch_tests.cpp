//
//  Copyright (C) 2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

TEST_CASE("Github #8424: direction on aromatic bonds in SMARTS") {
  SECTION("simplified") {
    auto m = "C/N=c1/[nH]cc(Br)nc1"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#6]/[#7]=[#6]1/[#7]:[#6]:[#6](-[#35]):[#7]:[#6]:1");
  }
  SECTION("as reported") {
    auto m =
        "CN1C(=O)CN(c2cc(C3CC3)cn3cc(CC(=O)N/N=c4\\cnc(Br)c[nH]4)nc23)C1=O"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    // should have slashes in both directions
    CHECK(smarts.find("/") != std::string::npos);
    CHECK(smarts.find("\\") != std::string::npos);
  }
}

TEST_CASE("repeated explicit H counts and charges") {
  SECTION("h counts") {
    std::vector<std::string> smartses = {
        "[N&H3&H0]",
        "[N&H0&H3]",
    };
    for (const auto &smarts : smartses) {
      INFO(smarts);
      auto m = v2::SmilesParse::MolFromSmarts(smarts);
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getNoImplicit());
      CHECK(m->getAtomWithIdx(0)->getNumExplicitHs() == 0);
    }
  }
  SECTION("charges") {
    std::vector<std::string> smartses = {"[N&+&+0]", "[N&+0&+2]"};
    for (const auto &smarts : smartses) {
      INFO(smarts);
      auto m = v2::SmilesParse::MolFromSmarts(smarts);
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(0)->getFormalCharge() == 0);
    }
  }
}

TEST_CASE("implicit Hs from SMILES should not make it into SMARTS") {
  SECTION("aromatic N") {
    auto m = "c1ccc[nH]1"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#6]1:[#6]:[#6]:[#6]:[#7]:1");
  }
  SECTION("chirality") {
    // as of this writing, we still keep Hs in the SMARTS for chiral centers.
    //   it's inconsistent to do so, but we have explicitly punted on fixing
    //   this for now.
    auto m = "C[C@H](N)F"_smiles;
    REQUIRE(m);
    auto smarts = MolToSmarts(*m);
    CHECK(smarts == "[#6]-[#6@H](-[#7])-[#9]");
  }
}

void checkMatches(const std::string &smarts, const std::string &smiles,
                  unsigned int nMatches, unsigned int lenFirst,
                  bool addHs = false) {
  // utility function that will find the matches between a smarts and smiles
  // if they match the expected values
  //  smarts : smarts string
  //  smiles : smiles string
  //  nMatches : expected number of matches
  //  lenFirst : length of the first match
  //
  // Return the list of all matches just in case want to do additional testing
  INFO(smarts + " " + smiles);
  auto matcher = v2::SmilesParse::MolFromSmarts(smarts);
  REQUIRE(matcher);
  // we will at the same time test the serialization:
  std::string pickle;
  MolPickler::pickleMol(*matcher, pickle);
  ROMol matcher2(pickle);

  auto mol = v2::SmilesParse::MolFromSmiles(smiles);
  REQUIRE(mol);
  if (addHs) {
    MolOps::addHs(*mol);
  }
  MolOps::findSSSR(*mol);

  MatchVectType mV;
  auto matches = SubstructMatch(*mol, *matcher, mV);
  CHECK(matches);
  CHECK(mV.size() == lenFirst);
  std::vector<MatchVectType> mVV;
  auto uniquify = true;
  auto matchCount = SubstructMatch(*mol, *matcher, mVV, uniquify);
  CHECK(matchCount == nMatches);
  CHECK(mVV[0].size() == lenFirst);

  matches = SubstructMatch(*mol, matcher2, mV);
  CHECK(matches);
  CHECK(mV.size() == lenFirst);
  matchCount = SubstructMatch(*mol, matcher2, mVV, true);
  CHECK(matchCount == nMatches);
  CHECK(mVV[0].size() == lenFirst);
}

TEST_CASE("k SMARTS extensions") {
  SECTION("parsing and writing") {
    auto q = "[k4]"_smarts;
    REQUIRE(q);
    auto smarts = MolToSmarts(*q);
    CHECK(smarts == "[k4]");
  }
  SECTION("matching") {
    auto m = "C1CC2N1CCCCCC2"_smiles;
    REQUIRE(m);
    auto k4 = "[k4]"_smarts;
    REQUIRE(k4);
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*m, *k4, matches));
    CHECK(matches.size() == 4);
    auto q4 = "[r4]"_smarts;
    REQUIRE(q4);
    CHECK(SubstructMatch(*m, *q4, matches));
    CHECK(matches.size() == 4);
    auto k8 = "[k8]"_smarts;
    REQUIRE(k8);
    CHECK(SubstructMatch(*m, *k8, matches));
    CHECK(matches.size() == 8);
    auto q8 = "[r8]"_smarts;
    REQUIRE(q8);
    CHECK(SubstructMatch(*m, *q8, matches));
    CHECK(matches.size() == 6);
  }
  SECTION("ranges") {
    std::string smiles = "C1CC2N1CCCCCC2C";
    std::vector<std::pair<std::string, size_t>> smartses = {
        {"[r{4-5}]", 4},  {"[k{4-5}]", 4},  {"[k{4-}]", 10}, {"[r{4-}]", 10},
        {"[k{-5}]", 4},   {"[k{8-}]", 8},   {"[r{8-}]", 6},  {"[r{4-8}]", 10},
        {"[!r{4-5}]", 7}, {"[!k{4-5}]", 7}, {"[!k{4-}]", 1}, {"[!r{4-}]", 1},
        {"[!k{-5}]", 7},  {"[!k{8-}]", 3},  {"[!r{8-}]", 5}, {"[!r{4-8}]", 1},
        {"[k]", 10},      {"[r]", 10},      {"[!k]", 1},     {"[!r]", 1},
    };
    for (const auto &[sma, val] : smartses) {
      checkMatches(sma, smiles, val, 1);
    }
  }
}

TEST_CASE("@{n} SMARTS bond ring count") {
  SECTION("parsing and writing") {
    auto q = "*@{0}*"_smarts;
    REQUIRE(q);
    CHECK(MolToSmarts(*q) == "*@{0}*");

    auto q2 = "*!@*"_smarts;
    REQUIRE(q2);
    CHECK(MolToSmarts(*q2) == "*!@*");

    auto q3 = "*@{2}*"_smarts;
    REQUIRE(q3);
    CHECK(MolToSmarts(*q3) == "*@{2}*");

    auto q4 = "*@{2-}*"_smarts;
    REQUIRE(q4);
    CHECK(MolToSmarts(*q4) == "*@{2-}*");

    auto q5 = "*@{-2}*"_smarts;
    REQUIRE(q5);
    CHECK(MolToSmarts(*q5) == "*@{-2}*");

    auto q6 = "*@{1-3}*"_smarts;
    REQUIRE(q6);
    CHECK(MolToSmarts(*q6) == "*@{1-3}*");
  }

  SECTION("matching") {
    auto naph = "c1ccc2ccccc2c1"_smiles;
    REQUIRE(naph);
    auto q = "*@{2}*"_smarts;
    REQUIRE(q);
    std::vector<MatchVectType> matches;
    CHECK(SubstructMatch(*naph, *q, matches));
    CHECK(matches.size() == 1);

    auto bip = "c1ccc2c(c1)c1ccccc12"_smiles;
    REQUIRE(bip);
    auto q2 = "*@{2}*"_smarts;
    REQUIRE(q2);
    matches.clear();
    CHECK(SubstructMatch(*bip, *q2, matches));
    CHECK(matches.size() == 2);

    auto cub = "C12C3C4C1C5C2C3C45"_smiles;
    REQUIRE(cub);
    auto q3 = "*@{2}*"_smarts;
    REQUIRE(q3);
    matches.clear();
    CHECK(SubstructMatch(*cub, *q3, matches));
    CHECK(matches.size() == 12);

    auto cyc = "C1CCCCC1"_smiles;
    REQUIRE(cyc);
    auto q4 = "*@{1}*"_smarts;
    REQUIRE(q4);
    matches.clear();
    CHECK(SubstructMatch(*cyc, *q4, matches));
    CHECK(matches.size() == 6);

    auto eth = "CC"_smiles;
    REQUIRE(eth);
    auto q5 = "*@{0}*"_smarts;
    REQUIRE(q5);
    matches.clear();
    CHECK(SubstructMatch(*eth, *q5, matches));
    CHECK(matches.size() == 1);
  }

  SECTION("ranges") {
    std::string smiles = "c1ccc2ccccc2c1";
    std::vector<std::pair<std::string, size_t>> smartses = {
        {"*@{1}*", 10},  {"*@{2}*", 1},   {"*@{2-}*", 1},  {"*@{1-2}*", 11},
        {"*@{-1}*", 10}, {"*@{-2}*", 11}, {"*!@{2}*", 10},
    };
    for (const auto &[sma, val] : smartses) {
      checkMatches(sma, smiles, val, 2);
    }
  }

  SECTION("composite bond range round-trip") {
    // @{1-2} combined with another bond primitive must not trigger an invalid
    // static_cast. The composite-bond path in _recurseBondSmarts passes the
    // raw QUERYBOND_QUERY * to getBondSmartsSimple, which does its own
    // dynamic_cast to BOND_RANGE_QUERY.
    auto q = "*@{1-2}&-*"_smarts;
    REQUIRE(q);
    CHECK(MolToSmarts(*q) == "*@{1-2}&-*");

    // Test less and greater forms in composite bonds too
    auto q2 = "*@{-2}&-*"_smarts;
    REQUIRE(q2);
    CHECK(MolToSmarts(*q2) == "*@{-2}&-*");

    auto q3 = "*@{3-}&-*"_smarts;
    REQUIRE(q3);
    CHECK(MolToSmarts(*q3) == "*@{3-}&-*");
  }

  SECTION("bond range query-query matching") {
    // QueryMatch does subset matching: pattern (q1) matches target (q2)
    // when the set of bonds matched by q1 is a subset of those matched by q2.
    SubstructMatchParameters params;
    params.useQueryQueryMatches = true;

    // Same range matches itself
    auto q1 = "*@{1-2}*"_smarts;
    REQUIRE(q1);
    auto q2 = "*@{1-2}*"_smarts;
    REQUIRE(q2);
    auto qb1 = static_cast<QueryBond *>(q1->getBondWithIdx(0));
    auto qb2 = static_cast<QueryBond *>(q2->getBondWithIdx(0));
    CHECK(qb1->QueryMatch(qb2));

    // Disjoint ranges do not match
    auto q3 = "*@{3-4}*"_smarts;
    REQUIRE(q3);
    auto qb3 = static_cast<QueryBond *>(q3->getBondWithIdx(0));
    CHECK(!qb1->QueryMatch(qb3));

    // Equality @{2} is subset of range @{1-3}
    auto qEq = "*@{2}*"_smarts;
    REQUIRE(qEq);
    auto qRange = "*@{1-3}*"_smarts;
    REQUIRE(qRange);
    auto qbEq = static_cast<QueryBond *>(qEq->getBondWithIdx(0));
    auto qbRange = static_cast<QueryBond *>(qRange->getBondWithIdx(0));
    CHECK(qbEq->QueryMatch(qbRange)); // pattern @{2} ⊆ target @{1-3} -> match
    CHECK(!qbRange->QueryMatch(qbEq)); // pattern @{1-3} ⊄ target @{2} -> no match

    // Range @{2-3} is subset of @{1-3}
    auto qSub = "*@{2-3}*"_smarts;
    REQUIRE(qSub);
    auto qbSub = static_cast<QueryBond *>(qSub->getBondWithIdx(0));
    CHECK(qbSub->QueryMatch(qbRange)); // pattern @{2-3} ⊆ target @{1-3} -> match
    CHECK(!qbRange->QueryMatch(qbSub)); // pattern @{1-3} ⊄ target @{2-3} -> no match

    // Equality @{4} is not subset of @{1-3}
    auto qEq4 = "*@{4}*"_smarts;
    REQUIRE(qEq4);
    auto qbEq4 = static_cast<QueryBond *>(qEq4->getBondWithIdx(0));
    CHECK(!qbEq4->QueryMatch(qbRange)); // pattern @{4} ⊄ target @{1-3}

    // Negation: @{2} vs !@{1-3} -> intervals overlap, negations differ -> no match
    auto qNegRange = "*!@{1-3}*"_smarts;
    REQUIRE(qNegRange);
    auto qbNegRange = static_cast<QueryBond *>(qNegRange->getBondWithIdx(0));
    CHECK(!qbEq->QueryMatch(qbNegRange)); // @{2} overlaps !@{1-3} -> no match

    // Negation: @{4} vs !@{1-3} -> intervals disjoint, negations differ -> match
    CHECK(qbEq4->QueryMatch(qbNegRange)); // @{4} disjoint from !@{1-3} -> match

    // Negation: !@{1-3} vs !@{1-2} -> both negated, target excluded [1,2] ⊆
    // pattern excluded [1,3] -> match
    auto qNegRange2 = "*!@{1-2}*"_smarts;
    REQUIRE(qNegRange2);
    auto qbNegRange2 = static_cast<QueryBond *>(qNegRange2->getBondWithIdx(0));
    CHECK(qbNegRange->QueryMatch(qbNegRange2)); // !@{1-3} ⊆ !@{1-2}

    // Negation: !@{1-2} vs !@{1-3} -> both negated, target excluded [1,3] ⊄
    // pattern excluded [1,2] -> no match
    CHECK(!qbNegRange2->QueryMatch(qbNegRange)); // !@{1-2} ⊄ !@{1-3}

    // --- Open-ended ranges ---
    // @{-N} is LessEqualQuery(N) with interval [INT_MIN, N]
    // @{N-} is GreaterEqualQuery(N) with interval [N, INT_MAX]

    // Case 1: both positive, LessEqualQuery (@{-N})
    auto qLE3 = "*@{-3}*"_smarts; // LessEqualQuery(3) → [MIN, 3]
    REQUIRE(qLE3);
    auto qLE5 = "*@{-5}*"_smarts; // LessEqualQuery(5) → [MIN, 5]
    REQUIRE(qLE5);
    auto qbLE3 = static_cast<QueryBond *>(qLE3->getBondWithIdx(0));
    auto qbLE5 = static_cast<QueryBond *>(qLE5->getBondWithIdx(0));
    auto qEq3 = "*@{3}*"_smarts;
    REQUIRE(qEq3);
    auto qbEq3 = static_cast<QueryBond *>(qEq3->getBondWithIdx(0));
    CHECK(qbEq3->QueryMatch(qbLE3)); // @{3} ⊆ @{-3} [3,3] ⊆ [MIN,3]
    CHECK(qbLE3->QueryMatch(qbLE5)); // @{-3} ⊆ @{-5} [MIN,3] ⊆ [MIN,5]
    CHECK(!qbLE5->QueryMatch(qbLE3)); // @{-5} ⊄ @{-3} [MIN,5] ⊄ [MIN,3]

    // Case 1: both positive, GreaterEqualQuery (@{N-})
    auto qGE2 = "*@{2-}*"_smarts; // GreaterEqualQuery(2) → [2, MAX]
    REQUIRE(qGE2);
    auto qGE3 = "*@{3-}*"_smarts; // GreaterEqualQuery(3) → [3, MAX]
    REQUIRE(qGE3);
    auto qbGE2 = static_cast<QueryBond *>(qGE2->getBondWithIdx(0));
    auto qbGE3 = static_cast<QueryBond *>(qGE3->getBondWithIdx(0));
    CHECK(qbGE3->QueryMatch(qbGE2)); // @{3-} ⊆ @{2-} [3,MAX] ⊆ [2,MAX]
    CHECK(!qbGE2->QueryMatch(qbGE3)); // @{2-} ⊄ @{3-} [2,MAX] ⊄ [3,MAX]

    // Case 1: cross-type LessEqual vs GreaterEqual
    CHECK(!qbLE3->QueryMatch(qbGE2)); // @{-3} ⊄ @{2-} [MIN,3] ⊄ [2,MAX]
    CHECK(!qbGE2->QueryMatch(qbLE3)); // @{2-} ⊄ @{-3} [2,MAX] ⊄ [MIN,3]

    // Case 2: both negated, open-ended
    auto qNegLE3 = "*!@{-3}*"_smarts; // !LessEqualQuery(3), excludes [MIN,3]
    REQUIRE(qNegLE3);
    auto qNegLE5 = "*!@{-5}*"_smarts; // !LessEqualQuery(5), excludes [MIN,5]
    REQUIRE(qNegLE5);
    auto qbNegLE3 = static_cast<QueryBond *>(qNegLE3->getBondWithIdx(0));
    auto qbNegLE5 = static_cast<QueryBond *>(qNegLE5->getBondWithIdx(0));
    CHECK(!qbNegLE3->QueryMatch(qbNegLE5)); // !@{-3} ⊄ !@{-5} target excl [MIN,5] ⊄ pattern excl [MIN,3]
    CHECK(qbNegLE5->QueryMatch(qbNegLE3)); // !@{-5} ⊆ !@{-3} target excl [MIN,3] ⊆ pattern excl [MIN,5]

    // Case 3: pattern negated, target positive (complement must fit)
    CHECK(qbNegLE3->QueryMatch(qbGE2)); // !@{-3} ⊆ @{2-} complement [4,MAX] ⊆ [2,MAX]
    CHECK(!qbNegLE3->QueryMatch(qbLE5)); // !@{-3} ⊄ @{-5} complement [4,MAX] ⊄ [MIN,5]
    auto qNegGE2 = "*!@{2-}*"_smarts; // !GreaterEqualQuery(2), excludes [2,MAX]
    REQUIRE(qNegGE2);
    auto qbNegGE2 = static_cast<QueryBond *>(qNegGE2->getBondWithIdx(0));
    CHECK(qbNegGE2->QueryMatch(qbLE3)); // !@{2-} ⊆ @{-3} complement [MIN,1] ⊆ [MIN,3]
    auto qLE4 = "*@{-4}*"_smarts; // LessEqualQuery(4) → [MIN, 4]
    REQUIRE(qLE4);
    auto qbLE4 = static_cast<QueryBond *>(qLE4->getBondWithIdx(0));
    CHECK(qbNegGE2->QueryMatch(qbLE4)); // !@{2-} ⊆ @{-4} complement [MIN,1] ⊆ [MIN,4]

    // Case 4: pattern positive, target negated (disjoint)
    CHECK(!qbGE2->QueryMatch(qbNegLE5)); // @{2-} NOT disjoint from !@{-5} [2,MAX] overlaps [MIN,5]
    CHECK(!qbLE3->QueryMatch(qbNegGE2)); // @{-3} NOT disjoint from !@{2-} [MIN,3] overlaps [2,MAX]
    auto qNegGE1 = "*!@{1-}*"_smarts; // !GreaterEqualQuery(1), excludes [1,MAX]
    REQUIRE(qNegGE1);
    auto qbNegGE1 = static_cast<QueryBond *>(qNegGE1->getBondWithIdx(0));
    CHECK(!qbGE2->QueryMatch(qbNegGE1)); // @{2-} NOT disjoint from !@{1-} [2,MAX] overlaps [1,MAX]
    CHECK(!qbLE3->QueryMatch(qbNegLE5)); // @{-3} NOT disjoint from !@{-5} [MIN,3] overlaps [MIN,5]

    // Substruct match with query-query matching
    {
      auto matches = SubstructMatch(*qRange, *qEq, params);
      CHECK(matches.size() == 1); // pattern @{2} ⊆ target @{1-3}
      matches = SubstructMatch(*qEq, *qRange, params);
      CHECK(matches.empty()); // pattern @{1-3} ⊄ target @{2}
    }
  }

  SECTION("atom range query-query matching") {
    SubstructMatchParameters params;
    params.useQueryQueryMatches = true;

    // Same range matches itself
    auto q1 = "[R{1-2}]"_smarts;
    REQUIRE(q1);
    auto q2 = "[R{1-2}]"_smarts;
    REQUIRE(q2);
    auto qa1 = static_cast<QueryAtom *>(q1->getAtomWithIdx(0));
    auto qa2 = static_cast<QueryAtom *>(q2->getAtomWithIdx(0));
    CHECK(qa1->QueryMatch(qa2));

    // Disjoint ranges do not match
    auto q3 = "[R{3-4}]"_smarts;
    REQUIRE(q3);
    auto qa3 = static_cast<QueryAtom *>(q3->getAtomWithIdx(0));
    CHECK(!qa1->QueryMatch(qa3));

    // Equality R2 is subset of range R{1-3}
    auto qEq = "[R2]"_smarts;
    REQUIRE(qEq);
    auto qRange = "[R{1-3}]"_smarts;
    REQUIRE(qRange);
    auto qaEq = static_cast<QueryAtom *>(qEq->getAtomWithIdx(0));
    auto qaRange = static_cast<QueryAtom *>(qRange->getAtomWithIdx(0));
    CHECK(qaEq->QueryMatch(qaRange)); // pattern R2 ⊆ target R{1-3}
    CHECK(!qaRange->QueryMatch(qaEq)); // pattern R{1-3} ⊄ target R2

    // Ring size ranges
    auto q8 = "[r{5-6}]"_smarts;
    REQUIRE(q8);
    auto q9 = "[r{5-6}]"_smarts;
    REQUIRE(q9);
    auto qa8 = static_cast<QueryAtom *>(q8->getAtomWithIdx(0));
    auto qa9 = static_cast<QueryAtom *>(q9->getAtomWithIdx(0));
    CHECK(qa8->QueryMatch(qa9));

    auto q10 = "[r{4-5}]"_smarts;
    REQUIRE(q10);
    auto qa10 = static_cast<QueryAtom *>(q10->getAtomWithIdx(0));
    CHECK(!qa8->QueryMatch(qa10)); // r{5-6} ⊄ r{4-5}

    // Negation: R2 vs !R{1-3} -> intervals overlap, negations differ -> no match
    auto qNegRange = "[!R{1-3}]"_smarts;
    REQUIRE(qNegRange);
    auto qaNegRange = static_cast<QueryAtom *>(qNegRange->getAtomWithIdx(0));
    CHECK(!qaEq->QueryMatch(qaNegRange)); // R2 overlaps !R{1-3} -> no match

    // Negation: R4 vs !R{1-3} -> intervals disjoint, negations differ -> match
    auto qEq4 = "[R4]"_smarts;
    REQUIRE(qEq4);
    auto qaEq4 = static_cast<QueryAtom *>(qEq4->getAtomWithIdx(0));
    CHECK(qaEq4->QueryMatch(qaNegRange)); // R4 disjoint from !R{1-3} -> match

    // --- Open-ended ranges ---
    // R{-N} is LessEqualQuery(N) with interval [INT_MIN, N]
    // R{N-} is GreaterEqualQuery(N) with interval [N, INT_MAX]

    // Case 1: both positive, LessEqualQuery (R{-N})
    auto qLE3 = "[R{-3}]"_smarts; // LessEqualQuery(3) → [MIN, 3]
    REQUIRE(qLE3);
    auto qLE5 = "[R{-5}]"_smarts; // LessEqualQuery(5) → [MIN, 5]
    REQUIRE(qLE5);
    auto qaLE3 = static_cast<QueryAtom *>(qLE3->getAtomWithIdx(0));
    auto qaLE5 = static_cast<QueryAtom *>(qLE5->getAtomWithIdx(0));
    CHECK(qaLE3->QueryMatch(qaLE5)); // R{-3} ⊆ R{-5} [MIN,3] ⊆ [MIN,5]
    CHECK(!qaLE5->QueryMatch(qaLE3)); // R{-5} ⊄ R{-3} [MIN,5] ⊄ [MIN,3]

    // Case 1: both positive, GreaterEqualQuery (R{N-})
    auto qGE2 = "[R{2-}]"_smarts; // GreaterEqualQuery(2) → [2, MAX]
    REQUIRE(qGE2);
    auto qGE3 = "[R{3-}]"_smarts; // GreaterEqualQuery(3) → [3, MAX]
    REQUIRE(qGE3);
    auto qaGE2 = static_cast<QueryAtom *>(qGE2->getAtomWithIdx(0));
    auto qaGE3 = static_cast<QueryAtom *>(qGE3->getAtomWithIdx(0));
    CHECK(qaGE3->QueryMatch(qaGE2)); // R{3-} ⊆ R{2-} [3,MAX] ⊆ [2,MAX]
    CHECK(!qaGE2->QueryMatch(qaGE3)); // R{2-} ⊄ R{3-} [2,MAX] ⊄ [3,MAX]

    // Case 1: cross-type LessEqual vs GreaterEqual
    CHECK(!qaLE3->QueryMatch(qaGE2)); // R{-3} ⊄ R{2-} [MIN,3] ⊄ [2,MAX]
    CHECK(!qaGE2->QueryMatch(qaLE3)); // R{2-} ⊄ R{-3} [2,MAX] ⊄ [MIN,3]

    // Case 2: both negated
    auto qNegR13 = "[!R{1-3}]"_smarts;
    REQUIRE(qNegR13);
    auto qNegR12 = "[!R{1-2}]"_smarts;
    REQUIRE(qNegR12);
    auto qaNegR13 = static_cast<QueryAtom *>(qNegR13->getAtomWithIdx(0));
    auto qaNegR12 = static_cast<QueryAtom *>(qNegR12->getAtomWithIdx(0));
    CHECK(qaNegR13->QueryMatch(qaNegR12)); // !R{1-3} ⊆ !R{1-2}
    CHECK(!qaNegR12->QueryMatch(qaNegR13)); // !R{1-2} ⊄ !R{1-3}

    // Case 2: both negated, open-ended
    auto qNegLE3 = "[!R{-3}]"_smarts; // !LessEqualQuery(3), excludes [MIN,3]
    REQUIRE(qNegLE3);
    auto qNegLE5 = "[!R{-5}]"_smarts; // !LessEqualQuery(5), excludes [MIN,5]
    REQUIRE(qNegLE5);
    auto qaNegLE3 = static_cast<QueryAtom *>(qNegLE3->getAtomWithIdx(0));
    auto qaNegLE5 = static_cast<QueryAtom *>(qNegLE5->getAtomWithIdx(0));
    CHECK(!qaNegLE3->QueryMatch(qaNegLE5)); // !R{-3} ⊄ !R{-5} target excl [MIN,5] ⊄ pattern excl [MIN,3]
    CHECK(qaNegLE5->QueryMatch(qaNegLE3)); // !R{-5} ⊆ !R{-3} target excl [MIN,3] ⊆ pattern excl [MIN,5]

    // Case 3: pattern negated, target positive (complement must fit)
    auto qNegR1 = "[!R1]"_smarts;
    REQUIRE(qNegR1);
    auto qaNegR1 = static_cast<QueryAtom *>(qNegR1->getAtomWithIdx(0));
    CHECK(!qaNegR1->QueryMatch(qa3)); // !R1 ⊄ R{3-4} two pieces can't fit
    CHECK(qaNegLE3->QueryMatch(qaGE2)); // !R{-3} ⊆ R{2-} complement [4,MAX] ⊆ [2,MAX]
    CHECK(!qaNegLE3->QueryMatch(qaLE5)); // !R{-3} ⊄ R{-5} complement [4,MAX] ⊄ [MIN,5]
    auto qNegGE2 = "[!R{2-}]"_smarts; // !GreaterEqualQuery(2), excludes [2,MAX]
    REQUIRE(qNegGE2);
    auto qaNegGE2 = static_cast<QueryAtom *>(qNegGE2->getAtomWithIdx(0));
    CHECK(qaNegGE2->QueryMatch(qaLE3)); // !R{2-} ⊆ R{-3} complement [MIN,1] ⊆ [MIN,3]
    auto qLE4 = "[R{-4}]"_smarts; // LessEqualQuery(4) → [MIN, 4]
    REQUIRE(qLE4);
    auto qaLE4 = static_cast<QueryAtom *>(qLE4->getAtomWithIdx(0));
    CHECK(qaNegGE2->QueryMatch(qaLE4)); // !R{2-} ⊆ R{-4} complement [MIN,1] ⊆ [MIN,4]

    // Case 4: pattern positive, target negated (disjoint)
    CHECK(qa3->QueryMatch(qaNegR1)); // R{3-4} disjoint from !R1
    auto qR1 = "[R1]"_smarts;
    REQUIRE(qR1);
    auto qaR1 = static_cast<QueryAtom *>(qR1->getAtomWithIdx(0));
    CHECK(!qaR1->QueryMatch(qaNegR1)); // R1 NOT disjoint from !R1
    CHECK(!qaR1->QueryMatch(qaNegLE3)); // R1 NOT disjoint from !R{-3} [1,1] overlaps [MIN,3]
    CHECK(!qaLE3->QueryMatch(qaNegGE2)); // R{-3} NOT disjoint from !R{2-} [MIN,3] overlaps [2,MAX]
    CHECK(!qaGE2->QueryMatch(qaNegLE3)); // R{2-} NOT disjoint from !R{-3} [2,MAX] overlaps [MIN,3]
    auto qNegLE2 = "[!R{-2}]"_smarts; // !LessEqualQuery(2), excludes [MIN,2]
    REQUIRE(qNegLE2);
    auto qaNegLE2 = static_cast<QueryAtom *>(qNegLE2->getAtomWithIdx(0));
    CHECK(!qaGE2->QueryMatch(qaNegLE2)); // R{2-} NOT disjoint from !R{-2} touch at 2
    CHECK(!qaLE3->QueryMatch(qaNegLE5)); // R{-3} NOT disjoint from !R{-5} [MIN,3] overlaps [MIN,5]

    // Substruct match with query-query matching
    {
      auto matches = SubstructMatch(*qRange, *qEq, params);
      CHECK(matches.size() == 1);
      matches = SubstructMatch(*qEq, *qRange, params);
      CHECK(matches.empty());
    }
  }
}
