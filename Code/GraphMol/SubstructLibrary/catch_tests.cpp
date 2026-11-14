//
//  Copyright (C) 2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/PatternFactory.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>

using namespace RDKit;

TEST_CASE("querying with a molbundle") {
  std::vector<std::string> libSmiles = {"CCCC", "CCOC", "CCNC"};
  SubstructLibrary ssslib;
  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SubstructMatchParameters params;
  SECTION("basics") {
    std::vector<std::string> qSmiles = {"CCC", "COC", "CNC"};
    MolBundle bundle;
    for (const auto &smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, params, numThreads) == 3);
      CHECK(ssslib.hasMatch(bundle, params, numThreads));
      auto libMatches = ssslib.getMatches(bundle, params, numThreads);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
  }
  SECTION("some bundle members don't match") {
    std::vector<std::string> qSmiles = {"CCC", "CSC", "CNC"};
    MolBundle bundle;
    for (const auto &smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, params, numThreads) == 2);
      CHECK(ssslib.hasMatch(bundle, params, numThreads));
      auto libMatches = ssslib.getMatches(bundle, params, numThreads);
      CHECK(libMatches.size() == 2);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 2});
    }
  }
  SECTION("no dupes please") {
    std::vector<std::string> qSmiles = {"CCC", "CNC", "CC"};
    MolBundle bundle;
    for (const auto &smi : qSmiles) {
      boost::shared_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bundle.addMol(mol);
    }
    for (auto numThreads : std::vector<int>{1, -1}) {
      CHECK(ssslib.countMatches(bundle, params, numThreads) == 3);
      CHECK(ssslib.hasMatch(bundle, params, numThreads));
      auto libMatches = ssslib.getMatches(bundle, params, numThreads);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
  }
}

TEST_CASE("using modified query parameters") {
  std::vector<std::string> libSmiles = {"C[C@H](F)Cl", "C[C@@H](F)Cl",
                                        "C[CH](F)Cl"};
  SubstructLibrary ssslib;
  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("basics") {
    auto qm = "C[C@@H](F)Cl"_smiles;
    {  // by default stereo is not used
      SubstructMatchParameters params;
      CHECK(ssslib.countMatches(*qm, params) == 3);
      CHECK(ssslib.hasMatch(*qm, params));
      auto libMatches = ssslib.getMatches(*qm, params);
      CHECK(libMatches.size() == 3);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    }
    {  // use stereo
      SubstructMatchParameters params;
      params.useChirality = true;
      CHECK(ssslib.countMatches(*qm, params) == 1);
      CHECK(ssslib.hasMatch(*qm, params));
      auto libMatches = ssslib.getMatches(*qm, params);
      CHECK(libMatches.size() == 1);
      std::sort(libMatches.begin(), libMatches.end());
      CHECK(libMatches == std::vector<unsigned int>{1});
    }
  }
}

TEST_CASE("searchOrder") {
  std::vector<std::string> libSmiles = {"CCCOC", "CCCCOCC", "CCOC", "COC",
                                        "CCCCCOC"};
  boost::shared_ptr<MolHolder> mholder(new MolHolder());
  boost::shared_ptr<PatternHolder> fpholder(new PatternHolder());

  SubstructLibrary ssslib(mholder, fpholder);

  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("basics") {
    auto qm = "COC"_smiles;
    for (auto numThreads : std::vector<int>{1, -1}) {
      {  // default search order
        ssslib.setSearchOrder(std::vector<unsigned int>());
        unsigned int maxResults = 2;
        SubstructMatchParameters params;
        CHECK(ssslib.countMatches(*qm, params, numThreads) == 5);
        CHECK(ssslib.hasMatch(*qm, params));
        auto libMatches =
            ssslib.getMatches(*qm, params, numThreads, maxResults);
        CHECK(libMatches.size() == 2);
        CHECK(libMatches == std::vector<unsigned int>{0, 1});
      }
      {  // specified search order
        unsigned int maxResults = 2;
        std::vector<unsigned int> searchOrder = {3, 2, 0, 1, 4};
        ssslib.setSearchOrder(searchOrder);
        SubstructMatchParameters params;
        CHECK(ssslib.countMatches(*qm, params, numThreads) == 5);
        CHECK(ssslib.hasMatch(*qm, params));
        auto libMatches =
            ssslib.getMatches(*qm, params, numThreads, maxResults);
        CHECK(libMatches.size() == 2);
        CHECK(libMatches == std::vector<unsigned int>{3, 2});
      }
      {  // search order not the full length
        unsigned int maxResults = 2;
        std::vector<unsigned int> searchOrder = {3, 2, 0};
        ssslib.setSearchOrder(searchOrder);
        SubstructMatchParameters params;
        CHECK(ssslib.countMatches(*qm, params, numThreads) == 3);
        CHECK(ssslib.hasMatch(*qm, params));
        auto libMatches =
            ssslib.getMatches(*qm, params, numThreads, maxResults);
        CHECK(libMatches.size() == 2);
        CHECK(libMatches == std::vector<unsigned int>{3, 2});
      }
    }
  }

  SECTION("overly long searchOrder") {
    auto qm = "COC"_smiles;
    unsigned int maxResults = 2;
    int numThreads = -1;
    std::vector<unsigned int> searchOrder = {3, 2, 0, 1, 4, 2, 1, 4, 0};
    ssslib.setSearchOrder(searchOrder);
    SubstructMatchParameters params;
    CHECK(ssslib.countMatches(*qm, params, numThreads) == 5);
    CHECK(ssslib.hasMatch(*qm, params));
    auto libMatches = ssslib.getMatches(*qm, params, numThreads, maxResults);
    CHECK(libMatches.size() == 2);
    CHECK(libMatches == std::vector<unsigned int>{3, 2});
  }
  SECTION("bogus indices in searchOrder") {
    auto qm = "COC"_smiles;
    std::vector<unsigned int> searchOrder = {3, 12, 0};
    CHECK_THROWS_AS(ssslib.setSearchOrder(searchOrder), IndexErrorException);
  }
#ifdef RDK_USE_BOOST_SERIALIZATION
  SECTION("serialization") {
    std::vector<unsigned int> searchOrder = {3, 0, 1};
    ssslib.setSearchOrder(searchOrder);
    std::string pickle = ssslib.Serialize();
    SubstructLibrary serialized;
    serialized.initFromString(pickle);
    CHECK(serialized.size() == ssslib.size());
    CHECK(serialized.getSearchOrder() == ssslib.getSearchOrder());
  }
#endif
}

void setSearchSmallestFirst(SubstructLibrary &ssslib) {
  const auto holder = ssslib.getMolHolder();
  std::vector<unsigned int> searchOrder(holder->size());
  std::iota(searchOrder.begin(), searchOrder.end(), 0);
  std::stable_sort(searchOrder.begin(), searchOrder.end(),
                   [holder](unsigned int i1, unsigned int i2) {
                     return holder->getMol(i1)->getNumAtoms() <
                            holder->getMol(i2)->getNumAtoms();
                   });
  ssslib.setSearchOrder(searchOrder);
}

TEST_CASE("searchOrderFunctionDemo") {
  std::vector<std::string> libSmiles = {"CCCOC", "CCCCOCC", "CCOC", "COC",
                                        "CCCCCOC"};
  boost::shared_ptr<MolHolder> mholder(new MolHolder());
  boost::shared_ptr<PatternHolder> fpholder(new PatternHolder());

  SubstructLibrary ssslib(mholder, fpholder);

  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("basics") {
    auto qm = "COC"_smiles;
    setSearchSmallestFirst(ssslib);
    auto libMatches = ssslib.getMatches(*qm);
    CHECK(libMatches.size() == 5);
    CHECK(libMatches == std::vector<unsigned int>{3, 2, 0, 1, 4});
  }
}

#ifdef RDK_USE_BOOST_SERIALIZATION
TEST_CASE("ExtendedQueryMol") {
  std::vector<std::string> libSmiles = {"COCC=O", "COOCC=O", "COOOCC=O",
                                        "COOOOCC=O"};
  boost::shared_ptr<MolHolder> mholder(new MolHolder());
  boost::shared_ptr<PatternHolder> fpholder(new TautomerPatternHolder());

  SubstructLibrary ssslib(mholder, fpholder);

  for (const auto &smi : libSmiles) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ssslib.addMol(*mol);
  }
  SECTION("mol") {
    auto mol = "COCC"_smiles;
    auto xqm = GeneralizedSubstruct::createExtendedQueryMol(*mol);

    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.size() == 1);
    CHECK(libMatches == std::vector<unsigned int>{0});
    CHECK(ssslib.countMatches(xqm) == 1);
    CHECK(ssslib.hasMatch(xqm));
  }
  SECTION("tautomer query") {
    auto mol = "COCC"_smiles;
    auto xqm = GeneralizedSubstruct::ExtendedQueryMol(
        std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol)));
    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.size() == 1);
    CHECK(libMatches == std::vector<unsigned int>{0});
    CHECK(ssslib.countMatches(xqm) == 1);
    CHECK(ssslib.hasMatch(xqm));
  }
  SECTION("mol bundle") {
    auto mol = "COCC |LN:1:1.3|"_smiles;
    auto xqm = GeneralizedSubstruct::createExtendedQueryMol(*mol);

    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.size() == 3);
    CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    CHECK(ssslib.countMatches(xqm) == 3);
    CHECK(ssslib.hasMatch(xqm));
  }
  SECTION("tautomer bundle") {
    auto mol = "COC=CO |LN:1:1.3|"_smiles;
    auto xqm = GeneralizedSubstruct::createExtendedQueryMol(*mol);

    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.size() == 3);
    CHECK(libMatches == std::vector<unsigned int>{0, 1, 2});
    CHECK(ssslib.countMatches(xqm) == 3);
    CHECK(ssslib.hasMatch(xqm));
  }
  SECTION("mol no match") {
    auto mol = "CNCC"_smiles;
    auto xqm = GeneralizedSubstruct::createExtendedQueryMol(*mol);

    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.empty());
    CHECK(ssslib.countMatches(xqm) == 0);
    CHECK(!ssslib.hasMatch(xqm));
  }
  SECTION("tautomer query no match") {
    auto mol = "CNCC"_smiles;
    auto xqm = GeneralizedSubstruct::ExtendedQueryMol(
        std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*mol)));
    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.empty());
    CHECK(ssslib.countMatches(xqm) == 0);
    CHECK(!ssslib.hasMatch(xqm));
  }
  SECTION("mol bundle no match") {
    auto mol = "CNCC |LN:1:1.3|"_smiles;
    auto xqm = GeneralizedSubstruct::createExtendedQueryMol(*mol);

    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.empty());
    CHECK(ssslib.countMatches(xqm) == 0);
    CHECK(!ssslib.hasMatch(xqm));
  }
  SECTION("tautomer bundle no match") {
    auto mol = "CNC=CO |LN:1:1.3|"_smiles;
    auto xqm = GeneralizedSubstruct::createExtendedQueryMol(*mol);

    auto libMatches = ssslib.getMatches(xqm);
    CHECK(libMatches.empty());
    CHECK(ssslib.countMatches(xqm) == 0);
    CHECK(!ssslib.hasMatch(xqm));
  }
}
#endif
TEST_CASE("github #W9835: ring query inside a recursive SMARTS and TrustedSmilesMolHolder") {
  // query_needs_rings() decides whether the holder's molecules need
  // symmetrizeSSSR() before matching.
  //  Previously this only supported RecursiveDescent not AtomAnd/AtomOr/AtomNot
  const std::vector<std::string> libSmiles = {"CCc1ccccc1", "CCCC",
                                              "O=Cc1ccccc1C(=O)O"};
  const std::vector<std::pair<std::string, std::string>> queries = {
      {"top-level ring query", "[c;r6]"},
      {"recursive query is the root", "[$([c;r6])]"},
      {"ring query under AtomAnd", "[c;$([c;r6])]"},
      {"ring query under AtomOr", "[c,$([c;r6])]"},
      {"ring query under AtomNot", "[!$([c;r6])]"},
      {"ring query two recursions deep", "[c;$([c;$([c;r6])])]"},
      {"recursive query, no ring query", "[C;$(CC)]"},
      {"ring count query under AtomAnd", "[C;$([C;R0])]"},
  };

  for (const auto &[label, smarts] : queries) {
    INFO(label << ": " << smarts);
    std::unique_ptr<ROMol> query(SmartsToMol(smarts));
    REQUIRE(query);

    // ground truth from fully sanitized molecules
    std::vector<unsigned int> expected;
    for (unsigned int i = 0; i < libSmiles.size(); ++i) {
      std::unique_ptr<ROMol> mol(SmilesToMol(libSmiles[i]));
      REQUIRE(mol);
      if (SubstructMatch(*mol, *query).size()) {
        expected.push_back(i);
      }
    }

    auto molHolder = boost::make_shared<CachedTrustedSmilesMolHolder>();
    auto patternHolder = boost::make_shared<PatternHolder>();
    SubstructLibrary lib(molHolder, patternHolder);
    for (const auto &smi : libSmiles) {
      std::unique_ptr<RWMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      lib.addMol(*mol);
    }
    // threw before the fix rather than returning a wrong answer
    CHECK(lib.getMatches(*query) == expected);
  }
}

TEST_CASE("CachedTrustedSmilesMolHolder should assign radicals") {
  // The holder parses with sanitize=false, and that "light" sanitization does
  // not assign radical electrons, so a radical query matched nothing until
  // getMol() called assignRadicals().  There is no SMARTS primitive for
  // radicals, so the query is built directly.
  //
  // The store is plain SMILES, so a radical count is INFERRED from valence
  // rather than read back.  Every species below therefore has the
  // valence-implied count; the last SECTION pins what that excludes.
  auto radicalQuery = [](int nRadicals) {
    auto *qa = new QueryAtom();
    qa->setQuery(makeAtomNumRadicalElectronsQuery(nRadicals));
    auto *q = new RWMol();
    q->addAtom(qa, true, true);
    return std::unique_ptr<ROMol>(static_cast<ROMol *>(q));
  };
  auto totalRadicals = [](const ROMol &m) {
    unsigned int n = 0;
    for (const auto atom : m.atoms()) {
      n += atom->getNumRadicalElectrons();
    }
    return n;
  };

  const std::vector<std::string> libSmiles = {
      "[CH3]",  // 1 radical electron
      "CC",     // none
      "[NH2]",  // 1
      "[O]",    // 2
  };

  SECTION("radical queries match through the holder") {
    for (int nRadicals : {0, 1, 2}) {
      INFO("querying for " << nRadicals << " radical electron(s)");
      auto query = radicalQuery(nRadicals);

      // ground truth from fully sanitized molecules
      std::vector<unsigned int> expected;
      for (unsigned int i = 0; i < libSmiles.size(); ++i) {
        std::unique_ptr<ROMol> mol(SmilesToMol(libSmiles[i]));
        REQUIRE(mol);
        if (SubstructMatch(*mol, *query).size()) {
          expected.push_back(i);
        }
      }

      auto molHolder = boost::make_shared<CachedSmilesMolHolder>();
      SubstructLibrary lib(molHolder);
      for (const auto &smi : libSmiles) {
        molHolder->addSmiles(smi);
      }
      CHECK(lib.getMatches(*query) == expected);
    }
  }

  SECTION("per-atom radical counts survive the holder") {
    //  addSmiles() stores the string verbatim, so atom order is preserved and
    //  a per-atom comparison is meaningful.  addMol() would canonicalise.
    auto molHolder = boost::make_shared<CachedSmilesMolHolder>();
    for (const auto &smi : libSmiles) {
      molHolder->addSmiles(smi);
    }
    for (unsigned int i = 0; i < libSmiles.size(); ++i) {
      INFO(libSmiles[i]);
      std::unique_ptr<ROMol> ref(SmilesToMol(libSmiles[i]));
      REQUIRE(ref);
      auto held = molHolder->getMol(i);
      REQUIRE(held);
      REQUIRE(held->getNumAtoms() == ref->getNumAtoms());
      for (unsigned int a = 0; a < ref->getNumAtoms(); ++a) {
        CHECK(held->getAtomWithIdx(a)->getNumRadicalElectrons() ==
              ref->getAtomWithIdx(a)->getNumRadicalElectrons());
      }
    }
  }

  SECTION("addMol keeps the radical count") {
    //  addMol() stores MolToSmiles(), which reorders: the COUNT has to survive
    //  even where the index does not.
    auto molHolder = boost::make_shared<CachedSmilesMolHolder>();
    SubstructLibrary lib(molHolder);
    for (const auto &smi : libSmiles) {
      std::unique_ptr<RWMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      lib.addMol(*mol);
    }
    for (unsigned int i = 0; i < libSmiles.size(); ++i) {
      INFO(libSmiles[i]);
      std::unique_ptr<ROMol> ref(SmilesToMol(libSmiles[i]));
      REQUIRE(ref);
      auto held = molHolder->getMol(i);
      REQUIRE(held);
      CHECK(totalRadicals(*held) == totalRadicals(*ref));
    }
  }

  SECTION("radicals are inferred, not stored") {
    //  This will be reevaluated whenever we support trusted CXSmiles
    //   For now use CachedSmilesMolHolder
    const std::vector<std::pair<std::string, unsigned int>> inferred = {
        {"CC[CH2] |^1:2|", 1},  // agrees with the inference, survives
        {"CC[CH2] |^2:2|", 1},  // asked for 2, inference wins
        {"CC[CH3] |^1:2|", 0},  // asked for 1, inference wins
        {"[CH2] |^1:0|", 2},    // asked for 1, inference wins
    };
    auto molHolder = boost::make_shared<CachedSmilesMolHolder>();
    for (const auto &[smi, expected] : inferred) {
      INFO(smi);
      auto idx = molHolder->addSmiles(smi);
      auto held = molHolder->getMol(idx);
      REQUIRE(held);
      CHECK(totalRadicals(*held) == expected);
    }
  }
}
