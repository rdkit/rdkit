//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <RDGeneral/test.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <algorithm>
#include <set>
#include <sstream>

#include <GraphMol/ChemReactions/Enumerate/CartesianProduct.h>
#include <GraphMol/ChemReactions/Enumerate/EvenSamplePairs.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSample.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSampleAllBBs.h>
#include <GraphMol/ChemReactions/Enumerate/Enumerate.h>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

using namespace RDKit;

#ifdef RDK_USE_BOOST_SERIALIZATION
// for each starting point check to see that the archive
//  starts at the same point
void pickleTest(EnumerationStrategyBase &en, size_t len) {
  boost::shared_ptr<EnumerationStrategyBase> base(en.copy());
  TEST_ASSERT(std::string(base->type()) == std::string(en.type()));

  for (size_t i = 0; i < len; ++i) {
    std::stringstream ss;
    {
      boost::archive::text_oarchive ar(ss);
      ar & base;
    }
    boost::shared_ptr<EnumerationStrategyBase> copy;
    {
      boost::archive::text_iarchive ar(ss);
      ar & copy;
    }
    TEST_ASSERT(std::string(base->type()) == std::string(copy->type()));
    TEST_ASSERT(base->next() == copy->next());
    TEST_ASSERT(base->getPosition() == en.next());
  }
}
#endif

void testSamplers() {
  EnumerationTypes::BBS bbs;
  bbs.resize(3);
  for (int i = 0; i < 10; ++i) {
    bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")));
  }

  for (int i = 0; i < 5; ++i) {
    bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")));
  }

  for (int i = 0; i < 6; ++i) {
    bbs[2].push_back(
        boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")));
  }

  ChemicalReaction rxn;
  CartesianProductStrategy cart;
  cart.initialize(rxn, bbs);
  RandomSampleStrategy rand;
  rand.initialize(rxn, bbs);
  RandomSampleAllBBsStrategy randBBs;
  randBBs.initialize(rxn, bbs);
  EvenSamplePairsStrategy even;
  even.initialize(rxn, bbs);
  std::vector<boost::shared_ptr<EnumerationStrategyBase>> enumerators;
  enumerators.emplace_back(cart.copy());
  enumerators.emplace_back(rand.copy());
  enumerators.emplace_back(randBBs.copy());
  enumerators.emplace_back(even.copy());
#ifdef RDK_USE_BOOST_SERIALIZATION
  for (auto &enumerator : enumerators) {
    TEST_ASSERT(enumerator->getNumPermutations() == 10 * 5 * 6);
    pickleTest(*enumerator, 10 * 5 * 6);
  }
#endif
  // for(auto&& i: enumerators) {
  //  TEST_ASSERT(i->getNumPermutations() == 10*5*6);
  //}
}

void testEvenSamplers() {
  EnumerationTypes::BBS bbs;
  bbs.resize(3);
  boost::uint64_t R1 = 600;
  boost::uint64_t R2 = 50;
  boost::uint64_t R3 = 1000;

  boost::shared_ptr<ROMol> m(SmilesToMol("C=CCN=C=S"));
  boost::shared_ptr<ROMol> m2(SmilesToMol("NCc1ncc(Cl)cc1Br"));
  boost::shared_ptr<ROMol> m3(SmilesToMol("NCCCc1ncc(Cl)cc1Br"));

  for (unsigned long i = 0; i < R1; ++i) {
    bbs[0].push_back(m);
  }

  for (unsigned long i = 0; i < R2; ++i) {
    bbs[1].push_back(m2);
  }

  for (unsigned long i = 0; i < R3; ++i) {
    bbs[2].push_back(m3);
  }

  ChemicalReaction rxn;
  EvenSamplePairsStrategy even;
  even.initialize(rxn, bbs);
  std::cout << even.getNumPermutations() << " " << R1 * R2 * R3 << std::endl;
  TEST_ASSERT(even.getNumPermutations() == R1 * R2 * R3);

  for (size_t i = 0; i < 5000; ++i) {
    even.next();
  }
  even.stats();
}

const char *smiresults[] = {
    "C=CCNC(=S)NCc1ncc(Cl)cc1Br",   "CC=CCNC(=S)NCc1ncc(Cl)cc1Br",
    "C=CCNC(=S)NCCc1ncc(Cl)cc1Br",  "CC=CCNC(=S)NCCc1ncc(Cl)cc1Br",
    "C=CCNC(=S)NCCCc1ncc(Cl)cc1Br", "CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br"};

namespace {
ChemicalReaction *makeEnumerationCacheRxn() {
  return RxnSmartsToChemicalReaction(
      "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);"
      "!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:"
      "2]");
}

EnumerationTypes::BBS makeEnumerationCacheBbs() {
  EnumerationTypes::BBS bbs;
  bbs.resize(2);

  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")));
  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CC=CCN=C=S")));

  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCCc1ncc(Cl)cc1Br")));
  bbs[1].push_back(
      boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")));
  return bbs;
}

EnumerationTypes::BBS makeEnumerationCacheBbsWithNonMatches() {
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbs();
  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CC")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("Clc1ccccc1")));
  return bbs;
}

EnumerationTypes::BBS makeEnumerationCacheBbsWithProtectedAtoms() {
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbsWithNonMatches();
  bbs[0][0]->getAtomWithIdx(0)->setProp(common_properties::_protected, 1);
  bbs[1][0]->getAtomWithIdx(0)->setProp(common_properties::_protected, 1);
  return bbs;
}

std::vector<std::multiset<std::string>> collectReagentSmiles(
    const EnumerationTypes::BBS &bbs) {
  std::vector<std::multiset<std::string>> results;
  results.resize(bbs.size());
  for (size_t i = 0; i < bbs.size(); ++i) {
    for (const auto &mol : bbs[i]) {
      results[i].insert(MolToSmiles(*mol, true));
    }
  }
  return results;
}

std::multiset<std::string> collectEnumerationSmiles(EnumerateLibrary &en) {
  std::multiset<std::string> results;
  for (; (bool)en;) {
    std::vector<std::vector<std::string>> smiles = en.nextSmiles();
    for (const auto &productSet : smiles) {
      TEST_ASSERT(productSet.size() == 1);
      results.insert(productSet[0]);
    }
  }
  return results;
}

std::vector<std::string> collectEnumerationSequence(EnumerateLibrary &en) {
  std::vector<std::string> results;
  for (; (bool)en;) {
    std::vector<std::vector<std::string>> smiles = en.nextSmiles();
    for (const auto &productSet : smiles) {
      TEST_ASSERT(productSet.size() == 1);
      results.push_back(productSet[0]);
    }
  }
  return results;
}

std::set<std::string> collectUniqueSmiles(const std::vector<std::string> &seq) {
  return std::set<std::string>(seq.begin(), seq.end());
}

ChemicalReaction *makeSymmetricDedupeRxn() {
  return RxnSmartsToChemicalReaction("[N;H2:1]>>[N:1]C(=O)C");
}

EnumerationTypes::BBS makeSymmetricDedupBbs() {
  EnumerationTypes::BBS bbs;
  bbs.resize(1);

  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCCN")));
  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CCN")));
  return bbs;
}

std::multiset<std::string> collectDirectProductSmiles(
    const ChemicalReaction &rxn, const EnumerationTypes::BBS &bbs) {
  std::multiset<std::string> results;
  for (size_t i = 0; i < bbs[0].size(); ++i) {
    for (size_t j = 0; j < bbs[1].size(); ++j) {
      MOL_SPTR_VECT reactants;
      reactants.push_back(bbs[0][i]);
      reactants.push_back(bbs[1][j]);
      std::vector<MOL_SPTR_VECT> products = rxn.runReactants(reactants);
      TEST_ASSERT(products.size() == 1);
      TEST_ASSERT(products[0].size() == 1);
      results.insert(MolToSmiles(*products[0][0], true));
    }
  }
  return results;
}
}  // namespace

void testEnumerationCacheOutputUnchanged() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbs();

  EnumerateLibrary cached(*rxn, bbs);
  EnumerateLibrary viaLibrary(*rxn, bbs);

  // drive the non-cached baseline from the library's filtered reagents so the
  // comparison stays valid even if removeNonmatchingReagents drops a reagent
  const std::multiset<std::string> expected =
      collectDirectProductSmiles(*rxn, cached.getReagents());
  const std::multiset<std::string> cachedResults = collectEnumerationSmiles(cached);
  const std::multiset<std::string> viaLibraryResults =
      collectEnumerationSmiles(viaLibrary);

  TEST_ASSERT(cachedResults == expected);
  TEST_ASSERT(viaLibraryResults == expected);
}

void testPrefilterResultsUnchanged() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbsWithNonMatches();

  EnumerationParams params;
  EnumerationTypes::BBS baseline = removeNonmatchingReagents(*rxn, bbs, params);
  EnumerateLibrary library(*rxn, bbs, params);

  TEST_ASSERT(collectReagentSmiles(library.getReagents()) ==
              collectReagentSmiles(baseline));

  params.reagentMaxMatchCount = 1;
  baseline = removeNonmatchingReagents(*rxn, bbs, params);
  EnumerateLibrary cappedLibrary(*rxn, bbs, params);
  TEST_ASSERT(collectReagentSmiles(cappedLibrary.getReagents()) ==
              collectReagentSmiles(baseline));
}

void testPrefilterPrimesSharedCache() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbsWithNonMatches();

  EnumerateLibrary library(*rxn, bbs);
  const EnumerationTypes::BBS &filtered = library.getReagents();
  size_t expectedCacheSize = 0;
  for (const auto &pool : filtered) {
    expectedCacheSize += pool.size();
  }

  TEST_ASSERT(library.getMatchCacheSize() == expectedCacheSize);
}

void testPrefilterFullEnumerationOutputUnchanged() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbsWithNonMatches();

  EnumerateLibrary library(*rxn, bbs);
  const std::multiset<std::string> expected =
      collectDirectProductSmiles(*rxn, library.getReagents());
  const std::multiset<std::string> actual = collectEnumerationSmiles(library);

  TEST_ASSERT(actual == expected);
}

void testPrefilterProtectedAtomFallback() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbsWithProtectedAtoms();

  EnumerationParams params;
  EnumerationTypes::BBS baseline = removeNonmatchingReagents(*rxn, bbs, params);
  EnumerateLibrary library(*rxn, bbs, params);

  // removal behavior must be identical to the 3-arg baseline even when some
  // reagents carry a _protected atom (those take the countMatches fallback)
  TEST_ASSERT(collectReagentSmiles(library.getReagents()) ==
              collectReagentSmiles(baseline));

  // reagents carrying a _protected atom are not primed into the shared cache;
  // only the protected-free survivors take the fast path
  size_t nonProtectedSurvivors = 0;
  for (const auto &pool : library.getReagents()) {
    for (const auto &mol : pool) {
      bool hasProtected = false;
      for (const auto atom : mol->atoms()) {
        if (atom->hasProp(common_properties::_protected)) {
          hasProtected = true;
          break;
        }
      }
      if (!hasProtected) {
        ++nonProtectedSurvivors;
      }
    }
  }
  TEST_ASSERT(library.getMatchCacheSize() == nonProtectedSurvivors);
}

void testEnumerationCacheReuseAcrossReset() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbs();

  EnumerateLibrary en(*rxn, bbs);
  const std::vector<std::string> firstPass = collectEnumerationSequence(en);

  en.reset();

  const std::vector<std::string> secondPass = collectEnumerationSequence(en);
  TEST_ASSERT(firstPass == secondPass);
}

void testDedupOutputUnchangedWhenDisabled() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeSymmetricDedupeRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeSymmetricDedupBbs();

  EnumerationParams params;
  params.dedupeSymmetricMatches = false;

  EnumerateLibrary explicitOff(*rxn, bbs, params);
  EnumerateLibrary defaultOff(*rxn, bbs);

  const std::vector<std::string> explicitOffSequence =
      collectEnumerationSequence(explicitOff);
  const std::vector<std::string> defaultOffSequence =
      collectEnumerationSequence(defaultOff);

  TEST_ASSERT(explicitOffSequence == defaultOffSequence);
  TEST_ASSERT(collectUniqueSmiles(defaultOffSequence).size() <
              defaultOffSequence.size());
}

void testDedupCollapsesSymmetricReagent() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeSymmetricDedupeRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeSymmetricDedupBbs();

  EnumerationParams offParams;
  offParams.dedupeSymmetricMatches = false;
  EnumerateLibrary offLibrary(*rxn, bbs, offParams);
  const std::vector<std::string> offSequence =
      collectEnumerationSequence(offLibrary);

  EnumerationParams onParams;
  onParams.dedupeSymmetricMatches = true;
  EnumerateLibrary onLibrary(*rxn, bbs, onParams);
  const std::vector<std::string> onSequence =
      collectEnumerationSequence(onLibrary);

  TEST_ASSERT(offSequence.size() > onSequence.size());
  TEST_ASSERT(collectUniqueSmiles(offSequence) ==
              collectUniqueSmiles(onSequence));
}

#ifdef RDK_USE_BOOST_SERIALIZATION
void testDedupFlagSerialized() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeSymmetricDedupeRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeSymmetricDedupBbs();

  EnumerationParams params;
  params.dedupeSymmetricMatches = true;
  EnumerateLibrary en(*rxn, bbs, params);
  const std::vector<std::string> baseline = collectEnumerationSequence(en);

  en.reset();
  std::string serialized = en.Serialize();
  EnumerateLibrary reloaded(serialized);

  TEST_ASSERT(reloaded.getDedupeSymmetricMatches());
  const std::vector<std::string> reloadedSequence =
      collectEnumerationSequence(reloaded);
  TEST_ASSERT(reloadedSequence == baseline);
}
#else
void testDedupFlagSerialized() {}
#endif

#ifdef RDK_USE_BOOST_SERIALIZATION
void testEnumerationCacheNotSerialized() {
  boost::scoped_ptr<ChemicalReaction> rxn(makeEnumerationCacheRxn());
  rxn->initReactantMatchers();
  EnumerationTypes::BBS bbs = makeEnumerationCacheBbs();

  EnumerateLibrary en(*rxn, bbs);
  const std::vector<std::string> baseline = collectEnumerationSequence(en);

  en.reset();
  std::string serialized = en.Serialize();
  EnumerateLibrary reloaded(serialized);

  const std::vector<std::string> reloadedSequence = collectEnumerationSequence(reloaded);
  TEST_ASSERT(reloadedSequence == baseline);
}
#else
void testEnumerationCacheNotSerialized() {}
#endif

void testEnumerations() {
  EnumerationTypes::BBS bbs;
  bbs.resize(2);

  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")));
  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CC=CCN=C=S")));

  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCCc1ncc(Cl)cc1Br")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")));

  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
      "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);"
      "!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:"
      "2]");

  {
    EnumerateLibrary en(*rxn, bbs);
    size_t i = 0;
    for (; (bool)en; ++i) {
      std::vector<std::vector<std::string>> res = en.nextSmiles();
      TEST_ASSERT(res.size() == 1);
      TEST_ASSERT(res[0].size() == 1);
      TEST_ASSERT(res[0][0] == smiresults[i]);
      TEST_ASSERT(i <= 6);
    }
    TEST_ASSERT(i == 6);
    // tests reset
    en.resetState();
    i = 0;
    for (; (bool)en; ++i) {
      std::vector<std::vector<std::string>> res = en.nextSmiles();
      TEST_ASSERT(res.size() == 1);
      TEST_ASSERT(res[0].size() == 1);
      TEST_ASSERT(res[0][0] == smiresults[i]);
      TEST_ASSERT(i <= 6);
    }
    TEST_ASSERT(i == 6);
  }

#ifdef RDK_USE_BOOST_SERIALIZATION
  {
    boost::shared_ptr<EnumerateLibrary> en(
        new EnumerateLibrary(*rxn, bbs, RandomSampleStrategy()));

    std::vector<std::vector<std::vector<std::string>>> smir;
    for (size_t j = 0; j < 10; ++j) {
      std::vector<std::vector<std::string>> smiles = en->nextSmiles();
      smir.push_back(smiles);
    }

    en->resetState();

    for (size_t i = 0; i < 1000; ++i) {
      // pickle and unpickle
      std::stringstream ss;
      {
        boost::archive::text_oarchive ar(ss);
        ar & en;
      }
      boost::shared_ptr<EnumerateLibrary> copy;
      {
        boost::archive::text_iarchive ar(ss);
        ar & copy;
      }

      for (size_t j = 0; j < 10; ++j) {
        TEST_ASSERT(en->nextSmiles() == copy->nextSmiles());
      }

      copy->resetState();
      for (size_t j = 0; j < 10; ++j) {
        TEST_ASSERT(smir[j] == copy->nextSmiles());
      }
    }
  }
#endif
  delete rxn;
}

const char *rxndata =
    "$RXN\nUntitled Document-1\n  ChemDraw10291618492D\n\n  3  1\n$MOL\n\n\n\n "
    " 2  1  0  0  0  0  0  0  0  0999 V2000\n    0.4125    0.0000    0.0000 N  "
    " 0  0  0  0  0  0  0  0  0  3  0  0\n   -0.4125    0.0000    0.0000 R2  0 "
    " 0  0  0  0  0  0  0  0  2  0  0\n  1  2  1  0        0\nM  "
    "END\n$MOL\n\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   -0.4125    "
    "0.0000    0.0000 R1  0  0  0  0  0  0  0  0  0  1  0  0\n    0.4125    "
    "0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0     "
    "   0\nM  END\n$MOL\n\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    "
    "0.4125    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  5  0  0\n   "
    "-0.4125    0.0000    0.0000 R4  0  0  0  0  0  0  0  0  0  4  0  0\n  1  "
    "2  1  0        0\nM  END\n$MOL\n\n\n\n 14 15  0  0  0  0  0  0  0  0999 "
    "V2000\n    0.5072   -0.5166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0 "
    " 0\n    0.5072    0.3084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n    1.2949   -0.7616    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n    1.7817   -0.0880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n    1.2967    0.5794    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n    1.5558   -1.5443    0.0000 R1  0  0  0  0  0  0  0  0  0  1  0  "
    "0\n   -0.2073    0.7208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n   -0.9218    0.3083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n   -0.9217   -0.5167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n   -0.2073   -0.9292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  "
    "0\n   -1.6362    0.7208    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  "
    "0\n    1.5452    1.3661    0.0000 N   0  0  0  0  0  0  0  0  0  5  0  "
    "0\n    2.3507    1.5443    0.0000 R4  0  0  0  0  0  0  0  0  0  4  0  "
    "0\n   -2.3507    0.3083    0.0000 R2  0  0  0  0  0  0  0  0  0  2  0  "
    "0\n  1  2  2  0        0\n  1  3  1  0        0\n  3  4  1  0        0\n  "
    "4  5  1  0        0\n  5  2  1  0        0\n  3  6  1  0        0\n  2  7 "
    " 1  0        0\n  7  8  2  0        0\n  8  9  1  0        0\n  9 10  2  "
    "0        0\n 10  1  1  0        0\n  8 11  1  0        0\n 12 13  1  0    "
    "    0\n 11 14  1  0        0\n 12  5  1  0        0\nM  END\n";

void testInsaneEnumerations() {
  EnumerationTypes::BBS bbs;
  bbs.resize(3);

  ChemicalReaction *rxn2 = RxnBlockToChemicalReaction(rxndata);
  // RxnOps::sanitizeRxn(*rxn2, MolOps::AdjustQueryParameters());
  MatchVectType tvect;

  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CCNCC")));
  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCC")));
  std::cerr << "0,0 "
            << (int)SubstructMatch(*bbs[0][0].get(),
                                   *rxn2->getReactants()[0].get(), tvect)
            << std::endl;
  std::cerr << "0,1 "
            << (int)SubstructMatch(*bbs[0][1].get(),
                                   *rxn2->getReactants()[0].get(), tvect)
            << std::endl;

  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("ClC1CCC1")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("ClC1CCC1Cl")));
  std::cerr << "1,0 "
            << (int)SubstructMatch(*bbs[1][0].get(),
                                   *rxn2->getReactants()[1].get(), tvect)
            << std::endl;
  std::cerr << "1,1 "
            << (int)SubstructMatch(*bbs[1][1].get(),
                                   *rxn2->getReactants()[1].get(), tvect)
            << std::endl;

  bbs[2].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CCNCC")));
  bbs[2].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCC")));
  std::cerr << "2,0 "
            << (int)SubstructMatch(*bbs[2][0].get(),
                                   *rxn2->getReactants()[2].get(), tvect)
            << std::endl;
  std::cerr << "2,1 "
            << (int)SubstructMatch(*bbs[2][1].get(),
                                   *rxn2->getReactants()[2].get(), tvect)
            << std::endl;

  {
    ChemicalReaction *rxn = RxnBlockToChemicalReaction(rxndata);
    RxnOps::sanitizeRxn(*rxn, MolOps::AdjustQueryParameters());
    EnumerationParams ThereCanBeOnlyOne;
    ThereCanBeOnlyOne.reagentMaxMatchCount = 1;
    EnumerationTypes::BBS bbs2 =
        removeNonmatchingReagents(*rxn, bbs, ThereCanBeOnlyOne);
    TEST_ASSERT(bbs2[0].size() == 1);
    TEST_ASSERT(bbs2[1].size() == 1);
    TEST_ASSERT(bbs2[2].size() == 1);

    delete rxn;
  }
  delete rxn2;
}

#ifdef RDK_USE_BOOST_SERIALIZATION
void testGithub1657() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github #1657: EnumerateLibrary with "
                          "initFromString called twice doesn't clear the "
                          "reaction"
                       << std::endl;
  EnumerationTypes::BBS bbs;
  bbs.resize(2);

  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("C=CCN=C=S")));
  bbs[0].push_back(boost::shared_ptr<ROMol>(SmilesToMol("CC=CCN=C=S")));

  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCc1ncc(Cl)cc1Br")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCCc1ncc(Cl)cc1Br")));
  bbs[1].push_back(boost::shared_ptr<ROMol>(SmilesToMol("NCCCc1ncc(Cl)cc1Br")));

  ChemicalReaction *rxn = RxnSmartsToChemicalReaction(
      "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);"
      "!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:"
      "2]");

  {  // we'll also test that initFromString() works at all
    EnumerateLibrary en1(*rxn, bbs);
    std::stringstream sstr;
    en1.toStream(sstr);
    EnumerateLibrary en2;
    en2.initFromString(sstr.str());

    EnumerateLibrary en3;
    en3.initFromString(sstr.str());
    bool ok = false;
    try {
      en3.initFromString(sstr.str());
    } catch (const ValueErrorException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  delete rxn;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}
#else
void testGithub1657() {}
#endif

int main() {
  RDLog::InitLogs();
#if 1
  testPrefilterResultsUnchanged();
  testPrefilterPrimesSharedCache();
  testPrefilterFullEnumerationOutputUnchanged();
  testPrefilterProtectedAtomFallback();
  testEnumerationCacheOutputUnchanged();
  testEnumerationCacheReuseAcrossReset();
  testDedupOutputUnchangedWhenDisabled();
  testDedupCollapsesSymmetricReagent();
  testDedupFlagSerialized();
  testEnumerationCacheNotSerialized();
  testSamplers();
  testEvenSamplers();
  testEnumerations();
  testInsaneEnumerations();
#endif
  testGithub1657();
}
