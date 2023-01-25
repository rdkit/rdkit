//  Copyright (C) 2020 Gareth Jones, Glysade LLC
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
#include <string>
#include <vector>
#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#define private public
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/RGroupDecomposition/RGroupDecompData.h>
#include <GraphMol/RGroupDecomposition/RGroupGa.h>

using namespace RDKit;

void testCoresLabelledProperly() {
  // Tests for an error in RGroupDecompositionParameters::prepareCore where
  // the same unindexed label could be mistakenly given to multiple rgroups
  // in different cores

  // See https://github.com/rdkit/rdkit/pull/3565

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test cores labelled properly" << std::endl;

  std::vector<std::string> coreSmi = {"N1([*:1])CCN([*:2])CC1",
                                      "C1(O[*:1])CCC(O[*:2])CC1",
                                      "C1([*:1])CCC([*:2])CC1"};

  std::vector<ROMOL_SPTR> cores;
  for (const auto &s : coreSmi) {
    cores.emplace_back(SmartsToMol(s));
  }

  RGroupDecompositionParameters params;
  params.alignment = RGroupCoreAlignment::None;
  params.scoreMethod = FingerprintVariance;
  RGroupDecomposition decomposition(cores, params);

  auto data = decomposition.data;
  std::set<int> rlabels;

  for (const auto &core : data->cores) {
    auto mol = core.second.core;
    for (const auto atom : mol->atoms()) {
      if (atom->hasProp(RLABEL)) {
        int rlabel = atom->getProp<int>(RLABEL);
        if (rlabel < 0) {
          TEST_ASSERT(rlabels.find(rlabel) == rlabels.end());
          rlabels.insert(rlabel);
        }
      }
    }
  }
}

std::pair<int, RData> makeRData(int attachment, std::vector<int> attachments,
                                const std::string &smiles) {
  auto rData = boost::make_shared<RGroupData>();
  auto mol = SmilesToMol(smiles);
  auto frags = MolOps::getMolFrags(*mol);
  for (auto &frag : frags) {
    rData->add(frag, attachments);
  }
  delete mol;
  std::pair<int, RData> pair(attachment, rData);
  return pair;
}

std::pair<int, RData> makeRData(int attachment, const std::string &smiles) {
  std::vector<int> attachments{attachment};
  return makeRData(attachment, attachments, smiles);
}

void testRingMatching3Score() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test scoring function for RingMatching3- see GitHub ##3924"
      << std::endl;

  R_DECOMP decomp1Mol1 = {makeRData(-4, "*[H]"), makeRData(-2, "*[H]"),
                          makeRData(-1, "*[H]"), makeRData(1, "*C([H])([H])C")};
  R_DECOMP decomp1Mol2 = {makeRData(-4, "*[H]"), makeRData(-3, "*[H]"),
                          makeRData(-2, "*[H]"), makeRData(-1, "*[H]"),
                          makeRData(1, "*C([H])([H])I")};
  R_DECOMP decomp1Mol3 = {makeRData(-4, "*[H]"), makeRData(-2, "*[H]"),
                          makeRData(-1, "*[H]"), makeRData(1, "*C([H])([H])F")};
  RGroupMatch match1Mol1(0, 0, decomp1Mol1, nullptr);
  RGroupMatch match1Mol2(0, 0, decomp1Mol2, nullptr);
  RGroupMatch match1Mol3(0, 0, decomp1Mol3, nullptr);
  std::set<int> labels{-4, -3, -2, -1, 1};
  std::vector<RGroupMatch> matches1Mol1{match1Mol1};
  std::vector<RGroupMatch> matches1Mol2{match1Mol2};
  std::vector<RGroupMatch> matches1Mol3{match1Mol3};
  std::vector<std::vector<RGroupMatch>> allMatches1 = {
      matches1Mol1, matches1Mol2, matches1Mol3};
  std::vector<size_t> permutation{0, 0, 0};

  R_DECOMP decomp2Mol1 = {makeRData(-4, "*[H]"), makeRData(-3, "*C([H])([H])C"),
                          makeRData(-2, "*[H]"), makeRData(-1, "*[H]")};
  R_DECOMP decomp2Mol2 = {makeRData(-4, "*[H]"), makeRData(-3, "*[H]"),
                          makeRData(-2, "*[H]"), makeRData(-1, "*[H]"),
                          makeRData(1, "*C([H])([H])I")};
  R_DECOMP decomp2Mol3 = {makeRData(-4, "*[H]"), makeRData(-2, "*[H]"),
                          makeRData(-1, "*[H]"), makeRData(1, "*C([H])([H])F")};
  RGroupMatch match2Mol1(0, 1, decomp2Mol1, nullptr);
  RGroupMatch match2Mol2(0, 0, decomp2Mol2, nullptr);
  RGroupMatch match2Mol3(0, 0, decomp2Mol3, nullptr);
  std::vector<RGroupMatch> matches2Mol1{match2Mol1};
  std::vector<RGroupMatch> matches2Mol2{match2Mol2};
  std::vector<RGroupMatch> matches2Mol3{match2Mol3};
  std::vector<std::vector<RGroupMatch>> allMatches2 = {
      matches2Mol1, matches2Mol2, matches2Mol3};

  RGroupScorer scorer;
  auto test1 = scorer.matchScore(permutation, allMatches1, labels);
  auto test2 = scorer.matchScore(permutation, allMatches2, labels);

  // expect test1 to have better score than test2 since all halogens are on R1

  TEST_ASSERT(test1 > test2);

  auto testFp1 = fingerprintVarianceScore(permutation, allMatches1, labels);
  auto testFp2 = fingerprintVarianceScore(permutation, allMatches2, labels);

  TEST_ASSERT(testFp1 > testFp2);
}

void testGithub3746() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test GA falls over to exhaustive on simple system"
                       << std::endl;
  const std::vector<ROMOL_SPTR> cores{"c1([*:1])c([*:2])c([*:3])ccc1"_smiles,
                                      "c1([*:1])c([*:2])c([*:3])cnc1"_smiles};
  const std::vector<const char *> smilesData{"c1(CO)cc(CN)ccc1",
                                             "c1(CO)cc(CN)cnc1"};

  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = true;
  params.removeHydrogensPostMatch = true;
  params.matchingStrategy = GA;
  params.removeHydrogensPostMatch = true;

  RGroupDecomposition decomposition(cores, params);
  size_t i = 0;
  for (const auto &smi : smilesData) {
    ROMOL_SPTR mol(static_cast<ROMol *>(SmilesToMol(smi)));
    TEST_ASSERT(decomposition.add(*mol) == static_cast<int>(i++));
  }

  auto data = decomposition.data;
  RGroupGa ga(*data);
  auto numberPermutations = ga.numberPermutations();

  TEST_ASSERT(numberPermutations == 4);
  // criteria for exhaustive search instead of GA
  TEST_ASSERT(numberPermutations < ga.getPopsize() * 100);
}

int main() {
  RDLog::InitLogs();
  boost::logging::disable_logs("rdApp.debug");

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group Decomposition Internals\n";

  testRingMatching3Score();
  testGithub3746();
  testCoresLabelledProperly();

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
