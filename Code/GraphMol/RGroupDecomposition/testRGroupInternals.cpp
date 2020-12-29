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

int main() {
  RDLog::InitLogs();
  boost::logging::disable_logs("rdApp.debug");

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group Decomposition Internals\n";

  testCoresLabelledProperly();

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
