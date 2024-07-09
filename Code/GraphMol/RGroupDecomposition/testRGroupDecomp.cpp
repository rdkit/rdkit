//  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
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
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/RGroupDecomposition/RGroupDecompData.h>
#include <GraphMol/RGroupDecomposition/RGroupUtils.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/Exceptions.h>
#include <boost/tokenizer.hpp>
#include <regex>

// #define DEBUG

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

using namespace RDKit;

#ifdef DEBUG
const bool DOASSERT = false;
#else
const bool DOASSERT = true;
#endif

typedef std::vector<std::unique_ptr<ROMol>> UMOLS;
#define UPTR(m) std::unique_ptr<ROMol>(m)

void CHECK_RGROUP(RGroupRows::const_iterator &it, const std::string &expected,
                  ROMol *mol = nullptr, bool doassert = DOASSERT) {
  std::ostringstream str;
  int i = 0;
  std::unique_ptr<ROMol> res;
  for (auto rgroups = it->begin(); rgroups != it->end(); ++rgroups, ++i) {
    if (i) {
      str << " ";
      if (mol) {
        res = molzip(*res, *rgroups->second.get());
      }
    } else if (mol) {
      res = std::unique_ptr<ROMol>(new ROMol(*rgroups->second.get()));
    }
    // rlabel:smiles
    str << rgroups->first << ":" << MolToSmiles(*rgroups->second.get());
  }
  std::string result = str.str();

  if (expected != result) {
    std::cerr << "Expected: '" << expected << "'" << std::endl;
    std::cerr << "Got:      '" << result << "'" << std::endl;
  }

  if (doassert) {
    TEST_ASSERT(result == expected)
    if (mol) {
      auto smi1 = MolToSmiles(*res);
      auto smi2 = MolToSmiles(*mol);
      TEST_ASSERT(smi1 == smi2)
    }
  }
}

void DUMP_RGROUP(RGroupRows::const_iterator &it, std::string &result) {
  std::ostringstream str;

  for (const auto &rgroups : *it) {
    // rlabel:smiles
    str << rgroups.first << ":" << MolToSmiles(*rgroups.second.get(), true)
        << " ";
  }
  std::cerr << str.str() << std::endl;
  result = str.str();
}

const char *symdata[5] = {"c1(Cl)ccccc1", "c1c(Cl)cccc1", "c1cccc(Cl)c1",
                          "c1cc(Cl)ccc1", "c1ccc(Cl)cc1"};

void testSymmetryMatching(RGroupScore scoreMethod = Match) {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "test rgroup decomp symmetry matching with score method "
      << scoreMethod << std::endl;

  UMOLS mols;
  RWMol *core = SmilesToMol("c1ccccc1");
  RGroupDecompositionParameters params;
  params.scoreMethod = scoreMethod;
  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 5; ++i) {
    ROMol *mol = SmilesToMol(symdata[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    mols.push_back(UPTR(mol));
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, "Core:c1ccc([*:1])cc1 R1:Cl[*:1]", mols[i].get());
  }
  delete core;
}

void testGaSymmetryMatching(RGroupScore scoreMethod) {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "test rgroup decomp symmetry matching using GA with scoring method "
      << scoreMethod << std::endl;

  UMOLS mols;
  RWMol *core = SmilesToMol("c1ccccc1");
  RGroupDecompositionParameters params;
  params.matchingStrategy = GA;
  params.scoreMethod = scoreMethod;
  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 5; ++i) {
    ROMol *mol = SmilesToMol(symdata[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    mols.push_back(UPTR(mol));
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, "Core:c1ccc([*:1])cc1 R1:Cl[*:1]", mols[i].get());
  }
  delete core;
}

void testGaBatch() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "test rgroup decomp symmetry matching using GA with parallel runs"
      << std::endl;

  UMOLS mols;
  RWMol *core = SmilesToMol("c1ccccc1");
  RGroupDecompositionParameters params;
  params.matchingStrategy = GA;
  params.scoreMethod = FingerprintVariance;
  params.gaNumberRuns = 3;
  params.gaParallelRuns = true;

  std::stringstream sstrm;
  rdWarningLog->SetTee(sstrm);
  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 5; ++i) {
    ROMol *mol = SmilesToMol(symdata[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    mols.push_back(UPTR(mol));
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  bool isParallelGaEnabled =
      (sstrm.str().find("This RDKit build does not enable GA parallel runs") ==
       std::string::npos);
#ifdef RDK_TEST_MULTITHREADED
  TEST_ASSERT(isParallelGaEnabled);
#else
  TEST_ASSERT(!isParallelGaEnabled);
#endif
  rdWarningLog->ClearTee();

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
    CHECK_RGROUP(it, "Core:c1ccc([*:1])cc1 R1:Cl[*:1]", mols[i].get());
  }
  delete core;
}

const char *matchRGroupOnlyData[] = {
    "c1(Cl)ccccc1", "c1c(Cl)cccc1",    "c1cc(Cl)ccc1",
    "c1ccc(Cl)cc1", "c1c(Cl)cccc(I)1",
};

void testRGroupOnlyMatching() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp rgroup only matching"
                       << std::endl;

  UMOLS mols;
  RWMol *core = SmilesToMol("c1ccccc1[1*]");
  RGroupDecompositionParameters params;
  params.labels = IsotopeLabels;
  params.onlyMatchAtRGroups = true;

  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 5; ++i) {
    ROMol *mol = SmilesToMol(matchRGroupOnlyData[i]);
    int res = decomp.add(*mol);
    if (i < 4) {
      TEST_ASSERT(res == i);
      mols.push_back(UPTR(mol));
    } else {
      TEST_ASSERT(res == -1);
      delete mol;
    }
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, "Core:c1ccc([*:1])cc1 R1:Cl[*:1]", mols[i].get());
  }
  delete core;
}

const char *ringData[3] = {"c1cocc1", "c1c[nH]cc1", "c1cscc1"};

const char *ringDataRes[3] = {"Core:c1ccoc1", "Core:c1cc[nH]c1",
                              "Core:c1ccsc1"};

void testRingMatching() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp ring matching" << std::endl;

  UMOLS mols;
  RWMol *core = SmilesToMol("c1ccc[1*]1");
  RGroupDecompositionParameters params;
  params.labels = IsotopeLabels;

  auto exceptionThrown = false;
  try {
    RGroupDecomposition decompError(*core, params);
  } catch (ValueErrorException &) {
    exceptionThrown = true;
  }
  TEST_ASSERT(exceptionThrown);

  params.allowNonTerminalRGroups = true;
  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 3; ++i) {
    ROMol *mol = SmilesToMol(ringData[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    mols.push_back(UPTR(mol));
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  auto cols = decomp.getRGroupsAsColumns();
  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    // Ring rgroups not supported by molzip yet.
    CHECK_RGROUP(it, ringDataRes[i]);
  }
  delete core;
}

const char *ringData2[3] = {"c1cocc1CCl", "c1c[nH]cc1CI", "c1cscc1CF"};

const char *ringDataRes2[3] = {"Core:c1cc(C[*:2])co1 R2:Cl[*:2]",
                               "Core:c1cc(C[*:2])c[nH]1 R2:I[*:2]",
                               "Core:c1cc(C[*:2])cs1 R2:F[*:2]"};

void testRingMatching2() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp full ring dummy core"
                       << std::endl;

  RWMol *core = SmartsToMol("*1***[*:1]1C[*:2]");
  RGroupDecompositionParameters params;
  params.allowNonTerminalRGroups = true;

  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 3; ++i) {
    ROMol *mol = SmilesToMol(ringData2[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    delete mol;
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, ringDataRes2[i]);
  }
  delete core;
}

const char *ringData3[3] = {"c1cocc1CCl", "c1c[nH]cc1CI", "c1cscc1CF"};

const char *ringDataRes3[3] = {"Core:c1cc([*:1])co1 R1:ClC[*:1]",
                               "Core:c1cc([*:1])c[nH]1 R1:IC[*:1]",
                               "Core:c1cc([*:1])cs1 R1:FC[*:1]"};

void testRingMatching3() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp full ring dummy core"
                       << std::endl;

  RWMol *core = SmartsToMol("*1***[*:1]1");
  // RWMol *core = SmartsToMol("*1****1");

  std::vector<RGroupScore> matchtypes{Match, FingerprintVariance};
  for (auto match : matchtypes) {
    RGroupDecompositionParameters params;
    // This test is currently failing using the default scoring method (the
    // halogens are not all in the same group)
    params.scoreMethod = match;
    params.allowNonTerminalRGroups = true;

    RGroupDecomposition decomp(*core, params);
    for (int i = 0; i < 3; ++i) {
      ROMol *mol = SmilesToMol(ringData3[i]);
      int res = decomp.add(*mol);
      delete mol;
      TEST_ASSERT(res == i);
    }

    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    // All Cl's should be labeled with the same rgroup
    int i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
         ++it, ++i) {
      CHECK_RGROUP(it, ringDataRes3[i]);
    }
  }
  delete core;
}

const char *coreSmi[] = {
    "C1CCNC(Cl)CC1", "C1CC(Cl)NCCC1", "C1CCNC(I)CC1", "C1CC(I)NCCC1",

    "C1CCSC(Cl)CC1", "C1CC(Cl)SCCC1", "C1CCSC(I)CC1", "C1CC(I)SCCC1",

    "C1CCOC(Cl)CC1", "C1CC(Cl)OCCC1", "C1CCOC(I)CC1", "C1CC(I)OCCC1"};

const char *coreSmiRes[] = {
    "Core:C1CCNC([*:1])CC1 R1:Cl[*:1]", "Core:C1CCNC([*:1])CC1 R1:Cl[*:1]",
    "Core:C1CCNC([*:1])CC1 R1:I[*:1]",  "Core:C1CCNC([*:1])CC1 R1:I[*:1]",
    "Core:C1CCSC([*:1])CC1 R1:Cl[*:1]", "Core:C1CCSC([*:1])CC1 R1:Cl[*:1]",
    "Core:C1CCSC([*:1])CC1 R1:I[*:1]",  "Core:C1CCSC([*:1])CC1 R1:I[*:1]",
    "Core:C1CCOC([*:1])CC1 R1:Cl[*:1]", "Core:C1CCOC([*:1])CC1 R1:Cl[*:1]",
    "Core:C1CCOC([*:1])CC1 R1:I[*:1]",  "Core:C1CCOC([*:1])CC1 R1:I[*:1]"};

void testMultiCore() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test multi core" << std::endl;
  std::vector<ROMOL_SPTR> cores;
  cores.emplace_back(SmartsToMol("C1CCNCCC1"));
  cores.emplace_back(SmilesToMol("C1CCOCCC1"));
  cores.emplace_back(SmilesToMol("C1CCSCCC1"));
  UMOLS mols;
  RGroupDecomposition decomp(cores);
  for (unsigned int i = 0; i < sizeof(coreSmi) / sizeof(const char *); ++i) {
    ROMol *mol = SmilesToMol(coreSmi[i]);
    unsigned int res = decomp.add(*mol);
    mols.push_back(UPTR(mol));
    TEST_ASSERT(res == i);
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    // molzip doesn't support double attachments yet (it probably should)
    CHECK_RGROUP(it, coreSmiRes[i]);
  }
}

void testGithub1550() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "test Github #1550: Kekulization error from R-group decomposition"
      << std::endl;

  RWMol *core = SmilesToMol("O=c1oc2ccccc2cc1");
  RGroupDecompositionParameters params;

  RGroupDecomposition decomp(*core, params);
  const char *smilesData[3] = {"O=c1cc(Cn2ccnc2)c2ccc(Oc3ccccc3)cc2o1",
                               "O=c1oc2ccccc2c(Cn2ccnc2)c1-c1ccccc1",
                               "COc1ccc2c(Cn3cncn3)cc(=O)oc2c1"};
  for (int i = 0; i < 3; ++i) {
    ROMol *mol = SmilesToMol(smilesData[i]);
    int res = decomp.add(*mol);
    delete mol;
    TEST_ASSERT(res == i);
  }

  decomp.process();
  RGroupColumns groups = decomp.getRGroupsAsColumns();

  RWMol *coreRes = (RWMol *)groups["Core"][0].get();
  TEST_ASSERT(coreRes->getNumAtoms() == 14);
  MolOps::Kekulize(*coreRes);
  RWMol *rg2 = (RWMol *)groups["R2"][0].get();
  TEST_ASSERT(rg2->getNumAtoms() == 7);
  MolOps::Kekulize(*rg2);

  delete core;
}

void testRemoveHs() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test remove sidechain Hs" << std::endl;

  RWMol *core = SmilesToMol("O=c1oc2ccccc2cc1");

  {
    RGroupDecompositionParameters params;
    RGroupDecomposition decomp(*core, params);
    const char *smilesData[3] = {"O=c1cc(Cn2ccnc2)c2ccc(Oc3ccccc3)cc2o1",
                                 "O=c1oc2ccccc2c(Cn2ccnc2)c1-c1ccccc1",
                                 "COc1ccc2c(Cn3cncn3)cc(=O)oc2c1"};
    for (int i = 0; i < 3; ++i) {
      ROMol *mol = SmilesToMol(smilesData[i]);
      int res = decomp.add(*mol);
      delete mol;
      TEST_ASSERT(res == i);
    }

    decomp.process();
    RGroupColumns groups = decomp.getRGroupsAsColumns();
    RWMol *rg2 = (RWMol *)groups["R2"][0].get();
    TEST_ASSERT(rg2->getNumAtoms() == 7);
  }
  {
    RGroupDecompositionParameters params;
    params.removeHydrogensPostMatch = false;
    RGroupDecomposition decomp(*core, params);
    const char *smilesData[3] = {"O=c1cc(Cn2ccnc2)c2ccc(Oc3ccccc3)cc2o1",
                                 "O=c1oc2ccccc2c(Cn2ccnc2)c1-c1ccccc1",
                                 "COc1ccc2c(Cn3cncn3)cc(=O)oc2c1"};
    for (int i = 0; i < 3; ++i) {
      ROMol *mol = SmilesToMol(smilesData[i]);
      int res = decomp.add(*mol);
      delete mol;
      TEST_ASSERT(res == i);
    }

    decomp.process();
    RGroupColumns groups = decomp.getRGroupsAsColumns();
    RWMol *rg2 = (RWMol *)groups["R2"][0].get();
    TEST_ASSERT(rg2->getNumAtoms() == 12);
  }
  delete core;
}

void testGitHubIssue1705() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "test grouping substituents in chunks as large as possible"
      << std::endl;
#if 1
  {
    RWMol *core = SmilesToMol("Oc1ccccc1");
    RGroupDecompositionParameters params;

    RGroupDecomposition decomp(*core, params);
    const char *smilesData[5] = {"Oc1ccccc1", "Oc1c(F)cccc1", "Oc1ccccc1F",
                                 "Oc1c(F)cc(N)cc1", "Oc1ccccc1Cl"};
    for (int i = 0; i < 5; ++i) {
      ROMol *mol = SmilesToMol(smilesData[i]);
      int res = decomp.add(*mol);
      delete mol;
      TEST_ASSERT(res == i);
    }

    decomp.process();
    std::stringstream ss;
    RGroupColumns groups = decomp.getRGroupsAsColumns();
    for (auto &column : groups) {
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    delete core;
    std::string expected = R"RES(Rgroup===Core
Oc1ccc([*:1])cc1[*:2]
Oc1ccc([*:1])cc1[*:2]
Oc1ccc([*:1])cc1[*:2]
Oc1ccc([*:1])cc1[*:2]
Oc1ccc([*:1])cc1[*:2]
Rgroup===R1
[H][*:1]
[H][*:1]
[H][*:1]
N[*:1]
[H][*:1]
Rgroup===R2
[H][*:2]
F[*:2]
F[*:2]
F[*:2]
Cl[*:2]
)RES";
#ifdef DEBUG
    if (ss.str() != expected) {
      std::cerr << __LINE__ << " ERROR got\n"
                << ss.str() << "\nexpected\n"
                << expected << std::endl;
    }
#else
    TEST_ASSERT(ss.str() == expected);
#endif
  }
#endif
  // std::cerr<<"n\n\n\n\n\n--------------------------------------------------------------\n\n\n\n\n";
  {
    RWMol *core = SmilesToMol("Cc1ccccc1");
    RGroupDecompositionParameters params;

    RGroupDecomposition decomp(*core, params);
    std::vector<std::string> smilesData = {"c1ccccc1C", "Fc1ccccc1C",
                                           "c1cccc(F)c1C", "Fc1cccc(F)c1C"};
    for (const auto &smi : smilesData) {
      ROMol *mol = SmilesToMol(smi);
      decomp.add(*mol);
      delete mol;
    }

    decomp.process();
    std::stringstream ss;
    RGroupColumns groups = decomp.getRGroupsAsColumns();
    for (auto &column : groups) {
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    delete core;
    std::string expected = R"RES(Rgroup===Core
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Rgroup===R1
[H][*:1]
[H][*:1]
[H][*:1]
F[*:1]
Rgroup===R2
[H][*:2]
F[*:2]
F[*:2]
F[*:2]
)RES";
#ifdef DEBUG
    if (ss.str() != expected) {
      std::cerr << __LINE__ << " ERROR got\n"
                << ss.str() << "\nexpected\n"
                << expected << std::endl;
    }
#else
    TEST_ASSERT(ss.str() == expected);
#endif
  }
}

void testMatchOnlyAtRgroupHs() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test matching only rgroups but allows Hs"
                       << std::endl;

  RWMol *core = SmilesToMol("*OCC");
  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp(*core, params);
  const char *smilesData[2] = {"OCC", "COCC"};
  for (auto &i : smilesData) {
    ROMol *mol = SmilesToMol(i);
    decomp.add(*mol);
    delete mol;
  }
  decomp.process();

  std::stringstream ss;
  RGroupColumns groups = decomp.getRGroupsAsColumns();
  for (auto &column : groups) {
    ss << "Rgroup===" << column.first << std::endl;
    for (auto &rgroup : column.second) {
      ss << MolToSmiles(*rgroup) << std::endl;
    }
  }
  delete core;
  TEST_ASSERT(
      ss.str() ==
      "Rgroup===Core\nCCO[*:1]\nCCO[*:1]\nRgroup===R1\n[H][*:1]\nC[*:1]\n");
}

void testGithub2332() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test github #2332: RGroupDecomposition: addHs() "
                          "call should set coords "
                       << std::endl;
  auto core = "*OCC"_smiles;
  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp(*core, params);
  std::string chains[2] = {
      R"CTAB(
  Mrv1810 03291913362D          

  4  3  0  0  0  0            999 V2000
    2.0625   -0.7145    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    1.2375   -0.7145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
M  END
)CTAB",
      R"CTAB(
  Mrv1810 03291913362D          

  3  2  0  0  0  0            999 V2000
    1.2375   -0.7145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
)CTAB"};
  for (const auto &chain : chains) {
    ROMol *mol = MolBlockToMol(chain);
    decomp.add(*mol);
    delete mol;
  }
  decomp.process();

  std::stringstream ss;
  RGroupColumns groups = decomp.getRGroupsAsColumns();
  auto &r1 = groups["R1"];
  TEST_ASSERT(r1.size() == 2);
  TEST_ASSERT(r1[1]->getAtomWithIdx(0)->getAtomicNum() == 1);
  auto conf = r1[1]->getConformer();
  TEST_ASSERT(!feq(conf.getAtomPos(0).x, 0.0));
  TEST_ASSERT(!feq(conf.getAtomPos(0).y, 0.0));
  TEST_ASSERT(feq(conf.getAtomPos(0).z, 0.0));
}

void testSDFGRoupMultiCoreNoneShouldMatch() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "testSDFGRoupMultiCoreNoneShouldMatch" << std::endl;
  std::string sdcores = R"CTAB(
  Mrv1813 05061918272D          

 13 14  0  0  0  0            999 V2000
   -1.1505    0.0026    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1505   -0.8225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4360   -1.2350    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2784   -0.8225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2784    0.0026    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4360    0.4151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9354    0.2575    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4202   -0.4099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9354   -1.0775    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    0.9907   -1.2333    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4360    1.2373    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2784    1.6497    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -3.2452   -0.4098    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  9  2  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  4 10  1  0  0  0  0
  5  6  1  0  0  0  0
  6 11  1  0  0  0  0
  7  8  1  0  0  0  0
  8 13  1  0  0  0  0
  8  9  1  0  0  0  0
 11 12  1  0  0  0  0
M  RGP  3  10   1  12   2  13   3
M  END
$$$$

  Mrv1813 05061918272D          

 13 14  0  0  0  0            999 V2000
    6.9524    0.1684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9524   -0.6567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6668   -1.0692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3813   -0.6567    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3813    0.1684    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.6668    0.5809    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1674    0.4233    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    5.6827   -0.2441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1674   -0.9117    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    9.0935   -1.0675    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    7.6668    1.4031    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    8.3813    1.8155    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    4.8576   -0.2440    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  9  2  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  4 10  1  0  0  0  0
  5  6  1  0  0  0  0
  6 11  1  0  0  0  0
  7  8  1  0  0  0  0
  8 13  1  0  0  0  0
  8  9  1  0  0  0  0
 11 12  1  0  0  0  0
M  RGP  3  10   1  12   2  13   3
M  END
$$$$)CTAB";
  std::string sdmols = R"CTAB(
  Mrv1813 05061918322D          

 15 17  0  0  0  0            999 V2000
    0.1742    0.6899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8886    0.2774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8886   -0.5476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1742   -0.9601    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1742   -1.7851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8886   -2.1976    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8886   -3.0226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1742   -3.4351    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5403   -3.0226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3249   -3.2775    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8099   -2.6101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3249   -1.9426    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5403   -2.1976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5403   -0.5476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5403    0.2774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1 15  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  4 14  1  0  0  0  0
  5  6  1  0  0  0  0
  5 13  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
  9 13  1  0  0  0  0
 10 11  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  1  0  0  0  0
 14 15  1  0  0  0  0
M  END
$$$$

  Mrv1813 05061918322D          

 14 15  0  0  0  0            999 V2000
    6.4368    0.3002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7223   -0.1123    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.7223   -0.9373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4368   -1.3498    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.4368   -2.1748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7223   -2.5873    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.0078   -2.1748    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2232   -2.4297    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    3.7383   -1.7623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2232   -1.0949    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9683   -0.3102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1613   -0.1387    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5203    0.3029    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0078   -1.3498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7 14  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  1  0  0  0  0
 10 14  1  0  0  0  0
 11 13  1  0  0  0  0
 11 12  1  0  0  0  0
M  END
$$$$

  Mrv1813 05061918322D          

 14 15  0  0  0  0            999 V2000
    0.8289   -7.9643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1144   -8.3768    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1144   -9.2018    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8289   -9.6143    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8289  -10.4393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1144  -10.8518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000  -10.4393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3847  -10.6942    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8696  -10.0268    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3847   -9.3593    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6396   -8.5747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4466   -8.4032    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0876   -7.9616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000   -9.6143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7 14  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  1  0  0  0  0
 10 14  1  0  0  0  0
 11 13  1  0  0  0  0
 11 12  1  0  0  0  0
M  END
$$$$

  Mrv1813 05061918322D          

 12 13  0  0  0  0            999 V2000
    5.3295   -8.1871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5844   -7.4025    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.0995   -6.7351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5844   -6.0676    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.3690   -6.3226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0835   -5.9101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0835   -5.0851    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.7980   -6.3226    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.7980   -7.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5124   -7.5601    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0835   -7.5601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3690   -7.1476    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2 12  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5 12  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  6  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
 11 12  1  0  0  0  0
M  END
$$$$)CTAB";

  {
    RGroupDecompositionParameters params;
    params.onlyMatchAtRGroups = true;
    params.removeHydrogensPostMatch = true;
    std::vector<ROMOL_SPTR> cores;

    {
      SDMolSupplier sdsup;
      sdsup.setData(sdcores);
      while (!sdsup.atEnd()) {
        cores.emplace_back(sdsup.next());
      }
    }

    RGroupDecomposition decomp(cores, params);

    {
      SDMolSupplier sdsup;
      sdsup.setData(sdmols);

      while (!sdsup.atEnd()) {
        ROMol *mol = sdsup.next();
        TEST_ASSERT(mol);
        int addedIndex = decomp.add(*mol);
        TEST_ASSERT(addedIndex == -1);  // none should match
        delete mol;
      }
    }
  }
  {
    RGroupDecompositionParameters params;
    params.onlyMatchAtRGroups = false;
    params.removeHydrogensPostMatch = true;
    std::vector<ROMOL_SPTR> cores;

    {
      SDMolSupplier sdsup;
      sdsup.setData(sdcores);
      while (!sdsup.atEnd()) {
        cores.emplace_back(sdsup.next());
      }
    }

    RGroupDecomposition decomp(cores, params);

    {
      SDMolSupplier sdsup;
      sdsup.setData(sdmols);

      while (!sdsup.atEnd()) {
        ROMol *mol = sdsup.next();
        TEST_ASSERT(mol);
        decomp.add(*mol);
        delete mol;
      }
    }

    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();

    /*
     * Working on issue Github5613 these core smiles for the 3rd target no
     * longer works, however the replacement is for the same core structure
     *
    const char *expected[4] = {
        "Core:C1C([*:6])C2C(N([*:2])[*:4])NC([*:1])NC2N1[*:5] "
        "R1:[H][*:1] R2:C(CC[*:2])CC[*:4] R4:C(CC[*:2])CC[*:4] R5:[H][*:5] "
        "R6:[H][*:6]",
        "Core:C1SC2NC([*:1])NC(N([*:2])[*:4])C2C1[*:6] "
        "R1:[H][*:1] R2:C[*:2] R4:[H][*:4] R6:CC(C)[*:6]",
        "Core:C1SC2CC([*:1])NC(N([*:2])[*:4])C2C1[*:6] "
        "R1:[H][*:1] R2:C[*:2] R4:[H][*:4] R6:CC(C)[*:6]",
        "Core:C1C2C(C(N([*:2])[*:4])NC1[*:1])N([*:6])CN2[*:5] R1:O[*:1] "
        "R2:[H][*:2] R4:[H][*:4] R5:C[*:5] R6:[H][*:6]"};
     */
    const char *expected[4] = {
        "Core:C1C([*:6])C2C(N([*:2])[*:4])NC([*:1])NC2N1[*:5] "
        "R1:[H][*:1] R2:C(CC[*:2])CC[*:4] R4:C(CC[*:2])CC[*:4] R5:[H][*:5] "
        "R6:[H][*:6]",
        "Core:C1SC2NC([*:1])NC(N([*:2])[*:4])C2C1[*:6] "
        "R1:[H][*:1] R2:C[*:2] R4:[H][*:4] R6:CC(C)[*:6]",
        "Core:C1C2SCC([*:6])C2C(N([*:2])[*:4])NC1[*:1] "
        "R1:[H][*:1] R2:C[*:2] R4:[H][*:4] R6:CC(C)[*:6]",
        "Core:C1C2C(C(N([*:2])[*:4])NC1[*:1])N([*:6])CN2[*:5] R1:O[*:1] "
        "R2:[H][*:2] R4:[H][*:4] R5:C[*:5] R6:[H][*:6]"};
    int i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
         ++it, ++i) {
      TEST_ASSERT(i < 4);
      // molzip doesn't support double attachment points yet
      CHECK_RGROUP(it, expected[i]);
    }
  }
}

void testRowColumnAlignmentProblem() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test a problem with row-column alignment"
                       << std::endl;
  std::vector<std::string> csmiles = {"c1c([*:1])cncn1", "c1c([*:1])cccn1"};
  std::vector<ROMOL_SPTR> cores;
  for (auto smi : csmiles) {
    cores.emplace_back(SmilesToMol(smi));
  }

  std::vector<std::string> msmiles = {"c1c(F)cccn1", "c1c(F)cncn1",
                                      "c1c(Cl)cccn1"};
  std::vector<std::unique_ptr<RWMol>> mols;
  for (auto smi : msmiles) {
    mols.emplace_back(SmilesToMol(smi));
  }

  {
    RGroupDecomposition decomp(cores);
    for (const auto &mol : mols) {
      TEST_ASSERT(decomp.add(*mol) >= 0);
    }
    decomp.process();

    auto rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == mols.size());
    // dump rgroups
    const char *expected[] = {"Core:c1cncc([*:1])c1 R1:F[*:1]",
                              "Core:c1ncc([*:1])cn1 R1:F[*:1]",
                              "Core:c1cncc([*:1])c1 R1:Cl[*:1]"};

    int i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
         ++it, ++i) {
      CHECK_RGROUP(it, expected[i], mols[i].get());
    }

    for (const auto &row : rows) {
      TEST_ASSERT(row.count("Core") == 1);
      TEST_ASSERT(row.count("R1") == 1);
    }
    TEST_ASSERT(rows[0].count("R2") == 0);
    TEST_ASSERT(rows[2].count("R2") == 0);
    TEST_ASSERT(rows[1].count("R2") == 0);

    auto cols = decomp.getRGroupsAsColumns();
    auto &core = cols["Core"];
    TEST_ASSERT(core.size() == 3);
    auto &R1 = cols["R1"];
    TEST_ASSERT(R1.size() == 3);
    for (const auto &rg : R1) {
      TEST_ASSERT(rg);
      TEST_ASSERT(rg->getNumAtoms());
    }
    auto &R2 = cols["R2"];
    TEST_ASSERT(R2.size() == 0);
  }
}

void testSymmetryIssues() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group symmetry issues\n";

  auto m1 = "c1c(F)cccn1"_smiles;
  auto m2 = "c1c(Cl)c(C)c(Br)cn1"_smiles;
  auto m3 = "c1c(O)cccn1"_smiles;
  auto m4 = "c1c(Br)c(C)c(F)cn1"_smiles;
  auto core = "c1c([*:1])c([*:2])ccn1"_smiles;
  {
    RGroupDecomposition decomp(*core);
    decomp.add(*m1);
    decomp.add(*m2);
    decomp.add(*m3);
    decomp.add(*m4);
    decomp.process();
    std::stringstream ss;
    auto groups = decomp.getRGroupsAsColumns();
    std::set<std::string> r_labels;
    for (auto &column : groups) {
      r_labels.insert(column.first);
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    // We want three groups added, fluorine as R1 and bromine as R3

    TEST_ASSERT(r_labels == std::set<std::string>({"Core", "R1", "R2", "R3"}));
    TEST_ASSERT(groups.size() == 4);

    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
c1ncc([*:3])c([*:2])c1[*:1]
c1ncc([*:3])c([*:2])c1[*:1]
c1ncc([*:3])c([*:2])c1[*:1]
c1ncc([*:3])c([*:2])c1[*:1]
Rgroup===R1
F[*:1]
Cl[*:1]
O[*:1]
F[*:1]
Rgroup===R2
[H][*:2]
C[*:2]
[H][*:2]
C[*:2]
Rgroup===R3
[H][*:3]
Br[*:3]
[H][*:3]
Br[*:3]
)RES");
  }
  {
    // repeat that without symmetrization (testing #3224)
    // Still three groups added, but bromine and fluorine
    // are not aligned between R1 and R3
    RGroupDecompositionParameters ps;
    ps.matchingStrategy = RDKit::NoSymmetrization;
    RGroupDecomposition decomp(*core, ps);
    decomp.add(*m1);
    decomp.add(*m2);
    decomp.add(*m3);
    decomp.add(*m4);
    decomp.process();
    std::stringstream ss;
    auto groups = decomp.getRGroupsAsColumns();
    std::set<std::string> r_labels;
    for (auto &column : groups) {
      r_labels.insert(column.first);
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    TEST_ASSERT(r_labels == std::set<std::string>({"Core", "R1", "R2", "R3"}));
    TEST_ASSERT(groups.size() == 4);

    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
c1ncc([*:3])c([*:2])c1[*:1]
c1ncc([*:3])c([*:2])c1[*:1]
c1ncc([*:3])c([*:2])c1[*:1]
c1ncc([*:3])c([*:2])c1[*:1]
Rgroup===R1
F[*:1]
Cl[*:1]
O[*:1]
Br[*:1]
Rgroup===R2
[H][*:2]
C[*:2]
[H][*:2]
C[*:2]
Rgroup===R3
[H][*:3]
Br[*:3]
[H][*:3]
F[*:3]
)RES");
  }
}

void testSymmetryPerformance() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group symmetry performance\n";
  boost::logging::disable_logs("rdApp.warning");

  std::string smis =
      R"DATA(CN(C)Cc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CNc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
Cc1cc2cc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NC4CCN(C)CC4)c([N+](=O)[O-])c3)ccc2[nH]1
Cc1cc2cc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)ccc2[nH]1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccccc1Cl
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(Cl)c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(Cl)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc([N+](=O)[O-])c1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(CO)c1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2ccccc2Cl)cc1[N+](=O)[O-]
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(Cl)c1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(Cl)cc1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc(Cl)c2)cc1[N+](=O)[O-]
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2ccc(Cl)cc2)cc1[N+](=O)[O-]
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc3c2ccn3C)cc1[N+](=O)[O-]
CC(=O)Nc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
Nc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
Nc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
COc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CN(C)c1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
N#Cc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
Cc1nc2ccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCC4CCOCC4)c([N+](=O)[O-])c3)cc2s1
Cc1nc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)ccc2s1
Cc1nc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN(C)C)c([N+](=O)[O-])c3)ccc2s1
CN(C)C(=O)CCc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)C(=O)Cc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)CCCc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)CCc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)C(=O)c1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)Cc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc(N3CCOCC3)c2)cc1[N+](=O)[O-]
Cc1nc(C)c(-c2cccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)c2)s1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cc(Cl)cc(Cl)c3)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(Cl)c3)cc2[N+](=O)[O-])CC1
CN(C)CCOc1ccc(-c2ccc(Cl)cc2)c(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NC5CCN(C)CC5)c([N+](=O)[O-])c4)c(Oc4cccc(Cl)c4)c3)CC2)c1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)OC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN(C)CCOc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccc(N)c(Cl)c3)cc2[N+](=O)[O-])CC1
CC(C)N1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Br)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CCCC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
Cc1n[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)c12
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(F)c3F)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(Br)c3)cc2[N+](=O)[O-])CC1
CCN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN1C(C)(C)CC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1(C)C
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NC5CCN(C6CCOCC6)CC5)c([N+](=O)[O-])c4)c(Oc4cccc(F)c4F)c3)CC2)=C(c2ccc(Cl)cc2)C1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cc(F)c4[nH]ccc4c3)cc2[N+](=O)[O-])CC1
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCCCN5CCOCC5)c([N+](=O)[O-])c4)c(Oc4cccc(F)c4F)c3)CC2)=C(c2ccc(Cl)cc2)C1
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCCCN5CCOCC5)c([N+](=O)[O-])c4)c(Oc4ccc(N)c(Cl)c4)c3)CC2)=C(c2ccc(Cl)cc2)C1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3cc(CCN4CCCC4)ccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(Cl)c1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(Cl)c3Cl)cc2[N+](=O)[O-])CC1
Cc1n[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NC4CCN(C)CC4)c([N+](=O)[O-])c3)c12
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CCCCC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(C(F)(F)F)c3)cc2[N+](=O)[O-])CC1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc3c2CCC(=O)N3)cc1[N+](=O)[O-]
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2S(=O)(=O)C(F)(F)F)CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cc(Cl)ccc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccc(F)cc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)C5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
Cc1c[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)c12
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(C(F)(F)F)c3Cl)cc2[N+](=O)[O-])CC1
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NC5CCN(C6CC6)CC5)c([N+](=O)[O-])c4)c(Oc4ccccc4Cl)c3)CC2)=C(c2ccc(Cl)cc2)C1
Cc1c[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NC4CCN(C)CC4)c([N+](=O)[O-])c3)c12
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCCCN5CCOCC5)c([N+](=O)[O-])c4)c(Oc4cc(Cl)ccc4Cl)c3)CC2)=C(c2ccc(Cl)cc2)C1
Cn1ccc2c(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)cccc21
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(N2CCOCC2)c1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2ccc3[nH]cc(CCC(=O)N4CCOCC4)c3c2)cc1[N+](=O)[O-]
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(OCc2ccccc2)c1
N#Cc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc2[nH]cc(CCC(=O)N3CCOCC3)c2c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc2[nH]cc(CCCN3CCOCC3)c2c1
CN(C)Cc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(-n2ccnc2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)cc1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc([N+](=O)[O-])c1
CCN(Cc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1)C(=O)OC(C)(C)C
CCN(Cc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1)C(=O)OC(C)(C)C
CCNCc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
CCNCc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CC(=O)Nc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
CC(C)(C)OC(=O)Nc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccccc1-c1ccccc1
CC(C)(C)OC(=O)Nc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(-c2ccccc2)c1
CN(C)CCc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(OCc2ccccc2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(N2CCOCC2)c1
Cc1nc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCC4CCOCC4)c([N+](=O)[O-])c3)ccc2s1
CC(C)(C)OC(=O)N1CCN(c2cccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCC4CCOCC4)c([N+](=O)[O-])c3)c2)CC1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc(OCc3ccccc3)c2)cc1[N+](=O)[O-]
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(OCc2ccccc2)c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(OCCN2CCOCC2)cc1
O=C1CCc2c(cccc2Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)N1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(OCc2ccccc2)cc1
CC(C)(C)OC(=O)N1CCN(c2ccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)cc2)CC1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(-c2ccncc2)c1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(-c2ccncc2)cc1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(-c2cccnc2)cc1
CN(C)C(=O)COc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
Cn1cnc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)ccc21)DATA";
  std::vector<ROMOL_SPTR> ms;
  boost::char_separator<char> sep(" \n");
  tokenizer tokens(smis, sep);
  for (tokenizer::iterator token = tokens.begin(); token != tokens.end();
       ++token) {
    std::string smi = *token;
    RWMol *m = SmilesToMol(smi);
    TEST_ASSERT(m);
    ms.push_back(ROMOL_SPTR(m));
  }
  auto core = "O=C(NS(=O)(=O)c1ccccc1)c1ccccc1Oc1ccccc1"_smiles;

  {
    std::cerr << "iterative" << std::endl;
    RGroupDecompositionParameters ps = RGroupDecompositionParameters();
    ps.timeout = 0.1;
    RGroupDecomposition decomp(*core, ps);
    bool ok = false;
    int ndone = -1;
    try {
      for (auto m : ms) {
        decomp.add(*m);
        ++ndone;
      }
    } catch (const std::runtime_error &) {
      ok = true;
    }
    TEST_ASSERT(ok);
    TEST_ASSERT(ndone >= 0);
  }
  {
    RGroupDecompositionParameters ps = RGroupDecompositionParameters();
    ps.timeout = 2.0;
    std::cerr << "bulk" << std::endl;
    std::vector<ROMOL_SPTR> cores;
    cores.push_back(ROMOL_SPTR(new ROMol(*core)));
    RGroupRows rows;
    bool ok = false;
    try {
      auto res = RGroupDecompose(cores, ms, rows, nullptr, ps);
      RDUNUSED_PARAM(res);
    } catch (const std::runtime_error &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    RGroupDecompositionParameters ps = RGroupDecompositionParameters();
#ifdef NDEBUG
    ps.timeout = 5.0;
#else
    ps.timeout = 25.0;
#endif
    ps.matchingStrategy = RDKit::NoSymmetrization;
    std::cerr << "bulk, no symmetry" << std::endl;
    std::vector<ROMOL_SPTR> cores;
    cores.push_back(ROMOL_SPTR(new ROMol(*core)));
    RGroupRows rows;
    bool ok = true;
    try {
      auto res = RGroupDecompose(cores, ms, rows, nullptr, ps);
      RDUNUSED_PARAM(res);
    } catch (const std::runtime_error &) {
      ok = false;
    }
    TEST_ASSERT(ok);
  }
  boost::logging::enable_logs("rdApp.warning");
}

void testScorePermutations() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing permutation scoring function\n";

  {
    auto core = "Cc1ccccc1"_smiles;
    std::vector<RWMOL_SPTR> mols{"c1ccccc1C"_smiles, "Fc1ccccc1C"_smiles,
                                 "c1cccc(F)c1C"_smiles, "Fc1cccc(F)c1C"_smiles};
    RGroupDecomposition decomp(*core);
    for (auto &m : mols) {
      decomp.add(*m);
    }
    decomp.process();
    std::stringstream ss;
    auto groups = decomp.getRGroupsAsColumns();
    std::set<std::string> r_labels;
    for (auto &column : groups) {
      r_labels.insert(column.first);
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    TEST_ASSERT(r_labels == std::set<std::string>({"Core", "R1", "R2"}));
    TEST_ASSERT(groups.size() == 3);
    std::string expected = R"RES(Rgroup===Core
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Rgroup===R1
[H][*:1]
[H][*:1]
[H][*:1]
F[*:1]
Rgroup===R2
[H][*:2]
F[*:2]
F[*:2]
F[*:2]
)RES";
#ifdef DEBUG
    if (ss.str() != expected) {
      std::cerr << __LINE__ << " ERROR got\n"
                << ss.str() << "\nexpected\n"
                << expected << std::endl;
    }
#else
    TEST_ASSERT(ss.str() == expected);
#endif
  }
  {
    auto core = "Cc1ccccc1"_smiles;
    std::vector<RWMOL_SPTR> mols{
        "c1(Cl)cccc(Cl)c1C"_smiles, "Fc1cccc(Cl)c1C"_smiles,
        "Clc1cccc(F)c1C"_smiles, "Fc1cccc(F)c1C"_smiles};
    RGroupDecomposition decomp(*core);
    for (auto &m : mols) {
      decomp.add(*m);
    }
    decomp.process();
    std::stringstream ss;
    auto groups = decomp.getRGroupsAsColumns();
    std::set<std::string> r_labels;
    for (auto &column : groups) {
      r_labels.insert(column.first);
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    TEST_ASSERT(r_labels == std::set<std::string>({"Core", "R1", "R2"}));
    TEST_ASSERT(groups.size() == 3);
    std::string expected = R"RES(Rgroup===Core
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Cc1c([*:1])cccc1[*:2]
Rgroup===R1
Cl[*:1]
Cl[*:1]
Cl[*:1]
F[*:1]
Rgroup===R2
Cl[*:2]
F[*:2]
F[*:2]
F[*:2]
)RES";
#ifdef DEBUG
    if (ss.str() != expected) {
      std::cerr << __LINE__ << " ERROR got\n"
                << ss.str() << "\nexpected\n"
                << expected << std::endl;
    }
#else
    TEST_ASSERT(ss.str() == expected);
#endif
  }
  {
    auto core = "O1C([*:1])([*:2])CCC1"_smiles;
    std::vector<RWMOL_SPTR> mols{"NC1CCCO1"_smiles,
                                 "OC1CCCO1"_smiles,
                                 "SC1CCCO1"_smiles,
                                 "O1C2(CN2)CCC1"_smiles,
                                 "OC1(P)CCCO1"_smiles,
                                 "NC1(O)CCCO1"_smiles,
                                 "CC1(N)CCCO1"_smiles,
                                 "CCC1(C)CCCO1"_smiles,
                                 "CCC1(C(C)C)CCCO1"_smiles,
                                 "CCC1(Cc2ccccc2)CCCO1"_smiles,
                                 "OCC1(Cc2ccccc2)CCCO1"_smiles,
                                 "CCCC1(CO)CCCO1"_smiles,
                                 "CCCCC1(CC(C)C)CCCO1"_smiles};
    RGroupDecompositionParameters params;
    params.removeHydrogensPostMatch = true;
    params.onlyMatchAtRGroups = true;
    RGroupDecomposition decomp(*core, params);
    for (auto &m : mols) {
      auto res = decomp.add(*m);
      TEST_ASSERT(res != -1);
    }
    decomp.process();
    std::stringstream ss;
    auto groups = decomp.getRGroupsAsColumns();
    std::set<std::string> r_labels;
    for (auto &column : groups) {
      r_labels.insert(column.first);
      ss << "Rgroup===" << column.first << std::endl;
      for (auto &rgroup : column.second) {
        ss << MolToSmiles(*rgroup) << std::endl;
      }
    }
    TEST_ASSERT(r_labels == std::set<std::string>({"Core", "R1", "R2"}));
    TEST_ASSERT(groups.size() == 3);
    std::string expected = R"RES(Rgroup===Core
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
C1COC([*:1])([*:2])C1
Rgroup===R1
N[*:1]
O[*:1]
S[*:1]
C(N[*:1])[*:2]
P[*:1]
O[*:1]
N[*:1]
C[*:1]
CC(C)[*:1]
c1ccc(C[*:1])cc1
c1ccc(C[*:1])cc1
OC[*:1]
CC(C)C[*:1]
Rgroup===R2
[H][*:2]
[H][*:2]
[H][*:2]
C(N[*:1])[*:2]
O[*:2]
N[*:2]
C[*:2]
CC[*:2]
CC[*:2]
CC[*:2]
OC[*:2]
CCC[*:2]
CCCC[*:2]
)RES";

    if (ss.str() != expected) {
      std::cerr << __LINE__ << " ERROR got\n"
                << ss.str() << "\nexpected\n"
                << expected << std::endl;
    }

    TEST_ASSERT(ss.str() == expected);
  }
}

void testMultiCorePreLabelled() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test multi core pre-labelled" << std::endl;

  struct MultiCoreRGD {
    static void test(
        const std::vector<ROMOL_SPTR> &cores,
        RGroupDecompositionParameters &params,
        const std::vector<std::string> &expectedLabels,
        const std::vector<std::string> &expectedRows,
        const std::vector<std::vector<std::string>> &expectedItems) {
      std::vector<ROMOL_SPTR> mols{"CNC(=O)C1=CN=CN1CC"_smiles,
                                   "Fc1ccc2ccc(Br)nc2n1"_smiles};
      params.removeHydrogensPostMatch = true;
      params.onlyMatchAtRGroups = true;
      RGroupDecomposition decomp(cores, params);
      unsigned int i = 0;
      for (const auto &m : mols) {
        unsigned int res = decomp.add(*m);
        TEST_ASSERT(res == i++);
      }
      decomp.process();
      RGroupRows rows = decomp.getRGroupsAsRows();
      i = 0;
      for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
           ++it, ++i) {
        CHECK_RGROUP(it, expectedRows[i]);
      }
      RGroupColumns groups = decomp.getRGroupsAsColumns();
      i = 0;
      TEST_ASSERT(groups.size() <= expectedLabels.size());
      for (const auto &pair : groups) {
#ifdef DEBUG
        if (pair.first != expectedLabels[i]) {
          std::cerr << __LINE__ << " ERROR: Expected " << expectedLabels[i]
                    << ", got " << pair.first << std::endl;
        }
#else
        TEST_ASSERT(pair.first == expectedLabels[i]);
#endif
        unsigned int j = 0;
        for (const auto &item : pair.second) {
#ifdef DEBUG
          if (expectedItems[i][j] != MolToSmiles(*item)) {
            std::cerr << __LINE__ << " ERROR: Expected " << expectedItems[i][j]
                      << ", got " << MolToSmiles(*item) << std::endl;
          }
#else
          TEST_ASSERT(expectedItems[i][j] == MolToSmiles(*item));
#endif
          ++j;
        }
        ++i;
      }
    }
  };

  std::string sdcores = R"CTAB(
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
    1.1100   -1.3431    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.5225   -0.6286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9705   -0.0156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2168   -0.3511    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3029   -1.1716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1419    0.7914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5289    1.3431    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9266    1.0463    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4976    0.0613    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  1  0
  1  5  2  0
  3  6  1  0
  6  7  2  0
  6  8  1  0
  4  9  1  0
M  RGP  2   8   1   9   2
V    8 *
V    9 *
M  END
$$$$

     RDKit          2D

 12 13  0  0  0  0  0  0  0  0999 V2000
   -6.5623    0.3977    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -5.8478   -0.0147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1333    0.3977    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4188   -0.0147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4188   -0.8397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1333   -1.2522    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8478   -0.8397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7044   -1.2522    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7044    0.3977    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9899   -0.0147    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9899   -0.8397    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2754    0.3978    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  2  3  2  0
  2  7  1  0
  9 10  2  0
 10 11  1  0
  8 11  2  0
  8  5  1  0
  4  9  1  0
 10 12  1  0
  1  2  1  0
M  RGP  2   1   2  12   1
V    1 *
V   12 *
M  END
$$$$
)CTAB";
  std::vector<ROMOL_SPTR> cores;
  SDMolSupplier sdsup;
  sdsup.setData(sdcores);
  while (!sdsup.atEnd()) {
    cores.emplace_back(sdsup.next());
  }

  std::vector<std::string> expectedRows{
      "Core:O=C(c1cncn1[*:2])[*:1] R1:CN[*:1] R2:CC[*:2]",
      "Core:c1cc2ccc([*:2])nc2nc1[*:1] R1:F[*:1] R2:Br[*:2]"};

  std::vector<std::vector<std::string>> expectedItems{
      {"O=C(c1cncn1[*:2])[*:1]", "c1cc2ccc([*:2])nc2nc1[*:1]"},
      {"CN[*:1]", "F[*:1]"},
      {"CC[*:2]", "Br[*:2]"},
  };

  std::vector<std::string> expectedLabels{"Core", "R1", "R2"};

  RGroupDecompositionParameters params;

  // test pre-labelled with MDL R-group labels, autodetect
  params.labels = AutoDetect;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with MDL R-group labels, no autodetect
  params.labels = MDLRGroupLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with MDL R-group labels, autodetect, no MCS alignment
  params.labels = AutoDetect;
  params.alignment = NoAlignment;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);

  // Reading from a MDL molblock also sets isotopic labels, so no need
  // to set them again; we only clear MDL R-group labels
  for (auto &core : cores) {
    for (auto a : core->atoms()) {
      if (a->hasProp(common_properties::_MolFileRLabel)) {
        a->clearProp(common_properties::_MolFileRLabel);
      }
    }
  }
  // test pre-labelled with isotopic labels, autodetect
  params.labels = AutoDetect;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with isotopic labels, no autodetect
  params.labels = IsotopeLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with isotopic labels, autodetect, no MCS alignment
  params.labels = AutoDetect;
  params.alignment = NoAlignment;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);

  for (auto &core : cores) {
    for (auto a : core->atoms()) {
      auto iso = a->getIsotope();
      if (iso) {
        a->setAtomMapNum(iso);
        a->setIsotope(0);
      }
    }
  }
  // test pre-labelled with atom map labels, autodetect
  params.labels = AutoDetect;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with atom map labels, no autodetect
  params.labels = AtomMapLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with atom map labels, autodetect, no MCS alignment
  params.labels = AutoDetect;
  params.alignment = NoAlignment;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);

  for (auto &core : cores) {
    for (auto a : core->atoms()) {
      if (a->getAtomMapNum()) {
        a->setAtomMapNum(0);
      }
    }
  }
  // test pre-labelled with dummy atom labels, autodetect

  expectedRows = {"Core:O=C(c1cncn1[*:2])[*:1] R1:CN[*:1] R2:CC[*:2]",
                  "Core:c1cc2ccc([*:2])nc2nc1[*:1] R1:Br[*:1] R2:F[*:2]"};

  expectedItems = {{"O=C(c1cncn1[*:2])[*:1]", "c1cc2ccc([*:2])nc2nc1[*:1]"},
                   {"CN[*:1]", "Br[*:1]"},
                   {"CC[*:2]", "F[*:2]"}};

  params.labels = AutoDetect;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
  // test pre-labelled with dummy atom labels, no autodetect
  params.labels = DummyAtomLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRows,
                     expectedItems);
}

void testCoreWithRGroupAdjQuery() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test core with query atom adjacent to R-group"
                       << std::endl;
  std::string sdcore_query = R"CTAB(
     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
   -3.6689   -0.8582    0.0000 R#  0  0  0  0  0  1  0  0  0  0  0  0
   -2.2421   -1.3211    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1279   -0.3169    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4403    1.1502    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3261    2.1543    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1007    1.6914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4132    0.2243    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2989   -0.7798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8400   -0.2386    0.0000 Q   0  0  0  0  0  0  0  0  0  0  0  0
    3.1525   -1.7057    0.0000 R#  0  0  0  0  0  1  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  7  9  1  0
  9 10  1  0
  8  3  1  0
M  RGP  2   1   1  10   2
M  END
$$$$
)CTAB";
  std::string sdcore_noquery =
      std::regex_replace(sdcore_query, std::regex("Q  "), "O  ");
  auto mol = "CNc1cccc(c1)OC1CCC1"_smiles;
  for (const auto &sdcore : {sdcore_query, sdcore_noquery}) {
    SDMolSupplier sdsup;
    sdsup.setData(sdcore);
    ROMOL_SPTR core(sdsup.next());
    RGroupDecompositionParameters params;
    params.removeHydrogensPostMatch = true;
    params.onlyMatchAtRGroups = true;
    RGroupDecomposition decomp(*core, params);
    TEST_ASSERT(decomp.add(*mol) == 0);
    TEST_ASSERT(decomp.process());
    RGroupColumns groups = decomp.getRGroupsAsColumns();
    TEST_ASSERT(groups.size() == 3);
    TEST_ASSERT(groups.find("R1") != groups.end());
    TEST_ASSERT(groups.find("R2") != groups.end());
    TEST_ASSERT(MolToSmiles(*groups.at("R1")[0]) == "C[*:1]");
    TEST_ASSERT(MolToSmiles(*groups.at("R2")[0]) == "C1CC([*:2])C1");

    auto rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    RGroupRows::const_iterator it = rows.begin();
    std::string expected(
        "Core:c1cc(N[*:1])cc(O[*:2])c1 R1:C[*:1] R2:C1CC([*:2])C1");
    CHECK_RGROUP(it, expected);
  }
}

void testMultipleCoreRelabellingIssues() {
  // This test fixes 2 issues with relabelling groups
  // Firstly, a new R group which appeared in a later core could have it's label
  // assigned to an unindexed group in a previous core
  // Secondly, a user defined r group which is not part of the decomposition
  // could have it's index assigned to an unindexed group.

  // See https://github.com/rdkit/rdkit/pull/3565

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test relabelling issues in multiple core decomp"
                       << std::endl;

  std::vector<std::shared_ptr<ROMol>> molecules;
  {
    std::fstream fh;
    std::string rdBase(getenv("RDBASE"));
    fh.open(rdBase + "/Docs/Notebooks/compounds.txt", std::ios::in);
    std::string line;
    getline(fh, line);

    while (getline(fh, line)) {
      int pos = line.find_last_of("\t");
      auto smiles = line.substr(pos + 1);
      std::shared_ptr<ROMol> mol(SmilesToMol(smiles));
      molecules.push_back(mol);
      if (molecules.size() == 30) {
        break;
      }
    }
  }

  std::vector<std::string> smi{
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])S2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])CS2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])CC2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])CO2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])C2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])C2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)=C([*:3])S2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])S2",
      "O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])O2",
      "O=C1C([*:2])([*:1])C([*:6])([*:5])N1"};
  std::vector<RGroupScore> matchtypes{Match, FingerprintVariance};
  for (auto match : matchtypes) {
    std::vector<ROMOL_SPTR> cores;
    for (const auto &s : smi) {
      cores.emplace_back(SmartsToMol(s));
    }

    RGroupDecompositionParameters params;
    params.scoreMethod = match;
    RGroupDecomposition decomposition(cores, params);
    for (auto &mol : molecules) {
      decomposition.add(*mol);
    }

    decomposition.process();
    const auto &columns = decomposition.getRGroupsAsColumns();
    TEST_ASSERT(columns.size() == 7u);
    for (auto &col : columns) {
      TEST_ASSERT(30U == col.second.size());
    }
  }
}

void testUnprocessedMapping() {
  // Tests a bug that results in an unprocessed mapping Invariant violation
  // The cause of the error is an rgroup mistakenly identified as containing
  // only hydrogens in a multicore decomp

  // See https://github.com/rdkit/rdkit/pull/3565

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test unprocessed mapping error in multiple core decomp" << std::endl;
  std::vector<std::string> structureSmi = {
      "Cn1nc(-c2ccccc2)cc1OC1CCC(OC2CCN(C(=O)OC(C)(C)C)CC2)CC1",
      "Cc1cccc(N2CCN(c3ncnc(Nc4ccc(S(C)(=O)=O)nc4C)c3F)[C@@H](C)C2)c1",
      "CC(C)(C)OC(=O)N1CCC(OCC2CCCCC2COc2ccc(Br)nn2)CC1",
      "CC(C)(C)OC(=O)N1CCC(CO[C@H]2CC[C@@H](c3ccc(S(C)(=O)=O)nc3)CC2)CC1",
      "CCCCCCCC/"
      "C=C\\CCCCCCCC(=O)Oc1ccc2c(c1)CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12"};
  std::vector<std::string> coreSmi = {"N1([*:1])CCN([*:2])CC1",
                                      "C1(O[*:1])CCC(O[*:2])CC1",
                                      "C1([*:1])CCC([*:2])CC1"};

  std::vector<RGroupScore> matchtypes{Match, FingerprintVariance};
  for (auto match : matchtypes) {
    std::vector<ROMOL_SPTR> cores;
    for (const auto &s : coreSmi) {
      cores.emplace_back(SmartsToMol(s));
    }

    RGroupDecompositionParameters params;
    params.scoreMethod = match;
    RGroupDecomposition decomposition(cores, params);
    for (auto &smi : structureSmi) {
      auto mol = SmilesToMol(smi);
      decomposition.add(*mol);
      delete mol;
    }

    auto result = decomposition.processAndScore();
    TEST_ASSERT(result.success);
    TEST_ASSERT(result.score != -1.0);
  }
}

void testGeminalRGroups() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test core with geminal R-groups" << std::endl;
  std::string core_ctab = R"CTAB(
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
   -0.6026    1.2267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171    0.8142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171   -0.0108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6026   -0.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1118   -0.0108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1118    0.8142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0714    1.7839    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0506    1.8398    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  2  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  0  0  0  0
  1  8  1  0  0  0  0
M  RGP  2   7   5   8   6
M  END
)CTAB";
  ROMOL_SPTR core(MolBlockToMol(core_ctab));
  const std::vector<const char *> smilesData{"C1CCCCC12CC2", "C1CCCCC1(C)C",
                                             "C1CCCCC1(Cl)Br"};

  // the test should yield the same results irrespective of the permutations
  // across the two parameters
  for (auto matchAtRGroup = 0; matchAtRGroup < 2; ++matchAtRGroup) {
    for (auto mdlRGroupLabels = 0; mdlRGroupLabels < 2; ++mdlRGroupLabels) {
      RGroupDecompositionParameters params;
      if (matchAtRGroup) {
        params.labels = MDLRGroupLabels;
      }
      if (mdlRGroupLabels) {
        params.labels = MDLRGroupLabels;
      }
      RGroupDecomposition decomp(*core, params);
      for (const auto &smi : smilesData) {
        ROMol *mol = SmilesToMol(smi);
        TEST_ASSERT(decomp.add(*mol) != -1);
        delete mol;
      }
      TEST_ASSERT(decomp.process());
      auto rows = decomp.getRGroupsAsRows();
      const std::vector<const char *> res{
          "Core:C1CCC([*:5])([*:6])CC1 R5:C(C[*:6])[*:5] R6:C(C[*:6])[*:5]",
          "Core:C1CCC([*:5])([*:6])CC1 R5:C[*:5] R6:C[*:6]",
          "Core:C1CCC([*:5])([*:6])CC1 R5:Br[*:5] R6:Cl[*:6]"};
      TEST_ASSERT(rows.size() == res.size());
      size_t i = 0;
      for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
           ++it) {
        CHECK_RGROUP(it, res.at(i++));
      }
    }
  }
}

void testMatchOnAnyAtom() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test core with matching at any atom" << std::endl;
  std::string core_ctab = R"CTAB(
  Mrv2008 01192109352D

 12 12  0  0  0  0            999 V2000
    3.7389   -3.2028    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5640   -3.2028    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    4.9764   -2.4884    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    4.5640   -1.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7389   -1.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3265   -2.4884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5015   -2.4884    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    3.3265   -3.9172    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    5.8014   -2.4884    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    2.0890   -1.7739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5015   -1.0595    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2640   -1.7739    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  5  6  8  0  0  0  0
  4  5  8  0  0  0  0
  3  4  8  0  0  0  0
  2  3  8  0  0  0  0
  1  2  8  0  0  0  0
  1  6  8  0  0  0  0
  6  7  1  0  0  0  0
  1  8  1  0  0  0  0
  3  9  1  0  0  0  0
  7 10  1  0  0  0  0
 10 11  2  0  0  0  0
 10 12  1  0  0  0  0
M  RGP  3   8   2   9   3  12   1
M  END
)CTAB";
  ROMOL_SPTR core(MolBlockToMol(core_ctab));
  const std::vector<const char *> smilesData{""};
  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = false;

  RGroupDecomposition decomp(*core, params);
  auto mol = "O=C(NC1CCN(c2ncccc2[N+](=O)[O-])CC1)c1cc(Cl)c(Cl)[nH]1"_smiles;
  int res = decomp.add(*mol);
  TEST_ASSERT(res == 0);

  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  RGroupRows::const_iterator it = rows.begin();
  std::string expected(
      "Core:O=C(NC1CCN([*:3])CC1)[*:1] R1:Clc1cc([*:1])[nH]c1Cl "
      "R3:O=[N+]([O-])c1cccnc1[*:3]");
  CHECK_RGROUP(it, expected);

  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp2(*core, params);
  res = decomp2.add(*mol);
  TEST_ASSERT(res == 0);

  decomp2.process();
  rows = decomp2.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  it = rows.begin();
  CHECK_RGROUP(it, expected);
}

void testNoAlignmentAndSymmetry() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test NoAlignment with symmetric groups" << std::endl;
  const std::vector<ROMOL_SPTR> cores{"c([*:1])1c([*:2])c([*:3])ccc1"_smiles,
                                      "c([*:3])1c([*:2])c([*:1])cnc1"_smiles};
  const std::vector<const char *> smilesData{"c1(CO)c(F)c(CN)ccc1",
                                             "c1(CO)c(Cl)c(CN)cnc1"};

  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = true;
  params.removeHydrogensPostMatch = true;
  params.alignment = NoAlignment;
  RGroupDecomposition decomp(cores, params);
  size_t i = 0;
  for (const auto &smi : smilesData) {
    ROMOL_SPTR mol(static_cast<ROMol *>(SmilesToMol(smi)));
    TEST_ASSERT(decomp.add(*mol) == static_cast<int>(i++));
  }
  TEST_ASSERT(decomp.process());
  auto rows = decomp.getRGroupsAsRows();
  const std::vector<const char *> res{
      "Core:c1cc([*:1])c([*:2])c([*:3])c1 R1:NC[*:1] R2:F[*:2] R3:OC[*:3]",
      "Core:c1ncc([*:3])c([*:2])c1[*:1] R1:NC[*:1] R2:Cl[*:2] R3:OC[*:3]"};
  TEST_ASSERT(rows.size() == res.size());
  i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
    CHECK_RGROUP(it, res.at(i++));
  }
}

void testSingleAtomBridge() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test single atom bridge between 2 user r groups"
                       << std::endl;

  auto core = "C1([*:1])C([*:2])CC1"_smiles;
  RGroupDecompositionParameters params;
  RGroupDecomposition decomp(*core, params);
  auto mol = "C1CC2NC12"_smiles;
  params.onlyMatchAtRGroups = true;
  auto res = decomp.add(*mol);
  TEST_ASSERT(res == 0);
  TEST_ASSERT(decomp.process());
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  const std::string expected(
      "Core:C1CC([*:2])C1[*:1] R1:N([*:1])[*:2]"
      " R2:N([*:1])[*:2]");
  RGroupRows::const_iterator it = rows.begin();
  CHECK_RGROUP(it, expected);

  core = "C1([*:1])CCC1"_smiles;
  RGroupDecomposition decomp3(*core, params);
  res = decomp3.add(*mol);
  TEST_ASSERT(res == -1);

  params.onlyMatchAtRGroups = false;
  RGroupDecomposition decomp2(*core, params);
  res = decomp2.add(*mol);
  TEST_ASSERT(res == 0);
  TEST_ASSERT(decomp2.process());
  rows = decomp2.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  it = rows.begin();
  CHECK_RGROUP(it, expected);

  params.onlyMatchAtRGroups = true;
  params.allowNonTerminalRGroups = true;
  core = "C1([*:1])[*:2]CC1"_smiles;
  RGroupDecomposition decomp4(*core, params);
  res = decomp4.add(*mol);
  TEST_ASSERT(res == 0);
  TEST_ASSERT(rows.size() == 1)
  it = rows.begin();
  CHECK_RGROUP(it, expected);

  // Now that issue 4505 (Cores with query atoms may fail to R-group-decompose
  // molecules) is resolved this will match with no substitution for R3
  core = "C1([*:1])C([*:2])*([*:3])C1"_smiles;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp5(*core, params);
  mol = "C1OC2NC12"_smiles;
  res = decomp5.add(*mol);
  TEST_ASSERT(res == 0);
  TEST_ASSERT(rows.size() == 1)
  it = rows.begin();
  CHECK_RGROUP(it, expected);
}

void testAddedRGroupsHaveCoords() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test added R groups have non-zero coords"
                       << std::endl;
  auto core = R"CTAB(
     RDKit          2D

 13 14  0  0  0  0  0  0  0  0999 V2000
   -3.9955    2.1866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7099    1.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7099    0.9490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9955    0.5365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2810    0.9490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2810    1.7741    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4244    2.1866    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9956   -0.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2811   -0.7009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1948   -1.5214    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3878   -1.6929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9753   -0.9785    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5274   -0.3654    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  4  8  1  0  0  0  0
  8  9  1  0  0  0  0
 10 11  1  0  0  0  0
  9 10  1  0  0  0  0
  9 13  1  0  0  0  0
 11 12  1  0  0  0  0
 12 13  1  0  0  0  0
  2  7  1  0  0  0  0
M  RGP  1   7   1
M  END
)CTAB"_ctab;
  TEST_ASSERT(core);
  auto mol = "COc1ccc(CC2CCNC2)cc1NC(C)=O"_smiles;
  RGroupDecompositionParameters params;
  RGroupDecomposition decomp(*core, params);
  TEST_ASSERT(decomp.add(*mol) == 0);
  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  for (const auto &row : rows) {
    auto it = row.find("Core");
    TEST_ASSERT(it != row.end());
    const auto &rgdCore = it->second;
    TEST_ASSERT(rgdCore->getNumConformers() == 1);
    size_t r2Num = 0;
    for (const auto atom : rgdCore->atoms()) {
      // test that R2 has non-zero coords and a sensible bond length
      // to its neighboring atom
      if (atom->getAtomicNum() == 0 && atom->getAtomMapNum() == 2) {
        ++r2Num;
        auto &r2Coord = rgdCore->getConformer().getAtomPos(atom->getIdx());
        TEST_ASSERT(fabs(r2Coord.x) > 1e-4);
        TEST_ASSERT(fabs(r2Coord.y) > 1e-4);
        for (const auto &nbri :
             boost::make_iterator_range(rgdCore->getAtomNeighbors(atom))) {
          const auto nbr = (*rgdCore)[nbri];
          const auto bond =
              rgdCore->getBondBetweenAtoms(nbr->getIdx(), atom->getIdx());
          TEST_ASSERT(bond);
          auto &nbrCoord = rgdCore->getConformer().getAtomPos(nbr->getIdx());
          auto bondLen = (nbrCoord - r2Coord).length();
          TEST_ASSERT(fabs(bondLen - 1.0) < 0.1);
        }
        TEST_ASSERT(atom->getDegree() == 1);
      }
    }
    TEST_ASSERT(r2Num == 1);
  }
}

void testUserMatchTypes() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test user rgroup label specification and matching"
                       << std::endl;

  struct TestMatchType {
    static void test(RWMol &core, RWMol &mol,
                     RGroupDecompositionParameters &parameters,
                     std::string &expected) {
      RGroupDecomposition decomp(core, parameters);
      auto res = decomp.add(mol);
      TEST_ASSERT(res == 0);
      TEST_ASSERT(decomp.process());
      auto rows = decomp.getRGroupsAsRows();
      TEST_ASSERT(rows.size() == 1)
      RGroupRows::const_iterator it = rows.begin();
      CHECK_RGROUP(it, expected);
    }
  };

  std::vector<RGroupScore> matchtype{Match, FingerprintVariance};
  for (auto match : matchtype) {
    auto mol = "C1CCCCC1(N)(O)"_smiles;
    auto core = "C1CCCCC1[*:1]"_smiles;
    core = "C1CCCCC1[*:1]"_smarts;
    RGroupDecompositionParameters params;
    params.onlyMatchAtRGroups = true;
    params.scoreMethod = match;
    RGroupDecomposition decomp(*core, params);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == -1);

    params.onlyMatchAtRGroups = false;
    std::string expected("Core:C1CCC([*:1])([*:2])CC1 R1:O[*:1] R2:N[*:2]");
    TestMatchType::test(*core, *mol, params, expected);
    core = "C1CCCCC1([*:1])([*:2])"_smiles;
    TestMatchType::test(*core, *mol, params, expected);
    core = "C1CCCC[*:2]1[*:1]"_smiles;
    params.allowNonTerminalRGroups = true;
    TestMatchType::test(*core, *mol, params, expected);
  }
}

void testUnlabelledRGroupsOnAromaticNitrogen() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test unlabelled R groups on aromatic nitrogens"
                       << std::endl;

  auto core = "c1ccc(-c2cccc3[nH]ncc23)nc1"_smiles;
  auto mol1 = "c1ccc(-c2cccc3n(C)ncc23)nc1"_smiles;
  auto mol2 = "c1ccc(-c2cccc3[nH]ncc23)[n+](CC)c1"_smiles;
  RGroupDecompositionParameters params;
  RGroupDecomposition decomp(*core, params);
  TEST_ASSERT(decomp.add(*mol1) == 0);
  TEST_ASSERT(decomp.add(*mol2) == 1);
  TEST_ASSERT(decomp.process());
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 2);
  size_t i = 0;
  std::vector<std::string> expected{
      "Core:c1ccc(-c2cccc3c2cnn3[*:2])nc1 R2:C[*:2]",
      "Core:c1cc[n+]([*:1])c(-c2cccc3c2cnn3[*:2])c1 R1:CC[*:1] R2:[H][*:2]",
  };
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
    CHECK_RGROUP(it, expected.at(i++));
  }
}

void testAddHsDoesNotFail() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that AddHs does not fail" << std::endl;
  auto core = "[1*]c1ccc([2*])cn1"_smiles;
  auto mol = "Fc1ccc([2*])cn1"_smiles;
  RGroupDecompositionParameters params;
  RGroupDecomposition decomp(*core, params);
  TEST_ASSERT(decomp.add(*mol) == 0);
}

void testNoTempLabels() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that temp labels are removed from results"
                       << std::endl;
  auto core = "[1*]c1ccccc1"_smiles;
  auto mol = "Cc1ccccc1"_smiles;
  RGroupDecompositionParameters params;
  RGroupDecomposition decomp(*core, params);
  TEST_ASSERT(decomp.add(*mol) == 0);
  TEST_ASSERT(decomp.process());
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1);
  const auto &res = rows.front();
  for (const auto &pair : res) {
    for (const auto a : pair.second->atoms()) {
      for (const auto &propName : a->getPropList()) {
        TEST_ASSERT(propName.find("label") == std::string::npos);
      }
    }
  }
}

void testNoSideChains() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test the no side chain return code" << std::endl;
  auto core = "[H]C([H])([H])[H]"_smarts;
  auto mol = "C"_smiles;
  RGroupDecompositionParameters params;
  RGroupDecomposition decomp(*core, params);
  TEST_ASSERT(decomp.add(*mol) == -2);
}

void testDoNotAddUnnecessaryRLabels() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that R labels are not added if not necessary"
                       << std::endl;
  {
    auto m1 = "c1c(F)cccn1"_smiles;
    auto m2 = "c1c(Cl)c(C)ccn1"_smiles;
    auto m3 = "c1c(O)cccn1"_smiles;
    auto m4 = "c1cc(C)c(F)cn1"_smiles;
    auto core = "c1c([*:1])c([*:2])ccn1"_smiles;
    // the test runs twice, with and without symmetrization
    // results should be the same as do not depend on permutations
    for (unsigned int i = 0; i < 2; ++i) {
      RGroupDecompositionParameters ps;
      if (i) {
        ps.matchingStrategy = RDKit::NoSymmetrization;
      }
      RGroupDecomposition decomp(*core, ps);
      TEST_ASSERT(decomp.add(*m1) == 0);
      TEST_ASSERT(decomp.add(*m2) == 1);
      TEST_ASSERT(decomp.add(*m3) == 2);
      TEST_ASSERT(decomp.add(*m4) == 3);
      decomp.process();
      std::stringstream ss;
      auto groups = decomp.getRGroupsAsColumns();
      std::set<std::string> r_labels;
      for (auto &column : groups) {
        r_labels.insert(column.first);
        ss << "Rgroup===" << column.first << std::endl;
        for (auto &rgroup : column.second) {
          ss << MolToSmiles(*rgroup) << std::endl;
        }
      }
      // We only want two groups added

      TEST_ASSERT(r_labels == std::set<std::string>({"Core", "R1", "R2"}));
      TEST_ASSERT(groups.size() == 3);

      TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
c1cc([*:2])c([*:1])cn1
c1cc([*:2])c([*:1])cn1
c1cc([*:2])c([*:1])cn1
c1cc([*:2])c([*:1])cn1
Rgroup===R1
F[*:1]
Cl[*:1]
O[*:1]
F[*:1]
Rgroup===R2
[H][*:2]
C[*:2]
[H][*:2]
C[*:2]
)RES");
    }
  }
  {
    auto cores = std::vector<ROMOL_SPTR>{"[*:1]c1ccccc1"_smiles,
                                         "[*:1]c1c([*:2])ccc(F)c1"_smiles};
    auto m1 = "Cc1c(Cl)ccc(F)c1"_smiles;
    auto m2 = "Cc1c(Br)cccc1"_smiles;
    size_t i;
    RGroupRows rows;
    RGroupDecomposition decomp(cores);
    TEST_ASSERT(decomp.add(*m1) == 0);
    decomp.process();
    rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1);
    std::vector<std::string> expected1{
        "Core:Fc1ccc([*:2])c([*:1])c1 R1:C[*:1] R2:Cl[*:2]"};
    i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
      CHECK_RGROUP(it, expected1.at(i++));
    }
    TEST_ASSERT(decomp.add(*m2) == 1);
    decomp.process();
    rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 2);
    std::vector<std::string> expected2{
        "Core:Fc1ccc([*:2])c([*:1])c1 R1:C[*:1] R2:Cl[*:2]",
        "Core:c1ccc([*:3])c([*:1])c1 R1:C[*:1] R3:Br[*:3]"};
    i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
      CHECK_RGROUP(it, expected2.at(i++));
    }
  }
  {
    auto cores = std::vector<ROMOL_SPTR>{"[*:1]c1c([*:2])ccc(F)c1"_smiles,
                                         "[*:1]c1ccccc1"_smiles,
                                         "[*:1]c1c([*:2])cccc1"_smiles};
    auto mols = std::vector<ROMOL_SPTR>{
        "Cc1c(Cl)ccc(F)c1"_smiles, "Cc1c(Cl)ccc(N)c1"_smiles,
        "Cc1ccccc1"_smiles, "Cc1c(Br)cccc1"_smiles, "Cc1cc(F)ccc1"_smiles};
    RGroupRows rows;
    RGroupDecomposition decomp(cores);
    int n = 0;
    for (const auto &m : mols) {
      TEST_ASSERT(decomp.add(*m) == n++);
    }
    decomp.process();
    rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 5);
    std::vector<std::string> expected{
        "Core:Fc1ccc([*:2])c([*:1])c1 R1:C[*:1] R2:Cl[*:2]",
        "Core:c1cc([*:2])c([*:1])cc1[*:3] R1:C[*:1] R2:Cl[*:2] R3:N[*:3]",
        "Core:c1cc([*:1])cc([*:3])c1 R1:C[*:1] R3:[H][*:3]",
        "Core:c1cc([*:2])c([*:1])cc1[*:3] R1:C[*:1] R2:Br[*:2] R3:[H][*:3]",
        "Core:Fc1ccc([*:2])c([*:1])c1 R1:C[*:1] R2:[H][*:2]"};
    size_t i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
      CHECK_RGROUP(it, expected.at(i++));
    }
  }
}

void testCoreWithAlsRecords() {
  auto core = R"CTAB(
  Mrv2008 11112113312D

  6  6  6  0  0  0            999 V2000
  -13.7277    2.6107    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -14.4421    2.1982    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -14.4421    1.3732    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -13.7277    0.9607    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -13.0132    1.3732    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -13.0132    2.1982    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  1 F    2   6   7
  2 F    2   6   7
  3 F    2   6   7
  4 F    2   6   7
  5 F    2   6   7
  6 F    2   6   7
M  ALS   1  2 F C   N   
M  ALS   2  2 F C   N   
M  ALS   3  2 F C   N   
M  ALS   4  2 F C   N   
M  ALS   5  2 F C   N   
M  ALS   6  2 F C   N   
M  END
)CTAB"_ctab;
  TEST_ASSERT(core);
  std::string sma = MolToSmarts(*core);
  TEST_ASSERT(sma == "[#6,#7]1:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]:[#6,#7]:1");

  auto structure = "ClC1=CN=C(C=C1)N1CCCC1"_smiles;
  RGroupDecomposition decomp(*core);
  TEST_ASSERT(decomp.add(*structure) == 0);
  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  auto core_out = rows[0]["Core"];
  auto core_mol_block = MolToMolBlock(*rows[0]["Core"]);
  auto pos = core_mol_block.find("ALS");
  TEST_ASSERT(pos == std::string::npos);
  std::string expected(
      "Core:c1cc([*:2])ncc1[*:1] R1:Cl[*:1] R2:C1CCN([*:2])C1");
  RGroupRows::const_iterator it = rows.begin();
  CHECK_RGROUP(it, expected);
}

void testAlignOutputCoreToMolecule() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that output core is aligned to input molecule"
                       << std::endl;
  struct Helper {
    static RDGeom::Point3D findPointForAtomNumber(const ROMol &mol,
                                                  int atomNumber) {
      for (const auto atom : mol.atoms()) {
        if (atom->getAtomicNum() == atomNumber) {
          return mol.getConformer().getAtomPos(atom->getIdx());
        }
      }
      throw std::runtime_error("Can't find atom in molecule");
    }
  };

  auto core = R"CTAB(
  Mrv2008 11162112382D

 15 16  0  0  0  0            999 V2000
  -13.4365    2.6419    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.0116    1.9347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.1867    1.9491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.7867    2.6706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.2116    3.3778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.0365    3.3634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.9618    2.6850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5619    3.4066    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7370    3.4210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5369    1.9779    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3121    4.1282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3370    2.6994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.5121    2.6850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0873    3.3922    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4872    4.1138    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  3  4  2  0  0  0  0
  5  6  2  0  0  0  0
  4  5  1  0  0  0  0
  2  3  1  0  0  0  0
  1  6  1  0  0  0  0
  4  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  7 10  2  0  0  0  0
 12 13  1  0  0  0  0
 13 14  2  0  0  0  0
 14 15  1  0  0  0  0
 11 15  2  0  0  0  0
 11  9  1  0  0  0  0
  9 12  2  0  0  0  0
M  END
      )CTAB"_ctab;
  auto mol = R"CTAB(
  -OEChem-03051316302D
  
 19 21  0     0  0  0  0  0  0999 V2000
   -1.7593    5.0073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6211    4.4999    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6183    3.4999    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7450    3.0022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8744    3.5045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8860    4.5096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7423    2.0022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6070    1.4999    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8750    1.5045    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8723    0.5045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0013   -1.0057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8697   -1.5068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7364   -1.0079    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7377   -0.0022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6056    0.5056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4735   -0.0022    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4735   -1.0079    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6056   -1.5057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  4  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 15  2  0  0  0  0
 10 11  1  0  0  0  0
 11 12  2  0  0  0  0
 12 13  1  0  0  0  0
 13 14  2  0  0  0  0
 14 19  1  0  0  0  0
 14 15  1  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 17 18  1  0  0  0  0
 18 19  2  0  0  0  0
M  END
)CTAB"_ctab;

  RGroupRows rows;
  RGroupDecomposition decomp(*core);
  TEST_ASSERT(decomp.add(*mol) == 0);
  decomp.process();
  rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  auto coreOut = rows[0]["Core"];

  for (int atomNumber = 7; atomNumber <= 8; atomNumber++) {
    const auto &coreInPoint = Helper::findPointForAtomNumber(*core, atomNumber);
    const auto &molInPoint = Helper::findPointForAtomNumber(*mol, atomNumber);
    const auto &coreOutPoint =
        Helper::findPointForAtomNumber(*coreOut, atomNumber);
    TEST_ASSERT(fabs(coreInPoint.x - molInPoint.x) > 0.25);
    TEST_ASSERT(fabs(coreOutPoint.x - molInPoint.x) < 1e-10);
    TEST_ASSERT(fabs(coreInPoint.y - molInPoint.y) > 0.25);
    TEST_ASSERT(fabs(coreOutPoint.y - molInPoint.y) < 1e-10);
  }
}

void testWildcardInInput() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that dummy atom in input molecule is handled correctly"
      << std::endl;

  auto core = R"CTAB(
Mrv2008 12012115162D          

  6  6  6  0  0  0            999 V2000
  -21.0938  -16.9652    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -21.8082  -17.3777    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -21.8082  -18.2027    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -21.0938  -18.6152    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -20.3793  -18.2027    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  -20.3793  -17.3777    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  6  0  0  0  0
  2  3  7  0  0  0  0
  3  4  6  0  0  0  0
  4  5  7  0  0  0  0
  5  6  6  0  0  0  0
  1  6  7  0  0  0  0
  1 F    2   6   7
  2 F    2   6   7
  3 F    2   6   7
  4 F    2   6   7
  5 F    2   6   7
  6 F    2   6   7
M  ALS   1  2 F C   N   
M  ALS   2  2 F C   N   
M  ALS   3  2 F C   N   
M  ALS   4  2 F C   N   
M  ALS   5  2 F C   N   
M  ALS   6  2 F C   N   
M  END
)CTAB"_ctab;

  auto structure = "CC1CCN(C1)C1=CC(O*)=C(Cl)C=C1C#N"_smiles;
  RGroupDecomposition decomp(*core);
  TEST_ASSERT(decomp.add(*structure) == 0);
  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  RGroupRows::const_iterator it = rows.begin();
  std::string expected(
      "Core:c1c([*:2])c([*:1])cc([*:4])c1[*:3] R1:Cl[*:1] R2:*O[*:2] "
      "R3:CC1CCN([*:3])C1 R4:N#C[*:4]");
  CHECK_RGROUP(it, expected);

  structure = "CC1CCN(C1)C1=CC([*:2])=C(Cl)C=C1C#N"_smiles;
  RGroupDecomposition decomp2(*core);
  TEST_ASSERT(decomp2.add(*structure) == 0);
  decomp2.process();
  rows = decomp2.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  it = rows.begin();
  expected =
      "Core:c1c([*:2])c([*:1])cc([*:4])c1[*:3] R1:Cl[*:1] R2:*[*:2] "
      "R3:CC1CCN([*:3])C1 R4:N#C[*:4]";
  CHECK_RGROUP(it, expected);
}

void testDoNotChooseUnrelatedCores() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that later cores with more R-groups\n"
                       << "are only chosen if superstructures of earlier\n"
                       << "cores" << std::endl;
  {
    // 1st test, two unrelated cores:
    // 1) 5 terms, 1 R-group
    // 2) 6 terms, 3 R-groups
    // dataset molecule can fit core 1 adding 1 non-user defined R label,
    // and core 2 with no addition of R labels
    ROMOL_SPTR core5Terms1RGroup = "[*:1]C1COCN1"_smiles;
    ROMOL_SPTR core6Terms3RGroups = "[*:1]c1ccc([*:2])c([*:3])c1"_smiles;
    std::vector<ROMOL_SPTR> cores{core5Terms1RGroup, core6Terms3RGroups};
    auto m = "Cc1cc(ccc1F)C1NC(CO1)c1cccs1"_smiles;
    // repeat the test twice, with cores in opposite orders
    // in both cases core ordering should be honored and the first
    // core (either 5-term of 6-term) should be chosen, even when adding
    // R labels is required, as the 2 cores are not structurally related
    for (unsigned int i = 0; i < 2; ++i) {
      std::vector<ROMOL_SPTR> orderedCores{cores[i], cores[1 - i]};
      RGroupDecomposition decomp(orderedCores);
      TEST_ASSERT(decomp.add(*m) == 0);
      TEST_ASSERT(decomp.process());
      auto cols = decomp.getRGroupsAsColumns();
      const auto &core = cols["Core"];
      TEST_ASSERT(core.size() == 1);
      TEST_ASSERT(
          core.front()->getRingInfo()->atomRings().front().size() ==
          orderedCores.front()->getRingInfo()->atomRings().front().size());
    }
  }
  {
    // 2nd test: two related cores:
    // 1) 5 terms, 2 R-groups
    // 2) 5 terms, 3 R-groups
    // dataset molecule 1 has 1 substituent, fits both cores
    // dataset molecule 2 has 2 substituents, fits both cores
    // dataset molecule 3 has 3 substituents, fits core 2 with no need to add
    // R labels
    ROMOL_SPTR core5Terms2RGroups = "[*:1]C1COC([*:2])N1"_smiles;
    ROMOL_SPTR core5Terms3RGroups = "[*:1]C1C([*:2])OC([*:3])N1"_smiles;
    std::vector<ROMOL_SPTR> cores{core5Terms2RGroups, core5Terms3RGroups};
    std::vector<ROMOL_SPTR> mols{"CC1NCCO1"_smiles, "CC1NC(F)CO1"_smiles,
                                 "CC1NC(F)C(Cl)O1"_smiles};
    // repeat the test twice, with cores in opposite orders
    // Molecules (1) and (2) should always pick the first core
    // in the order provided, though both cores could fit
    // Molecule (3) should always pick the more specific core 2
    for (unsigned int i = 0; i < 2; ++i) {
      std::vector<ROMOL_SPTR> orderedCores{cores[i], cores[1 - i]};
      RGroupDecomposition decomp(orderedCores);
      int j = 0;
      for (const auto &m : mols) {
        TEST_ASSERT(decomp.add(*m) == j++);
      }
      TEST_ASSERT(decomp.process());
      auto cols = decomp.getRGroupsAsColumns();
      const auto &core = cols["Core"];
      TEST_ASSERT(core.size() == 3);
      TEST_ASSERT(MolToSmiles(*core.at(0)) ==
                  MolToSmiles(*orderedCores.front()));
      TEST_ASSERT(MolToSmiles(*core.at(1)) ==
                  MolToSmiles(*orderedCores.front()));
      TEST_ASSERT(MolToSmiles(*core.at(2)) == MolToSmiles(*core5Terms3RGroups));
    }
  }
}

void atomDegreePreconditionBug() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that we don't get a bad atom degree precondition violation when "
         "the input structure has dummy atoms"
      << std::endl;

  auto structure = R"CTAB(
     RDKit          2D

 12 12  0  0  0  0  0  0  0  0999 V2000
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000   -2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -5.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  4  8  1  0
  8  9  2  0
  8 10  1  0
 10 11  1  0
  7  2  1  0
  6 12  1  0
M  RGP  1  12   3
M  END

)CTAB"_ctab;

  auto core = "[#6]1:[#7]:[#6]:[#6]:[#6]:[#7]:1"_smarts;
  RGroupDecomposition decomp(*core);
  TEST_ASSERT(decomp.add(*structure) == 0);
  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  RGroupRows::const_iterator it = rows.begin();
  std::string expected(
      "Core:c1c([*:1])nc([*:3])nc1[*:2] R1:COC(=O)[*:1] R2:C[*:2] R3:*[*:3]");
  // Check R3 atom labelling
  auto r3 = rows[0]["R3"];
  TEST_ASSERT(r3->getNumAtoms() == 2)
  TEST_ASSERT(r3->getAtomWithIdx(0)->hasProp(common_properties::dummyLabel));
  TEST_ASSERT(r3->getAtomWithIdx(1)->hasProp(common_properties::dummyLabel));
  TEST_ASSERT(r3->getAtomWithIdx(0)->getProp<std::string>(
                  common_properties::dummyLabel) == "*");
  TEST_ASSERT(r3->getAtomWithIdx(1)->getProp<std::string>(
                  common_properties::dummyLabel) == "R3");
  CHECK_RGROUP(it, expected);
}

void testGithub5222() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that Github5222 is fixed" << std::endl;

  auto core = R"CTAB(
  ChemDraw04112214222D

  6  6  3  0  0  0  0  0  0  0999 V2000
   -0.7145    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -0.4125    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.4125    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.8250    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  2 F    2   6   7
  4 F    2   6   7
  6 F    2   6   7
M  ALS   2  2 F C   N   
M  ALS   4  2 F C   N   
M  ALS   6  2 F C   N   
M  END
)CTAB"_ctab;
  std::vector<std::string> smiArray(10, "COc1ccccc1");
  smiArray.push_back("COc1ccncn1");
  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  RGroupDecomposition decomp(*core, params);
  for (const auto &smiles : smiArray) {
    ROMol *mol = SmilesToMol(smiles);
    int res = decomp.add(*mol);
    TEST_ASSERT(res >= 0);
    delete mol;
  }

  decomp.process();
  std::cerr << "Best mapping" << std::endl;
  RGroupRows rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 11);
  for (const auto &row : rows) {
    TEST_ASSERT(row.size() == 2);
    TEST_ASSERT(row.count("Core") == 1);
    TEST_ASSERT(row.count("R1") == 1);
    auto mol = row.at("R1");
    auto groupSmiles = MolToSmiles(*mol);
    TEST_ASSERT(groupSmiles == "CO[*:1]");
  }
}

void testGithub5569() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that Github5269 is fixed" << std::endl;
  auto core = R"CTAB(
ChemDraw09152209202D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0001    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0001   -0.8250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  4  1  0      
  4  5  1  0      
  5  6  1  0      
  6  1  1  0      
M  END
)CTAB"_ctab;

  auto test = R"CTAB(
     RDKit          2D

 20 21  0  0  0  0  0  0  0  0999 V2000
   -2.6437    1.7625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3582    1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3582    0.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6438    0.1125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9293    0.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9293    1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2148    1.7625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5003    1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2141    1.7625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9286    1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3692    2.3459    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.7975    2.3459    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6438   -0.7125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3582   -1.1250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3582   -1.9500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6438   -2.3625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9293   -1.9500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9293   -1.1250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2271   -2.9459    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0604   -2.9459    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  1  1  0
  6  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
  9 11  1  0
  9 12  1  0
  4 13  1  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 16 17  1  0
 17 18  1  0
 18 13  1  0
 16 19  1  0
 16 20  1  0
M  END
)CTAB"_ctab;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  RGroupDecomposition decomp(*core, params);
  decomp.add(*test);

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  auto r2 = rows[0]["R2"];
  auto match = std::find_if(r2->atoms().begin(), r2->atoms().end(),
                            [](Atom *a) { return a->getAtomicNum() == 0; });
  auto dummy = *match;
  int neighborIndex = *r2->getAtomNeighbors(dummy).first;
  auto conf = r2->getConformer();
  auto p1 = conf.getAtomPos(dummy->getIdx());
  auto p2 = conf.getAtomPos(neighborIndex);
  auto length = (p1 - p2).length();
  TEST_ASSERT(fabs(length - 1.0) < 0.25);
}

void testMultipleGroupsToUnlabelledCoreAtomGithub5573() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that Github5573 is fixed" << std::endl;
  auto core = R"CTAB(
  Mrv2008 09172211422D          

  8  8  0  0  0  0            999 V2000
   -7.2098   -3.1928    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9243   -3.6053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9243   -4.4304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2098   -4.8429    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4954   -4.4304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4954   -3.6053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7809   -4.8429    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0664   -4.4304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  5  7  1  0  0  0  0
  7  8  1  0  0  0  0
M  END
)CTAB"_ctab;

  auto test = R"CTAB(
  Mrv2008 09172211422D          

 10 10  0  0  0  0            999 V2000
   -2.6437    1.7625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3582    1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3582    0.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6438    0.1125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9293    0.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9293    1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1730    1.4791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5085    2.1612    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2148    0.1125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5004    0.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  1  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
  5  9  1  0  0  0  0
  9 10  1  0  0  0  0
M  END
)CTAB"_ctab;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  RGroupDecomposition decomp(*core, params);
  auto result = decomp.add(*test);
  TEST_ASSERT(result == 0);
  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1)
  auto row = rows[0];
  std::string expected("Core:COC1CCC([*:1])([*:2])CN1 R1:C[*:1] R2:C[*:2]");
  RGroupRows::const_iterator it = rows.begin();
  CHECK_RGROUP(it, expected);
}

void testGithub4505() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test GitHub 4505 is fixed" << std::endl;
  {
    // this is the first example from issue 4505- I have changed it so that the
    // molecule is not an exact match to the core (if no sidechains are found
    // the decomposition fails)
    auto core = "[n]1([*:1])cc([*:2])ccc1"_smarts;
    auto mol = "n1cc(OC)ccc1"_smiles;
    RGroupDecompositionParameters params;
    params.removeAllHydrogenRGroups = false;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected("Core:c1cncc([*:2])c1 R2:CO[*:2]");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
  {
    // This is the second example from issue 5505.  The core is not listed in
    // the example code so this is my guess
    auto core = "[#6]([*:1])([*:2])~1~[#6]~[#6]~[#6]~[#6]~[#6]~1"_smarts;
    auto mol = "C=1(C)CCCCC1"_smiles;
    RGroupDecompositionParameters params;
    params.removeAllHydrogenRGroups = false;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected("Core:C1=C([*:1])CCCC1 R1:C[*:1]");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
}

void testMultipleGroupsToUnlabelledCoreAtom() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test unlabeled core atom issues" << std::endl;
  {
    // Sulfonamide example
    auto core = "[#6]-1-[#6]-[#6]-[#6]-[#7]-[#16]-1"_smarts;
    auto mol = "O=S1(=O)CCCCN1"_smiles;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected("Core:C1CCS(=[*:1])(=[*:2])NC1 R1:O=[*:1] R2:O=[*:2]");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
  {
    // Smiles/smarts version of github 5573
    auto core = "[#6]-[#8]-[#6]-1-[#6]-[#6]-[#6]-[#6]-[#7]-1"_smarts;
    auto mol = "COC1CCC(C)(C)CN1"_smiles;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected("Core:COC1CCC([*:1])([*:2])CN1 R1:C[*:1] R2:C[*:2]");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
  {
    // Check that sidechains cluster properly
    auto core = "[#6]-[#8]-[#6]-1-[#6]-[#6]-[#6]-[#6]-[#7]-1"_smarts;
    std::vector<std::string> smilesVec{
        "COC1CCC(C)(C)CN1", "COC1CCC(CC)(C)CN1", "COC1CCC(C)(COC)CN1",
        "COC1CCC(C)(OC)CN1", "COC1CCC(CCC)(C)CN1"};
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.matchingStrategy = Exhaustive;
    params.scoreMethod = FingerprintVariance;
    RGroupDecomposition decomp(*core, params);
    for (auto smiles : smilesVec) {
      auto mol = SmilesToMol(smiles);
      auto result = decomp.add(*mol);
      TEST_ASSERT(result > -1);
      delete mol;
    }
    decomp.process();
    auto rows = decomp.getRGroupsAsRows();
    std::vector<std::string> expected{
        "Core:COC1CCC([*:1])([*:2])CN1 R1:C[*:1] R2:C[*:2]",
        "Core:COC1CCC([*:1])([*:2])CN1 R1:CC[*:1] R2:C[*:2]",
        "Core:COC1CCC([*:1])([*:2])CN1 R1:COC[*:1] R2:C[*:2]",
        "Core:COC1CCC([*:1])([*:2])CN1 R1:CO[*:1] R2:C[*:2]",
        "Core:COC1CCC([*:1])([*:2])CN1 R1:CCC[*:1] R2:C[*:2]"};
    TEST_ASSERT(rows.size() == expected.size());
    int i = 0;
    for (auto row = rows.cbegin(); row != rows.cend(); ++row, ++i) {
      CHECK_RGROUP(row, expected[i]);
    }
  }
  {
    // Check core with terminal wildcard - dummy atom labels allowed
    auto core = "[*]-[#8]-[#6]-1-[#6]-[#6]-[#6]-[#6]-[#7]-1"_smarts;
    auto mol = "COC1CCC(C)(C)CN1"_smiles;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.labels = DummyAtomLabels;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected(
        "Core:C1CC([*:2])([*:3])CNC1O[*:1] R1:C[*:1] R2:C[*:2] R3:C[*:3]");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
    // Check core with terminal wildcard - dummy atom labels not allowed
    params.labels = IsotopeLabels;
    RGroupDecomposition decomp2(*core, params);
    result = decomp2.add(*mol);
    TEST_ASSERT(result == 0)
    decomp2.process();
    rows = decomp2.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    row = rows[0];
    std::string expected2("Core:COC1CCC([*:1])([*:2])CN1 R1:C[*:1] R2:C[*:2]");
    it = rows.begin();
    CHECK_RGROUP(it, expected2);
  }
  {
    // Check core with wildcard in ring
    auto core = "[#6]-[#8]-[#6]-1-[#7]-[#6]-[#6]-[#6]-[*]-1"_smarts;
    std::vector<std::string> smilesVec{"COC1CCC(C)(C)CN1", "COC1NCC(C)(C)CO1",
                                       "COC1NCC(C)(C)CN1"};
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    params.matchingStrategy = Exhaustive;
    params.scoreMethod = FingerprintVariance;
    RGroupDecomposition decomp(*core, params);
    for (auto smiles : smilesVec) {
      auto mol = SmilesToMol(smiles);
      auto result = decomp.add(*mol);
      TEST_ASSERT(result > -1);
      delete mol;
    }
    decomp.process();
    auto rows = decomp.getRGroupsAsRows();
    std::vector<std::string> expected{
        "Core:COC1CCC([*:1])([*:2])CN1 R1:C[*:1] R2:C[*:2]",
        "Core:COC1NCC([*:1])([*:2])CO1 R1:C[*:1] R2:C[*:2]",
        "Core:COC1NCC([*:1])([*:2])CN1 R1:C[*:1] R2:C[*:2]"};
    TEST_ASSERT(rows.size() == expected.size());
    int i = 0;
    for (auto row = rows.cbegin(); row != rows.cend(); ++row, ++i) {
      CHECK_RGROUP(row, expected[i]);
    }
  }
}

void testGithub5613() {
  {
    // Original issue from 5613
    auto core = "[*:1]C(=O)NC1CCN([*:3])C1"_smarts;
    auto mol =
        "Cc1c(c(c([nH]1)C(=O)N[C@H]2CCN(C2)c3cc(nc(n3)Cl)C(=O)O)Cl)Cl"_smiles;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = false;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected(
        "Core:O=C(N[C@H]1CCN([*:3])C1)[*:1] "
        "R1:Cc1[nH]c([*:1])c(Cl)c1Cl "
        "R3:O=C(O)c1cc([*:3])nc(Cl)n1");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
  {
    auto core = "[1*]C(=O)*C1~C~C~*([3*])~*~C~1[2*]"_smarts;
    auto mol =
        "CO[C@H]1CN(c2nc(-c3cnc(OCCN4CCN(C)CC4)cn3)c(C(=O)O)s2)CC[C@H]1NC(=O)c1[nH]c(C)c(Cl)c1Cl"_smiles;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = false;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected(
        "Core:O=C(N[C@@H]1CCN([*:3])C[C@@H]1[*:2])[*:1] "
        "R1:Cc1[nH]c([*:1])c(Cl)c1Cl R2:CO[*:2] "
        "R3:CN1CCN(CCOc2cnc(-c3nc([*:3])sc3C(=O)O)cn2)CC1");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
  {
    auto core = "C(=O)*C1~C~C~*~*~C~1"_smarts;
    auto mol =
        "CO[C@H]1CN(c2nc(-c3cnc(OCCN4CCN(C)CC4)cn3)c(C(=O)O)s2)CC[C@H]1NC(=O)c1[nH]c(C)c(Cl)c1Cl"_smiles;
    RGroupDecompositionParameters params;
    params.allowMultipleRGroupsOnUnlabelled = true;
    RGroupDecomposition decomp(*core, params);
    auto result = decomp.add(*mol);
    TEST_ASSERT(result == 0)
    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();
    TEST_ASSERT(rows.size() == 1)
    auto row = rows[0];
    std::string expected(
        "Core:O=C(N[C@@H]1CCN([*:1])C[C@@H]1[*:2])[*:3] "
        "R1:CN1CCN(CCOc2cnc(-c3nc([*:1])sc3C(=O)O)cn2)CC1 "
        "R2:CO[*:2] R3:Cc1[nH]c([*:3])c(Cl)c1Cl");
    RGroupRows::const_iterator it = rows.begin();
    CHECK_RGROUP(it, expected);
  }
}

void testGitHub5631() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that Github563 (proper placement of core R groups) is fixed"
      << std::endl;
  auto core = R"CTAB(
     RDKit          2D

 12 12  0  0  0  0  0  0  0  0999 V2000
   -3.4154    3.2137    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5903    3.2137    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1778    3.9282    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5903    4.6427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4154    4.6427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8278    3.9282    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6528    3.9282    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8278    2.4993    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -1.3528    3.9282    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   -5.0654    4.6427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6528    5.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8904    4.6427    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  5  6  8  0
  4  5  8  0
  3  4  8  0
  2  3  8  0
  1  2  8  0
  1  6  8  0
  6  7  1  0
  1  8  1  0
  3  9  1  0
  7 10  1  0
 10 11  2  0
 10 12  1  0
M  ISO  3   8   2   9   3  12   1
M  RGP  3   8   2   9   3  12   1
M  END
)CTAB"_ctab;

  auto mol1 = R"CTAB(
  -OEChem-02051811442D

 25 27  0     0  0  0  0  0  0999 V2000
   -0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    2.0104    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7350    2.0001    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.6070    1.5001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4790    2.0002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4790    3.0002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6159    3.5053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7439    3.0052    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3465    3.4976    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.3495    4.4976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4849    5.0002    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.2170    4.9951    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1281    4.5828    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8007    5.3248    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3049    6.1952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3214    5.9901    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.7160    7.1067    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    7.7947    5.2157    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    1.7328   -0.0038    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    1.7313   -1.0038    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5995    0.4950    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  1  6  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  4  7  1  0  0  0  0
  7 12  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  1  0  0  0  0
 11 12  1  0  0  0  0
 10 13  1  0  0  0  0
 13 14  1  0  0  0  0
 14 15  2  0  0  0  0
 14 16  1  0  0  0  0
 16 20  1  0  0  0  0
 16 17  2  0  0  0  0
 17 18  1  0  0  0  0
 18 19  2  0  0  0  0
 19 20  1  0  0  0  0
 19 21  1  0  0  0  0
 18 22  1  0  0  0  0
  3 23  1  0  0  0  0
 23 24  2  0  0  0  0
 23 25  1  0  0  0  0
M  CHG  2  23   1  25  -1
M  END
$$$$
)CTAB"_ctab;

  auto mol2 = R"CTAB(
  -OEChem-02051811442D

 45 49  0     1  0  0  0  0  0999 V2000
   -2.9472   -3.8597    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8765   -2.8622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6440   -2.2211    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2684   -1.2927    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2691   -1.3606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0260   -2.3354    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6270   -0.5940    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9700    0.3454    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6420   -0.7667    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    1.1236   -1.3417    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.0104    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
    0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    0.4975    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    1.2077   -0.4429    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8525    0.6702    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1954    1.6096    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    3.0104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8098    3.5999    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5017    4.5529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4983    4.5516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8122    3.6017    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.0846    5.3617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0793    5.2591    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6762    6.2745    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0895    5.3618    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0889    5.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6768    6.0627    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2755    6.9786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2762    7.0868    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6781    6.2789    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8669    7.7850    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8610    7.6760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4524    8.4824    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4465    8.3734    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0457    9.1803    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0449    9.0708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.4449    8.1542    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8558    7.3461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8566    7.4557    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4389    8.0452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7990   -0.4450    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -4.6141   -2.4639    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  6  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  1  1  0  0  0
 10 16  1  0  0  0  0
 10 12  1  0  0  0  0
 12 13  1  0  0  0  0
 13 14  1  0  0  0  0
 14 15  1  0  0  0  0
 15 16  1  0  0  0  0
 16 17  1  1  0  0  0
 16 18  1  0  0  0  0
 18 19  1  0  0  0  0
 14 20  1  0  0  0  0
 20 24  1  0  0  0  0
 20 21  2  0  0  0  0
 21 22  1  0  0  0  0
 22 23  2  0  0  0  0
 23 24  1  0  0  0  0
 23 25  1  0  0  0  0
 25 26  2  0  0  0  0
 25 27  1  0  0  0  0
 22 28  1  0  0  0  0
 28 33  2  0  0  0  0
 28 29  1  0  0  0  0
 29 30  2  0  0  0  0
 30 31  1  0  0  0  0
 31 32  2  0  0  0  0
 32 33  1  0  0  0  0
 31 34  1  0  0  0  0
 34 35  1  0  0  0  0
 35 36  1  0  0  0  0
 36 37  1  0  0  0  0
 37 42  1  0  0  0  0
 37 38  1  0  0  0  0
 38 39  1  0  0  0  0
 39 40  1  0  0  0  0
 40 41  1  0  0  0  0
 41 42  1  0  0  0  0
 40 43  1  0  0  0  0
  4 44  1  0  0  0  0
  3 45  1  0  0  0  0
M  END
)CTAB"_ctab;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = false;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp(*core, params);
  auto result = decomp.add(*mol1);
  TEST_ASSERT(result == 0);
  result = decomp.add(*mol2);
  TEST_ASSERT(result == 1);
  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 2) {
    auto coreRgd = rows[0]["Core"];
    auto match = std::find_if(
        coreRgd->atoms().begin(), coreRgd->atoms().end(), [](Atom *a) {
          return a->getAtomicNum() == 0 && a->getAtomMapNum() == 2;
        });
    TEST_ASSERT(match != coreRgd->atoms().end());
    auto dummy = *match;
    auto neighbor = *coreRgd->atomNeighbors(dummy).begin();
    auto &conf = coreRgd->getConformer();
    auto &dummyPoint = conf.getAtomPos(dummy->getIdx());
    auto &neighborPoint = conf.getAtomPos(neighbor->getIdx());
    auto length = (dummyPoint - neighborPoint).length();
    TEST_ASSERT(fabs(length - 1.0) < 0.005);
    // R2 dummy should be directly above neighbor
    TEST_ASSERT(fabs(dummyPoint.x - neighborPoint.x) < 0.05);
  }

  {
    auto coreRgd = rows[1]["Core"];
    auto match = std::find_if(
        coreRgd->atoms().begin(), coreRgd->atoms().end(), [](const Atom *a) {
          return a->getAtomicNum() == 0 && a->getAtomMapNum() == 2;
        });
    TEST_ASSERT(match != coreRgd->atoms().end());
    auto dummy = *match;
    auto &conf = coreRgd->getConformer();
    auto &dummyPoint = conf.getAtomPos(dummy->getIdx());
    // R2 dummy should be over input chiral oxygen, which is first oxygen of
    // degree 2 in input mol block
    auto &inputPoint = mol2->getConformer(0).getAtomPos(15);
    TEST_ASSERT(fabs(dummyPoint.x - inputPoint.x) < 0.05);
    TEST_ASSERT(fabs(dummyPoint.y - inputPoint.y) < 0.05);
  }
}

void testRGroupCoordinatesAddedToCore() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that coordinates for R groups are properly added to core when the core has coordinates and the target does not"
      << std::endl;
  auto core = R"CTAB(ACS Document 1996
  ChemDraw05202112262D

 17 17  0  0  0  0  0  0  0  0999 V2000
    1.9511    0.6607    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2362    1.0726    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.5214    1.4844    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8244    0.3577    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2375   -0.3563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8256   -1.0712    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0006   -1.0719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4113   -1.7868    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.3578    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0006    0.3570    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2375   -0.3585    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0625   -0.3592    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.6481    1.7875    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.2387   -1.7853    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -1.2363   -1.7875    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4137    1.0711    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0625   -0.3556    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  2  0
  2  4  1  0
  4 10  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  7  9  1  0
  9 10  2  0
  9 11  1  0
 11 12  3  0
  2 13  1  0
  6 14  1  0
  8 15  1  0
 10 16  1  0
  5 17  1  0
M  END
)CTAB"_ctab;
  auto mol =
      "Brc1cc(Br)c(Oc2ccc(cc2C#N)S(=O)(=O)Nc3ncc(Br)s3)c(c1)c4ccccc4"_smiles;
  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = false;
  params.onlyMatchAtRGroups = false;
  RGroupDecomposition decomp(*core, params);
  auto result = decomp.add(*mol);
  TEST_ASSERT(result == 0);
  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  auto coreRgd = rows[0]["Core"];
  auto numberGroups = 0;
  for (const auto atom : coreRgd->atoms()) {
    if (int rGroupNum = atom->getAtomMapNum(); rGroupNum > 0) {
      auto coreAtoms = core->atoms();
      auto originalAtom = std::find_if(
          coreAtoms.begin(), coreAtoms.end(), [rGroupNum](const auto &a) {
            return static_cast<int>(a->getIsotope()) == rGroupNum;
          });
      TEST_ASSERT(originalAtom != coreAtoms.end());
      const auto &originalPoint =
          core->getConformer(0).getAtomPos((*originalAtom)->getIdx());
      const auto &outputPoint =
          coreRgd->getConformer(0).getAtomPos(atom->getIdx());
      TEST_ASSERT(originalPoint.x == outputPoint.x);
      TEST_ASSERT(originalPoint.y == outputPoint.y);
      TEST_ASSERT(originalPoint.z == outputPoint.z);
      numberGroups++;
    }
  }
  TEST_ASSERT(numberGroups == 2);
}

void testStereoGroupsPreserved() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that stereo group information is copied from input structure to core and R groups"
      << std::endl;
  auto core = R"CTAB(
  ChemDraw02132309392D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 10 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -1.071690 1.031297 0.000000 0
M  V30 2 C -1.071690 0.206260 0.000000 0
M  V30 3 C -0.357230 -0.206259 0.000000 0
M  V30 4 C 0.357230 0.206260 0.000000 0
M  V30 5 C 0.357230 1.031297 0.000000 0
M  V30 6 C -0.357230 1.443816 0.000000 0
M  V30 7 C -0.357230 -1.031297 0.000000 0
M  V30 8 N 0.357230 -1.443816 0.000000 0
M  V30 9 R1 -1.071690 -1.443816 0.000000 0
M  V30 10 Cl 1.071690 1.443816 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 3 7
M  V30 8 1 7 8 CFG=1
M  V30 9 1 7 9
M  V30 10 1 5 10
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;

  auto mol1 = R"CTAB(
  ChemDraw02032311272D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 17 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 0.348623 1.031250 0.000000 0
M  V30 2 C 0.348623 0.206250 0.000000 0
M  V30 3 C 1.063094 -0.206250 0.000000 0
M  V30 4 C 1.777565 0.206250 0.000000 0
M  V30 5 C 1.777565 1.031250 0.000000 0
M  V30 6 C 1.063094 1.443750 0.000000 0
M  V30 7 C 1.063094 -1.031250 0.000000 0
M  V30 8 N 1.777565 -1.443750 0.000000 0
M  V30 9 C 0.348623 -1.443750 0.000000 0
M  V30 10 C -0.365848 -1.031250 0.000000 0
M  V30 11 C -0.452084 -0.210769 0.000000 0
M  V30 12 C -1.259056 -0.039242 0.000000 0
M  V30 13 C -1.671556 -0.753713 0.000000 0
M  V30 14 C -1.119523 -1.366808 0.000000 0
M  V30 15 O -2.492036 -0.839949 0.000000 0
M  V30 16 Cl 2.492036 1.443750 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 3 7
M  V30 8 1 7 8 CFG=1
M  V30 9 1 7 9
M  V30 10 1 9 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 14 1 13 14
M  V30 15 1 14 10
M  V30 16 1 13 15 CFG=3
M  V30 17 1 5 16
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(2 7 13)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;

  auto mol2 = R"CTAB(
  ChemDraw02032311382D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 18 19 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.348623 1.388485 0.000000 0
M  V30 2 C 0.348623 0.563485 0.000000 0
M  V30 3 C 1.063094 0.150985 0.000000 0
M  V30 4 C 1.777565 0.563485 0.000000 0
M  V30 5 C 1.777565 1.388485 0.000000 0
M  V30 6 C 1.063094 1.800985 0.000000 0
M  V30 7 C 1.063094 -0.674015 0.000000 0
M  V30 8 N 1.777565 -1.086515 0.000000 0
M  V30 9 C 0.348623 -1.086515 0.000000 0
M  V30 10 C -0.365848 -0.674015 0.000000 0
M  V30 11 C -0.452084 0.146466 0.000000 0
M  V30 12 C -1.259056 0.317993 0.000000 0
M  V30 13 C -1.671556 -0.396478 0.000000 0
M  V30 14 C -1.119523 -1.009572 0.000000 0
M  V30 15 O -2.492036 -0.482714 0.000000 0
M  V30 16 Cl 2.492036 1.800986 0.000000 0
M  V30 17 O -0.063877 -1.800986 0.000000 0
M  V30 18 C 0.761123 -1.800986 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 3 7
M  V30 8 1 7 8 CFG=1
M  V30 9 1 7 9
M  V30 10 1 9 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 14 1 13 14
M  V30 15 1 14 10
M  V30 16 1 13 15 CFG=3
M  V30 17 1 5 16
M  V30 18 1 9 17 CFG=1
M  V30 19 1 9 18 CFG=3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(2 7 9)
M  V30 MDLV30/STEREL1 ATOMS=(1 13)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  params.onlyMatchAtRGroups = false;
  RGroupDecomposition decomp(*core, params);
  auto result = decomp.add(*mol1);
  TEST_ASSERT(result == 0);
  result = decomp.add(*mol2);
  TEST_ASSERT(result == 1);
  decomp.process();

  RGroupRows rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 2);
  TEST_ASSERT(rows[0].size() == 2);
  TEST_ASSERT(rows[1].size() == 2);

  auto core1 = rows[0]["Core"];
  TEST_ASSERT(core1->getStereoGroups().size() == 1);
  TEST_ASSERT(core1->getStereoGroups()[0].getGroupType() ==
              StereoGroupType::STEREO_ABSOLUTE);
  auto r1 = rows[0]["R1"];
  TEST_ASSERT(r1->getStereoGroups().size() == 1);
  TEST_ASSERT(r1->getStereoGroups()[0].getGroupType() ==
              StereoGroupType::STEREO_ABSOLUTE);

  auto core2 = rows[1]["Core"];
  TEST_ASSERT(core2->getStereoGroups().size() == 1);
  TEST_ASSERT(core2->getStereoGroups()[0].getGroupType() ==
              StereoGroupType::STEREO_ABSOLUTE);
  auto r2 = rows[1]["R1"];
  TEST_ASSERT(r2->getStereoGroups().size() == 2);
  TEST_ASSERT(r2->getStereoGroups()[0].getGroupType() ==
              StereoGroupType::STEREO_ABSOLUTE);
  TEST_ASSERT(r2->getStereoGroups()[1].getGroupType() ==
              StereoGroupType::STEREO_OR);
}

void testEnumeratedCore() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that enumerated cores behave properly"
                       << std::endl;

  auto core = R"CTAB(
  Mrv2008 08242317002D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.4167 11.04 0 0
M  V30 2 C -5.7503 10.27 0 0
M  V30 3 C -5.7503 8.6883 0 0
M  V30 4 C -4.4167 7.9183 0 0
M  V30 5 C -3.083 8.6883 0 0
M  V30 6 C -3.083 10.2283 0 0
M  V30 7 * -5.0835 10.655 0 0
M  V30 8 F -5.0835 12.965 0 0
M  V30 9 * -3.083 9.4583 0 0
M  V30 10 Cl -1.928 11.4589 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(2 1 2) ATTACH=ANY
M  V30 8 1 9 10 ENDPTS=(2 5 6) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;

  auto mol1 = "CC1=CC=C(F)C(Cl)=C1"_smiles;
  auto mol2 = "CCC1=C(F)C=CC(Cl)=C1"_smiles;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  params.onlyMatchAtRGroups = false;
  params.doEnumeration = true;

  const char *expected[] = {"Core:Fc1ccc([*:2])cc1Cl R2:C[*:2]",
                            "Core:Fc1ccc(Cl)cc1[*:1] R1:CC[*:1]"};

  RGroupDecomposition decomp(*core, params);
  const auto add11 = decomp.add(*mol1);
  TEST_ASSERT(add11 == 0);
  const auto add12 = decomp.add(*mol2);
  TEST_ASSERT(add12 == 1);
  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    TEST_ASSERT(i < 2);
    CHECK_RGROUP(it, expected[i]);
  }
}

void testTautomerCore() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Test that cores handled as tautomers behave properly"
                       << std::endl;

  const auto core1 = "Oc1ccccn1"_smiles;
  const auto core2 = "O=C1NC=CC=C1"_smiles;
  const auto mol1 = "Cc1cnc(O)cc1Cl"_smiles;
  const auto mol2 = "CC1=CNC(=O)C=C1F"_smiles;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  params.onlyMatchAtRGroups = false;
  params.doTautomers = true;

  const char *expected1[] = {
      "Core:Oc1cc([*:2])c([*:1])cn1 R1:C[*:1] R2:Cl[*:2]",
      "Core:O=c1cc([*:2])c([*:1])c[nH]1 R1:C[*:1] R2:F[*:2]"};
  const char *expected2[] = {
      "Core:Oc1cc([*:1])c([*:2])cn1 R1:Cl[*:1] R2:C[*:2]",
      "Core:O=c1cc([*:1])c([*:2])c[nH]1 R1:F[*:1] R2:C[*:2]"};

  RGroupDecomposition decomp1(*core1, params);
  const auto add11 = decomp1.add(*mol1);
  TEST_ASSERT(add11 == 0);
  const auto add12 = decomp1.add(*mol2);
  TEST_ASSERT(add12 == 1);
  decomp1.process();
  auto rows = decomp1.getRGroupsAsRows();
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    TEST_ASSERT(i < 2);
    CHECK_RGROUP(it, expected1[i]);
  }

  RGroupDecomposition decomp2(*core2, params);
  const auto add21 = decomp2.add(*mol1);
  TEST_ASSERT(add21 == 0);
  const auto add22 = decomp2.add(*mol2);
  TEST_ASSERT(add22 == 1);
  decomp2.process();
  rows = decomp2.getRGroupsAsRows();
  i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    TEST_ASSERT(i < 2);
    CHECK_RGROUP(it, expected2[i]);
  }

  auto core3 = R"CTAB(
  Mrv2008 08072313382D          

  9  9  0  0  0  0            999 V2000
    5.9823    5.0875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9823    4.2625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2679    3.8500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5534    4.2625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5534    5.0875    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.2679    5.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2679    6.3250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.6968    3.8500    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    5.2679    3.0250    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
  1  6  1  0  0  0  0
  2  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  RGP  2   8   1   9   2
M  END
)CTAB"_ctab;
  auto smiles = MolToSmiles(*core3);
  RGroupDecomposition decomp3(*core3, params);
  const auto add31 = decomp3.add(*mol1);
  TEST_ASSERT(add31 == 0);
  const auto add32 = decomp3.add(*mol2);
  TEST_ASSERT(add32 == 1);
  decomp3.process();
  rows = decomp3.getRGroupsAsRows();
  i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    TEST_ASSERT(i < 2);
    CHECK_RGROUP(it, expected2[i]);
  }
}

void testStereoBondBug() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that stereo bonds adjacent to the core or attachment atoms are handled correctly"
      << std::endl;

  const auto core = R"CTAB(ACS Document 1996
  ChemDraw10242316092D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.7145    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0      
  2  3  1  0      
  3  4  2  0      
  4  5  1  0      
  5  6  2  0      
  6  1  1  0      
M  END
)CTAB"_ctab;
  const auto mol = "C/C=C/C1=CC=CC=C1"_smiles;
  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  params.onlyMatchAtRGroups = false;
  params.doEnumeration = false;
  RGroupDecomposition decomp(*core, params);
  const auto add1 = decomp.add(*mol);
  TEST_ASSERT(add1 == 0);
  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  auto r1 = rows[0]["R1"];
  // Check to see that Stereo bond is present and defined
  bool foundStereo = false;
  for (const auto bond : r1->bonds()) {
    if (bond->getStereo() > Bond::STEREOANY) {
      TEST_ASSERT(!foundStereo);
      foundStereo = true;
      TEST_ASSERT(bond->getStereoAtoms().size() == 2);
    }
  }
  TEST_ASSERT(foundStereo);

  const auto core2 = R"CTAB(
  ChemDraw10242316432D

  7  7  0  0  0  0  0  0  0  0999 V2000
   -1.0717    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0717   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3572    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3572    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0717    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0        0
  2  3  1  0        0
  3  4  2  0        0
  4  5  1  0        0
  5  6  2  0        0
  6  1  1  0        0
  5  7  1  0        0
M  END
)CTAB"_ctab;
  RGroupDecomposition decomp2(*core2, params);
  const auto add2 = decomp2.add(*mol);
  TEST_ASSERT(add2 == 0);
  decomp2.process();
  rows = decomp2.getRGroupsAsRows();
  r1 = rows[0]["R1"];
  // Check to see that Stereo bond is not present
  foundStereo = false;
  for (const auto bond : r1->bonds()) {
    if (bond->getStereo() > Bond::STEREOANY) {
      foundStereo = true;
    }
  }
  TEST_ASSERT(!foundStereo);

  const auto core3 = R"CTAB(
  ChemDraw10252316142D

  7  7  0  0  0  0  0  0  0  0999 V2000
   -0.7145    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0        0
  2  3  1  0        0
  3  4  1  0        0
  4  5  2  0        0
  5  6  1  0        0
  6  1  1  0        0
  6  7  2  0        0
M  END
)CTAB"_ctab;
  const auto mol3 = "C/C=C1N=CCC=C/1"_smiles;
  RGroupDecomposition decomp3(*core3, params);
  const auto add3 = decomp3.add(*mol3);
  TEST_ASSERT(add3 == 0);
  decomp3.process();
  rows = decomp3.getRGroupsAsRows();
  const auto c1 = rows[0]["Core"];
  // Check to see that Stereo bond is not present
  foundStereo = false;
  for (const auto bond : c1->bonds()) {
    if (bond->getStereo() > Bond::STEREOANY) {
      TEST_ASSERT(!foundStereo);
      foundStereo = true;
      TEST_ASSERT(bond->getStereoAtoms().size() == 2);
    }
  }
  TEST_ASSERT(foundStereo);
}

void testNotEnumeratedCore() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test that enumerated setting for non enumerated cores behaves properly"
      << std::endl;

  const auto core = "C1CCCCC1"_smarts;
  const auto mol = "C1CCCCC1C"_smiles;

  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  params.onlyMatchAtRGroups = false;
  params.doEnumeration = true;
  params.doTautomers = false;

  const char *expected = "Core:C1CCC([*:1])CC1 R1:C[*:1]";

  RGroupDecomposition decomp(*core, params);
  const auto add11 = decomp.add(*mol);
  TEST_ASSERT(add11 == 0);
  decomp.process();
  auto rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1);
  RGroupRows::const_iterator it = rows.begin();
  CHECK_RGROUP(it, expected);
}

void testRgroupDecompZipping() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog)
      << "Test we can reconstruct rgroup decomps that break rings" << std::endl;

  const auto core = "N1OCC1"_smiles;
  const auto mol = "C1CC2ONC12"_smiles;
  RGroupDecompositionParameters params;
  params.matchingStrategy = GreedyChunks;
  params.allowMultipleRGroupsOnUnlabelled = true;
  params.onlyMatchAtRGroups = false;
  params.doEnumeration = true;
  params.doTautomers = false;
  RGroupDecomposition decomp(*core, params);
  const auto add11 = decomp.add(*mol);
  TEST_ASSERT(add11 == 0);
  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  TEST_ASSERT(rows.size() == 1);
  RGroupRows::const_iterator it = rows.begin();
  std::vector<ROMOL_SPTR> mols;
  for (auto rgroups = it->begin(); rgroups != it->end(); ++rgroups) {
    mols.push_back(rgroups->second);
  }
  auto res = molzip(mols);
  TEST_ASSERT(MolToSmiles(*res) == "C1CC2ONC12")

  for (RGroupRow &rgroup : rows) {
    res = molzip(rgroup);
    TEST_ASSERT(MolToSmiles(*res) == "C1CC2ONC12")
  }
}

int main() {
  RDLog::InitLogs();
  boost::logging::disable_logs("rdApp.debug");
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group Decomposition \n";

#if 1
  testSymmetryMatching(FingerprintVariance);
  testSymmetryMatching();
  testRGroupOnlyMatching();
  testRingMatching();
  testRingMatching3();
  testMultiCore();
  testGithub1550();
  testRemoveHs();

  testMatchOnlyAtRgroupHs();
  testRingMatching2();
  testGitHubIssue1705();
  testGithub2332();
  testSDFGRoupMultiCoreNoneShouldMatch();
  testRowColumnAlignmentProblem();
  testSymmetryIssues();
  testMultipleCoreRelabellingIssues();

  testGaSymmetryMatching(FingerprintVariance);
  testGaSymmetryMatching(Match);
  testGaBatch();

  testUnprocessedMapping();
  testSingleAtomBridge();
#endif
  testSymmetryPerformance();
  testScorePermutations();
  testMultiCorePreLabelled();
  testCoreWithRGroupAdjQuery();
  testGeminalRGroups();
  testMatchOnAnyAtom();
  testNoAlignmentAndSymmetry();
  testAddedRGroupsHaveCoords();
  testUserMatchTypes();
  testUnlabelledRGroupsOnAromaticNitrogen();
  testAddHsDoesNotFail();
  testNoTempLabels();
  testNoSideChains();
  testDoNotAddUnnecessaryRLabels();
  testCoreWithAlsRecords();
  testAlignOutputCoreToMolecule();
  testWildcardInInput();
  testDoNotChooseUnrelatedCores();
  atomDegreePreconditionBug();
  testGithub5222();
  testGithub5569();
  testGithub4505();
  testMultipleGroupsToUnlabelledCoreAtomGithub5573();
  testMultipleGroupsToUnlabelledCoreAtom();
  testGitHub5631();
  testGithub5613();
  testRGroupCoordinatesAddedToCore();
  testStereoGroupsPreserved();
  testTautomerCore();
  testEnumeratedCore();
  testStereoBondBug();
  testNotEnumeratedCore();
  testRgroupDecompZipping();
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
