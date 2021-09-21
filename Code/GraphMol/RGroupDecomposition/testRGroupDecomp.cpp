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
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/RGroupDecomposition/RGroupDecompData.h>
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
    str << rgroups->first << ":" << MolToSmiles(*rgroups->second.get(), true);
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
    str << rgroups.first << "\t" << MolToSmiles(*rgroups.second.get(), true)
        << "\t";
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
  bool isParallelGaEnabled = (sstrm.str().find("This RDKit build does not enable GA parallel runs") == std::string::npos);
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
Oc1ccc([*:2])cc1[*:1]
Oc1ccc([*:2])cc1[*:1]
Oc1ccc([*:2])cc1[*:1]
Oc1ccc([*:2])cc1[*:1]
Oc1ccc([*:2])cc1[*:1]
Rgroup===R1
[H][*:1]
F[*:1]
F[*:1]
F[*:1]
Cl[*:1]
Rgroup===R2
[H][*:2]
[H][*:2]
[H][*:2]
N[*:2]
[H][*:2]
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
F[*:1]
F[*:1]
F[*:1]
Rgroup===R2
[H][*:2]
[H][*:2]
[H][*:2]
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

      int idx = 0;
      while (!sdsup.atEnd()) {
        ROMol *mol = sdsup.next();
        TEST_ASSERT(mol);
        int addedIndex = decomp.add(*mol);
        TEST_ASSERT(addedIndex == -1);  // none should match
        ++idx;
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

      int idx = 0;
      while (!sdsup.atEnd()) {
        ROMol *mol = sdsup.next();
        TEST_ASSERT(mol);
        decomp.add(*mol);
        ++idx;
        delete mol;
      }
    }

    decomp.process();
    RGroupRows rows = decomp.getRGroupsAsRows();

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
  {  // repeat that without symmetrization (testing #3224)
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
    ps.timeout = 1.0;
#else
    ps.timeout = 5.0;
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
F[*:1]
F[*:1]
F[*:1]
Rgroup===R2
[H][*:2]
[H][*:2]
[H][*:2]
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
      if (molecules.size() == 30) break;
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

  core = "C1([*:1])C([*:2])*([*:3])C1"_smiles;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp5(*core, params);
  mol = "C1OC2NC12"_smiles;
  res = decomp5.add(*mol);
  TEST_ASSERT(res == -1);
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
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
