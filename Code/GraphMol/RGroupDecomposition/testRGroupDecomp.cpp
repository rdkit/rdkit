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
#include <RDBoost/test.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/Exceptions.h>

using namespace RDKit;

void CHECK_RGROUP(RGroupRows::const_iterator &it, std::string expected,
                  bool doassert = true) {
  std::ostringstream str;
  int i = 0;

  for (auto rgroups = it->begin(); rgroups != it->end(); ++rgroups, ++i) {
    if (i) str << " ";
    // rlabel:smiles
    str << rgroups->first << ":" << MolToSmiles(*rgroups->second.get(), true);
  }
  std::string result = str.str();

  if (expected != result) {
    std::cerr << "Expected: " << expected << std::endl;
    std::cerr << "Got:      " << result << std::endl;
  }

  if (doassert) TEST_ASSERT(result == expected);
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

const char *symdata[5] = {"c1(Cl)ccccc1", "c1c(Cl)cccc1", "c1c(Cl)cccc1",
                          "c1cc(Cl)ccc1", "c1ccc(Cl)cc1"};

void testSymmetryMatching() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp symmetry matching" << std::endl;

  RWMol *core = SmilesToMol("c1ccccc1");
  RGroupDecomposition decomp(*core);
  for (int i = 0; i < 5; ++i) {
    ROMol *mol = SmilesToMol(symdata[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    delete mol;
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();

  std::ostringstream str;

  // All Cl's should be labeled with the same rgroup
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it) {
    CHECK_RGROUP(it, "Core:c1ccc([*:1])cc1 R1:Cl[*:1]");
  }
  delete core;
}

const char *matchRGroupOnlyData[5] = {
    "c1(Cl)ccccc1", "c1c(Cl)cccc1",    "c1cc(Cl)ccc1",
    "c1ccc(Cl)cc1", "c1c(Cl)cccc(I)1",
};

void testRGroupOnlyMatching() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp rgroup only matching"
                       << std::endl;

  RWMol *core = SmilesToMol("c1ccccc1[1*]");
  RGroupDecompositionParameters params(IsotopeLabels);
  params.onlyMatchAtRGroups = true;

  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 5; ++i) {
    ROMol *mol = SmilesToMol(matchRGroupOnlyData[i]);
    int res = decomp.add(*mol);
    if (i < 4) {
      TEST_ASSERT(res == i);
    } else {
      TEST_ASSERT(res == -1);
    }
    delete mol;
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  std::ostringstream str;

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, "Core:c1ccc([*:1])cc1 R1:Cl[*:1]");
  }
  delete core;
}

const char *ringData[3] = {"c1cocc1", "c1c[nH]cc1", "c1cscc1"};

const char *ringDataRes[3] = {"Core:c1cc:[*:1]:c1 R1:o(:[*:1]):[*:1]",
                              "Core:c1cc:[*:1]:c1 R1:[H]n(:[*:1]):[*:1]",
                              "Core:c1cc:[*:1]:c1 R1:s(:[*:1]):[*:1]"};

void testRingMatching() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp ring matching" << std::endl;

  RWMol *core = SmilesToMol("c1ccc[1*]1");
  RGroupDecompositionParameters params(IsotopeLabels);

  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 3; ++i) {
    ROMol *mol = SmilesToMol(ringData[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    delete mol;
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  std::ostringstream str;

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, ringDataRes[i]);
  }
  delete core;
}

const char *ringData2[3] = {"c1cocc1CCl", "c1c[nH]cc1CI", "c1cscc1CF"};

const char *ringDataRes2[3] = {
    "Core:*1**[*:1](C[*:2])*1 R1:[H]c1oc([H])c([*:1])c1[H] R2:Cl[*:2]",
    "Core:*1**[*:1](C[*:2])*1 R1:[H]c1c([*:1])c([H])n([H])c1[H] "
    "R2:I[*:2]",
    "Core:*1**[*:1](C[*:2])*1 R1:[H]c1sc([H])c([*:1])c1[H] R2:F[*:2]"};

void testRingMatching2() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp full ring dummy core"
                       << std::endl;

  RWMol *core = SmartsToMol("*1***[*:1]1C[*:2]");
  RGroupDecompositionParameters params;

  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 3; ++i) {
    ROMol *mol = SmilesToMol(ringData2[i]);
    int res = decomp.add(*mol);
    TEST_ASSERT(res == i);
    delete mol;
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  std::ostringstream str;

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, ringDataRes2[i]);
  }
  delete core;
}

const char *ringData3[3] = {"c1cocc1CCl", "c1c[nH]cc1CI", "c1cscc1CF"};

const char *ringDataRes3[3] = {
    "Core:c1co([*:2])cc1[*:1] R1:[H]C([H])(Cl)[*:1]",
    "Core:c1cn([*:2])cc1[*:1] R1:[H]C([H])(I)[*:1] R2:[H][*:2]",
    "Core:c1cs([*:2])cc1[*:1] R1:[H]C([H])(F)[*:1]"};

void testRingMatching3() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test rgroup decomp full ring dummy core"
                       << std::endl;

  RWMol *core = SmartsToMol("*1***[*:1]1");
  RGroupDecompositionParameters params;

  RGroupDecomposition decomp(*core, params);
  for (int i = 0; i < 3; ++i) {
    ROMol *mol = SmilesToMol(ringData3[i]);
    int res = decomp.add(*mol);
    delete mol;
    TEST_ASSERT(res == i);
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  std::ostringstream str;

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
    CHECK_RGROUP(it, ringDataRes3[i]);
  }
  delete core;
}

const char *coreSmi[] = {
    "C1CCNC(Cl)CC1", "C1CC(Cl)NCCC1", "C1CCNC(I)CC1", "C1CC(I)NCCC1",

    "C1CCSC(Cl)CC1", "C1CC(Cl)SCCC1", "C1CCSC(I)CC1", "C1CC(I)SCCC1",

    "C1CCOC(Cl)CC1", "C1CC(Cl)OCCC1", "C1CCOC(I)CC1", "C1CC(I)OCCC1"};

const char *coreSmiRes[] = {"Core:C1CCC([*:2])N([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCC([*:2])N([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCC([*:2])N([*:1])CC1 R1:I[*:1].[H][*:1]",
                            "Core:C1CCC([*:2])N([*:1])CC1 R1:I[*:1].[H][*:1]",
                            "Core:C1CCSC([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCSC([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCSC([*:1])CC1 R1:I[*:1].[H][*:1]",
                            "Core:C1CCSC([*:1])CC1 R1:I[*:1].[H][*:1]",
                            "Core:C1CCOC([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCOC([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCOC([*:1])CC1 R1:I[*:1].[H][*:1]",
                            "Core:C1CCOC([*:1])CC1 R1:I[*:1].[H][*:1]"};

void testMultiCore() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test multi core" << std::endl;
  std::vector<ROMOL_SPTR> cores;
  cores.push_back(ROMOL_SPTR(SmartsToMol("C1CCNCCC1")));
  cores.push_back(ROMOL_SPTR(SmilesToMol("C1CCOCCC1")));
  cores.push_back(ROMOL_SPTR(SmilesToMol("C1CCSCCC1")));

  RGroupDecomposition decomp(cores);
  for (unsigned int i = 0; i < sizeof(coreSmi) / sizeof(const char *); ++i) {
    ROMol *mol = SmilesToMol(coreSmi[i]);
    unsigned int res = decomp.add(*mol);
    delete mol;
    TEST_ASSERT(res == i);
  }

  decomp.process();
  RGroupRows rows = decomp.getRGroupsAsRows();
  std::ostringstream str;

  // All Cl's should be labeled with the same rgroup
  int i = 0;
  for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
       ++it, ++i) {
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
  TEST_ASSERT(rg2->getNumAtoms() == 12);
  MolOps::Kekulize(*rg2);
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
    TEST_ASSERT(rg2->getNumAtoms() == 12);
  }
  {
    RGroupDecompositionParameters params;
    params.removeHydrogensPostMatch = true;
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
}

void testGitHubIssue1705() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test preferring grouping non hydrogens over hydrogens if possible" << std::endl;

  RWMol *core = SmilesToMol("Oc1ccccc1");
  RGroupDecompositionParameters params;

  RGroupDecomposition decomp(*core, params);
  const char *smilesData[5] = {"Oc1ccccc1","Oc1c(F)cccc1","Oc1ccccc1F","Oc1c(F)cc(N)cc1","Oc1ccccc1Cl"};
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
    for (auto &rgroup : column.second ) {
      ss << MolToSmiles(*rgroup) << std::endl;
    }
  }

  TEST_ASSERT(ss.str() == "Rgroup===Core\nOc1ccc([*:2])cc1[*:1]\nOc1ccc([*:2])cc1[*:1]\nOc1ccc([*:2])cc1[*:1]\nOc1ccc([*:2])cc1[*:1]\nOc1ccc([*:2])cc1[*:1]\nRgroup===R1\n[H][*:1]\nF[*:1]\nF[*:1]\nF[*:1]\nCl[*:1]\nRgroup===R2\n[H][*:2]\n[H][*:2]\n[H][*:2]\n[H]N([H])[*:2]\n[H][*:2]\n");


}

void testMatchOnlyAtRgroupHs() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "test matching only rgroups but allows Hs" << std::endl;

  RWMol *core = SmilesToMol("*OCC");
  RGroupDecompositionParameters params;
  params.onlyMatchAtRGroups = true;
  RGroupDecomposition decomp(*core, params);
  const char *smilesData[2] = {"OCC","COCC"};
  for (int i = 0; i < 2; ++i) {
    ROMol *mol = SmilesToMol(smilesData[i]);
    decomp.add(*mol);
  }
  decomp.process();

  std::stringstream ss;
  RGroupColumns groups = decomp.getRGroupsAsColumns();
  for (auto &column : groups) {
    ss << "Rgroup===" << column.first << std::endl;
    for (auto &rgroup : column.second ) {
      ss << MolToSmiles(*rgroup) << std::endl;
    }
  }
  std::cerr << ss.str() << std::endl;
  TEST_ASSERT(ss.str() == "Rgroup===Core\nCCO[*:1]\nCCO[*:1]\nRgroup===R1\n[H][*:1]\n[H]C([H])([H])[*:1]\n");
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group Decomposition \n";

#if 1
  testSymmetryMatching();
  testRGroupOnlyMatching();
  testRingMatching();
  testRingMatching2();
  testRingMatching3();
  testMultiCore();
#endif
  testGithub1550();
  testRemoveHs();

  testGitHubIssue1705();
  testMatchOnlyAtRgroupHs();
  
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
