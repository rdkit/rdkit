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
    std::cerr << "Expected: '" << expected << "'" << std::endl;
    std::cerr << "Got:      '" << result << "'" << std::endl;
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
                              "Core:c1cc:[*:1]:c1 R1:[nH](:[*:1]):[*:1]",
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
    "Core:*1**[*:1](C[*:2])*1 R1:c1cc([*:1])co1 R2:Cl[*:2]",
    "Core:*1**[*:1](C[*:2])*1 R1:c1cc([*:1])c[nH]1 "
    "R2:I[*:2]",
    "Core:*1**[*:1](C[*:2])*1 R1:c1cc([*:1])cs1 R2:F[*:2]"};

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
    "Core:c1co([*:2])cc1[*:1] R1:ClC[*:1]",
    "Core:c1cn([*:2])cc1[*:1] R1:IC[*:1] R2:[H][*:2]",
    "Core:c1cs([*:2])cc1[*:1] R1:FC[*:1]"};

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

const char *coreSmiRes[] = {"Core:C1CCNC([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCNC([*:1])CC1 R1:Cl[*:1].[H][*:1]",
                            "Core:C1CCNC([*:1])CC1 R1:I[*:1].[H][*:1]",
                            "Core:C1CCNC([*:1])CC1 R1:I[*:1].[H][*:1]",
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
      << "test preferring grouping non hydrogens over hydrogens if possible"
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
    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
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
)RES");
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
    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
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
)RES" || ss.str() == R"RES(Rgroup===Core
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
)RES");
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
  for (int i = 0; i < 2; ++i) {
    ROMol *mol = SmilesToMol(smilesData[i]);
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
  TEST_ASSERT(ss.str() ==
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
  for (int i = 0; i < 2; ++i) {
    ROMol *mol = MolBlockToMol(chains[i]);
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
        cores.push_back(ROMOL_SPTR(sdsup.next()));
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
        cores.push_back(ROMOL_SPTR(sdsup.next()));
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
        "Core:N1C(N([*:2])[*:4])C2C(NC1[*:1])[*:5]C([*:3])[*:6]2 "
        "R2:C(CC[*:2])CC[*:4] R4:C(CC[*:2])CC[*:4] R5:N([*:5])[*:5] "
        "R6:C([*:6])[*:6]",
	"Core:N1C(N([*:2])[*:4])C2C(NC1[*:1])[*:5]C([*:3])[*:6]2 "
	"R2:C[*:2] R4:[H][*:4] R5:S([*:5])[*:5] R6:CC(C)C([*:6])[*:6]",
	"Core:C1C([*:1])NC(N([*:2])[*:4])C2C1[*:5]C([*:3])[*:6]2 "
	"R2:C[*:2] R4:[H][*:4] R5:S([*:5])[*:5] R6:CC(C)C([*:6])[*:6]",
	"Core:C1C([*:1])NC(N([*:2])[*:4])C2C1[*:5]C([*:3])[*:6]2 "
	"R2:[H][*:2] R4:[H][*:4] R5:CN([*:5])[*:5] R6:N([*:6])[*:6]"
    };

    int i = 0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end();
         ++it, ++i) {
      TEST_ASSERT(i < 4);
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
    cores.push_back(ROMOL_SPTR(SmilesToMol(smi)));
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

    int i=0;
    for (RGroupRows::const_iterator it = rows.begin(); it != rows.end(); ++it, ++i) {
      CHECK_RGROUP(it, expected[i]);
    }

    for (const auto row : rows) {
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
    for (const auto rg : R1) {
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
  BOOST_LOG(rdInfoLog) << "Testing R-Group symmetry issues \n";

  auto m1 = "c1c(F)cccn1"_smiles; 
  auto m2 = "c1c(Cl)c(C)ccn1"_smiles;
  auto m3 = "c1c(O)cccn1"_smiles;
  auto m4 = "c1cc(C)c(F)cn1"_smiles;
  auto core = "c1c([*:1])c([*:2])ccn1"_smiles;
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
int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group Decomposition \n";

#if 1
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
#endif
  testRowColumnAlignmentProblem();
  testSymmetryIssues();
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
