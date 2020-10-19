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
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

using namespace RDKit;

void CHECK_RGROUP(RGroupRows::const_iterator &it, const std::string &expected,
                  bool doassert = true) {
  std::ostringstream str;
  int i = 0;

  for (auto rgroups = it->begin(); rgroups != it->end(); ++rgroups, ++i) {
    if (i) {
      str << " ";
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
    TEST_ASSERT(result == expected);
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
  RGroupDecompositionParameters params;
  params.labels = IsotopeLabels;
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
  RGroupDecompositionParameters params;
  params.labels = IsotopeLabels;

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
  cores.emplace_back(SmartsToMol("C1CCNCCC1"));
  cores.emplace_back(SmilesToMol("C1CCOCCC1"));
  cores.emplace_back(SmilesToMol("C1CCSCCC1"));

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
        "Core:N1C(N([*:2])[*:4])C2C(NC1[*:1])[*:5]C([*:3])[*:6]2 "
        "R2:C(CC[*:2])CC[*:4] R4:C(CC[*:2])CC[*:4] R5:N([*:5])[*:5] "
        "R6:C([*:6])[*:6]",
        "Core:N1C(N([*:2])[*:4])C2C(NC1[*:1])[*:5]C([*:3])[*:6]2 "
        "R2:C[*:2] R4:[H][*:4] R5:S([*:5])[*:5] R6:CC(C)C([*:6])[*:6]",
        "Core:C1C([*:1])NC(N([*:2])[*:4])C2C1[*:5]C([*:3])[*:6]2 "
        "R2:C[*:2] R4:[H][*:4] R5:S([*:5])[*:5] R6:CC(C)C([*:6])[*:6]",
        "Core:C1C([*:1])NC(N([*:2])[*:4])C2C1[*:5]C([*:3])[*:6]2 "
        "R2:[H][*:2] R4:[H][*:4] R5:CN([*:5])[*:5] R6:N([*:6])[*:6]"};

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
      CHECK_RGROUP(it, expected[i]);
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
  BOOST_LOG(rdInfoLog) << "Testing R-Group symmetry issues \n";

  auto m1 = "c1c(F)cccn1"_smiles;
  auto m2 = "c1c(Cl)c(C)ccn1"_smiles;
  auto m3 = "c1c(O)cccn1"_smiles;
  auto m4 = "c1cc(C)c(F)cn1"_smiles;
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
  {  // repeat that without symmetrization (testing #3224)
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
[H][*:1]
Rgroup===R2
[H][*:2]
C[*:2]
[H][*:2]
C[*:2]
Rgroup===R3
[H][*:3]
[H][*:3]
[H][*:3]
F[*:3]
)RES");
  }
}

void testSymmetryPerformance() {
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing R-Group symmetry issues \n";
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
    ps.timeout = 1.0;
    RGroupDecomposition decomp(*core, ps);
    bool ok = false;
    try {
      size_t ndone = 0;
      for (auto m : ms) {
        decomp.add(*m);
        ++ndone;
      }
    } catch (const std::runtime_error &) {
      ok = true;
    }
    TEST_ASSERT(ok);
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
    ps.timeout = 1.0;
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
  BOOST_LOG(rdInfoLog) << "Testing permutation scoring function \n";

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
    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
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
    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
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
)RES");
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
    TEST_ASSERT(ss.str() == R"RES(Rgroup===Core
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
O[*:1]
O[*:1]
C[*:1]
C[*:1]
CC(C)[*:1]
c1ccc(C[*:1])cc1
c1ccc(C[*:1])cc1
CCC[*:1]
CC(C)C[*:1]
Rgroup===R2
[H][*:2]
[H][*:2]
[H][*:2]
C(N[*:1])[*:2]
P[*:2]
N[*:2]
N[*:2]
CC[*:2]
CC[*:2]
CC[*:2]
OC[*:2]
OC[*:2]
CCCC[*:2]
)RES");
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
        CHECK_RGROUP(it, expectedRows[i] /*, false*/);
      }
      RGroupColumns groups = decomp.getRGroupsAsColumns();
      i = 0;
      TEST_ASSERT(groups.size() == 3);
      for (const auto &pair : groups) {
        /*
        if (pair.first != expectedLabels[i]) {
          std::cerr << "ERROR: Expected " << expectedLabels[i] << ", got "
                    << pair.first << std::endl;
        }
        */
        TEST_ASSERT(pair.first == expectedLabels[i]);
        unsigned int j = 0;
        for (const auto &item : pair.second) {
          /*
          if (expectedItems[i][j] != MolToSmiles(*item)) {
            std::cerr << "ERROR: Expected " << expectedItems[i][j] << ", got "
                      << MolToSmiles(*item) << std::endl;
          }
          */
          TEST_ASSERT(expectedItems[i][j] == MolToSmiles(*item));
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
  std::vector<std::string> expectedRowsAutodetect{
      "Core:O=C(c1cncn1[*:2])[*:1] R1:CN[*:1] R2:CC[*:2]",
      "Core:*1:*c2c(*c([*:2])c[*:1]2)nc1[*:3] R1:c(:[*:1]):[*:1] R2:Br[*:2]"};
  std::vector<std::vector<std::string>> expectedItemsAutodetect{
      {"O=C(c1cncn1[*:2])[*:1]", "*1:*c2c(*c([*:2])c[*:1]2)nc1[*:3]"},
      {"CN[*:1]", "c(:[*:1]):[*:1]"},
      {"CC[*:2]", "Br[*:2]"}};
  std::vector<std::string> expectedRowsNoAutodetect{
      "Core:O=C(c1cncn1[*:2])[*:1] R1:CN[*:1] R2:CC[*:2]",
      "Core:*1:*c2*cc([*:2])*c2nc1[*:1] R1:F[*:1] R2:Br[*:2]"};
  std::vector<std::vector<std::string>> expectedItemsNoAutodetect{
      {"O=C(c1cncn1[*:2])[*:1]", "*1:*c2*cc([*:2])*c2nc1[*:1]"},
      {"CN[*:1]", "F[*:1]"},
      {"CC[*:2]", "Br[*:2]"}};
  std::vector<std::string> expectedLabels{"Core", "R1", "R2"};
  RGroupDecompositionParameters params;

  // test pre-labelled with MDL R-group labels, autodetect
  params.labels = AutoDetect;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsAutodetect,
                     expectedItemsAutodetect);
  // test pre-labelled with MDL R-group labels, no autodetect
  params.labels = MDLRGroupLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsNoAutodetect,
                     expectedItemsNoAutodetect);
  // test pre-labelled with MDL R-group labels, autodetect, no MCS alignment
  params.labels = AutoDetect;
  params.alignment = NoAlignment;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsNoAutodetect,
                     expectedItemsNoAutodetect);

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
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsAutodetect,
                     expectedItemsAutodetect);
  // test pre-labelled with isotopic labels, no autodetect
  params.labels = IsotopeLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsNoAutodetect,
                     expectedItemsNoAutodetect);
  // test pre-labelled with isotopic labels, autodetect, no MCS alignment
  params.labels = AutoDetect;
  params.alignment = NoAlignment;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsNoAutodetect,
                     expectedItemsNoAutodetect);

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
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsAutodetect,
                     expectedItemsAutodetect);
  // test pre-labelled with atom map labels, no autodetect
  params.labels = AtomMapLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsNoAutodetect,
                     expectedItemsNoAutodetect);
  // test pre-labelled with atom map labels, autodetect, no MCS alignment
  params.labels = AutoDetect;
  params.alignment = NoAlignment;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsNoAutodetect,
                     expectedItemsNoAutodetect);

  for (auto &core : cores) {
    for (auto a : core->atoms()) {
      if (a->getAtomMapNum()) {
        a->setAtomMapNum(0);
      }
    }
  }
  // test pre-labelled with dummy atom labels, autodetect
  expectedRowsAutodetect = std::vector<std::string>{
      "Core:O=C(c1cncn1[*:2])[*:1] R1:CN[*:1] R2:CC[*:2]",
      "Core:c1c([*:2])[*:3]c2nc([*:6])[*:5]:[*:4]c2[*:1]1 R1:c(:[*:1]):[*:1] "
      "R2:Br[*:2]"};
  expectedItemsAutodetect = std::vector<std::vector<std::string>>{
      {"O=C(c1cncn1[*:2])[*:1]",
       "c1c([*:2])[*:3]c2nc([*:6])[*:5]:[*:4]c2[*:1]1"},
      {"CN[*:1]", "c(:[*:1]):[*:1]"},
      {"CC[*:2]", "Br[*:2]"}};
  params.labels = AutoDetect;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsAutodetect,
                     expectedItemsAutodetect);
  // test pre-labelled with dummy atom labels, no autodetect
  // in this case there is no difference from autodetect as the RGD code
  // cannot tell the difference between query atoms and dummy R-groups
  params.labels = DummyAtomLabels | RelabelDuplicateLabels;
  params.alignment = MCS;
  MultiCoreRGD::test(cores, params, expectedLabels, expectedRowsAutodetect,
                     expectedItemsAutodetect);
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
  testRowColumnAlignmentProblem();
  testSymmetryIssues();
#endif
  testSymmetryPerformance();
  testScorePermutations();
  testMultiCorePreLabelled();
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  return 0;
}
