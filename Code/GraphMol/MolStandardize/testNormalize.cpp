//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogUtils.h>
#include "Normalize.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolStandardize/MolStandardize.h>

#include <iostream>
#include <fstream>

using namespace RDKit;
using namespace MolStandardize;

void test1() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test1" << std::endl;
  std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7;

  Normalizer normalizer;

  // Test sulfoxide normalization.
  smi1 = "CS(C)=O";
  std::shared_ptr<ROMol> m1(SmilesToMol(smi1));
  ROMOL_SPTR normalized(normalizer.normalize(*m1));
  TEST_ASSERT(MolToSmiles(*normalized) == "C[S+](C)[O-]");

  // Test sulfone
  smi2 = "C[S+2]([O-])([O-])O";
  std::shared_ptr<ROMol> m2(SmilesToMol(smi2));
  ROMOL_SPTR normalized2(normalizer.normalize(*m2));
  TEST_ASSERT(MolToSmiles(*normalized2) == "CS(=O)(=O)O");

  // Test 1,3-separated charges are recombined.
  smi3 = "CC([O-])=[N+](C)C";
  std::shared_ptr<ROMol> m3(SmilesToMol(smi3));
  ROMOL_SPTR normalized3(normalizer.normalize(*m3));
  TEST_ASSERT(MolToSmiles(*normalized3) == "CC(=O)N(C)C");

  // Test 1,3-separated charges are recombined.
  smi4 = "C[n+]1ccccc1[O-]";
  std::shared_ptr<ROMol> m4(SmilesToMol(smi4));
  ROMOL_SPTR normalized4(normalizer.normalize(*m4));
  TEST_ASSERT(MolToSmiles(*normalized4) == "Cn1ccccc1=O");

  // Test a case where 1,3-separated charges should not be recombined.
  smi5 = "CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]";
  std::shared_ptr<ROMol> m5(SmilesToMol(smi5));
  ROMOL_SPTR normalized5(normalizer.normalize(*m5));
  TEST_ASSERT(MolToSmiles(*normalized5) ==
              "CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]");

  // Test 1,5-separated charges are recombined.
  smi6 = R"(C[N+](C)=C\C=C\[O-])";
  std::shared_ptr<ROMol> m6(SmilesToMol(smi6));
  ROMOL_SPTR normalized6(normalizer.normalize(*m6));
  TEST_ASSERT(MolToSmiles(*normalized6) == "CN(C)C=CC=O");

  // Test a case where 1,5-separated charges should not be recombined.
  smi7 = "C[N+]1=C2C=[N+]([O-])C=CN2CCC1";
  std::shared_ptr<ROMol> m7(SmilesToMol(smi7));
  ROMOL_SPTR normalized7(normalizer.normalize(*m7));
  TEST_ASSERT(MolToSmiles(*normalized7) == "C[N+]1=C2C=[N+]([O-])C=CN2CCC1");

  // Failed on 1k normalize test sanitizeMol step
  std::string smi8 = "O=c1cc([O-])[n+](C2OC(CO)C(O)C2O)c2sccn12";
  std::shared_ptr<ROMol> m8(SmilesToMol(smi8));
  ROMOL_SPTR normalized8(normalizer.normalize(*m8));
  TEST_ASSERT(MolToSmiles(*normalized8) ==
              "O=c1cc([O-])[n+](C2OC(CO)C(O)C2O)c2sccn12");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  test1();
  return 0;
}
