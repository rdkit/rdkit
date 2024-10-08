//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <future>
#include <thread>
#endif

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit::v2;

TEST_CASE("v2 basics") {
  {
    auto mol = SmilesParse::MolFromSmiles("CCO[H]");
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
  }
  {
    SmilesParse::SmilesParserParams ps;
    ps.removeHs = false;
    auto mol = SmilesParse::MolFromSmiles("CCO[H]", ps);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  {
    auto mol = SmilesParse::MolFromSmarts("[H]CC[R]");
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 4);
  }
  {
    SmilesParse::SmartsParserParams ps;
    ps.mergeHs = true;
    auto mol = SmilesParse::MolFromSmarts("[H]CC[R]", ps);
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 3);
  }
  {
    auto atm = SmilesParse::AtomFromSmiles("C");
    REQUIRE(atm);
  }
  {
    auto bnd = SmilesParse::BondFromSmiles("-");
    REQUIRE(bnd);
  }
  {
    auto atm = SmilesParse::AtomFromSmarts("[R]");
    REQUIRE(atm);
  }
  {
    auto bnd = SmilesParse::BondFromSmarts("@");
    REQUIRE(bnd);
  }
}
TEST_CASE("handling of aromatic Al in SMILES") {
  SECTION("basics") {
    auto mol = SmilesParse::MolFromSmiles("[Al+]1cccccccccc1");
    REQUIRE(mol);
    auto smi = RDKit::MolToSmiles(*mol);
    CHECK(smi.find("Al") != std::string::npos);
  }
}

#ifdef RDKit_BUILD_THREADSAFE_SSS
TEST_CASE("multithreaded") {
  // SMILES from Regress/Data/znp.50k.smi
  std::vector<std::string> smileses = {
      "CC(=O)OC(C)(C)CCC1OC(C)(C)OC1(C)C1CCC2(O)C3=CC(=O)C4CC5OC(C)(C)OC5CC4(C)C3CCC12C ZINC70701530",
      "O=C([O-])CC[NH2+]CCCC[NH2+]CCC(=O)[O-] ZINC01628630",
      "CCOC(=O)C1CCC[NH+](Cc2c([O-])ccc3c2O/C(=C\\c2cc(Br)cc4c2OCOC4)C3=O)C1 ZINC71389266",
      "COc1ccc(C(C)OCS(=O)(=O)c2ccccc2)cc1 ZINC03847541",
      "CC(=O)OC1CC(O)C2(C)C(CC(OC(C)=O)C3(C)C2CCC2(C)C(c4ccoc4)C(=O)C4OC423)C1(C)C ZINC67911113",
      "CCCC[NH+](CCCC)Cc1c([O-])cc(C)c2c1O/C(=C\\C=C\\c1ccccc1)C2=O ZINC41584184",
      "CC([NH2+]CCC(=O)[O-])C(O)c1ccccc1 ZINC02124320",
      "O=C(OCc1ccccc1)C(F)Cl ZINC02386460",
      "CC(=O)NCC(C)CCC(=O)OC1CC2C3CC=C4CC(OC(C)=O)CCC4(C)C3CCC2(C)C1NC(C)=O ZINC70698628",
      "CCCC[NH+](C)Cc1c([O-])cc(C)c2c1O/C(=C\\c1cc(Cl)cc3c1OCOC3)C2=O ZINC41584052",
      "Cc1cc2c(c3c1C(=O)/C(=C/c1ccc(Br)cc1)O3)CN(CCCc1ccccc1)CO2 ZINC41648570",
      "CSc1nc(N(C)c2ccccc2)n2nccc2n1 ZINC71404846",
      "C=C1CCCC2(C)CC3OC(=O)C(C[NH+](Cc4ccc(OC)c(OC)c4)CC(C)O)C3CC12 ZINC08790343",
      "O=C(NCCc1ccccc1F)c1cc2c3ccccc3[nH]c2c(-c2ccc3c(c2)OCO3)n1 ZINC08790829",
      "Cc1cc([O-])c(C[NH+]2CCN(c3ccc(F)cc3)CC2)c2c1C(=O)/C(=C/c1ccc(Cl)cc1Cl)O2 ZINC41585528",
      "C=C(C)CCCC1(C)CCCC(C)O1 ZINC02132870",
      "OCC1OC(C[NH2+]C2CCOCC2)C(N2CCCC2)C1O ZINC20465048",
      "CC(=O)Oc1ccc2c(c1)CCC1C2CCC2(C)C1CCC21OCCO1 ZINC04026970",
      "CC(=O)c1cccc(OCC(=O)Nc2ccc(C)cc2)c1 ZINC03499981",
      "CC(C)(C)C(=O)COc1ccc2c(-c3cc4ccccc4oc3=O)cc(=O)oc2c1 ZINC02138521",
      "CCOC(=O)C(Cc1ccc(O)cc1)OCC ZINC14455512",
      "CCCN1CCCN(C)CCN(C)CCC[NH+](CCC)CC1 ZINC68604820",
      "CCOC(=O)N1CCN(Cc2c(O)ccc3c(=O)c(Oc4ccc(OC)cc4)coc23)CC1 ZINC41586146",
      "Cc1ccc(NC(=O)NC2C=C(C(=O)NC3CCCNC3=O)CC(O)C2O)cc1 ZINC03840146",
      "CC1CC(C)C[NH+](Cc2c([O-])ccc3c2O/C(=C\\C=C\\c2ccccc2)C3=O)C1 ZINC41583797",
      "CSCCC(NC(=O)C(Cc1ccccc1)NC(=O)N1CC(=O)Nc2ccccc21)C(=O)[O-] ZINC12889537",
      "CCC(C)C(NC(=O)C1(c2ccccc2)CCN(C(=O)OC(C)(C)C)CC1)C(=O)NC(CC(C)C)C(=O)[O-] ZINC12887538",
      "COc1ccc(C(=O)NCc2cc(=O)oc3cc(C)c(Cl)c(C)c23)cc1 ZINC41285557",
      "C[N+](C)=CC1=C(O)/C(=C/c2ccccc2)CCC1 ZINC35614474",
      "COc1ccc(OC)c(-c2c(C)c3ccc(OC)cc3oc2=O)c1 ZINC04085533",
      "c1ccc(COC2COC3C([NH2+]Cc4ccc(-c5ccccc5)cc4)COC23)cc1 ZINC05397166",
      "NC1=C(C(=O)Nc2ccccc2)C(c2ccc(Br)cc2)C2Oc3ccccc3C2O1 ZINC03851792",
      "O=C(Nc1nccs1)C1C[NH+]2CCC1CC2Cn1cc(CCO)nn1 ZINC08298377",
      "C=C(C)COc1ccc2oc(-c3ccc(OC)cc3)c(C(=O)NCc3ccccc3)c2c1 ZINC15880519",
      "CCC(C)C(NC(=O)N1CCc2cc(OC)c(OC)cc2C1)C(=O)Nc1nccs1 ZINC12896765",
      "COc1ccc(S(=O)(=O)NCC2CC3CC[NH+]2CC3C[NH+]2CCN(c3ccccc3)CC2)cc1 ZINC08636242",
      "CCOC(=O)C1CCCN(C(=O)Cc2c(C)c3ccc(O)c(O)c3oc2=O)C1 ZINC12894801",
      "CC(C)c1ccc2c(c1)C(O)CC1C(C)(CO)CCCC21C ZINC13302244",
      "CCc1c(O)c(=O)ccn1C1OC(CO)C(O)C1O ZINC03847032",
      "CC1=C(C(=O)OC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(O)C2O)C(C)(C)CCC1 ZINC35465935",
      "Cc1cc([O-])c(C[NH+]2CCCC2)c2c1C(=O)/C(=C/c1ccc3c(c1)OCO3)O2 ZINC41513428",
      "Cc1c(C)c2ccc(OC3OC(CO)C(O)C(O)C3O)cc2oc1=O ZINC04064983",
      "COc1ccc(/C=C2\\Oc3c(C[NH+]4CCN(c5ccccc5)CC4)c([O-])ccc3C2=O)c(OC)c1OC ZINC41585163",
      "C[NH+](C)Cc1ccccc1-c1ccc2n(c1=O)CC1CN(C(=O)c3ccco3)CC2C1 ZINC04237216",
      "CCOC(=O)C1=CCCN(C)C1 ZINC31555018",
      "C=CCN(C)C(=O)C(C)C1CCC2(C)Cc3sc(NC(=O)CC(C)(C)C)nc3C(C)C2C1O ZINC03840722",
      "COC(=O)CNC(=O)N1CCc2[nH]c[nH+]c2C1c1ccccc1OC ZINC20763158",
      "CC(=O)NCC1CC2CC[NH+]1CC2c1cc(-c2ccccc2Br)nc(C)n1 ZINC08623367",
      "CC12C=CC(=O)C=C1CCC1C2CCC2(C)C(O)CCC12 ZINC04073880",
      "[N-]=[N+]=CC(=O)OCC([NH3+])C(=O)[O-] ZINC03977741"};
  std::vector<std::string> csmileses;
  for (const auto &smi : smileses) {
    auto m = SmilesParse::MolFromSmiles(smi);
    REQUIRE(m);
    csmileses.push_back(MolToSmiles(*m));
  }
  // quick expansion to get to 800 structures
  for (auto i = 0u; i < 4; ++i) {
    smileses.insert(smileses.end(), smileses.begin(), smileses.end());
    csmileses.insert(csmileses.end(), csmileses.begin(), csmileses.end());
  }
  size_t nThreads = 4;
  std::vector<std::thread> threads;
  for (auto i = 0u; i < nThreads; ++i) {
    auto func = [&smileses, &csmileses](unsigned) {
      for (auto j = 0u; j < smileses.size(); j++) {
        auto m = SmilesParse::MolFromSmiles(smileses[j]);
        REQUIRE(m);
        CHECK(MolToSmiles(*m) == csmileses[j]);
      }
    };
    threads.emplace_back(func, i);
  }
  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}
#endif
