//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Tautomer.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;
using namespace MolStandardize;

void testEnumerator() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing tautomer enumeration" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string tautomerFile =
      rdbase + "/Data/MolStandardize/tautomerTransforms.in";
  auto tautparams = std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(tautomerFile));
  unsigned int ntautomers = tautparams->getNumTautomers();
  TEST_ASSERT(ntautomers == 34);

  TautomerCatalog tautcat(tautparams.get());
  TautomerEnumerator te;

  // Enumerate 1,3 keto/enol tautomer.
  std::string smi1 = "C1(=CCCCC1)O";
  std::shared_ptr<ROMol> m1(SmilesToMol(smi1));
  std::vector<ROMOL_SPTR> res = te.enumerate(*m1, &tautcat);
  std::vector<std::string> ans = {"O=C1CCCCC1", "OC1=CCCCC1"};
  TEST_ASSERT(res.size() == ans.size());
  for (size_t i = 0; i < res.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res[i]) == ans[i]);
  }

  // Enumerate 1,3 keto/enol tautomer.
  std::string smi2 = "C1(CCCCC1)=O";
  std::shared_ptr<ROMol> m2(SmilesToMol(smi2));
  std::vector<ROMOL_SPTR> res2 = te.enumerate(*m2, &tautcat);
  std::vector<std::string> ans2 = {"O=C1CCCCC1", "OC1=CCCCC1"};
  TEST_ASSERT(res2.size() == ans2.size());
  for (size_t i = 0; i < res2.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res2[i]) == ans2[i]);
  }

  // Enumerate acetophenone keto/enol tautomer.
  std::string smi3 = "C(=C)(O)C1=CC=CC=C1";
  std::shared_ptr<ROMol> m3(SmilesToMol(smi3));
  std::vector<ROMOL_SPTR> res3 = te.enumerate(*m3, &tautcat);
  std::vector<std::string> ans3 = {"C=C(O)c1ccccc1", "CC(=O)c1ccccc1"};
  TEST_ASSERT(res3.size() == ans3.size());
  for (size_t i = 0; i < res3.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res3[i]) == ans3[i]);
  }

  // Enumerate acetone keto/enol tautomer.
  std::string smi4 = "CC(C)=O";
  std::shared_ptr<ROMol> m4(SmilesToMol(smi4));
  std::vector<ROMOL_SPTR> res4 = te.enumerate(*m4, &tautcat);
  std::vector<std::string> ans4 = {"C=C(C)O", "CC(C)=O"};
  TEST_ASSERT(res4.size() == ans4.size());
  for (size_t i = 0; i < res4.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res4[i]) == ans4[i]);
  }

  // keto/enol tautomer
  std::string smi5 = "OC(C)=C(C)C";
  std::shared_ptr<ROMol> m5(SmilesToMol(smi5));
  std::vector<ROMOL_SPTR> res5 = te.enumerate(*m5, &tautcat);
  std::vector<std::string> ans5 = {"C=C(O)C(C)C", "CC(=O)C(C)C", "CC(C)=C(C)O"};
  TEST_ASSERT(res5.size() == ans5.size());
  for (size_t i = 0; i < res5.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res5[i]) == ans5[i]);
  }

  // 1-phenyl-2-propanone enol/keto
  std::string smi6 = "c1(ccccc1)CC(=O)C";
  std::shared_ptr<ROMol> m6(SmilesToMol(smi6));
  std::vector<ROMOL_SPTR> res6 = te.enumerate(*m6, &tautcat);
  std::vector<std::string> ans6 = {"C=C(O)Cc1ccccc1", "CC(=O)Cc1ccccc1",
                                   "CC(O)=Cc1ccccc1"};
  TEST_ASSERT(res6.size() == ans6.size());
  for (size_t i = 0; i < res6.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res6[i]) == ans6[i]);
  }

  // 1,5 keto/enol tautomer
  std::string smi7 = "Oc1nccc2cc[nH]c(=N)c12";
  std::shared_ptr<ROMol> m7(SmilesToMol(smi7));
  std::vector<ROMOL_SPTR> res7 = te.enumerate(*m7, &tautcat);
  std::vector<std::string> ans7 = {
      "N=c1[nH]ccc2cc[nH]c(=O)c12", "N=c1[nH]ccc2ccnc(O)c12",
      "N=c1nccc2cc[nH]c(O)c1-2",    "Nc1[nH]ccc2ccnc(=O)c1-2",
      "Nc1nccc2cc[nH]c(=O)c12",     "Nc1nccc2ccnc(O)c12"};
  TEST_ASSERT(res7.size() == ans7.size());
  for (size_t i = 0; i < res7.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res7[i]) == ans7[i]);
  }

  // 1,5 keto/enol tautomer
  std::string smi8 = "C1(C=CCCC1)=O";
  std::shared_ptr<ROMol> m8(SmilesToMol(smi8));
  std::vector<ROMOL_SPTR> res8 = te.enumerate(*m8, &tautcat);
  std::vector<std::string> ans8 = {"O=C1C=CCCC1", "O=C1CC=CCC1", "OC1=CC=CCC1",
                                   "OC1=CCC=CC1", "OC1=CCCC=C1"};
  TEST_ASSERT(res8.size() == ans8.size());
  for (size_t i = 0; i < res8.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res8[i]) == ans8[i]);
  }

  // 1,5 keto/enol tautomer
  std::string smi9 = "C1(=CC=CCC1)O";
  std::shared_ptr<ROMol> m9(SmilesToMol(smi9));
  std::vector<ROMOL_SPTR> res9 = te.enumerate(*m9, &tautcat);
  std::vector<std::string> ans9 = {"O=C1C=CCCC1", "O=C1CC=CCC1", "OC1=CC=CCC1",
                                   "OC1=CCC=CC1", "OC1=CCCC=C1"};
  TEST_ASSERT(res9.size() == ans9.size());
  for (size_t i = 0; i < res9.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res9[i]) == ans9[i]);
  }

  // aliphatic imine tautomer
  std::string smi10 = "C1(CCCCC1)=N";
  std::shared_ptr<ROMol> m10(SmilesToMol(smi10));
  std::vector<ROMOL_SPTR> res10 = te.enumerate(*m10, &tautcat);
  std::vector<std::string> ans10 = {"N=C1CCCCC1", "NC1=CCCCC1"};
  TEST_ASSERT(res10.size() == ans10.size());
  for (size_t i = 0; i < res10.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res10[i]) == ans10[i]);
  }

  // aliphatic imine tautomer
  std::string smi11 = "C1(=CCCCC1)N";
  std::shared_ptr<ROMol> m11(SmilesToMol(smi11));
  std::vector<ROMOL_SPTR> res11 = te.enumerate(*m11, &tautcat);
  std::vector<std::string> ans11 = {"N=C1CCCCC1", "NC1=CCCCC1"};
  TEST_ASSERT(res11.size() == ans11.size());
  for (size_t i = 0; i < res11.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res11[i]) == ans11[i]);
  }

  // special imine tautomer
  std::string smi12 = "C1(C=CC=CN1)=CC";
  std::shared_ptr<ROMol> m12(SmilesToMol(smi12));
  std::vector<ROMOL_SPTR> res12 = te.enumerate(*m12, &tautcat);
  std::vector<std::string> ans12 = {"CC=C1C=CC=CN1", "CC=C1C=CCC=N1",
                                    "CCc1ccccn1"};
  TEST_ASSERT(res12.size() == ans12.size());
  for (size_t i = 0; i < res12.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res12[i]) == ans12[i]);
  }

  // special imine tautomer
  std::string smi13 = "C1(=NC=CC=C1)CC";
  std::shared_ptr<ROMol> m13(SmilesToMol(smi13));
  std::vector<ROMOL_SPTR> res13 = te.enumerate(*m13, &tautcat);
  std::vector<std::string> ans13 = {"CC=C1C=CC=CN1", "CC=C1C=CCC=N1",
                                    "CCc1ccccn1"};
  TEST_ASSERT(res13.size() == ans13.size());
  for (size_t i = 0; i < res13.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res13[i]) == ans13[i]);
  }

  // 1,3 aromatic heteroatom H shift
  std::string smi14 = "O=c1cccc[nH]1";
  std::shared_ptr<ROMol> m14(SmilesToMol(smi14));
  std::vector<ROMOL_SPTR> res14 = te.enumerate(*m14, &tautcat);
  std::vector<std::string> ans14 = {"O=c1cccc[nH]1", "Oc1ccccn1"};
  TEST_ASSERT(res14.size() == ans14.size());
  for (size_t i = 0; i < res14.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res14[i]) == ans14[i]);
  }

  // 1,3 aromatic heteroatom H shift
  std::string smi15 = "Oc1ccccn1";
  std::shared_ptr<ROMol> m15(SmilesToMol(smi15));
  std::vector<ROMOL_SPTR> res15 = te.enumerate(*m15, &tautcat);
  std::vector<std::string> ans15 = {"O=c1cccc[nH]1", "Oc1ccccn1"};
  TEST_ASSERT(res15.size() == ans15.size());
  for (size_t i = 0; i < res15.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res15[i]) == ans15[i]);
  }

  // 1,3 aromatic heteroatom H shift
  std::string smi16 = "Oc1ncc[nH]1";
  std::shared_ptr<ROMol> m16(SmilesToMol(smi16));
  std::vector<ROMOL_SPTR> res16 = te.enumerate(*m16, &tautcat);
  std::vector<std::string> ans16 = {"O=c1[nH]cc[nH]1", "Oc1ncc[nH]1"};
  TEST_ASSERT(res16.size() == ans16.size());
  for (size_t i = 0; i < res16.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res16[i]) == ans16[i]);
  }

  // 1,3 heteroatom H shift
  std::string smi17 = "OC(C)=NC";
  std::shared_ptr<ROMol> m17(SmilesToMol(smi17));
  std::vector<ROMOL_SPTR> res17 = te.enumerate(*m17, &tautcat);
  std::vector<std::string> ans17 = {"C=C(O)NC", "CN=C(C)O", "CNC(C)=O"};
  TEST_ASSERT(res17.size() == ans17.size());
  for (size_t i = 0; i < res17.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res17[i]) == ans17[i]);
  }

  // 1,3 heteroatom H shift
  std::string smi18 = "CNC(C)=O";
  std::shared_ptr<ROMol> m18(SmilesToMol(smi18));
  std::vector<ROMOL_SPTR> res18 = te.enumerate(*m18, &tautcat);
  std::vector<std::string> ans18 = {"C=C(O)NC", "CN=C(C)O", "CNC(C)=O"};
  TEST_ASSERT(res18.size() == ans18.size());
  for (size_t i = 0; i < res18.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res18[i]) == ans18[i]);
  }

  // 1,3 heteroatom H shift
  std::string smi19 = "S=C(N)N";
  std::shared_ptr<ROMol> m19(SmilesToMol(smi19));
  std::vector<ROMOL_SPTR> res19 = te.enumerate(*m19, &tautcat);
  std::vector<std::string> ans19 = {"N=C(N)S", "NC(N)=S"};
  TEST_ASSERT(res19.size() == ans19.size());
  for (size_t i = 0; i < res19.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res19[i]) == ans19[i]);
  }

  // 1,3 heteroatom H shift
  std::string smi20 = "SC(N)=N";
  std::shared_ptr<ROMol> m20(SmilesToMol(smi20));
  std::vector<ROMOL_SPTR> res20 = te.enumerate(*m20, &tautcat);
  std::vector<std::string> ans20 = {"N=C(N)S", "NC(N)=S"};
  TEST_ASSERT(res20.size() == ans20.size());
  for (size_t i = 0; i < res20.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res20[i]) == ans20[i]);
  }

  // 1,3 heteroatom H shift
  std::string smi21 = "N=c1[nH]ccn(C)1";
  std::shared_ptr<ROMol> m21(SmilesToMol(smi21));
  std::vector<ROMOL_SPTR> res21 = te.enumerate(*m21, &tautcat);
  std::vector<std::string> ans21 = {"Cn1cc[nH]c1=N", "Cn1ccnc1N"};
  TEST_ASSERT(res21.size() == ans21.size());
  for (size_t i = 0; i < res21.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res21[i]) == ans21[i]);
  }

  // 1,3 heteroatom H shift
  std::string smi22 = "CN=c1[nH]cncc1";
  std::shared_ptr<ROMol> m22(SmilesToMol(smi22));
  std::vector<ROMOL_SPTR> res22 = te.enumerate(*m22, &tautcat);
  std::vector<std::string> ans22 = {"CN=c1cc[nH]cn1", "CN=c1ccnc[nH]1",
                                    "CNc1ccncn1"};
  TEST_ASSERT(res22.size() == ans22.size());
  for (size_t i = 0; i < res22.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res22[i]) == ans22[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi23 = "Oc1cccc2ccncc12";
  std::shared_ptr<ROMol> m23(SmilesToMol(smi23));
  std::vector<ROMOL_SPTR> res23 = te.enumerate(*m23, &tautcat);
  std::vector<std::string> ans23 = {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"};
  TEST_ASSERT(res23.size() == ans23.size());
  for (size_t i = 0; i < res23.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res23[i]) == ans23[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi24 = "O=c1cccc2cc[nH]cc1-2";
  std::shared_ptr<ROMol> m24(SmilesToMol(smi24));
  std::vector<ROMOL_SPTR> res24 = te.enumerate(*m24, &tautcat);
  std::vector<std::string> ans24 = {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"};
  TEST_ASSERT(res24.size() == ans24.size());
  for (size_t i = 0; i < res24.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res24[i]) == ans24[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi25 = "Cc1n[nH]c2ncnn12";
  std::shared_ptr<ROMol> m25(SmilesToMol(smi25));
  std::vector<ROMOL_SPTR> res25 = te.enumerate(*m25, &tautcat);
  std::vector<std::string> ans25 = {"C=C1NN=C2N=CNN12", "C=C1NN=C2NC=NN12",
                                    "C=C1NNc2ncnn21",   "Cc1n[nH]c2ncnn12",
                                    "Cc1nnc2[nH]cnn12", "Cc1nnc2nc[nH]n12"};
  TEST_ASSERT(res25.size() == ans25.size());
  for (size_t i = 0; i < res25.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res25[i]) == ans25[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi26 = "Cc1nnc2nc[nH]n12";
  std::shared_ptr<ROMol> m26(SmilesToMol(smi26));
  std::vector<ROMOL_SPTR> res26 = te.enumerate(*m26, &tautcat);
  std::vector<std::string> ans26 = {"C=C1NN=C2N=CNN12", "C=C1NN=C2NC=NN12",
                                    "C=C1NNc2ncnn21",   "Cc1n[nH]c2ncnn12",
                                    "Cc1nnc2[nH]cnn12", "Cc1nnc2nc[nH]n12"};
  TEST_ASSERT(res26.size() == ans26.size());
  for (size_t i = 0; i < res26.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res26[i]) == ans26[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi29 = "Oc1ccncc1";
  std::shared_ptr<ROMol> m29(SmilesToMol(smi29));
  std::vector<ROMOL_SPTR> res29 = te.enumerate(*m29, &tautcat);
  std::vector<std::string> ans29 = {"O=c1cc[nH]cc1", "Oc1ccncc1"};
  TEST_ASSERT(res29.size() == ans29.size());
  for (size_t i = 0; i < res29.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res29[i]) == ans29[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi27 = "Oc1c(cccc3)c3nc2ccncc12";
  std::shared_ptr<ROMol> m27(SmilesToMol(smi27));
  std::vector<ROMOL_SPTR> res27 = te.enumerate(*m27, &tautcat);
  std::vector<std::string> ans27 = {"O=c1c2c[nH]ccc-2nc2ccccc12",
                                    "O=c1c2ccccc2[nH]c2ccncc12",
                                    "Oc1c2ccccc2nc2ccncc12"};
  TEST_ASSERT(res27.size() == ans27.size());
  for (size_t i = 0; i < res27.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res27[i]) == ans27[i]);
  }

  // 1,3 and 1,5 aromatic heteroatom H shift
  std::string smi28 = "Oc1ncncc1";
  std::shared_ptr<ROMol> m28(SmilesToMol(smi28));
  std::vector<ROMOL_SPTR> res28 = te.enumerate(*m28, &tautcat);
  std::vector<std::string> ans28 = {"O=c1cc[nH]cn1", "O=c1ccnc[nH]1",
                                    "Oc1ccncn1"};
  TEST_ASSERT(res28.size() == ans28.size());
  for (size_t i = 0; i < res28.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res28[i]) == ans28[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi30 = "C2(=C1C(=NC=N1)[NH]C(=N2)N)O";
  std::shared_ptr<ROMol> m30(SmilesToMol(smi30));
  std::vector<ROMOL_SPTR> res30 = te.enumerate(*m30, &tautcat);
  std::vector<std::string> ans30 = {
      "N=c1[nH]c(=O)c2[nH]cnc2[nH]1", "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
      "N=c1[nH]c2ncnc-2c(O)[nH]1",    "N=c1nc(O)c2[nH]cnc2[nH]1",
      "N=c1nc(O)c2nc[nH]c2[nH]1",     "N=c1nc2[nH]cnc2c(O)[nH]1",
      "N=c1nc2nc[nH]c2c(O)[nH]1",     "Nc1nc(=O)c2[nH]cnc2[nH]1",
      "Nc1nc(=O)c2nc[nH]c2[nH]1",     "Nc1nc(O)c2[nH]cnc2n1",
      "Nc1nc(O)c2nc[nH]c2n1",         "Nc1nc(O)c2ncnc-2[nH]1",
      "Nc1nc2[nH]cnc2c(=O)[nH]1",     "Nc1nc2nc[nH]c2c(=O)[nH]1",
      "Nc1nc2ncnc-2c(O)[nH]1"};
  TEST_ASSERT(res30.size() == ans30.size());
  for (size_t i = 0; i < res30.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res30[i]) == ans30[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi31 = "C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O";
  std::shared_ptr<ROMol> m31(SmilesToMol(smi31));
  std::vector<ROMOL_SPTR> res31 = te.enumerate(*m31, &tautcat);
  std::vector<std::string> ans31 = {
      "N=c1[nH]c(=O)c2[nH]cnc2[nH]1", "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
      "N=c1[nH]c2ncnc-2c(O)[nH]1",    "N=c1nc(O)c2[nH]cnc2[nH]1",
      "N=c1nc(O)c2nc[nH]c2[nH]1",     "N=c1nc2[nH]cnc2c(O)[nH]1",
      "N=c1nc2nc[nH]c2c(O)[nH]1",     "Nc1nc(=O)c2[nH]cnc2[nH]1",
      "Nc1nc(=O)c2nc[nH]c2[nH]1",     "Nc1nc(O)c2[nH]cnc2n1",
      "Nc1nc(O)c2nc[nH]c2n1",         "Nc1nc(O)c2ncnc-2[nH]1",
      "Nc1nc2[nH]cnc2c(=O)[nH]1",     "Nc1nc2nc[nH]c2c(=O)[nH]1",
      "Nc1nc2ncnc-2c(O)[nH]1"};
  TEST_ASSERT(res31.size() == ans31.size());
  for (size_t i = 0; i < res31.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res31[i]) == ans31[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi32 = "Oc1n(C)ncc1";
  std::shared_ptr<ROMol> m32(SmilesToMol(smi32));
  std::vector<ROMOL_SPTR> res32 = te.enumerate(*m32, &tautcat);
  std::vector<std::string> ans32 = {"CN1N=CCC1=O", "Cn1[nH]ccc1=O",
                                    "Cn1nccc1O"};
  TEST_ASSERT(res32.size() == ans32.size());
  for (size_t i = 0; i < res32.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res32[i]) == ans32[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi33 = "O=c1nc2[nH]ccn2cc1";
  std::shared_ptr<ROMol> m33(SmilesToMol(smi33));
  std::vector<ROMOL_SPTR> res33 = te.enumerate(*m33, &tautcat);
  std::vector<std::string> ans33 = {"O=c1ccn2cc[nH]c2n1", "O=c1ccn2ccnc2[nH]1",
                                    "Oc1ccn2ccnc2n1"};
  TEST_ASSERT(res33.size() == ans33.size());
  for (size_t i = 0; i < res33.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res33[i]) == ans33[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi34 = "N=c1nc[nH]cc1";
  std::shared_ptr<ROMol> m34(SmilesToMol(smi34));
  std::vector<ROMOL_SPTR> res34 = te.enumerate(*m34, &tautcat);
  std::vector<std::string> ans34 = {"N=c1cc[nH]cn1", "N=c1ccnc[nH]1",
                                    "Nc1ccncn1"};
  TEST_ASSERT(res34.size() == ans34.size());
  for (size_t i = 0; i < res34.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res34[i]) == ans34[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi35 = "N=c(c1)ccn2cc[nH]c12";
  std::shared_ptr<ROMol> m35(SmilesToMol(smi35));
  std::vector<ROMOL_SPTR> res35 = te.enumerate(*m35, &tautcat);
  std::vector<std::string> ans35 = {"N=c1ccn2cc[nH]c2c1", "Nc1ccn2ccnc2c1"};
  TEST_ASSERT(res35.size() == ans35.size());
  for (size_t i = 0; i < res35.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res35[i]) == ans35[i]);
  }

  // 1,5 aromatic heteroatom H shift
  std::string smi36 = "CN=c1nc[nH]cc1";
  std::shared_ptr<ROMol> m36(SmilesToMol(smi36));
  std::vector<ROMOL_SPTR> res36 = te.enumerate(*m36, &tautcat);
  std::vector<std::string> ans36 = {"CN=c1cc[nH]cn1", "CN=c1ccnc[nH]1",
                                    "CNc1ccncn1"};
  TEST_ASSERT(res36.size() == ans36.size());
  for (size_t i = 0; i < res36.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res36[i]) == ans36[i]);
  }

  // 1,7 aromatic heteroatom H shift
  std::string smi37 = "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1";
  std::shared_ptr<ROMol> m37(SmilesToMol(smi37));
  std::vector<ROMOL_SPTR> res37 = te.enumerate(*m37, &tautcat);
  std::vector<std::string> ans37 = {"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
                                    "c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2",
                                    "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2"};
  TEST_ASSERT(res37.size() == ans37.size());
  for (size_t i = 0; i < res37.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res37[i]) == ans37[i]);
  }

  // 1,7 aromatic heteroatom H shift
  std::string smi38 = "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2";
  std::shared_ptr<ROMol> m38(SmilesToMol(smi38));
  std::vector<ROMOL_SPTR> res38 = te.enumerate(*m38, &tautcat);
  std::vector<std::string> ans38 = {"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
                                    "c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2",
                                    "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2"};
  TEST_ASSERT(res38.size() == ans38.size());
  for (size_t i = 0; i < res38.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res38[i]) == ans38[i]);
  }

  // 1,9 aromatic heteroatom H shift
  std::string smi39 = "CNc1ccnc2ncnn21";
  std::shared_ptr<ROMol> m39(SmilesToMol(smi39));
  std::vector<ROMOL_SPTR> res39 = te.enumerate(*m39, &tautcat);
  std::vector<std::string> ans39 = {"CN=c1cc[nH]c2ncnn12",
                                    "CN=c1ccnc2[nH]cnn12",
                                    "CN=c1ccnc2nc[nH]n12", "CNc1ccnc2ncnn12"};
  TEST_ASSERT(res39.size() == ans39.size());
  for (size_t i = 0; i < res39.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res39[i]) == ans39[i]);
  }

  // 1,9 aromatic heteroatom H shift
  std::string smi40 = "CN=c1ccnc2nc[nH]n21";
  std::shared_ptr<ROMol> m40(SmilesToMol(smi40));
  std::vector<ROMOL_SPTR> res40 = te.enumerate(*m40, &tautcat);
  std::vector<std::string> ans40 = {"CN=c1cc[nH]c2ncnn12",
                                    "CN=c1ccnc2[nH]cnn12",
                                    "CN=c1ccnc2nc[nH]n12", "CNc1ccnc2ncnn12"};
  TEST_ASSERT(res40.size() == ans40.size());
  for (size_t i = 0; i < res40.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res40[i]) == ans40[i]);
  }

  // 1,11 aromatic heteroatom H shift
  std::string smi41 = "Nc1ccc(C=C2C=CC(=O)C=C2)cc1";
  std::shared_ptr<ROMol> m41(SmilesToMol(smi41));
  std::vector<ROMOL_SPTR> res41 = te.enumerate(*m41, &tautcat);
  std::vector<std::string> ans41 = {
      "N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1", "N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
      "N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1", "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"};
  TEST_ASSERT(res41.size() == ans41.size());
  for (size_t i = 0; i < res41.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res41[i]) == ans41[i]);
  }

  // 1,11 aromatic heteroatom H shift
  std::string smi42 = "N=C1C=CC(=Cc2ccc(O)cc2)C=C1";
  std::shared_ptr<ROMol> m42(SmilesToMol(smi42));
  std::vector<ROMOL_SPTR> res42 = te.enumerate(*m42, &tautcat);
  std::vector<std::string> ans42 = {
      "N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1", "N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
      "N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1", "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"};
  TEST_ASSERT(res42.size() == ans42.size());
  for (size_t i = 0; i < res42.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res42[i]) == ans42[i]);
  }

  // heterocyclic tautomer
  std::string smi43 = "n1ccc2ccc[nH]c12";
  std::shared_ptr<ROMol> m43(SmilesToMol(smi43));
  std::vector<ROMOL_SPTR> res43 = te.enumerate(*m43, &tautcat);
  std::vector<std::string> ans43 = {"c1c[nH]c2nccc-2c1", "c1cnc2[nH]ccc2c1"};
  TEST_ASSERT(res43.size() == ans43.size());
  for (size_t i = 0; i < res43.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res43[i]) == ans43[i]);
  }

  // heterocyclic tautomer
  std::string smi44 = "c1cc(=O)[nH]c2nccn12";
  std::shared_ptr<ROMol> m44(SmilesToMol(smi44));
  std::vector<ROMOL_SPTR> res44 = te.enumerate(*m44, &tautcat);
  std::vector<std::string> ans44 = {"O=c1ccn2cc[nH]c2n1", "O=c1ccn2ccnc2[nH]1",
                                    "Oc1ccn2ccnc2n1"};
  TEST_ASSERT(res44.size() == ans44.size());
  for (size_t i = 0; i < res44.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res44[i]) == ans44[i]);
  }

  // heterocyclic tautomer
  std::string smi45 = "c1cnc2c[nH]ccc12";
  std::shared_ptr<ROMol> m45(SmilesToMol(smi45));
  std::vector<ROMOL_SPTR> res45 = te.enumerate(*m45, &tautcat);
  std::vector<std::string> ans45 = {"c1cc2cc[nH]c2cn1", "c1cc2cc[nH]cc-2n1"};
  TEST_ASSERT(res45.size() == ans45.size());
  for (size_t i = 0; i < res45.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res45[i]) == ans45[i]);
  }

  // heterocyclic tautomer
  std::string smi46 = "n1ccc2c[nH]ccc12";
  std::shared_ptr<ROMol> m46(SmilesToMol(smi46));
  std::vector<ROMOL_SPTR> res46 = te.enumerate(*m46, &tautcat);
  std::vector<std::string> ans46 = {"c1cc2[nH]ccc2cn1", "c1cc2c[nH]ccc-2n1"};
  TEST_ASSERT(res46.size() == ans46.size());
  for (size_t i = 0; i < res46.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res46[i]) == ans46[i]);
  }

  // heterocyclic tautomer
  std::string smi47 = "c1cnc2ccc[nH]c12";
  std::shared_ptr<ROMol> m47(SmilesToMol(smi47));
  std::vector<ROMOL_SPTR> res47 = te.enumerate(*m47, &tautcat);
  std::vector<std::string> ans47 = {"c1c[nH]c2ccnc-2c1", "c1cnc2cc[nH]c2c1"};
  TEST_ASSERT(res47.size() == ans47.size());
  for (size_t i = 0; i < res47.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res47[i]) == ans47[i]);
  }

  // furanone tautomer
  std::string smi48 = "C1=CC=C(O1)O";
  std::shared_ptr<ROMol> m48(SmilesToMol(smi48));
  std::vector<ROMOL_SPTR> res48 = te.enumerate(*m48, &tautcat);
  std::vector<std::string> ans48 = {"O=C1CC=CO1", "Oc1ccco1"};
  TEST_ASSERT(res48.size() == ans48.size());
  for (size_t i = 0; i < res48.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res48[i]) == ans48[i]);
  }

  // furanone tautomer
  std::string smi49 = "O=C1CC=CO1";
  std::shared_ptr<ROMol> m49(SmilesToMol(smi49));
  std::vector<ROMOL_SPTR> res49 = te.enumerate(*m49, &tautcat);
  std::vector<std::string> ans49 = {"O=C1CC=CO1", "Oc1ccco1"};
  TEST_ASSERT(res49.size() == ans49.size());
  for (size_t i = 0; i < res49.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res49[i]) == ans49[i]);
  }

  // keten/ynol tautomer
  std::string smi50 = "CC=C=O";
  std::shared_ptr<ROMol> m50(SmilesToMol(smi50));
  std::vector<ROMOL_SPTR> res50 = te.enumerate(*m50, &tautcat);
  std::vector<std::string> ans50 = {"CC#CO", "CC=C=O"};
  TEST_ASSERT(res50.size() == ans50.size());
  for (size_t i = 0; i < res50.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res50[i]) == ans50[i]);
  }

  // keten/ynol tautomer
  std::string smi51 = "CC#CO";
  std::shared_ptr<ROMol> m51(SmilesToMol(smi51));
  std::vector<ROMOL_SPTR> res51 = te.enumerate(*m51, &tautcat);
  std::vector<std::string> ans51 = {"CC#CO", "CC=C=O"};
  TEST_ASSERT(res51.size() == ans51.size());
  for (size_t i = 0; i < res51.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res51[i]) == ans51[i]);
  }

  // ionic nitro/aci-nitro tautomer
  std::string smi52 = "C([N+](=O)[O-])C";
  std::shared_ptr<ROMol> m52(SmilesToMol(smi52));
  std::vector<ROMOL_SPTR> res52 = te.enumerate(*m52, &tautcat);
  std::vector<std::string> ans52 = {"CC=[N+]([O-])O", "CC[N+](=O)[O-]"};
  TEST_ASSERT(res52.size() == ans52.size());
  for (size_t i = 0; i < res52.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res52[i]) == ans52[i]);
  }

  // ionic nitro/aci-nitro tautomer
  std::string smi53 = "C(=[N+](O)[O-])C";
  std::shared_ptr<ROMol> m53(SmilesToMol(smi53));
  std::vector<ROMOL_SPTR> res53 = te.enumerate(*m53, &tautcat);
  std::vector<std::string> ans53 = {"CC=[N+]([O-])O", "CC[N+](=O)[O-]"};
  TEST_ASSERT(res53.size() == ans53.size());
  for (size_t i = 0; i < res53.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res53[i]) == ans53[i]);
  }

  // oxim nitroso tautomer
  std::string smi54 = "CC(C)=NO";
  std::shared_ptr<ROMol> m54(SmilesToMol(smi54));
  std::vector<ROMOL_SPTR> res54 = te.enumerate(*m54, &tautcat);
  std::vector<std::string> ans54 = {"C=C(C)NO", "CC(C)=NO", "CC(C)N=O"};
  TEST_ASSERT(res54.size() == ans54.size());
  for (size_t i = 0; i < res54.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res54[i]) == ans54[i]);
  }

  // oxim nitroso tautomer
  std::string smi55 = "CC(C)N=O";
  std::shared_ptr<ROMol> m55(SmilesToMol(smi55));
  std::vector<ROMOL_SPTR> res55 = te.enumerate(*m55, &tautcat);
  std::vector<std::string> ans55 = {"C=C(C)NO", "CC(C)=NO", "CC(C)N=O"};
  TEST_ASSERT(res55.size() == ans55.size());
  for (size_t i = 0; i < res55.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res55[i]) == ans55[i]);
  }

  // oxim/nitroso tautomer via phenol
  std::string smi56 = "O=Nc1ccc(O)cc1";
  std::shared_ptr<ROMol> m56(SmilesToMol(smi56));
  std::vector<ROMOL_SPTR> res56 = te.enumerate(*m56, &tautcat);
  std::vector<std::string> ans56 = {"O=C1C=CC(=NO)C=C1", "O=NC1C=CC(=O)C=C1",
                                    "O=Nc1ccc(O)cc1"};
  TEST_ASSERT(res56.size() == ans56.size());
  for (size_t i = 0; i < res56.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res56[i]) == ans56[i]);
  }

  // oxim/nitroso tautomer via phenol
  std::string smi57 = "O=C1C=CC(=NO)C=C1";
  std::shared_ptr<ROMol> m57(SmilesToMol(smi57));
  std::vector<ROMOL_SPTR> res57 = te.enumerate(*m57, &tautcat);
  std::vector<std::string> ans57 = {"O=C1C=CC(=NO)C=C1", "O=NC1C=CC(=O)C=C1",
                                    "O=Nc1ccc(O)cc1"};
  TEST_ASSERT(res57.size() == ans57.size());
  for (size_t i = 0; i < res57.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res57[i]) == ans57[i]);
  }

  // cyano/iso-cyanic acid tautomer
  std::string smi58 = "C(#N)O";
  std::shared_ptr<ROMol> m58(SmilesToMol(smi58));
  std::vector<ROMOL_SPTR> res58 = te.enumerate(*m58, &tautcat);
  std::vector<std::string> ans58 = {"N#CO", "N=C=O"};
  TEST_ASSERT(res58.size() == ans58.size());
  for (size_t i = 0; i < res58.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res58[i]) == ans58[i]);
  }

  // cyano/iso-cyanic acid tautomer
  std::string smi59 = "C(=N)=O";
  std::shared_ptr<ROMol> m59(SmilesToMol(smi59));
  std::vector<ROMOL_SPTR> res59 = te.enumerate(*m59, &tautcat);
  std::vector<std::string> ans59 = {"N#CO", "N=C=O"};
  TEST_ASSERT(res59.size() == ans59.size());
  for (size_t i = 0; i < res59.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res59[i]) == ans59[i]);
  }

  // isocyanide tautomer
  std::string smi60 = "C#N";
  std::shared_ptr<ROMol> m60(SmilesToMol(smi60));
  std::vector<ROMOL_SPTR> res60 = te.enumerate(*m60, &tautcat);
  std::vector<std::string> ans60 = {"C#N", "[C-]#[NH+]"};
  TEST_ASSERT(res60.size() == ans60.size());
  for (size_t i = 0; i < res60.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res60[i]) == ans60[i]);
  }

  // isocyanide tautomer
  std::string smi61 = "[C-]#[NH+]";
  std::shared_ptr<ROMol> m61(SmilesToMol(smi61));
  std::vector<ROMOL_SPTR> res61 = te.enumerate(*m61, &tautcat);
  std::vector<std::string> ans61 = {"C#N", "[C-]#[NH+]"};
  TEST_ASSERT(res61.size() == ans61.size());
  for (size_t i = 0; i < res61.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res61[i]) == ans61[i]);
  }

  // phosphonic acid tautomer
  std::string smi62 = "[PH](=O)(O)(O)";
  std::shared_ptr<ROMol> m62(SmilesToMol(smi62));
  std::vector<ROMOL_SPTR> res62 = te.enumerate(*m62, &tautcat);
  std::vector<std::string> ans62 = {"O=[PH](O)O", "OP(O)O"};
  TEST_ASSERT(res62.size() == ans62.size());
  for (size_t i = 0; i < res62.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res62[i]) == ans62[i]);
  }

  // phosphonic acid tautomer
  std::string smi63 = "P(O)(O)O";
  std::shared_ptr<ROMol> m63(SmilesToMol(smi63));
  std::vector<ROMOL_SPTR> res63 = te.enumerate(*m63, &tautcat);
  std::vector<std::string> ans63 = {"O=[PH](O)O", "OP(O)O"};
  TEST_ASSERT(res63.size() == ans63.size());
  for (size_t i = 0; i < res63.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res63[i]) == ans63[i]);
  }

  // Remove stereochemistry from mobile double bonds
  std::string smi64 = "c1(ccccc1)/C=C(/O)\\C";
  std::shared_ptr<ROMol> m64(SmilesToMol(smi64));
  std::vector<ROMOL_SPTR> res64 = te.enumerate(*m64, &tautcat);
  std::vector<std::string> ans64 = {"CC(O)=Cc1ccccc1", "C=C(O)Cc1ccccc1",
                                    "CC(=O)Cc1ccccc1"};
  TEST_ASSERT(res64.size() == ans64.size());
  for (size_t i = 0; i < res64.size(); ++i) {
    // std::cout << MolToSmiles(*res64[i]) << ", " << ans64[i] << std::endl;
    TEST_ASSERT(MolToSmiles(*res64[i]) == ans64[i]);
  }

  // Remove stereochemistry from mobile double bonds
  std::string smi65 = "C/C=C/C(C)=O";
  std::shared_ptr<ROMol> m65(SmilesToMol(smi65));
  std::vector<ROMOL_SPTR> res65 = te.enumerate(*m65, &tautcat);
  std::vector<std::string> ans65 = {"CC=CC(C)=O", "C=C(O)C=CC", "C=CC=C(C)O",
                                    "C=CCC(=C)O", "C=CCC(C)=O"};
  TEST_ASSERT(res65.size() == ans65.size());
  for (size_t i = 0; i < res65.size(); ++i) {
    TEST_ASSERT(MolToSmiles(*res65[i]) == ans65[i]);
  }

  // Remove stereochemistry from mobile double bonds
  std::string smi66 = "C/C=C\\C(C)=O";
  std::shared_ptr<ROMol> m66(SmilesToMol(smi66));
  std::vector<ROMOL_SPTR> res66 = te.enumerate(*m66, &tautcat);
  std::vector<std::string> ans66 = {"C=C(O)C=CC", "C=CC=C(C)O", "CC=CC(C)=O",
                                    "C=CCC(=C)O", "C=CCC(C)=O"};
  TEST_ASSERT(res66.size() == ans66.size());

  std::vector<std::string> sm66;
  for (const auto r : res66) {
    sm66.push_back(MolToSmiles(*r));
  }
  // sort both for alphabetical order
  std::sort(sm66.begin(), sm66.end());
  std::sort(ans66.begin(), ans66.end());
  TEST_ASSERT(sm66 == ans66);

  // Gaunine tautomers
  std::string smi67 = "N1C(N)=NC=2N=CNC2C1=O";
  std::shared_ptr<ROMol> m67(SmilesToMol(smi67));
  std::vector<ROMOL_SPTR> res67 = te.enumerate(*m67, &tautcat);
  std::vector<std::string> ans67 = {
      "N=c1[nH]c(=O)c2[nH]cnc2[nH]1", "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
      "N=c1[nH]c2ncnc-2c(O)[nH]1",    "N=c1nc(O)c2[nH]cnc2[nH]1",
      "N=c1nc(O)c2nc[nH]c2[nH]1",     "N=c1nc2[nH]cnc2c(O)[nH]1",
      "N=c1nc2nc[nH]c2c(O)[nH]1",     "Nc1nc(=O)c2[nH]cnc2[nH]1",
      "Nc1nc(=O)c2nc[nH]c2[nH]1",     "Nc1nc(O)c2[nH]cnc2n1",
      "Nc1nc(O)c2nc[nH]c2n1",         "Nc1nc(O)c2ncnc-2[nH]1",
      "Nc1nc2[nH]cnc2c(=O)[nH]1",     "Nc1nc2nc[nH]c2c(=O)[nH]1",
      "Nc1nc2ncnc-2c(O)[nH]1"};
  TEST_ASSERT(res67.size() == ans67.size());
  std::vector<std::string> sm67;
  for (const auto r : res67) {
    sm67.push_back(MolToSmiles(*r));
  }
  // sort both for alphabetical order
  std::sort(sm67.begin(), sm67.end());
  std::sort(ans67.begin(), ans67.end());
  TEST_ASSERT(sm67 == ans67);

  // Test a structure with hundreds of tautomers.
  std::string smi68 = "[H][C](CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O";
  std::shared_ptr<ROMol> m68(SmilesToMol(smi68));
  std::vector<ROMOL_SPTR> res68 = te.enumerate(*m68, &tautcat);
  TEST_ASSERT(res68.size() == 375);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testCanonicalize() {
  std::string rdbase = getenv("RDBASE");
  std::string tautomerFile =
      rdbase + "/Data/MolStandardize/tautomerTransforms.in";
  auto *tautparams = new TautomerCatalogParams(tautomerFile);
  unsigned int ntautomers = tautparams->getNumTautomers();
  TEST_ASSERT(ntautomers == 34);

  TautomerCatalog tautcat(tautparams);
  // TautomerCanonicalizer tc;

  // Enumerate 1,3 keto/enol tautomer.
  std::string smi1 = "C1(=CCCCC1)O";
  std::shared_ptr<ROMol> m1(SmilesToMol(smi1));
  // ROMol *res = tc.canonicalize(*m1, &tautcat);
  //    TEST_ASSERT(MolToSmiles(*res) == "O=C1CCCCC1");
}

int main() {
  testEnumerator();
  // testCanonicalize();
  return 0;
}
