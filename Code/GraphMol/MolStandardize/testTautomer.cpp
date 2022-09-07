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
#include <chrono>
#include <ctime>

using namespace RDKit;
using namespace MolStandardize;

void testEnumerator() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing tautomer enumeration" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string tautomerFile =
      rdbase + "/Code/GraphMol/MolStandardize/test_data/tautomerTransforms.in";
  auto tautparams = std::unique_ptr<TautomerCatalogParams>(
      new TautomerCatalogParams(tautomerFile));

  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 36);

  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));

  std::function<void(const std::string &, const std::vector<std::string> &)>
      checkAns([te](const std::string &smi,
                    const std::vector<std::string> &ans) {
        ROMOL_SPTR m(SmilesToMol(smi));
        TautomerEnumeratorResult res = te.enumerate(*m);
        TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
        size_t sz = std::max(res.size(), ans.size());
        bool exceedingTautomer = false;
        bool missingTautomer = false;
        bool wrongTautomer = false;
        for (size_t i = 0; i < sz; ++i) {
          if (i >= res.size()) {
            missingTautomer = true;
            std::cerr << "missingTautomer, ans " << ans[i] << std::endl;
          } else if (i >= ans.size()) {
            exceedingTautomer = true;
            std::cerr << "exceedingTautomer, taut " << MolToSmiles(*res.at(i))
                      << std::endl;
          } else if (MolToSmiles(*res[i]) != ans[i]) {
            wrongTautomer = true;
            std::cerr << "wrongTautomer, taut " << MolToSmiles(*res[i])
                      << ", ans " << ans[i] << std::endl;
          }
        }
        TEST_ASSERT(!(missingTautomer || exceedingTautomer || wrongTautomer));
      });

  // Enumerate 1,3 keto/enol tautomer.
  checkAns("C1(=CCCCC1)O", {"O=C1CCCCC1", "OC1=CCCCC1"});

  // Enumerate 1,3 keto/enol tautomer.
  checkAns("C1(=CCCCC1)O", {"O=C1CCCCC1", "OC1=CCCCC1"});

  // Enumerate 1,3 keto/enol tautomer.
  checkAns("C1(CCCCC1)=O", {"O=C1CCCCC1", "OC1=CCCCC1"});

  // Enumerate acetophenone keto/enol tautomer.
  checkAns("C(=C)(O)C1=CC=CC=C1", {"C=C(O)c1ccccc1", "CC(=O)c1ccccc1"});

  // Enumerate acetone keto/enol tautomer.
  checkAns("CC(C)=O", {"C=C(C)O", "CC(C)=O"});

  // keto/enol tautomer
  checkAns("OC(C)=C(C)C", {"C=C(O)C(C)C", "CC(=O)C(C)C", "CC(C)=C(C)O"});

  // 1-phenyl-2-propanone enol/keto
  checkAns("c1(ccccc1)CC(=O)C",
           {"C=C(O)Cc1ccccc1", "CC(=O)Cc1ccccc1", "CC(O)=Cc1ccccc1"});

  // 1,5 keto/enol tautomer
  checkAns("Oc1nccc2cc[nH]c(=N)c12",
           {"N=c1[nH]ccc2cc[nH]c(=O)c12", "N=c1[nH]ccc2ccnc(O)c12",
            "N=c1nccc2cc[nH]c(O)c1-2", "Nc1[nH]ccc2ccnc(=O)c1-2",
            "Nc1nccc2cc[nH]c(=O)c12", "Nc1nccc2ccnc(O)c12"});

  // 1,5 keto/enol tautomer
  checkAns("C1(C=CCCC1)=O", {"O=C1C=CCCC1", "O=C1CC=CCC1", "OC1=CC=CCC1",
                             "OC1=CCC=CC1", "OC1=CCCC=C1"});

  // 1,5 keto/enol tautomer
  checkAns("C1(=CC=CCC1)O", {"O=C1C=CCCC1", "O=C1CC=CCC1", "OC1=CC=CCC1",
                             "OC1=CCC=CC1", "OC1=CCCC=C1"});

  // aliphatic imine tautomer
  checkAns("C1(CCCCC1)=N", {"N=C1CCCCC1", "NC1=CCCCC1"});

  // aliphatic imine tautomer
  checkAns("C1(=CCCCC1)N", {"N=C1CCCCC1", "NC1=CCCCC1"});

  // special imine tautomer
  checkAns("C1(C=CC=CN1)=CC", {"CC=C1C=CC=CN1", "CC=C1C=CCC=N1", "CCc1ccccn1"});

  // special imine tautomer
  checkAns("C1(=NC=CC=C1)CC", {"CC=C1C=CC=CN1", "CC=C1C=CCC=N1", "CCc1ccccn1"});

  // 1,3 aromatic heteroatom H shift
  checkAns("O=c1cccc[nH]1", {"O=c1cccc[nH]1", "Oc1ccccn1"});

  // 1,3 aromatic heteroatom H shift
  checkAns("Oc1ccccn1", {"O=c1cccc[nH]1", "Oc1ccccn1"});

  // 1,3 aromatic heteroatom H shift
  checkAns("Oc1ncc[nH]1", {"O=c1[nH]cc[nH]1", "Oc1ncc[nH]1"});

  // 1,3 heteroatom H shift
  checkAns("OC(C)=NC", {"C=C(O)NC", "CN=C(C)O", "CNC(C)=O"});

  // 1,3 heteroatom H shift
  checkAns("CNC(C)=O", {"C=C(O)NC", "CN=C(C)O", "CNC(C)=O"});

  // 1,3 heteroatom H shift
  checkAns("S=C(N)N", {"N=C(N)S", "NC(N)=S"});

  // 1,3 heteroatom H shift
  checkAns("SC(N)=N", {"N=C(N)S", "NC(N)=S"});

  // 1,3 heteroatom H shift
  checkAns("N=c1[nH]ccn(C)1", {"Cn1cc[nH]c1=N", "Cn1ccnc1N"});

  // 1,3 heteroatom H shift
  checkAns("CN=c1[nH]cncc1",
           {"CN=c1cc[nH]cn1", "CN=c1ccnc[nH]1", "CNc1ccncn1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("Oc1cccc2ccncc12", {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"});

  // 1,5 aromatic heteroatom H shift
  checkAns("O=c1cccc2cc[nH]cc1-2", {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"});

  // 1,5 aromatic heteroatom H shift
  checkAns("Cc1n[nH]c2ncnn12",
           {"C=C1NN=C2N=CNN12", "C=C1NN=C2NC=NN12", "C=C1NNc2ncnn21",
            "Cc1n[nH]c2ncnn12", "Cc1nnc2[nH]cnn12", "Cc1nnc2nc[nH]n12"});

  // 1,5 aromatic heteroatom H shift
  checkAns("Cc1nnc2nc[nH]n12",
           {"C=C1NN=C2N=CNN12", "C=C1NN=C2NC=NN12", "C=C1NNc2ncnn21",
            "Cc1n[nH]c2ncnn12", "Cc1nnc2[nH]cnn12", "Cc1nnc2nc[nH]n12"});

  // 1,5 aromatic heteroatom H shift
  checkAns("Oc1ccncc1", {"O=c1cc[nH]cc1", "Oc1ccncc1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("Oc1c(cccc3)c3nc2ccncc12",
           {"O=c1c2c[nH]ccc-2nc2ccccc12", "O=c1c2ccccc2[nH]c2ccncc12",
            "Oc1c2ccccc2nc2ccncc12"});

  // 1,3 and 1,5 aromatic heteroatom H shift
  checkAns("Oc1ncncc1", {"O=c1cc[nH]cn1", "O=c1ccnc[nH]1", "Oc1ccncn1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("C2(=C1C(=NC=N1)[NH]C(=N2)N)O",
           {"N=c1[nH]c(=O)c2[nH]cnc2[nH]1", "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
            "N=c1[nH]c2ncnc-2c(O)[nH]1", "N=c1nc(O)c2[nH]cnc2[nH]1",
            "N=c1nc(O)c2nc[nH]c2[nH]1", "N=c1nc2[nH]cnc2c(O)[nH]1",
            "N=c1nc2nc[nH]c2c(O)[nH]1", "Nc1nc(=O)c2[nH]cnc2[nH]1",
            "Nc1nc(=O)c2nc[nH]c2[nH]1", "Nc1nc(O)c2[nH]cnc2n1",
            "Nc1nc(O)c2nc[nH]c2n1", "Nc1nc(O)c2ncnc-2[nH]1",
            "Nc1nc2[nH]cnc2c(=O)[nH]1", "Nc1nc2nc[nH]c2c(=O)[nH]1",
            "Nc1nc2ncnc-2c(O)[nH]1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O",
           {"N=c1[nH]c(=O)c2[nH]cnc2[nH]1", "N=c1[nH]c(=O)c2nc[nH]c2[nH]1",
            "N=c1[nH]c2ncnc-2c(O)[nH]1", "N=c1nc(O)c2[nH]cnc2[nH]1",
            "N=c1nc(O)c2nc[nH]c2[nH]1", "N=c1nc2[nH]cnc2c(O)[nH]1",
            "N=c1nc2nc[nH]c2c(O)[nH]1", "Nc1nc(=O)c2[nH]cnc2[nH]1",
            "Nc1nc(=O)c2nc[nH]c2[nH]1", "Nc1nc(O)c2[nH]cnc2n1",
            "Nc1nc(O)c2nc[nH]c2n1", "Nc1nc(O)c2ncnc-2[nH]1",
            "Nc1nc2[nH]cnc2c(=O)[nH]1", "Nc1nc2nc[nH]c2c(=O)[nH]1",
            "Nc1nc2ncnc-2c(O)[nH]1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("Oc1n(C)ncc1", {"CN1N=CCC1=O", "Cn1[nH]ccc1=O", "Cn1nccc1O"});

  // 1,5 aromatic heteroatom H shift
  checkAns("O=c1nc2[nH]ccn2cc1",
           {"O=c1ccn2cc[nH]c2n1", "O=c1ccn2ccnc2[nH]1", "Oc1ccn2ccnc2n1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("N=c1nc[nH]cc1", {"N=c1cc[nH]cn1", "N=c1ccnc[nH]1", "Nc1ccncn1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("N=c(c1)ccn2cc[nH]c12", {"N=c1ccn2cc[nH]c2c1", "Nc1ccn2ccnc2c1"});

  // 1,5 aromatic heteroatom H shift
  checkAns("CN=c1nc[nH]cc1",
           {"CN=c1cc[nH]cn1", "CN=c1ccnc[nH]1", "CNc1ccncn1"});

  // 1,7 aromatic heteroatom H shift
  checkAns("c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
           {"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
            "c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2",
            "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2"});

  // 1,7 aromatic heteroatom H shift
  checkAns("c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2",
           {"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
            "c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2",
            "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2"});

  // 1,9 aromatic heteroatom H shift
  checkAns("CNc1ccnc2ncnn21", {"CN=c1cc[nH]c2ncnn12", "CN=c1ccnc2[nH]cnn12",
                               "CN=c1ccnc2nc[nH]n12", "CNc1ccnc2ncnn12"});

  // 1,9 aromatic heteroatom H shift
  checkAns("CN=c1ccnc2nc[nH]n21", {"CN=c1cc[nH]c2ncnn12", "CN=c1ccnc2[nH]cnn12",
                                   "CN=c1ccnc2nc[nH]n12", "CNc1ccnc2ncnn12"});

  // 1,11 aromatic heteroatom H shift
  checkAns("Nc1ccc(C=C2C=CC(=O)C=C2)cc1",
           {"N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1", "N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
            "N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1", "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"});

  // 1,11 aromatic heteroatom H shift
  checkAns("N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
           {"N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1", "N=C1C=CC(=Cc2ccc(O)cc2)C=C1",
            "N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1", "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"});

  // heterocyclic tautomer
  checkAns("n1ccc2ccc[nH]c12", {"c1c[nH]c2nccc-2c1", "c1cnc2[nH]ccc2c1"});

  // heterocyclic tautomer
  checkAns("c1cc(=O)[nH]c2nccn12",
           {"O=c1ccn2cc[nH]c2n1", "O=c1ccn2ccnc2[nH]1", "Oc1ccn2ccnc2n1"});

  // heterocyclic tautomer
  checkAns("c1cnc2c[nH]ccc12", {"c1cc2cc[nH]c2cn1", "c1cc2cc[nH]cc-2n1"});

  // heterocyclic tautomer
  checkAns("n1ccc2c[nH]ccc12", {"c1cc2[nH]ccc2cn1", "c1cc2c[nH]ccc-2n1"});

  // heterocyclic tautomer
  checkAns("c1cnc2ccc[nH]c12", {"c1c[nH]c2ccnc-2c1", "c1cnc2cc[nH]c2c1"});

  // furanone tautomer
  checkAns("C1=CC=C(O1)O", {"O=C1CC=CO1", "Oc1ccco1"});

  // furanone tautomer
  checkAns("O=C1CC=CO1", {"O=C1CC=CO1", "Oc1ccco1"});

  // keten/ynol tautomer
  checkAns("CC=C=O", {"CC#CO", "CC=C=O"});

  // keten/ynol tautomer
  checkAns("CC#CO", {"CC#CO", "CC=C=O"});

  // ionic nitro/aci-nitro tautomer
  checkAns("C([N+](=O)[O-])C", {"CC=[N+]([O-])O", "CC[N+](=O)[O-]"});

  // ionic nitro/aci-nitro tautomer
  checkAns("C(=[N+](O)[O-])C", {"CC=[N+]([O-])O", "CC[N+](=O)[O-]"});

  // oxim nitroso tautomer
  checkAns("CC(C)=NO", {"C=C(C)NO", "CC(C)=NO", "CC(C)N=O"});

  // oxim nitroso tautomer
  checkAns("CC(C)N=O", {"C=C(C)NO", "CC(C)=NO", "CC(C)N=O"});

  // oxim/nitroso tautomer via phenol
  checkAns("O=Nc1ccc(O)cc1",
           {"O=C1C=CC(=NO)C=C1", "O=NC1C=CC(=O)C=C1", "O=Nc1ccc(O)cc1"});

  // oxim/nitroso tautomer via phenol
  checkAns("O=C1C=CC(=NO)C=C1",
           {"O=C1C=CC(=NO)C=C1", "O=NC1C=CC(=O)C=C1", "O=Nc1ccc(O)cc1"});

  // cyano/iso-cyanic acid tautomer
  checkAns("C(#N)O", {"N#CO", "N=C=O"});

  // cyano/iso-cyanic acid tautomer
  checkAns("C(=N)=O", {"N#CO", "N=C=O"});

  // formamidinesulfinic acid tautomer
  checkAns("NC(N)=S(=O)=O",
           {"N=C(N)S(=O)O", "N=C(N)[SH](=O)=O", "NC(N)=S(=O)=O"});

  // formamidinesulfinic acid tautomer
  checkAns("NC(=N)S(=O)O",
           {"N=C(N)S(=O)O", "N=C(N)[SH](=O)=O", "NC(N)=S(=O)=O"});

  // formamidinesulfonic acid tautomer
  checkAns("NC(=N)S(=O)(=O)O", {"N=C(N)S(=O)(=O)O"});

  // isocyanide tautomer
  checkAns("C#N", {"C#N", "[C-]#[NH+]"});

  // isocyanide tautomer
  checkAns("[C-]#[NH+]", {"C#N", "[C-]#[NH+]"});

  // phosphonic acid tautomer
  checkAns("[PH](=O)(O)(O)", {"O=[PH](O)O", "OP(O)O"});

  // phosphonic acid tautomer
  checkAns("P(O)(O)O", {"O=[PH](O)O", "OP(O)O"});

  // Remove stereochemistry from mobile double bonds
  checkAns("c1(ccccc1)/C=C(/O)\\C",
           {"C=C(O)Cc1ccccc1", "CC(=O)Cc1ccccc1", "CC(O)=Cc1ccccc1"});

  // Remove stereochemistry from mobile double bonds
  checkAns("C/C=C/C(C)=O", {"C=C(O)C=CC", "C=CC=C(C)O", "C=CCC(=C)O",
                            "C=CCC(C)=O", "CC=CC(C)=O"});

  // Remove stereochemistry from mobile double bonds
  std::string smi66 = "C/C=C\\C(C)=O";
  ROMOL_SPTR m66(SmilesToMol(smi66));
  TautomerEnumeratorResult res66 = te.enumerate(*m66);
  std::vector<std::string> ans66 = {"C=C(O)C=CC", "C=CC=C(C)O", "C=CCC(=C)O",
                                    "C=CCC(C)=O", "CC=CC(C)=O"};
  TEST_ASSERT(res66.size() == ans66.size());
  TEST_ASSERT(res66.status() == TautomerEnumeratorStatus::Completed);

  std::vector<std::string> sm66;
  for (const auto &r : res66) {
    sm66.push_back(MolToSmiles(*r));
  }
  // sort both for alphabetical order
  std::sort(sm66.begin(), sm66.end());
  std::sort(ans66.begin(), ans66.end());
  TEST_ASSERT(sm66 == ans66);

  // Gaunine tautomers
  std::string smi67 = "N1C(N)=NC=2N=CNC2C1=O";
  ROMOL_SPTR m67(SmilesToMol(smi67));
  TautomerEnumeratorResult res67 = te.enumerate(*m67);
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
  TEST_ASSERT(res67.status() == TautomerEnumeratorStatus::Completed);
  std::vector<std::string> sm67;
  for (const auto &r : res67) {
    sm67.push_back(MolToSmiles(*r));
  }
  // sort both by alphabetical order
  std::sort(sm67.begin(), sm67.end());
  std::sort(ans67.begin(), ans67.end());
  TEST_ASSERT(sm67 == ans67);

  // Test a structure with hundreds of tautomers.
  std::string smi68 = "[H][C](CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O";
  ROMOL_SPTR m68(SmilesToMol(smi68));
  TautomerEnumeratorResult res68 = te.enumerate(*m68);
  // the maxTransforms limit is hit before the maxTautomers one
  TEST_ASSERT(res68.size() == 292);
  TEST_ASSERT(res68.status() == TautomerEnumeratorStatus::MaxTransformsReached);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testEnumeratorParams() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing TautomerEnumerator params"
      << std::endl;

  // Test a structure with hundreds of tautomers.
  std::string smi68 = "C(CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O";
  ROMOL_SPTR m68(SmilesToMol(smi68));

  {
    TautomerEnumerator te;
    TautomerEnumeratorResult res68 = te.enumerate(*m68);
    TEST_ASSERT(res68.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res68.size() == 72);
  }
  {  // test v1 of the tautomerization parameters
    std::unique_ptr<TautomerEnumerator> te(getV1TautomerEnumerator());
    TautomerEnumeratorResult res68 = te->enumerate(*m68);
    TEST_ASSERT(res68.status() ==
                TautomerEnumeratorStatus::MaxTransformsReached);
    TEST_ASSERT(res68.size() == 292);
  }
  {
    CleanupParameters params;
    params.maxTautomers = 50;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res68 = te.enumerate(*m68);
    TEST_ASSERT(res68.size() == 50);
    TEST_ASSERT(res68.status() ==
                TautomerEnumeratorStatus::MaxTautomersReached);
  }
  std::string sAlaSmi = "C[C@H](N)C(=O)O";
  ROMOL_SPTR sAla(SmilesToMol(sAlaSmi));
  {
    // test remove (S)-Ala stereochemistry
    TEST_ASSERT(sAla->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(sAla->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = true;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*sAla);
    for (const auto &taut : res) {
      TEST_ASSERT(taut->getAtomWithIdx(1)->getChiralTag() ==
                  Atom::CHI_UNSPECIFIED);
      TEST_ASSERT(
          !taut->getAtomWithIdx(1)->hasProp(common_properties::_CIPCode));
    }
  }
  {
    // test retain (S)-Ala stereochemistry
    TEST_ASSERT(sAla->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(sAla->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*sAla);
    for (const auto &taut : res) {
      const auto tautAtom = taut->getAtomWithIdx(1);
      if (tautAtom->getHybridization() == Atom::SP3) {
        TEST_ASSERT(tautAtom->hasProp(common_properties::_CIPCode));
        TEST_ASSERT(
            tautAtom->getProp<std::string>(common_properties::_CIPCode) == "S");
        TEST_ASSERT(tautAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
      } else {
        TEST_ASSERT(!tautAtom->hasProp(common_properties::_CIPCode));
        TEST_ASSERT(tautAtom->getChiralTag() == Atom::CHI_UNSPECIFIED);
      }
    }
  }
  std::string eEnolSmi = "C/C=C/O";
  ROMOL_SPTR eEnol(SmilesToMol(eEnolSmi));
  TEST_ASSERT(eEnol->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
  {
    // test remove enol E stereochemistry
    CleanupParameters params;
    params.tautomerRemoveBondStereo = true;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*eEnol);
    for (const auto &taut : res) {
      TEST_ASSERT(taut->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);
    }
  }
  {
    // test retain enol E stereochemistry
    CleanupParameters params;
    params.tautomerRemoveBondStereo = false;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*eEnol);
    for (const auto &taut : res) {
      if (taut->getBondWithIdx(1)->getBondType() == Bond::DOUBLE) {
        TEST_ASSERT(taut->getBondWithIdx(1)->getStereo() == Bond::STEREOE);
      }
    }
  }
  ROMOL_SPTR zEnol = "C/C=C\\O"_smiles;
  TEST_ASSERT(zEnol->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
  {
    // test remove enol Z stereochemistry
    CleanupParameters params;
    params.tautomerRemoveBondStereo = true;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*zEnol);
    for (const auto &taut : res) {
      TEST_ASSERT(taut->getBondWithIdx(1)->getStereo() == Bond::STEREONONE);
    }
  }
  {
    // test retain enol Z stereochemistry
    CleanupParameters params;
    params.tautomerRemoveBondStereo = false;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*zEnol);
    for (const auto &taut : res) {
      if (taut->getBondWithIdx(1)->getBondType() == Bond::DOUBLE) {
        TEST_ASSERT(taut->getBondWithIdx(1)->getStereo() == Bond::STEREOZ);
      }
    }
  }
  ROMOL_SPTR chembl2024142 =
      "[2H]C1=C(C(=C2C(=C1[2H])C(=O)C(=C(C2=O)C([2H])([2H])[2H])C/C=C(\\C)/CC([2H])([2H])/C=C(/CC/C=C(\\C)/CCC=C(C)C)\\C([2H])([2H])[2H])[2H])[2H]"_smiles;
  MolOps::RemoveHsParameters hparams;
  hparams.removeAndTrackIsotopes = true;
  chembl2024142.reset(MolOps::removeHs(*chembl2024142, hparams));
  TEST_ASSERT(chembl2024142->getAtomWithIdx(12)->hasProp(
      common_properties::_isotopicHs));
  {
    // test remove isotopic Hs involved in tautomerism
    CleanupParameters params;
    params.tautomerRemoveIsotopicHs = true;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*chembl2024142);
    for (const auto &taut : res) {
      const auto tautAtom = taut->getAtomWithIdx(12);
      TEST_ASSERT(!tautAtom->hasProp(common_properties::_isotopicHs));
    }
  }
  {
    // test retain isotopic Hs involved in tautomerism
    CleanupParameters params;
    params.tautomerRemoveIsotopicHs = false;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*chembl2024142);
    for (const auto &taut : res) {
      const auto tautAtom = taut->getAtomWithIdx(12);
      TEST_ASSERT(tautAtom->hasProp(common_properties::_isotopicHs));
    }
  }
  ROMOL_SPTR enolexample = "[2H]OC=C"_smiles;
  enolexample.reset(MolOps::removeHs(*enolexample, hparams));
  TEST_ASSERT(
      enolexample->getAtomWithIdx(0)->hasProp(common_properties::_isotopicHs));
  {
    CleanupParameters params;
    params.tautomerRemoveIsotopicHs = true;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*enolexample);
    for (const auto &taut : res) {
      const auto tautAtom = taut->getAtomWithIdx(0);
      TEST_ASSERT(!(tautAtom->hasProp(common_properties::_isotopicHs) &&
                    !tautAtom->getTotalNumHs()));
    }
  }
  {
    CleanupParameters params;
    params.tautomerRemoveIsotopicHs = false;
    TautomerEnumerator te(params);
    TautomerEnumeratorResult res = te.enumerate(*enolexample);
    for (const auto &taut : res) {
      const auto tautAtom = taut->getAtomWithIdx(0);
      TEST_ASSERT(!(tautAtom->hasProp(common_properties::_isotopicHs) &&
                    !tautAtom->getTotalNumHs()));
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testEnumeratorCallback() {
  class MyTautomerEnumeratorCallback : public TautomerEnumeratorCallback {
   public:
    MyTautomerEnumeratorCallback(double timeoutMs)
        : d_timeoutMs(timeoutMs), d_start(std::chrono::system_clock::now()) {}
    bool operator()(const ROMol &, const TautomerEnumeratorResult &) override {
      double elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
                             std::chrono::system_clock::now() - d_start)
                             .count();
      return (elapsedMs < d_timeoutMs);
    }

   private:
    double d_timeoutMs;
    std::chrono::time_point<std::chrono::system_clock> d_start;
  };

  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing TautomerEnumerator callback"
      << std::endl;

  // Test a structure with hundreds of tautomers.
  std::string smi68 = "[H][C](CO)(NC(=O)C1=C(O)C(O)=CC=C1)C(O)=O";
  ROMOL_SPTR m68(SmilesToMol(smi68));

  CleanupParameters params;
  params.maxTransforms = 10000;
  params.maxTautomers = 10000;
  params.tautomerTransformData =
      MolStandardize::defaults::defaultTautomerTransformsv1;
  {
    TautomerEnumerator te(params);
    te.setCallback(new MyTautomerEnumeratorCallback(50.0));
    TautomerEnumeratorResult res68 = te.enumerate(*m68);
    // either the enumeration was canceled due to timeout
    // or it has completed very quickly
    bool hasReachedTimeout =
        (res68.size() < 375 &&
         res68.status() == TautomerEnumeratorStatus::Canceled);
    bool hasCompleted = (res68.size() == 375 &&
                         res68.status() == TautomerEnumeratorStatus::Completed);
    if (hasReachedTimeout) {
      std::cerr << "Enumeration was canceled due to timeout (50 ms)"
                << std::endl;
    }
    if (hasCompleted) {
      std::cerr << "Enumeration has completed" << std::endl;
    }
    TEST_ASSERT(hasReachedTimeout || hasCompleted);
    TEST_ASSERT(hasReachedTimeout ^ hasCompleted);
  }
  {
    TautomerEnumerator te(params);
    te.setCallback(new MyTautomerEnumeratorCallback(10000.0));
    TautomerEnumeratorResult res68 = te.enumerate(*m68);
    std::cerr << res68.size() << std::endl;
    // either the enumeration completed
    // or it ran very slowly and was canceled due to timeout
    bool hasReachedTimeout =
        (res68.size() < 375 &&
         res68.status() == TautomerEnumeratorStatus::Canceled);
    bool hasCompleted = (res68.size() == 375 &&
                         res68.status() == TautomerEnumeratorStatus::Completed);
    if (hasReachedTimeout) {
      std::cerr << "Enumeration was canceled due to timeout (10 s)"
                << std::endl;
    }
    if (hasCompleted) {
      std::cerr << "Enumeration has completed" << std::endl;
    }
    TEST_ASSERT(hasReachedTimeout || hasCompleted);
    TEST_ASSERT(hasReachedTimeout ^ hasCompleted);
  }
  {
    // GitHub #4736
    TautomerEnumerator te(params);
    te.setCallback(new MyTautomerEnumeratorCallback(50.0));
    TautomerEnumeratorResult res68 = te.enumerate(*m68);
    TautomerEnumerator teCopy(te);
    TautomerEnumeratorResult res68Copy = teCopy.enumerate(*m68);
    TEST_ASSERT(res68.status() == res68Copy.status());
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

// tests from the molvs repo:
// https://github.com/mcs07/MolVS/blob/456f2fe723acfedbf634a8bcfe943b83ea7d4c20/tests/test_tautomer.py
std::vector<std::pair<std::string, std::string>> canonTautomerData{
    {"C1(=CCCCC1)O", "O=C1CCCCC1"},
    {"C1(CCCCC1)=O", "O=C1CCCCC1"},
    {"C(=C)(O)C1=CC=CC=C1", "CC(=O)c1ccccc1"},
    {"CC(C)=O", "CC(C)=O"},
    {"OC(C)=C(C)C", "CC(=O)C(C)C"},
    {"c1(ccccc1)CC(=O)C", "CC(=O)Cc1ccccc1"},
    {"Oc1nccc2cc[nH]c(=N)c12", "Nc1nccc2cc[nH]c(=O)c12"},
    {"C1(C=CCCC1)=O", "O=C1C=CCCC1"},
    {"C1(CCCCC1)=N", "N=C1CCCCC1"},
    {"C1(=CCCCC1)N", "N=C1CCCCC1"},
    {"C1(C=CC=CN1)=CC", "CCc1ccccn1"},
    {"C1(=NC=CC=C1)CC", "CCc1ccccn1"},
    {"O=c1cccc[nH]1", "O=c1cccc[nH]1"},
    {"Oc1ccccn1", "O=c1cccc[nH]1"},
    {"Oc1ncc[nH]1", "O=c1[nH]cc[nH]1"},
    {"OC(C)=NC", "CNC(C)=O"},
    {"CNC(C)=O", "CNC(C)=O"},
    {"S=C(N)N", "NC(N)=S"},
    {"SC(N)=N", "NC(N)=S"},
    {"N=c1[nH]ccn(C)1", "Cn1ccnc1N"},
    {"CN=c1[nH]cncc1", "CNc1ccncn1"},
    {"Oc1cccc2ccncc12", "Oc1cccc2ccncc12"},
    {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"},
    {"Cc1n[nH]c2ncnn12", "Cc1n[nH]c2ncnn12"},
    {"Cc1nnc2nc[nH]n12", "Cc1n[nH]c2ncnn12"},
    {"Oc1ccncc1", "O=c1cc[nH]cc1"},
    {"Oc1c(cccc3)c3nc2ccncc12", "O=c1c2ccccc2[nH]c2ccncc12"},
    {"Oc1ncncc1", "O=c1cc[nH]cn1"},
    {"C2(=C1C(=NC=N1)[NH]C(=N2)N)O", "Nc1nc(=O)c2[nH]cnc2[nH]1"},
    {"C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O", "Nc1nc(=O)c2[nH]cnc2[nH]1"},
    {"Oc1n(C)ncc1", "Cn1[nH]ccc1=O"},
    {"O=c1nc2[nH]ccn2cc1", "O=c1ccn2cc[nH]c2n1"},
    {"N=c1nc[nH]cc1", "Nc1ccncn1"},
    {"N=c(c1)ccn2cc[nH]c12", "Nc1ccn2ccnc2c1"},
    {"CN=c1nc[nH]cc1", "CNc1ccncn1"},
    {"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
     "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1"},
    {"c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2",
     "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1"},
    {"CNc1ccnc2ncnn21", "CNc1ccnc2ncnn12"},
    {"CN=c1ccnc2nc[nH]n21", "CNc1ccnc2ncnn12"},
    {"Nc1ccc(C=C2C=CC(=O)C=C2)cc1", "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"},
    {"N=C1C=CC(=Cc2ccc(O)cc2)C=C1", "Nc1ccc(C=C2C=CC(=O)C=C2)cc1"},
    {"n1ccc2ccc[nH]c12", "c1cnc2[nH]ccc2c1"},
    {"c1cc(=O)[nH]c2nccn12", "O=c1ccn2cc[nH]c2n1"},
    {"c1cnc2c[nH]ccc12", "c1cc2cc[nH]c2cn1"},
    {"n1ccc2c[nH]ccc12", "c1cc2[nH]ccc2cn1"},
    {"c1cnc2ccc[nH]c12", "c1cnc2cc[nH]c2c1"},
    {"C1=CC=C(O1)O", "Oc1ccco1"},
    {"O=C1CC=CO1", "Oc1ccco1"},
    {"CC=C=O", "CC=C=O"},
    {"CC#CO", "CC=C=O"},
    {"C([N+](=O)[O-])C", "CC[N+](=O)[O-]"},
    {"C(=[N+](O)[O-])C", "CC[N+](=O)[O-]"},
    {"CC(C)=NO", "CC(C)=NO"},
    {"CC(C)N=O", "CC(C)=NO"},
    {"O=Nc1ccc(O)cc1", "O=Nc1ccc(O)cc1"},
    {"O=C1C=CC(=NO)C=C1", "O=Nc1ccc(O)cc1"},
    {"C(#N)O", "N=C=O"},
    {"C(=N)=O", "N=C=O"},
    {"N=C(N)S(=O)O", "N=C(N)S(=O)O"},
    {"C#N", "C#N"},
    {"[C-]#[NH+]", "C#N"},
    {"[PH](=O)(O)(O)", "O=[PH](O)O"},
    {"P(O)(O)O", "O=[PH](O)O"}};

void testCanonicalize() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing tautomer canonicalization"
      << std::endl;

  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));

  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);

  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));

  for (const auto &itm : canonTautomerData) {
    std::unique_ptr<ROMol> mol{SmilesToMol(itm.first)};
    TEST_ASSERT(mol);
    std::unique_ptr<ROMol> res{te.canonicalize(*mol)};
    TEST_ASSERT(res);
    TEST_ASSERT(MolToSmiles(*res) == itm.second);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testPickCanonical() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing pickCanonical"
                       << std::endl;

  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));

  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);

  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));

  for (const auto &itm : canonTautomerData) {
    std::unique_ptr<ROMol> mol{SmilesToMol(itm.first)};
    TEST_ASSERT(mol);
    auto tautRes = te.enumerate(*mol);
    std::unique_ptr<ROMol> res{te.pickCanonical(tautRes)};
    TEST_ASSERT(res);
    // std::cerr << itm.first<<" -> "<<MolToSmiles(*res)<<std::endl;
    TEST_ASSERT(MolToSmiles(*res) == itm.second);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testCustomScoreFunc() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing custom scoring functions"
      << std::endl;

  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));

  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);

  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));

  // silly examples just using the scoreRings() function
  std::vector<std::pair<std::string, std::string>> subsetTautomerData{
      {"C1(=CCCCC1)O", "O=C1CCCCC1"},
      {"C1(CCCCC1)=O", "O=C1CCCCC1"},
      {"C(=C)(O)C1=CC=CC=C1", "C=C(O)c1ccccc1"},
      {"CC(C)=O", "C=C(C)O"},
      {"OC(C)=C(C)C", "C=C(O)C(C)C"},
  };
  for (const auto &itm : subsetTautomerData) {
    std::unique_ptr<ROMol> mol{SmilesToMol(itm.first)};
    TEST_ASSERT(mol);
    {
      // this uses the non-templated pickCanonical() function
      std::unique_ptr<ROMol> res{
          te.canonicalize(*mol, [](const ROMol &m) -> int {
            return MolStandardize::TautomerScoringFunctions::scoreRings(m);
          })};
      TEST_ASSERT(res);
      TEST_ASSERT(MolToSmiles(*res) == itm.second);
    }
    {
      // this uses the non-templated pickCanonical() overload
      auto tautRes = te.enumerate(*mol);
      std::unique_ptr<ROMol> res{
          te.pickCanonical(tautRes, [](const ROMol &m) -> int {
            return MolStandardize::TautomerScoringFunctions::scoreRings(m);
          })};
      TEST_ASSERT(res);
      TEST_ASSERT(MolToSmiles(*res) == itm.second);
    }
    {
      // this tests the templated pickCanonical() overload on a std::vector
      auto tautRes = te.enumerate(*mol);
      std::unique_ptr<ROMol> res{
          te.pickCanonical(tautRes.tautomers(), [](const ROMol &m) -> int {
            return MolStandardize::TautomerScoringFunctions::scoreRings(m);
          })};
      TEST_ASSERT(res);
      TEST_ASSERT(MolToSmiles(*res) == itm.second);
    }
    {
      // this tests the templated pickCanonical() overload
      // with a different iterable container
      auto tautRes = te.enumerate(*mol);
      std::set<ROMOL_SPTR> tautomerSet(tautRes.begin(), tautRes.end());
      std::unique_ptr<ROMol> res{
          te.pickCanonical(tautomerSet, [](const ROMol &m) -> int {
            return MolStandardize::TautomerScoringFunctions::scoreRings(m);
          })};
      TEST_ASSERT(res);
      TEST_ASSERT(MolToSmiles(*res) == itm.second);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testEnumerationProblems() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing tautomer enumeration problems"
      << std::endl;

  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));

  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);

  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));
#if 1
  {  // from the discussion of #2908
    auto mol = "O=C(C1=C[NH+]=CC=C1)[O-]"_smiles;
    auto tautRes = te.enumerate(*mol);
    TEST_ASSERT(tautRes.size() == 1);
  }
#endif
  {  // one of the examples from the tautobase paper
    auto m =
        "[S:1]=[c:2]1[nH+:3][c:5]([NH2:9])[nH:8][c:7]2[c:4]1[n:6][nH:10][n:11]2"_smiles;
    TEST_ASSERT(m);

    auto tautRes = te.enumerate(*m);
    // for (auto taut : tauts) {
    //   std::cerr << MolToSmiles(*taut) << std::endl;
    // }
    TEST_ASSERT(tautRes.size() == 12);
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testPickCanonical2() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing pickCanonical"
                       << std::endl;

  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));
  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);

  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));
  {
    auto mol = "CN=c1nc[nH]cc1"_smiles;
    TEST_ASSERT(mol);
    auto tautRes = te.enumerate(*mol);
    for (const auto &taut : tautRes) {
      std::cerr << MolToSmiles(*taut) << std::endl;
    }
    std::unique_ptr<ROMol> canon{te.pickCanonical(tautRes)};
    std::cerr << "res: " << MolToSmiles(*canon) << std::endl;
  }
  {
    auto mol = "CN=c1[nH]cccc1"_smiles;
    TEST_ASSERT(mol);
    auto tautRes = te.enumerate(*mol);
    for (const auto &taut : tautRes) {
      std::cerr << MolToSmiles(*taut) << std::endl;
    }
    std::unique_ptr<ROMol> canon{te.pickCanonical(tautRes)};
    std::cerr << "res: " << MolToSmiles(*canon) << std::endl;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testEnumerateDetails() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing getting details back "
         "from tautomer enumeration"
      << std::endl;
  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));
  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);
  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));
  {
    auto mol = "c1ccccc1CN=c1[nH]cccc1"_smiles;
    TEST_ASSERT(mol);

    auto tautRes = te.enumerate(*mol);
    TEST_ASSERT(tautRes.size() == 2);
    TEST_ASSERT(tautRes.modifiedAtoms().count() == 2);
    TEST_ASSERT(tautRes.modifiedBonds().count() == 7);
    TEST_ASSERT(tautRes.modifiedAtoms().test(7));
    TEST_ASSERT(tautRes.modifiedAtoms().test(9));
    TEST_ASSERT(!tautRes.modifiedBonds().test(0));
    TEST_ASSERT(tautRes.modifiedBonds().test(7));
    TEST_ASSERT(tautRes.modifiedBonds().test(8));
    TEST_ASSERT(tautRes.modifiedBonds().test(14));
  }
  {
    // test the deprecated form
    auto mol = "c1ccccc1CN=c1[nH]cccc1"_smiles;
    TEST_ASSERT(mol);
    boost::dynamic_bitset<> atomsModified(mol->getNumAtoms());
    boost::dynamic_bitset<> bondsModified(mol->getNumBonds());

#if defined(_MSC_VER)
#pragma warning(suppress : 4996)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    auto tauts = te.enumerate(*mol, &atomsModified, &bondsModified);
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    TEST_ASSERT(tauts.size() == 2);
    TEST_ASSERT(atomsModified.count() == 2);
    TEST_ASSERT(bondsModified.count() == 7);
    TEST_ASSERT(atomsModified[7]);
    TEST_ASSERT(atomsModified[9]);
    TEST_ASSERT(!bondsModified[0]);
    TEST_ASSERT(bondsModified[7]);
    TEST_ASSERT(bondsModified[8]);
    TEST_ASSERT(bondsModified[14]);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub2990() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github #2990: "
                          "Tautomer enumeration "
                          "should remove stereo in all tautomers"
                       << std::endl;
  auto tautparams =
      std::unique_ptr<TautomerCatalogParams>(new TautomerCatalogParams(""));

  unsigned int ntransforms = tautparams->getTransforms().size();
  TEST_ASSERT(ntransforms == 37);
  TautomerEnumerator te(new TautomerCatalog(tautparams.get()));
  {
    // atom stereo
    auto mol = "COC(=O)[C@@H](N)CO"_smiles;
    TEST_ASSERT(mol);
    auto res = te.enumerate(*mol);
    for (const auto &taut : res) {
      auto smi = MolToSmiles(*taut);
      // std::cerr << smi << std::endl;
      TEST_ASSERT(smi.find("@H") == std::string::npos);
    }
  }
  {
    // atom stereo, atoms not in the tautomer zone are still ok
    auto mol = "C[C@](Cl)(F)COC(=O)[C@@H](N)CO"_smiles;
    TEST_ASSERT(mol);
    auto res = te.enumerate(*mol);
    for (const auto &taut : res) {
      auto smi = MolToSmiles(*taut);
      // std::cerr << smi << std::endl;
      TEST_ASSERT(smi.find("@H") == std::string::npos);
      TEST_ASSERT(smi.find("@]") != std::string::npos);
    }
  }
  {
    // bond stereo
    auto mol = "C/C=C/C/N=c1/[nH]cccc1"_smiles;
    TEST_ASSERT(mol);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondDir() !=
                Bond::BondDir::NONE);
    TEST_ASSERT(mol->getBondBetweenAtoms(2, 3)->getBondDir() !=
                Bond::BondDir::NONE);
    TEST_ASSERT(mol->getBondBetweenAtoms(3, 4)->getBondDir() !=
                Bond::BondDir::NONE);
    TEST_ASSERT(mol->getBondBetweenAtoms(5, 6)->getBondDir() !=
                Bond::BondDir::NONE);
    TEST_ASSERT(mol->getBondBetweenAtoms(1, 2)->getStereo() >
                Bond::BondStereo::STEREOANY);
    TEST_ASSERT(mol->getBondBetweenAtoms(4, 5)->getStereo() >
                Bond::BondStereo::STEREOANY);

    auto res = te.enumerate(*mol);
    for (const auto &taut : res) {
      TEST_ASSERT(taut->getBondBetweenAtoms(0, 1)->getBondDir() !=
                  Bond::BondDir::NONE);
      TEST_ASSERT(taut->getBondBetweenAtoms(2, 3)->getBondDir() !=
                  Bond::BondDir::NONE);
      TEST_ASSERT(taut->getBondBetweenAtoms(3, 4)->getBondDir() ==
                  Bond::BondDir::NONE);
      TEST_ASSERT(taut->getBondBetweenAtoms(5, 6)->getBondDir() ==
                  Bond::BondDir::NONE);
      TEST_ASSERT(taut->getBondBetweenAtoms(1, 2)->getStereo() >
                  Bond::BondStereo::STEREOANY);
      TEST_ASSERT(taut->getBondBetweenAtoms(4, 5)->getStereo() ==
                  Bond::BondStereo::STEREONONE);
    }
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testPickCanonicalCIPChangeOnChiralCenter() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n testPickCanonicalCIPChangeOnChiralCenter"
      << std::endl;

  struct CanonicalTaut {
    static ROMOL_SPTR get(const TautomerEnumeratorResult &res) {
      std::vector<int> scores;
      scores.reserve(res.size());
      std::transform(res.begin(), res.end(), std::back_inserter(scores),
                     [](const ROMOL_SPTR &m) {
                       return TautomerScoringFunctions::scoreTautomer(*m);
                     });
      std::vector<size_t> indices(res.size());
      std::iota(indices.begin(), indices.end(), 0);
      int bestIdx =
          *std::max_element(indices.begin(), indices.end(),
                            [scores](const size_t &a, const size_t &b) {
                              if (scores.at(a) != scores.at(b)) {
                                return (scores.at(a) < scores.at(b));
                              }
                              return (a < b);
                            });
      TEST_ASSERT(*std::max_element(scores.begin(), scores.end()) ==
                  scores.at(bestIdx));
      return res.at(bestIdx);
    }
  };

  auto mol = "CC\\C=C(/O)[C@@H](C)C(C)=O"_smiles;
  TEST_ASSERT(mol.get());
  TEST_ASSERT(mol->getAtomWithIdx(5)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(mol->getAtomWithIdx(5)->getProp<std::string>(
                  common_properties::_CIPCode) == "R");
  {
    // here the chirality disappears as the chiral center is itself involved in
    // tautomerism
    TautomerEnumerator te;
    ROMOL_SPTR canTaut(te.canonicalize(*mol));
    TEST_ASSERT(canTaut.get());
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(
        !canTaut->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(MolToSmiles(*canTaut) == "CCCC(=O)C(C)C(C)=O");
  }
  {
    // here the chirality stays even if the chiral center is itself involved in
    // tautomerism because of the tautomerRemoveSp3Stereo parameter being set to
    // false
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    TautomerEnumerator te(params);
    ROMOL_SPTR canTaut(te.canonicalize(*mol));
    TEST_ASSERT(canTaut.get());
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    TEST_ASSERT(MolToSmiles(*canTaut) == "CCCC(=O)[C@@H](C)C(C)=O");
  }
  {
    // here the chirality disappears as the chiral center is itself involved in
    // tautomerism; the reassignStereo setting has no influence
    TautomerEnumerator te;
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 8);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(
        !bestTaut->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)C(C)C(C)=O");
  }
  {
    // here the chirality disappears as the chiral center is itself involved in
    // tautomerism; the reassignStereo setting has no influence
    CleanupParameters params;
    params.tautomerReassignStereo = false;
    TautomerEnumerator te(params);
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 8);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_UNSPECIFIED);
    TEST_ASSERT(
        !bestTaut->getAtomWithIdx(5)->hasProp(common_properties::_CIPCode));
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)C(C)C(C)=O");
  }
  {
    // here the chirality stays even if the chiral center is itself involved in
    // tautomerism because of the tautomerRemoveSp3Stereo parameter being set to
    // false. As reassignStereo by default is true, the CIP code has  been
    // recomputed and therefore it is now S (correct)
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    TautomerEnumerator te(params);
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 8);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)[C@@H](C)C(C)=O");
  }
  {
    // here the chirality stays even if the chiral center is itself involved in
    // tautomerism because of the tautomerRemoveSp3Stereo parameter being set to
    // false. As reassignStereo is false, the CIP code has not been recomputed
    // and therefore it is still R (incorrect)
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    params.tautomerReassignStereo = false;
    TautomerEnumerator te(params);
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 8);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "R");
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)[C@@H](C)C(C)=O");
  }

  mol = "CC\\C=C(/O)[C@@](CC)(C)C(C)=O"_smiles;
  TEST_ASSERT(mol.get());
  TEST_ASSERT(mol->getAtomWithIdx(5)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW);
  TEST_ASSERT(mol->getAtomWithIdx(5)->getProp<std::string>(
                  common_properties::_CIPCode) == "S");
  // here the chirality stays no matter how tautomerRemoveSp3Stereo
  // is set as the chiral center is not involved in tautomerism
  {
    TautomerEnumerator te;
    ROMOL_SPTR canTaut(te.canonicalize(*mol));
    TEST_ASSERT(canTaut.get());
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "R");
    TEST_ASSERT(MolToSmiles(*canTaut) == "CCCC(=O)[C@](C)(CC)C(C)=O");
  }
  {
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    TautomerEnumerator te(params);
    ROMOL_SPTR canTaut(te.canonicalize(*mol));
    TEST_ASSERT(canTaut.get());
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(canTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "R");
    TEST_ASSERT(MolToSmiles(*canTaut) == "CCCC(=O)[C@](C)(CC)C(C)=O");
  }
  {
    // as reassignStereo by default is true, the CIP code has been recomputed
    // and therefore it is now R (correct)
    TautomerEnumerator te;
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 4);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "R");
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)[C@](C)(CC)C(C)=O");
  }
  {
    // as reassignStereo is false, the CIP code has not been recomputed
    // and therefore it is still S (incorrect)
    CleanupParameters params;
    params.tautomerReassignStereo = false;
    TautomerEnumerator te(params);
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 4);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)[C@](C)(CC)C(C)=O");
  }
  {
    // as reassignStereo by default is true, the CIP code has  been recomputed
    // and therefore it is now R (correct)
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    TautomerEnumerator te(params);
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 4);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "R");
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)[C@](C)(CC)C(C)=O");
  }
  {
    // here the chirality stays even if the tautomerRemoveSp3Stereo parameter
    // is set to false as the chiral center is not involved in tautomerism.
    // As reassignStereo is false, the CIP code has not been recomputed
    // and therefore it is still S (incorrect)
    CleanupParameters params;
    params.tautomerRemoveSp3Stereo = false;
    params.tautomerReassignStereo = false;
    TautomerEnumerator te(params);
    auto res = te.enumerate(*mol);
    TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
    TEST_ASSERT(res.size() == 4);
    ROMOL_SPTR bestTaut = CanonicalTaut::get(res);
    TEST_ASSERT(bestTaut.get());
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(bestTaut->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::_CIPCode) == "S");
    TEST_ASSERT(MolToSmiles(*bestTaut) == "CCCC(=O)[C@](C)(CC)C(C)=O");
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testTautomerEnumeratorResult_const_iterator() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n testTautomerEnumeratorResult_const_iterator"
      << std::endl;
  // CHEMBL3480964
  RWMOL_SPTR mol = "Cc1nnc(NC(=O)N2CCN(Cc3ccc(F)cc3)C(=O)C2)s1"_smiles;
  TautomerEnumerator te;
  auto res = te.enumerate(*mol);
  TEST_ASSERT(res.status() == TautomerEnumeratorStatus::Completed);
  TEST_ASSERT(res.size() == 6);
  auto it = res.begin();
  auto it2 = res.begin();
  // Test semantic requirements of bidirectional_iterator
  // https://en.cppreference.com/w/cpp/iterator/bidirectional_iterator
  TEST_ASSERT(it == it2);
  TEST_ASSERT(it++ == it2);
  TEST_ASSERT(it == ++it2);
  TEST_ASSERT(it == it2);
  TEST_ASSERT(it-- == it2);
  TEST_ASSERT(it == --it2);
  TEST_ASSERT(it == it2);
  ++it;
  ++it2;
  TEST_ASSERT(++(--it) == it2);
  TEST_ASSERT(--(++it) == it2);
  TEST_ASSERT(std::addressof(--it) == std::addressof(it));
  ++it;
  TEST_ASSERT(it == it2);
  it--;
  --it2;
  TEST_ASSERT(it == it2);
  TEST_ASSERT(*it == res[0]);
  TEST_ASSERT(*it++ == res.at(0));
  TEST_ASSERT(*it == res[1]);
  TEST_ASSERT(*++it == res.at(2));
  TEST_ASSERT(*it == res[2]);
  ++it;
  TEST_ASSERT(*it == res[3]);
  ++it;
  TEST_ASSERT(*it == res[4]);
  it++;
  TEST_ASSERT(*it == res[5]);
  TEST_ASSERT(*it-- == res.at(5));
  TEST_ASSERT(*it == res[4]);
  TEST_ASSERT(*--it == res.at(3));
  TEST_ASSERT(*it == res[3]);
  --it;
  TEST_ASSERT(*it == res[2]);
  --it;
  TEST_ASSERT(*it == res[1]);
  it--;
  TEST_ASSERT(*it == res[0]);
  std::ptrdiff_t i = 0;
  for (auto t : res) {
    TEST_ASSERT(t == res[i++]);
  }
  i = 0;
  for (auto it = res.begin(); it != res.end(); ++it) {
    TEST_ASSERT(std::distance(res.begin(), it) == i);
    TEST_ASSERT(*it == res[i]);
    TEST_ASSERT(it->getNumAtoms() == res[i++]->getNumAtoms());
  }
  i = res.size();
  for (auto it = res.end(); it != res.begin();) {
    TEST_ASSERT(std::distance(res.begin(), it) == i);
    TEST_ASSERT(*--it == res[--i]);
    TEST_ASSERT(it->getNumAtoms() == res[i]->getNumAtoms());
  }
  i = 0;
  for (const auto &pair : res.smilesTautomerMap()) {
    TEST_ASSERT(pair.first == MolToSmiles(*res[i]));
    TEST_ASSERT(pair.second.tautomer == res[i++]);
  }
  i = 0;
  for (auto it = res.smilesTautomerMap().begin();
       it != res.smilesTautomerMap().end(); ++it) {
    TEST_ASSERT(std::distance(res.smilesTautomerMap().begin(), it) == i);
    TEST_ASSERT(it->first == MolToSmiles(*res[i]));
    TEST_ASSERT(it->second.tautomer == res[i++]);
  }
  i = res.smilesTautomerMap().size();
  for (auto it = res.smilesTautomerMap().end();
       it != res.smilesTautomerMap().begin();) {
    TEST_ASSERT(std::distance(res.smilesTautomerMap().begin(), it) == i);
    TEST_ASSERT((--it)->first == MolToSmiles(*res[--i]));
    TEST_ASSERT(it->second.tautomer == res[i]);
  }
}

void testGithub3430() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n testGithub3430"
                       << std::endl;
  // The "guanidine terminal=N" rule should not apply to aromatic C
  // as this balances the "aromatic C = exocyclic N" rule with no net
  // effect on the score
  std::vector<ROMOL_SPTR> mols{"Cc1ccc(NC(=O)N=c2[nH]c(C)cn2C)nc1"_smiles,
                               "CCCCC(=O)N=c1nc(C)c2ncn(C)c2[nH]1"_smiles,
                               "c12ccccc1[nH]c(=N)[nH]2"_smiles};
  for (auto mol : mols) {
    TEST_ASSERT(mol);
    TautomerEnumerator te;
    auto res = te.enumerate(*mol);
    std::vector<int> scores;
    scores.reserve(res.size());
    std::transform(res.begin(), res.end(), std::back_inserter(scores),
                   [](const ROMOL_SPTR &m) {
                     return TautomerScoringFunctions::scoreTautomer(*m);
                   });

    std::sort(scores.begin(), scores.end(), std::greater<int>());
    TEST_ASSERT(scores[1] < scores[0]);
  }
}

void testGithub3755() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n testGithub3755"
                       << std::endl;
  // hydrates, aminals and hemiaminals should be scored lower than
  // carboxylic acids, amides, amidines, and guanidines
  std::vector<std::pair<std::string, std::string>> orig_vs_expected{
      {"OC(=O)C(N)CO", "NC(CO)C(=O)O"}, {"C([C@@H](C(=O)O)N)O", "NC(CO)C(=O)O"},
      {"OC(=O)C(N)CN", "NCC(N)C(=O)O"}, {"NC(=O)C(N)CO", "NC(=O)C(N)CO"},
      {"NC(=N)C(N)CO", "N=C(N)C(N)CO"}, {"NC(=N)NC(N)CO", "N=C(N)NC(N)CO"}};
  TautomerEnumerator te;
  for (const auto &pair : orig_vs_expected) {
    ROMOL_SPTR orig(SmilesToMol(pair.first));
    TEST_ASSERT(orig);
    ROMOL_SPTR canonical(te.canonicalize(*orig));
    TEST_ASSERT(MolToSmiles(*canonical) == pair.second);
  }
}

int main() {
  RDLog::InitLogs();
#if 1
  testEnumerator();
  testEnumeratorParams();
  testEnumeratorCallback();
  testCanonicalize();
  testPickCanonical();
  testCustomScoreFunc();
  testEnumerationProblems();
#endif
  testPickCanonical2();
  testEnumerateDetails();
  testGithub2990();
  testPickCanonicalCIPChangeOnChiralCenter();
  testTautomerEnumeratorResult_const_iterator();
  testGithub3430();
  testGithub3755();
  return 0;
}
