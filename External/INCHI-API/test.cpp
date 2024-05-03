#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include "inchi.h"

using namespace RDKit;

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock(const std::vector<ROMol *> &mols, unsigned int count,
              unsigned int idx, const std::vector<std::string> &inchis,
              const std::vector<std::string> &keys) {
  for (unsigned int j = 0; j < 200; j++) {
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) {
        continue;
      }
      ROMol *mol = mols[i];
      ExtraInchiReturnValues tmp;
      std::string inchi = MolToInchi(*mol, tmp);
      TEST_ASSERT(inchi == inchis[i]);
      std::string key = InchiToInchiKey(inchi);
      TEST_ASSERT(key == keys[i]);
      std::string key2 = MolToInchiKey(*mol);
      TEST_ASSERT(key2 == keys[i]);

      ROMol *mol2 = InchiToMol(inchi, tmp);
      TEST_ASSERT(mol2);
      ExtraInchiReturnValues tmp2;
      std::string inchi2 = MolToInchi(*mol2, tmp2);
      TEST_ASSERT(inchi == inchi2);
      delete mol2;

      std::string mol_block = MolToMolBlock(*mol);
      ExtraInchiReturnValues tmp3;
      std::string inchi3 = MolBlockToInchi(mol_block, tmp3);
      TEST_ASSERT(inchi == inchi3);
    }
  }
};
}  // namespace

#include <thread>
#include <future>
void testMultiThread() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  std::cerr << "reading molecules" << std::endl;
  std::vector<ROMol *> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = suppl.next();
    TEST_ASSERT(mol != nullptr);
    mols.push_back(mol);
  }
  TEST_ASSERT(mols.size() == 100);

  std::cerr << "generating reference data" << std::endl;
  std::vector<std::string> inchis;
  std::vector<std::string> keys;
  for (const auto *mol : mols) {
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*mol, tmp);
    std::string key = InchiToInchiKey(inchi);
    inchis.push_back(inchi);
    keys.push_back(key);
  }

  std::vector<std::future<void>> tg;
  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(std::async(std::launch::async, runblock, std::ref(mols),
                               count, i, std::ref(inchis), std::ref(keys)));
  }
  for (auto &fut : tg) {
    fut.get();
  }

  for (auto &mol : mols) {
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMultiThread2() {
  const char *smiles_init[] = {
      "C#Cc1ccc2c(c1)C(C(=O)Nc1nncn1C1CC1)CCO2",
      "C=C(C(=O)C1CCOc2ccc(Cl)cc21)c1cncn1C1CC1",
      "C=CC(=O)N(C(=O)[C@@H]1CCOc2ccc(Cl)cc21)c1nncn1C1CC1",
      "C=CC(=O)N[C@@]1(C(=O)Nc2nncn2C2CC2)CCOc2ccc(Cl)cc21",
      "CC#CC(=O)N(C(=O)[C@@H]1CCOc2ccc(Cl)cc21)c1nncn1C1CC1",
      "CC(C(=O)Nc1nncn1C1CC1)c1cccc(Cl)c1",
      "CC(C)c1scnc1NC(=O)C1CCOc2ccc(Cl)cc21",
      "CC(C)n1c(NC(=O)C2CCOc3ccccc32)nc2ccccc21",
  };
  const char *inchikeys_init[] = {
      "FSPOGANMWXTKGJ-UHFFFAOYSA-N", "CIBQANUIVWKYMR-UHFFFAOYSA-N",
      "GBPCXKPRLGZQQB-CYBMUJFWSA-N", "IOOIGOWFYCNXKX-SFHVURJKSA-N",
      "XSYMIYIQKSFGOG-CQSZACIVSA-N", "BLJFTMZJZUGEGJ-UHFFFAOYSA-N",
      "CEVCYTUHJSCLCQ-UHFFFAOYSA-N", "FHPRHOONYCYTIX-UHFFFAOYSA-N",
  };

  auto lambda = [](const std::string &smiles,
                   const std::string &expected_inchikey) {
    ROMol *mol = SmilesToMol(smiles);

    std::string key = MolToInchiKey(*mol);
    TEST_ASSERT(key == expected_inchikey);

    delete mol;
  };

  const std::vector<std::string> smiles(smiles_init, std::end(smiles_init));
  const std::vector<std::string> inchikeys(inchikeys_init,
                                           std::end(inchikeys_init));

  std::vector<std::future<void>> futures;
  for (size_t i = 0; i < smiles.size(); ++i) {
    futures.emplace_back(std::async(std::launch::async, lambda,
                                    std::ref(smiles[i]),
                                    std::ref(inchikeys[i])));
  }
  for (auto &fut : futures) {
    fut.get();
  }
}
#else
void testMultiThread() {
  {}
}
void testMultiThread2() {
  {}
}
#endif

void testGithubIssue3() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github issue 3: bad inchis when mol has no stereoinfo"
      << std::endl;
  {
    std::string fName = getenv("RDBASE");
    fName += "/External/INCHI-API/test_data/github3.mol";
    auto *m = static_cast<ROMol *>(MolFileToMol(fName));
    TEST_ASSERT(m);
    std::string smi = MolToSmiles(*m, true);
    TEST_ASSERT(smi == "CNC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO");

    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi ==
                "InChI=1S/C7H17NO5/c1-8-2-4(10)6(12)7(13)5(11)3-9/"
                "h4-13H,2-3H2,1H3/t4-,5+,6+,7+/m0/s1");

    // blow out the stereo information with a copy:
    auto *m2 = new RWMol(*m);
    m2->clearComputedProps();
    MolOps::sanitizeMol(*m2);

    inchi = MolToInchi(*m2, tmp);
    TEST_ASSERT(inchi ==
                "InChI=1S/C7H17NO5/c1-8-2-4(10)6(12)7(13)5(11)3-9/"
                "h4-13H,2-3H2,1H3/t4-,5+,6+,7+/m0/s1");

    delete m;
    delete m2;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue8() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing a consequence of the fix for github issue "
                          "8: bad mols from some inchis"
                       << std::endl;
  {
    std::string fName = getenv("RDBASE");
    fName += "/External/INCHI-API/test_data/github8_extra.mol";
    auto *m = static_cast<ROMol *>(MolFileToMol(fName, true, false));
    TEST_ASSERT(m);

    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi ==
                "InChI=1S/C20H13IN2O3/"
                "c21-15-8-14-18(9-16(15)23)26-17-7-10(22)5-6-13(17)19(14)11-3-"
                "1-2-4-12(11)20(24)25/h1-9,23H,22H2,(H,24,25)/b23-16+/i21-2");

    ExtraInchiReturnValues tmp2;
    ROMol *m2 = InchiToMol(inchi, tmp2);
    TEST_ASSERT(m2);
    std::string smi = MolToSmiles(*m2, true);
    TEST_ASSERT(smi ==
                "[H]/N=c1\\cc2oc3cc(N)ccc3c(-c3ccccc3C(=O)O)c-2cc1[125I]");

    inchi = MolToInchi(*m2, tmp2);
    TEST_ASSERT(inchi ==
                "InChI=1S/C20H13IN2O3/"
                "c21-15-8-14-18(9-16(15)23)26-17-7-10(22)5-6-13(17)19(14)11-3-"
                "1-2-4-12(11)20(24)25/h1-9,23H,22H2,(H,24,25)/b23-16+/i21-2");
    delete m;
    delete m2;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue40() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github issue 40: bad MWs from inchis"
                       << std::endl;
  {
    ExtraInchiReturnValues tmp;
    std::string inchi =
        "InChI=1S/C10H9N3O/c1-7-11-10(14)9(13-12-7)8-5-3-2-4-6-8/"
        "h2-6H,1H3,(H,11,12,14)";
    ROMol *m = InchiToMol(inchi, tmp);
    TEST_ASSERT(m);

    double mw = Descriptors::calcAMW(*m);
    TEST_ASSERT(feq(mw, 187.202));

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue67() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github issue 67: seg fault from inchi"
                       << std::endl;
  {
    ExtraInchiReturnValues tmp;
    std::string inchi =
        "InChI=1S/C18H17N3/"
        "c19-18(20)21-14-17-12-10-16(11-13-17)9-5-4-8-15-6-2-1-3-7-15/"
        "h1-3,6-13H,14H2,(H4,19,20,21)/b9-8+";
    ROMol *m = InchiToMol(inchi, tmp);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 21);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue68() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github issue 68: hang while reading InChI"
                       << std::endl;
  {
    ExtraInchiReturnValues tmp;
    std::string inchi =
        "InChI=1S/C11H20NSSi2.Li/c1-14(2)8-9-15(3,4)12(14)10-11-6-5-7-13-11;/"
        "h6-7H,8-10H2,1-4H3;/q-1;+1";
    BOOST_LOG(rdInfoLog) << "  parse 1:" << std::endl;
    ROMol *m = InchiToMol(inchi, tmp);
    BOOST_LOG(rdInfoLog) << "  done" << std::endl;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 16);

    delete m;
  }
  {
    ExtraInchiReturnValues tmp;
    std::string inchi =
        "InChI=1S/C12H22NSSi2.Li/"
        "c1-15(2)10-11-16(3,4)13(15)8-7-12-6-5-9-14-12;/"
        "h6,9H,7-8,10-11H2,1-4H3;/q-1;+1";
    BOOST_LOG(rdInfoLog) << "  parse 2:" << std::endl;
    ROMol *m = InchiToMol(inchi, tmp);
    BOOST_LOG(rdInfoLog) << "  done" << std::endl;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 17);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue296() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github issue 296: problems with chiral S and inchi"
      << std::endl;
  {
    std::string smiles = "C[S@@](=O)C(C)(C)C";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/C5H12OS/c1-5(2,3)7(4)6/h1-4H3/t7-/m1/s1");
    delete m;
  }

  {
    std::string fName = getenv("RDBASE");
    fName += "/External/INCHI-API/test_data/github296.mol";
    auto *m = static_cast<ROMol *>(MolFileToMol(fName));
    TEST_ASSERT(m);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/C5H12OS/c1-5(2,3)7(4)6/h1-4H3/t7-/m1/s1");
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue437() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github issue 437: problems with InChI and "
                          "ring stereochemistry"
                       << std::endl;

  // start with some general chirality checks
  {
    std::string smiles = "C[C@@H](F)Cl";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi == "InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3/t2-/m1/s1");

    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    // std::cerr<<" smi1: "<<smi1<<std::endl;
    // std::cerr<<" smi2: "<<smi2<<std::endl;
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "O=[S@](Cl)F";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi == "InChI=1S/ClFOS/c1-4(2)3/t4-/m1/s1");
    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    // std::cerr<<" smi1: "<<smi1<<std::endl;
    // std::cerr<<" smi2: "<<smi2<<std::endl;
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "C[C@@H](F)[C@H](C)[C@H](C)F";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C6H12F2/c1-4(5(2)7)6(3)8/h4-6H,1-3H3/t4-,5+,6-");

    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "C[C@@H](F)[C@H](C)[C@@H](C)F";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C6H12F2/c1-4(5(2)7)6(3)8/h4-6H,1-3H3/t5-,6-/m1/s1");

    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    // std::cerr<<" smi1: "<<smi1<<std::endl;
    // std::cerr<<" smi2: "<<smi2<<std::endl;
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "C[C@H](F)[C@H](C)[C@H](C)F";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C6H12F2/c1-4(5(2)7)6(3)8/h4-6H,1-3H3/t5-,6-/m0/s1");

    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "C[C@H](F)[C@H](C)[C@@H](C)F";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C6H12F2/c1-4(5(2)7)6(3)8/h4-6H,1-3H3/t4-,5-,6+");

    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }

  // now the rings
  {
    std::string smiles = "C[C@@H]1CC[C@@H](O)C[C@@H]1F";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C7H13FO/c1-5-2-3-6(9)4-7(5)8/h5-7,9H,2-4H2,1H3/"
                "t5-,6-,7+/m1/s1");

    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }

  {
    std::string smiles = "C[C@@]1(F)CC[C@@](O)(F)CC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(
        inchi ==
        "InChI=1S/C7H12F2O/c1-6(8)2-4-7(9,10)5-3-6/h10H,2-5H2,1H3/t6-,7+");
    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "C[C@H]1CC[C@H](O)CC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C7H14O/c1-6-2-4-7(8)5-3-6/h6-8H,2-5H2,1H3/t6-,7-");
    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    std::string smiles = "C[C@H]1CC[C@](O)(F)CC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    // std::cerr<<" result: "<<inchi<<std::endl;
    TEST_ASSERT(inchi ==
                "InChI=1S/C7H13FO/c1-6-2-4-7(8,9)5-3-6/h6,9H,2-5H2,1H3/t6-,7+");
    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  {
    // this was the original molecule in the bug report
    std::string smiles =
        "OC(=NCCc1ccccc1)[C@H]1CC[C@H](Cn2c(O)nc3ccccc3c2=O)CC1";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    std::string smi1 = MolToSmiles(*m, true);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi ==
                "InChI=1S/C24H27N3O3/"
                "c28-22(25-15-14-17-6-2-1-3-7-17)19-12-10-18(11-13-19)16-27-23("
                "29)20-8-4-5-9-21(20)26-24(27)30/"
                "h1-9,18-19H,10-16H2,(H,25,28)(H,26,30)/t18-,19-");
    delete m;
    m = InchiToMol(inchi, tmp);
    std::string smi2 = MolToSmiles(*m, true);
    delete m;
    TEST_ASSERT(smi1 == smi2);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue562() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github issue 562: InChI radicals not properly converted"
      << std::endl;

  {
    std::string inchi = "InChI=1S/HO/h1H";
    ExtraInchiReturnValues tmp;
    ROMol *m = InchiToMol(inchi, tmp);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumExplicitHs() == 1);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNoImplicit() == true);

    std::string oinchi = MolToInchi(*m, tmp);
    TEST_ASSERT(oinchi == inchi);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue614() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github issue 614: segfault from MolToInchi "
                          "when bad bond stereochem info is present"
                       << std::endl;

  {
    std::string smiles = "C/C=C/C";
    RWMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    m->removeAtom(3);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    std::cerr << " inchi: " << inchi << std::endl;
    TEST_ASSERT(inchi == "InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithubIssue1572() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github issue 1572: Incorrect conversion from InChI"
      << std::endl;

  {
    std::string inchi =
        "InChI=1S/C16H10N6O2/"
        "c23-21(15-9-17-11-5-1-3-7-13(11)19-15)22(24)16-10-18-12-6-2-4-8-14(12)"
        "20-16/h1-10H";
    ExtraInchiReturnValues tmp;
    RWMol *m = InchiToMol(inchi, tmp, false);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 24);
    m->updatePropertyCache(false);
    TEST_ASSERT(m->getAtomWithIdx(20)->getAtomicNum() == 7);
    TEST_ASSERT(m->getAtomWithIdx(21)->getAtomicNum() == 7);
    TEST_ASSERT(m->getAtomWithIdx(22)->getAtomicNum() == 8);
    TEST_ASSERT(m->getAtomWithIdx(23)->getAtomicNum() == 8);
    TEST_ASSERT(m->getAtomWithIdx(20)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(21)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(22)->getFormalCharge() == 0);
    TEST_ASSERT(m->getAtomWithIdx(23)->getFormalCharge() == 0);

    // and the molecule should be sanitizable:
    MolOps::sanitizeMol(*m);
    TEST_ASSERT(m->getAtomWithIdx(20)->getFormalCharge() == 1);
    TEST_ASSERT(m->getAtomWithIdx(21)->getFormalCharge() == 1);
    TEST_ASSERT(m->getAtomWithIdx(22)->getFormalCharge() == -1);
    TEST_ASSERT(m->getAtomWithIdx(23)->getFormalCharge() == -1);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMolBlockToInchi() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing mol block to InChI" << std::endl;

  {
    std::string molb = R"MOL(
  Mrv1824 02111920092D          

  6  6  0  0  0  0            999 V2000
   -5.5134    3.5259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    3.1134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5134    1.8759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    3.1134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  2  0  0  0  0
M  END
    )MOL";
    ExtraInchiReturnValues tmp;
    std::string inchi = MolBlockToInchi(molb, tmp);
    TEST_ASSERT(inchi == "InChI=1S/C5H8O/c1-2-4-6-5-3-1/h1-2H,3-5H2");
  }
  {
    std::string molb = R"MOL(BDBM163075
     RDKit          2D

 27 30  0  0  0  0  0  0  0  0999 V2000
    1.6146   -5.5162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9260   -4.0489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0484    1.6535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2594    1.8470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4379    3.0237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1670    0.4398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6489    0.4769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3781    3.0608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9796   -3.6396    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9460    3.1800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6752    0.5961    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1571    0.3205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8863    2.9045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0946   -2.6362    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4878   -3.7959    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5729    2.1226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7839    1.3780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0647    1.9663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2758    1.5343    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0711   -1.5187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6916    0.9089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6816   -0.1486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4207   -1.6751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1734    0.0078    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1997    1.0652    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3020   -0.4613    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8110   -3.0456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2 27  1  0
  3  5  2  0
  3  6  1  0
  4  7  2  0
  4  8  1  0
  5 10  1  0
  6 11  2  0
  7 12  1  0
  8 13  2  0
  9 14  2  0
  9 15  1  0
 10 18  2  0
 11 18  1  0
 12 19  2  0
 13 19  1  0
 14 20  1  0
 15 27  1  0
 16 18  1  0
 16 21  1  0
 17 19  1  0
 17 24  1  0
 20 22  2  0
 20 23  1  0
 21 25  2  0
 21 26  1  0
 22 24  1  0
 22 25  1  0
 23 26  2  0
 23 27  1  0
M  CHG  2  15  -1  27   1
M  END
)MOL";

    ExtraInchiReturnValues tmp;
    std::string inchi = MolBlockToInchi(molb, tmp);
    TEST_ASSERT(inchi ==
                "InChI=1S/C23H23N3O/"
                "c1-2-27-15-9-14-20-22(24-17-19-12-7-4-8-13-19)25-21(26-23(20)"
                "27)16-18-10-5-3-6-11-18/h3-15H,2,16-17H2,1H3,(H,24,25,26)");
    TEST_ASSERT(tmp.messagePtr ==
                "Charges were rearranged; Accepted unusual valence(s): O(4)");
  }

  {
    std::string molb = R"MOL(
  Mrv1824 02121905282D          

 10 11  0  0  0  0            999 V2000
   -4.6875   -1.1393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4020   -1.5518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4020   -2.3768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6875   -2.7893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9730   -2.3768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9730   -1.5518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2586   -2.7893    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5441   -1.5518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5441   -2.3768    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9608   -0.9684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  1  6  1  0  0  0  0
  7  9  1  0  0  0  0
  6  8  1  0  0  0  0
  7  5  1  0  0  0  0
  8 10  1  0  0  0  0
  8  9  2  0  0  0  0
M  END
    )MOL";
    {
      ExtraInchiReturnValues tmp;
      std::string inchi = MolBlockToInchi(molb, tmp);
      TEST_ASSERT(
          inchi ==
          "InChI=1S/C8H8N2/c1-6-7-4-2-3-5-8(7)10-9-6/h2-5H,1H3,(H,9,10)");
    }
    {
      ExtraInchiReturnValues tmp;
      std::string inchi = MolBlockToInchi(molb, tmp, "/FixedH");
      TEST_ASSERT(
          inchi ==
          "InChI=1/C8H8N2/c1-6-7-4-2-3-5-8(7)10-9-6/h2-5H,1H3,(H,9,10)/f/h10H");
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub3365() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github #3365: problems with high radical counts" << std::endl;

  {
    auto m = "[C]"_smiles;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 4);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/C");
  }
  {
    auto m = "[CH]"_smiles;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 3);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/CH/h1H");
  }
  {
    auto m = "[CH2]"_smiles;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 2);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/CH2/h1H2");
  }
  {
    auto m = "[CH3]"_smiles;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/CH3/h1H3");
  }
  {
    auto m = "C[SH](C)=O"_smiles;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/C2H7OS/c1-4(2)3/h4H,1-2H3");
  }
  {
    auto m = "C[SH](C)(O)O"_smiles;
    TEST_ASSERT(m);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 1);
    ExtraInchiReturnValues tmp;
    std::string inchi = MolToInchi(*m, tmp);
    TEST_ASSERT(inchi == "InChI=1S/C2H9O2S/c1-5(2,3)4/h3-5H,1-2H3");
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub3645() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github #3645: Seg fault when parsing InChI"
                       << std::endl;

  {
    std::string inchi =
        "InChI=1S/C9H22O3SSi/c1-5-9-13(14,10-6-2,11-7-3)12-8-4/h5-9H2,1-4H3";
    ExtraInchiReturnValues tmp;
    bool ok = false;
    try {
      std::unique_ptr<ROMol> m{InchiToMol(inchi, tmp)};
    } catch (const MolSanitizeException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  {
    std::string inchi =
        "InChI=1S/C8H20O3SSi/c1-5-9-12(13,8-4,10-6-2)11-7-3/h5-8H2,1-4H3";
    ExtraInchiReturnValues tmp;
    bool ok = false;
    try {
      std::unique_ptr<ROMol> m{InchiToMol(inchi, tmp)};
    } catch (const MolSanitizeException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test_clean_up_on_kekulization_error() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "To INCHI conversion leaks on kekulization failure"
                       << std::endl;

  {
    // This will fail kekulization because of the ambiguous valence
    // on the N atom
    SmilesParserParams params;
    params.sanitize = false;

    std::unique_ptr<ROMol> mol{SmilesToMol("c1ccnc1", params)};
    TEST_ASSERT(mol);

    ExtraInchiReturnValues tmp;
    bool ok = true;
    try {
      auto inchi = MolToInchi(*mol, tmp);
    } catch (const std::exception &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub6172() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github #6172: bogus bond type when parsing InChI"
      << std::endl;

  {
    std::string inchi =
        "InChI=1S/2C18H36N3OP.2CH3.2ClH.Sn/c2*22-23(19-16-10-4-1-5-11-16,20-17-12-6-2-7-13-17)21-18-14-8-3-9-15-18;;;;;/h2*16-18H,1-15H2,(H3,19,20,21,22);2*1H3;2*1H;/q;;2*-1;;;+4/p-2";
    ExtraInchiReturnValues tmp;
    std::unique_ptr<ROMol> m{InchiToMol(inchi, tmp)};
    TEST_ASSERT(!m);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub5311() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing github #5311" << std::endl;

  {
    auto m = "O[PH2]=O"_smiles;
    TEST_ASSERT(m);
    ExtraInchiReturnValues tmp;
    auto inchi = MolToInchi(*m, tmp);
    std::cerr << "!!! " << inchi << std::endl;
    TEST_ASSERT(inchi == "InChI=1S/H3O2P/c1-3-2/h3H2,(H,1,2)");
    BOOST_LOG(rdInfoLog) << "done" << std::endl;
  }
}
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();

  testGithubIssue3();
  testGithubIssue8();
  testGithubIssue40();
  testGithubIssue67();
  testGithubIssue68();
  testGithubIssue296();
  testMultiThread();
  testMultiThread2();
  testGithubIssue437();
  testGithubIssue614();
  testGithubIssue1572();
  testMolBlockToInchi();
  testGithubIssue562();
  testGithub3365();
  testGithub3645();
  test_clean_up_on_kekulization_error();
  testGithub6172();
  testGithub5311();
}
