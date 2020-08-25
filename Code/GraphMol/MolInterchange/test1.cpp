//
//  Copyright (C) 2018 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <clocale>
#include <cstdlib>
#include <chrono>

#include <string>
#include <fstream>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "test1: basics" << std::endl;

#if 1
  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/MolInterchange/test_data/test1.json";
    std::ifstream inStream(fName);
    if (!inStream || (inStream.bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }

    auto mols = MolInterchange::JSONDataStreamToMols(&inStream);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    // m->debugMol(std::cerr);
    TEST_ASSERT(m->getNumAtoms() == 15);
    TEST_ASSERT(m->getNumBonds() == 15);
    TEST_ASSERT(m->getAtomWithIdx(0)->getIsAromatic());
    TEST_ASSERT(m->getAtomWithIdx(13)->getFormalCharge() == 1);
    TEST_ASSERT(m->getAtomWithIdx(12)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CCW);
    TEST_ASSERT(m->getBondBetweenAtoms(10, 11));
    TEST_ASSERT(m->getBondBetweenAtoms(10, 11)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(10, 11)->getStereo() == Bond::STEREOCIS);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getIsAromatic());
    TEST_ASSERT(m->getNumConformers() == 0);
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                std::string("example 1"));
    TEST_ASSERT(m->hasProp("prop1"));
    TEST_ASSERT(m->getProp<int>("prop1") == 1);
    TEST_ASSERT(m->hasProp("prop2"));
    TEST_ASSERT(feq(m->getProp<double>("prop2"), 3.14));
    TEST_ASSERT(m->hasProp("prop3"));
    TEST_ASSERT(m->getProp<std::string>("prop3") == "foo");
    TEST_ASSERT(m->getRingInfo()->isInitialized());
    TEST_ASSERT(m->getRingInfo()->atomRings().size() == 1);
    TEST_ASSERT(m->getRingInfo()->atomRings()[0].size() == 6);
  }
  {
    std::string fName =
        rdbase + "/Code/GraphMol/MolInterchange/test_data/test2.json";
    std::ifstream inStream(fName);
    if (!inStream || (inStream.bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }

    auto mols = MolInterchange::JSONDataStreamToMols(&inStream);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getAtomWithIdx(1)->getChiralTag() ==
                Atom::CHI_TETRAHEDRAL_CW);
    TEST_ASSERT(m->getNumConformers() == 2);
    TEST_ASSERT(!m->getConformer(0).is3D());
    TEST_ASSERT(m->getConformer(1).is3D());
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                std::string("example 2"));
    TEST_ASSERT(m->getAtomWithIdx(1)->getIsotope() == 0);
    TEST_ASSERT(m->getAtomWithIdx(2)->getIsotope() == 35);
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge));
    TEST_ASSERT(feq(m->getAtomWithIdx(0)->getProp<double>(
                        common_properties::_GasteigerCharge),
                    -0.352));
  }

  {
    std::string json =
        "{\"commonchem\": {\"version\": 10 },"
        " \"defaults\": {\"atom\": {\"chg\": 0, \"impHs\": 0, "
        "\"stereo\": \"unspecified\", \"nrad\": 0, \"z\": 6}, "
        "\"bond\": {\"bo\": 1, \"stereo\": \"unspecified\", "
        "\"stereoAtoms\": []}}, \"molecules\": [{\"name\": \"no name\", "
        "\"atoms\": [{\"z\": 6, \"impHs\": 2}, {\"z\": 8}, {\"z\": 26}], "
        "\"bonds\": [{\"atoms\": [0, 1], \"bo\": 2}, {\"atoms\": [1, 2], "
        "\"bo\": 0}], \"extensions\": [{\"formatVersion\": 1, "
        "\"name\": \"rdkitRepresentation\", \"formatVersion\": 1,"
        "\"toolkitVersion\": \"2018.03.1.dev1\", "
        "\"aromaticAtoms\": [], \"aromaticBonds\": [], \"cipRanks\": [0, 1, "
        "2], \"cipCodes\": [], \"atomRings\": []}]}]}";
    auto mols = MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::ZERO);
  }
#endif

  {
    std::string json =
        "{\"commonchem\":{\"version\":10 },"
        "\"defaults\":{\"atom\":{\"z\":6,\"impHs\":3,\"chg\":0,\"nRad\":0,"
        "\"isotope\":0,"
        "\"stereo\":\"unspecified\"},\"bond\":{\"bo\":1,\"stereo\":"
        "\"unspecified\"}},"
        "\"molecules\":[{\"name\":\"mol1 "
        "name\",\"atoms\":[{},{}],\"bonds\":[{\"bo\":1, \"atoms\":[0, 1]}]}]}";
    auto mols = MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2)
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1));
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void roundtripSmi(const char *smi) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smi));
  TEST_ASSERT(mol);
  mol->setProp("_Name", "test mol");
  auto json = MolInterchange::MolToJSONData(*mol);
  std::cerr << json << std::endl;
  std::string smi1 = MolToSmiles(*mol);
  auto newMols = MolInterchange::JSONDataToMols(json);
  TEST_ASSERT(newMols.size() == 1);
  std::string smi2 = MolToSmiles(*newMols[0]);
  if (smi1 != smi2) {
    mol->debugMol(std::cerr);
    newMols[0]->debugMol(std::cerr);
    std::cerr << "smi1: " << smi1 << std::endl;
    std::cerr << "smi2: " << smi2 << std::endl;
  }
  TEST_ASSERT(smi1 == smi2);
}

void test2() {
  BOOST_LOG(rdInfoLog) << "test2: basic writing" << std::endl;
#if 1
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("CC"));
    TEST_ASSERT(mol);
    mol->setProp("_Name", "mol1 name");
    auto json = MolInterchange::MolToJSONData(*mol);
    std::cerr << json << std::endl;
  }
#endif
  roundtripSmi("F[C@@](Cl)(O)C");
  roundtripSmi("c1ccccc1");
  roundtripSmi("CCC1=C(N)C=C(C)N=C1");
#if 0
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("F[C@](Cl)(O)/C=C/C"));
    TEST_ASSERT(mol);
    mol->setProp("_Name", "test mol");
    mol->getBondBetweenAtoms(4, 5)->setStereo(Bond::STEREOTRANS);
    auto json = MolInterchange::MolToJSONData(*mol, "test2 mol2");
    std::cerr << json << std::endl;
    std::string smi1 = MolToSmiles(*mol);
    auto newMols = MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(newMols.size() == 1);
    mol->debugMol(std::cerr);
    newMols[0]->debugMol(std::cerr);
    std::string smi2 = MolToSmiles(*newMols[0]);
    std::cerr << "smi1: " << smi1 << std::endl;
    std::cerr << "smi2: " << smi2 << std::endl;
    TEST_ASSERT(smi1 == smi2);
  }
#endif
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test3() {
  BOOST_LOG(rdInfoLog) << "test3: writing conformers" << std::endl;
  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/MolInterchange/test_data/test2.json";
    std::ifstream inStream(fName);
    auto mols = MolInterchange::JSONDataStreamToMols(&inStream);
    TEST_ASSERT(mols.size() == 1);
    TEST_ASSERT(mols[0]->getNumConformers() == 2);
    TEST_ASSERT(!mols[0]->getConformer(0).is3D());
    TEST_ASSERT(mols[0]->getConformer(1).is3D());
    std::string json = MolInterchange::MolToJSONData(*mols[0]);
    std::cerr << json << std::endl;
    TEST_ASSERT(json.find("conformers") != std::string::npos);
    auto newMols = MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(newMols.size() == 1);
    TEST_ASSERT(newMols[0]->getNumConformers() == 2);
    TEST_ASSERT(!newMols[0]->getConformer(0).is3D());
    TEST_ASSERT(newMols[0]->getConformer(1).is3D());
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test4() {
  BOOST_LOG(rdInfoLog) << "test4: writing properties" << std::endl;
#if 1
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("CC"));
    TEST_ASSERT(mol);
    mol->setProp("foo_string", "bar");
    mol->setProp("foo_int", 1);
    mol->setProp("foo_double", 1.2);
    auto json = MolInterchange::MolToJSONData(*mol);
    std::cerr << json << std::endl;
    TEST_ASSERT(json.find("foo_string") != std::string::npos);
    TEST_ASSERT(json.find("foo_int") != std::string::npos);
    TEST_ASSERT(json.find("foo_double") != std::string::npos);
    auto newMols = MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(newMols.size() == 1);
    TEST_ASSERT(newMols[0]->hasProp("foo_string"));
    TEST_ASSERT(newMols[0]->getProp<std::string>("foo_string") == "bar");
    TEST_ASSERT(newMols[0]->hasProp("foo_int"));
    TEST_ASSERT(newMols[0]->getProp<int>("foo_int") == 1);
    TEST_ASSERT(newMols[0]->hasProp("foo_double"));
    TEST_ASSERT(newMols[0]->getProp<double>("foo_double") == 1.2);
  }
#endif
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test5() {
  BOOST_LOG(rdInfoLog) << "test5: writing partial charges" << std::endl;

  {
    std::unique_ptr<RWMol> mol(SmilesToMol("CO"));
    TEST_ASSERT(mol);
    mol->getAtomWithIdx(0)->setProp(common_properties::_GasteigerCharge, 0.5);
    mol->getAtomWithIdx(1)->setProp(common_properties::_GasteigerCharge, -0.5);

    auto json = MolInterchange::MolToJSONData(*mol);
    std::cerr << json << std::endl;
    TEST_ASSERT(json.find("partialCharges") != std::string::npos);
    auto newMols = MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(newMols.size() == 1);
    TEST_ASSERT(newMols[0]->getAtomWithIdx(0)->hasProp(
        common_properties::_GasteigerCharge));
    TEST_ASSERT(feq(newMols[0]->getAtomWithIdx(0)->getProp<double>(
                        common_properties::_GasteigerCharge),
                    0.5));
    TEST_ASSERT(newMols[0]->getAtomWithIdx(1)->hasProp(
        common_properties::_GasteigerCharge));
    TEST_ASSERT(feq(newMols[0]->getAtomWithIdx(1)->getProp<double>(
                        common_properties::_GasteigerCharge),
                    -0.5));
  }
}

void benchmarking() {
  BOOST_LOG(rdInfoLog) << "benchmarking performance" << std::endl;
  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/MolInterchange/test_data/znp.50k.smi";
    SmilesMolSupplier suppl(fName);
    std::vector<RWMol *> mols;
    auto smir_t1 = std::chrono::system_clock::now();
    while (mols.size() < 20000) {
      mols.push_back(static_cast<RWMol *>(suppl.next()));
    }
    auto smir_t2 = std::chrono::system_clock::now();
    std::cerr << "construction of " << mols.size() << " took "
              << std::chrono::duration<double>(smir_t2 - smir_t1).count()
              << std::endl;
    for (auto &m : mols) {
      MolOps::Kekulize(*m);
    }
    auto jsonw_t1 = std::chrono::system_clock::now();
    auto json = MolInterchange::MolsToJSONData(mols);
    auto jsonw_t2 = std::chrono::system_clock::now();
    std::cerr << "json generation took "
              << std::chrono::duration<double>(jsonw_t2 - jsonw_t1).count()
              << std::endl;

    auto jsonr_t1 = std::chrono::system_clock::now();
    auto newms = MolInterchange::JSONDataToMols(json);
    auto jsonr_t2 = std::chrono::system_clock::now();
    std::cerr << "json parsing took "
              << std::chrono::duration<double>(jsonr_t2 - jsonr_t1).count()
              << std::endl;
    newms.clear();

    auto pklw_t1 = std::chrono::system_clock::now();
    std::vector<std::string> pkls;
    pkls.reserve(mols.size());
    for (const auto &mol : mols) {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      pkls.push_back(pkl);
    }
    auto pklw_t2 = std::chrono::system_clock::now();
    std::cerr << "pickle generation took "
              << std::chrono::duration<double>(pklw_t2 - pklw_t1).count()
              << std::endl;

    auto pklr_t1 = std::chrono::system_clock::now();
    for (const auto &pkl : pkls) {
      ROMol m;
      MolPickler::molFromPickle(pkl, m);
    }
    auto pklr_t2 = std::chrono::system_clock::now();
    std::cerr << "pickle parsing took "
              << std::chrono::duration<double>(pklr_t2 - pklr_t1).count()
              << std::endl;

    for (auto &m : mols) {
      delete m;
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test6() {
  BOOST_LOG(rdInfoLog) << "testing parse options" << std::endl;
  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MolInterchange/test_data/test3.json";
  std::ifstream inStream(fName);
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  const std::string jsond(std::istreambuf_iterator<char>(inStream), {});
  {
    MolInterchange::JSONParseParameters ps;
    auto mols = MolInterchange::JSONDataToMols(jsond, ps);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getNumConformers() == 2);
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                std::string("example 2"));
    TEST_ASSERT(m->hasProp("prop3"));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge));
  }
  {
    MolInterchange::JSONParseParameters ps;
    ps.parseConformers = false;
    auto mols = MolInterchange::JSONDataToMols(jsond, ps);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getNumConformers() == 0);
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                std::string("example 2"));
    TEST_ASSERT(m->hasProp("prop3"));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge));
  }
  {
    MolInterchange::JSONParseParameters ps;
    ps.parseConformers = false;
    ps.parseProperties = false;
    auto mols = MolInterchange::JSONDataToMols(jsond, ps);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getNumConformers() == 0);
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                std::string("example 2"));  // we always parse the name
    TEST_ASSERT(!m->hasProp("prop3"));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge));
  }
  {
    MolInterchange::JSONParseParameters ps;
    ps.parseProperties = false;
    auto mols = MolInterchange::JSONDataToMols(jsond, ps);
    TEST_ASSERT(mols.size() == 1);
    auto m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 6);
    TEST_ASSERT(m->getNumBonds() == 5);
    TEST_ASSERT(m->getNumConformers() == 2);
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) ==
                std::string("example 2"));  // we always parse the name
    TEST_ASSERT(!m->hasProp("prop3"));
    TEST_ASSERT(
        m->getAtomWithIdx(0)->hasProp(common_properties::_GasteigerCharge));
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub2046() {
  BOOST_LOG(rdInfoLog) << "testing github #2046: CIPRank values from "
                          "JSONDataToMols are not unsigned"
                       << std::endl;
  roundtripSmi("C1CCO[C@H]1F");
  {
    std::unique_ptr<RWMol> mol(SmilesToMol("C1CCO[C@H]1F"));
    TEST_ASSERT(mol);
    mol->setProp("_Name", "mol1 name");
    auto jsond = MolInterchange::MolToJSONData(*mol);
    auto mols = MolInterchange::JSONDataToMols(jsond);
    TEST_ASSERT(mols[0]->getAtomWithIdx(3)->getProp<unsigned int>(
                    RDKit::common_properties::_CIPRank) > 0);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testEitherStereo() {
  BOOST_LOG(rdInfoLog) << "testing 'either' stereochemistry" << std::endl;
  {
    auto mol = "CC=CC/C=C/C"_smiles;
    TEST_ASSERT(mol);
    mol->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    auto jsond = MolInterchange::MolToJSONData(*mol);
    auto mols = MolInterchange::JSONDataToMols(jsond);
    TEST_ASSERT(mols[0]->getBondWithIdx(1)->getStereo() == Bond::STEREOANY);
  }
  {
    auto mol = "CC=CC"_smiles;
    TEST_ASSERT(mol);
    mol->getBondWithIdx(1)->setStereo(Bond::STEREOANY);
    auto jsond = MolInterchange::MolToJSONData(*mol);
    auto mols = MolInterchange::JSONDataToMols(jsond);
    TEST_ASSERT(mols[0]->getBondWithIdx(1)->getStereo() == Bond::STEREOANY);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void RunTests() {
#if 1
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
#endif
  testGithub2046();
  testEitherStereo();

  // benchmarking();
  // test2();
}

#if 0
// must be in German Locale for test...
void testLocaleSwitcher() {
  float d = -1.0;
  char buffer[1024];
  sprintf(buffer, "%0.2f", d);
  if (std::string(buffer) != "-1,00") {
    BOOST_LOG(rdInfoLog) << " ---- no German locale support (skipping) ---- "
                         << std::endl;
    return;
  }

  {
    RDKit::Utils::LocaleSwitcher ls;
    sprintf(buffer, "%0.2f", d);
    CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
    // test locale switcher recursion
    {
      RDKit::Utils::LocaleSwitcher ls;
      sprintf(buffer, "%0.2f", d);
      CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
    }
    // should still be in the "C" variant
    sprintf(buffer, "%0.2f", d);
    CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
  }

  // Should be back in German Locale
  sprintf(buffer, "%0.2f", d);
  CHECK_INVARIANT(std::string(buffer) == "-1,00", "Locale Switcher Fail");
}

#ifdef RDK_TEST_MULTITHREADED

#include <RDGeneral/BoostStartInclude.h>
#include <boost/thread.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace {
void runblock() { testLocaleSwitcher(); }
}
void testMultiThreadedSwitcher() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading Locale Switching"
                        << std::endl;

  boost::thread_group tg;
  unsigned int count = 100;
  for (unsigned int i = 0; i < count; ++i) {
    tg.add_thread(new boost::thread(runblock));
  }
  tg.join_all();
  BOOST_LOG(rdErrorLog) << "    Test multithreading (Done)" << std::endl;
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
}

#else

void testMultiThreadedSwitcher() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << " ---- Multithreaded tests disabled ---- "
                       << std::endl;
}
#endif
#endif

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  RunTests();  // run with C locale

  return 0;
}
