//
//  Copyright (C) 2018 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <clocale>
#include <cstdlib>

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

    std::vector<boost::shared_ptr<RWMol>> mols =
        MolInterchange::JSONDataStreamToMols(&inStream);
    TEST_ASSERT(mols.size() == 1);
    RWMol *m = mols[0].get();
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

    std::vector<boost::shared_ptr<RWMol>> mols =
        MolInterchange::JSONDataStreamToMols(&inStream);
    TEST_ASSERT(mols.size() == 1);
    RWMol *m = mols[0].get();
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
  }

  {
    std::string json =
        "{\"moljson-header\": {\"version\": 10, \"name\": \"example "
        "molecules\"}, \"atomDefaults\": {\"chg\": 0, \"impHs\": 0, "
        "\"stereo\": \"unspecified\", \"nrad\": 0, \"Z\": 6}, "
        "\"bondDefaults\": {\"bo\": 1, \"stereo\": \"unspecified\", "
        "\"stereoAtoms\": []}, \"molecules\": [{\"name\": \"no name\", "
        "\"atoms\": [{\"Z\": 6, \"impHs\": 2}, {\"Z\": 8}, {\"Z\": 26}], "
        "\"bonds\": [{\"atoms\": [0, 1], \"bo\": 2}, {\"atoms\": [1, 2], "
        "\"bo\": 0}], \"representations\": [{\"format_version\": 1, "
        "\"toolkit\": \"RDKit\", \"toolkit_version\": \"2018.03.1.dev1\", "
        "\"aromaticAtoms\": [], \"aromaticBonds\": [], \"cipRanks\": [0, 1, "
        "2], \"cipCodes\": [], \"atomRings\": []}]}]}";
    std::vector<boost::shared_ptr<RWMol>> mols =
        MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(mols.size() == 1);
    RWMol *m = mols[0].get();
    TEST_ASSERT(m);
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::ZERO);
  }
#endif

  {
    std::string json =
        "{\"moljson-header\":{\"version\":10,\"name\":\"test2 mols\"},"
        "\"atomDefaults\":{\"Z\":6,\"impHs\":3,\"chg\":0,\"nRad\":0,"
        "\"isotope\":0,"
        "\"stereo\":\"unspecified\"},\"bondDefaults\":{\"bo\":1,\"stereo\":"
        "\"unspecified\"},"
        "\"molecules\":[{\"name\":\"mol1 "
        "name\",\"atoms\":[{},{}],\"bonds\":[{\"bo\":1, \"atoms\":[0, 1]}]}]}";
    std::vector<boost::shared_ptr<RWMol>> mols =
        MolInterchange::JSONDataToMols(json);
    TEST_ASSERT(mols.size() == 1);
    RWMol *m = mols[0].get();
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
  auto json = MolInterchange::MolToJSONData(*mol, "test2 mol2");
  std::cerr << json << std::endl;
  std::string smi1 = MolToSmiles(*mol);
  auto newMols = MolInterchange::JSONDataToMols(json);
  TEST_ASSERT(newMols.size() == 1);
  // mol->debugMol(std::cerr);
  // newMols[0]->debugMol(std::cerr);
  std::string smi2 = MolToSmiles(*newMols[0]);
  if (smi1 != smi2) {
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
    auto json = MolInterchange::MolToJSONData(*mol, "test2 mols");
    std::cerr << json << std::endl;
  }
#endif
  roundtripSmi("F[C@@](Cl)(O)C");
  roundtripSmi("c1ccccc1");
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

void RunTests() {
#if 1
  // test1();
  test2();
#endif
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
