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

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

#if 0
#include <boost/timer/timer.hpp>
void test2() {
  BOOST_LOG(rdInfoLog) << "test2: timing" << std::endl;

  std::string rdbase = getenv("RDBASE");
  {
    std::string fName =
        rdbase + "/Code/GraphMol/MolInterchange/test_data/10000mols.json";
    std::ifstream inStream(fName);
    if (!inStream || (inStream.bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fName;
      throw BadFileException(errout.str());
    }
    std::string data((std::istreambuf_iterator<char>(inStream)),
                     std::istreambuf_iterator<char>());

    std::vector<boost::shared_ptr<RWMol>> mols;
    std::cerr << "from json" << std::endl;
    {
      boost::timer::auto_cpu_timer t;
      mols = MolInterchange::JSONDataToMols(data);
    }
    std::cerr << "done" << std::endl;
    TEST_ASSERT(mols.size() == 10000);
#if 0
    std::vector<std::string> smis;
    for (auto &mol : mols) {
      smis.push_back(mol->getProp<std::string>("_Name"));
    }

    std::vector<boost::shared_ptr<RWMol>> mols2;
    std::cerr << "from smiles" << std::endl;
    {
      boost::timer::auto_cpu_timer t;
      for (const auto &smi : smis) {
        mols2.push_back(boost::shared_ptr<RWMol>(SmilesToMol(smi)));
      }
    }
    std::cerr << "done" << std::endl;
    mols2.clear();
    std::cerr << "from smiles no sanit" << std::endl;
    {
      boost::timer::auto_cpu_timer t;
      for (const auto &smi : smis) {
        mols2.push_back(
            boost::shared_ptr<RWMol>(SmilesToMol(smi, false, false)));
      }
    }
#elif 1
    std::vector<std::string> pkls;
    for (auto &mol : mols) {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      pkls.push_back(pkl);
    }

    std::vector<boost::shared_ptr<RWMol>> mols2;
    std::cerr << "from pickle" << std::endl;
    {
      boost::timer::auto_cpu_timer t;
      for (const auto &pkl : pkls) {
        mols2.push_back(boost::shared_ptr<RWMol>(new RWMol(pkl)));
      }
    }
    std::cerr << "done" << std::endl;

#endif
    std::cerr << "done" << std::endl;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
#endif

void RunTests() {
#if 1
  test1();
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
