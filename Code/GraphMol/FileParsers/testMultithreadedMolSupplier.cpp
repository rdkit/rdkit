//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/ConcurrentQueue.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/test.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>

#include "MultithreadedMolSupplier.h"
#include "MultithreadedSDMolSupplier.h"
#include "MultithreadedSmilesMolSupplier.h"

namespace io = boost::iostreams;
using namespace RDKit;
using namespace std::chrono;

// thread safe printing for debugging
// Usage Example: PrintThread{} << "something";
struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

void testSmiConcurrent(std::string path, std::string delimiter,
                       int smilesColumn, int nameColumn, bool titleLine,
                       bool sanitize, unsigned int numWriterThreads,
                       size_t sizeInputQueue, size_t sizeOutputQueue,
                       unsigned int expectedResult) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  unsigned int nMols = 0;
  boost::dynamic_bitset<> bitVector(expectedResult);
  MultithreadedSmilesMolSupplier sup(fname, delimiter, smilesColumn, nameColumn,
                                     titleLine, sanitize, numWriterThreads,
                                     sizeInputQueue, sizeOutputQueue);
  // we have not called the next method yet
  TEST_ASSERT(sup.getLastRecordId() == 0);
  // initially no bit is set in the bitVector, sanity check
  TEST_ASSERT(!bitVector.any());

  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    unsigned int id = sup.getLastRecordId();
    bitVector[id - 1] = 1;
    if (mol != nullptr) {
      ++nMols;
    }
    delete mol;
  }
  // if all bits are set then we have seen possible ids
  TEST_ASSERT(bitVector.all());
  TEST_ASSERT(nMols == expectedResult);
}

void testSDConcurrent(std::string path, bool sanitize, bool removeHs,
                      bool strictParsing, unsigned int numWriterThreads,
                      size_t sizeInputQueue, size_t sizeOutputQueue,
                      unsigned int expectedResult) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  unsigned int nMols = 0;
  boost::dynamic_bitset<> bitVector(expectedResult);
  MultithreadedSDMolSupplier sup(fname, sanitize, removeHs, strictParsing,
                                 numWriterThreads, sizeInputQueue,
                                 sizeOutputQueue);
  // we have not called the next method yet
  TEST_ASSERT(sup.getLastRecordId() == 0);
  // initially no bit is set in the bitVector, sanity check
  TEST_ASSERT(!bitVector.any());
  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      unsigned int id = sup.getLastRecordId();
      bitVector[id - 1] = 1;
      ++nMols;
    }
    delete mol;
  }
  // if all bits are set then we have seen possible ids
  TEST_ASSERT(bitVector.all());
  TEST_ASSERT(nMols == expectedResult);
}

void testSmiOld(std::string path, std::string delimiter, int smilesColumn,
                int nameColumn, bool titleLine, bool sanitize,
                unsigned int expectedResult) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  unsigned int numMols = 0;
  SmilesMolSupplier sup(fname, delimiter, smilesColumn, nameColumn, titleLine,
                        sanitize);
  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      ++numMols;
    }
    delete mol;
  }
  TEST_ASSERT(numMols == expectedResult);
}

void testSDOld(std::string path, bool sanitize, bool removeHs,
               bool strictParsing, unsigned int expectedResult) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  unsigned int numMols = 0;
  SDMolSupplier sup(fname, sanitize, removeHs, strictParsing);
  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      ++numMols;
    }
    delete mol;
  }
  TEST_ASSERT(numMols == expectedResult);
}

void testSmiProperties() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  std::vector<std::string> nameVector, tpsaVector;
  std::string tempStr;
  SmilesMolSupplier sup(fname, ",", 1, 0, true);
  MultithreadedSmilesMolSupplier multiSup(fname, ",", 1, 0, true);
  while (!multiSup.atEnd()) {
    ROMol* mol = multiSup.next();
    if (mol != nullptr) {
      mol->getProp(common_properties::_Name, tempStr);
      nameVector.push_back(tempStr);
      mol->getProp("TPSA", tempStr);
      tpsaVector.push_back(tempStr);
    }
  }

  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      mol->getProp(common_properties::_Name, tempStr);
      TEST_ASSERT(std::find(nameVector.begin(), nameVector.end(), tempStr) !=
                  nameVector.end());
      mol->getProp("TPSA", tempStr);
      TEST_ASSERT(std::find(tpsaVector.begin(), tpsaVector.end(), tempStr) !=
                  tpsaVector.end());
    }
  }
}

void testSmiCorrectness() {
  /*
          TEST CORRECTNESS
  */
  std::string path = "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  unsigned int expectedResult = 10;
  testSmiConcurrent(path, ",", 1, 0, false, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  expectedResult = 10;
  testSmiConcurrent(path, ",", 1, 0, true, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/first_200.tpsa.csv";
  expectedResult = 200;
  testSmiConcurrent(path, ",", 0, -1, true, true, 2, 5, 5, expectedResult);

  path = "/Regress/Data/znp.50k.smi";
  expectedResult = 50000;
  testSmiConcurrent(path, " \t", 0, 1, false, true, 4, 1000, 100,
                    expectedResult);
  /*
          TEST PROPERTIES
  */
  testSmiProperties();
}

void testSDProperties() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  std::vector<std::string> nameVector;
  std::string tempStr;
  SDMolSupplier sup(fname, false);
  MultithreadedSDMolSupplier multiSup(fname, false);

  while (!multiSup.atEnd()) {
    ROMol* mol = multiSup.next();
    if (mol != nullptr) {
      TEST_ASSERT(mol->hasProp(common_properties::_Name));
      mol->getProp(common_properties::_Name, tempStr);
      nameVector.push_back(tempStr);
      TEST_ASSERT(mol->hasProp("NCI_AIDS_Antiviral_Screen_Conclusion"));
      TEST_ASSERT(mol->hasProp("CAS_RN"));
      TEST_ASSERT(mol->hasProp("NSC"));
    }
  }

  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      mol->getProp(common_properties::_Name, tempStr);
      TEST_ASSERT(std::find(nameVector.begin(), nameVector.end(), tempStr) !=
                  nameVector.end());
    }
  }
}

void testSDCorrectness() {
  /*
          TEST CORRECTNESS
  */

  std::string path = "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  unsigned int expectedResult = 16;
  testSDConcurrent(path, false, true, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/outNCI_few_molsupplier.sdf";
  expectedResult = 16;
  testSDConcurrent(path, true, true, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/esters_end.sdf";
  expectedResult = 6;
  testSDConcurrent(path, true, true, true, 2, 5, 5, expectedResult);

  // strict parsing results in reading 0 out of 2 records
  path = "/Code/GraphMol/FileParsers/test_data/strictLax1.sdf";
  expectedResult = 0;
  testSDConcurrent(path, true, true, true, 2, 5, 5, expectedResult);

  expectedResult = 2;
  testSDConcurrent(path, true, true, false, 2, 5, 5, expectedResult);

  path = "/Regress/Data/O2.sdf";
  expectedResult = 1;
  testSDConcurrent(path, true, true, false, 2, 5, 5, expectedResult);

  path = "/Regress/Data/mols1000.sdf";
  expectedResult = 1000;
  testSDConcurrent(path, false, true, true, 2, 5, 5, expectedResult);

  /*
          TEST PROPERTIES
  */
  testSDProperties();
}

void testPerformance() {
  /*
          TEST PERFORMANCE
  */
  std::string path = "/Regress/Data/znp.50k.smi";
  unsigned int expectedResult = 50000;

  auto start = high_resolution_clock::now();
  testSmiOld(path, " \t", 0, 1, false, true, expectedResult);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
  std::cout << "Duration for SmilesMolSupplier: " << duration.count()
            << " (milliseconds) \n";

  unsigned int maxThreadCount = std::thread::hardware_concurrency() + 4;
  std::cout << "Maximum Threads available: " << maxThreadCount << "\n";
  for (unsigned int i = maxThreadCount; i >= 1; --i) {
    start = high_resolution_clock::now();
    testSmiConcurrent(path, " \t", 0, 1, false, true, i, 1000, 100,
                      expectedResult);
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration for testSmiConcurent with " << i
              << "  writer threads: " << duration.count()
              << " (milliseconds) \n";
  }

  path = "/Regress/Data/chembl26_very_active.sdf";
  expectedResult = 35767;

  start = high_resolution_clock::now();
  testSDOld(path, false, true, false, expectedResult);
  stop = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(stop - start);
  std::cout << "Duration for SDMolSupplier: " << duration.count()
            << " (milliseconds) \n";

  for (unsigned int i = maxThreadCount; i >= 1; --i) {
    start = high_resolution_clock::now();
    testSDConcurrent(path, false, true, true, i, 1000, 4000, expectedResult);
    stop = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration for testSDConcurent with " << i
              << "  writer threads: " << duration.count()
              << " (milliseconds) \n";
  }
}

int main() {
  RDLog::InitLogs();

#ifdef RDK_TEST_MULTITHREADED

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSmiCorrectness();
  BOOST_LOG(rdErrorLog) << "Finished: testSmiCorrectness()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSDCorrectness();
  BOOST_LOG(rdErrorLog) << "Finished: testSDCorrectness()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  /*
    BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
    testPerfomance();
    BOOST_LOG(rdErrorLog) << "Finished: testPerformance()\n";
    BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  */

#endif

  return 0;
}
