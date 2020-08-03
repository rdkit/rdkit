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
#include "MultithreadedSmilesMolSupplier.h"

namespace io = boost::iostreams;
using namespace RDKit;
using namespace std::chrono;

struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

void testSmiConcurrent(std::string path, std::string delimiter,
                       int smilesColumn, int nameColumn, bool titleLine,
                       bool sanitize, int numWriterThreads,
                       size_t sizeInputQueue, size_t sizeOutputQueue,
                       unsigned int expectedResult) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  std::atomic<unsigned int> nMols(0);
  MultithreadedSmilesMolSupplier sup(fname, delimiter, smilesColumn, nameColumn,
                                     titleLine, sanitize, numWriterThreads,
                                     sizeInputQueue, sizeOutputQueue);
  std::string last_smiles;
  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      ++nMols;
    }
    delete mol;
  }
  TEST_ASSERT(sup.getLastRecordId() == expectedResult);
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

void testSmiCorrectness() {
  /*
          TEST CORRECTNESS
  */
  std::string path = "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  unsigned int expectedResult = 10;
  testSmiConcurrent(path, ",", 1, 0, false, true, 2, 5, 5, expectedResult);
  testSmiConcurrent(path, ",", 1, 0, false, true, 3, 5, 5, expectedResult);
  testSmiConcurrent(path, ",", 1, 0, false, true, 4, 5, 10, expectedResult);
  testSmiConcurrent(path, ",", 1, 0, false, true, 1, 10, 10, expectedResult);
  testSmiConcurrent(path, ",", 1, 0, false, true, 3, 1, 1, expectedResult);
}

void testSmiPerformance() {
  /*
          TEST PERFORMANCE
  */
  std::string path = "/Regress/Data/znp.50k.smi";
  unsigned int expectedResult = 50000;

  auto start = high_resolution_clock::now();
  testSmiOld(path, " \t", 0, 1, false, true, expectedResult);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  std::cout << "Duration for SmilesMolSupplier: " << duration.count()
            << " (seconds) \n";

  unsigned int maxThreadCount = std::thread::hardware_concurrency() + 4;
  std::cout << "Maximum Threads available: " << maxThreadCount << "\n";
  for (unsigned int i = maxThreadCount; i >= 1; --i) {
    start = high_resolution_clock::now();
    testSmiConcurrent(path, " \t", 0, 1, false, true, i, 1000, 100,
                      expectedResult);
    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "Duration for testSmiConcurent with " << i
              << "  writer threads: " << duration.count() << " (seconds) \n";
  }
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSmiCorrectness();
  BOOST_LOG(rdErrorLog) << "Finished: testSmiCorrectness()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  /*
    BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
    testSmiPerformance();
    BOOST_LOG(rdErrorLog) << "Finished: testSmiPerformance()\n";
    BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";
  */
  return 0;
}
