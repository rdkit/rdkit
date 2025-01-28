//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <chrono>
#include <memory>

#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/test.h>
#include <RDStreams/streams.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "MultithreadedSDMolSupplier.h"
#include "MultithreadedSmilesMolSupplier.h"

namespace io = boost::iostreams;
using namespace RDKit;
using namespace std::chrono;

// thread safe printing for debugging
// Usage Example: PrintThread{} << "something";
struct PrintThread : public std::stringstream {
  inline static std::mutex cout_mutex;
  ~PrintThread() override {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

void testSmiConcurrent(std::istream *strm, bool takeOwnership,
                       std::string delimiter, int smilesColumn, int nameColumn,
                       bool titleLine, bool sanitize,
                       unsigned int numWriterThreads, size_t sizeInputQueue,
                       size_t sizeOutputQueue, unsigned int expectedResult,
                       bool extras = false) {
  unsigned int nMols = 0;
  boost::dynamic_bitset<> bitVector(expectedResult);
  MultithreadedSmilesMolSupplier sup(
      strm, takeOwnership, delimiter, smilesColumn, nameColumn, titleLine,
      sanitize, numWriterThreads, sizeInputQueue, sizeOutputQueue);
  // we have not called the next method yet
  TEST_ASSERT(sup.getLastRecordId() == 0);
  // initially no bit is set in the bitVector, sanity check
  TEST_ASSERT(!bitVector.any());

  while (!sup.atEnd()) {
    ROMol *mol = sup.next();
    if (mol) {
      unsigned int id = sup.getLastRecordId();
      bitVector[id - 1] = 1;
      if (extras) {
        std::unique_ptr<ExplicitBitVect> fp(
            MorganFingerprints::getFingerprintAsBitVect(*mol, 2, 2048));
      }
      ++nMols;
    }
    delete mol;
  }
  // if all bits are set then we have seen possible ids
  TEST_ASSERT(bitVector.all());
  TEST_ASSERT(nMols == expectedResult);
}

void testSmiConcurrent(std::string path, std::string delimiter,
                       int smilesColumn, int nameColumn, bool titleLine,
                       bool sanitize, unsigned int numWriterThreads,
                       size_t sizeInputQueue, size_t sizeOutputQueue,
                       unsigned int expectedResult, bool extras = false) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  std::istream *strm = new std::ifstream(fname.c_str());
  testSmiConcurrent(strm, true, delimiter, smilesColumn, nameColumn, titleLine,
                    sanitize, numWriterThreads, sizeInputQueue, sizeOutputQueue,
                    expectedResult, extras);
}

void testSmiOld(std::string path, std::string delimiter, int smilesColumn,
                int nameColumn, bool titleLine, bool sanitize,
                unsigned int expectedResult, bool extras = false) {
  unsigned int numMols = 0;
  SmilesMolSupplier sup(path, delimiter, smilesColumn, nameColumn, titleLine,
                        sanitize);
  while (!sup.atEnd()) {
    ROMol *mol = sup.next();
    if (mol) {
      if (extras) {
        std::unique_ptr<ExplicitBitVect> fp(
            MorganFingerprints::getFingerprintAsBitVect(*mol, 2, 2048));
      }
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
    std::unique_ptr<ROMol> mol{multiSup.next()};
    if (mol != nullptr) {
      mol->getProp(common_properties::_Name, tempStr);
      nameVector.push_back(tempStr);
      mol->getProp("TPSA", tempStr);
      tpsaVector.push_back(tempStr);
    }
  }

  while (!sup.atEnd()) {
    std::unique_ptr<ROMol> mol{sup.next()};
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
  std::string path;
  unsigned int expectedResult;
  std::string rdbase = getenv("RDBASE");
  path = "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  expectedResult = 10;
  testSmiConcurrent(path, ",", 1, 0, false, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/fewSmi.2.csv";
  expectedResult = 10;
  testSmiConcurrent(path, ",", 1, 0, true, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/first_200.tpsa.csv";
  expectedResult = 200;
  testSmiConcurrent(path, ",", 0, -1, true, true, 2, 5, 5, expectedResult);

#ifdef RDK_USE_BOOST_IOSTREAMS

  path = rdbase + "/Regress/Data/znp.50k.smi.gz";
  std::istream *strm = new gzstream(path);
  expectedResult = 50000;
  testSmiConcurrent(strm, true, " \t", 0, 1, false, false, 3, 1000, 100,
                    expectedResult);
#endif
  /*

     TEST PROPERTIES

  */
  testSmiProperties();
}

void testSDConcurrent(std::istream *strm, bool takeOwnership, bool sanitize,
                      bool removeHs, bool strictParsing,
                      unsigned int numWriterThreads, size_t sizeInputQueue,
                      size_t sizeOutputQueue, unsigned int expectedResult,
                      bool extras = false) {
  unsigned int nMols = 0;
  boost::dynamic_bitset<> bitVector(expectedResult);
  MultithreadedSDMolSupplier sup(strm, takeOwnership, sanitize, removeHs,
                                 strictParsing, numWriterThreads,
                                 sizeInputQueue, sizeOutputQueue);
  // we have not called the next method yet
  TEST_ASSERT(sup.getLastRecordId() == 0);
  // initially no bit is set in the bitVector, sanity check
  TEST_ASSERT(!bitVector.any());
  while (!sup.atEnd()) {
    ROMol *mol = sup.next();
    if (mol) {
      unsigned int id = sup.getLastRecordId();
      bitVector[id - 1] = 1;
      ++nMols;
      if (extras) {
        std::unique_ptr<ExplicitBitVect> fp(
            MorganFingerprints::getFingerprintAsBitVect(*mol, 2, 2048));
      }
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
                      unsigned int expectedResult, bool extras = false) {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + path;
  std::istream *strm = new std::ifstream(fname.c_str());
  testSDConcurrent(strm, true, sanitize, removeHs, strictParsing,
                   numWriterThreads, sizeInputQueue, sizeOutputQueue,
                   expectedResult, extras);
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
    std::unique_ptr<ROMol> mol{multiSup.next()};
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
    std::unique_ptr<ROMol> mol{sup.next()};
    if (mol != nullptr) {
      mol->getProp(common_properties::_Name, tempStr);
      TEST_ASSERT(std::find(nameVector.begin(), nameVector.end(), tempStr) !=
                  nameVector.end());
    }
  }
}

void testSDOld(std::string path, bool sanitize, bool removeHs,
               bool strictParsing, unsigned int expectedResult,
               bool extras = false) {
  unsigned int numMols = 0;
  SDMolSupplier sup(path, sanitize, removeHs, strictParsing);
  while (!sup.atEnd()) {
    ROMol *mol = sup.next();
    if (mol) {
      ++numMols;
      if (extras) {
        std::unique_ptr<ExplicitBitVect> fp(
            MorganFingerprints::getFingerprintAsBitVect(*mol, 2, 2048));
      }
    }
    delete mol;
  }
  TEST_ASSERT(numMols == expectedResult);
}

void testSDCorrectness() {
  /*
          TEST CORRECTNESS
  */
  std::string rdbase = getenv("RDBASE");
  std::string path;
  unsigned int expectedResult;
  path = "/Code/GraphMol/FileParsers/test_data/NCI_aids_few.sdf";
  expectedResult = 16;
  testSDConcurrent(path, false, true, true, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/esters_end.sdf";
  expectedResult = 6;
  testSDConcurrent(path, true, true, true, 2, 5, 5, expectedResult);

  // strict parsing results in reading 0 out of 2 records
  path = "/Code/GraphMol/FileParsers/test_data/strictLax1.sdf";
  expectedResult = 0;
  testSDConcurrent(path, true, true, true, 2, 5, 5, expectedResult);

  expectedResult = 2;
  testSDConcurrent(path, true, true, false, 2, 5, 5, expectedResult);

  path = "/Code/GraphMol/FileParsers/test_data/O2.sdf";
  expectedResult = 1;
  testSDConcurrent(path, true, true, false, 2, 5, 5, expectedResult);

#ifdef RDK_USE_BOOST_IOSTREAMS

  path = rdbase + "/Regress/Data/mols.1000.sdf.gz";
  std::istream *strm = new gzstream(path);
  expectedResult = 1000;
  testSDConcurrent(strm, true, false, true, true, 2, 5, 5, expectedResult);

#endif
  /*
     TEST PROPERTIES
  */
  testSDProperties();
}

void testPerformance() {
  /*
     TEST PERFORMANCE

     NOTE: Only use this method when you have
     extracted the files znp.50k.smi.gz and chembl26_very_active.sdf.gz in the
     $RDBASE/Regress/Data directory.
  */
  std::string rdbase = getenv("RDBASE");
  unsigned int maxThreadCount =
      std::min(std::thread::hardware_concurrency() + 4, 4u);
#if 1

  {
    std::string path = "/Regress/Data/znp.50k.smi";
    std::string gzpath = "/Regress/Data/znp.50k.smi.gz";
    unsigned int expectedResult = 50000;
    auto start = high_resolution_clock::now();
    // NOTE: have to use path instead of stream, since the tellg()
    //       method, which is used in implementation of the supplier
    // 			 fails for this file.

    std::string delim = " \t";
    int smilesColumn = 0;
    int nameColumn = 1;
    bool titleLine = false;
    bool sanitize = true;
    bool extras = true;
    testSmiOld(rdbase + path, delim, smilesColumn, nameColumn, titleLine,
               sanitize, expectedResult, extras);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration for SmilesMolSupplier: " << duration.count()
              << " (milliseconds) \n";

#if RDK_USE_BOOST_IOSTREAMS
    for (unsigned int i = maxThreadCount; i >= 1; --i) {
      std::istream *strm = new gzstream(rdbase + gzpath);
      start = high_resolution_clock::now();
      bool takeOwnership = true;
      testSmiConcurrent(strm, takeOwnership, delim, smilesColumn, nameColumn,
                        titleLine, sanitize, i, 1000, 100, expectedResult,
                        extras);
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      std::cout << "Duration for testSmiConcurent with " << i
                << "  writer threads: " << duration.count()
                << " (milliseconds) \n";
    }
#endif
  }
#endif

  {
    std::string delim = " \t";
    bool sanitize = true;
    bool removeHs = true;
    bool strictParsing = false;
    bool extras = true;

    std::string path = "/Regress/Data/chembl26_very_active.sdf";
    std::string gzpath = "/Regress/Data/chembl26_very_active.sdf.gz";
    unsigned int expectedResult = 35767;
    auto start = high_resolution_clock::now();
    // NOTE: have to use path instead of stream, since the tellg()
    //       method, which is used in implementation of the supplier
    // 			 fails for this file.
    testSDOld(rdbase + path, sanitize, removeHs, strictParsing, expectedResult,
              extras);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Duration for SDMolSupplier: " << duration.count()
              << " (milliseconds) \n";

#if RDK_USE_BOOST_IOSTREAMS
    for (unsigned int i = maxThreadCount; i >= 1; --i) {
      std::istream *strm = new gzstream(rdbase + gzpath);
      bool takeOwnership = true;
      start = high_resolution_clock::now();
      testSDConcurrent(strm, takeOwnership, sanitize, removeHs, strictParsing,
                       i, 1000, 4000, expectedResult, extras);
      stop = high_resolution_clock::now();
      duration = duration_cast<milliseconds>(stop - start);
      std::cout << "Duration for testSDConcurent with " << i
                << "  writer threads: " << duration.count()
                << " (milliseconds) \n";
    }
#endif
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
    testPerformance();
    BOOST_LOG(rdErrorLog) << "Finished: testPerformance()\n";
    BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";
  */

#endif

  return 0;
}
