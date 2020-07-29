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
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "MultithreadedMolSupplier.h"
#include "MultithreadedSmilesMolSupplier.h"

namespace io = boost::iostreams;
using namespace RDKit;

void testSmiGeneral(int numWriterThreads, size_t sizeInputQueue,
                    size_t sizeOutputQueue) {
  std::string mname;
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  unsigned int nMols = 0;
  MultithreadedSmilesMolSupplier sup(fname, ",", 1, 0, true, true,
                                     numWriterThreads, sizeInputQueue,
                                     sizeOutputQueue);
  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol != nullptr) {
      ++nMols;
    }
    delete mol;
  }
  TEST_ASSERT(nMols == 10);
}

void testSmi() {
  //! try to construct a multithreaded smiles mol supplier with different
  //! parameters
  testSmiGeneral(2, 5, 5);
  testSmiGeneral(5, 5, 5);
  testSmiGeneral(2, 5, 10);
  testSmiGeneral(2, 10, 5);

  //! an extreme case
  testSmiGeneral(2, 1, 1);
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSmi();
  BOOST_LOG(rdErrorLog) << "Finished: testSmi()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  return 0;
}
