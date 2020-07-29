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

void func(MultithreadedSmilesMolSupplier& sup, unsigned int& i) {
  while (!sup.atEnd()) {
    ROMol* mol = sup.next();
    if (mol == nullptr) {
      break;
    }
    delete mol;
    i++;
  }
}

void testSmi() {
  std::string mname;
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/FileParsers/test_data/fewSmi.csv";
  unsigned int i;
  //! create the mol supplier with the usual syntax
  MultithreadedSmilesMolSupplier sup(fname, ",", 1, 0, true);

  sup.startThreads();
  //! run your method in a seperate thread
  std::thread th(func, std::ref(sup), std::ref(i));
  sup.joinReaderAndWriter();
  //! join the thread after the reader and writer method
  th.join();

  std::cout << "i = " << i << std::endl;
  TEST_ASSERT(i == 10);
}

int main() {
  RDLog::InitLogs();

  BOOST_LOG(rdErrorLog) << "\n-----------------------------------------\n";
  testSmi();
  BOOST_LOG(rdErrorLog) << "Finished: testSmi()\n";
  BOOST_LOG(rdErrorLog) << "-----------------------------------------\n\n";

  return 0;
}
