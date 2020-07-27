//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MultithreadedMolSupplier.h"

#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <chrono>
#include <future>
#include <iomanip>
#include <mutex>
#include <sstream>

namespace RDKit {

//! method for thread-safe printing, only for debugging
struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

void MultithreadedMolSupplier::inputProducer() {
  std::string record;
  unsigned int lineNum;
  while (extractNextRecord(record, lineNum)) {
    auto r = std::tuple<std::string, unsigned int>{record, lineNum};
    /*
PrintThread{}
<< "Reader id: 0, Pushing elements into the input queue, lineNum: "
<< lineNum << std::endl;
    */
    d_inputQueue->push(r);
  }
  // PrintThread{} << "Exit Input Producer, done reading" << std::endl;
}

void MultithreadedMolSupplier::inputConsumer(size_t id) {
  std::tuple<std::string, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    ROMol* mol = processMoleculeRecord(std::get<0>(r), std::get<1>(r));
    /*
                PrintThread{} << "Input Consumer id: " << id
                  << ", Pushing elements into the output queue, lineNum: "
                  << std::get<1>(r) << std::endl;
                */
    d_outputQueue->push(mol);
  }
  // PrintThread{} << "Exit Input Consumer: " << id
  //             << ", done processing molecules" << std::endl;
}

ROMol* MultithreadedMolSupplier::next() {
  ROMol* mol;
  // std::cout << "Size of output queue: " << d_outputQueue->size() <<
  // std::endl;
  if (d_outputQueue->pop(mol)) {
    //  PrintThread{} << "Popping elements from the output queue: "
    //              << MolToSmiles(*mol) << std::endl;
    return mol;
  }
  // PrintThread{} << "We have a nullptr" << std::endl;
  return nullptr;
}

void MultithreadedMolSupplier::joinReaderAndWriters() {
  reader.join();
  d_inputQueue->setDone();
  std::for_each(writers.begin(), writers.end(),
                std::mem_fn(&std::thread::join));
  d_outputQueue->setDone();
}

void MultithreadedMolSupplier::startThreads() {
  //! start reader threads
  reader = std::thread(&MultithreadedMolSupplier::inputProducer, this);
  //! start writer threads
  for (int i = 0; i < d_numWriterThreads; i++) {
    writers.push_back(
        std::thread(&MultithreadedMolSupplier::inputConsumer, this, i));
  }
}

bool MultithreadedMolSupplier::atEnd() { return d_outputQueue->getDone(); }

void MultithreadedMolSupplier::reset() { ; }

}  // namespace RDKit
