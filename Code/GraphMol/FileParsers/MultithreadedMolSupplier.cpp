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
    PrintThread{}
        << "Reader id: 0, Pushing elements into the input queue, lineNum: "
        << lineNum << std::endl;
    d_inputQueue->push(r);
  }
  PrintThread{} << "Exit Input Producer, done reading" << std::endl;
}

void MultithreadedMolSupplier::inputConsumer(size_t id) {
  std::tuple<std::string, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    ROMol* mol = processMoleculeRecord(std::get<0>(r), std::get<1>(r));
    std::shared_ptr<ROMol> shared_mol(mol);
    PrintThread{} << "Consumer id: " << id
                  << ", Pushing elements into the output queue, lineNum: "
                  << std::get<1>(r) << std::endl;
    //! possible deadlock here if the output queue is full
    d_outputQueue->push(shared_mol);
  }
  PrintThread{} << "Exit Input Consumer: " << id
                << ", done processing molecules" << std::endl;
}

ROMol* MultithreadedMolSupplier::next() {
  std::cout << "Next method is being called" << std::endl;
  std::shared_ptr<ROMol> mol(nullptr);
  if (d_outputQueue->pop(mol)) {
    PrintThread{} << "Popping elements from the output queue"
                  << MolToSmiles(*mol.get()) << std::endl;
    return mol.get();
  }
  PrintThread{} << "We have a nullptr" << std::endl;
  return nullptr;
}

void MultithreadedMolSupplier::startThreads() {
  std::vector<std::thread> writers(d_numWriterThreads);

  std::thread reader(&MultithreadedMolSupplier::inputProducer, this);
  for (int i = 0; i < d_numWriterThreads; i++) {
    writers[i] = std::thread(&MultithreadedMolSupplier::inputConsumer, this, i);
  }
  reader.join();
  //! done reading elements from the file
  d_inputQueue->setDone();

  std::for_each(writers.begin(), writers.end(),
                std::mem_fn(&std::thread::join));
  //! done wrting elements from the input queue
  d_outputQueue->setDone();
}

bool MultithreadedMolSupplier::atEnd() { return d_outputQueue->isEmpty(); }

void MultithreadedMolSupplier::reset() { ; }

}  // namespace RDKit
