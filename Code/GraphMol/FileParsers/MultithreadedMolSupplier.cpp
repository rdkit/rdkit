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

namespace RDKit {

//! method for thread-safe printing, only for debugging
struct PrintThread : public std::stringstream {
  static inline std::mutex cout_mutex;
  ~PrintThread() {
    std::lock_guard<std::mutex> l{cout_mutex};
    std::cout << rdbuf();
  }
};

void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum;
  while (extractNextRecord(record, lineNum)) {
    auto r = std::tuple<std::string, int>{record, lineNum};
    PrintThread{}
        << "Reader id: 0, Pushing elements into the input queue, lineNum: "
        << lineNum << std::endl;
    d_inputQueue->push(r);
  }
  PrintThread{} << "Exit Reader, done reading" << std::endl;
}

void MultithreadedMolSupplier::writer(size_t id) {
  std::tuple<std::string, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    ROMol* mol = processMoleculeRecord(std::get<0>(r), std::get<1>(r));
    PrintThread{} << "---------------------------------------------------------"
                     "-------------"
                  << "Writer id: " << id
                  << ", Pushing elements into the output queue, lineNum: "
                  << std::get<1>(r) << std::endl;
    d_outputQueue->push(mol);
  }
  PrintThread{} << "-----------------------------------------------------------"
                   "-----------"
                << "Exit Writer: " << id << ", done processing molecules"
                << std::endl;
}

ROMol* MultithreadedMolSupplier::next() {
  ROMol* mol;
  if (d_outputQueue->pop(mol)) {
    PrintThread{} << "---------------------------------------------------------"
                     "-------------"
                  << "Next method, popped mol: " << MolToSmiles(*mol)
                  << std::endl;

    return mol;
  }
  return nullptr;
}

void MultithreadedMolSupplier::joinReaderAndWriter() {
  readerThread.join();
  d_inputQueue->setDone();
  for (auto& thread : writerThreads) {
    thread.join();
  }
  d_outputQueue->setDone();
}

void MultithreadedMolSupplier::startThreads() {
  //! run the reader function in a seperate thread
  readerThread = std::thread(&MultithreadedMolSupplier::reader, this);
  //! run the writer function in seperate threads
  for (int i = 0; i < d_numWriterThreads; i++) {
    writerThreads.push_back(
        std::thread(&MultithreadedMolSupplier::writer, this, i));
  }
}

bool MultithreadedMolSupplier::atEnd() { return d_outputQueue->getDone(); }

void MultithreadedMolSupplier::reset() { ; }

}  // namespace RDKit
