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
namespace RDKit {

void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum;
  while (extractNextRecord(record, lineNum)) {
    auto r = std::tuple<std::string, int>{record, lineNum};
    d_inputQueue->push(r);
  }
  d_inputQueue->setDone();
}

void MultithreadedMolSupplier::writer() {
  std::tuple<std::string, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    ROMol* mol = processMoleculeRecord(std::get<0>(r), std::get<1>(r));
    d_outputQueue->push(mol);
  }
  if (threadCounter != d_numWriterThreads) {
    ++threadCounter;
  } else {
    d_outputQueue->setDone();
  }
}

ROMol* MultithreadedMolSupplier::next() {
  ROMol* mol;
  if (d_outputQueue->pop(mol)) {
    return mol;
  }
  return nullptr;
}

void MultithreadedMolSupplier::endThreads() {
  readerThread.join();
  for (auto& thread : writerThreads) {
    thread.join();
  }
}

void MultithreadedMolSupplier::startThreads() {
  //! run the reader function in a seperate thread
  readerThread = std::thread(&MultithreadedMolSupplier::reader, this);
  //! run the writer function in seperate threads
  for (int i = 0; i < d_numWriterThreads; i++) {
    writerThreads.push_back(
        std::thread(&MultithreadedMolSupplier::writer, this));
  }
}

bool MultithreadedMolSupplier::atEnd() {
  return (d_outputQueue->isEmpty() && d_outputQueue->getDone());
}

void MultithreadedMolSupplier::reset() { ; }

}  // namespace RDKit
