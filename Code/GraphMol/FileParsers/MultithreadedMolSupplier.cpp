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
MultithreadedMolSupplier::~MultithreadedMolSupplier() {
  endThreads();
  // rethrow the exceptions once the threads have been executed
  for (auto const& e : d_exceptions) {
    try {
      if (e != nullptr) {
        std::rethrow_exception(e);
      }
    } catch (std::exception const& ex) {
      BOOST_LOG(rdErrorLog) << "ERROR: " << ex.what() << "\n";
    }
  }
  // exception_ptr is a smart shared pointer so we only
  // need to clear the vector
  d_exceptions.clear();
  // destroy all objects in the input queue
  d_inputQueue->clear();
  // delete the pointer to the input queue
  delete d_inputQueue;
  std::tuple<ROMol*, std::string, unsigned int> r;
  while (d_outputQueue->pop(r)) {
    ROMol* m = std::get<0>(r);
    delete m;
  }
  // destroy all objects in the output queue
  d_outputQueue->clear();
  // delete the pointer to the output queue
  delete d_outputQueue;
}

void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum, index;
  while (extractNextRecord(record, lineNum, index)) {
    auto r = std::tuple<std::string, unsigned int, unsigned int>{
        record, lineNum, index};
    d_inputQueue->push(r);
  }
  d_inputQueue->setDone();
}

void MultithreadedMolSupplier::writer() {
  std::tuple<std::string, unsigned int, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    try {
      ROMol* mol = processMoleculeRecord(std::get<0>(r), std::get<1>(r));
      auto temp = std::tuple<ROMol*, std::string, unsigned int>{
          mol, std::get<0>(r), std::get<2>(r)};
      d_outputQueue->push(temp);
    } catch (...) {
      std::lock_guard<std::mutex> lock(d_mutex);
      d_exceptions.push_back(std::current_exception());
      // fill the queue wih a null value
      auto nullValue = std::tuple<ROMol*, std::string, unsigned int>{
          nullptr, std::get<0>(r), std::get<2>(r)};
      d_outputQueue->push(nullValue);
    }
  }

  if (d_threadCounter != d_numWriterThreads) {
    ++d_threadCounter;
  } else {
    d_outputQueue->setDone();
  }
}

ROMol* MultithreadedMolSupplier::next() {
  std::tuple<ROMol*, std::string, unsigned int> r;
  if (d_outputQueue->pop(r)) {
    ROMol* mol = std::get<0>(r);
    d_lastItemText = std::get<1>(r);
    d_lastRecordId = std::get<2>(r);
    return mol;
  }
  return nullptr;
}

void MultithreadedMolSupplier::endThreads() {
  d_readerThread.join();
  for (auto& thread : d_writerThreads) {
    thread.join();
  }
}

void MultithreadedMolSupplier::startThreads() {
  // run the reader function in a seperate thread
  d_readerThread = std::thread(&MultithreadedMolSupplier::reader, this);
  // run the writer function in seperate threads
  for (unsigned int i = 0; i < d_numWriterThreads; i++) {
    d_writerThreads.emplace_back(
        std::thread(&MultithreadedMolSupplier::writer, this));
  }
}

bool MultithreadedMolSupplier::atEnd() {
  return (d_outputQueue->isEmpty() && d_outputQueue->getDone());
}

unsigned int MultithreadedMolSupplier::getLastRecordId() const {
  return d_lastRecordId;
}

std::string MultithreadedMolSupplier::getLastItemText() const {
  return d_lastItemText;
}

void MultithreadedMolSupplier::reset() {
  UNDER_CONSTRUCTION("reset() not supported for MultithreadedMolSupplier();");
}

}  // namespace RDKit
