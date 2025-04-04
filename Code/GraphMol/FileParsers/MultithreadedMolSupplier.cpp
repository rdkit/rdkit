#ifdef RDK_BUILD_THREADSAFE_SSS
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

#include <RDGeneral/RDLog.h>

namespace RDKit {

namespace v2 {
namespace FileParsers {

MultithreadedMolSupplier::~MultithreadedMolSupplier() {
  endThreads();
  // destroy all objects in the input queue
  d_inputQueue->clear();
  if (df_started) {
    std::tuple<RWMol *, std::string, unsigned int> r;
    while (d_outputQueue->pop(r)) {
      RWMol *m = std::get<0>(r);
      delete m;
    }
  }
  // destroy all objects in the output queue
  d_outputQueue->clear();
}

void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum, index;
  while (extractNextRecord(record, lineNum, index)) {
    if (readCallback) {
      try {
        record = readCallback(record, index);
      } catch (std::exception &e) {
        BOOST_LOG(rdErrorLog)
            << "Read callback exception: " << e.what() << std::endl;
      }
    }
    auto r = std::make_tuple(record, lineNum, index);
    d_inputQueue->push(r);
  }
  d_inputQueue->setDone();
}

void MultithreadedMolSupplier::writer() {
  std::tuple<std::string, unsigned int, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    try {
      std::unique_ptr<RWMol> mol(
          processMoleculeRecord(std::get<0>(r), std::get<1>(r)));
      if (mol && writeCallback) {
        writeCallback(*mol, std::get<0>(r), std::get<2>(r));
      }
      auto temp = std::tuple<RWMol *, std::string, unsigned int>{
          mol.release(), std::get<0>(r), std::get<2>(r)};
      d_outputQueue->push(temp);
    } catch (...) {
      // fill the queue wih a null value
      auto nullValue = std::tuple<RWMol *, std::string, unsigned int>{
          nullptr, std::get<0>(r), std::get<2>(r)};
      d_outputQueue->push(nullValue);
    }
  }
  if (d_threadCounter != d_params.numWriterThreads) {
    ++d_threadCounter;
  } else {
    d_outputQueue->setDone();
  }
}

std::unique_ptr<RWMol> MultithreadedMolSupplier::next() {
  if (!df_started) {
    df_started = true;
    startThreads();
  }
  std::tuple<RWMol *, std::string, unsigned int> r;
  if (d_outputQueue->pop(r)) {
    d_lastItemText = std::get<1>(r);
    d_lastRecordId = std::get<2>(r);
    std::unique_ptr<RWMol> res{std::get<0>(r)};
    if (res && nextCallback) {
      try {
        nextCallback(*res, *this);
      } catch (...) {
        // Ignore exception and proceed with mol as is.
      }
    }
    return res;
  }
  return nullptr;
}

void MultithreadedMolSupplier::endThreads() {
  if (!df_started) {
    return;
  }
  d_readerThread.join();
  for (auto &thread : d_writerThreads) {
    thread.join();
  }
}

void MultithreadedMolSupplier::startThreads() {
  // run the reader function in a seperate thread
  d_readerThread = std::thread(&MultithreadedMolSupplier::reader, this);
  // run the writer function in seperate threads
  for (unsigned int i = 0; i < d_params.numWriterThreads; i++) {
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
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
