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

void MultithreadedMolSupplier::close() {
  df_forceStop = true;
  d_outputQueue->setDone();

  if (df_started) {
    // Clear the queues until they are empty
    //  d_inputQueue->clear is not thread-safe
    std::tuple<std::string, unsigned int, unsigned int> r;
    while (d_inputQueue->pop(r)) {
    }
    // clear the output queues, they might be full
    //  and blocking the writer threads, note
    //  that while ending threads the writers may
    //  put a few more items back in the queue
    std::tuple<RWMol *, std::string, unsigned int> mol_r;
    while (d_outputQueue->pop(mol_r)) {
      delete std::get<0>(mol_r);
    }
  }

  endThreads();

  // notify the queue again that it is done in case
  //  anyone is waiting on it
  d_outputQueue->setDone();

  // destroy all objects in the input and output queues
  //  and anything missed put in the queues while
  //  the threads were endings
  if (df_started) {
    d_inputQueue->clear();
  }

  if (d_outputQueue) {
    // destroy all objects in the output queue
    std::tuple<RWMol *, std::string, unsigned int> r;
    while (d_outputQueue->pop(r)) {
      delete std::get<0>(r);
    }
  }

  // close external streams if any
  //  destructors are called child to parent, however the threads
  //  need to be ended before shutting down streams, so override this
  //  in the child class.
  closeStreams();
  df_started = false;
}

void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum, index;
  while (!df_forceStop && extractNextRecord(record, lineNum, index)) {
    if (readCallback) {
      try {
        record = readCallback(record, index);
      } catch (std::exception &e) {
        BOOST_LOG(rdErrorLog)
            << "Read callback exception: " << e.what() << std::endl;
      }
    }
    auto r = std::make_tuple(record, lineNum, index);
    if (!df_forceStop) {
      d_inputQueue->push(r);
    }
  }
  d_inputQueue->setDone();
}

void MultithreadedMolSupplier::writer() {
  std::tuple<std::string, unsigned int, unsigned int> r;
  while (!df_forceStop && d_inputQueue->pop(r)) {
    try {
      std::unique_ptr<RWMol> mol(
          processMoleculeRecord(std::get<0>(r), std::get<1>(r)));
      if (!df_forceStop && mol && writeCallback) {
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

  // we need a lock here otherwise two threads
  //  can increment d_threadCounter even though it's
  //  atomic.
  d_threadCounterMutex.lock();
  if (d_threadCounter < d_params.numWriterThreads) {
    ++d_threadCounter;
    d_threadCounterMutex.unlock();
  } else {
    // Here we need to unlock the threadCounterMutex before we setDone on the
    //  outputQueue.  This causes a notification to the queue which may actually
    //  have elements in it.  This notification may unblock the queue which
    //  allows waiting threads to get their last attempt at adding to it
    //  which will end up here and deadlock.
    d_threadCounterMutex.unlock();
    d_outputQueue->setDone();
  }
}

std::unique_ptr<RWMol> MultithreadedMolSupplier::next() {
  if (!df_started) {
    df_started = true;
    startThreads();
  }
  std::tuple<RWMol *, std::string, unsigned int> r;
  while (!d_outputQueue->pop(r)) {
    if (df_forceStop) {
      throw std::runtime_error(
          "Multhreded supplier closed while waiting for a mol");
    } else if (d_outputQueue->getDone()) {
      df_eofHitOnRead = true;
      return nullptr;
    }
  }
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

// this calls joins on the reader and writer threads
//  and waits until completion.  To actually force a stop
//  call close which handles the input and output queues
void MultithreadedMolSupplier::endThreads() {
  if (!df_started) {
    return;
  }

  // stop the writers before stopping the readers
  //  otherwise there might be a deadlock
  for (auto &thread : d_writerThreads) {
    thread.join();
  }
  d_readerThread.join();
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
  return getEOFHitOnRead() ||
         (d_inputQueue->isEmpty() && d_inputQueue->getDone() &&
          d_outputQueue->isEmpty() && d_outputQueue->getDone());
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
