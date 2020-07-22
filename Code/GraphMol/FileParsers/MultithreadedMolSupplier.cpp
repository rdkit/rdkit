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

MultithreadedMolSupplier::MultithreadedMolSupplier() { init(); }

MultithreadedMolSupplier::MultithreadedMolSupplier(std::istream* inStream,
                                                   bool takeOwnership,
                                                   int numWriterThreads,
                                                   size_t sizeInputQueue,
                                                   size_t sizeOutputQueue) {
  dp_inStream = inStream;
  df_owner = takeOwnership;
  if (numWriterThreads == -1) {
    d_numWriterThreads = (int)getNumThreadsToUse(numWriterThreads);
  } else {
    d_numWriterThreads =
        std::min(numWriterThreads, (int)getNumThreadsToUse(numWriterThreads));
  }
  d_sizeInputQueue = sizeInputQueue;
  d_sizeOutputQueue = sizeOutputQueue;
  d_inputQueue = new ConcurrentQueue<Record>(d_sizeInputQueue);
  d_outputQueue =
      new ConcurrentQueue<std::shared_ptr<ROMol>>(d_sizeOutputQueue);
}

MultithreadedMolSupplier::MultithreadedMolSupplier(const std::string fileName,
                                                   bool takeOwnership,
                                                   int numWriterThreads,
                                                   size_t sizeInputQueue,
                                                   size_t sizeOutputQueue) {
  dp_inStream = openAndCheckStream(fileName);
  CHECK_INVARIANT(dp_inStream, "bad instream");
  CHECK_INVARIANT(!(dp_inStream->eof()), "early EOF");
  df_owner = takeOwnership;
  if (numWriterThreads == -1) {
    d_numWriterThreads = (int)getNumThreadsToUse(numWriterThreads);
  } else {
    d_numWriterThreads =
        std::min(numWriterThreads, (int)getNumThreadsToUse(numWriterThreads));
  }

  d_sizeInputQueue = sizeInputQueue;
  d_sizeOutputQueue = sizeOutputQueue;
  d_inputQueue = new ConcurrentQueue<Record>(d_sizeInputQueue);
  d_outputQueue =
      new ConcurrentQueue<std::shared_ptr<ROMol>>(d_sizeOutputQueue);
}

void MultithreadedMolSupplier::init() {
  dp_inStream = nullptr;
  df_owner = true;
  d_numWriterThreads = 1;
  d_sizeInputQueue = 10;
  d_sizeOutputQueue = 10;
  d_inputQueue = new ConcurrentQueue<Record>(d_sizeInputQueue);
  d_outputQueue =
      new ConcurrentQueue<std::shared_ptr<ROMol>>(d_sizeOutputQueue);
}

void MultithreadedMolSupplier::inputProducer() {
  std::string record;
  unsigned int lineNum;
  while (extractNextRecord(record, lineNum)) {
    struct Record r;
    r.record = record;
    r.lineNum = lineNum;
    d_inputQueue->push(r);
  }
}

void MultithreadedMolSupplier::inputConsumer() {
  struct Record r while (d_inputQueue->pop(r)) {
    ROMol* mol = processMoleculeRecord(r.record, r.lineNum);
    std::shared_ptr<ROMol> shared_mol(mol);
    d_outputQueue->push(shared_mol);
  }
}

ROMol* MultithreadedMolSupplier::next() {
  std::shared_ptr<ROMol> mol(nullptr);
  if (d_outputQueue->pop(mol)) {
    return mol;
  }
  return nullptr;
}

void MultithreadedMolSupplier::startThreads() {
  std::vector<std::thread> writers(d_numWriterThreads);
  std::thread reader = std::thread(inputProducer);
  for (int i = 0; i < d_numWriterThreads; i++) {
    writers[i] = std::thread(inputConsumers);
  }
  reader.join();
  reader.setDone();
  std::for_each(writers.begin(), writers.end(),
                std::mem_fn(&std::thread::join));
  writers.setDone();
}
