//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_THREADSAFE_SSS
#ifndef MULTITHREADED_MOL_SUPPLIER
#define MULTITHREADED_MOL_SUPPLIER

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/ConcurrentQueue.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/StreamOps.h>

#include <boost/tokenizer.hpp>
#include <cstdlib>
#include <sstream>

#include "FileParsers.h"
#include "MolSupplier.h"

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
class RDKIT_FILEPARSERS_EXPORT MultithreadedMolSupplier : public MolSupplier {
 public:
  MultithreadedMolSupplier();
  MultithreadedMolSupplier(std::istream *inStream, bool takeOwnership,
                           int numWriterThreads, size_t sizeInputQueue,
                           size_t sizeOutputQueue);
  MultithreadedMolSupplier(std::string fileName, bool takeOwnership,
                           int numWriterThreads, size_t sizeInputQueue,
                           size_t sizeOutputQueue);
  //! intialize data memebers
  void init();

  //! reads lines from input stream to populate the input queue
  void inputProducer();

  //! parses lines from the input queue converting them to ROMol objects
  //! populating the output queue
  void inputConsumer();

  virtual void reset();
  virtual ROMol *next();
  virtual bool atEnd();
  virtual bool extractNextRecord(std::string &record) = 0;
  virtual ROMol *processMoleculeRecord(const std::string &record) = 0;

 private:
  const int d_numReaderThread = 1;  // fix number of reader threads to 1
  int d_numWriterThreads;           // number of writer threads
  size_t d_sizeInputQueue;          // size of input queue
  size_t d_sizeOutputQueue;         // size of output queue
  Concurrent Queue<std::string> *d_inputQueue;  // concurrent input queue
  Concurrent Queue<ROMol> *d_outputQueue;       // concurrent output queue
};
}  // namespace RDKit
#endif
#endif
