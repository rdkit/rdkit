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
#include <RDGeneral/RDThreads.h>
#include <RDGeneral/StreamOps.h>

#include <atomic>
#include <boost/tokenizer.hpp>

#include "FileParsers.h"
#include "MolSupplier.h"

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
class RDKIT_FILEPARSERS_EXPORT MultithreadedMolSupplier : public MolSupplier {
  //! this is an abstract base class to concurrently supply molecules one at a
  //! time
 public:
  MultithreadedMolSupplier(){};
  virtual ~MultithreadedMolSupplier();
  //! pop elements from the output queue
  ROMol *next();
  //! returns true when all records have been read from the supplier
  bool atEnd();

  //! included for the interface, always returns false
  bool getEOFHitOnRead() const { return false; }

  //! returns the record id of the last extracted item
  //! Note: d_LastRecordId = 0, initially therefore the value 0 is returned
  //! if and only if the function is called before extracting the first
  //! record
  unsigned int getLastRecordId() const;
  //! returns the text block for the last extracted item
  std::string getLastItemText() const;

 protected:
  //! starts reader and writer threads
  void startThreads();

 private:
  //! reads lines from input stream to populate the input queue
  void reader();
  //! parses lines from the input queue converting them to ROMol objects
  //! populating the output queue
  void writer();
  //! finalizes the reader and writer threads
  void endThreads();
  //! disable automatic copy constructors and assignment operators
  //! for this class and its subclasses.  They will likely be
  //! carrying around stream pointers and copying those is a recipe
  //! for disaster.
  MultithreadedMolSupplier(const MultithreadedMolSupplier &);
  MultithreadedMolSupplier &operator=(const MultithreadedMolSupplier &);
  //! not yet implemented
  virtual void reset();
  virtual void init() = 0;
  virtual bool getEnd() const = 0;
  //! extracts next record from the input file or stream
  virtual bool extractNextRecord(std::string &record, unsigned int &lineNum,
                                 unsigned int &index) = 0;
  //! processes the record into an ROMol object
  virtual ROMol *processMoleculeRecord(const std::string &record,
                                       unsigned int lineNum) = 0;

 private:
  std::atomic<unsigned int> d_threadCounter{1};  //! thread counter
  std::vector<std::thread> d_writerThreads;      //! vector writer threads
  std::thread d_readerThread;                    //! single reader thread

 protected:
  unsigned int d_lastRecordId = 0;           //! stores last extracted record id
  std::string d_lastItemText;                //! stores last extracted record
  const unsigned int d_numReaderThread = 1;  //! number of reader thread
  unsigned int d_numWriterThreads;           //! number of writer threads
  size_t d_sizeInputQueue;                   //! size of input queue
  size_t d_sizeOutputQueue;                  //! size of output queue

  ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>
      *d_inputQueue;  //! concurrent input queue
  ConcurrentQueue<std::tuple<ROMol *, std::string, unsigned int>>
      *d_outputQueue;  //! concurrent output queue
};
}  // namespace RDKit
#endif
#endif
