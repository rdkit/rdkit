//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_BUILD_THREADSAFE_SSS
#ifndef MULTITHREADED_MOL_SUPPLIER
#define MULTITHREADED_MOL_SUPPLIER

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/ConcurrentQueue.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/RDThreads.h>
#include <RDGeneral/StreamOps.h>

#include <functional>
#include <atomic>
#include <boost/tokenizer.hpp>

#include "FileParsers.h"
#include "MolSupplier.h"

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
namespace v2 {
namespace FileParsers {
class RDKIT_FILEPARSERS_EXPORT MultithreadedMolSupplier : public MolSupplier {
  //! this is an abstract base class to concurrently supply molecules one at a
  //! time
 public:
  using readCallBackFn_t =
      std::function<std::string(const std::string &, unsigned int)>;
  using nextCallBackFn_t =
      std::function<void(RWMol &, const MultithreadedMolSupplier &)>;
  using writeCallBackFn_t =
      std::function<void(RWMol &, const std::string &, unsigned int)>;

  using inputQueue_t =
      ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>;
  using outputQueue_t =
      ConcurrentQueue<std::tuple<RWMol *, std::string, unsigned int>>;

  struct Parameters {
    unsigned int numWriterThreads = 1;
    size_t sizeInputQueue = 5;
    size_t sizeOutputQueue = 5;
  };

  MultithreadedMolSupplier() {}

  // Derived classes MUST have a destructor that calls close
  //  to properly end threads while the instance is alive
  ~MultithreadedMolSupplier() override { close(); }

  //! shut down the supplier
  void close() final;

  //! pop elements from the output queue
  std::unique_ptr<RWMol> next() final;

  //! returns true when all records have been read from the supplier
  bool atEnd() final;

  //! included for the interface, always returns false
  bool getEOFHitOnRead() const { return false; }

  //! returns the record id of the last extracted item
  //! Note: d_LastRecordId = 0, initially therefore the value 0 is returned
  //! if and only if the function is called before extracting the first
  //! record
  unsigned int getLastRecordId() const;

  //! returns the text block for the last extracted item
  std::string getLastItemText() const;

  //! sets the callback to be applied to molecules before they are returned by
  ///! the next() function
  /*!
    \param cb: a function that takes a reference to an RWMol and a const
    reference to the MultithreadedMolSupplier. This can modify the molecule in
    place

   */
  void setNextCallback(nextCallBackFn_t cb) { nextCallback = cb; }

  //! sets the callback to be applied to molecules after they are processed, but
  ///! before they are written to the output queue
  /*!
    \param cb: a function that takes a reference to an RWMol, a const reference
    to the string record, and an unsigned int record id. This can modify the
    molecule in place
  */
  void setWriteCallback(writeCallBackFn_t cb) { writeCallback = cb; }

  //! sets the callback to be applied to input text records before they are
  ///! added to the input queue
  /*!
    \param cb: a function that takes a const reference to the string record and
    an unsigned int record id and returns the modified string record
  */
  void setReadCallback(readCallBackFn_t cb) { readCallback = cb; }

  //! not yet implemented
  void init() final{};

  //! not yet implemented
  void reset() final;

 protected:
  virtual bool getEnd() const = 0;

  void initFromSettings(bool takeOwnership, const Parameters &params);

  //! extracts next record from the input file or stream
  virtual bool extractNextRecord(std::string &record, unsigned int &lineNum,
                                 unsigned int &index) = 0;

  //! processes the record into an RWMol object
  virtual std::unique_ptr<RWMol> processMoleculeRecord(
      const std::string &record, unsigned int lineNum) = 0;

  //!< stores last extracted record id
  std::atomic<unsigned int> d_lastRecordId = 0;

  //!< concurrent input queue
  std::unique_ptr<inputQueue_t> d_inputQueue;

  //!< concurrent output queue
  std::unique_ptr<outputQueue_t> d_outputQueue;

  Parameters d_params;

 private:
  //! Close down any external streams
  void closeStreams();

  //! starts reader and writer threads
  void startThreads();

  //! finalizes the reader and writer threads
  void endThreads();

  //! reads lines from input stream to populate the input queue
  void reader();

  //! parses lines from the input queue converting them to RWMol objects
  //! populating the output queue
  void writer();

  //! disable automatic copy constructors and assignment operators
  //! for this class and its subclasses.  They will likely be
  //! carrying around stream pointers and copying those is a recipe
  //! for disaster.
  MultithreadedMolSupplier(const MultithreadedMolSupplier &) = delete;

  MultithreadedMolSupplier &operator=(const MultithreadedMolSupplier &) =
      delete;

  std::atomic<bool> df_started = false;
  std::atomic<bool> df_forceStop = false;

  std::mutex d_threadCounterMutex;
  std::atomic<unsigned int> d_threadCounter{1};  //!< thread counter
  std::vector<std::thread> d_writerThreads;      //!< vector writer threads
  std::thread d_readerThread;                    //!< single reader thread

  std::string d_lastItemText;  //!< stores last extracted record

  readCallBackFn_t readCallback = nullptr;
  nextCallBackFn_t nextCallback = nullptr;
  writeCallBackFn_t writeCallback = nullptr;
};
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
#endif
