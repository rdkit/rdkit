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
  struct Parameters {
    unsigned int numWriterThreads = 1;
    size_t sizeInputQueue = 5;
    size_t sizeOutputQueue = 5;
  };

  MultithreadedMolSupplier() {}
  ~MultithreadedMolSupplier() override;
  //! pop elements from the output queue
  std::unique_ptr<RWMol> next() override;

  //! returns true when all records have been read from the supplier
  bool atEnd() override;

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
  template <typename T>
  void setNextCallback(T cb) {
    nextCallback = cb;
  }
  //! sets the callback to be applied to molecules after they are processed, but
  ///! before they are written to the output queue
  /*!
    \param cb: a function that takes a reference to an RWMol, a const reference
    to the string record, and an unsigned int record id. This can modify the
    molecule in place
  */
  template <typename T>
  void setWriteCallback(T cb) {
    writeCallback = cb;
  }
  //! sets the callback to be applied to input text records before they are
  ///! added to the input queue
  /*!
    \param cb: a function that takes a const reference to the string record and
    an unsigned int record id and returns the modified string record
  */
  template <typename T>
  void setReadCallback(T cb) {
    readCallback = cb;
  }

 protected:
  //! starts reader and writer threads
  void startThreads();

 private:
  //! reads lines from input stream to populate the input queue
  void reader();
  //! parses lines from the input queue converting them to RWMol objects
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
  void reset() override;
  void init() override = 0;
  virtual bool getEnd() const = 0;
  //! extracts next record from the input file or stream
  virtual bool extractNextRecord(std::string &record, unsigned int &lineNum,
                                 unsigned int &index) = 0;
  //! processes the record into an RWMol object
  virtual RWMol *processMoleculeRecord(const std::string &record,
                                       unsigned int lineNum) = 0;

  std::atomic<unsigned int> d_threadCounter{1};  //!< thread counter
  std::vector<std::thread> d_writerThreads;      //!< vector writer threads
  std::thread d_readerThread;                    //!< single reader thread

 protected:
  std::atomic<bool> df_started = false;
  std::atomic<unsigned int> d_lastRecordId =
      0;                       //!< stores last extracted record id
  std::string d_lastItemText;  //!< stores last extracted record
  const unsigned int d_numReaderThread = 1;  //!< number of reader thread

  std::unique_ptr<
      ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>>
      d_inputQueue;  //!< concurrent input queue
  std::unique_ptr<
      ConcurrentQueue<std::tuple<RWMol *, std::string, unsigned int>>>
      d_outputQueue;  //!< concurrent output queue
  Parameters d_params;
  std::function<void(RWMol &, const MultithreadedMolSupplier &)> nextCallback =
      nullptr;
  std::function<void(RWMol &, const std::string &, unsigned int)>
      writeCallback = nullptr;
  std::function<std::string(const std::string &, unsigned int)> readCallback =
      nullptr;
};
}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
class RDKIT_FILEPARSERS_EXPORT MultithreadedMolSupplier : public MolSupplier {
  //! this is an abstract base class to concurrently supply molecules one at a
  //! time
 public:
  using ContainedType = v2::FileParsers::MultithreadedMolSupplier;
  MultithreadedMolSupplier() {}

  //! included for the interface, always returns false
  bool getEOFHitOnRead() const {
    if (dp_supplier) {
      return static_cast<ContainedType *>(dp_supplier.get())->getEOFHitOnRead();
    }
    return false;
  }

  //! returns the record id of the last extracted item
  //! Note: d_LastRecordId = 0, initially therefore the value 0 is returned
  //! if and only if the function is called before extracting the first
  //! record
  unsigned int getLastRecordId() const {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->getLastRecordId();
  }
  //! returns the text block for the last extracted item
  std::string getLastItemText() const {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->getLastItemText();
  }
};
}  // namespace v1
}  // namespace RDKit
#endif
#endif
