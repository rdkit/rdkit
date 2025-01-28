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
#ifndef MULTITHREADED_SMILES_MOL_SUPPLIER
#define MULTITHREADED_SMILES_MOL_SUPPLIER
#include "MultithreadedMolSupplier.h"
namespace RDKit {
namespace v2 {
namespace FileParsers {
//! This class is still a bit experimental and the public API may change
//! in future releases.
class RDKIT_FILEPARSERS_EXPORT MultithreadedSmilesMolSupplier
    : public MultithreadedMolSupplier {
 public:
  explicit MultithreadedSmilesMolSupplier(
      const std::string &fileName, const Parameters &params = Parameters(),
      const SmilesMolSupplierParams &parseParams = SmilesMolSupplierParams());

  explicit MultithreadedSmilesMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const Parameters &params = Parameters(),
      const SmilesMolSupplierParams &parseParams = SmilesMolSupplierParams());
  MultithreadedSmilesMolSupplier();
  ~MultithreadedSmilesMolSupplier() override;

  void init() override {}
  //! returns df_end
  bool getEnd() const override;
  //! reads and processes the title line
  void processTitleLine();
  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum,
                         unsigned int &index) override;
  //! parses the record and returns the resulting molecule
  RWMol *processMoleculeRecord(const std::string &record,
                               unsigned int lineNum) override;

 private:
  void initFromSettings(
      bool takeOwnership, const Parameters &params,
      const SmilesMolSupplierParams &parseParams = SmilesMolSupplierParams());

 private:
  bool df_end = false;                 //!< have we reached the end of the file?
  int d_line = 0;                      //!< line number we are currently on
  STR_VECT d_props;                    //!< vector of property names
  unsigned int d_currentRecordId = 1;  //!< current record id
  SmilesMolSupplierParams d_parseParams;
};
}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
class RDKIT_FILEPARSERS_EXPORT MultithreadedSmilesMolSupplier
    : public MolSupplier {
  //! this is an abstract base class to concurrently supply molecules one at a
  //! time
 public:
  using ContainedType = v2::FileParsers::MultithreadedSmilesMolSupplier;
  MultithreadedSmilesMolSupplier() {}
  explicit MultithreadedSmilesMolSupplier(
      const std::string &fileName, const std::string &delimiter = " \t",
      int smilesColumn = 0, int nameColumn = 1, bool titleLine = true,
      bool sanitize = true, unsigned int numWriterThreads = 1,
      size_t sizeInputQueue = 5, size_t sizeOutputQueue = 5) {
    v2::FileParsers::MultithreadedSmilesMolSupplier::Parameters params;
    params.numWriterThreads = numWriterThreads;
    params.sizeInputQueue = sizeInputQueue;
    params.sizeOutputQueue = sizeOutputQueue;
    v2::FileParsers::SmilesMolSupplierParams parseParams;
    parseParams.delimiter = delimiter;
    parseParams.smilesColumn = smilesColumn;
    parseParams.nameColumn = nameColumn;
    parseParams.titleLine = titleLine;
    parseParams.parseParameters.sanitize = sanitize;

    dp_supplier.reset(new v2::FileParsers::MultithreadedSmilesMolSupplier(
        fileName, params, parseParams));
  }

  explicit MultithreadedSmilesMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const std::string &delimiter = " \t", int smilesColumn = 0,
      int nameColumn = 1, bool titleLine = true, bool sanitize = true,
      unsigned int numWriterThreads = 1, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5) {
    v2::FileParsers::MultithreadedSmilesMolSupplier::Parameters params;
    params.numWriterThreads = numWriterThreads;
    params.sizeInputQueue = sizeInputQueue;
    params.sizeOutputQueue = sizeOutputQueue;
    v2::FileParsers::SmilesMolSupplierParams parseParams;
    parseParams.delimiter = delimiter;
    parseParams.smilesColumn = smilesColumn;
    parseParams.nameColumn = nameColumn;
    parseParams.titleLine = titleLine;
    parseParams.parseParameters.sanitize = sanitize;

    dp_supplier.reset(new v2::FileParsers::MultithreadedSmilesMolSupplier(
        inStream, takeOwnership, params, parseParams));
  }

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
