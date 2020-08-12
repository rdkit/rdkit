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
#ifndef MULTITHREADED_SD_MOL_SUPPLIER
#define MULTITHREADED_SD_MOL_SUPPLIER
#include "MultithreadedMolSupplier.h"

namespace RDKit {
class RDKIT_FILEPARSERS_EXPORT MultithreadedSDMolSupplier
    : public MultithreadedMolSupplier {
 public:
  explicit MultithreadedSDMolSupplier(
      const std::string &fileName, bool sanitize = true, bool removeHs = true,
      bool strictParsing = true, unsigned int numWriterThreads = 2,
      size_t sizeInputQueue = 5, size_t sizeOutputQueue = 5);

  explicit MultithreadedSDMolSupplier(
      std::istream *inStream, bool takeOwnership = true, bool sanitize = true,
      bool removeHs = true, bool strictParsing = true,
      unsigned int numWriterThreads = 2, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5);

  MultithreadedSDMolSupplier();
  ~MultithreadedSDMolSupplier();
  void init(){};

  //! initialize data members
  void _init(bool takeOwnership, bool sanitize, bool removeHs,
             bool strictParsing, unsigned int numWriterThreads,
             size_t sizeInputQueue, size_t sizeOutputQueue);
  void checkForEnd();
  bool getEnd() const;
  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum,
                         unsigned int &index);
  //! parses the record and returns the resulting molecule
  ROMol *processMoleculeRecord(const std::string &record, unsigned int lineNum);

 private:
  bool df_end = false;  //! have we reached the end of the file?
  int d_line = 0;       //! line number we are currently on
  bool df_sanitize = true, df_removeHs = true, df_strictParsing = true;
  bool df_processPropertyLists = true;
  bool df_eofHitOnRead = false;
  unsigned int d_currentRecordId = 1;  //! current record id
};
}  // namespace RDKit
#endif
#endif
