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

//! This class is still a bit experimental and the public API may change
//! in future releases.
class RDKIT_FILEPARSERS_EXPORT MultithreadedSDMolSupplier
    : public MultithreadedMolSupplier {
 public:
  explicit MultithreadedSDMolSupplier(
      const std::string &fileName, bool sanitize = true, bool removeHs = true,
      bool strictParsing = true, unsigned int numWriterThreads = 1,
      size_t sizeInputQueue = 5, size_t sizeOutputQueue = 5);

  explicit MultithreadedSDMolSupplier(
      std::istream *inStream, bool takeOwnership = true, bool sanitize = true,
      bool removeHs = true, bool strictParsing = true,
      unsigned int numWriterThreads = 1, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5);

  MultithreadedSDMolSupplier();
  ~MultithreadedSDMolSupplier();
  void init(){};

  void checkForEnd();
  bool getEnd() const;
  void setProcessPropertyLists(bool val) { df_processPropertyLists = val; }
  bool getProcessPropertyLists() const { return df_processPropertyLists; }
  bool getEOFHitOnRead() const { return df_eofHitOnRead; }

  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum,
                         unsigned int &index);
  void readMolProps(ROMol *mol, std::istringstream &inStream);
  //! parses the record and returns the resulting molecule
  ROMol *processMoleculeRecord(const std::string &record, unsigned int lineNum);

 private:
  void initFromSettings(bool takeOwnership, bool sanitize, bool removeHs,
                        bool strictParsing, unsigned int numWriterThreads,
                        size_t sizeInputQueue, size_t sizeOutputQueue);

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
