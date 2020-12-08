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
#ifndef MULTITHREADED_SMILES_MOL_SUPPLIER
#define MULTITHREADED_SMILES_MOL_SUPPLIER
#include "MultithreadedMolSupplier.h"
namespace RDKit {
//! This class is still a bit experimental and the public API may change
//! in future releases.    
class RDKIT_FILEPARSERS_EXPORT MultithreadedSmilesMolSupplier
    : public MultithreadedMolSupplier {
 public:
  explicit MultithreadedSmilesMolSupplier(
      const std::string &fileName, const std::string &delimiter = " \t",
      int smilesColumn = 0, int nameColumn = 1, bool titleLine = true,
      bool sanitize = true, unsigned int numWriterThreads = 1,
      size_t sizeInputQueue = 5, size_t sizeOutputQueue = 5);

  explicit MultithreadedSmilesMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const std::string &delimiter = " \t", int smilesColumn = 0,
      int nameColumn = 1, bool titleLine = true, bool sanitize = true,
      unsigned int numWriterThreads = 1, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5);
  MultithreadedSmilesMolSupplier();
  ~MultithreadedSmilesMolSupplier();

  void init(){};
  //! returns df_end
  bool getEnd() const;
  //! reads and processes the title line
  void processTitleLine();
  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum,
                         unsigned int &index);
  //! parses the record and returns the resulting molecule
  ROMol *processMoleculeRecord(const std::string &record, unsigned int lineNum);

 private:
  void initFromSettings(bool takeOwnership, const std::string &delimiter,
                        int smilesColumn, int nameColumn, bool titleLine,
                        bool sanitize, unsigned int numWriterThreads,
                        size_t sizeInputQueue, size_t sizeOutputQueue);

 private:
  bool df_end = false;      //! have we reached the end of the file?
  int d_line = 0;           //! line number we are currently on
  std::string d_delim;      //! the delimiter string
  bool df_sanitize = true;  //! sanitize molecules before returning them?
  STR_VECT d_props;         //! vector of property names
  bool df_title = true;     //! do we have a title line?
  int d_smi = 0;            //! column id for the smile string
  int d_name = 1;           //! column id for the name
  unsigned int d_currentRecordId = 1;  //! current record id
};
}  // namespace RDKit
#endif
#endif
