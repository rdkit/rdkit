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
class RDKIT_FILEPARSERS_EXPORT MultithreadedSmilesMolSupplier
    : public MultithreadedMolSupplier {
 public:
  explicit MultithreadedSmilesMolSupplier(
      const std::string &fileName, const std::string &delimiter = " \t",
      int smilesColumn = 0, int nameColumn = 1, bool titleLine = true,
      bool sanitize = true, int numWriterThreads = 2, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5);
  explicit MultithreadedSmilesMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const std::string &delimiter = " \t", int smilesColumn = 0,
      int nameColumn = 1, bool titleLine = true, bool sanitize = true,
      int numWriterThreads = 2, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5);
  MultithreadedSmilesMolSupplier();
  ~MultithreadedSmilesMolSupplier();

  //! initialize data members
  void init();
  //! Returns the position of the beginning of the next
  //! non-comment line in the input stream. -1 is returned if
  //! no line could be read;
  long int skipComments();
  //! checks if there is a line to be read from the file
  void checkForEnd();
  //! returns df_end
  bool getEnd();
  //! get next line
  std::string nextLine();
  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum);
  //! parses the record and returns the resulting molecule
  ROMol *processMoleculeRecord(const std::string &record, unsigned int lineNum);

 private:
  bool df_end = false;      // have we reached the end of the file?
  int d_next = 0;           // the  molecule we are ready to read
  int d_line = 0;           // line number we are currently on
  std::string d_delim;      // the delimiter string
  bool df_sanitize = true;  // sanitize molecules before returning them?
  STR_VECT d_props;         // vector of property names
  bool df_title = true;     // do we have a title line?
  int d_smi = 0;            // column id for the smile string
  int d_name = 1;           // column id for the name
};
}  // namespace RDKit
#endif
#endif
