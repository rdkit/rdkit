//
//  Copyright (C) 2024 greg landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_MOLSUPPLIER_v1_H
#define RD_MOLSUPPLIER_v1_H

namespace RDKit {
inline namespace v1 {
/*!
//
//  Here are a couple of ways one can interact with MolSuppliers:
//
//  1) Lazy (ForwardIterator):
//     while(!supplier.atEnd()){
//       ROMol *mol = supplier.next();
//       if(mol){
//           do something;
//       }
//     }
//  2) Random Access:
//     for(int i=0;i<supplier.length();i++){
//       ROMol *mol = supplier[i];
//       if(mol){
//           do something;
//       }
//     }
//
//
*/
class RDKIT_FILEPARSERS_EXPORT MolSupplier {
  // this is an abstract base class to supply molecules one at a time
 public:
  MolSupplier() {}
  virtual ~MolSupplier() {}
  void init() {
    if (dp_supplier) {
      dp_supplier->init();
    }
  }
  void reset() {
    if (dp_supplier) {
      dp_supplier->reset();
    }
  }

  bool atEnd() {
    if (dp_supplier) {
      return dp_supplier->atEnd();
    }
    return true;
  }
  ROMol *next() {
    PRECONDITION(dp_supplier, "no supplier");
    return dp_supplier->next().release();
  }

  virtual void close() {
    if (dp_supplier) {
      dp_supplier->close();
    }
  }

 private:
  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  MolSupplier(const MolSupplier &);
  MolSupplier &operator=(const MolSupplier &);

 protected:
  std::unique_ptr<v2::FileParsers::MolSupplier> dp_supplier;
};

// \brief a supplier from an SD file that only reads forward:
class RDKIT_FILEPARSERS_EXPORT ForwardSDMolSupplier : public MolSupplier {
  /*************************************************************************
   * A lazy mol supplier from a SD file.
   *  - When new molecules are read using "next" their positions in the file are
   *noted.
   ***********************************************************************************/
 public:
  using ContainedType = v2::FileParsers::ForwardSDMolSupplier;
  ForwardSDMolSupplier() {}

  explicit ForwardSDMolSupplier(std::istream *inStream,
                                bool takeOwnership = true, bool sanitize = true,
                                bool removeHs = true,
                                bool strictParsing = false) {
    v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    dp_supplier.reset(new v2::FileParsers::ForwardSDMolSupplier(
        inStream, takeOwnership, params));
  };

  ~ForwardSDMolSupplier() override {}

  void setProcessPropertyLists(bool val) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())
        ->setProcessPropertyLists(val);
  }
  bool getProcessPropertyLists() const {
    if (dp_supplier) {
      return static_cast<ContainedType *>(dp_supplier.get())
          ->getProcessPropertyLists();
    }
    return false;
  }

  bool getEOFHitOnRead() const {
    if (dp_supplier) {
      return static_cast<ContainedType *>(dp_supplier.get())->getEOFHitOnRead();
    }
    return false;
  }
};

// \brief a lazy supplier from an SD file
class RDKIT_FILEPARSERS_EXPORT SDMolSupplier : public ForwardSDMolSupplier {
  /*************************************************************************
   * A lazy mol supplier from a SD file.
   *  - When new molecules are read using "next" their positions in the file are
   *noted.
   *  - A call to the "length" will automatically parse the entire file and
   *cache all the mol
   *    block positions
   *  - [] operator is used to access a molecule at "idx", calling next
   *following this will result
   *    in the next molecule after "idx"
   ***********************************************************************************/

 public:
  using ContainedType = v2::FileParsers::SDMolSupplier;
  SDMolSupplier() { dp_supplier.reset(new ContainedType()); }

  /*!
   *   \param fileName - the name of the SD file
   *   \param sanitize - if true sanitize the molecule before returning it
   *   \param removeHs - if true remove Hs from the molecule before returning it
   *                     (triggers sanitization)
   *   \param strictParsing - if set to false, the parser is more lax about
   * correctness
   *                          of the contents.
   */
  explicit SDMolSupplier(const std::string &fileName, bool sanitize = true,
                         bool removeHs = true, bool strictParsing = true) {
    v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    dp_supplier.reset(new v2::FileParsers::SDMolSupplier(fileName, params));
  }

  explicit SDMolSupplier(std::istream *inStream, bool takeOwnership = true,
                         bool sanitize = true, bool removeHs = true,
                         bool strictParsing = true) {
    v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    dp_supplier.reset(
        new v2::FileParsers::SDMolSupplier(inStream, takeOwnership, params));
  }

  void moveTo(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())->moveTo(idx);
  }
  ROMol *operator[](unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())
        ->
        operator[](idx)
        .release();
  }
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->getItemText(idx);
  }
  unsigned int length() {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->length();
  }
  void setData(const std::string &text, bool sanitize = true,
               bool removeHs = true) {
    PRECONDITION(dp_supplier, "no supplier");
    v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    static_cast<ContainedType *>(dp_supplier.get())->setData(text, params);
  }
  void setData(const std::string &text, bool sanitize, bool removeHs,
               bool strictParsing) {
    v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    static_cast<ContainedType *>(dp_supplier.get())->setData(text, params);
  }
  /*! Resets our internal state and sets the indices of molecules in the stream.
   *  The client should be *very* careful about calling this method, as it's
   *trivial
   *  to end up with a completely useless supplier.
   *
   *   \param locs - the vector of stream positions.
   *
   *  Note that this can be used not only to make reading selected molecules
   *from a
   *  large SD file much faster, but it can also allow subsetting an SD file or
   *  rearranging the order of the molecules.
   */
  void setStreamIndices(const std::vector<std::streampos> &locs) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())->setStreamIndices(locs);
  }
};

//! lazy file parser for Smiles tables
class RDKIT_FILEPARSERS_EXPORT SmilesMolSupplier : public MolSupplier {
  /**************************************************************************
   * Lazy file parser for Smiles table file, similar to the lazy SD
   * file parser above
   * - As an when new molecules are read using "next" their
   *    positions in the file are noted.
   *  - A call to the "length" will automatically parse the entire
   *    file and cache all the mol block positions
   *  - [] operator is used to access a molecule at "idx", calling
   *    next following this will result in the next molecule after
   *    "idx"
   ***************************************************************************/
 public:
  using ContainedType = v2::FileParsers::SmilesMolSupplier;
  /*!
   *   \param fileName - the name of smiles table file
   *   \param delimiter - delimiting characters between records on a each
   *     line NOTE that this is not a string, the tokenizer looks for
   *     the individual characters in delimiter, not the full string
   *     itself.  So the default delimiter: " \t", means " " or "\t".
   *   \param smilesColumn - column number for the SMILES string (defaults
   *     to the first column)
   *   \param nameColumn - column number for the molecule name (defaults to
   *     the second column) If set to -1 we assume that no name is
   *     available for the molecule and the name is defaulted to the
   *     smiles string
   *   \param titleLine - if true, the first line is assumed to list the
   *     names of properties in order separated by 'delimiter'. It is
   *     also assume that the 'SMILES' column and the 'name' column
   *     are not specified here if false - no title line is assumed
   *     and the properties are recorded as the "columnX" where "X" is
   *     the column number
   *   \param sanitize - if true sanitize the molecule before returning it
   */
  explicit SmilesMolSupplier(const std::string &fileName,
                             const std::string &delimiter = " \t",
                             int smilesColumn = 0, int nameColumn = 1,
                             bool titleLine = true, bool sanitize = true) {
    v2::FileParsers::SmilesMolSupplierParams params;
    params.delimiter = delimiter;
    params.smilesColumn = smilesColumn;
    params.nameColumn = nameColumn;
    params.titleLine = titleLine;
    params.parseParameters.sanitize = sanitize;
    dp_supplier.reset(new v2::FileParsers::SmilesMolSupplier(fileName, params));
  }
  explicit SmilesMolSupplier(std::istream *inStream, bool takeOwnership = true,
                             const std::string &delimiter = " \t",
                             int smilesColumn = 0, int nameColumn = 1,
                             bool titleLine = true, bool sanitize = true) {
    v2::FileParsers::SmilesMolSupplierParams params;
    params.delimiter = delimiter;
    params.smilesColumn = smilesColumn;
    params.nameColumn = nameColumn;
    params.titleLine = titleLine;
    params.parseParameters.sanitize = sanitize;
    dp_supplier.reset(new v2::FileParsers::SmilesMolSupplier(
        inStream, takeOwnership, params));
  }
  SmilesMolSupplier() { dp_supplier.reset(new ContainedType()); }

  void setData(const std::string &text, const std::string &delimiter = " ",
               int smilesColumn = 0, int nameColumn = 1, bool titleLine = true,
               bool sanitize = true) {
    PRECONDITION(dp_supplier, "no supplier");
    v2::FileParsers::SmilesMolSupplierParams params;
    params.delimiter = delimiter;
    params.smilesColumn = smilesColumn;
    params.nameColumn = nameColumn;
    params.titleLine = titleLine;
    params.parseParameters.sanitize = sanitize;
    static_cast<ContainedType *>(dp_supplier.get())->setData(text, params);
  }
  void moveTo(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())->moveTo(idx);
  }
  ROMol *operator[](unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())
        ->
        operator[](idx)
        .release();
  }
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->getItemText(idx);
  }
  unsigned int length() {
    PRECONDITION(dp_supplier, "no supplier")
    return static_cast<ContainedType *>(dp_supplier.get())->length();
  }
};

//! lazy file parser for TDT files
class RDKIT_FILEPARSERS_EXPORT TDTMolSupplier : public MolSupplier {
  /**************************************************************************
   * Lazy file parser for TDT files, similar to the lazy SD
   * file parser above
   * - As an when new molecules are read using "next" their
   *    positions in the file are noted.
   *  - A call to the "length" will automatically parse the entire
   *    file and cache all the mol block positions
   *  - [] operator is used to access a molecule at "idx", calling
   *    next following this will result in the next molecule after
   *    "idx"
   ***************************************************************************/
 public:
  using ContainedType = v2::FileParsers::TDTMolSupplier;
  /*!
   *   \param fileName - the name of the TDT file
   *   \param nameRecord - property name for the molecule name.
   *     If empty (the default), the name defaults to be empty
   *   \param confId2D - if >=0 and 2D coordinates are provided, the 2D
   *                   structure (depiction) in the input will be read into the
   *                   corresponding conformer id.
   *   \param confId3D - if >=0 and 3D coordinates are provided, the 3D
   *                   structure (depiction) in the input will be read into the
   *                   corresponding conformer id.
   *   \param sanitize - if true sanitize the molecule before returning it
   */
  explicit TDTMolSupplier(const std::string &fileName,
                          const std::string &nameRecord = "", int confId2D = -1,
                          int confId3D = 0, bool sanitize = true) {
    v2::FileParsers::TDTMolSupplierParams params;
    params.nameRecord = nameRecord;
    params.confId2D = confId2D;
    params.confId3D = confId3D;
    params.parseParameters.sanitize = sanitize;
    dp_supplier.reset(new v2::FileParsers::TDTMolSupplier(fileName, params));
  }
  explicit TDTMolSupplier(std::istream *inStream, bool takeOwnership = true,
                          const std::string &nameRecord = "", int confId2D = -1,
                          int confId3D = 0, bool sanitize = true) {
    v2::FileParsers::TDTMolSupplierParams params;
    params.nameRecord = nameRecord;
    params.confId2D = confId2D;
    params.confId3D = confId3D;
    params.parseParameters.sanitize = sanitize;
    dp_supplier.reset(
        new v2::FileParsers::TDTMolSupplier(inStream, takeOwnership, params));
  }
  TDTMolSupplier() { dp_supplier.reset(new ContainedType()); }
  void setData(const std::string &text, const std::string &nameRecord = "",
               int confId2D = -1, int confId3D = 0, bool sanitize = true) {
    PRECONDITION(dp_supplier, "no supplier");
    v2::FileParsers::TDTMolSupplierParams params;
    params.nameRecord = nameRecord;
    params.confId2D = confId2D;
    params.confId3D = confId3D;
    params.parseParameters.sanitize = sanitize;
    static_cast<ContainedType *>(dp_supplier.get())->setData(text, params);
  }
  void moveTo(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())->moveTo(idx);
  }
  ROMol *operator[](unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())
        ->
        operator[](idx)
        .release();
  }
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->getItemText(idx);
  }
  unsigned int length() {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->length();
  }
};

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
//! lazy file parser for MAE files
class RDKIT_FILEPARSERS_EXPORT MaeMolSupplier : public MolSupplier {
  /**
   * Due to maeparser's shared_ptr<istream> Reader interface, MaeMolSupplier
   * always requires taking ownership of the istream ptr, as the shared ptr will
   * always clear it upon destruction.
   */

 public:
  using ContainedType = v2::FileParsers::MaeMolSupplier;
  MaeMolSupplier() { dp_supplier.reset(new ContainedType()); }

  explicit MaeMolSupplier(std::shared_ptr<std::istream> inStream,
                          bool sanitize = true, bool removeHs = true) {
    v2::FileParsers::MaeMolSupplierParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    dp_supplier.reset(new ContainedType(inStream, params));
  }

  explicit MaeMolSupplier(std::istream *inStream, bool takeOwnership = true,
                          bool sanitize = true, bool removeHs = true) {
    v2::FileParsers::MaeMolSupplierParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    dp_supplier.reset(new ContainedType(inStream, takeOwnership, params));
  }

  explicit MaeMolSupplier(const std::string &fname, bool sanitize = true,
                          bool removeHs = true) {
    v2::FileParsers::MaeMolSupplierParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    dp_supplier.reset(new ContainedType(fname, params));
  }
  void moveTo(unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())->moveTo(idx);
  }
  RWMol *operator[](unsigned int idx) {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())
        ->
        operator[](idx)
        .release();
  }
  unsigned int length() {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())->length();
  }

  void setData(const std::string &text, bool sanitize = true,
               bool removeHs = true) {
    PRECONDITION(dp_supplier, "no supplier");
    v2::FileParsers::MaeMolSupplierParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    static_cast<ContainedType *>(dp_supplier.get())->setData(text, params);
  }
};
#endif  // RDK_BUILD_MAEPARSER_SUPPORT

#if 0

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
  ~MultithreadedSDMolSupplier() override;
  void init() override {}

  void checkForEnd();
  bool getEnd() const override;
  void setProcessPropertyLists(bool val) { df_processPropertyLists = val; }
  bool getProcessPropertyLists() const { return df_processPropertyLists; }
  bool getEOFHitOnRead() const { return df_eofHitOnRead; }

  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum,
                         unsigned int &index) override;
  void readMolProps(RWMol *mol, std::istringstream &inStream);
  //! parses the record and returns the resulting molecule
  RWMol *processMoleculeRecord(const std::string &record,
                               unsigned int lineNum) override;

 private:
  void initFromSettings(bool takeOwnership, bool sanitize, bool removeHs,
                        bool strictParsing, unsigned int numWriterThreads,
                        size_t sizeInputQueue, size_t sizeOutputQueue);

 private:
  bool df_end = false;  //!< have we reached the end of the file?
  int d_line = 0;       //!< line number we are currently on
  bool df_sanitize = true, df_removeHs = true, df_strictParsing = true;
  bool df_processPropertyLists = true;
  bool df_eofHitOnRead = false;
  unsigned int d_currentRecordId = 1;  //!< current record id
};

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
  ROMol *processMoleculeRecord(const std::string &record,
                               unsigned int lineNum) override;

 private:
  void initFromSettings(bool takeOwnership, const std::string &delimiter,
                        int smilesColumn, int nameColumn, bool titleLine,
                        bool sanitize, unsigned int numWriterThreads,
                        size_t sizeInputQueue, size_t sizeOutputQueue);

 private:
  bool df_end = false;      //!< have we reached the end of the file?
  int d_line = 0;           //!< line number we are currently on
  std::string d_delim;      //!< the delimiter string
  bool df_sanitize = true;  //!< sanitize molecules before returning them?
  STR_VECT d_props;         //!< vector of property names
  bool df_title = true;     //!< do we have a title line?
  int d_smi = 0;            //!< column id for the smile string
  int d_name = 1;           //!< column id for the name
  unsigned int d_currentRecordId = 1;  //!< current record id
};

#endif
}  // namespace v1
}  // namespace RDKit

#endif
