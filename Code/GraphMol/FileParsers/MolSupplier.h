//
//  Copyright (C) 2002-2019 greg landrum, Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOLSUPPLIER_H
#define RD_MOLSUPPLIER_H

#include <RDGeneral/types.h>

#include <string>
#include <list>
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <GraphMol/ROMol.h>
#include <RDGeneral/BadFileException.h>

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
namespace schrodinger {
namespace mae {
class Reader;
class Block;
}  // namespace mae
}  // namespace schrodinger
#endif  // RDK_BUILD_MAEPARSER_SUPPORT

namespace RDKit {
RDKIT_FILEPARSERS_EXPORT std::string strip(const std::string &orig);

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
  MolSupplier(){};
  virtual ~MolSupplier(){};
  virtual void init() = 0;
  virtual void reset() = 0;
  virtual bool atEnd() = 0;
  virtual ROMol *next() = 0;

 private:
  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  MolSupplier(const MolSupplier &);
  MolSupplier &operator=(const MolSupplier &);

 protected:
  // stream to read the molecules from:
  std::istream *dp_inStream = nullptr;
  // do we own dp_inStream?
  bool df_owner = false;
  // opens a stream for reading and verifies that it can be read from.
  // if not it throws an exception
  // the caller owns the resulting stream
  std::istream *openAndCheckStream(const std::string &filename) {
    // FIX: this binary mode of opening file is here because of a bug in
    // VC++ 6.0
    // the function "tellg" does not work correctly if we do not open it this
    // way
    //   Jan 2009: Confirmed that this is still the case in visual studio 2008
    std::ifstream *strm =
        new std::ifstream(filename.c_str(), std::ios_base::binary);
    if ((!(*strm)) || strm->bad()) {
      std::ostringstream errout;
      errout << "Bad input file " << filename;
      delete strm;
      throw BadFileException(errout.str());
    }

    strm->peek();
    if (strm->bad() || strm->eof()) {
      std::ostringstream errout;
      errout << "Invalid input file " << filename;
      delete strm;
      throw BadFileException(errout.str());
    }
    return static_cast<std::istream *>(strm);
  }
};

// \brief a supplier from an SD file that only reads forward:
class RDKIT_FILEPARSERS_EXPORT ForwardSDMolSupplier : public MolSupplier {
  /*************************************************************************
   * A lazy mol supplier from a SD file.
   *  - When new molecules are read using "next" their positions in the file are
   *noted.
   ***********************************************************************************/
 public:
  ForwardSDMolSupplier() { init(); };

  explicit ForwardSDMolSupplier(std::istream *inStream,
                                bool takeOwnership = true, bool sanitize = true,
                                bool removeHs = true,
                                bool strictParsing = false);

  virtual ~ForwardSDMolSupplier() {
    if (df_owner && dp_inStream) {
      delete dp_inStream;
      df_owner = false;
      dp_inStream = nullptr;
    }
  };

  virtual void init();
  virtual void reset();
  virtual ROMol *next();
  virtual bool atEnd();

  void setProcessPropertyLists(bool val) { df_processPropertyLists = val; }
  bool getProcessPropertyLists() const { return df_processPropertyLists; }

  bool getEOFHitOnRead() const { return df_eofHitOnRead; }

 protected:
  virtual void checkForEnd();
  ROMol *_next();
  virtual void readMolProps(ROMol *);
  bool df_end = false;
  int d_line = 0;  // line number we are currently on
  bool df_sanitize = true, df_removeHs = true, df_strictParsing = true;
  bool df_processPropertyLists = true;
  bool df_eofHitOnRead = false;
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
  SDMolSupplier() { init(); };

  /*!
   *   \param fileName - the name of the SD file
   *   \param sanitize - if true sanitize the molecule before returning it
   *   \param removeHs - if true remove Hs from the molecule before returning it
   *                     (triggers sanitization)
   *   \param strictParsing - if not set, the parser is more lax about
   * correctness
   *                          of the contents.
   */
  explicit SDMolSupplier(const std::string &fileName, bool sanitize = true,
                         bool removeHs = true, bool strictParsing = true);

  explicit SDMolSupplier(std::istream *inStream, bool takeOwnership = true,
                         bool sanitize = true, bool removeHs = true,
                         bool strictParsing = true);

  ~SDMolSupplier(){};
  void init();
  void reset();
  ROMol *next();
  bool atEnd();
  void moveTo(unsigned int idx);
  ROMol *operator[](unsigned int idx);
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx);
  unsigned int length();
  void setData(const std::string &text, bool sanitize = true,
               bool removeHs = true);
  void setData(const std::string &text, bool sanitize, bool removeHs,
               bool strictParsing);

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
  void setStreamIndices(const std::vector<std::streampos> &locs);

 private:
  void checkForEnd();
  void setDataCommon(const std::string &text, bool sanitize, bool removeHs);
  int d_len = 0;   // total number of mol blocks in the file (initialized to -1)
  int d_last = 0;  // the molecule we are ready to read
  std::vector<std::streampos> d_molpos;
};

//! lazy file parser for Smiles tables
class RDKIT_FILEPARSERS_EXPORT SmilesMolSupplier : public MolSupplier {
  /**************************************************************************
   * Lazy file parser for Smiles table file, similar to the lazy SD
   * file parser above
   * - As an when new molecules are read using "next" their
   *    positions in the file are noted.
   *  - A call to the "length" will autamatically parse the entire
   *    file and cache all the mol block positions
   *  - [] operator is used to access a molecule at "idx", calling
   *    next following this will result in the next molecule after
   *    "idx"
   ***************************************************************************/
 public:
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
                             bool titleLine = true, bool sanitize = true);
  SmilesMolSupplier();
  explicit SmilesMolSupplier(std::istream *inStream, bool takeOwnership = true,
                             const std::string &delimiter = " \t",
                             int smilesColumn = 0, int nameColumn = 1,
                             bool titleLine = true, bool sanitize = true);

  ~SmilesMolSupplier();
  void setData(const std::string &text, const std::string &delimiter = " ",
               int smilesColumn = 0, int nameColumn = 1, bool titleLine = true,
               bool sanitize = true);
  void init();
  void reset();
  ROMol *next();
  bool atEnd();
  void moveTo(unsigned int idx);
  ROMol *operator[](unsigned int idx);
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx);
  unsigned int length();

 private:
  ROMol *processLine(std::string inLine);
  void processTitleLine();
  std::string nextLine();
  long int skipComments();
  void checkForEnd();

  bool df_end = false;  // have we reached the end of the file?
  int d_len = 0;        // total number of smiles in the file
  int d_next = 0;       // the  molecule we are ready to read
  int d_line = 0;       // line number we are currently on
  std::vector<std::streampos>
      d_molpos;  // vector of positions in the file for molecules
  std::vector<int> d_lineNums;
  std::string d_delim;      // the delimiter string
  bool df_sanitize = true;  // sanitize molecules before returning them?
  STR_VECT d_props;         // vector of property names
  bool df_title = true;     // do we have a title line?
  int d_smi = 0;            // column id for the smile string
  int d_name = 1;           // column id for the name
};

//! lazy file parser for TDT files
class RDKIT_FILEPARSERS_EXPORT TDTMolSupplier : public MolSupplier {
  /**************************************************************************
   * Lazy file parser for TDT files, similar to the lazy SD
   * file parser above
   * - As an when new molecules are read using "next" their
   *    positions in the file are noted.
   *  - A call to the "length" will autamatically parse the entire
   *    file and cache all the mol block positions
   *  - [] operator is used to access a molecule at "idx", calling
   *    next following this will result in the next molecule after
   *    "idx"
   ***************************************************************************/
 public:
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
                          int confId3D = 0, bool sanitize = true);
  explicit TDTMolSupplier(std::istream *inStream, bool takeOwnership = true,
                          const std::string &nameRecord = "", int confId2D = -1,
                          int confId3D = 0, bool sanitize = true);
  TDTMolSupplier();
  ~TDTMolSupplier();
  void setData(const std::string &text, const std::string &nameRecord = "",
               int confId2D = -1, int confId3D = 0, bool sanitize = true);
  void init();
  void reset();
  ROMol *next();
  bool atEnd();
  void moveTo(unsigned int idx);
  ROMol *operator[](unsigned int idx);
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx);
  unsigned int length();

 private:
  bool advanceToNextRecord();
  void checkForEnd();
  ROMol *parseMol(std::string inLine);

  bool df_end = false;  // have we reached the end of the file?
  int d_len = 0;        // total number of mols in the file
  int d_last = 0;       // the molecule we are ready to read
  int d_line = 0;       // line number we are currently on
  int d_confId2D = -1;  // id to use for 2D conformers
  int d_confId3D = 0;   // id to use for 3D conformers
  std::vector<std::streampos>
      d_molpos;             // vector of positions in the file for molecules
  bool df_sanitize = true;  // sanitize molecules before returning them?
  std::string d_nameProp =
      "";  // local storage for the property providing mol names
};

//! lazy file parser for PDB files
class RDKIT_FILEPARSERS_EXPORT PDBMolSupplier : public MolSupplier {
 public:
  explicit PDBMolSupplier(std::istream *inStream, bool takeOwnership = true,
                          bool sanitize = true, bool removeHs = true,
                          unsigned int flavor = 0,
                          bool proximityBonding = true);
  explicit PDBMolSupplier(const std::string &fname, bool sanitize = true,
                          bool removeHs = true, unsigned int flavor = 0,
                          bool proximityBonding = true);

  virtual ~PDBMolSupplier() {
    if (df_owner && dp_inStream) delete dp_inStream;
  };

  virtual void init();
  virtual void reset();
  virtual ROMol *next();
  virtual bool atEnd();

 protected:
  bool df_sanitize, df_removeHs, df_proximityBonding;
  unsigned int d_flavor;
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
  MaeMolSupplier() { init(); };

  explicit MaeMolSupplier(std::shared_ptr<std::istream> inStream,
                          bool sanitize = true, bool removeHs = true);

  explicit MaeMolSupplier(std::istream *inStream, bool takeOwnership = true,
                          bool sanitize = true, bool removeHs = true);

  explicit MaeMolSupplier(const std::string &fname, bool sanitize = true,
                          bool removeHs = true);

  virtual ~MaeMolSupplier(){
      // The dp_sInStream shared_ptr will take care of cleaning up.
  };

  virtual void init();
  virtual void reset();
  virtual ROMol *next();
  virtual bool atEnd();

 protected:
  bool df_sanitize, df_removeHs;
  std::shared_ptr<schrodinger::mae::Reader> d_reader;
  std::shared_ptr<schrodinger::mae::Block> d_next_struct;
  std::shared_ptr<std::istream> dp_sInStream;
  std::string d_stored_exc;
};
#endif  // RDK_BUILD_MAEPARSER_SUPPORT
}  // namespace RDKit

#endif
