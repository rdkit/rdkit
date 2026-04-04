//
//  Copyright (C) 2002-2024 greg landrum and other RDKit contributors
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
#include <string_view>
#include <list>
#include <memory>
#include <vector>
#include <fstream>
#include <GraphMol/ROMol.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#endif

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

namespace v2 {
namespace FileParsers {
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
  virtual void init() = 0;
  virtual void reset() = 0;
  virtual bool atEnd() = 0;
  virtual std::unique_ptr<RWMol> next() = 0;

  virtual void close() {
    if (df_owner) {
      delete dp_inStream;
      df_owner = false;
    }
    dp_inStream = nullptr;
  }

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

template <typename Supplier>
struct ForwardSupplierIter {
  using iterator_category = std::input_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = std::shared_ptr<RWMol>;

  Supplier *supplier = nullptr;
  std::optional<value_type> current;
  ForwardSupplierIter() = default;
  explicit ForwardSupplierIter(Supplier *supplier)
      : supplier(supplier), current(supplier->nextShared()) {}
  value_type operator*() const { return current.value(); }
  ForwardSupplierIter &operator++() {
    if (supplier->atEnd()) {
      current.reset();
      return *this;
    }
    current = supplier->nextShared();
    if (!current.value() && supplier->atEnd()) {
      current.reset();
    }
    return *this;
  }
  ForwardSupplierIter operator++(int) {
    if (supplier->atEnd()) {
      current.reset();
      return *this;
    }
    ForwardSupplierIter tmp = *this;
    ++(*this);
    return tmp;
  }
  bool operator==(const ForwardSupplierIter &other) const {
    return !current.has_value() && !other.current.has_value();
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
  ForwardSDMolSupplier() { init(); }

  explicit ForwardSDMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const MolFileParserParams &params = MolFileParserParams());

  ~ForwardSDMolSupplier() override { close(); }

  void init() override;
  void reset() override;
  std::unique_ptr<RWMol> next() override;
  std::shared_ptr<RWMol> nextShared() {
    return std::shared_ptr<RWMol>(this->next());
  };
  bool atEnd() override;

  void setProcessPropertyLists(bool val) { df_processPropertyLists = val; }
  bool getProcessPropertyLists() const { return df_processPropertyLists; }

  bool getEOFHitOnRead() const { return df_eofHitOnRead; }

  ForwardSupplierIter<ForwardSDMolSupplier> begin() {
    if (d_line) {
      throw ValueErrorException(
          "Cannot create an iterator for a ForwardSDMolSupplier that has already "
          "been read from.");
    }
    return ForwardSupplierIter(this);
  }
  ForwardSupplierIter<ForwardSDMolSupplier> end() {
    return ForwardSupplierIter<ForwardSDMolSupplier>();
  }

 protected:
  virtual void checkForEnd();
  std::unique_ptr<RWMol> _next();
  virtual void readMolProps(ROMol &);
  bool df_end = false;
  int d_line = 0;  // line number we are currently on
  MolFileParserParams d_params;
  bool df_processPropertyLists = true;
  bool df_eofHitOnRead = false;
};
// clang-format off
static_assert(
    std::ranges::input_range<ForwardSDMolSupplier>
  );
// clang-format on

template <typename Supplier>
struct RandomAccessSupplierIter {
  using iterator_category = std::random_access_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = std::shared_ptr<RWMol>;

  Supplier *supplier = nullptr;
  size_t current_idx = 0;
  RandomAccessSupplierIter() = default;
  explicit RandomAccessSupplierIter(Supplier *supplier)
      : supplier(supplier), current_idx(0) {}
  RandomAccessSupplierIter(Supplier *supplier, size_t idx)
      : supplier(supplier), current_idx(idx) {}
  value_type operator*() const { return supplier->getShared(current_idx); }
  value_type operator[](size_t idx) const { return supplier->getShared(idx); }
  RandomAccessSupplierIter &operator++() {
    ++current_idx;
    return *this;
  }
  RandomAccessSupplierIter operator++(int) {
    RandomAccessSupplierIter tmp(this->supplier, current_idx);
    ++(*this);
    return tmp;
  }
  RandomAccessSupplierIter &operator--() {
    --current_idx;
    return *this;
  }
  RandomAccessSupplierIter operator--(int) {
    RandomAccessSupplierIter tmp(this->supplier, current_idx);
    --(*this);
    return tmp;
  }
  RandomAccessSupplierIter &operator+=(difference_type n) {
    current_idx += n;
    return *this;
  }
  RandomAccessSupplierIter &operator-=(difference_type n) {
    current_idx -= n;
    return *this;
  }
  RandomAccessSupplierIter operator+(difference_type n) const {
    return RandomAccessSupplierIter(this->supplier, current_idx + n);
  }
  RandomAccessSupplierIter operator-(difference_type n) const {
    return RandomAccessSupplierIter(this->supplier, current_idx - n);
  }
  difference_type operator-(const RandomAccessSupplierIter &other) const {
    return static_cast<difference_type>(current_idx) -
           static_cast<difference_type>(other.current_idx);
  }
  friend RandomAccessSupplierIter operator+(
      difference_type n, const RandomAccessSupplierIter &it) {
    return RandomAccessSupplierIter(it.supplier, it.current_idx + n);
  }
  auto operator<=>(const RandomAccessSupplierIter &other) const {
    return current_idx <=> other.current_idx;
  }
  bool operator==(const RandomAccessSupplierIter &other) const {
    return this->supplier == other.supplier &&
           this->current_idx == other.current_idx;
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
  SDMolSupplier() { init(); }

  /*!
   *   \param fileName - the name of the SD file
   *   \param sanitize - if true sanitize the molecule before returning it
   *   \param removeHs - if true remove Hs from the molecule before returning it
   *                     (triggers sanitization)
   *   \param strictParsing - if set to false, the parser is more lax about
   * correctness
   *                          of the contents.
   */
  explicit SDMolSupplier(
      const std::string &fileName,
      const MolFileParserParams &params = MolFileParserParams());

  explicit SDMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const MolFileParserParams &params = MolFileParserParams());

  ~SDMolSupplier() override { close(); }
  void init() override;
  void reset() override;
  std::unique_ptr<RWMol> next() override;
  bool atEnd() override;
  void moveTo(unsigned int idx);
  std::unique_ptr<RWMol> operator[](unsigned int idx);
  ///! \brief returns a shared pointer to the molecule at the given index.
  std::shared_ptr<RWMol> getShared(unsigned int idx);
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx);
  unsigned int length();
  void setData(const std::string &text);
  void setData(const std::string &text, const MolFileParserParams &params);

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

  RandomAccessSupplierIter<SDMolSupplier> begin() {
    return RandomAccessSupplierIter(this);
  }
  RandomAccessSupplierIter<SDMolSupplier> end() {
    return RandomAccessSupplierIter(this, length());
  }

  void setCaching(bool val) { d_cacheMolecules = val; }
  bool getCaching() const { return d_cacheMolecules; }

 private:
  void checkForEnd() override;
  void peekCheckForEnd(char *bufPtr, char *bufEnd, std::streampos molStartPos);
  void buildIndexTo(unsigned int targetIdx);
  int d_len = 0;   // total number of mol blocks in the file (initialized to -1)
  int d_last = 0;  // the molecule we are ready to read
  std::vector<std::streampos> d_molpos;
  bool d_cacheMolecules = false;
  std::vector<std::optional<std::shared_ptr<RWMol>>> d_molCache;
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::mutex d_readMutex;
  std::mutex d_cacheMutex;
#endif
};
// clang-format off
static_assert(
    std::ranges::random_access_range<SDMolSupplier>
  );
// clang-format on

struct SmilesMolSupplierParams {
  std::string delimiter = " \t";
  int smilesColumn = 0;
  int nameColumn = 1;
  bool titleLine = true;
  v2::SmilesParse::SmilesParserParams parseParameters = {
      true,   // sanitize
      false,  // allowCXSMILES
      true,   // strictCXSMILES
      false,  // parseName
      true,   // removeHs
      false,  // skipCleanup
      false,  // debugParse
      {}      // replacements
  };
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
  explicit SmilesMolSupplier(
      const std::string &fileName,
      const SmilesMolSupplierParams &params = SmilesMolSupplierParams());
  SmilesMolSupplier();
  explicit SmilesMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const SmilesMolSupplierParams &params = SmilesMolSupplierParams());

  ~SmilesMolSupplier() override { close(); }
  void setData(const std::string &text, const SmilesMolSupplierParams &params =
                                            SmilesMolSupplierParams());
  void init() override;
  void reset() override;
  std::unique_ptr<RWMol> next() override;
  bool atEnd() override;
  void moveTo(unsigned int idx);
  std::unique_ptr<RWMol> operator[](unsigned int idx);
  std::shared_ptr<RWMol> getShared(unsigned int idx);
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx);
  unsigned int length();

  RandomAccessSupplierIter<SmilesMolSupplier> begin() {
    return RandomAccessSupplierIter(this);
  }
  RandomAccessSupplierIter<SmilesMolSupplier> end() {
    return RandomAccessSupplierIter(this, length());
  }

  void setCaching(bool val) { d_cacheMolecules = val; }
  bool getCaching() const { return d_cacheMolecules; }

 private:
  std::unique_ptr<RWMol> processLine(std::string inLine);
  void processTitleLine();
  std::string nextLine();
  long int skipComments();
  void checkForEnd();

  bool df_end = false;  // have we reached the end of the file?
  long d_len = 0;       // total number of smiles in the file
  long d_next = 0;      // the  molecule we are ready to read
  size_t d_line = 0;    // line number we are currently on
  SmilesMolSupplierParams d_params;
  std::vector<std::streampos>
      d_molpos;  // vector of positions in the file for molecules
  std::vector<int> d_lineNums;
  STR_VECT d_props;  // vector of property names
  bool d_cacheMolecules = false;
  std::vector<std::optional<std::shared_ptr<RWMol>>> d_molCache;
#ifdef RDK_BUILD_THREADSAFE_SSS
  std::mutex d_readMutex;
  std::mutex d_cacheMutex;
#endif
};

struct TDTMolSupplierParams {
  std::string nameRecord = "";
  int confId2D = -1;
  int confId3D = -1;
  v2::SmilesParse::SmilesParserParams parseParameters = {
      true,   // sanitize
      false,  // allowCXSMILES
      true,   // strictCXSMILES
      false,  // parseName
      true,   // removeHs
      false,  // skipCleanup
      false,  // debugParse
      {}      // replacements
  };
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
  explicit TDTMolSupplier(
      const std::string &fileName,
      const TDTMolSupplierParams &params = TDTMolSupplierParams());
  explicit TDTMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const TDTMolSupplierParams &params = TDTMolSupplierParams());
  TDTMolSupplier();
  ~TDTMolSupplier() override { close(); }
  void setData(const std::string &text,
               const TDTMolSupplierParams &params = TDTMolSupplierParams());
  void init() override;
  void reset() override;
  std::unique_ptr<RWMol> next() override;
  bool atEnd() override;
  void moveTo(unsigned int idx);
  std::unique_ptr<RWMol> operator[](unsigned int idx);
  /*! \brief returns the text block for a particular item
   *
   *  \param idx - which item to return
   */
  std::string getItemText(unsigned int idx);
  unsigned int length();

 private:
  bool advanceToNextRecord();
  void checkForEnd();
  std::unique_ptr<RWMol> parseMol(std::string inLine);

  bool df_end = false;  // have we reached the end of the file?
  int d_len = 0;        // total number of mols in the file
  int d_last = 0;       // the molecule we are ready to read
  int d_line = 0;       // line number we are currently on
  std::vector<std::streampos>
      d_molpos;  // vector of positions in the file for molecules
  TDTMolSupplierParams d_params;
};

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
struct MaeMolSupplierParams {
  bool sanitize = true;
  bool removeHs = true;
};
//! lazy file parser for MAE files
class RDKIT_FILEPARSERS_EXPORT MaeMolSupplier : public MolSupplier {
  /**
   * Due to maeparser's shared_ptr<istream> Reader interface, MaeMolSupplier
   * always requires taking ownership of the istream ptr, as the shared ptr will
   * always clear it upon destruction.
   */

 public:
  MaeMolSupplier() {}

  explicit MaeMolSupplier(
      std::shared_ptr<std::istream> inStream,
      const MaeMolSupplierParams &params = MaeMolSupplierParams());

  explicit MaeMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const MaeMolSupplierParams &params = MaeMolSupplierParams());

  explicit MaeMolSupplier(
      const std::string &fname,
      const MaeMolSupplierParams &params = MaeMolSupplierParams());

  ~MaeMolSupplier() override {}

  void init() override;
  void reset() override;
  std::unique_ptr<RWMol> next() override;
  bool atEnd() override;
  void moveTo(unsigned int idx);
  std::unique_ptr<RWMol> operator[](unsigned int idx);
  unsigned int length();

  void close() override { dp_sInStream.reset(); }

  void setData(const std::string &text,
               const MaeMolSupplierParams &params = MaeMolSupplierParams());

 private:
  void moveToNextBlock();

 protected:
  MaeMolSupplierParams d_params;
  std::shared_ptr<schrodinger::mae::Reader> d_reader;
  std::shared_ptr<schrodinger::mae::Block> d_next_struct;
  std::shared_ptr<std::istream> dp_sInStream;
  std::string d_stored_exc;
  unsigned d_position;
  unsigned d_length;
};
#endif  // RDK_BUILD_MAEPARSER_SUPPORT

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit

#include "MolSupplier.v1API.h"

#endif
