//
// Copyright (c) 2016 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FPBREADER_H_DEC2015
#define RD_FPBREADER_H_DEC2015
/*! \file FPBReader.h

  \brief contains a simple class for reading and searching FPB files

  \b Note that this functionality is experimental and the API may change
     in future releases.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <RDGeneral/BadFileException.h>
#include <DataStructs/ExplicitBitVect.h>

#include <cstdint>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace RDKit {
namespace detail {
struct FPBReader_impl;
}

//! class for reading and searching FPB files
/*!
  basic usage:
  \code
  FPBReader reader("foo.fpb");
  reader.init();
  boost::shared_ptr<ExplicitBitVect> ebv = reader.getFP(95);
  std::vector<std::pair<double, unsigned int> > nbrs =
      reader.getTanimotoNeighbors(*ebv.get(), 0.70);
  \endcode

  \b Note: this functionality is experimental and the API may change
     in future releases.

  <b>Note on thread safety</b>
  Operations that involve reading from the FPB file are not thread safe.
  This means that the \c init() method is not thread safe and none of the
  search operations are thread safe when an \c FPBReader is initialized in
  \c lazyRead mode.

*/
class RDKIT_DATASTRUCTS_EXPORT FPBReader {
 public:
  FPBReader()
      
        {};
  //! ctor for reading from a named file
  /*!
  \param fname the name of the file to reads
  \param lazyRead if set to \c false all fingerprints from the file will be read
  into memory when \c init() is called.
  */
  FPBReader(const char *fname, bool lazyRead = false) {
    _initFromFilename(fname, lazyRead);
  };
  //! \overload
  FPBReader(const std::string &fname, bool lazyRead = false) {
    _initFromFilename(fname.c_str(), lazyRead);
  };
  //! ctor for reading from an open istream
  /*!
  \param inStream the stream to read from
  \param takeOwnership if set, we will take over ownership of the stream pointer
  \param lazyRead if set to \c false all fingerprints from the file will be read
  into memory when \c init() is called.

  Some additional notes:
    - if \c lazyRead is set, \c inStream must support the \c seekg() and \c
  tellg() operations.

  */
  FPBReader(std::istream *inStream, bool takeOwnership = true,
            bool lazyRead = false)
      : dp_istrm(inStream),
        dp_impl(nullptr),
        df_owner(takeOwnership),
        df_init(false),
        df_lazyRead(lazyRead){};
  ~FPBReader() {
    destroy();
    if (df_owner) delete dp_istrm;
    dp_istrm = nullptr;
    df_init = false;
  };

  //! Read the data from the file and initialize internal data structures
  /*!
  This must be called before most of the other methods of this class.

  Some notes:
  \li if \c lazyRead is not set, all fingerprints will be read into memory. This
  can require substantial amounts of memory for large files.
  \li For large files, this can take a long time.
  \li If \c lazyRead and \c takeOwnership are both \c false it is safe to close
  and delete inStream after calling \c init()
  */
  void init();
  //! cleanup
  /*!
  Cleans up whatever memory was allocated during init()
  */
  void cleanup() {
    if (!df_init) return;
    destroy();
    df_init = false;
  };
  //! returns the requested fingerprint as an \c ExplicitBitVect
  boost::shared_ptr<ExplicitBitVect> getFP(unsigned int idx) const;
  //! returns the requested fingerprint as an array of bytes
  boost::shared_array<std::uint8_t> getBytes(unsigned int idx) const;

  //! returns the id of the requested fingerprint
  std::string getId(unsigned int idx) const;
  //! returns the fingerprint and id of the requested fingerprint
  std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> operator[](
      unsigned int idx) const {
    return std::make_pair(getFP(idx), getId(idx));
  };

  //! returns beginning and end indices of fingerprints having on-bit counts
  //! within the range (including end points)
  std::pair<unsigned int, unsigned int> getFPIdsInCountRange(
      unsigned int minCount, unsigned int maxCount);

  //! returns the number of fingerprints
  unsigned int length() const;
  //! returns the number of bits in our fingerprints
  unsigned int nBits() const;

  //! returns the tanimoto similarity between the specified fingerprint and the
  //! provided fingerprint
  double getTanimoto(unsigned int idx, const std::uint8_t *bv) const;
  //! \overload
  double getTanimoto(unsigned int idx,
                     boost::shared_array<std::uint8_t> bv) const {
    return getTanimoto(idx, bv.get());
  };
  //! \overload
  double getTanimoto(unsigned int idx, const ExplicitBitVect &ebv) const;

  //! returns tanimoto neighbors that are within a similarity threshold
  /*!
  The result vector of (similarity,index) pairs is sorted in order
  of decreasing similarity

    \param bv the query fingerprint
    \param threshold the minimum similarity to return
    \param usePopcountScreen if this is true (the default) the popcount of the
           neighbors will be used to reduce the number of calculations that need
           to be done

  */
  std::vector<std::pair<double, unsigned int>> getTanimotoNeighbors(
      const std::uint8_t *bv, double threshold = 0.7,
      bool usePopcountScreen = true) const;
  //! \overload
  std::vector<std::pair<double, unsigned int>> getTanimotoNeighbors(
      boost::shared_array<std::uint8_t> bv, double threshold = 0.7,
      bool usePopcountScreen = true) const {
    return getTanimotoNeighbors(bv.get(), threshold, usePopcountScreen);
  };
  //! \overload
  std::vector<std::pair<double, unsigned int>> getTanimotoNeighbors(
      const ExplicitBitVect &ebv, double threshold = 0.7,
      bool usePopcountScreen = true) const;

  //! returns the Tversky similarity between the specified fingerprint and the
  //! provided fingerprint
  /*!

    \param idx the fingerprint to compare to
    \param bv the query fingerprint
    \param ca the Tversky a coefficient
    \param cb the Tversky a coefficient

   */
  double getTversky(unsigned int idx, const std::uint8_t *bv, double ca,
                    double cb) const;
  //! \overload
  double getTversky(unsigned int idx, boost::shared_array<std::uint8_t> bv,
                    double ca, double cb) const {
    return getTversky(idx, bv.get(), ca, cb);
  };
  //! \overload
  double getTversky(unsigned int idx, const ExplicitBitVect &ebv, double ca,
                    double cb) const;

  //! returns Tversky neighbors that are within a similarity threshold
  /*!
  The result vector of (similarity,index) pairs is sorted in order
  of decreasing similarity

    \param bv the query fingerprint
    \param ca the Tversky a coefficient
    \param cb the Tversky a coefficient
    \param threshold the minimum similarity to return
    \param usePopcountScreen if this is true (the default) the popcount of the
           neighbors will be used to reduce the number of calculations that need
           to be done

  */
  std::vector<std::pair<double, unsigned int>> getTverskyNeighbors(
      const std::uint8_t *bv, double ca, double cb, double threshold = 0.7,
      bool usePopcountScreen = true) const;
  //! \overload
  std::vector<std::pair<double, unsigned int>> getTverskyNeighbors(
      boost::shared_array<std::uint8_t> bv, double ca, double cb,
      double threshold = 0.7, bool usePopcountScreen = true) const {
    return getTverskyNeighbors(bv.get(), ca, cb, threshold, usePopcountScreen);
  };
  //! \overload
  std::vector<std::pair<double, unsigned int>> getTverskyNeighbors(
      const ExplicitBitVect &ebv, double ca, double cb, double threshold = 0.7,
      bool usePopcountScreen = true) const;

  //! returns indices of all fingerprints that completely contain this one
  /*! (i.e. where all the bits set in the query are also set in the db
   molecule)
   */
  std::vector<unsigned int> getContainingNeighbors(
      const std::uint8_t *bv) const;
  //! \overload
  std::vector<unsigned int> getContainingNeighbors(
      boost::shared_array<std::uint8_t> bv) const {
    return getContainingNeighbors(bv.get());
  };
  //! \overload
  std::vector<unsigned int> getContainingNeighbors(
      const ExplicitBitVect &ebv) const;

 private:
  std::istream *dp_istrm{nullptr};
  detail::FPBReader_impl *dp_impl{nullptr};  // implementation details
  bool df_owner{false};
  bool df_init{false};
  bool df_lazyRead{false};

  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  FPBReader(const FPBReader &);
  FPBReader &operator=(const FPBReader &);
  void destroy();
  void _initFromFilename(const char *fname, bool lazyRead) {
    std::istream *tmpStream = static_cast<std::istream *>(
        new std::ifstream(fname, std::ios_base::binary));
    if (!(*tmpStream) || (tmpStream->bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fname;
      delete tmpStream;
      throw BadFileException(errout.str());
    }
    dp_istrm = tmpStream;
    dp_impl = nullptr;
    df_owner = true;
    df_init = false;
    df_lazyRead = lazyRead;
  }
};
}  // namespace RDKit
#endif
