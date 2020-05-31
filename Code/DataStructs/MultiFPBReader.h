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
#ifndef RD_MULTIFPBREADER_H_APR2016
#define RD_MULTIFPBREADER_H_APR2016
/*! \file MultiFPBReader.h

  \brief contains a class for reading and searching collections of FPB files

  \b Note that this functionality is experimental and the API may change
     in future releases.
*/

#include <RDGeneral/Exceptions.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/FPBReader.h>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>

namespace RDKit {

//! class for reading and searching multiple FPB files
/*!
  basic usage:
  \code
  FPBReader r1("foo1.fpb"),r2("foo2.fpb");
  std::vector<FPBReader *> readers;
  readers.append(&r1);
  readers.append(&r2);
  MultiFPBReader fpbs(readers);
  fpbs.init();
  boost::shared_ptr<ExplicitBitVect> ebv = fpbs.getReader(0)->getFP(95);
  std::vector<boost::tuple<double,unsigned int, unsigned int> > nbrs =
      fpbs.getTanimotoNeighbors(*ebv.get(), 0.70);
  \endcode

  \b Note: this functionality is experimental and the API may change
     in future releases.

  <b>Note on thread safety</b>
  Operations that involve reading from FPB files are not thread safe.
  This means that the \c init() method is not thread safe and none of the
  search operations are thread safe when an \c FPBReader is initialized in
  \c lazyRead mode.

*/
class RDKIT_DATASTRUCTS_EXPORT MultiFPBReader {
 public:
  typedef boost::tuple<double, unsigned int, unsigned int> ResultTuple;
  MultiFPBReader()
       {};

  /*!
    \param initOnSearch: if this is true, the \c init() method on child readers
    will not be called until the first search is done. This is useful with large
    FPB readers.
  */
  MultiFPBReader(bool initOnSearch)
      : df_init(false),
        df_initOnSearch(initOnSearch),
        df_takeOwnership(false){};
  /*!
    \param readers: the set of FPBReader objects to use.
    \param takeOwnership: if true, we own the memory for the FPBReaders
    \param initOnSearch: if this is true, the \c init() method on child readers
    will not be called until the first search is done. This is useful with large
    FPB readers.
  */
  MultiFPBReader(std::vector<FPBReader *> &readers, bool takeOwnership = false,
                 bool initOnSearch = false);

  ~MultiFPBReader() {
    df_init = false;
    if (df_takeOwnership) {
      BOOST_FOREACH (FPBReader *rdr, d_readers) { delete rdr; };
      d_readers.clear();
    }
  };

  //! Read the data from the file and initialize internal data structures
  /*!
  This must be called before most of the other methods of this class.
  It calls the \c init() method on each of the child FPBReaders

  */
  void init();

  //! returns the number of readers
  unsigned int length() const { return d_readers.size(); };
  //! returns the number of bits in our fingerprints (all readers are expected
  //! to have the same length)
  unsigned int nBits() const;

  //! returns a particular reader
  /*!

    \param which: the reader to return

  */
  FPBReader *getReader(unsigned int which);

  //! adds a new FPBReader to our list
  /*!

    This does no error checking on the reader, so be careful.

    If \c takeOwnership is \c true then we will take ownership of the memory.

    \param rdr: the reader to add. If we have already been initialized, the
    reader's \c init() method will be called

    \returns a count of the current number of readers
  */
  unsigned int addReader(FPBReader *rdr) {
    PRECONDITION(rdr, "no reader provided");
    d_readers.push_back(rdr);
    if (df_init) rdr->init();
    return d_readers.size();
  };

  //! returns tanimoto neighbors that are within a similarity threshold
  /*!
  The result vector of (similarity,index,reader) tuples is sorted in order
  of decreasing similarity

    \param bv the query fingerprint
    \param threshold the minimum similarity to return
    \param numThreads  Sets the number of threads to use (more than one thread
    will only be used if the RDKit was build with multithread support) If set to
    zero, the max supported by the system will be used.

  */
  std::vector<ResultTuple> getTanimotoNeighbors(const std::uint8_t *bv,
                                                double threshold = 0.7,
                                                int numThreads = 1) const;
  //! \overload
  std::vector<ResultTuple> getTanimotoNeighbors(
      boost::shared_array<std::uint8_t> bv, double threshold = 0.7,
      int numThreads = 1) const {
    return getTanimotoNeighbors(bv.get(), threshold, numThreads);
  };
  //! \overload
  std::vector<ResultTuple> getTanimotoNeighbors(const ExplicitBitVect &ebv,
                                                double threshold = 0.7,
                                                int numThreads = 1) const;

  //! returns Tversky neighbors that are within a similarity threshold
  /*!
  The result vector of (similarity,index) pairs is sorted in order
  of decreasing similarity

    \param bv the query fingerprint
    \param ca the Tversky a coefficient
    \param cb the Tversky a coefficient
    \param threshold the minimum similarity to return
    \param numThreads  Sets the number of threads to use (more than one thread
    will only be used if the RDKit was build with multithread support) If set to
    zero, the max supported by the system will be used.

  */
  std::vector<ResultTuple> getTverskyNeighbors(const std::uint8_t *bv,
                                               double ca, double cb,
                                               double threshold = 0.7,
                                               int numThreads = 1) const;
  //! \overload
  std::vector<ResultTuple> getTverskyNeighbors(
      boost::shared_array<std::uint8_t> bv, double ca, double cb,
      double threshold = 0.7, int numThreads = 1) const {
    return getTverskyNeighbors(bv.get(), ca, cb, threshold, numThreads);
  };
  //! \overload
  std::vector<ResultTuple> getTverskyNeighbors(const ExplicitBitVect &ebv,
                                               double ca, double cb,
                                               double threshold = 0.7,
                                               int numThreads = 1) const;

  //! returns indices of all fingerprints that completely contain this one
  /*! (i.e. where all the bits set in the query are also set in the db
   molecule)
   */
  std::vector<std::pair<unsigned int, unsigned int>> getContainingNeighbors(
      const std::uint8_t *bv, int numThreads = 1) const;
  //! \overload
  std::vector<std::pair<unsigned int, unsigned int>> getContainingNeighbors(
      boost::shared_array<std::uint8_t> bv, int numThreads = 1) const {
    return getContainingNeighbors(bv.get(), numThreads);
  };
  //! \overload
  std::vector<std::pair<unsigned int, unsigned int>> getContainingNeighbors(
      const ExplicitBitVect &ebv, int numThreads = 1) const;

 private:
  std::vector<FPBReader *> d_readers;
  bool df_init{false}, df_initOnSearch{false}, df_takeOwnership{false};

  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  MultiFPBReader(const MultiFPBReader &);
  MultiFPBReader &operator=(const MultiFPBReader &);
};
}  // namespace RDKit
#endif
