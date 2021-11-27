//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef RDK_SUBSTRUCT_LIBRARY
#define RDK_SUBSTRUCT_LIBRARY
#include <utility>

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>

#include <algorithm>
#include <string>
#include <boost/lexical_cast.hpp>

namespace RDKit {

RDKIT_SUBSTRUCTLIBRARY_EXPORT bool SubstructLibraryCanSerialize();

//! Base class API for holding molecules to substructure search.
/*!
  This is an API that hides the implementation details used for
  indexing molecules for substructure searching.  It simply
  provides an API for adding and getting molecules from a set.
 */
class RDKIT_SUBSTRUCTLIBRARY_EXPORT MolHolderBase {
 public:
  virtual ~MolHolderBase() {}

  //! Add a new molecule to the substructure search library
  //!  Returns the molecules index in the library
  virtual unsigned int addMol(const ROMol &m) = 0;

  // implementations should throw IndexError on out of range
  virtual boost::shared_ptr<ROMol> getMol(unsigned int) const = 0;

  //! Get the current library size
  virtual unsigned int size() const = 0;
};

//! Concrete class that holds molecules in memory
/*!
    This is currently one of the faster implementations.
    However it is very memory intensive.
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT MolHolder : public MolHolderBase {
  std::vector<boost::shared_ptr<ROMol>> mols;

 public:
  MolHolder() : MolHolderBase(), mols() {}

  unsigned int addMol(const ROMol &m) override {
    mols.push_back(boost::make_shared<ROMol>(m));
    return size() - 1;
  }

  boost::shared_ptr<ROMol> getMol(unsigned int idx) const override {
    if (idx >= mols.size()) {
      throw IndexErrorException(idx);
    }
    return mols[idx];
  }

  unsigned int size() const override {
    return rdcast<unsigned int>(mols.size());
  }

  std::vector<boost::shared_ptr<ROMol>> &getMols() { return mols; }
  const std::vector<boost::shared_ptr<ROMol>> &getMols() const { return mols; }
};

//! Concrete class that holds binary cached molecules in memory
/*!
  This implementation uses quite a bit less memory than the
  non cached implementation.  However, due to the reduced speed
  it should be used in conjunction with a pattern fingerprinter.

  See RDKit::FPHolder
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT CachedMolHolder : public MolHolderBase {
  std::vector<std::string> mols;

 public:
  CachedMolHolder() : MolHolderBase(), mols() {}

  unsigned int addMol(const ROMol &m) override {
    mols.emplace_back();
    MolPickler::pickleMol(m, mols.back());
    return size() - 1;
  }

  //! Adds a pickled binary molecule, no validity checking of the input
  //!  is done.
  unsigned int addBinary(const std::string &pickle) {
    mols.push_back(pickle);
    return size() - 1;
  }

  boost::shared_ptr<ROMol> getMol(unsigned int idx) const override {
    if (idx >= mols.size()) {
      throw IndexErrorException(idx);
    }
    boost::shared_ptr<ROMol> mol(new ROMol);
    MolPickler::molFromPickle(mols[idx], mol.get());
    return mol;
  }

  unsigned int size() const override {
    return rdcast<unsigned int>(mols.size());
  }

  std::vector<std::string> &getMols() { return mols; }
  const std::vector<std::string> &getMols() const { return mols; }
};

//! Concrete class that holds smiles strings in memory
/*!
    This implementation uses quite a bit less memory than the
    cached binary or uncached implementation.  However, due to the
    reduced speed it should be used in conjunction with a pattern
    fingerprinter.

    See RDKit::FPHolder
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT CachedSmilesMolHolder
    : public MolHolderBase {
  std::vector<std::string> mols;

 public:
  CachedSmilesMolHolder() : MolHolderBase(), mols() {}

  unsigned int addMol(const ROMol &m) override {
    bool doIsomericSmiles = true;
    mols.push_back(MolToSmiles(m, doIsomericSmiles));
    return size() - 1;
  }

  //! Add a smiles to the dataset, no validation is done
  //! to the inputs.
  unsigned int addSmiles(const std::string &smiles) {
    mols.push_back(smiles);
    return size() - 1;
  }

  boost::shared_ptr<ROMol> getMol(unsigned int idx) const override {
    if (idx >= mols.size()) {
      throw IndexErrorException(idx);
    }

    boost::shared_ptr<ROMol> mol(SmilesToMol(mols[idx]));
    return mol;
  }

  unsigned int size() const override {
    return rdcast<unsigned int>(mols.size());
  }

  std::vector<std::string> &getMols() { return mols; }
  const std::vector<std::string> &getMols() const { return mols; }
};

//! Concrete class that holds trusted smiles strings in memory
/*!
    A trusted smiles is essentially a smiles string that
    RDKit has generated.  This indicates that fewer
    sanitization steps are required.  See
    http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html

    This implementation uses quite a bit less memory than the
    cached binary or uncached implementation.  However, due to the
    reduced speed it should be used in conjunction with a pattern
    fingerprinter.

    See RDKit::FPHolder
*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT CachedTrustedSmilesMolHolder
    : public MolHolderBase {
  std::vector<std::string> mols;

 public:
  CachedTrustedSmilesMolHolder() : MolHolderBase(), mols() {}

  unsigned int addMol(const ROMol &m) override {
    bool doIsomericSmiles = true;
    mols.push_back(MolToSmiles(m, doIsomericSmiles));
    return size() - 1;
  }

  //! Add a smiles to the dataset, no validation is done
  //! to the inputs.
  unsigned int addSmiles(const std::string &smiles) {
    mols.push_back(smiles);
    return size() - 1;
  }

  boost::shared_ptr<ROMol> getMol(unsigned int idx) const override {
    if (idx >= mols.size()) {
      throw IndexErrorException(idx);
    }

    RWMol *m = SmilesToMol(mols[idx], 0, false);
    if (m) {
      m->updatePropertyCache();
    }
    return boost::shared_ptr<ROMol>(m);
  }

  unsigned int size() const override {
    return rdcast<unsigned int>(mols.size());
  }

  std::vector<std::string> &getMols() { return mols; }
  const std::vector<std::string> &getMols() const { return mols; }
};

//! Base FPI for the fingerprinter used to rule out impossible matches
class RDKIT_SUBSTRUCTLIBRARY_EXPORT FPHolderBase {
  std::vector<ExplicitBitVect *> fps;

 public:
  virtual ~FPHolderBase() {
    for (size_t i = 0; i < fps.size(); ++i) {
      delete fps[i];
    }
  }

  virtual unsigned int size() const { return rdcast<unsigned int>(fps.size()); }

  //! Adds a molecule to the fingerprinter
  unsigned int addMol(const ROMol &m) {
    fps.push_back(makeFingerprint(m));
    return rdcast<unsigned int>(fps.size() - 1);
  }

  //! Adds a raw bit vector pointer to the fingerprinter, which takes ownership
  //! PLEASE NOTE: make sure that the passed ExplicitBitVect
  //! is compatible with the one generated by makeFingerprint()
  unsigned int addFingerprint(ExplicitBitVect *v) {
    fps.push_back(v);
    return rdcast<unsigned int>(fps.size() - 1);
  }

  //! Adds a raw bit vector to the fingerprinter
  //! PLEASE NOTE: make sure that the passed ExplicitBitVect
  //! is compatible with the one generated by makeFingerprint()
  unsigned int addFingerprint(const ExplicitBitVect &v) {
    return addFingerprint(new ExplicitBitVect(v));
  }

  //! Return false if a substructure search can never match the molecule
  bool passesFilter(unsigned int idx, const ExplicitBitVect &query) const {
    if (idx >= fps.size()) {
      throw IndexErrorException(idx);
    }

    return AllProbeBitsMatch(query, *fps[idx]);
  }

  //! Get the bit vector at the specified index (throws IndexError if out of
  //! range)
  const ExplicitBitVect &getFingerprint(unsigned int idx) const {
    if (idx >= fps.size()) {
      throw IndexErrorException(idx);
    }
    return *fps[idx];
  }

  //! make the query vector
  //!  Caller owns the vector!
  virtual ExplicitBitVect *makeFingerprint(const ROMol &m) const = 0;

  std::vector<ExplicitBitVect *> &getFingerprints() { return fps; }
  const std::vector<ExplicitBitVect *> &getFingerprints() const { return fps; }
};

//! Uses the pattern fingerprinter with a user-defined number of bits (default:
//! 2048) to rule out matches
class RDKIT_SUBSTRUCTLIBRARY_EXPORT PatternHolder : public FPHolderBase {
  unsigned int numBits;

 public:
  PatternHolder() : FPHolderBase(), numBits(defaultNumBits()) {}
  PatternHolder(unsigned int numBits) : FPHolderBase(), numBits(numBits) {}
  //! Caller owns the vector!
  ExplicitBitVect *makeFingerprint(const ROMol &m) const override {
    return PatternFingerprintMol(m, numBits);
  }
  const unsigned int &getNumBits() const { return numBits; };
  unsigned int &getNumBits() { return numBits; };
  static unsigned int defaultNumBits() {
    static const unsigned int DEFAULT_NUM_BITS = 2048;
    return DEFAULT_NUM_BITS;
  };
};

class RDKIT_SUBSTRUCTLIBRARY_EXPORT TautomerPatternHolder
    : public PatternHolder {
 public:
  TautomerPatternHolder() : PatternHolder() {}
  TautomerPatternHolder(unsigned int numBits) : PatternHolder(numBits) {}
  ExplicitBitVect *makeFingerprint(const ROMol &m) const override {
    std::vector<unsigned int> *atomCounts = nullptr;
    ExplicitBitVect *setOnlyBits = nullptr;
    const bool tautomericFingerprint = true;
    return PatternFingerprintMol(m, getNumBits(), atomCounts, setOnlyBits,
                                 tautomericFingerprint);
  }
};

class RDKIT_SUBSTRUCTLIBRARY_EXPORT KeyHolderBase {
 public:
  virtual ~KeyHolderBase() {}

  //! Add a key to the database getting it from the molecule
  virtual unsigned int addMol(const ROMol &m) = 0;

  //! Add a key to the database, this needs to be in the same order
  //!  as the molecule, no validation is done
  virtual unsigned int addKey(const std::string &) = 0;

  // !get the key at the requested index
  // implementations should throw IndexError on out of range
  virtual const std::string &getKey(unsigned int) const = 0;

  // !get keys from a bunch of indices
  virtual std::vector<std::string> getKeys(
      const std::vector<unsigned int> &indices) const = 0;
  //! Get the current keeyholder size
  virtual unsigned int size() const = 0;
};

class RDKIT_SUBSTRUCTLIBRARY_EXPORT KeyFromPropHolder : public KeyHolderBase {
  std::string propname;
  std::vector<std::string> keys;
  const std::string empty_string = {};

 public:
  KeyFromPropHolder(const std::string &propname = "_Name")
      : propname(propname) {}

  std::string &getPropName() { return propname; }
  const std::string &getPropName() const { return propname; }

  std::vector<std::string> &getKeys() { return keys; }
  const std::vector<std::string> &getKeys() const { return keys; }

  unsigned int addMol(const ROMol &m) override {
    std::string key;
    if (m.getPropIfPresent(propname, key)) {
      keys.push_back(std::move(key));
    } else {
      // XXX is this a warning? it could be verbose.  Should we push back the
      // string repr of the
      //  numeric index?
      const static std::string prefix("LIBIDX-");
      keys.emplace_back(prefix + boost::lexical_cast<std::string>(keys.size()));
    }
    return keys.size() - 1u;
  };

  unsigned int addKey(const std::string &key) override {
    keys.push_back(key);
    return keys.size() - 1u;
  }

  const std::string &getKey(unsigned int idx) const override {
    if (idx >= keys.size()) {
      throw IndexErrorException(idx);
    }
    return keys[idx];
  }

  std::vector<std::string> getKeys(
      const std::vector<unsigned int> &indices) const override {
    std::vector<std::string> res;
    std::transform(indices.begin(), indices.end(), std::back_inserter(res),
                   [=](unsigned idx) { return keys.at(idx); });
    return res;
  }
  unsigned int size() const override { return keys.size(); }
};

//! Substructure Search a library of molecules
/*!  This class allows for multithreaded substructure searches of
     large datasets.

     The implementations can use fingerprints to speed up searches
     and have molecules cached as binary forms to reduce memory
     usage.

     basic usage:
     \code
     SubstructLibrary lib;
     lib.addMol(mol);
     std::vector<unsigned int> results = lib.getMatches(query);
     for(std::vector<unsigned int>::const_iterator matchIndex=results.begin();
             matchIndex != results.end();
             ++matchIndex) {
       boost::shared_ptr<ROMol> match = lib.getMol(*matchIndex);
     }
     \endcode

     Using different mol holders and pattern fingerprints.

     \code
     boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder = \
        boost::make_shared<CachedTrustedSmilesMolHolder>();
     boost::shared_ptr<PatternHolder> patternHolder = \
        boost::make_shared<PatternHolder>();

     SubstructLibrary lib(molHolder, patternHolder);
     lib.addMol(mol);
     \endcode

     Cached molecule holders create molecules on demand.  There are currently
     three styles of cached molecules.

       CachedMolHolder: stores molecules in the rdkit binary format.
       CachedSmilesMolHolder: stores molecules in smiles format.
       CachedTrustedSmilesMolHolder: stores molecules in smiles format.

     The CachedTrustedSmilesMolHolder is made to add molecules from
     a trusted source.  This makes the basic assumption that RDKit was
     used to sanitize and canonicalize the smiles string.  In practice
     this is considerably faster than using arbitrary smiles strings since
     certain assumptions can be made.  Molecules generated from trusted
     smiles do not have ring information (although this is created
     in the molecule being searched if necessary).

     When loading from external data, as opposed to using the "addMol" API,
     care must be taken to ensure that the pattern fingerprints and smiles
     are synchronized.

     Each pattern holder has an API point for making its fingerprint.  This
     is useful to ensure that the pattern stored in the database will be
     compatible with the patterns made when analyzing queries.

     \code
     boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder = \
         boost::make_shared<CachedTrustedSmilesMolHolder>();
     boost::shared_ptr<PatternHolder> patternHolder =    \
         boost::make_shared<PatternHolder>();

     // the PatternHolder instance is able to make fingerprints.
     //  These, of course, can be read from a file.  For demonstration
     //   purposes we construct them here.
     const std::string trustedSmiles = "c1ccccc1";
     ROMol *m = SmilesToMol(trustedSmiles);
     const ExplicitBitVect *bitVector = patternHolder->makeFingerprint(*m);

     // The trusted smiles and bitVector can be read from any source.
     //  This is the fastest way to load a substruct library.
     molHolder->addSmiles( trustedSmiles );
     patternHolder->addFingerprint( *bitVector );
     SubstructLibrary lib(molHolder, patternHolder);
     delete m;
     delete bitVector;
     \endcode

     Finally, using the KeyFromPropHolder will store user ids or keys.
     By default, it uses RDKit's default _Name prop, but can be changed
     to any property.

     \code
     boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder = \
         boost::make_shared<CachedTrustedSmilesMolHolder>();
     boost::shared_ptr<KeyFromPropHolder> keyHolder = \
         boost::make_shared<KeyFromPropHolder>();
     SubstructLibrary lib(molHolder, keyHolder);
     ...

     You can get the keys in multiple through the use of the keyholder
     auto key = lib.getKeys().getKey(idx);
     auto keys = lib.getKeys().getKeys(lib.GetMatch(query));
     \endcode

*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT SubstructLibrary {
  boost::shared_ptr<MolHolderBase> molholder;
  boost::shared_ptr<FPHolderBase> fpholder;
  boost::shared_ptr<KeyHolderBase> keyholder;

  MolHolderBase *mols;  // used for a small optimization
  FPHolderBase *fps{nullptr};
  bool is_tautomerquery = false;
  std::vector<unsigned int> searchOrder;

 public:
  SubstructLibrary()
      : molholder(new MolHolder),
        fpholder(),
        keyholder(),
        mols(molholder.get()) {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules)
      : molholder(std::move(molecules)),
        fpholder(),
        keyholder(),
        mols(molholder.get()),
        fps(nullptr) {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<FPHolderBase> fingerprints)
      : molholder(std::move(molecules)),
        fpholder(std::move(fingerprints)),
        keyholder(),
        mols(molholder.get()),
        fps(fpholder.get()) {
    if (fpholder.get() &&
        dynamic_cast<TautomerPatternHolder *>(fpholder.get()) != nullptr) {
      is_tautomerquery = true;
    }
  }

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<KeyHolderBase> keys)
      : molholder(std::move(molecules)),
        fpholder(),
        keyholder(std::move(keys)),
        mols(molholder.get()),
        fps(nullptr) {
    if (fpholder.get() &&
        dynamic_cast<TautomerPatternHolder *>(fpholder.get()) != nullptr) {
      is_tautomerquery = true;
    }
  }

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<FPHolderBase> fingerprints,
                   boost::shared_ptr<KeyHolderBase> keys)
      : molholder(std::move(molecules)),
        fpholder(std::move(fingerprints)),
        keyholder(std::move(keys)),
        mols(molholder.get()),
        fps(fpholder.get()) {
    if (fpholder.get() &&
        dynamic_cast<TautomerPatternHolder *>(fpholder.get()) != nullptr) {
      is_tautomerquery = true;
    }
  }

  SubstructLibrary(const std::string &pickle)
      : molholder(new MolHolder),
        fpholder(),
        mols(molholder.get()),
        fps(nullptr) {
    initFromString(pickle);
    if (fpholder.get() &&
        dynamic_cast<TautomerPatternHolder *>(fpholder.get()) != nullptr) {
      is_tautomerquery = true;
    }
  }

  //! Get the underlying molecule holder implementation
  boost::shared_ptr<MolHolderBase> &getMolHolder() { return molholder; }

  const boost::shared_ptr<MolHolderBase> &getMolHolder() const {
    return molholder;
  }

  //! Get the underlying molecule holder implementation
  boost::shared_ptr<FPHolderBase> &getFpHolder() { return fpholder; }

  //! Get the underlying molecule holder implementation
  const boost::shared_ptr<FPHolderBase> &getFpHolder() const {
    return fpholder;
  }

  //! Get the underlying molecule holder implementation
  boost::shared_ptr<KeyHolderBase> &getKeyHolder() { return keyholder; }

  //! Get the underlying molecule holder implementation
  const boost::shared_ptr<KeyHolderBase> &getKeyHolder() const {
    return keyholder;
  }

  const MolHolderBase &getMolecules() const {
    PRECONDITION(mols, "Molecule holder NULL in SubstructLibrary");
    return *mols;
  }

  //! Get the underlying fingerprint implementation.
  /*! Throws a value error if no fingerprints have been set */
  FPHolderBase &getFingerprints() {
    if (!fps) {
      throw ValueErrorException("Substruct Library does not have fingerprints");
    }
    return *fps;
  }

  const FPHolderBase &getFingerprints() const {
    if (!fps) {
      throw ValueErrorException("Substruct Library does not have fingerprints");
    }
    return *fps;
  }

  //! Get the underlying key holder implementation.
  /*! Throws a value error if no keyholder have been set */
  KeyHolderBase &getKeys() {
    if (!keyholder.get()) {
      throw ValueErrorException("Substruct Library does not have fingerprints");
    }
    return *keyholder.get();
  }

  //! Get the underlying key holder implementation.
  /*! Throws a value error if no keyholder have been set */
  const KeyHolderBase &getKeys() const {
    if (!keyholder.get()) {
      throw ValueErrorException("Substruct Library does not have fingerprints");
    }
    return *keyholder.get();
  }

  //! Add a molecule to the library
  /*!
    \param mol Molecule to add

    returns index for the molecule in the library
  */
  unsigned int addMol(const ROMol &mol);

  //! Get the matching indices for the query
  /*!
    \param query       Query or Tautomer Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
                              [default true]
    \param useChirality  use atomic CIP codes as part of the comparison
                         [default true]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
                                 [default false]
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all
                       [default -1]
  */
  template <class Query>
  std::vector<unsigned int> getMatches(const Query &query,
                                       bool recursionPossible = true,
                                       bool useChirality = true,
                                       bool useQueryQueryMatches = false,
                                       int numThreads = -1,
                                       int maxResults = -1) const {
    SubstructMatchParameters params;
    params.recursionPossible = recursionPossible;
    params.useChirality = useChirality;
    params.useQueryQueryMatches = useQueryQueryMatches;
    return getMatches(query, 0, size(), params, numThreads, maxResults);
  }
  //! overload
  template <class Query>
  std::vector<unsigned int> getMatches(const Query &query,
                                       const SubstructMatchParameters &params,
                                       int numThreads = -1,
                                       int maxResults = -1) const {
    return getMatches(query, 0, size(), params, numThreads, maxResults);
  }
  //! Get the matching indices for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (non-inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
                       [default true]
    \param useChirality  use atomic CIP codes as part of the comparison
                       [default true]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
                                [default false]
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all
                       [default -1]
  */
  template <class Query>
  std::vector<unsigned int> getMatches(
      const Query &query, unsigned int startIdx, unsigned int endIdx,
      bool recursionPossible = true, bool useChirality = true,
      bool useQueryQueryMatches = false, int numThreads = -1,
      int maxResults = -1) const {
    SubstructMatchParameters params;
    params.recursionPossible = recursionPossible;
    params.useChirality = useChirality;
    params.useQueryQueryMatches = useQueryQueryMatches;
    return getMatches(query, startIdx, endIdx, params, numThreads, maxResults);
  };
  //! overload
  std::vector<unsigned int> getMatches(const ROMol &query,
                                       unsigned int startIdx,
                                       unsigned int endIdx,
                                       const SubstructMatchParameters &params,
                                       int numThreads = -1,
                                       int maxResults = -1) const;
  //! overload
  std::vector<unsigned int> getMatches(const MolBundle &query,
                                       unsigned int startIdx,
                                       unsigned int endIdx,
                                       const SubstructMatchParameters &params,
                                       int numThreads = -1,
                                       int maxResults = -1) const;
  //! overload
  std::vector<unsigned int> getMatches(const TautomerQuery &query,
                                       unsigned int startIdx,
                                       unsigned int endIdx,
                                       const SubstructMatchParameters &params,
                                       int numThreads = -1,
                                       int maxResults = -1) const;

  //! Return the number of matches for the query
  /*!
    \param query       Molecule or Tautomer Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
                              [default true]
    \param useChirality  use atomic CIP codes as part of the comparison
                         [default true]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
                                 [default false]
    \param numThreads  If -1 use all available processors [default -1]
  */
  template <class Query>
  unsigned int countMatches(const Query &query, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1) const {
    SubstructMatchParameters params;
    params.recursionPossible = recursionPossible;
    params.useChirality = useChirality;
    params.useQueryQueryMatches = useQueryQueryMatches;
    return countMatches(query, 0, size(), params, numThreads);
  }
  //! overload
  template <class Query>
  unsigned int countMatches(const Query &query,
                            const SubstructMatchParameters &params,
                            int numThreads = -1) const {
    return countMatches(query, 0, size(), params, numThreads);
  }

  //! Return the number of matches for the query

  //! Return the number of matches for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (non-inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
                              [default true]
    \param useChirality  use atomic CIP codes as part of the comparison
                         [default true]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
                                 [default false]
    \param numThreads  If -1 use all available processors [default -1]
  */
  template <class Query>
  unsigned int countMatches(const Query &query, unsigned int startIdx,
                            unsigned int endIdx, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1) const {
    SubstructMatchParameters params;
    params.recursionPossible = recursionPossible;
    params.useChirality = useChirality;
    params.useQueryQueryMatches = useQueryQueryMatches;
    return countMatches(query, startIdx, endIdx, params, numThreads);
  };

  //! overload
  unsigned int countMatches(const ROMol &query, unsigned int startIdx,
                            unsigned int endIdx,
                            const SubstructMatchParameters &params,
                            int numThreads = -1) const;
  //! overload
  unsigned int countMatches(const TautomerQuery &query, unsigned int startIdx,
                            unsigned int endIdx,
                            const SubstructMatchParameters &params,
                            int numThreads = -1) const;
  //! overload
  unsigned int countMatches(const MolBundle &query, unsigned int startIdx,
                            unsigned int endIdx,
                            const SubstructMatchParameters &params,
                            int numThreads = -1) const;

  //! Returns true if any match exists for the query
  /*!
    \param query       Molecule or Tautomer Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
                              [default true]
    \param useChirality  use atomic CIP codes as part of the comparison
                         [default true]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries
                                 will be used as part of the matching
                                 [default false]
    \param numThreads  If -1 use all available processors [default -1]
  */
  template <class Query>
  bool hasMatch(const Query &query, bool recursionPossible = true,
                bool useChirality = true, bool useQueryQueryMatches = false,
                int numThreads = -1) const {
    SubstructMatchParameters params;
    params.recursionPossible = recursionPossible;
    params.useChirality = useChirality;
    params.useQueryQueryMatches = useQueryQueryMatches;
    return hasMatch(query, 0, size(), params, numThreads);
  }
  //! overload
  template <class Query>
  bool hasMatch(const Query &query, const SubstructMatchParameters &params,
                int numThreads = -1) const {
    return hasMatch(query, 0, size(), params, numThreads);
  }
  //! Returns true if any match exists for the query between the specified
  //! indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are
    allowed [default true] \param useChirality  use atomic CIP codes as part
    of the comparison [default true] \param useQueryQueryMatches  if set, the
    contents of atom and bond queries will be used as part of the matching
                                 [default false]
    \param numThreads  If -1 use all available processors [default -1]
  */
  template <class Query>
  bool hasMatch(const Query &query, unsigned int startIdx, unsigned int endIdx,
                bool recursionPossible = true, bool useChirality = true,
                bool useQueryQueryMatches = false, int numThreads = -1) const {
    SubstructMatchParameters params;
    params.recursionPossible = recursionPossible;
    params.useChirality = useChirality;
    params.useQueryQueryMatches = useQueryQueryMatches;
    return hasMatch(query, startIdx, endIdx, params, numThreads);
  };
  //! overload
  bool hasMatch(const ROMol &query, unsigned int startIdx, unsigned int endIdx,
                const SubstructMatchParameters &params,
                int numThreads = -1) const;
  //! overload
  bool hasMatch(const TautomerQuery &query, unsigned int startIdx,
                unsigned int endIdx, const SubstructMatchParameters &params,
                int numThreads = -1) const;
  //! overload
  bool hasMatch(const MolBundle &query, unsigned int startIdx,
                unsigned int endIdx, const SubstructMatchParameters &params,
                int numThreads = -1) const;
  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library (n.b. could contain
    null)
  */
  boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    // expects implementation to throw IndexError if out of range
    PRECONDITION(mols, "molholder is null in SubstructLibrary");
    return mols->getMol(idx);
  }

  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library (n.b. could contain
    null)
  */
  boost::shared_ptr<ROMol> operator[](unsigned int idx) {
    // expects implementation to throw IndexError if out of range
    PRECONDITION(mols, "molholder is null in SubstructLibrary");
    return mols->getMol(idx);
  }

  //! return the number of molecules in the library
  unsigned int size() const {
    PRECONDITION(mols, "molholder is null in SubstructLibrary");
    return rdcast<unsigned int>(molholder->size());
  }

  //! does error checking
  void setSearchOrder(const std::vector<unsigned int> &order) {
    for (const auto idx : order) {
      if (idx >= mols->size()) {
        throw IndexErrorException(idx);
      }
    }
    searchOrder = order;
  }

  const std::vector<unsigned int> &getSearchOrder() const {
    return searchOrder;
  }

  std::vector<unsigned int> &getSearchOrder() { return searchOrder; }
  //! access required for serialization
  void resetHolders() {
    is_tautomerquery = false;
    mols = molholder.get();
    fps = fpholder.get();
    if (fps && dynamic_cast<TautomerPatternHolder *>(fps) != nullptr) {
      is_tautomerquery = true;
    }
  }

  //! serializes (pickles) to a stream
  void toStream(std::ostream &ss) const;
  //! returns a string with a serialized (pickled) representation
  std::string Serialize() const;
  //! initializes from a stream pickle
  void initFromStream(std::istream &ss);
  //! initializes from a string pickle
  void initFromString(const std::string &text);
};
}  // namespace RDKit

#include "SubstructLibrarySerialization.h"
#endif
