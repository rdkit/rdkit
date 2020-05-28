//  Copyright (c) 2017-2019, Novartis Institutes for BioMedical Research Inc.
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
#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/MolOps.h>

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

  virtual unsigned int addMol(const ROMol &m) {
    mols.push_back(boost::make_shared<ROMol>(m));
    return size() - 1;
  }

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= mols.size()) throw IndexErrorException(idx);
    return mols[idx];
  }

  virtual unsigned int size() const {
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

  virtual unsigned int addMol(const ROMol &m) {
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

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= mols.size()) throw IndexErrorException(idx);
    boost::shared_ptr<ROMol> mol(new ROMol);
    MolPickler::molFromPickle(mols[idx], mol.get());
    return mol;
  }

  virtual unsigned int size() const {
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

  virtual unsigned int addMol(const ROMol &m) {
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

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= mols.size()) throw IndexErrorException(idx);

    boost::shared_ptr<ROMol> mol(SmilesToMol(mols[idx]));
    return mol;
  }

  virtual unsigned int size() const {
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

  virtual unsigned int addMol(const ROMol &m) {
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

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    if (idx >= mols.size()) throw IndexErrorException(idx);

    RWMol *m = SmilesToMol(mols[idx], 0, false);
    m->updatePropertyCache();
    return boost::shared_ptr<ROMol>(m);
  }

  virtual unsigned int size() const {
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
    for (size_t i = 0; i < fps.size(); ++i) delete fps[i];
  }

  virtual unsigned int size() const {
    return rdcast<unsigned int>(fps.size());
  }
  
  //! Adds a molecule to the fingerprinter
  unsigned int addMol(const ROMol &m) {
    fps.push_back(makeFingerprint(m));
    return rdcast<unsigned int>(fps.size() - 1);
  }

  //! Adds a raw bit vector to the fingerprinter
  unsigned int addFingerprint(const ExplicitBitVect &v) {
    fps.push_back(new ExplicitBitVect(v));
    return rdcast<unsigned int>(fps.size() - 1);
  }

  //! Return false if a substructure search can never match the molecule
  bool passesFilter(unsigned int idx, const ExplicitBitVect &query) const {
    if (idx >= fps.size()) throw IndexErrorException(idx);

    return AllProbeBitsMatch(query, *fps[idx]);
  }

  //! Get the bit vector at the specified index (throws IndexError if out of
  //! range)
  const ExplicitBitVect &getFingerprint(unsigned int idx) const {
    if (idx >= fps.size()) throw IndexErrorException(idx);
    return *fps[idx];
  }

  //! make the query vector
  //!  Caller owns the vector!
  virtual ExplicitBitVect *makeFingerprint(const ROMol &m) const = 0;

  std::vector<ExplicitBitVect *> &getFingerprints() { return fps; }
  const std::vector<ExplicitBitVect *> &getFingerprints() const { return fps; }
};

//! Uses the pattern fingerprinter to rule out matches
class RDKIT_SUBSTRUCTLIBRARY_EXPORT PatternHolder : public FPHolderBase {
 public:
  //! Caller owns the vector!
  virtual ExplicitBitVect *makeFingerprint(const ROMol &m) const {
    return PatternFingerprintMol(m, 2048);
  }
};

//! Substructure Search a library of molecules
/*!  This class allows for multithreaded substructure searches os
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

*/
class RDKIT_SUBSTRUCTLIBRARY_EXPORT SubstructLibrary {
  boost::shared_ptr<MolHolderBase> molholder;
  boost::shared_ptr<FPHolderBase> fpholder;
  MolHolderBase *mols;  // used for a small optimization
  FPHolderBase *fps{nullptr};

 public:
  SubstructLibrary()
      : molholder(new MolHolder),
        fpholder(),
        mols(molholder.get())
        {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules)
      : molholder(molecules), fpholder(), mols(molholder.get()), fps(nullptr) {}

  SubstructLibrary(boost::shared_ptr<MolHolderBase> molecules,
                   boost::shared_ptr<FPHolderBase> fingerprints)
      : molholder(molecules),
        fpholder(fingerprints),
        mols(molholder.get()),
        fps(fpholder.get()) {}

  SubstructLibrary(const std::string &pickle)
      : molholder(new MolHolder),
        fpholder(),
        mols(molholder.get()),
        fps(nullptr) {
    initFromString(pickle);
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

  const MolHolderBase &getMolecules() const {
    PRECONDITION(mols, "Molecule holder NULL in SubstructLibrary");
    return *mols;
  }

  //! Get the underlying fingerprint implementation.
  /*! Throws a value error if no fingerprints have been set */
  FPHolderBase &getFingerprints() {
    if (!fps)
      throw ValueErrorException("Substruct Library does not have fingerprints");
    return *fps;
  }

  const FPHolderBase &getFingerprints() const {
    if (!fps)
      throw ValueErrorException("Substruct Library does not have fingerprints");
    return *fps;
  }

  //! Add a molecule to the library
  /*!
    \param mol Molecule to add

    returns index for the molecule in the library
  */
  unsigned int addMol(const ROMol &mol);

  //! Get the matching indices for the query
  /*!
    \param query       Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all [default
    -1]
  */
  std::vector<unsigned int> getMatches(const ROMol &query,
                                       bool recursionPossible = true,
                                       bool useChirality = true,
                                       bool useQueryQueryMatches = false,
                                       int numThreads = -1,
                                       int maxResults = -1);
  //! Get the matching indices for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (non-inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all [default
    -1]
  */
  std::vector<unsigned int> getMatches(
      const ROMol &query, unsigned int startIdx, unsigned int endIdx,
      bool recursionPossible = true, bool useChirality = true,
      bool useQueryQueryMatches = false, int numThreads = -1,
      int maxResults = -1);

  //! Return the number of matches for the query
  /*!
    \param query       Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  unsigned int countMatches(const ROMol &query, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1);
  //! Return the number of matches for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (non-inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  unsigned int countMatches(const ROMol &query, unsigned int startIdx,
                            unsigned int endIdx, bool recursionPossible = true,
                            bool useChirality = true,
                            bool useQueryQueryMatches = false,
                            int numThreads = -1);

  //! Returns true if any match exists for the query
  /*!
    \param query       Query to match against molecules
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  bool hasMatch(const ROMol &query, bool recursionPossible = true,
                bool useChirality = true, bool useQueryQueryMatches = false,
                int numThreads = -1);
  //! Returns true if any match exists for the query between the specified
  //! indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (inclusive) of the search.
    \param recursionPossible  flags whether or not recursive matches are allowed
    [ default true ]
    \param useChirality  use atomic CIP codes as part of the comparison [
    default true ]
    \param useQueryQueryMatches  if set, the contents of atom and bond queries [
    default false ]
                                 will be used as part of the matching
    \param numThreads  If -1 use all available processors [default -1]
  */
  bool hasMatch(const ROMol &query, unsigned int startIdx, unsigned int endIdx,
                bool recursionPossible = true, bool useChirality = true,
                bool useQueryQueryMatches = false, int numThreads = -1);

  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library
  */
  boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    // expects implementation to throw IndexError if out of range
    PRECONDITION(mols, "molholder is null in SubstructLibrary");
    return mols->getMol(idx);
  }

  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library
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

  //! access required for serialization
  void resetHolders() {
    mols = molholder.get();
    fps = fpholder.get();
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
