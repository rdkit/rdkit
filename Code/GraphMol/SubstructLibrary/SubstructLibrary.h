//  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
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
#ifndef RDKIT_SUBSTRUCT_LIBRARY
#define RDKIT_SUBSTRUCT_LIBRARY

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>

namespace RDKit {

//! Base class API for holding molecules so substructure search.
/*!
  This is an API that hides the implementation details used for
  indexing molecules for substructure searching.  It simply
  provides an API for adding and getting molecules from a set.
 */
class MolHolderBase {
public:
  virtual ~MolHolderBase() {}

  //! Add a new molecule to the substructure search library
  //!  Returns the molecules index in the library
  virtual unsigned int addMol( const ROMol &m ) = 0;
  
  virtual boost::shared_ptr<ROMol> getMol(unsigned int) const = 0;
     // throws IndexError?

  //! Get the current library size
  virtual unsigned int size() const = 0;
};

//! Concrete class that holds molecules in memory
/*!
    This is currently one of the faster implementations.
    However it is very memory intensive.
*/
class MolHolder: public MolHolderBase {
  std::vector<boost::shared_ptr<ROMol> > mols;
public:
  MolHolder() : MolHolderBase(), mols() {}

  virtual unsigned int addMol( const ROMol &m ) {
    mols.push_back( boost::make_shared<ROMol>(m) );
    return size()-1;
  }
    
  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    // throws IndexError?
    return mols[idx];
  }
  
  virtual unsigned int size() const { return rdcast<unsigned int>(mols.size()); }

};

//! Concrete class that holds binary cached molecules in memory
/*!
  This implementation uses quite a bit less memory than the
  non cached implementation.  However, due to the reduced speed
  it should be used in conjunction with a pattern fingerprinter.

  See RDKit::FPHolder
*/
class CachedMolHolder: public MolHolderBase {
  std::vector<std::string> mols;
public:
  CachedMolHolder() : MolHolderBase(), mols() {}

  virtual unsigned int addMol( const ROMol &m ) {
    mols.push_back( std::string() );
    MolPickler::pickleMol(m, mols.back());
    return size()-1;
  }

  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    // throws IndexError?
    boost::shared_ptr<ROMol> mol(new ROMol);
    MolPickler::molFromPickle(mols[idx], mol.get());
    return mol;
  }
  
  virtual unsigned int size() const { return rdcast<unsigned int>(mols.size()); }
};

//! Concrete class that holds smiles strings in memory
/*!
    This implementation uses quite a bit less memory than the
    cached binary or uncached implementation.  However, due to the
    reduced speed it should be used in conjunction with a pattern
    fingerprinter.

    See RDKit::FPHolder
*/
class CachedSmilesMolHolder: public MolHolderBase {
  std::vector<std::string> mols;
public:
  CachedSmilesMolHolder() : MolHolderBase(), mols() {}

  virtual unsigned int addMol( const ROMol &m ) {

    bool doIsomericSmiles = true;
    mols.push_back(MolToSmiles(m, doIsomericSmiles));
    return size()-1;
  }
  
  virtual boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    // throws IndexError?
    boost::shared_ptr<ROMol> mol(SmilesToMol(mols[idx]));
    return mol;
  }
  
  virtual unsigned int size() const { return rdcast<unsigned int>(mols.size()); }

};


//! Base FPI for the fingerprinter used to rule out impossible matches
class FPHolderBase {
  std::vector<ExplicitBitVect *> fps;

public:
 virtual ~FPHolderBase() {
   for(size_t i=0;i<fps.size();++i)
     delete fps[i];
   
 }

 //! Add a molecule to the fingerprinter
 unsigned int addMol( const ROMol &m) {
   fps.push_back(getQueryBits(m));
   return rdcast<unsigned int>(fps.size() - 1);
 }

 //! Return false if a substructure search can never match the molecule
 bool passesFilter(unsigned int idx, const ExplicitBitVect &query) const {
   PRECONDITION(idx < fps.size(), "idx out of range in PatternHolder::passesFilter");
   return AllProbeBitsMatch(query, *fps[idx]);
 }
 
 virtual ExplicitBitVect *getQueryBits(const ROMol &m) const = 0;
 
};

//! Uses the pattern fingerprinter to rule out matches
class PatternHolder : public FPHolderBase {
public:
  virtual ExplicitBitVect *getQueryBits(const ROMol &m) const {
    return PatternFingerprintMol(m, 2048);
  }
};

//! Substtructure Search a library of molecules
/*!  This class allows for multithreaded substructure searches os
     large datasets.

     The implementations can use fingerprints to speed up searches
     and have molecules cached as binary forms to reduce memory
     usage.
*/
class SubstructLibrary {
  boost::shared_ptr<MolHolderBase> molholder;
  boost::shared_ptr<FPHolderBase> fpholder;
  MolHolderBase *mols; // used for a small optimization
  FPHolderBase *fps;
public:
  SubstructLibrary() : molholder( new MolHolder ), fpholder() {}
  
  SubstructLibrary( boost::shared_ptr<MolHolderBase> molecules ) :
       molholder(molecules), fpholder(), mols(molholder.get()), fps(0) {
  }

  SubstructLibrary( boost::shared_ptr<MolHolderBase>  molecules,
                    boost::shared_ptr<FPHolderBase> fingerprints ) :
       molholder(molecules), fpholder(fingerprints),
       mols(molholder.get()), fps(fpholder.get()) {
  }
      

  //!Get the underlying molecule holder implementation
        MolHolderBase & getMolHolder() { return *mols; }
  const MolHolderBase & getMolecules() const { return *mols; }

  //!Get the underlying fingerprint implementation.
  /*! Throws a value error if no fingerprints have been set */
  FPHolderBase & getFingerprints() {
    PRECONDITION(fps, "No fingerprints set in library");
    return *fps;
  }
  const FPHolderBase & getFingerprints() const {
    PRECONDITION(fps, "No fingerprints set in library");
    return *fps;
  }

  //!Add a molecule to the library
  /*!
    \param mol Molecule to add
  */
  unsigned int addMol(const ROMol &mol);

  //!Get the matching indices for the query
  /*!
    \param query       Query to match against molecules
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all [default -1]    
  */
  std::vector<unsigned int> getMatches(const ROMol &query, int numThreads=-1,
                                       int maxResults=-1);
  //!Get the matching indices for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (inclusive) of the search.
    \param numThreads  If -1 use all available processors [default -1]
    \param maxResults  Maximum results to return, -1 means return all [default -1]
  */  
  std::vector<unsigned int> getMatches(const ROMol &query,
                                       unsigned int startIdx, unsigned int endIdx,
                                       int numThreads=-1,
                                       int maxResults=-1);

  //!Return the number of matches for the query
  /*!
    \param query       Query to match against molecules
    \param numThreads  If -1 use all available processors [default -1]
  */  
  unsigned int countMatches(const ROMol &query, int numThreads=-1);
  //!Return the number of matches for the query between the given indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (inclusive) of the search.
    \param numThreads  If -1 use all available processors [default -1]
  */    
  unsigned int countMatches(const ROMol &query,
                            unsigned int startIdx, unsigned int endIdx, int numThreads=-1);

  //! Returns true if any match exists for the query
  /*!
    \param query       Query to match against molecules
    \param numThreads  If -1 use all available processors [default -1]
  */  
  bool hasMatch(const ROMol &query, int numThreads=-1);
  //! Returns true if any match exists for the query between the specified indices
  /*!
    \param query       Query to match against molecules
    \param startIdx    Start index of the search
    \param endIdx      Ending idx (inclusive) of the search.
    \param numThreads  If -1 use all available processors [default -1]
  */      
  bool hasMatch(const ROMol &query,
                unsigned int startIdx, unsigned int endIdx, int numThreads=-1);


  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library
  */
  boost::shared_ptr<ROMol> getMol(unsigned int idx) const {
    return mols->getMol(idx);
  }

  //! Returns the molecule at the given index
  /*!
    \param idx       Index of the molecule in the library
  */  
  boost::shared_ptr<ROMol> operator[] (unsigned int idx) {
    return mols->getMol(idx);
  }

  //! return the number of molecules in the library
  unsigned int size() const {
    return rdcast<unsigned int>(molholder->size());
  }
};
}

#endif
