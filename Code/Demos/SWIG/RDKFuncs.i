// $Id$
//
// Copyright (C) 2008-2010 Greg Landrum
// All Rights Reserved
//
%include "std_string.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "std_pair.i"
%include "boost_shared_ptr.i"

%{
#include <vector>
#include <list>
#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <RDGeneral/Invariant.h>  
#include <RDBoost/Exceptions.h>  
#include <RDGeneral/types.h>
#include <GraphMol/MolOps.h>
#include "RDKFuncs.h"
%}


%include "exception.i"
%exception {
  try {
    $action
  } catch (RDKit::MolSanitizeException &e) {
    std::string msg="Sanitization error: ";
    msg += e.message();
    SWIG_exception(SWIG_ValueError,msg.c_str());
  } catch (Invar::Invariant &e) {
    std::string msg="Invariant error: "+e.getMessage();
    SWIG_exception(SWIG_RuntimeError,msg.c_str());
  } catch (IndexErrorException &e) {
    SWIG_exception(SWIG_IndexError,"bad index");
  } catch (ValueErrorException &e) {
    SWIG_exception(SWIG_ValueError,"bad value");
  } catch (KeyErrorException &e) {
    SWIG_exception(SWIG_ValueError,"bad key");
  } catch (...) {
    SWIG_exception(SWIG_RuntimeError,"Unknown exception");
  }
 }


%apply int { boost::int32_t };
%apply unsigned int { boost::uint32_t };
%apply long long { boost::int64_t };
%apply unsigned long long { boost::uint64_t };

%shared_ptr(RDKit::ROMol)
%shared_ptr(RDKit::RWMol)
%shared_ptr(RDKit::Atom)
%shared_ptr(RDKit::Bond)
#if SWIGCSHARP
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(ROMol, boost::shared_ptr<RDKit::ROMol> )
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(ROMol_Vect, std::vector< boost::shared_ptr<RDKit::ROMol> >);
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(Int_Vect, std::vector<int>)
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(Int_Pair, std::pair<int,int> );
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(Match_Vect, std::vector< std::pair<int,int> >);
#endif

typedef std::vector<int> INT_VECT;
typedef std::vector<std::vector<int> > INT_VECT;
typedef std::list<int > INT_LIST;

%template(ROMol_Vect) std::vector< boost::shared_ptr<RDKit::ROMol> >;
%template(Int_Vect) std::vector<int>;
%rename(equals) std::vector<int>::operator==;
%template(Int_Pair) std::pair<int, int >;
%template(Match_Vect) std::vector<std::pair<int,int> >;
%template(ROMol_Vect_Vect) std::vector< std::vector< boost::shared_ptr<RDKit::ROMol> > >;
%template(Int_Vect_Vect) std::vector<std::vector<int> >;
%template(Match_Vect_Vect) std::vector<std::vector<std::pair<int,int> > >;
%template(Char_Vect) std::vector<char>;

%ignore getAllAtomsWithBookmark;
%ignore getAtomBookmarks;
%ignore getAllBondsWithBookmark;
%ignore getBondBookmarks;
%ignore getAtomNeighbors;
%ignore getAtomBonds;
%ignore getAtomPMap;
%ignore getBondPMap;
%ignore getVertices;
%ignore getEdges;
%ignore getTopology;
%ignore debugMol;
%ignore beginAtoms;
%ignore endAtoms;
%ignore beginAromaticAtoms;
%ignore endAromaticAtoms;
%ignore beginHeteros;
%ignore endHeteros;
%ignore beginQueryAtoms;
%ignore endQueryAtoms;
%ignore beginBonds;
%ignore endBonds;
%ignore getPropList;
%ignore getConformer;
%ignore addConformer;
%ignore beginConformers;
%ignore endConformers;
%ignore RDKit::ROMol::getAtomDegree(const Atom *) const;
%ignore RDKit::ROMol::setAtomBookmark(Atom *,int);
%ignore RDKit::ROMol::clearAtomBookmark(const int, const Atom *);
%ignore RDKit::ROMol::setBondBookmark(Bond *,int);
%ignore RDKit::ROMol::clearBondBookmark(int, const Bond *);
%ignore RDKit::ROMol::getAtomWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondBetweenAtoms(unsigned int,unsigned int) const ;

%include <GraphMol/ROMol.h>
%extend RDKit::ROMol {
  bool hasSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    return SubstructMatch(*($self),query,mv,true,useChirality);
  };

  std::vector<std::pair<int, int> >
  getSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    SubstructMatch(*($self),query,mv,true,useChirality);
    return mv;

  };

  std::vector< std::vector<std::pair<int, int> > >
  getSubstructMatches(RDKit::ROMol &query,bool uniquify=true,
                      bool useChirality=false){
    std::vector<RDKit::MatchVectType> mv;
    SubstructMatch(*($self),query,mv,uniquify,true,useChirality);
    return mv;
  };

  std::string getProp(const std::string &propName){
    std::string res;
    $self->getProp(propName,res);
    return res;
  }

  void setProp(const std::string &propName,const std::string &propVal){
    $self->setProp(propName,propVal);
  }
}

%{
#include <GraphMol/RWMol.h>
%}
%ignore insertMol;
%ignore addAtom;
%ignore removeAtom;
%ignore addBond;
%ignore removeBond;
%ignore createPartialBond;
%ignore finishPartialBond;
%ignore replaceAtom;
%ignore getLastAtom;
%ignore getActiveAtom;
%include <GraphMol/RWMol.h>


%ignore copy;
%ignore setOwningMol;
%ignore setIdx;
%ignore setDativeFlag;
%ignore clearDativeFlag;
//%ignore hasQuery;
//%ignore setQuery;
//%ignore getQuery;
//%ignore expandQuery;
%ignore getPropList;
%ignore getPerturbationOrder;
%ignore RDKit::Atom::Match(const Atom *) const;
%include <GraphMol/Atom.h>
%extend RDKit::Atom {
  std::string getProp(const std::string &propName){
    std::string res;
    $self->getProp(propName,res);
    return res;
  }

  void setProp(const std::string &propName,const std::string &propVal){
    $self->setProp(propName,propVal);
  }
}


%ignore setIsConjugated;
%ignore setOwningMol;
%ignore setBeginAtom;
%ignore setEndAtom;
%ignore getStereoAtoms;
%ignore RDKit::Bond::getValenceContrib(const Atom *) const;
%ignore RDKit::Bond::Match(const Bond *) const;
%include <GraphMol/Bond.h>
%extend RDKit::Bond {
  std::string getProp(const std::string &propName){
    std::string res;
    $self->getProp(propName,res);
    return res;
  }

  void setProp(const std::string &propName,const std::string &propVal){
    $self->setProp(propName,propVal);
  }
}

%ignore initialize;
%ignore RDKit::RingInfo::isInitialized;
%ignore addRing;
%ignore reset;
%ignore preallocate;
%include <GraphMol/RingInfo.h>

%ignore RDKit::ChemicalReactionException::ChemicalReactionException(std::string const);
%ignore RDKit::ChemicalReaction::validate;
%include <GraphMol/ChemReactions/Reaction.h>
%extend RDKit::ChemicalReaction {
  bool validateReaction() const {
    unsigned int nErr=0,nWarn=0;
    bool res=$self->validate(nErr,nWarn);
    return res;
  }
}

%{
#include <GraphMol/FileParsers/MolSupplier.h>
%}

%include <GraphMol/FileParsers/MolSupplier.h>

%ignore ToSVG;
%include "RDKFuncs.h"

// ------------------- EBV fingerprints
%{
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
%}
#if SWIGCSHARP
%csmethodmodifiers ExplicitBitVect::ToString() const "public override";
#endif
%ignore ExplicitBitVect::dp_bits;
%include <DataStructs/ExplicitBitVect.h>
%newobject RDKFingerprintMol;
%newobject LayeredFingerprintMol;
%include <GraphMol/Fingerprints/Fingerprints.h>

// ------------------- SIV fingerprints
%{
#include <boost/cstdint.hpp>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
%}

#if SWIGCSHARP
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(UInt_Pair, std::pair<unsigned int,int> );
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(UInt_Pair_Vect, std::vector< std::pair<unsigned int,int> >);
#endif
%template(UInt_Pair) std::pair<unsigned int, int >;
%template(UInt_Pair_Vect) std::vector<std::pair<unsigned int,int> >;

%ignore RDKit::SparseIntVect<boost::uint32_t>::getNonzeroElements const;
%include <DataStructs/SparseIntVect.h>
%template(SparseIntVectu32) RDKit::SparseIntVect<boost::uint32_t>;
%extend RDKit::SparseIntVect<boost::uint32_t> {
  std::vector<std::pair<unsigned int, int> >
  getNonzero() const{
    std::vector<std::pair<unsigned int, int> > res;
    for(std::map<boost::uint32_t,int>::const_iterator es=$self->getNonzeroElements().begin();
        es!=$self->getNonzeroElements().end();++es){
      res.push_back(std::make_pair((unsigned int)es->first,(int)es->second));
    }
    return res;
  }
}
%template(UInt32_Vect) std::vector<boost::uint32_t>;
%newobject RDKit::MorganFingerprints::getFingerprint;
%rename(getMorganFingerprint) RDKit::MorganFingerprints::getFingerprint;
%newobject RDKit::MorganFingerprints::getFingerprintAsBitVect;
%rename(getMorganFingerprintAsBitVect) RDKit::MorganFingerprints::getFingerprintAsBitVect;
%include <GraphMol/Fingerprints/MorganFingerprints.h>


%template(SparseIntVecti32) RDKit::SparseIntVect<boost::int32_t>;
%template(SparseIntVecti64) RDKit::SparseIntVect<boost::int64_t>;
%newobject getAtomPairFingerprint;
%newobject getHashedAtomPairFingerprint;
%newobject getHashedAtomPairFingerprintAsBitVect;
%newobject getTopologicalTorsionFingerprint;
%newobject getHashedTopologicalTorsionFingerprint;
%newobject getHashedTopologicalTorsionFingerprintAsBitVect;
%include <GraphMol/Fingerprints/AtomPairs.h>


%include <DataStructs/BitOps.h>
%template(TanimotoSimilarityEBV) TanimotoSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(DiceSimilarityEBV) DiceSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(DiceSimilaritySIVu32) RDKit::DiceSimilarity<boost::uint32_t>;
%template(DiceSimilaritySIVi32) RDKit::DiceSimilarity<boost::int32_t>;
%template(DiceSimilaritySIVi64) RDKit::DiceSimilarity<boost::int64_t>;
%template(TanimotoSimilaritySIVu32) RDKit::TanimotoSimilarity<boost::uint32_t>;
%template(TanimotoSimilaritySIVi32) RDKit::TanimotoSimilarity<boost::int32_t>;
%template(TanimotoSimilaritySIVi64) RDKit::TanimotoSimilarity<boost::int64_t>;

%extend RDKit::ROMol {
  double MolLogP(RDKit::ROMol &mol){
    double logp,mr;
    RDKit::Descriptors::CalcCrippenDescriptors(mol,logp,mr);
    return logp;
  }
  double MolMR(RDKit::ROMol &mol){
    double logp,mr;
    RDKit::Descriptors::CalcCrippenDescriptors(mol,logp,mr);
    return mr;
  }
}

%{
#include <GraphMol/ChemTransforms/ChemTransforms.h>
%}
%newobject deleteSubstructs;
%newobject replaceSidechains;
%newobject replaceCore;
%newobject MurckoDecompose;
%include <GraphMol/ChemTransforms/ChemTransforms.h>

%{
#include <GraphMol/MolOps.h>
#include <RDGeneral/types.h>
%}

%include <RDGeneral/types.h>
%include <GraphMol/MolOps.h>


typedef std::vector<int> PATH_TYPE;
typedef std::list<std::vector<int> > PATH_LIST;
%template(Int_Vect_List) std::list<std::vector<int> >;
%template(Int_Int_Vect_List_Map) std::map<int,std::list<std::vector<int> > >;

%{
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
%}
%include <GraphMol/Subgraphs/Subgraphs.h>
%inline %{
  std::vector<int> calcPathDiscriminators(RDKit::ROMol &mol,RDKit::PATH_TYPE &path){
    std::vector<int> res(3);
    RDKit::Subgraphs::DiscrimTuple tpl=RDKit::Subgraphs::calcPathDiscriminators(mol,path);
    res[0]=boost::get<0>(tpl);
    res[1]=boost::get<1>(tpl);
    res[2]=boost::get<2>(tpl);
    return res;
  }
%}

%ignore calcPathDiscriminators;
%newobject pathToSubmol;
%include <GraphMol/Subgraphs/SubgraphUtils.h>

