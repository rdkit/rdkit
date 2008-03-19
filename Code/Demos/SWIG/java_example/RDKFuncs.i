// $Id$
//
// Copyright (C) 2008 Greg Landrum
// All Rights Reserved
//
%module RDKFuncs
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"
%include "boost_shared_ptr.i"

%{
#include <vector>
#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Atom.h>
#include "RDKFuncs.h"

%}

SWIG_SHARED_PTR(ROMol, RDKit::ROMol)
SWIG_SHARED_PTR(Atom, RDKit::Atom)
SWIG_SHARED_PTR(Bond, RDKit::Bond)
%template(ROMol_Vect) std::vector< boost::shared_ptr<RDKit::ROMol> >;
%template(ROMol_Vect_Vect) std::vector< std::vector< boost::shared_ptr<RDKit::ROMol> > >;
%template(Int_Vect) std::vector<int>;
%template(Int_Vect_Vect) std::vector<std::vector<int> >;
%template(Int_Pair) std::pair<int, int >;
%template(Match_Vect) std::vector<std::pair<int,int> >;
%template(Match_Vect_Vect) std::vector<std::vector<std::pair<int,int> > >;


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
%ignore RDKit::ROMol::hasProp(std::string const) const ;
%ignore RDKit::ROMol::clearProp(std::string const) const ;
%ignore RDKit::ROMol::getAtomWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondWithIdx(unsigned int) const ;
%ignore RDKit::ROMol::getBondBetweenAtoms(unsigned int,unsigned int) const ;

%include <GraphMol/ROMol.h>
%extend RDKit::ROMol {
  bool hasSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    return SubstructMatch(*($self),query,mv,true,useChirality,false);
  };

  std::vector<std::pair<int, int> >
  getSubstructMatch(RDKit::ROMol &query,bool useChirality=false){
    RDKit::MatchVectType mv;
    SubstructMatch(*($self),query,mv,true,useChirality,false);
    return mv;
  };

  std::vector< std::vector<std::pair<int, int> > >
  getSubstructMatches(RDKit::ROMol &query,bool uniquify=true,
                      bool useChirality=false){
    std::vector<RDKit::MatchVectType> mv;
    SubstructMatch(*($self),query,mv,uniquify,true,useChirality,false);
    return mv;
  };

}

%ignore copy;
%ignore setOwningMol;
%ignore setIdx;
%ignore setFormalCharge;
%ignore setNoImplicit;
%ignore setNumExplicitHs;
%ignore setIsAromatic;
%ignore setMass;
%ignore setDativeFlag;
%ignore clearDativeFlag;
%ignore setChiralTag;
%ignore setHybridizationType;
%ignore hasQuery;
%ignore setQuery;
%ignore getQuery;
%ignore expandQuery;
%ignore getPropList;
%ignore getPerturbationOrder;
%ignore RDKit::Atom::Match(const Atom *) const;
%ignore clearProp(std::string const) const;
%ignore hasProp(std::string const) const;
%include <GraphMol/Atom.h>

%ignore setIsConjugated;
%ignore setOwningMol;
%ignore setBeginAtom;
%ignore setEndAtom;
%ignore setBondDir;
%ignore setStereo;
%ignore getStereoAtoms;
%ignore RDKit::Bond::getValenceContrib(const Atom *) const;
%ignore RDKit::Bond::Match(const Bond *) const;
%include <GraphMol/Bond.h>

%ignore initialize;
%ignore isInitialized;
%ignore addRing;
%ignore reset;
%ignore preallocate;
%include <GraphMol/RingInfo.h>

%ignore RDKit::ChemicalReactionException::ChemicalReactionException(std::string const);
%include <GraphMol/ChemReactions/Reaction.h>

%include "RDKFuncs.h"


