/*
* $Id: Atom.i 2519 2013-05-17 03:01:18Z glandrum $
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following
*       disclaimer in the documentation and/or other materials provided
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc.
*       nor the names of its contributors may be used to endorse or promote
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"

%{
#include <Query/QueryObjects.h>
#include <RDGeneral/types.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/Atom.h>
#include <GraphMol/ConjugHybrid.cpp>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/point.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

// from Python wrapper
std::string describeQueryHelper(const RDKit::Atom::QUERYATOM_QUERY *q, unsigned int depth) {
  std::string res;
  if (q) {
    for (unsigned int i = 0; i < depth; ++i) {
      res += "  ";
    }
    res += q->getFullDescription() + "\n";
    for (const auto &child :
         boost::make_iterator_range(q->beginChildren(), q->endChildren())) {
      res += describeQueryHelper(child.get(), depth + 1);
    }
  }
  return res;
}

%}

%ignore RDKit::Atom::Match(const Atom *) const;
%ignore RDKit::Atom::expandQuery;
%template(Bond_Vect) std::vector<RDKit::Bond*>;

%include "enums.swg"
#if swifjava
%javaconst(1);
#endif
%include <GraphMol/Atom.h>

%newobject RDKit::Atom::getProp;
%newobject RDKit::Atom::getBonds;

%extend RDKit::Atom {
  std::string getProp(const std::string key){
    std::string res;
    ($self)->getProp(key, res);
    return res;
  }

  /* Methods from ConjugHybrid.cpp */
  void markConjAtomBonds() {
    RDKit::markConjAtomBonds(($self));
  }
  int numBondsPlusLonePairs() {
    return RDKit::numBondsPlusLonePairs(($self));
  }
  bool atomHasConjugatedBond() {
    return RDKit::MolOps::atomHasConjugatedBond(($self));
  }
   /* From MolTransforms.h */
  void transformAtom(RDGeom::Transform3D &tform) {
  MolTransforms::transformAtom(($self), tform);
  }

  /* Based on Python wrappers and unit tests */
  bool IsInRing(){
    if(!($self)->getOwningMol().getRingInfo()->isInitialized()){
      RDKit::MolOps::findSSSR(($self)->getOwningMol());
    }
    return ($self)->getOwningMol().getRingInfo()->numAtomRings(($self)->getIdx())!=0;
  }

  bool IsInRingSize(int size){
    if(!($self)->getOwningMol().getRingInfo()->isInitialized()){
      RDKit::MolOps::findSSSR(($self)->getOwningMol());
    }
    return ($self)->getOwningMol().getRingInfo()->isAtomInRingOfSize(($self)->getIdx(),size);
  }

  std::vector<RDKit::Bond*> *getBonds() {
    std::vector<RDKit::Bond*> *bonds = new std::vector<RDKit::Bond*>;
    RDKit::ROMol *parent = &($self)->getOwningMol();
    RDKit::ROMol::OEDGE_ITER begin,end;
    boost::tie(begin,end) = parent->getAtomBonds(($self));
    while(begin!=end){
      RDKit::Bond *tmpB = (*parent)[*begin];
      bonds->push_back(tmpB);
      begin++;
    }
    return bonds;
  }

  // also matches ATOM_NULL_QUERY
  void setQuery(RDKit::ATOM_OR_QUERY *query) {
    $self->setQuery(query);
  }

  void setQuery(RDKit::ATOM_EQUALS_QUERY *query) {
    $self->setQuery(query);
  }

  // From Python Wrapper
  std::string AtomGetSmarts(bool doKekule=false, bool allHsExplicit=false,
                          bool isomericSmiles=true) {
	std::string res;
	if (($self)->hasQuery()) {
	  res = RDKit::SmartsWrite::GetAtomSmarts(static_cast<const RDKit::QueryAtom *>(($self)));
	} else {
	  // FIX: this should not be necessary
	  res = RDKit::SmilesWrite::GetAtomSmiles(($self), doKekule, nullptr, allHsExplicit,
									   isomericSmiles);
	}
	return res;
  }

  // from Python Wrapper
  std::string describeQuery() {
	std::string res = "";
	  res = describeQueryHelper(($self)->getQuery(), 0);
	return res;
  }



}
