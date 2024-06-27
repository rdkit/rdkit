//
// Created by Gareth Jones on 9/2/2020.
//
// Copyright 2020 Schrodinger, Inc
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.



%{
#include <GraphMol/SubstanceGroup.h>

const RDKit::SubstanceGroup *getSubstanceGroupWithIdx(RDKit::ROMol &mol, unsigned int idx) {
  auto &groups = RDKit::getSubstanceGroups(mol);
  return &(groups[idx]);
}

unsigned int getSubstanceGroupCount(RDKit::ROMol &mol) {
    return RDKit::getSubstanceGroups(mol).size();
}
%}


// Base class RDProps is wrapped with shared_ptr, so SubstanceGroup must be too.
%shared_ptr(RDKit::SubstanceGroup)
%ignore getSubstanceGroups;

RDKit::SubstanceGroup *getSubstanceGroupWithIdx(RDKit::ROMol &mol, unsigned int idx);
unsigned int getSubstanceGroupCount(RDKit::ROMol &mol);

%ignore RDKit::SubstanceGroup::getAtoms;
%rename(getAtoms) RDKit::SubstanceGroup::getSgAtoms;
%ignore RDKit::SubstanceGroup::getBonds;
%rename(getBonds) RDKit::SubstanceGroup::getSgBonds;
%ignore RDKit::SubstanceGroup::getBrackets;
%ignore RDKit::SubstanceGroup::getCStates;
%ignore RDKit::SubstanceGroup::getAttachPoints;

%include <GraphMol/SubstanceGroup.h>

%extend RDKit::SubstanceGroup {
  // Wrap getAtoms, getParentAtoms and getBonds to return vector<int> by value (not reference).
  // int instead of unsigned int as vector<unsigned int> is not wrapped in C#
  // Not using references to avoid a free() error from Java
  const std::vector<int> getSgAtoms() {
    auto atoms = $self->getAtoms();
    std::vector<int> atomIdxs(atoms.begin(), atoms.end());
    return atomIdxs;
  }
  const std::vector<int> getSgBonds() {
    auto bonds = $self->getBonds();
    std::vector<int> bondIdxs(bonds.begin(), bonds.end());
    return bondIdxs;
  }
  const std::vector<int> getSgParentAtoms() {
    auto atoms = $self->getParentAtoms();
    std::vector<int> atomIdxs(atoms.begin(), atoms.end());
    return atomIdxs;
  }

  // SWIG does not do well wrapping vectors of objects, so provide accessors
  const RDKit::SubstanceGroup::Bracket *getBracket(unsigned int idx)  {
    auto &brackets = $self->getBrackets();
    return &(brackets[idx]);
  }
  const RDKit::SubstanceGroup::CState *getCState(unsigned int idx)  {
    auto &cstates = $self->getCStates();
    return &(cstates[idx]);
  }
  const RDKit::SubstanceGroup::AttachPoint *getAttachPoint(unsigned int idx)  {
    auto &attachPoints = $self->getAttachPoints();
    return &(attachPoints[idx]);
  }
  size_t getBracketCount() {
    return $self->getBrackets().size();
  }
  size_t getCStateCount() {
    return $self->getCStates().size();
  }
  size_t getAttachPointCount() {
    return $self->getAttachPoints().size();
  }
}

%template(getStringProp) RDKit::SubstanceGroup::getProp<std::string>;
%template(getUIntProp) RDKit::SubstanceGroup::getProp<unsigned int>;
%template(getStringVectProp) RDKit::SubstanceGroup::getProp<RDKit::STR_VECT>;
%template(getUIntVectProp) RDKit::SubstanceGroup::getProp<RDKit::UINT_VECT>;

