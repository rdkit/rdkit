//
//  Copyright (C) 2018 T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file StereoGroup.h

  \brief Defines the class StereoGroup which stores relationships between
  the absolute configurations of atoms within a structure.

*/

#include <RDGeneral/export.h>
#ifndef RD_StereoGroup_092018
#define RD_StereoGroup_092018

#include <vector>

namespace RDKit {
class Atom;

// OR means that it is known to be one or the other, but not both
// AND means that it is known to be a mix.
enum class StereoGroupType {
  STEREO_ABSOLUTE = 0,
  STEREO_OR = 1,
  STEREO_AND = 2
};

//! StereoGroup is a collection of atoms with a known stereochemical
//! relationship
/*!
  Used to help represent a sample with unknown stereochemistry, or that is a mix
  of diastereomers.

 */
class RDKIT_GRAPHMOL_EXPORT StereoGroup {
private:
  StereoGroupType d_grouptype;
  std::vector<Atom *> d_atoms;

public:
  // Takes control of atoms if possible.
  StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms);
  StereoGroupType getGroupType() const;
  const std::vector<Atom *>& getAtoms() const;
};

void remove_groups_with_atom(const Atom* atom, std::vector<StereoGroup>& groups);
void remove_groups_with_atoms(const std::vector<Atom*>& atoms, std::vector<StereoGroup>& groups);


}  // namespace RDKit

#endif
