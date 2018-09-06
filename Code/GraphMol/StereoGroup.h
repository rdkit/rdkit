//
//  Copyright (C) 2003-2015 Greg Landrum and Rational Discovery LLC
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
#ifndef __RD_StereoGroup_H__
#define __RD_StereoGroup_H__

#include <vector>

namespace RDKit {
class Atom;

// OR means that it is known to be one or the other, but not both
// AND means that it is known to be a mix.
enum class StereoGroupType : int { STEREO_ABSOLUTE = 0, STEREO_OR = 1, STEREO_AND = 2 };

//! StereoGroup is a collection of atoms with a known stereochemical
//! relationship
/*!
  Used to help represent a sample with unknown stereochemistry, or that is a mix
  of diastereomers.

 */
class RDKIT_GRAPHMOL_EXPORT StereoGroup {
 public:
  StereoGroup(StereoGroupType grouptype, std::vector<Atom*>&& atoms)
      : grouptype(grouptype), atoms(atoms) {}

  const StereoGroupType grouptype;
  const std::vector<Atom*> atoms;
};

}  // namespace RDKit

#endif
