//
//  Copyright (C) 2018-2021 Greg Landrum and other RDKit contributors
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

#include <iostream>
#include <vector>

namespace RDKit {
class Atom;
class Bond;
class ROMol;

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
  StereoGroupType d_grouptype{StereoGroupType::STEREO_ABSOLUTE};
  std::vector<Atom *> d_atoms;
  std::vector<Bond *> d_bonds;

  // The group ID for AND/OR groups (it has no meaning in ABS groups).
  // 0 means no group ID is defined.
  unsigned d_readId = 0u;
  unsigned d_writeId = 0u;

 public:
  StereoGroup() {}
  // Takes control of atoms if possible.
  StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms,
              std::vector<Bond *> &&bonds, unsigned readId = 0);
  StereoGroup(StereoGroupType grouptype, const std::vector<Atom *> &atoms,
              const std::vector<Bond *> &bonds, unsigned readId = 0);

  StereoGroup(const StereoGroup &other) = default;
  StereoGroup &operator=(const StereoGroup &other) = default;
  StereoGroup(StereoGroup &&other) = default;
  StereoGroup &operator=(StereoGroup &&other) = default;

  StereoGroupType getGroupType() const;
  const std::vector<Atom *> &getAtoms() const;
  const std::vector<Bond *> &getBonds() const;

  unsigned getReadId() const { return d_readId; }
  unsigned getWriteId() const { return d_writeId; }
  void setWriteId(unsigned id) { d_writeId = id; }

  // Seems odd to have to define these, but otherwise the SWIG wrappers
  // won't build
  bool operator==(const StereoGroup &other) const {
    return (d_grouptype == other.d_grouptype) && (d_atoms == other.d_atoms) &&
           (d_bonds == other.d_bonds);
  }
  bool operator!=(const StereoGroup &other) const {
    return (d_grouptype != other.d_grouptype) || (d_atoms != other.d_atoms) ||
           (d_bonds != other.d_bonds);
  }
  friend RDKIT_GRAPHMOL_EXPORT void removeAtomFromGroups(
      const Atom *atom, std::vector<StereoGroup> &groups);
};
RDKIT_GRAPHMOL_EXPORT void removeAtomFromGroups(
    const Atom *atom, std::vector<StereoGroup> &groups);
RDKIT_GRAPHMOL_EXPORT void removeGroupsWithAtom(
    const Atom *atom, std::vector<StereoGroup> &groups);
RDKIT_GRAPHMOL_EXPORT void removeGroupsWithAtoms(
    const std::vector<Atom *> &atoms, std::vector<StereoGroup> &groups);

//! Assign Group output IDs to all AND and OR StereoGroups in the vector
//! that don't already have one. The IDs are assigned based on the order
//! of the groups.
RDKIT_GRAPHMOL_EXPORT void assignStereoGroupIds(
    std::vector<StereoGroup> &groups);

//! Copy StereoGroup "read" IDs to "write" IDs so that they will be preserved
//! when the mol is exported.
RDKIT_GRAPHMOL_EXPORT void forwardStereoGroupIds(ROMol &mol);

}  // namespace RDKit

//! allows StereoGroup objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::StereoGroup &stg);

#endif
