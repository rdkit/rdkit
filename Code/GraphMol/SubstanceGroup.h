//
//
//  Copyright (C) 2018-2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file SubstanceGroup.h

  \brief Defines the SubstanceGroup class

*/
#include <RDGeneral/export.h>
#ifndef _RD_SGROUP_H
#define _RD_SGROUP_H

#include <unordered_map>

#include <Geometry/point.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDProps.h>
#include <boost/smart_ptr.hpp>

namespace RDKit {
class ROMol;
class RWMol;
class Bond;
class Atom;

//! used to indicate errors from incorrect sgroup access
class RDKIT_GRAPHMOL_EXPORT SubstanceGroupException
    : public std::runtime_error {
 public:
  //! construct with an error message
  SubstanceGroupException(const char *msg) : std::runtime_error(msg){};
  //! construct with an error message
  SubstanceGroupException(const std::string &msg) : std::runtime_error(msg){};
};

//! The class for representing SubstanceGroups
/*!
  <b>Notes:</b>
  - These are inspired by the SGroups in the MDL formats
  - Implementation is based on 2010 MDL SD specification:
    http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf
  - See SGroups.md for further, more comprehensive notes.

*/

class RDKIT_GRAPHMOL_EXPORT SubstanceGroup : public RDProps {
 public:
  //! Bond type (see V3000 spec)
  enum class BondType {
    XBOND,  // External/Crossing bond
    CBOND,  // Internal/Contained bond
  };

  typedef std::array<RDGeom::Point3D, 3> Bracket;

  //! Data structure for SAP lines (see V3000 spec)
  //! lvIdx may not be set; this signaled with value -1
  struct AttachPoint {
    unsigned int aIdx;
    int lvIdx;
    std::string id;
    bool operator==(const AttachPoint &other) const {
      return aIdx == other.aIdx && lvIdx == other.lvIdx && id == other.id;
    }
  };

  //! See specification for V3000 CSTATE
  //! vector may or not be considered, depending on TYPE
  struct CState {
    unsigned int bondIdx;
    RDGeom::Point3D vector;
    bool operator==(const CState &other) const {
      // note that we ignore coordinates for this
      return bondIdx == other.bondIdx;
    }
  };

  //! No default constructor
  #ifndef SWIG
   // Unfortunately, SWIG generated wrapper code uses temporary variables that require a default ctor not be deleted.
   SubstanceGroup() = delete;
  #endif // !SWIG

  //! Main Constructor. Ownership is only set on this side of the relationship:
  //! mol->addSubstanceGroup(sgroup) still needs to be called to get ownership
  //! on the other side.
  SubstanceGroup(ROMol *owning_mol, const std::string &type);

  SubstanceGroup(const SubstanceGroup &other) = default;
  SubstanceGroup(SubstanceGroup &&other) = default;

  SubstanceGroup &operator=(const SubstanceGroup &other) = default;
  SubstanceGroup &operator=(SubstanceGroup &&other) = default;

  //! Destructor
  ~SubstanceGroup(){};

  //! returns whether or not this belongs to a molecule
  bool hasOwningMol() const { return dp_mol != nullptr; };

  //! Get the molecule that owns this instance
  ROMol &getOwningMol() const {
    PRECONDITION(dp_mol, "no owner");
    return *dp_mol;
  }

  //! get the index of this sgroup in dp_mol's sgroups vector
  //! (do not mistake this by the ID!)00
  unsigned int getIndexInMol() const;

  /* Atom and Bond methods */
  void addAtomWithIdx(unsigned int idx);
  void addParentAtomWithIdx(unsigned int idx);
  void addBondWithIdx(unsigned int idx);
  void addAtomWithBookmark(int mark);
  void addParentAtomWithBookmark(int mark);
  void addBondWithBookmark(int mark);

  void addBracket(const Bracket &bracket);
  void addCState(unsigned int bondIdx, const RDGeom::Point3D &vector);
  void addAttachPoint(unsigned int aIdx, int lvIdx, const std::string &idStr);

  BondType getBondType(unsigned int bondIdx) const;

  const std::vector<unsigned int> &getAtoms() const { return d_atoms; }
  const std::vector<unsigned int> &getParentAtoms() const { return d_patoms; }
  const std::vector<unsigned int> &getBonds() const { return d_bonds; }

  const std::vector<Bracket> &getBrackets() const { return d_brackets; }
  const std::vector<CState> &getCStates() const { return d_cstates; }
  const std::vector<AttachPoint> &getAttachPoints() const { return d_saps; }

  //! adjusts our atom IDs to reflect that an atom has been removed from the
  //! parent molecule
  //!   decrements all atom IDs that are higher than \c atomIdx
  //!   raises a \c SubstanceGroupException if \c atomIdx is actually part of
  //!   this substance group
  //! \returns whether or not anything was changed
  bool adjustToRemovedAtom(unsigned int atomIdx);

  //! \returns whether or not anything the specified atom is part of the
  //! definition of this substance group
  bool includesAtom(unsigned int atomIdx) const;

  //! adjusts our bond IDs to reflect that a bond has been removed from the
  //! parent molecule
  //!   decrements all bond IDs that are higher than \c bondIdx
  //!   raises a \c SubstanceGroupException if \c bondIdx is actually part of
  //!   this substance group
  //! \returns whether or not anything was changed
  bool adjustToRemovedBond(unsigned int bondIdx);

  //! \returns whether or not anything the specified bond is part of the
  //! definition of this substance group
  bool includesBond(unsigned int bondIdx) const;

  //! Set owning molecule
  //! This only updates atoms and bonds; parent sgroup has to be updated
  //! independently, since parent might not exist at the time this is
  //! called.
  void setOwningMol(ROMol *mol);

  bool operator==(const SubstanceGroup &other) const {
    // we ignore brackets and cstates, which involve coordinates
    return dp_mol == other.dp_mol && d_atoms == other.d_atoms &&
           d_patoms == other.d_patoms && d_bonds == other.d_bonds &&
           d_saps == other.d_saps;
  }

 private:
  ROMol *dp_mol = nullptr;  // owning molecule

  std::vector<unsigned int> d_atoms;
  std::vector<unsigned int> d_patoms;
  std::vector<unsigned int> d_bonds;

  std::vector<Bracket> d_brackets;
  std::vector<CState> d_cstates;
  std::vector<AttachPoint> d_saps;
};

namespace SubstanceGroupChecks {

const std::vector<std::string> sGroupTypes = {
    // polymer sgroups:
    "SRU", "MON", "COP", "CRO", "GRA", "MOD", "MER", "ANY",
    // formulations/mixtures:
    "COM", "MIX", "FOR",
    // other
    "SUP", "MUL", "DAT", "GEN"};

const std::vector<std::string> sGroupSubtypes = {"ALT", "RAN", "BLO"};
const std::vector<std::string> sGroupConnectTypes = {"HH", "HT", "EU"};

RDKIT_GRAPHMOL_EXPORT bool isValidType(const std::string &type);

RDKIT_GRAPHMOL_EXPORT bool isValidSubType(const std::string &type);

RDKIT_GRAPHMOL_EXPORT bool isValidConnectType(const std::string &type);

RDKIT_GRAPHMOL_EXPORT bool isSubstanceGroupIdFree(const ROMol &mol,
                                                  unsigned int id);

}  // namespace SubstanceGroupChecks

//! \name SubstanceGroups and molecules
//@{

RDKIT_GRAPHMOL_EXPORT std::vector<SubstanceGroup> &getSubstanceGroups(
    ROMol &mol);
RDKIT_GRAPHMOL_EXPORT const std::vector<SubstanceGroup> &getSubstanceGroups(
    const ROMol &mol);

//! Add a new SubstanceGroup. A copy is added, so we can be sure that no other
//! references to the SubstanceGroup exist.
/*!
  \param sgroup - SubstanceGroup to be added to the molecule.
*/
RDKIT_GRAPHMOL_EXPORT unsigned int addSubstanceGroup(ROMol &mol,
                                                     SubstanceGroup sgroup);

//! Removes SubstanceGroups which reference a particular atom index
/*!
  \param mol - molecule to be edited.
  \param idx - atom index
*/
RDKIT_GRAPHMOL_EXPORT void removeSubstanceGroupsReferencingAtom(
    RWMol &mol, unsigned int idx);
//! Removes SubstanceGroups which reference a particular bond index
/*!
  \param mol - molecule to be edited.
  \param idx - bond index
*/
RDKIT_GRAPHMOL_EXPORT void removeSubstanceGroupsReferencingBond(
    RWMol &mol, unsigned int idx);
//@}

}  // namespace RDKit

//! allows SubstanceGroup objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::SubstanceGroup &sg);
#endif
