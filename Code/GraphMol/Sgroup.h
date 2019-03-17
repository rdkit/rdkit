//
//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Sgroup.h

  \brief Defines the SGroup class

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
class Bond;
class Atom;

//! used to indicate errors from incorrect sgroup access
class RDKIT_GRAPHMOL_EXPORT SGroupException : public std::runtime_error {
 public:
  //! construct with an error message
  SGroupException(const char *msg) : std::runtime_error(msg){};
  //! construct with an error message
  SGroupException(const std::string &msg) : std::runtime_error(msg){};
};

//! The class for representing SGroups
/*!
  <b>Notes:</b>
  - Implementation is based on 2010 MDL SD specification:
    http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf
  - See SGroups.md for further, more comprehensive notes.

*/

class RDKIT_GRAPHMOL_EXPORT SGroup : public RDProps {
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
  SGroup() = delete;

  //! Main Constructor. Ownsership is only set on this side of the relationship:
  //! mol->addSGroup(sgroup) still needs to be called to get ownership on the
  //! other side.
  SGroup(ROMol *owning_mol, const std::string &type);

  SGroup(const SGroup &other) = default;
  SGroup(SGroup &&other) = default;

  SGroup &operator=(const SGroup &other) = default;
  SGroup &operator=(SGroup &&other) = default;

  //! Destructor
  ~SGroup(){};

  //! Get the molecule that owns this conformation
  ROMol &getOwningMol() const { return *dp_mol; }

  //! get the index of this sgroup in dp_mol's sgroups vector
  //! (do not mistake this by the ID!)
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

  //! Set owning molecule
  //! This only updates atoms and bonds; parent sgroup has to be updated
  //! independently, since parent might not exist at the time this is called.
  void setOwningMol(ROMol *mol);

  bool operator==(const SGroup &other) const {
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

namespace SGroupChecks {

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

RDKIT_GRAPHMOL_EXPORT bool isSGroupIdFree(const ROMol &mol, unsigned int id);

}  // namespace SGroupChecks

//! \name SGroups and molecules
//@{

RDKIT_GRAPHMOL_EXPORT std::vector<SGroup> &getSGroups(ROMol &mol);
RDKIT_GRAPHMOL_EXPORT const std::vector<SGroup> &getSGroups(const ROMol &mol);

//! Add a new SGroup. A copy is added, so we can be sure that no other
//! references to the SGroup exist.
/*!
  \param sgroup - SGroup to be added to the molecule.
*/
RDKIT_GRAPHMOL_EXPORT unsigned int addSGroup(ROMol &mol, SGroup sgroup);
//@}

}  // namespace RDKit

//! allows SGroup objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::SGroup &sg);
#endif
