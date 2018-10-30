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

class RDKIT_GRAPHMOL_EXPORT SGroup {
 public:
  friend class ROMol;

  enum class BondType {
    XBOND,  // External
    CBOND,  // Internal
  };

  typedef std::array<RDGeom::Point3D, 3> Bracket;

  struct AttachPoint {
    Atom *aAtom;
    Atom *lvAtom;
    std::string id;
  };

  struct CState {
    Bond *bond;
    boost::shared_ptr<RDGeom::Point3D> vector;
  };

  typedef std::vector<Atom *> ATOM_PTR_VECT;
  typedef ATOM_PTR_VECT::iterator ATOM_PTR_VECT_I;
  typedef ATOM_PTR_VECT::const_iterator ATOM_PTR_VECT_CI;
  typedef std::vector<Bond *> BOND_PTR_VECT;
  typedef BOND_PTR_VECT::iterator BOND_PTR_VECT_I;
  typedef BOND_PTR_VECT::const_iterator BOND_PTR_VECT_CI;

  //! No default constructor
  SGroup() = delete;

  //! Constructor. Ownsership is only set on this side of the
  //! relationship: mol->addSGroup(sgroup) still needs to be
  //! called to get ownership on the other side.
  SGroup(ROMol *owning_mol, const std::string &_type);

  //! Copy Constructor: initialize from a second sgroup.
  SGroup(const SGroup &other) = default;

  //! Destructor
  ~SGroup(){};

  //! Get the molecule that owns this conformation
  ROMol &getOwningMol() const { return *dp_mol; }

  //! get the ID of this sgroup
  inline unsigned int getId() const { return d_id; }

  //! get the COMPNO of this sgroup
  inline unsigned int getCompNo() const { return d_compno; }

  //! get the index of this sgroup in dp_mol's sgroups vector
  //! (do not mistake this by the ID!)
  unsigned int getIndexInMol() const;

  //! get the type of the SGroup
  inline const std::string &getType() const { return d_type; }

  //! check if SGroup has the given property set
  bool hasProp(const std::string &prop) const;

  //! get string properties
  const std::string &getProp(const std::string &prop) const;

  //! get parent SGroup
  inline SGroup *getParent() const { return d_parent; }

  //! set the ID of this sgroup
  void setId(unsigned int id);

  //! set the COMPNO of this sgroup
  void setCompNo(unsigned int compno) { d_compno = compno; }

  //! set the type of the SGroup
  inline void setType(const std::string &type) { d_type = type; }

  //! set string properties
  inline void setProp(const std::string &prop, const std::string &value) {
    d_prop[prop] = value;
  };

  //! set the parent SGroup
  inline void setParent(SGroup *parent) { d_parent = parent; }

  /* Atom and Bond methods */
  void addAtomWithIdx(unsigned int idx);
  void addPAtomWithIdx(unsigned int idx);
  void addBondWithIdx(unsigned int idx);
  void addAtomWithBookmark(int mark);
  void addPAtomWithBookmark(int mark);
  void addBondWithBookmark(int mark);
  void addBracket(const Bracket &bracket);
  void addCState(Bond *bond, RDGeom::Point3D *vector);
  void addAttachPoint(Atom *aAtomPtr, Atom *lvAtomPtr, std::string idStr);
  void addDataField(const std::string &data);

  BondType getBondType(Bond *bond) const;

  inline const ATOM_PTR_VECT &getAtoms() const { return d_atoms; };
  inline const ATOM_PTR_VECT &getPAtoms() const { return d_patoms; };
  inline const BOND_PTR_VECT &getBonds() const { return d_bonds; };
  inline const std::vector<Bracket> &getBrackets() const { return d_brackets; }
  inline const std::vector<CState> &getCStates() const { return d_cstates; }
  inline const std::vector<std::string> &getDataFields() const {
    return d_dataFields;
  }
  inline const std::vector<AttachPoint> &getAttachPoints() const {
    return d_saps;
  }
  inline const std::unordered_map<std::string, std::string> &getProps() const {
    return d_prop;
  }

 protected:
  //! Set owning moelcule
  //! This only updates atoms and bonds; parent sgroup has to be updated
  //! independently, since parent might not exist at the time this is called.
  void setOwningMol(ROMol *mol);

 private:
  //! Update pointers (atoms, patoms, bonds, parent SGroup) with the ones
  //! at matching indexes in other_mol. Must be called only after all of
  //! these objects exist in other_mol!
  void updateOwningMol(ROMol *other_mol);

  /* ID of the group. If not 0, must be unique */
  unsigned int d_id = 0;
  unsigned int d_compno = 0;

  std::string d_type;       // type of the sgroup
  ROMol *dp_mol = nullptr;  // owning molecule
  SGroup *d_parent = nullptr;

  ATOM_PTR_VECT d_atoms;
  ATOM_PTR_VECT d_patoms;
  BOND_PTR_VECT d_bonds;
  std::vector<Bracket> d_brackets;
  std::vector<CState> d_cstates;
  std::vector<std::string> d_dataFields;
  std::vector<AttachPoint> d_saps;

  std::unordered_map<std::string, std::string>
      d_prop;  // Storage for string properties
};

typedef boost::shared_ptr<SGroup> SGROUP_SPTR;

bool SGroupTypeOK(std::string typ);

bool SGroupSubTypeOK(std::string typ);

bool SGroupConnectTypeOK(std::string typ);

}  // namespace RDKit

//! allows SGroup objects to be dumped to streams
RDKIT_GRAPHMOL_EXPORT std::ostream &operator<<(std::ostream &target,
                                               const RDKit::SGroup &sg);
#endif
