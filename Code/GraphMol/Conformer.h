//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_CONFORMER_H
#define _RD_CONFORMER_H

#include <Geometry/point.h>
#include <RDGeneral/types.h>
#include <boost/smart_ptr.hpp>
#include <RDGeneral/RDProps.h>

namespace RDKit {
class ROMol;

//! used to indicate errors from incorrect conformer access
class RDKIT_GRAPHMOL_EXPORT ConformerException : public std::exception {
 public:
  //! construct with an error message
  ConformerException(const char *msg) : _msg(msg){};
  //! construct with an error message
  ConformerException(const std::string &msg) : _msg(msg){};
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); };
  ~ConformerException() noexcept {};

 private:
  std::string _msg;
};

//! The class for representing 2D or 3D conformation of a molecule
/*!
  This class contains
  - a pointer to the owing molecule
  - a vector of 3D points (positions of atoms)
*/
class RDKIT_GRAPHMOL_EXPORT Conformer : public RDProps {
 public:
  friend class ROMol;

  //! Constructor
  Conformer() { d_positions.clear(); };

  //! Constructor with number of atoms specified ID specification
  Conformer(unsigned int numAtoms) {
    if (numAtoms) {
      d_positions.resize(numAtoms, RDGeom::Point3D(0.0, 0.0, 0.0));
    } else {
      d_positions.resize(0);
      d_positions.clear();
    }
  };

  //! Copy Constructor: initialize from a second conformation.
  Conformer(const Conformer &other);
  Conformer &operator=(const Conformer &other);

  //! Destructor
  ~Conformer(){};

  //! Resize the conformer so that more atoms location can be added.
  //! Useful, for e.g., when adding hydrogens
  void resize(unsigned int size) { d_positions.resize(size); }

  //! Reserve more space for atom position
  void reserve(unsigned int size) { d_positions.reserve(size); }

  //! returns whether or not this instance belongs to a molecule
  bool hasOwningMol() const { return dp_mol != nullptr; };

  //! Get the molecule that owns this instance
  ROMol &getOwningMol() const {
    PRECONDITION(dp_mol, "no owner");
    return *dp_mol;
  }

  //! Get a const reference to the vector of atom positions
  const RDGeom::POINT3D_VECT &getPositions() const;

  //! Get a reference to the atom positions
  RDGeom::POINT3D_VECT &getPositions();

  //! Get the position of the specified atom
  const RDGeom::Point3D &getAtomPos(unsigned int atomId) const;
  //! overload
  template <class U>
  const RDGeom::Point3D &getAtomPos(U atomId) const {
    return getAtomPos(rdcast<unsigned int>(atomId));
  }

  //! Get the position of the specified atom
  RDGeom::Point3D &getAtomPos(unsigned int atomId);
  //! overload
  template <class U>
  RDGeom::Point3D &getAtomPos(U atomId) {
    return getAtomPos(rdcast<unsigned int>(atomId));
  }

  //! Set the position of the specified atom
  inline void setAtomPos(unsigned int atomId, const RDGeom::Point3D &position) {
    // RANGE_CHECK(0,atomId,d_positions.size()-1);
    if (atomId >= d_positions.size()) {
      d_positions.resize(atomId + 1, RDGeom::Point3D(0.0, 0.0, 0.0));
    }
    d_positions[atomId] = position;
  }
  //! overload
  template <class U>
  void setAtomPos(U atomId, const RDGeom::Point3D &position) {
    return setAtomPos(rdcast<unsigned int>(atomId), position);
  }
  //! get the ID of this conformer
  inline unsigned int getId() const { return d_id; }

  //! set the ID of this conformer
  inline void setId(unsigned int id) { d_id = id; }

  //! Get the number of atoms
  inline unsigned int getNumAtoms() const {
    return rdcast<unsigned int>(d_positions.size());
  }
  inline bool is3D() const { return df_is3D; }
  inline void set3D(bool v) { df_is3D = v; }

 protected:
  //! Set owning molecule
  void setOwningMol(ROMol *mol);

  //! Set owning molecule
  void setOwningMol(ROMol &mol);

 private:
  bool df_is3D{true};                // is this a 3D conformation?
  unsigned int d_id{0};              // id is the conformation
  ROMol *dp_mol{nullptr};            // owning molecule
  RDGeom::POINT3D_VECT d_positions;  // positions of the atoms
  void initFromOther(const Conformer &conf);
};

typedef boost::shared_ptr<Conformer> CONFORMER_SPTR;

//! Returns true if any of the z coords are non zero, false otherwise
/*!
  \param conf  Conformer object to analyze
*/
inline bool hasNonZeroZCoords(const Conformer &conf) {
  for (auto p : conf.getPositions()) {
    if (p.z != 0.0) return true;
  }
  return false;
}

}  // namespace RDKit

#endif
