//
//  Copyright (C) 2017-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file MolBundle.h

  \brief Defines a class for managing bundles of molecules

*/

#include <RDGeneral/export.h>
#ifndef RD_MOLBUNDLE_AUG2017
#define RD_MOLBUNDLE_AUG2017

/// Std stuff
#include <vector>

// boost stuff
#include <RDGeneral/BoostStartInclude.h>
#include <boost/smart_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>

// our stuff
#include <RDGeneral/Exceptions.h>

namespace RDKit {
class ROMol;

//! MolBundle contains a collection of related ROMols
/*!
  This is designed to allow handling things like enumerating link nodes,
  polymers, etc.
*/
class MolBundle : public RDProps {
 public:
  MolBundle() : RDProps(){};

  //! copy constructor
  MolBundle(const MolBundle &other) : RDProps(other) { d_mols = other.d_mols; };
  // FIX: need serialization/deserialization

  virtual ~MolBundle(){};

  //! returns our molecules
  virtual const std::vector<boost::shared_ptr<ROMol>> &getMols() const {
    return d_mols;
  };

  //! adds a new molecule and returns the total number of molecules
  virtual size_t addMol(boost::shared_ptr<ROMol> nmol) {
    PRECONDITION(nmol.get(), "bad mol pointer");
    d_mols.push_back(nmol);
    return (d_mols.size());
  }
  //! returns the number of molecules from the bundle
  virtual size_t size() const { return d_mols.size(); };
  //! returns a particular molecule in the bundle
  virtual const boost::shared_ptr<ROMol> getMol(size_t idx) const {
    if (idx >= d_mols.size()) throw IndexErrorException(static_cast<int>(idx));
    return d_mols[idx];
  };
  //! returns a particular molecule from the bundle
  virtual const boost::shared_ptr<ROMol> operator[](size_t idx) const {
    return getMol(idx);
  };

 protected:
  std::vector<boost::shared_ptr<ROMol>> d_mols;
};

//! FixedMolSizeMolBundle contains a collection of ROMols with the same
//! number of atoms and bonds.
/*!
  This is designed to allow handling things like enhanced stereochemistry,
  but can no doubt be (ab)used in other ways.

  Implementation note: at the moment this isn't taking advantage of the fact
  that the number of atoms and bonds remains constant. This may be used in the
  future to allow this to be more efficient.

*/
class FixedMolSizeMolBundle : public MolBundle {
 public:
  FixedMolSizeMolBundle() : MolBundle(){};

  //! copy constructor
  FixedMolSizeMolBundle(const FixedMolSizeMolBundle &other)
      : MolBundle(other){};

  ~FixedMolSizeMolBundle() override{};

  //! adds a new molecule and returns the total number of molecules
  //!  enforces that the new molecule has the same number of atoms and bonds
  //!  as the molecules that are already there.
  size_t addMol(boost::shared_ptr<ROMol> nmol) override {
    PRECONDITION(nmol.get(), "bad mol pointer");
    if (d_mols.size()) {
      if (nmol->getNumAtoms() != d_mols[0]->getNumAtoms())
        throw ValueErrorException(
            "all molecules in a bundle must have the same number of atoms");
      // REVIEW: should we allow different numbers of bonds?
      if (nmol->getNumBonds() != d_mols[0]->getNumBonds())
        throw ValueErrorException(
            "all molecules in a bundle must have the same number of bonds");
    }
    d_mols.push_back(nmol);
    return (d_mols.size());
  }
};

};  // namespace RDKit
#endif
