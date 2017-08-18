//
//  Copyright (C) 2017 Greg Landrum
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

#ifndef RD_MOLBUNDLE_AUG2017
#define RD_MOLBUNDLE_AUG2017

/// Std stuff
#include <vector>

// boost stuff
#include <RDGeneral/BoostStartInclude.h>
#include <boost/smart_ptr.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
class ROMol;

//! MolBundle contains (conceptually) a collection of ROMols with the same
//! topology
/*!
  This is designed for allowing handling of things like enhanced stereochemistr.


 */

class MolBundle : public RDProps {
 public:
  MolBundle() : RDProps(){};

  //! copy constructor with a twist
  /*!
    \param other     the molecule to be copied
  */
  MolBundle(const MolBundle &other) : RDProps(other) { d_mols = other.d_mols; };
  //! construct from a pickle string
  // MolBundle(const std::string &binStr);

  virtual ~MolBundle(){};

  virtual const std::vector<boost::shared_ptr<ROMol> > &getMols() const {
    return d_mols;
  };

  unsigned int addMol(boost::shared_ptr<ROMol> nmol) {
    if (d_mols.size()) {
      // FIX: verify size
    }
    d_mols.push_back(nmol);
    return (d_mols.size());
  }
  unsigned int size() const { return d_mols.size(); };
  const boost::shared_ptr<ROMol> getMol(size_t idx) const {
    PRECONDITION(idx < d_mols.size(), "bad index");
    return d_mols[idx];
  };
  const boost::shared_ptr<ROMol> operator[](size_t idx) const {
    return getMol(idx);
  };

 private:
  std::vector<boost::shared_ptr<ROMol> > d_mols;
};

};  // end of RDKit namespace
#endif
