//
//  Copyright (C) 2017-2023 Greg Landrum and other RDKit contributors
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

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

// our stuff
#include <RDGeneral/Exceptions.h>
#include <GraphMol/MolPickler.h>

namespace RDKit {
class ROMol;

inline bool MolBundleCanSerialize() {
#ifdef RDK_USE_BOOST_SERIALIZATION
  return true;
#else
  return false;
#endif
};

//! MolBundle contains a collection of related ROMols
/*!
  This is designed to allow handling things like enumerating link nodes,
  polymers, etc.
*/
class MolBundle : public RDProps {
 public:
  MolBundle() : RDProps() {}

  //! copy constructor
  MolBundle(const MolBundle &other) : RDProps(other) { d_mols = other.d_mols; }
  MolBundle(const std::string &pkl) { initFromString(pkl); }
  virtual ~MolBundle() {}

  MolBundle &operator=(const MolBundle &other) = default;

  //! returns our molecules
  virtual const std::vector<boost::shared_ptr<ROMol>> &getMols() const {
    return d_mols;
  }

  //! adds a new molecule and returns the total number of molecules
  virtual size_t addMol(boost::shared_ptr<ROMol> nmol) {
    PRECONDITION(nmol.get(), "bad mol pointer");
    d_mols.push_back(nmol);
    return (d_mols.size());
  }
  //! returns the number of molecules from the bundle
  virtual size_t size() const { return d_mols.size(); }
  //! returns whether or not the bundle is empty
  virtual bool empty() const { return d_mols.empty(); }

  //! returns a particular molecule in the bundle
  virtual const boost::shared_ptr<ROMol> getMol(size_t idx) const {
    if (idx >= d_mols.size()) {
      throw IndexErrorException(static_cast<int>(idx));
    }
    return d_mols[idx];
  }
  //! returns a particular molecule from the bundle
  virtual const boost::shared_ptr<ROMol> operator[](size_t idx) const {
    return getMol(idx);
  }

  //! serializes (pickles) to a stream
  void toStream(std::ostream &ss) const {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    boost::archive::text_oarchive ar(ss);
    ar << *this;
#endif
  };
  //! returns a string with a serialized (pickled) representation
  std::string serialize() const {
    std::stringstream ss;
    toStream(ss);
    return ss.str();
  };
  //! initializes from a stream pickle
  void initFromStream(std::istream &ss) {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    boost::archive::text_iarchive ar(ss);
    ar >> *this;
#endif
  };
  //! initializes from a string pickle
  void initFromString(const std::string &text) {
    std::stringstream ss(text);
    initFromStream(ss);
  };

#ifdef RDK_USE_BOOST_SERIALIZATION
  // FIX: we don't currently serialize properties
  template <class Archive>
  void save(Archive &ar, const unsigned int version) const {
    RDUNUSED_PARAM(version);
    std::vector<std::string> pkls;
    for (const auto &mol : d_mols) {
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      pkls.push_back(pkl);
    }
    ar << pkls;
  }

  template <class Archive>
  void load(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);

    std::vector<std::string> pkls;
    ar >> pkls;
    d_mols.clear();
    for (const auto &pkl : pkls) {
      d_mols.push_back(ROMOL_SPTR(new ROMol(pkl)));
    }
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

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
  FixedMolSizeMolBundle() : MolBundle() {}

  //! copy constructor
  FixedMolSizeMolBundle(const FixedMolSizeMolBundle &other)
      : MolBundle(other) {}

  ~FixedMolSizeMolBundle() override = default;

  //! adds a new molecule and returns the total number of molecules
  //!  enforces that the new molecule has the same number of atoms and bonds
  //!  as the molecules that are already there.
  size_t addMol(boost::shared_ptr<ROMol> nmol) override {
    PRECONDITION(nmol.get(), "bad mol pointer");
    if (d_mols.size()) {
      if (nmol->getNumAtoms() != d_mols[0]->getNumAtoms()) {
        throw ValueErrorException(
            "all molecules in a bundle must have the same number of atoms");
      }
      // REVIEW: should we allow different numbers of bonds?
      if (nmol->getNumBonds() != d_mols[0]->getNumBonds()) {
        throw ValueErrorException(
            "all molecules in a bundle must have the same number of bonds");
      }
    }
    d_mols.push_back(nmol);
    return (d_mols.size());
  }
};

};  // namespace RDKit
#endif
