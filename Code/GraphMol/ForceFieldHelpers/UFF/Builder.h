//
//  Copyright (C) 2004-2018 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_UFFBUILDER_H
#define RD_UFFBUILDER_H

#include <vector>
#include <string>
#include <boost/shared_array.hpp>
#include <boost/scoped_ptr.hpp>
#ifdef RDK_THREADSAFE_SSS
#include <mutex>
#endif
#include <boost/noncopyable.hpp>

namespace ForceFields {
class ForceField;
namespace UFF {
class AtomicParams;
}
}  // namespace ForceFields

namespace RDKit {
class ROMol;
namespace UFF {
typedef std::vector<const ForceFields::UFF::AtomicParams *> AtomicParamVect;

//! Builds and returns a UFF force field for a molecule
/*!

  \param mol        the molecule to use
  \param vdwThresh  the threshold to be used in adding van der Waals terms
                    to the force field. Any non-bonded contact whose current
                    distance is greater than \c vdwThresh * the minimum value
                    for that contact will not be included.
  \param confId     the optional conformer id, if this isn't provided, the
  molecule's
                    default confId will be used.
  \param ignoreInterfragInteractions if true, nonbonded terms will not be added
  between
                                     fragments

  \return the new force field. The client is responsible for free'ing this.
*/
RDKIT_FORCEFIELDHELPERS_EXPORT ForceFields::ForceField *constructForceField(
    ROMol &mol, double vdwThresh = 100.0, int confId = -1,
    bool ignoreInterfragInteractions = true);

//! Builds and returns a UFF force field for a molecule
/*!

  \param mol        the molecule to use
  \param params     a vector with pointers to the ForceFields::UFF::AtomicParams
                    structures to be used
  \param vdwThresh  the threshold to be used in adding van der Waals terms
                    to the force field. Any non-bonded contact whose current
                    distance is greater than \c vdwThresh * the minimum value
                    for that contact will not be included.
  \param confId     the optional conformer id, if this isn't provided, the
  molecule's
                    default confId will be used.
  \param ignoreInterfragInteractions if true, nonbonded terms will not be added
  between
                                     fragments

  \return the new force field. The client is responsible for free'ing this.
*/
RDKIT_FORCEFIELDHELPERS_EXPORT ForceFields::ForceField *constructForceField(
    ROMol &mol, const AtomicParamVect &params, double vdwThresh = 100.0,
    int confId = -1, bool ignoreInterfragInteractions = true);

namespace Tools {
class RDKIT_FORCEFIELDHELPERS_EXPORT DefaultTorsionBondSmarts
    : private boost::noncopyable {
 public:
  static const std::string &string() { return ds_string; }
  static const ROMol *query();

 private:
  DefaultTorsionBondSmarts() {}
  static void create();
  static const std::string ds_string;
  static boost::scoped_ptr<const ROMol> ds_instance;
#ifdef RDK_THREADSAFE_SSS
  static std::once_flag ds_flag;
#endif
};

enum { RELATION_1_2 = 0, RELATION_1_3 = 1, RELATION_1_4 = 2, RELATION_1_X = 3 };
// these functions are primarily exposed so they can be tested.
RDKIT_FORCEFIELDHELPERS_EXPORT unsigned int twoBitCellPos(unsigned int nAtoms,
                                                          int i, int j);
RDKIT_FORCEFIELDHELPERS_EXPORT void setTwoBitCell(
    boost::shared_array<std::uint8_t> &res, unsigned int pos,
    std::uint8_t value);
RDKIT_FORCEFIELDHELPERS_EXPORT std::uint8_t getTwoBitCell(
    boost::shared_array<std::uint8_t> &res, unsigned int pos);
RDKIT_FORCEFIELDHELPERS_EXPORT boost::shared_array<std::uint8_t>
buildNeighborMatrix(const ROMol &mol);
RDKIT_FORCEFIELDHELPERS_EXPORT void addBonds(const ROMol &mol,
                                             const AtomicParamVect &params,
                                             ForceFields::ForceField *field);
RDKIT_FORCEFIELDHELPERS_EXPORT void addAngles(const ROMol &mol,
                                              const AtomicParamVect &params,
                                              ForceFields::ForceField *field);
RDKIT_FORCEFIELDHELPERS_EXPORT void addNonbonded(
    const ROMol &mol, int confId, const AtomicParamVect &params,
    ForceFields::ForceField *field,
    boost::shared_array<std::uint8_t> neighborMatrix, double vdwThresh = 100.0,
    bool ignoreInterfragInteractions = true);
RDKIT_FORCEFIELDHELPERS_EXPORT void addTorsions(
    const ROMol &mol, const AtomicParamVect &params,
    ForceFields::ForceField *field,
    const std::string &torsionBondSmarts = DefaultTorsionBondSmarts::string());
RDKIT_FORCEFIELDHELPERS_EXPORT void addInversions(
    const ROMol &mol, const AtomicParamVect &params,
    ForceFields::ForceField *field);
}  // namespace Tools
}  // namespace UFF
}  // namespace RDKit

#endif
