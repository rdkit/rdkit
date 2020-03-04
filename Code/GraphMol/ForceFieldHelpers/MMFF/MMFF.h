//
//  Copyright (C) 2015-2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MMFFCONVENIENCE_H
#define RD_MMFFCONVENIENCE_H
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <RDGeneral/RDThreads.h>
#include "AtomTyper.h"
#include "Builder.h"

namespace RDKit {
class ROMol;
namespace MMFF {
//! Convenience function for optimizing a molecule using MMFF
/*
  \param mol        the molecule to use
  \param maxIters   the maximum number of force-field iterations
  \param mmffVariant the MMFF variant to use, should be "MMFF94" or "MMFF94S"
  \param nonBondedThresh  the threshold to be used in adding non-bonded terms
                          to the force field. Any non-bonded contact whose
  current
                    distance is greater than \c nonBondedThresh * the minimum
  value
                    for that contact will not be included.
  \param confId     the optional conformer id, if this isn't provided, the
  molecule's
                    default confId will be used.
  \param ignoreInterfragInteractions if true, nonbonded terms will not be added
  between
                                     fragments

  \return a pair with:
     first: -1 if parameters were missing, 0 if the optimization converged, 1 if
  more iterations are required.
     second: the energy
*/
std::pair<int, double> MMFFOptimizeMolecule(
    ROMol &mol, int maxIters = 1000, std::string mmffVariant = "MMFF94",
    double nonBondedThresh = 10.0, int confId = -1,
    bool ignoreInterfragInteractions = true) {
  std::pair<int, double> res = std::make_pair(-1, -1);
  MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);
  if (mmffMolProperties.isValid()) {
    ForceFields::ForceField *ff = MMFF::constructForceField(
        mol, nonBondedThresh, confId, ignoreInterfragInteractions);
    res = ForceFieldsHelper::OptimizeMolecule(*ff, maxIters);
    delete ff;
  }
  return res;
}

//! Convenience function for optimizing all of a molecule's conformations using
// MMFF
/*
  \param mol        the molecule to use
  \param res        vector of (needsMore,energy) pairs
  \param numThreads the number of simultaneous threads to use (only has an
                    effect if the RDKit is compiled with thread support).
                    If set to zero, the max supported by the system will be
  used.
  \param maxIters   the maximum number of force-field iterations
  \param mmffVariant the MMFF variant to use, should be "MMFF94" or "MMFF94S"
  \param nonBondedThresh  the threshold to be used in adding non-bonded terms
                          to the force field. Any non-bonded contact whose
  current
                    distance is greater than \c nonBondedThresh * the minimum
  value
                    for that contact will not be included.
  \param ignoreInterfragInteractions if true, nonbonded terms will not be added
  between
                                     fragments

*/
void MMFFOptimizeMoleculeConfs(ROMol &mol,
                               std::vector<std::pair<int, double>> &res,
                               int numThreads = 1, int maxIters = 1000,
                               std::string mmffVariant = "MMFF94",
                               double nonBondedThresh = 10.0,
                               bool ignoreInterfragInteractions = true) {
  MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);
  if (mmffMolProperties.isValid()) {
    ForceFields::ForceField *ff = MMFF::constructForceField(
        mol, &mmffMolProperties, nonBondedThresh, -1, ignoreInterfragInteractions);
    ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    delete ff;
  } else {
    res.resize(mol.getNumConformers());
    for (unsigned int i = 0; i < mol.getNumConformers(); ++i) {
      res[i] = std::make_pair(static_cast<int>(-1), static_cast<double>(-1));
    }
  }
}
}  // namespace MMFF
}  // end of namespace RDKit
#endif
