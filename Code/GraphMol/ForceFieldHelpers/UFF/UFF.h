//
//  Copyright (C) 2015-2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_UFFCONVENIENCE_H
#define RD_UFFCONVENIENCE_H
#include <ForceField/ForceField.h>
#include <RDGeneral/RDThreads.h>
#include "Builder.h"

namespace RDKit {
class ROMol;
namespace UFF {
//! Convenience function for optimizing a molecule using UFF
/*
  \param mol        the molecule to use
  \param maxIters   the maximum number of force-field iterations
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

  \return a pair with:
     first: 0 if the optimization converged, 1 if more iterations are required.
     second: the energy
*/
std::pair<int, double> UFFOptimizeMolecule(
    ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
    bool ignoreInterfragInteractions = true) {
  ForceFields::ForceField *ff = UFF::constructForceField(
      mol, vdwThresh, confId, ignoreInterfragInteractions);
  ff->initialize();
  int res = ff->minimize(maxIters);
  double e = ff->calcEnergy();
  delete ff;
  return std::make_pair(res, e);
}
#ifdef RDK_THREADSAFE_SSS
namespace detail {
void UFFOptimizeMoleculeConfsHelper_(ForceFields::ForceField ff, ROMol *mol,
                                     std::vector<std::pair<int, double>> *res,
                                     unsigned int threadIdx,
                                     unsigned int numThreads, int maxIters) {
  unsigned int i = 0;
  ff.positions().resize(mol->getNumAtoms());
  for (ROMol::ConformerIterator cit = mol->beginConformers();
       cit != mol->endConformers(); ++cit, ++i) {
    if (i % numThreads != threadIdx) continue;
    for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
      ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    }
    ff.initialize();
    int needsMore = ff.minimize(maxIters);
    double e = ff.calcEnergy();
    (*res)[i] = std::make_pair(needsMore, e);
  }
}
}  // end of detail namespace
#endif
//! Convenience function for optimizing all of a molecule's conformations using
// UFF
/*
  \param mol        the molecule to use
  \param res        vector of (needsMore,energy)
  \param numThreads the number of simultaneous threads to use (only has an
                    effect if the RDKit is compiled with thread support).
                    If set to zero, the max supported by the system will be
  used.
  \param maxIters   the maximum number of force-field iterations
  \param vdwThresh  the threshold to be used in adding van der Waals terms
                    to the force field. Any non-bonded contact whose current
                    distance is greater than \c vdwThresh * the minimum value
                    for that contact will not be included.
  \param ignoreInterfragInteractions if true, nonbonded terms will not be added
  between
                                     fragments

*/
void UFFOptimizeMoleculeConfs(ROMol &mol,
                              std::vector<std::pair<int, double>> &res,
                              int numThreads = 1, int maxIters = 1000,
                              double vdwThresh = 10.0,
                              bool ignoreInterfragInteractions = true) {
  res.resize(mol.getNumConformers());
  numThreads = getNumThreadsToUse(numThreads);
  if (numThreads == 1) {
    unsigned int i = 0;
    for (ROMol::ConformerIterator cit = mol.beginConformers();
         cit != mol.endConformers(); ++cit, ++i) {
      res[i] = UFFOptimizeMolecule(mol, maxIters, vdwThresh, (*cit)->getId(),
                                   ignoreInterfragInteractions);
    }
  }
#ifdef RDK_THREADSAFE_SSS
  else {
    ForceFields::ForceField *ff = UFF::constructForceField(
        mol, vdwThresh, -1, ignoreInterfragInteractions);
    std::vector<std::thread> tg;
    for (int ti = 0; ti < numThreads; ++ti) {
      tg.emplace_back(std::thread(detail::UFFOptimizeMoleculeConfsHelper_, *ff,
                                  &mol, &res, ti, numThreads, maxIters));
    }
    for (auto &thread : tg) {
      if (thread.joinable()) thread.join();
    }
    delete ff;
  }
#endif
}
}  // end of namespace UFF
}  // end of namespace RDKit
#endif
