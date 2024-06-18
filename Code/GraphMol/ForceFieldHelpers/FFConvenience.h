//
//  Copyright (C) 2019 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FFCONVENIENCE_H
#define RD_FFCONVENIENCE_H
#include <ForceField/ForceField.h>
#include <RDGeneral/RDThreads.h>

namespace RDKit {
class ROMol;
namespace ForceFieldsHelper {
namespace detail {
#ifdef RDK_BUILD_THREADSAFE_SSS
inline void OptimizeMoleculeConfsHelper_(
    ForceFields::ForceField ff, ROMol *mol,
    std::vector<std::pair<int, double>> *res, unsigned int threadIdx,
    unsigned int numThreads, int maxIters) {
  PRECONDITION(mol, "mol must not be nullptr");
  PRECONDITION(res, "res must not be nullptr");
  PRECONDITION(res->size() >= mol->getNumConformers(),
               "res->size() must be >= mol->getNumConformers()");
  unsigned int i = 0;
  ff.positions().resize(mol->getNumAtoms());
  for (ROMol::ConformerIterator cit = mol->beginConformers();
       cit != mol->endConformers(); ++cit, ++i) {
    if (i % numThreads != threadIdx) {
      continue;
    }
    for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
      ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    }
    ff.initialize();
    int needsMore = ff.minimize(maxIters);
    double e = ff.calcEnergy();
    (*res)[i] = std::make_pair(needsMore, e);
  }
}

inline void OptimizeMoleculeConfsMT(ROMol &mol,
                                    const ForceFields::ForceField &ff,
                                    std::vector<std::pair<int, double>> &res,
                                    int numThreads, int maxIters) {
  std::vector<std::thread> tg;
  for (int ti = 0; ti < numThreads; ++ti) {
    tg.emplace_back(std::thread(detail::OptimizeMoleculeConfsHelper_, ff, &mol,
                                &res, ti, numThreads, maxIters));
  }
  for (auto &thread : tg) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}
#endif

inline void OptimizeMoleculeConfsST(ROMol &mol, ForceFields::ForceField &ff,
                                    std::vector<std::pair<int, double>> &res,
                                    int maxIters) {
  PRECONDITION(res.size() >= mol.getNumConformers(),
               "res.size() must be >= mol.getNumConformers()");
  unsigned int i = 0;
  for (ROMol::ConformerIterator cit = mol.beginConformers();
       cit != mol.endConformers(); ++cit, ++i) {
    for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
      ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    }
    ff.initialize();
    int needsMore = ff.minimize(maxIters);
    double e = ff.calcEnergy();
    res[i] = std::make_pair(needsMore, e);
  }
}
}  // namespace detail

//! Convenience function for optimizing a molecule using a pre-generated
//! force-field
/*
  \param ff         the force-field
  \param res        vector of (needsMore,energy) pairs
  \param maxIters   the maximum number of force-field iterations

  \return a pair with:
     first: -1 if parameters were missing, 0 if the optimization converged, 1 if
  more iterations are required.
     second: the energy
*/
inline std::pair<int, double> OptimizeMolecule(ForceFields::ForceField &ff,
                                               int maxIters = 1000) {
  ff.initialize();
  int res = ff.minimize(maxIters);
  double e = ff.calcEnergy();
  return std::make_pair(res, e);
}

//! Convenience function for optimizing all of a molecule's conformations using
/// a pre-generated force-field
/*
  \param mol        the molecule to use
  \param ff         the force-field
  \param res        vector of (needsMore,energy) pairs
  \param numThreads the number of simultaneous threads to use (only has an
                    effect if the RDKit is compiled with thread support).
                    If set to zero, the max supported by the system will be
  used.
  \param maxIters   the maximum number of force-field iterations

*/
inline void OptimizeMoleculeConfs(ROMol &mol, ForceFields::ForceField &ff,
                                  std::vector<std::pair<int, double>> &res,
                                  int numThreads = 1, int maxIters = 1000) {
  res.resize(mol.getNumConformers());
  numThreads = getNumThreadsToUse(numThreads);
  if (numThreads == 1) {
    detail::OptimizeMoleculeConfsST(mol, ff, res, maxIters);
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    detail::OptimizeMoleculeConfsMT(mol, ff, res, numThreads, maxIters);
  }
#endif
}

//! Convenience Function for generating an empty force Field with just
/// the molecules' atoms position.
/*
  \param mol  The molecule which positions should be added to the force field.
  \param confId Id of the conformer which positions are used in the force field.

*/
inline std::unique_ptr<ForceFields::ForceField> createEmptyForceFieldForMol(
    ROMol &mol, int confId = -1) {
  auto res = std::make_unique<ForceFields::ForceField>();
  auto &conf = mol.getConformer(confId);
  for (auto &pt : conf.getPositions()) {
    res->positions().push_back(&pt);
  }
  return res;
}
}  // end of namespace ForceFieldsHelper
}  // end of namespace RDKit
#endif
