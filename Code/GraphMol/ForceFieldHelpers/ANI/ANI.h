#include <RDGeneral/export.h>
#ifndef RD_ANICONVENIENCE_H_
#define RD_ANICONVENIENCE_H_

#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include "Builder.h"

namespace RDKit {
class ROMol;
namespace ANI {
//! Convenience function for optimizing a molecule using ANI
/*
  \param mol          the molecule to use
  \param maxIters     the maximum number of force-field iterations
  \param model        Param dir, such as one in /Code/ForceField/ANI/Params/
  \param confId       the optional conformer id, if this isn't provided, the
  molecule's default confId will be used.

  \return a pair with:
     first: 0 if the optimization converged, 1 if more iterations are required.
     second: the energy
*/
std::pair<int, double> ANIOptimizeMolecule(ROMol &mol, std::string model,
                                           int confId = -1,
                                           int maxIters = 1000) {
  std::unique_ptr<ForceFields::ForceField> ff(
      ANI::constructForceField(mol, model, confId));
  std::pair<int, double> res =
      ForceFieldsHelper::OptimizeMolecule(*ff, maxIters);
  return res;
}

//! Convenience function for optimizing all of a molecule's conformations using
// ANI
/*
  \param mol          the molecule to use
  \param res          vector of (needsMore,energy)
  \param numThreads   the number of simultaneous threads to use (only has an
                      effect if the RDKit is compiled with thread support).
                      If set to zero, the max supported by the system will be
                      used.
  \param model        Param dir, such as one in /Code/ForceField/ANI/Params/
*/
void ANIOptimizeMoleculeConfs(ROMol &mol,
                              std::vector<std::pair<int, double>> &res,
                              std::string model,
                              int numThreads = 1, int maxIters = 1000) {
  std::unique_ptr<ForceFields::ForceField> ff(
      ANI::constructForceField(mol, model, -1));
  ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
}
}  // namespace ANI
}  // namespace RDKit

#endif