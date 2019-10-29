//
// Copyright (C) 2018 by Greg Landrum
//
#ifndef EHTTOOLS_H_20181226
#define EHTTOOLS_H_20181226
/*! \file

  \brief Contains an interface to the YaEHMOP extended Hueckel program.

  \b Note: This interface is experimental and may change from version to
  version.

*/

#include <string>
#include <memory>

namespace RDKit {
class ROMol;
namespace EHTTools {

struct EHTResults {
  unsigned int numAtoms;
  unsigned int numOrbitals;
  std::unique_ptr<double[]> overlapPopulationMatrix;
  std::unique_ptr<double[]> reducedOverlapPopulationMatrix;
  std::unique_ptr<double[]> chargeMatrix;
  std::unique_ptr<double[]> reducedChargeMatrix;
  std::unique_ptr<double[]> atomicCharges;
  double fermiEnergy;
  double totalEnergy;
  EHTResults() = default;
  EHTResults(const EHTResults &) = delete;
  EHTResults &operator=(const EHTResults &) = delete;
};

//! Runs an extended Hueckel calculation for a molecule
//!   The results are returned in the EHTResults structure
bool runMol(const ROMol &mol, EHTResults &results, int confId = -1);

}  // namespace EHTTools
}  // namespace RDKit
#endif
