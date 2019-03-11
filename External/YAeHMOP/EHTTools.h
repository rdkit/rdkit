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

namespace RDKit {
class ROMol;
namespace EHTTools {
extern const std::string _EHTCharge;        // used to hold partial charges
extern const std::string _EHTMullikenOP;    // used to hold overlap populations
extern const std::string _EHTChargeMatrix;  // used to hold charge matrix

//! Runs an extended Hueckel calculation for a molecule
//!   The results are stored as molecular properties
bool runMol(const ROMol &mol, int confId = -1);

}  // namespace EHTTools
}  // namespace RDKit
#endif
