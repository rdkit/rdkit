//
//  Copyright (C) 2005-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include "MolDescriptors.h"
#include <map>
#include <list>
#include <algorithm>
#include <sstream>

namespace RDKit {
namespace Descriptors {

const std::string amwVersion = "1.0.0";
double calcAMW(const ROMol &mol, bool onlyHeavy) {
  return MolOps::getAvgMolWt(mol, onlyHeavy);
}

const std::string NumHeavyAtomsVersion = "1.0.0";
unsigned int calcNumHeavyAtoms(const ROMol &mol) {
  return mol.getNumHeavyAtoms();
}

const std::string NumAtomsVersion = "1.0.0";
unsigned int calcNumAtoms(const ROMol &mol) {
  bool onlyExplicit = false;
  return mol.getNumAtoms(onlyExplicit);
}

const std::string exactmwVersion = "1.1.0";
double calcExactMW(const ROMol &mol, bool onlyHeavy) {
  return MolOps::getExactMolWt(mol, onlyHeavy);
}

static std::string _molFormulaVersion = "1.3.0";
std::string calcMolFormula(const ROMol &mol, bool separateIsotopes,
                           bool abbreviateHIsotopes) {
  return MolOps::getMolFormula(mol, separateIsotopes, abbreviateHIsotopes);
}

}  // end of namespace Descriptors
}  // end of namespace RDKit
