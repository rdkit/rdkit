//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include "StructChecker.h"
#include "Pattern.h"
#include "Utilites.h"

namespace RDKit {
namespace StructureCheck {
/*
 * Returns the total charge of all atoms in molecule.
 */
int TotalCharge(const ROMol &mol);

class RDKIT_STRUCTCHECKER_EXPORT ChargeFix {
  const StructCheckerOptions &Options;
  RWMol &Mol;
  std::vector<unsigned> BondColor;
  std::vector<unsigned> AtomColor;
  std::vector<double> AtompKaValue;
  std::vector<double> AtomOldpKaValue;

 public:
  ChargeFix(const StructCheckerOptions &op, RWMol &mol)
      : Options(op), Mol(mol) {
    resetColors();
    resetValues();
  }
  /*
   * Removes hydrogens from *mp until desired_charge is reached. The
   * positions for hydrogen removal are selected by "acidity" combined
   * with a refinement algorithm. It returns TRUE if molecule could be
   * neutralized and FALSE if any problem were encountered.
   * *ndeprot and *nrefine are set to the number of deprotonations
   * and refinement cycles performed.
   */
  bool rechargeMolecule(unsigned &ndeprot, unsigned &nrefine);

 private:  // internal helpers:
  bool setpKaValues();
  void decrementMarkedCharges();
  int markMostAcidicAtoms(double &pKa_value, double &gap);
  int refineAcidicAtoms(std::vector<unsigned> &numbering);
  void resetColors();
  void resetValues();
};
}  // namespace StructureCheck
}  // namespace RDKit
