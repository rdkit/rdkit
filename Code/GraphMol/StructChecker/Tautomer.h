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

namespace RDKit {
namespace StructureCheck {

class RDKIT_STRUCTCHECKER_EXPORT StructCheckTautomer {
  RWMol &Mol;
  const StructCheckerOptions &Options;

 public:
  StructCheckTautomer(RWMol &mol, const StructCheckerOptions &options)
      : Mol(mol), Options(options) {}
  bool applyTautomer(unsigned it);
};
};  // namespace StructureCheck
}  // namespace RDKit
