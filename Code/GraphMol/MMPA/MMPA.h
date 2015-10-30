//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"

namespace RDKit {

  namespace MMPA {
  
    bool fragmentMol(const ROMol &mol,
                     std::vector< std::pair<ROMOL_SPTR,ROMOL_SPTR> >& result,
                     unsigned int maxCuts=3,
                     const std::string& pattern="[#6+0;!$(*=,#[!#6])]!@!=!#[*]");
  }
} // namespace RDKit
