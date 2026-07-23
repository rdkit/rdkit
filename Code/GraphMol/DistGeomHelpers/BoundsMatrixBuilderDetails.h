//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDGeneral/Invariant.h"
#include <vector>
#include <optional>

#ifndef RD_BOUNDS_MATRIX_BUILDER_DETAILS_H
#define RD_BOUNDS_MATRIX_BUILDER_DETAILS_H

namespace RDKit {
namespace DGeomHelpers {
enum class TorsionType {
  CIS = 0,
  TRANS,
  FLEXIBLE,
  CUSTOM,
  NONE  // don't set the bound
};

struct TorsionValue {
  TorsionType type = TorsionType::NONE;
  std::optional<double> value = {};
  std::optional<double> extraDist = {};
  bool isForced = false;
};

//! A structure used to store 14 paths - cis/trans info
struct Path14Configuration {
  unsigned int bid1, bid2, bid3;
  unsigned int aid1, aid2, aid3, aid4;
  TorsionValue type;
};
using PATH14_VECT = std::vector<Path14Configuration>;
}  // namespace DGeomHelpers
}  // namespace RDKit
#endif