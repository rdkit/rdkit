//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MacroAtomInfo.h"

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <array>
#include <utility>

namespace RDKit {
namespace {
const std::array<std::pair<MonomerClass, const char *>, 4> monomerClassNames = {
    {
        {MonomerClass::AminoAcid, "AminoAcid"},
        {MonomerClass::NucleicAcid, "NucleicAcid"},
        {MonomerClass::Chemical, "Chemical"},
        {MonomerClass::Other, "Other"},
    }};
}  // namespace

const char *monomerClassToString(MonomerClass monomerClass) {
  for (const auto &[value, name] : monomerClassNames) {
    if (value == monomerClass) {
      return name;
    }
  }
  POSTCONDITION(false, "unknown monomer class");
  return "";
}

MonomerClass monomerClassFromString(const std::string &monomerClass) {
  for (const auto &[value, name] : monomerClassNames) {
    if (monomerClass == name) {
      return value;
    }
  }
  BOOST_LOG(rdWarningLog) << "unrecognized monomer class '" << monomerClass
                          << "'; treating it as Other" << std::endl;
  return MonomerClass::Other;
}

}  // namespace RDKit
