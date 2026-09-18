//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_CTRE_H
#define RDKIT_CTRE_H

#include <ctre.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace RDKit::ctre_utils {

template <ctll::fixed_string Pattern>
std::vector<std::string> split(std::string_view input) {
  std::vector<std::string> result;
  size_t previousEnd = 0;
  for (const auto &match : ctre::search_all<Pattern>(input)) {
    const auto delimiter = match.template get<0>().to_view();
    const auto delimiterStart =
        static_cast<size_t>(delimiter.data() - input.data());
    result.emplace_back(input.substr(previousEnd, delimiterStart - previousEnd));
    previousEnd = delimiterStart + delimiter.size();
  }
  if (previousEnd < input.size()) {
    result.emplace_back(input.substr(previousEnd));
  }
  return result;
}

}  // namespace RDKit::ctre_utils

#endif  // RDKIT_CTRE_H
