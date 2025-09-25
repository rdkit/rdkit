//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDMOL_THROW_H
#define RDMOL_THROW_H

#include <string>
#include <stdexcept>

[[noreturn]] inline void raiseNonImplementedFunction(const std::string& func) {
  std::string errorMessage = "This function is not implemented: " + func;
  throw std::runtime_error(errorMessage);
}

[[noreturn]] inline void raiseNonImplementedDetail(const std::string& detail) {
  std::string errorMessage = "This detail is not implemented: " + detail;
  throw std::runtime_error(errorMessage);
}

#endif // RDMOL_THROW_H
