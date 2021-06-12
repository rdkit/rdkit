//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <stdexcept>
#include <string>

#include <RDGeneral/export.h>

namespace RDKit {
namespace CIPLabeler {

class RDKIT_CIPLABELER_EXPORT TooManyNodesException
    : public std::runtime_error {
 public:
  TooManyNodesException(const std::string &msg) : std::runtime_error(msg) {}
};

}  // namespace CIPLabeler
}  // namespace RDKit
