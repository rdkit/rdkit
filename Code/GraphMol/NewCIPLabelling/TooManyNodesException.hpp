//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
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

namespace RDKit
{
namespace NewCIPLabelling
{

class TooManyNodesException : public std::runtime_error
{
  public:
    TooManyNodesException(const std::string& msg) : std::runtime_error(msg){};
};

} // namespace NewCIPLabelling
} // namespace RDKit