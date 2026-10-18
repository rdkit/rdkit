//
//  Copyright (C) 2021 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*
 * This file is just to make sure exceptions in this directory
 * get built into the library.
 */

#include "BadFileException.h"
#include "Exceptions.h"
#include "FileParseException.h"

IndexErrorException::~IndexErrorException() noexcept = default;
ValueErrorException::~ValueErrorException() noexcept = default;
KeyErrorException::~KeyErrorException() noexcept = default;

namespace RDKit {
BadFileException::~BadFileException() noexcept = default;
FileParseException::~FileParseException() noexcept = default;
}  // namespace RDKit
