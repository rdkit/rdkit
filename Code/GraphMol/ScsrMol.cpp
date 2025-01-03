//
//  Copyright (C) 2003-2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/tokenizer.hpp>

// our stuff
#include "RWMol.h"
#include "ScsrMol.h"
#include <RDGeneral/FileParseException.h>

namespace RDKit {
SCSRMol::SCSRMol() { mol = std::unique_ptr<ROMol>(new RWMol()); }

SCSRMol::SCSRMol(const SCSRMol &other) {
  mol = std::unique_ptr<ROMol>(other.mol.get());
  for (auto &oneTemplate : other.templates) {
    templates.push_back(std::make_unique<RWMol>(*oneTemplate));
  }
}
SCSRMol &SCSRMol::operator=(const SCSRMol &) {
  throw RDKit::FileParseException("Assignment not allowed for SCSRMol");
};

SCSRMol::SCSRMol(SCSRMol &&other) noexcept {
  mol = std::move(other.mol);
  for (auto &oneTemplate : other.templates) {
    templates.push_back(std::move(oneTemplate));
  }
}
SCSRMol &SCSRMol::operator=(SCSRMol &&other) noexcept {
  mol = std::move(other.mol);
  for (auto &oneTemplate : other.templates) {
    templates.push_back(std::move(oneTemplate));
  }
  return *this;
}

}  // namespace RDKit
