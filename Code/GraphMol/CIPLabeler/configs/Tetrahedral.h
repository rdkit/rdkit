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

#include "Configuration.h"

namespace RDKit {
namespace CIPLabeler {

class Tetrahedral : public Configuration {
 public:
  Tetrahedral() = delete;

  Tetrahedral(const CIPMol &mol, Atom *focus);

  void setPrimaryLabel(Descriptor desc) override;

  Descriptor label(const Rules &comp) override;

  Descriptor label(Node *node, Digraph &digraph, const Rules &comp) override;

 private:
  Descriptor label(Node *node, const Rules &comp) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
