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

class Sp2Bond : public Configuration {
private:
  Bond *bond;

public:
  static const int OPPOSITE = 0x1;
  static const int TOGETHER = 0x2;

  Sp2Bond();

  Sp2Bond(Bond *bond, std::vector<Atom *> &&foci,
          std::vector<Atom *> &&carriers, int cfg);

  void setPrimaryLabel(CIPMol *mol, Descriptor desc) override;

  Descriptor label(const SequenceRule *comp) override;

  Descriptor label(Node *root1, Digraph *digraph,
                   const SequenceRule *rules) override;
}; // namespace CIPLabeler

} // namespace CIPLabeler
} // namespace RDKit