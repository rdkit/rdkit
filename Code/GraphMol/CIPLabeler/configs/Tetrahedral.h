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
  static const int LEFT = 0x1;
  static const int RIGHT = 0x2;

  Tetrahedral();

  Tetrahedral(Atom *focus, std::vector<Atom *> &&carriers, int cfg);

  void setPrimaryLabel(CIPMol *mol, Descriptor desc) override;

  Descriptor label(const SequenceRule *comp) override;

  Descriptor label(Node *node, Digraph *digraph,
                   const SequenceRule *comp) override;

private:
  Descriptor label(Node *node, const SequenceRule *comp) const;
};

} // namespace CIPLabeler
} // namespace RDKit