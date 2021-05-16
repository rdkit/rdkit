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
 public:
  Sp2Bond() = delete;

  Sp2Bond(const CIPMol &mol, Bond *bond, Atom *startAtom, Atom *endAtom,
          Bond::BondStereo cfg);

  void setPrimaryLabel(Descriptor desc) override;

  Descriptor label(const Rules &comp) override;

  Descriptor label(Node *root1, Digraph &digraph, const Rules &comp) override;

 private:
  Bond *dp_bond;

  // bond->getStereo() can return both E/Z or CIS/TRANS,
  // so we cache CIS/TRANS we found.
  Bond::BondStereo d_cfg;

};  // namespace CIPLabeler

}  // namespace CIPLabeler
}  // namespace RDKit
