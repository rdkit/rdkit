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

#include "BaseMol.hpp"

namespace RDKit {

class ROMol;
class Bond;
class Atom;

namespace NewCIPLabelling {

using RdkA = Atom*;
using RdkB = Bond*;

class RDKitCipMol : public BaseMol<RdkA, RdkB> {
 private:
  ROMol* mol;

 public:
  explicit RDKitCipMol(ROMol* mol);

  int getNumAtoms() const override;

  int getNumBonds() const override;

  RdkA getAtom(int idx) const override;

  int getAtomIdx(RdkA atom) const override;

  std::vector<RdkA> atoms() const override;

  RdkB getBond(int idx) const override;

  int getBondIdx(RdkB bond) const override;

  std::vector<RdkB> getBonds(RdkA atom) const override;

  RdkA getOther(RdkB bond, RdkA atom) const override;

  RdkA getBeg(RdkB bond) const override;

  RdkA getEnd(RdkB bond) const override;

  bool isInRing(RdkB bond) const override;

  int getAtomicNum(RdkA atom) const override;

  int getNumHydrogens(RdkA atom) const override;

#if 0
***** NOT IMPLEMENTED YET *****
  Descriptor getAtomDescriptor(RdkA atom,
                               const std::string& key) const override;

  Descriptor getBondDescriptor(RdkB bond,
                               const std::string& key) const override;
#endif

  int getMassNum(RdkA atom) const override;

  int getCharge(RdkA atom) const override;

  int getBondOrder(RdkB bond) const override;

  void setAtomDescriptor(RdkA atom, const std::string& key,
                         Descriptor desc) override;

  void setBondDescriptor(RdkB bond, const std::string& key,
                         Descriptor desc) override;
};

}  // namespace NewCIPLabelling
}  // namespace RDKit