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

#include <string>

#include "Descriptor.hpp"
#include "Mancude.hpp"

namespace RDKit {
namespace NewCIPLabelling {

/**
 * Defines how we can access the properties and connections
 * of a molecule.
 *
 * @param <A> atom type
 * @param <B> bond type
 */
template <typename A, typename B> class BaseMol {
protected:
  std::vector<Mancude::Fraction> atomnums;

public:
  BaseMol() = default;

  virtual ~BaseMol() = default;

  Mancude::Fraction getFractionalAtomicNum(A atom) const {
    if (atomnums.empty())
      const_cast<BaseMol<A, B> *>(this)->atomnums =
          std::move(Mancude::calcFracAtomNums(this));
    return atomnums[getAtomIdx(atom)];
  }

  virtual int getNumAtoms() const = 0;

  virtual int getNumBonds() const = 0;

  virtual A getAtom(int idx) const = 0;

  virtual int getAtomIdx(A atom) const = 0;

  virtual std::vector<A> atoms() const = 0;

  virtual B getBond(int idx) const = 0;

  virtual int getBondIdx(B bond) const = 0;

  virtual std::vector<B> getBonds(A atom) const = 0;

  virtual A getOther(B bond, A atom) const = 0;

  virtual A getBeg(B bond) const = 0;

  virtual A getEnd(B bond) const = 0;

  virtual bool isInRing(B bond) const = 0;

  virtual int getAtomicNum(A atom) const = 0;

  virtual int getNumHydrogens(A atom) const = 0;

#if 0
***** NOT IMPLEMENTED YET *****
  virtual Descriptor getAtomDescriptor(A atom,
                                       const std::string& key) const = 0;

  virtual Descriptor getBondDescriptor(B bond,
                                       const std::string& key) const = 0;
#endif

  virtual int getMassNum(A atom) const = 0;

  virtual int getCharge(A atom) const = 0;

  virtual int getBondOrder(B bond) const = 0;

  virtual void setAtomDescriptor(A atom, const std::string &key,
                                 Descriptor desc) = 0;

  virtual void setBondDescriptor(B bond, const std::string &key,
                                 Descriptor desc) = 0;
};

} // namespace NewCIPLabelling
} // namespace RDKit