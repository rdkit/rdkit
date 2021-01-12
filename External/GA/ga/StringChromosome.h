//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef BINARYSTRINGCHROMOSOME_H_
#define BINARYSTRINGCHROMOSOME_H_

#include "StringChromosomeBase.h"
#include "BinaryStringChromosomePolicy.h"
#include "IntegerStringChromosomePolicy.h"
#include "../util/export.h"

namespace GapeGa {

template <typename T, typename ChromosomePolicy>
class StringChromosome : public StringChromosomeBase<T, ChromosomePolicy> {
 private:
  // we should never use the default implementation
  StringChromosome() = delete;
};

template <>
class GA_EXPORT StringChromosome<bool, BinaryStringChromosomePolicy>
    : public StringChromosomeBase<bool, BinaryStringChromosomePolicy> {
 private:
  StringChromosome(const StringChromosome &other) = delete;
  StringChromosome &operator=(const StringChromosome &other) = delete;

 public:
  StringChromosome(int length_, GarethUtil::RandomUtil &rng_,
                   BinaryStringChromosomePolicy &chromosomePolicy_)
      : StringChromosomeBase(length_, rng_, chromosomePolicy_) {
    ;
  }

  int decodeToInt(int start, int nBits) const;
  int decodeByte(int byteNo) const { return decodeToInt(byteNo * 8, 8); }
};

using IntegerStringChromosome =
    StringChromosome<int, IntegerStringChromosomePolicy>;

template <>
class GA_EXPORT StringChromosome<int, IntegerStringChromosomePolicy>
    : public StringChromosomeBase<int, IntegerStringChromosomePolicy> {
 private:
  StringChromosome(const StringChromosome &other) = delete;
  StringChromosome &operator=(const StringChromosome &other) = delete;

 public:
  StringChromosome(int length_, GarethUtil::RandomUtil &rng_,
                   IntegerStringChromosomePolicy &chromosomePolicy_)
      : StringChromosomeBase(length_, rng_, chromosomePolicy_) {
    ;
  }

  void fullMixing(const IntegerStringChromosome &parent2,
                  IntegerStringChromosome &child1,
                  IntegerStringChromosome &child2) const;

  void fullMixingAndCrossover(const IntegerStringChromosome &parent2,
                              IntegerStringChromosome &child1,
                              IntegerStringChromosome &child2) const;
};

using BinaryStringChromosome =
    StringChromosome<bool, BinaryStringChromosomePolicy>;

}  // namespace GapeGa

#endif /* BINARYSTRINGCHROMOSOME_H_ */
