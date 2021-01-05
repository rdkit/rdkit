//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "StringChromosome.h"

namespace GapeGa {

/**
 * A decode to integer function for binary string chromosomes
 *
 * @param start
 * @param nBits
 * @return
 */
int StringChromosome<bool, BinaryStringChromosomePolicy>::decodeToInt(
    int start, int nBits) const {
  assert(length >= (start + nBits));

  int mask = 1, result = 0;
  bool *ptr = string.get() + start;
  for (int i = 0; i < nBits; i++, mask <<= 1) {
    if (*ptr++ == true) result |= mask;
  }
  return result;
}

/**
 * Performs full mixing.  Assumes that the strings primarily comprise dummy
 * values. The children then comprise all non dummy values from the parents.
 * Where both parents have non-dummy values set at the same position child1
 * takes the value from parent1 and child2 takes the value from parent2.
 * @param parent2
 * @param child1
 * @param child2
 */
void StringChromosome<int, IntegerStringChromosomePolicy>::fullMixing(
    const IntegerStringChromosome &parent2, IntegerStringChromosome &child1,
    IntegerStringChromosome &child2) const {
  int *c1 = child1.getString();
  int *c2 = child2.getString();
  int *p1 = getString();
  int *p2 = parent2.getString();

  for (int i = 0; i < length; i++, c1++, c2++, p1++, p2++) {
    if (*p1 == -1 && *p2 == -1)
      *c1 = *c2 = -1;
    else if (*p1 != -1 && *p2 == -1)
      *c1 = *c2 = *p1;
    else if (*p1 == -1 && *p2 != -1)
      *c1 = *c2 = *p2;
    else {
      *c1 = *p1;
      *c2 = *p2;
    }
  }
}

/**
 * Assumes that the strings primarily comprise dummy values. Proceed as for
 * 1 point crossover on integer strings. However if the child is given a
 * dummy value when the other parent has a non-dummy value at a position,
 * then that value is copied to the child.
 * @param parent2
 * @param child1
 * @param child2
 */
void StringChromosome<int, IntegerStringChromosomePolicy>::fullMixingAndCrossover(
    const IntegerStringChromosome &parent2, IntegerStringChromosome &child1,
    IntegerStringChromosome &child2) const {
  // select crossover site
  int site = getRng().randomInt(0, length);

  bool switchFlag = false;
  if (getChromosomePolicy().isAllowSwitch()) switchFlag = getRng().randomBoolean();

  int *c1 = switchFlag ? child2.string.get() : child1.string.get();
  int *c2 = switchFlag ? child1.string.get() : child2.string.get();
  int *p1 = string.get();
  int *p2 = parent2.string.get();

  // create child before cross point
  int i = 0;
  for (; i < site; i++, p1++, p2++, c1++, c2++) {
    if (*p1 == -1 && *p2 != -1)
      *c1 = *c2 = *p2;
    else if (*p1 != -1 && *p2 == -1)
      *c1 = *c2 = *p1;
    else {
      *c2 = *p2;
      *c1 = *p1;
    }
  }

  // child after cross point
  for (; i < length; i++, p1++, p2++, c1++, c2++) {
    if (*p1 == -1 && *p2 != -1)
      *c1 = *c2 = *p2;
    else if (*p1 != -1 && *p2 == -1)
      *c1 = *c2 = *p1;
    else {
      *c2 = *p1;
      *c1 = *p2;
    }
  }
}


}  // namespace GapeGa
