/**
 *
 * StringChromosome.h
 *
 * A template class for GA array based (e.g. bool, int) chromosomes.
 *
 * A ChromosomePolicy is required.  It's contract includes a initialization,
 * mutate and allowSwitch methods.
 *
 * 4/13 Gareth Jones
 */

#ifndef STRINGCHROMOSOME_H_
#define STRINGCHROMOSOME_H_

#include "../util/RandomUtil.h"
#include <cassert>
#include <string>
#include <sstream>
#include <typeinfo>
#include <memory>

namespace GapeGa {

using namespace GarethUtil;

template <typename T, typename ChromosomePolicy>
class StringChromosomeBase {
 private:
  RandomUtil &rng;
  ChromosomePolicy &chromosomePolicy;
  StringChromosomeBase(const StringChromosomeBase &other);
  StringChromosomeBase &operator=(const StringChromosomeBase &other);

 protected:
  const int length;
  std::unique_ptr<T[]> string;

 public:
  StringChromosomeBase(int length_, RandomUtil &rng_,
                       ChromosomePolicy &chromosomePolicy_)
      : rng(rng_),
        chromosomePolicy(chromosomePolicy_),
        length(length_),
        string(new T[length_]) {
    ;
  }

  virtual ~StringChromosomeBase() {}

  void initialize();
  bool equals(const StringChromosomeBase &other) const;
  void copyGene(const StringChromosomeBase &other);
  void mutate();
  void twoPointCrossover(const StringChromosomeBase &parent2,
                         StringChromosomeBase &child1,
                         StringChromosomeBase &child2) const;
  void onePointCrossover(const StringChromosomeBase &parent2,
                         StringChromosomeBase &child1,
                         StringChromosomeBase &child2) const;
  void fullMixing(const StringChromosomeBase &parent2,
                  StringChromosomeBase &child1,
                  StringChromosomeBase &child2) const;
  void fullMixingAndCrossover(const StringChromosomeBase &parent2,
                              StringChromosomeBase &child1,
                              StringChromosomeBase &child2) const;
  std::string geneInfo() const;
  const T getValue(int pos) const;
  T *const getString() const;

};

/**
 * Initializes the chromosome to new values chosen by the policy
 */
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::initialize() {
  for (int i = 0; i < length; i++) {
    string[i] = chromosomePolicy.initialize(i);
  }
}

/**
 *
 * @param other
 * @return true if two chromosome strings are the same
 */
template <typename T, typename ChromosomePolicy>
bool StringChromosomeBase<T, ChromosomePolicy>::equals(
    const StringChromosomeBase &other) const {
  if (other.length != length) {
    return false;
  }
  T *ptr = string.get();
  T *otherPtr = other.string.get();
  for (int i = 0; i < length; i++, ptr++, otherPtr++) {
    if (*ptr != *otherPtr) {
      return false;
    }
  }
  return true;
}

/**
 * copies the gene to another chromosome
 *
 * @param other
 */
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::copyGene(
    const StringChromosomeBase &other) {
  assert(length == other.length);
  assert(typeid(other.string) == typeid(string));
  T *ptr = string.get(), *otherPtr = other.string.get();
  for (int i = 0; i < length; i++, ptr++, otherPtr++) {
    *ptr = *otherPtr;
  }
}

/**
 * Randomly mutates at least one position in the string
 */
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::mutate() {
  double pMutate = 1.0 / (length);
  bool mutated = false;

  T *ptr = string.get();
  for (int i = 0; i < length; i++, ptr++) {
    if (rng.normalRand() < pMutate) {
      T newVal = chromosomePolicy.mutate(i, *ptr);
      *ptr = newVal;
      mutated = true;
    }
  }

  if (!mutated) {
    mutate();
  }
}

/**
 * Performs two point (or circular) crossover
 *
 * @param parent2
 * @param child1
 * @param child2
 */
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::twoPointCrossover(
    const StringChromosomeBase &parent2, StringChromosomeBase &child1,
    StringChromosomeBase &child2) const {
  // choose cross point sites
  int site1 = rng.randomInt(0, length);
  int site2 = rng.randomInt(0, length - 1);
  if (site2 >= site1)
    site2++;
  else {
    int n = site1;
    site1 = site2;
    site2 = n;
  }

  bool switchFlag = false;
  if (chromosomePolicy.isAllowSwitch()) switchFlag = rng.randomBoolean();

  T *c1 = switchFlag ? child2.string : child1.string;
  T *c2 = switchFlag ? child1.string : child2.string;
  T *p1 = string;
  T *p2 = parent2.string;

  // create children before 1st cross point
  int pos = 0;
  for (; pos < site1; pos++, c1++, c2++, p1++, p2++) {
    *c1 = *p1;
    *c2 = *p2;
  }

  // children between cross points
  for (; pos < site2; pos++, c1++, c2++, p1++, p2++) {
    *c1 = *p2;
    *c2 = *p1;
  }

  // children after 2nd cross point
  for (; pos < length; pos++, c1++, c2++, p1++, p2++) {
    *c1 = *p1;
    *c2 = *p2;
  }
}

/**
 * Performs one point crossover
 *
 * @param parent2
 * @param child1
 * @param child2
 */
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::onePointCrossover(
    const StringChromosomeBase &parent2, StringChromosomeBase &child1,
    StringChromosomeBase &child2) const {
  // choose cross point site
  int site = rng.randomInt(0, length - 1);
  bool switchFlag = false;
  if (chromosomePolicy.isAllowSwitch()) switchFlag = rng.randomBoolean();

  T *c1 = switchFlag ? child2.string.get() : child1.string.get();
  T *c2 = switchFlag ? child1.string.get() : child2.string.get();
  T *p1 = string.get();
  T *p2 = parent2.string.get();

  // create children before cross point
  int pos = 0;
  for (; pos < site; pos++, c1++, c2++, p1++, p2++) {
    *c1 = *p1;
    *c2 = *p2;
  }

  // children after cross point
  for (; pos < length; pos++, c1++, c2++, p1++, p2++) {
    *c1 = *p2;
    *c2 = *p1;
  }
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
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::fullMixing(
    const StringChromosomeBase &parent2, StringChromosomeBase &child1,
    StringChromosomeBase &child2) const {
  T *c1 = child1.string;
  T *c2 = child2.string;
  T *p1 = string;
  T *p2 = parent2.string;

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
template <typename T, typename ChromosomePolicy>
void StringChromosomeBase<T, ChromosomePolicy>::fullMixingAndCrossover(
    const StringChromosomeBase &parent2, StringChromosomeBase &child1,
    StringChromosomeBase &child2) const {
  // select crossover site
  int site = rng.randomInt(0, length);

  bool switchFlag = false;
  if (chromosomePolicy.isAllowSwitch()) switchFlag = rng.randomBoolean();

  T *c1 = switchFlag ? child2.string.get() : child1.string.get();
  T *c2 = switchFlag ? child1.string.get() : child2.string.get();
  T *p1 = string.get();
  T *p2 = parent2.string.get();

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

/**
 *
 * @return  a string representation of the chromosome
 *
 */
template <typename T, typename ChromosomePolicy>
std::string StringChromosomeBase<T, ChromosomePolicy>::geneInfo() const {
  std::stringstream ss;
  for (int i = 0; i < length; i++) {
    if (i > 0) ss << ' ';
    ss << string[i];
  }
  return ss.str();
}

/**
 *
 * @param pos position on the chromosome
 * @return  the value on the chromosome at the position
 */
template <typename T, typename ChromosomePolicy>
const T StringChromosomeBase<T, ChromosomePolicy>::getValue(int pos) const {
  return string[pos];
}

template <typename T, typename ChromosomePolicy>
T *const StringChromosomeBase<T, ChromosomePolicy>::getString() const {
  return string.get();
}

}  // namespace GapeGa
#endif /* STRINGCHROMOSOME_H_ */
