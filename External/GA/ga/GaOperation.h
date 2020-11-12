///
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef GAOPERATION_H_
#define GAOPERATION_H_

#include <vector>

namespace GapeGa {

template <typename Chromosome>
class GaOperation {
 private:
  const size_t nParents, nChildren;
  const double weight;
  void (*opfunction)(const std::vector<std::shared_ptr<Chromosome>>& parents,
                     std::vector<std::shared_ptr<Chromosome>>& children);
  GaOperation(const GaOperation&) = delete;
  GaOperation& operator=(const GaOperation& other) = delete;

 public:
  GaOperation(int nParents_, int nChildren_, double weight_,
              void (*opfunction_)(
                  const std::vector<std::shared_ptr<Chromosome>>& parents,
                  std::vector<std::shared_ptr<Chromosome>>& children))
      : nParents(nParents_),
        nChildren(nChildren_),
        weight(weight_),
        opfunction(opfunction_) {}
  virtual ~GaOperation() {}

  size_t getnChildren() const { return nChildren; }

  size_t getnParents() const { return nParents; }

  double getWeight() const { return weight; }

  void (*getOpfunction())(
      const std::vector<std::shared_ptr<Chromosome>>& parents,
      std::vector<std::shared_ptr<Chromosome>>& children) {
    return opfunction;
  }
};

// Templates for string based chromosomes

template <typename T>
void mutateOperation(const std::vector<std::shared_ptr<T>>& parents,
                     std::vector<std::shared_ptr<T>>& children) {
  auto parent = parents[0];
  auto child = children[0];
  child->copyGene(*parent);
  child->mutate();
}

template <typename T>
void onePointCrossoverOperation(const std::vector<std::shared_ptr<T>>& parents,
                                std::vector<std::shared_ptr<T>>& children) {
  auto parent1 = parents[0];
  auto child1 = children[0];
  auto parent2 = parents[1];
  auto child2 = children[1];

  parent1->onePointCrossover(*parent2, *child1, *child2);
}

template <typename T>
void twoPointCrossoverOperation(const std::vector<std::shared_ptr<T>>& parents,
                                std::vector<std::shared_ptr<T>>& children) {
  auto parent1 = parents[0];
  auto child1 = children[0];
  auto parent2 = parents[1];
  auto child2 = children[1];

  parent1->twoPointCrossover(*parent2, *child1, *child2);
}

template <typename T>
void fullMixingOperation(const std::vector<std::shared_ptr<T>>& parents,
                                std::vector<std::shared_ptr<T>>& children) {
  auto parent1 = parents[0];
  auto child1 = children[0];
  auto parent2 = parents[1];
  auto child2 = children[1];

  parent1->fullMixing(*parent2, *child1, *child2);
}

template <typename T>
void fullMixingAndCrossoverOperation(const std::vector<std::shared_ptr<T>>& parents,
                                std::vector<std::shared_ptr<T>>& children) {
  auto parent1 = parents[0];
  auto child1 = children[0];
  auto parent2 = parents[1];
  auto child2 = children[1];

  parent1->fullMixingAndCrossover(*parent2, *child1, *child2);
}

}  // namespace GapeGa

#endif /* GAOPERATION_H_ */
