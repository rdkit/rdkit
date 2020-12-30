//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*
 * A template for providing a steady state no-dupliates GA with
 * linear-normalized parent selection.
 *
 * Requires valid Chromosome and PopulationPolicy classes.
 */

#ifndef LINKEDPOPLINEARSEL_H_
#define LINKEDPOPLINEARSEL_H_

#include <memory>
#include <map>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <cmath>
#include <boost/format.hpp>
#include <limits>

#include "GaOperation.h"
#include "../util/RandomUtil.h"

// Don't include Reporter for RDKit build as it does not play well with Windows
// dynamic DLL
#ifdef INCLUDE_REPORTER
#include "../util/Reporter.h"
#endif

namespace GapeGa {

using namespace GarethUtil;

template <typename Chromosome, typename PopulationPolicy>
class LinkedPopLinearSel {
 private:
  PopulationPolicy& populationPolicy;
  RandomUtil& rng;
  const size_t popsize;
  const std::vector<std::shared_ptr<GaOperation<Chromosome>>> operations;
  const double selectionPressure;

  std::multimap<double, std::shared_ptr<Chromosome>> population;
  std::vector<std::shared_ptr<Chromosome>> freeChromosomes;
  int nOperations = 0;
  int nNiche = 0;
  int nFail = 0;
  int nDuplicates = 0;
  int nAdded = 0;
  double totalScaledFitness, scaledFitnessStep;
  double totalOperatorWeights;

  static constexpr double SELECT_START = 5000.0;

  LinkedPopLinearSel(const LinkedPopLinearSel& other);
  LinkedPopLinearSel& operator=(const LinkedPopLinearSel& other);

  bool addToPopulation(std::shared_ptr<Chromosome>& chromosome);
  bool addToPopulation(std::multimap<double, std::shared_ptr<Chromosome>>& pop,
                       std::shared_ptr<Chromosome>& chromosome);

  std::shared_ptr<Chromosome>& selectParent();
  const typename std::multimap<double,
                               std::shared_ptr<Chromosome>>::const_iterator
  findExactMatch(Chromosome& c) const;
  double bestScore = -std::numeric_limits<double>::max();

 public:
  LinkedPopLinearSel(PopulationPolicy& populationPolicy_);

  virtual ~LinkedPopLinearSel(){};
  void create();
  void iterate();
  void rebuild();
  std::string info() const;
  const std::shared_ptr<Chromosome>& getBest() const;
  const std::vector<std::shared_ptr<Chromosome>> getTiedBest(
      double tolerance = 1e-6) const;
  std::string populationInfo() const;
  double getBestScore() const { return bestScore; }
};

/**
 * Constructor.  Initializes state from the policy.  Determines steady state
 * parameters
 *
 *
 * @param populationPolicy_
 */
template <typename Chromosome, typename PopulationPolicy>
LinkedPopLinearSel<Chromosome, PopulationPolicy>::LinkedPopLinearSel(
    PopulationPolicy& populationPolicy_)
    : populationPolicy(populationPolicy_),
      rng(populationPolicy.getRng()),
      popsize(populationPolicy.getPopsize()),
      operations(populationPolicy.getOperations()),
      selectionPressure(populationPolicy.getSelectionPressure()),
      population(),
      freeChromosomes() {
  totalOperatorWeights = 0;
  for (auto& operation : operations) {
    totalOperatorWeights += operation->getWeight();
  }

  // selection pressure is defined as the ratio between the fitness of the worst
  // individual to the best individual.

  if (selectionPressure < 1.0) {
    throw std::invalid_argument(
        (boost::format("LinkedPopLinearSel: Selection pressure should be "
                       "greater than 1.0 [current value is %|.2f|]") %
         selectionPressure)
            .str());
  }

  double endFitness = selectionPressure * SELECT_START;
  totalScaledFitness = popsize * ((SELECT_START + endFitness) / 2.0);
  scaledFitnessStep =
      (((2.0 * totalScaledFitness) / popsize) - (2.0 * SELECT_START)) /
      (popsize - 1);

  double predictTotalScaledFitness = .0, currentFitness = SELECT_START;
  for (size_t i = 0; i < popsize; i++) {
    predictTotalScaledFitness += currentFitness;
    currentFitness += scaledFitnessStep;
  }
#ifdef INCLUDE_REPORTER
  REPORT(Reporter::TRACE) << "totalScaledFitness " << totalScaledFitness
                          << " predictTotalScaledFitness "
                          << predictTotalScaledFitness;
#endif
  assert(abs(totalScaledFitness - predictTotalScaledFitness) < 1.0e-5);

  double predictEndFitness = SELECT_START + (popsize - 1.0) * scaledFitnessStep;
#ifdef INCLUDE_REPORTER
  REPORT(Reporter::TRACE) << "endFitness " << endFitness
                          << " predictEndFitness " << predictEndFitness;
#endif
  assert(abs(endFitness - predictEndFitness) < 1.0e-10);

  freeChromosomes.reserve(10);
}

/**
 * Builds the initial random population
 */
template <typename Chromosome, typename PopulationPolicy>
void LinkedPopLinearSel<Chromosome, PopulationPolicy>::create() {
  while (population.size() < popsize) {
    std::shared_ptr<Chromosome> chromosome =
        populationPolicy.createChromosome();
    chromosome->initialize();
    chromosome->score();
    addToPopulation(chromosome);
  }

#ifdef INCLUDE_REPORTER
  REPORT(Reporter::DETAIL) << "Created population: " << info();
#endif
}

/**
 * Rebuilds the population
 */
template <typename Chromosome, typename PopulationPolicy>
void LinkedPopLinearSel<Chromosome, PopulationPolicy>::rebuild() {
  std::multimap<double, std::shared_ptr<Chromosome>> newPopulation;
  for (auto& entry : population) {
    auto chromosome = entry.second;
    chromosome->score();
    addToPopulation(newPopulation, chromosome);
  }

  population.swap(newPopulation);
#ifdef INCLUDE_REPORTER
  REPORT(Reporter::DETAIL) << "Rebuilt population: " << info();
#endif
}

/**
 * Selects a parent using roulette wheel parent selection.
 *
 * @return
 */
template <typename Chromosome, typename PopulationPolicy>
std::shared_ptr<Chromosome>&
LinkedPopLinearSel<Chromosome, PopulationPolicy>::selectParent() {
  double val = rng.normalRand() * totalScaledFitness;
  double sum = SELECT_START, currentFitness = SELECT_START;
  auto iterator = population.begin();
  for (size_t i = 0; i < popsize; i++) {
    if (val <= sum) return iterator->second;
    currentFitness += scaledFitnessStep;
    sum += currentFitness;
    ++iterator;
  }
  // may get here because of numeric errors
  assert(abs(sum - totalScaledFitness) < (0.000001 * scaledFitnessStep));
  iterator = population.end();
  --iterator;
  return iterator->second;
}

/**
 * Applies a single operation.
 */
template <typename Chromosome, typename PopulationPolicy>
void LinkedPopLinearSel<Chromosome, PopulationPolicy>::iterate() {
  thread_local std::vector<std::shared_ptr<Chromosome>> parents, children;
  parents.clear();
  children.clear();

  // select an operator.
  double total = 0, val = rng.normalRand() * totalOperatorWeights;
  std::shared_ptr<GaOperation<Chromosome>> selectedOperation = nullptr;
  for (auto& operation : operations) {
    total += operation->getWeight();
    if (val <= total) {
      selectedOperation = operation;
      break;
    }
  }
  if (selectedOperation == nullptr) {
    // may get here because of numeric errors
    assert((total - totalOperatorWeights) < 0.000001);
    selectedOperation = operations.back();
  }

  for (size_t i = 0; i < selectedOperation->getnChildren(); i++) {
    std::shared_ptr<Chromosome> child = nullptr;
    if (freeChromosomes.empty()) {
      child = populationPolicy.createChromosome();
    } else {
      child = freeChromosomes.back();
      freeChromosomes.pop_back();
    }
    children.push_back(child);
  }

  for (size_t i = 0; i < selectedOperation->getnParents(); i++) {
    std::shared_ptr<Chromosome> parent = selectParent();
#ifdef INCLUDE_REPORTER
    REPORT(Reporter::TRACE) << "Parent " << i << ": " << parent->info();
#endif
    parents.push_back(parent);
  }

  (*selectedOperation->getOpfunction())(parents, children);

  int i = 0;
  for (auto& child : children) {
#ifdef INCLUDE_REPORTER
    REPORT(Reporter::TRACE) << "Child  " << i << ": " << child->info();
#endif
    child->score();
    addToPopulation(child);
    i++;
  }

  nOperations++;
#ifdef INCLUDE_REPORTER
  REPORT(Reporter::TRACE) << "Finished iteration " << nOperations;
#endif
}

/**
 * Adds a new individual to the population
 *
 * @param chromosome
 * @return
 */
template <typename Chromosome, typename PopulationPolicy>
bool LinkedPopLinearSel<Chromosome, PopulationPolicy>::addToPopulation(
    std::shared_ptr<Chromosome>& chromosome) {
  return addToPopulation(population, chromosome);
}

/**
 * Adds a new individual to the specified population.  A different population is
 * required during rebuild.
 *
 * @param pop
 * @param chromosome
 * @return
 */
template <typename Chromosome, typename PopulationPolicy>
bool LinkedPopLinearSel<Chromosome, PopulationPolicy>::addToPopulation(
    std::multimap<double, std::shared_ptr<Chromosome>>& pop,
    std::shared_ptr<Chromosome>& chromosome) {
  if (!chromosome->isOk()) {
#ifdef INCLUDE_REPORTER
    REPORT(Reporter::TRACE) << "Bad chromosome not adding to population";
#endif
    nFail++;
    freeChromosomes.push_back(chromosome);
    return false;
  }
  const typename std::multimap<
      double, std::shared_ptr<Chromosome>>::const_iterator match =
      findExactMatch(*chromosome);
  // don't add if we have an exact match
  if (match != population.end()) {
#ifdef INCLUDE_REPORTER
    REPORT(Reporter::TRACE) << "Found exact match not adding to population";
#endif
    nDuplicates++;
    freeChromosomes.push_back(chromosome);
    return false;
  }

  // TODO niche matching

  if (pop.size() >= popsize) {
    // remove the worst individual- this is at the head of the list
    auto worst = pop.begin()->second;
    pop.erase(pop.begin());
#ifdef INCLUDE_REPORTER
    REPORT(Reporter::TRACE) << "Removing the worst individual";
#endif
    freeChromosomes.push_back(worst);
  }

  // add the new one
  double fitness = chromosome->getFitness();
  std::pair<double, std::shared_ptr<Chromosome>> pair(fitness, chromosome);
  pop.insert(pair);
#ifdef INCLUDE_REPORTER
  REPORT(Reporter::TRACE) << "Inserted a new chromosome, fitness " << fitness
                          << " popsize " << pop.size() << " "
                          << chromosome->info();
#endif

  nAdded++;

  if (fitness > bestScore) {
    boost::format format = boost::format("Op %5d new best: ") % nOperations;
    bestScore = fitness;
#ifdef INCLUDE_REPORTER
    REPORT(Reporter::DETAIL) << format << getBest()->info();
#endif
  }
  return true;
}

/**
 * @return an informational string about the state of the population.
 */
template <typename Chromosome, typename PopulationPolicy>
std::string LinkedPopLinearSel<Chromosome, PopulationPolicy>::info() const {
  boost::format format =
      boost::format("Op %5d Fit %7.3f Added %5d Dups %5d Niche %5d Fail %4d") %
      nOperations % bestScore % nAdded % nDuplicates % nNiche % nFail;
  return format.str();
}

/**
 * Returns the best individual in the population
 *
 * @return
 */
template <typename Chromosome, typename PopulationPolicy>
const std::shared_ptr<Chromosome>&
LinkedPopLinearSel<Chromosome, PopulationPolicy>::getBest() const {
  auto iter = population.end();
  --iter;
  return iter->second;
}

template <typename Chromosome, typename PopulationPolicy>
const std::vector<std::shared_ptr<Chromosome>>
LinkedPopLinearSel<Chromosome, PopulationPolicy>::getTiedBest(
    double tolerance) const {
  std::vector<std::shared_ptr<Chromosome>> ties;

  auto iter = population.end();
  --iter;
  ties.push_back(iter->second);
  auto bestScore = iter->first;
  --iter;
  while (std::fabs(iter->first - bestScore) < tolerance) {
    ties.push_back(iter->second);
    if (iter == population.begin()) {
      break;
    }
    --iter;
  }

  return ties;
}

/**
 * Creates a string with summary information about the population
 * @return
 */
template <typename Chromosome, typename PopulationPolicy>
std::string LinkedPopLinearSel<Chromosome, PopulationPolicy>::populationInfo()
    const {
  std::string rtn;
  int no = 1;
  boost::format format = boost::format("%3d: ");

  for (auto iterator = population.begin(); iterator != population.end();
       ++iterator) {
    format % no;
    rtn += format.str() + iterator->second->info() + "\n";
    no++;
  }
  return rtn;
}

/**
 * @return true if this chromosome is already present in the population
 */
template <typename Chromosome, typename PopulationPolicy>
const typename std::multimap<double,
                             std::shared_ptr<Chromosome>>::const_iterator
LinkedPopLinearSel<Chromosome, PopulationPolicy>::findExactMatch(
    Chromosome& c) const {
  std::pair<std::multimap<char, int>::iterator,
            std::multimap<char, int>::iterator>
      ret;
  // any chromosome in the population with the same gene must have the same
  // fitness
  auto iterators = population.equal_range(c.getFitness());
  for (auto iterator = iterators.first; iterator != iterators.second;
       ++iterator) {
    const std::shared_ptr<Chromosome>& other = iterator->second;
    if (c.equals(*other)) return iterator;
  }
  return population.end();
}

}  // namespace GapeGa

#endif /* LINKEDPOPLINEARSEL_H_ */
