
#include <chrono>
#include "RGroupGa.h"
#include "RGroupDecompData.h"

namespace RDKit {

using namespace std;

RGroupDecompositionChromosome::RGroupDecompositionChromosome(RGroupGa& rGroupGa)
    : IntegerStringChromosome(rGroupGa.chromosomeLength(), rGroupGa.getRng(),
                              rGroupGa.getChromosomePolicy()),
      rGroupGa(rGroupGa) {
  permutation.reserve(rGroupGa.numberDecompositions());
}

std::string RGroupDecompositionChromosome::info() const {
  auto format = boost::format("Fit %7.3f : ") % fitness;
  return format.str() + geneInfo();
}

double RGroupDecompositionChromosome::score() {
  auto &rGroupData = rGroupGa.getRGroupData();
  RGroupScore scoreMethod =
      static_cast<RGroupScore>(rGroupData.params.scoreMethod);
  if (operationName != RgroupMutate) {
    decode();
  }
  if (scoreMethod == FingerprintVariance && labelsToVarianceData.size() > 0 && operationName == RgroupMutate) {
    fitness = fingerprintVarianceGroupScore(labelsToVarianceData);
    // assert(fitness == recalculateScore());
  } else {
    fitness = rGroupData.score(permutation, &labelsToVarianceData);
  }
  return fitness;
}

double RGroupDecompositionChromosome::recalculateScore() {
  std::cerr << "Recalculating score" << std::endl;
  auto &rGroupData = rGroupGa.getRGroupData();
  return rGroupData.score(permutation);
}

void RGroupDecompositionChromosome::decode() {
  auto values = getString();
  permutation.clear();
  auto &matches = rGroupGa.getRGroupData().matches;
  auto pos = 0;
  for (const auto &m : matches) {
    if (m.size() == 1) {
      permutation.push_back(0);
    } else {
      permutation.push_back(values[pos]);
      pos++;
    }
  }
}

void RGroupDecompositionChromosome::copyGene(
    const StringChromosomeBase<int, IntegerStringChromosomePolicy>& other) {
  StringChromosomeBase<int, IntegerStringChromosomePolicy>::copyGene(other);
  const auto& parent = static_cast<const RGroupDecompositionChromosome&>(other);
  copyVarianceData(parent.labelsToVarianceData, labelsToVarianceData);
}

RGroupGa::RGroupGa(const RGroupDecompData& rGroupData)
    : rGroupData(rGroupData),
      chromosomePolicy(getRng(), rGroupData.matches.size()) {
  setSelectionPressure(1.0001);
  // setSelectionPressure(1.001);

  const auto &matches = rGroupData.matches;
  auto pos = 0;
  for (auto m : matches) {
    if (m.size() == 1) continue;
    chromosomePolicy.setMax(pos, m.size());
    pos++;
  }
  chromLength = pos;
  numberDecomps = matches.size();

  // TODO refine these settings
  auto popsize = 100 + chromLength / 10;
  if (popsize > 200) popsize = 200;
  noIterations = 5000 + chromLength * 100;
  if (noIterations > 100000) noIterations = 100000;
  noIterations = 1000000;

  // profiler settings
  // if (popsize > 100) popsize = 100;
  // if (noIterations > 200) noIterations = 2000;

  setPopsize(popsize);
}

void RGroupGa::rGroupMutateOperation(
    const std::vector<std::shared_ptr<RGroupDecompositionChromosome>>& parents,
    std::vector<std::shared_ptr<RGroupDecompositionChromosome>>& children) {
  auto parent = parents[0];
  auto child = children[0];
  child->copyGene(*parent);
  // double pMutate = parent->getLength() > 100 ? .01 : -1.0;
  child->mutate();
  child->setOperationName(RgroupMutate);
  child->decode();

  auto& labelsToVarianceData = child->getLabelsToVarianceData();
  if (labelsToVarianceData.size() == 0) return;
#ifdef DEBUG
  std::cerr << "RGroup mutate start" << std::endl;
#endif
#ifdef DEBUG
  std::cerr << "Starting child score" << std::endl;
  fingerprintVarianceScore(labelsToVarianceData);
#endif

  auto &parentPermutation = parent->getPermutation();
  auto &childPermutation = child->getPermutation();
  const auto& rgroupData = parent->getRGroupGA().getRGroupData();
  const auto& matches = rgroupData.matches;
  const auto& labels = rgroupData.labels;

  for (auto pos = 0U; pos < parentPermutation.size(); pos++) {
    int parentValue = parentPermutation.at(pos);
    int childValue = childPermutation.at(pos);
    if (parentValue != childValue) {
      removeVarianceData(pos, parentValue, matches, labels,
                         labelsToVarianceData);
#ifdef DEBUG
      std::cerr << "After removing parent" << std::endl;
      fingerprintVarianceScore(labelsToVarianceData);
#endif
      addVarianceData(pos, childValue, matches, labels, labelsToVarianceData);
#ifdef DEBUG
      std::cerr << "After adding child" << std::endl;
      fingerprintVarianceScore(labelsToVarianceData);
#endif
    }
  }
#ifdef DEBUG
  std::cerr << "Final recalculating" << std::endl;
  child->recalculateScore();
  std::cerr << "RGroup mutate done" << std::endl;
#endif
}

void RGroupGa::rGroupCrossoverOperation(
    const std::vector<std::shared_ptr<RGroupDecompositionChromosome>>& parents,
    std::vector<std::shared_ptr<RGroupDecompositionChromosome>>& children) {
  auto parent1 = parents[0];
  auto child1 = children[0];
  auto parent2 = parents[1];
  auto child2 = children[1];

  child1->setOperationName(Crossover);
  child2->setOperationName(Crossover);
  clearVarianceData(child1->getLabelsToVarianceData());
  clearVarianceData(child2->getLabelsToVarianceData());

  parent1->twoPointCrossover(*parent2, *child1, *child2);
}

void RGroupGa::createOperations() {
  // bias to mutation as that operator is so efficient
  auto mutationOperation =
      make_shared<GaOperation<RGroupDecompositionChromosome>>(
          1, 1, 75.0, &rGroupMutateOperation);
  auto crossoverOperation =
      make_shared<GaOperation<RGroupDecompositionChromosome>>(
          2, 2, 25.0, &rGroupCrossoverOperation);
  operations.reserve(2);
  operations.push_back(mutationOperation);
  operations.push_back(crossoverOperation);
}

std::string timeInfo(const chrono::system_clock::time_point start) {
  auto now = chrono::high_resolution_clock::now();
  chrono::duration<double, milli> milliseconds = now - start;
  auto seconds = milliseconds.count() / 1000.0;
  auto format = boost::format("Time %7.2f") % seconds;
  return format.str();
}

vector<vector<size_t>> RGroupGa::run() {
  auto startTime = std::chrono::high_resolution_clock::now();
  createOperations();
  population = unique_ptr<RGroupGaPopulation>(new RGroupGaPopulation(*this));
  auto format = boost::format(
                    "Running GA number operations %5d population size %5d "
                    "chromosome length %5d %s") %
                noIterations % getPopsize() % chromLength % timeInfo(startTime);
  REPORT(Reporter::INFO) << format.str();
  population->create();
  double bestScore = population->getBestScore();
  REPORT(Reporter::INFO) << population->info();
  REPORT(Reporter::DEBUG) << population->populationInfo();

  int nOps = 0;
  int lastImprovementOp = 0;
  while (nOps < noIterations) {
    population->iterate();
    nOps++;
    if (nOps % 1000 == 0) {
      REPORT(Reporter::INFO)
          << population->info() << " " << timeInfo(startTime);
    }
    if (population->getBestScore() > bestScore) {
      bestScore = population->getBestScore();
      lastImprovementOp = nOps;
      auto format = boost::format("OP %5d Fit %7.3f %s") % nOps % bestScore %
                    timeInfo(startTime);
      REPORT(Reporter::INFO) << format.str();
    }
    if (nOps - lastImprovementOp > 5000) {
      REPORT(Reporter::INFO) << "Op " << nOps << " No improvement since "
                             << lastImprovementOp << " finishing..";
      break;
    }
  }
  const shared_ptr<RGroupDecompositionChromosome> best = population->getBest();
  REPORT(Reporter::INFO) << "Best solution " << best->info() << endl;
  REPORT(Reporter::DETAIL) << population->populationInfo();

  auto ties = population->getTiedBest();
  auto permutations =
      GarethUtil::mapToNewList<shared_ptr<RGroupDecompositionChromosome>,
                               vector<size_t>>(
          ties, [](const shared_ptr<RGroupDecompositionChromosome> c) {
            return c->getPermutation();
          });
  REPORT(Reporter::INFO) << "Execution " << timeInfo(startTime);
  return permutations;
}

shared_ptr<RGroupDecompositionChromosome> RGroupGa::createChromosome() {
  return make_shared<RGroupDecompositionChromosome>(*this);
}

void copyVarianceData(const map<int, shared_ptr<VarianceDataForLabel>>& from,
                      map<int, shared_ptr<VarianceDataForLabel>>& to) {
  for (auto it = from.cbegin(); it != from.end(); ++it) {
    auto df = to.find(it->first);
    if (df == to.end()) {
      to.emplace(it->first, make_shared<VarianceDataForLabel>(*it->second));
    } else {
      auto& fromData = it->second;
      auto& toData = df->second;
      toData->numberFingerprints = fromData->numberFingerprints;
      toData->bitCounts.assign(fromData->bitCounts.cbegin(),
                               fromData->bitCounts.cend());
    }
  }
}

void clearVarianceData(map<int, shared_ptr<VarianceDataForLabel>>& data) {
  for (auto it = data.begin(); it != data.end(); ++it) {
    auto d = it->second;
    d->numberFingerprints = 0;
    d->bitCounts.assign(d->bitCounts.size(), 0.0);
  }
}

}  // namespace RDKit