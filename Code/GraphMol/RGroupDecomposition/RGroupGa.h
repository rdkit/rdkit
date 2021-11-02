//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RGROUPDECOMPGA_H_
#define RGROUPDECOMPGA_H_

#include <vector>
#include <memory>
#include <map>
#include <chrono>

#include "../../../External/GA/ga/StringChromosome.h"
#include "../../../External/GA/ga/GaBase.h"
#include "../../../External/GA/ga/GaOperation.h"
#include "../../../External/GA/ga/LinkedPopLinearSel.h"
#include "../../../External/GA/ga/IntegerStringChromosomePolicy.h"
#include "RGroupFingerprintScore.h"
#include "RGroupScore.h"

namespace RDKit {

using namespace GapeGa;
using namespace std;

class RGroupDecompositionChromosome;
class RGroupGa;
struct RGroupDecompData;

typedef LinkedPopLinearSel<RGroupDecompositionChromosome, RGroupGa>
    RGroupGaPopulation;

typedef enum {
  RgroupMutate = 0x01,
  Crossover = 0x02,
  Create = 0x04,
} OperationName;

class RGroupDecompositionChromosome : public IntegerStringChromosome {
 public:
  RGroupDecompositionChromosome(RGroupGa& rGroupGa);

  double getFitness() const { return fitness; }

  OperationName getOperationName() const { return operationName; }

  void setOperationName(OperationName operationName) {
    this->operationName = operationName;
  }

  std::string info() const;

  double score();

  double recalculateScore();

  bool isOk() { return true; }

  void decode();

  void copyGene(const StringChromosomeBase& other) override;

  FingerprintVarianceScoreData& getFingerprintVarianceScoreData() {
    return fingerprintVarianceScoreData;
  }

  const vector<size_t>& getPermutation() const { return permutation; }

  const RGroupGa& getRGroupGA() const { return rGroupGa; }

 private:
  RGroupDecompositionChromosome(const RGroupDecompositionChromosome& other) =
      delete;
  RGroupDecompositionChromosome& operator=(
      const RGroupDecompositionChromosome& other) = delete;
  double fitness;
  FingerprintVarianceScoreData fingerprintVarianceScoreData;
  OperationName operationName = Create;
  vector<size_t> permutation;

  RGroupGa& rGroupGa;
};

struct GaResult {
  RGroupScorer rGroupScorer;

  GaResult(const double score, const vector<vector<size_t>>& permutations)
      : rGroupScorer(RGroupScorer(permutations, score)) {}
  GaResult(const GaResult& other) : rGroupScorer(other.rGroupScorer) {}

  GaResult() {}

  // Copy constructor required by MSVC for future<GaResult>
  GaResult& operator=(const GaResult& other);
};

class RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupGa : public GaBase {
 public:
  RGroupGa(const RGroupDecompData& rGroupData,
           const chrono::steady_clock::time_point* const t0 = nullptr);

  IntegerStringChromosomePolicy& getChromosomePolicy() {
    return chromosomePolicy;
  }

  int chromosomeLength() const { return chromLength; }

  int numberDecompositions() const { return numberDecomps; }

  GaResult run(int runNumber = 1);

  vector<GaResult> runBatch();

  shared_ptr<RGroupDecompositionChromosome> createChromosome();

  const RGroupDecompData& getRGroupData() const { return rGroupData; }

  const vector<shared_ptr<GaOperation<RGroupDecompositionChromosome>>>
  getOperations() const;

  unsigned int numberPermutations() const { return numPermutations; }

 private:
  RGroupGa(const RGroupGa& other) = delete;
  RGroupGa& operator=(const RGroupGa& other) = delete;
  const RGroupDecompData& rGroupData;

  IntegerStringChromosomePolicy chromosomePolicy;
  int numberOperations;
  int numberOperationsWithoutImprovement;
  int chromLength;
  unsigned int numPermutations;
  int numberDecomps;
  const chrono::steady_clock::time_point* const t0;

  void createOperations();

  static void rGroupMutateOperation(
      const std::vector<std::shared_ptr<RGroupDecompositionChromosome>>&
          parents,
      std::vector<std::shared_ptr<RGroupDecompositionChromosome>>& children);
  static void rGroupCrossoverOperation(
      const std::vector<std::shared_ptr<RGroupDecompositionChromosome>>&
          parents,
      std::vector<std::shared_ptr<RGroupDecompositionChromosome>>& children);
};

void copyVarianceData(const FingerprintVarianceScoreData& fromData,
                      FingerprintVarianceScoreData& toData);

void clearVarianceData(
    FingerprintVarianceScoreData& fingerprintVarianceScoreData);

}  // namespace RDKit

#endif  // RGROUPDECOMPGA_H_
