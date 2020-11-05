#ifndef RGROUPDECOMPGA_H_
#define RGROUPDECOMPGA_H_

#include <vector>
#include <memory>
#include <map>

#include "../../../External/GA/ga/StringChromosome.h"
#include "../../../External/GA/ga/GaBase.h"
#include "../../../External/GA/ga/GaOperation.h"
#include "../../../External/GA/ga/LinkedPopLinearSel.h"
#include "../../../External/GA/ga/IntegerStringChromosomePolicy.h"

namespace RDKit {

using namespace GapeGa;
using namespace std;

class RGroupDecompositionChromosome;
class RGroupGa;
struct RGroupDecompData;
struct VarianceDataForLabel;

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

  map<int, shared_ptr<VarianceDataForLabel>>& getLabelsToVarianceData() {
    return labelsToVarianceData;
  }

  const vector<size_t>& getPermutation() const { return permutation; }

  const RGroupGa& getRGroupGA() const { return rGroupGa; };

 private:
  RGroupDecompositionChromosome(const RGroupDecompositionChromosome& other) =
      delete;
  RGroupDecompositionChromosome& operator=(
      const RGroupDecompositionChromosome& other) = delete;
  double fitness;
  map<int, shared_ptr<VarianceDataForLabel>> labelsToVarianceData;
  OperationName operationName = Create;
  vector<size_t> permutation;

  RGroupGa& rGroupGa;
};

class RGroupGa : public GaBase {
 public:
  RGroupGa(const RGroupDecompData& rGroupData);

  IntegerStringChromosomePolicy& getChromosomePolicy() {
    return chromosomePolicy;
  }

  int chromosomeLength() const { return chromLength; }

  int numberDecompositions() const { return numberDecomps; }

  vector<vector<size_t>> run();

  shared_ptr<RGroupDecompositionChromosome> createChromosome();

  const RGroupDecompData& getRGroupData() const { return rGroupData; };

  const vector<shared_ptr<GaOperation<RGroupDecompositionChromosome>>>&
  getOperations() const {
    return operations;
  }

  double getBestScore() const { return population->getBest()->getFitness(); };

 private:
  RGroupGa(const RGroupGa& other) = delete;
  RGroupGa& operator=(const RGroupGa& other) = delete;
  const RGroupDecompData& rGroupData;
  IntegerStringChromosomePolicy chromosomePolicy;
  vector<shared_ptr<GaOperation<RGroupDecompositionChromosome>>> operations;
  unique_ptr<RGroupGaPopulation> population;
  int noIterations;
  int chromLength;
  int numberDecomps;

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

void copyVarianceData(
    const std::map<int, std::shared_ptr<VarianceDataForLabel>>& from,
    std::map<int, std::shared_ptr<VarianceDataForLabel>>& to);

void clearVarianceData(
    std::map<int, std::shared_ptr<VarianceDataForLabel>>& data);

}  // namespace RDKit

#endif  // RGROUPDECOMPGA_H_
