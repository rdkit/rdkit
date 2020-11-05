//
// Created by gareth on 11/4/20.
//

#ifndef RDKIT_RGROUPFINGERPRINTSCORE_H
#define RDKIT_RGROUPFINGERPRINTSCORE_H

#include "RGroupMatch.h"
#include <vector>

namespace RDKit {

struct VarianceDataForLabel {
  const int label;
  int numberFingerprints;
  std::vector<int> bitCounts;

  VarianceDataForLabel(const int &label, int numberFingerprints,
                       const std::vector<int> &bitCounts);
  VarianceDataForLabel(const int &label);
  VarianceDataForLabel(const VarianceDataForLabel& other) = default;
  VarianceDataForLabel& operator= (const VarianceDataForLabel& other) = delete;
  void addRgroupData(RGroupData *rgroupData);
  void removeRgroupData(RGroupData *rgroupData);
  double variance() const;
};

double fingerprintVarianceScore(
    const std::vector<size_t> &bitCount,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels,
    std::map<int, std::shared_ptr<VarianceDataForLabel>> *labelsToVarianceData =
    nullptr);
double fingerprintDistanceScore(
    const std::vector<size_t> &bitCount,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels);
double fingerprintVarianceGroupScore(
    const std::map<int, std::shared_ptr<VarianceDataForLabel>> &bitCountsByLabel);

void addVarianceData(int matchNumber, int permutationNumber,
                     const std::vector<std::vector<RGroupMatch>> &matches,
                     const std::set<int> &labels,
                     std::map<int, std::shared_ptr<VarianceDataForLabel>>
                     &labelsToVarianceData);
void removeVarianceData(int matchNumber, int permutationNumber,
                        const std::vector<std::vector<RGroupMatch>> &matches,
                        const std::set<int> &labels,
                        std::map<int, std::shared_ptr<VarianceDataForLabel>>
                        &labelsToVarianceData);


}

#endif  // RDKIT_RGROUPFINGERPRINTSCORE_H
