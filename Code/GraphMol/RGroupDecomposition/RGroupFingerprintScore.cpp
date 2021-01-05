//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RGroupFingerprintScore.h"
#include "GraphMol/Fingerprints/Fingerprints.h"
#include "GraphMol//Fingerprints/MorganFingerprints.h"
#include "../../../External/GA/util/Util.h"
#include <memory>
#include <vector>
#include <map>
#include <mutex>
#include <thread>

// #define DEBUG

namespace RDKit {

static const int fingerprintSize = 512;
static const bool useTopologicalFingerprints = false;

// Add fingerprint information to RGroupData
void addFingerprintToRGroupData(RGroupData *rgroupData) {
  if (rgroupData->fingerprint == nullptr) {
    RWMol mol(*rgroupData->combinedMol);
    for (auto atom : mol.atoms()) {
      // replace attachment atom by Boron
      // TODO- Handle multiple attachments differently?
      if (atom->getAtomicNum() == 0) {
        atom->setAtomicNum(5);
        if (atom->getIsotope() > 0) atom->setIsotope(0);
      }
    }
    try {
      MolOps::sanitizeMol(mol);
    } catch (MolSanitizeException &) {
      BOOST_LOG(rdWarningLog)
          << "Failed to sanitize RGroup fingerprint mol for "
          << rgroupData->smiles << std::endl;
    }
#ifdef DEBUG
    std::cerr << "Fingerprint mol smiles " << MolToSmiles(mol) << std::endl;
#endif
    auto fingerprint = useTopologicalFingerprints
                           ? RDKFingerprintMol(mol, 1, 7, fingerprintSize)
                           : MorganFingerprints::getFingerprintAsBitVect(
                                 mol, 2, fingerprintSize);
    fingerprint->getOnBits(rgroupData->fingerprintOnBits);
    rgroupData->fingerprint = std::unique_ptr<ExplicitBitVect>(fingerprint);

#ifdef DEBUG
    std::cerr << "Combined mol smiles " << MolToSmiles(*rgroupData->combinedMol)
              << std::endl;
#endif
  }
}

// calculate centroid of fingerprints
static std::vector<double> fingerprintCenter(
    const std::vector<ExplicitBitVect *> &fingerprints) {
  std::vector<double> center(fingerprintSize, 0.0);
  for (auto fingerprint : fingerprints) {
    for (int i = 0; i < fingerprintSize; i++) {
      if (fingerprint->operator[](i)) {
        center[i] += 1.0;
      }
    }
  }

  for (int i = 0; i < fingerprintSize; i++) {
    center[i] /= fingerprints.size();
  }

#ifdef DEBUG
  std::cerr << "Center " << GarethUtil::collectionToString(center, ",")
            << std::endl;
#endif
  return center;
}

// calculate Euclidean distance between center and fingerprint
static double euclideanDistance(const std::vector<double> &center,
                                const ExplicitBitVect *fingerprint) {
  double sqrSum = 0.0;
  for (int i = 0; i < fingerprintSize; i++) {
    double bit = fingerprint->operator[](i) ? 1.0 : 0.0;
    double d = center[i] - bit;
    sqrSum += d * d;
  }
  auto distance = sqrt(sqrSum);

#ifdef DEBUG
  std::cerr << "Distance " << distance << "sqrSum " << sqrSum << std::endl;
#endif
  return distance;
}

// TODO Profile fingerprintDistanceScore
// Fingerprint score based on distance to fingerprint centroid for rgroups at
// each label Quite slow
double fingerprintDistanceScore(
    const std::vector<size_t> &permutation,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels) {
#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Fingerprint Scoring permutation "
            << " num matches: " << matches.size() << std::endl;
#endif

  std::vector<double> labelScores;
  labelScores.reserve(labels.size());

  // For each label (group)
  for (int l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif
    std::vector<ExplicitBitVect *> fingerprints;
    fingerprints.reserve(permutation.size());

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg != matches[m][permutation[m]].rgroups.end()) {
        auto rgroupData = rg->second;
        if (rgroupData->fingerprint == nullptr) {
          addFingerprintToRGroupData(rgroupData.get());
        }
        auto fingerprint = rgroupData->fingerprint.get();
        fingerprints.push_back(fingerprint);
      }
    }

    auto center = fingerprintCenter(fingerprints);
    auto distanceSum =
        std::accumulate(fingerprints.begin(), fingerprints.end(), 0.0,
                        [&center](double ss, ExplicitBitVect *fp) {
                          return ss + euclideanDistance(center, fp);
                        });
    auto rmsDistance = sqrt(distanceSum / fingerprints.size());
#ifdef DEBUG
    std::cerr << "RMS distance to center " << rmsDistance
              << " sum of distances " << distanceSum << std::endl;
#endif

    labelScores.push_back(rmsDistance);
  }

  // arithmetic mean of scores for each label
  auto sum =
      std::accumulate(labelScores.begin(), labelScores.end(), 0.0,
                      [](double sum, double score) { return sum + score; });
  auto score = sum / labelScores.size();

  // want to minimize this score
#ifdef DEBUG
  std::cerr << "Fingerprint score " << score << std::endl;
#endif
  return -score;
}

// Adds or subtracts a molecule match to the rgroup fingerprint bit counts
// vectors
void modifyVarianceData(
    int matchNumber, int permutationNumber,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels,
    std::map<int, std::shared_ptr<VarianceDataForLabel>> &labelsToVarianceData,
    bool add) {
  // For each label (group)
  for (int l : labels) {
    auto match = matches[matchNumber][permutationNumber].rgroups;
    auto rg = match.find(l);
    if (rg != match.end()) {
      auto rgroupData = rg->second;
      std::shared_ptr<VarianceDataForLabel> variableDataForLabel;
      auto df = labelsToVarianceData.find(l);
      if (df == labelsToVarianceData.end()) {
        variableDataForLabel = std::make_shared<VarianceDataForLabel>(l);
        labelsToVarianceData.emplace(l, variableDataForLabel);
      } else {
        variableDataForLabel = df->second;
      }
      if (add) {
        variableDataForLabel->addRgroupData(rgroupData.get());
      } else {
        variableDataForLabel->removeRgroupData(rgroupData.get());
      }
    }
  }
}

// Adds a molecule match to the rgroup fingerprint bit counts
// vectors
void addVarianceData(int matchNumber, int permutationNumber,
                     const std::vector<std::vector<RGroupMatch>> &matches,
                     const std::set<int> &labels,
                     std::map<int, std::shared_ptr<VarianceDataForLabel>>
                         &labelsToVarianceData) {
  modifyVarianceData(matchNumber, permutationNumber, matches, labels,
                     labelsToVarianceData, true);
}

// Subtracts a molecule match from the rgroup fingerprint bit counts
// vectors
void removeVarianceData(int matchNumber, int permutationNumber,
                        const std::vector<std::vector<RGroupMatch>> &matches,
                        const std::set<int> &labels,
                        std::map<int, std::shared_ptr<VarianceDataForLabel>>
                            &labelsToVarianceData) {
  modifyVarianceData(matchNumber, permutationNumber, matches, labels,
                     labelsToVarianceData, false);
}

// fingerprint variance score
// The arithmetic mean of the mean fingerprint bit variances for the
// fingerprints at each rgroup position.
double fingerprintVarianceScore(
    const std::vector<size_t> &permutation,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels,
    std::map<int, std::shared_ptr<VarianceDataForLabel>>
        *labelsToVarianceData) {
#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Fingerprint Scoring permutation "
            << " num matches: " << matches.size() << std::endl;
#endif

  std::map<int, std::shared_ptr<VarianceDataForLabel>> bitCountsByLabel;
  if (!labelsToVarianceData) {
    labelsToVarianceData = &bitCountsByLabel;
  }

  // For each label (group)
  for (int l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif

    std::shared_ptr<VarianceDataForLabel> variableDataForLabel;
    auto d = labelsToVarianceData->find(l);
    if (d == labelsToVarianceData->end()) {
      variableDataForLabel = std::make_shared<VarianceDataForLabel>(l);
      labelsToVarianceData->emplace(l, variableDataForLabel);
    } else {
      variableDataForLabel = d->second;
    }

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      auto rg = matches[m][permutation[m]].rgroups.find(l);
      if (rg != matches[m][permutation[m]].rgroups.end()) {
        auto rgroupData = rg->second;
        variableDataForLabel->addRgroupData(rgroupData.get());
      }
    }
  }

  return fingerprintVarianceGroupScore(*labelsToVarianceData);
}

// calculates fingerprint variance score from rgroup bit counts
double fingerprintVarianceGroupScore(
    const std::map<int, std::shared_ptr<VarianceDataForLabel>>
        &bitCountsByLabel) {
  // arithmetic mean of scores for each label
#ifdef DEBUG
  std::cerr << "fingerprint variance score: ";
#endif
  auto sum = std::accumulate(
      bitCountsByLabel.cbegin(), bitCountsByLabel.cend(), 0.0,
      [](double sum,
         std::pair<int, std::shared_ptr<VarianceDataForLabel>> pair) {
        auto variance = pair.second->variance();
#ifdef DEBUG
        std::cerr << variance << ',';
#endif
        return sum + variance;
      });
  auto score = sum / bitCountsByLabel.size();
#ifdef DEBUG
  std::cerr << " sum " << sum << " score " << score << std::endl;
#endif
  // want to minimize this score
  return -score;
}

VarianceDataForLabel::VarianceDataForLabel(const int &label,
                                           int numberFingerprints,
                                           const std::vector<int> &bitCounts)
    : label(label),
      numberFingerprints(numberFingerprints),
      bitCounts(bitCounts) {}

VarianceDataForLabel::VarianceDataForLabel(const int &label) : label(label) {
  numberFingerprints = 0;
  bitCounts = std::vector<int>(fingerprintSize, 0.0);
}

static std::mutex groupMutex;

// add an rgroup structure to a bit counts array
void VarianceDataForLabel::addRgroupData(RGroupData *rgroupData) {
  if (rgroupData->fingerprint == nullptr) {
    const std::lock_guard<std::mutex> lock(groupMutex);
    if (rgroupData->fingerprint == nullptr) {
      addFingerprintToRGroupData(rgroupData);
    }
  }
  ++numberFingerprints;
  const auto &onBits = rgroupData->fingerprintOnBits;
  for (int b : onBits) {
    ++bitCounts[b];
  }
}

// remove an rgroup structure to a bit counts array
void VarianceDataForLabel::removeRgroupData(RGroupData *rgroupData) {
  if (rgroupData->fingerprint == nullptr) {
    addFingerprintToRGroupData(rgroupData);
  }
  --numberFingerprints;
  const auto &onBits = rgroupData->fingerprintOnBits;
  for (int b : onBits) {
    --bitCounts[b];
  }
}

// calculate the mean variance for a bit counts array
double VarianceDataForLabel::variance() const {
  auto lambda = [this](double sum, int bitCount) {
    if (bitCount == 0) return sum;
    // variance calculation because fingerprint is binary:
    // sum  == squared sum == bit count
    // ss = sqrSum - (sum * sum) / cnt;
    auto ss = bitCount - (bitCount * bitCount) / (double)numberFingerprints;
    double variancePerBit = ss / (double)numberFingerprints;
#ifdef DEBUG
    std::cerr << variancePerBit << ',';
#endif

    return sum + variancePerBit;
  };

#ifdef DEBUG
  std::cerr << "Bitcounts " << GarethUtil::collectionToString(bitCounts, ",")
            << std::endl;
  std::cerr << "Variance per bit ";
#endif
  auto totalVariance =
      std::accumulate(bitCounts.begin(), bitCounts.end(), 0.0, lambda);
#ifdef DEBUG
  std::cerr << std::endl;
#endif
  auto rmsVariance = sqrt(totalVariance);
#ifdef DEBUG
  std::cerr << "Total Variance " << totalVariance << " RMS Variance "
            << rmsVariance << std::endl;
#endif
  return rmsVariance;
}

}  // namespace RDKit
