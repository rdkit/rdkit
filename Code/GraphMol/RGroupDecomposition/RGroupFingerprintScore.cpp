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
#include <utility>
#include <vector>
#include <map>
#ifdef RDK_THREADSAFE_SSS
#include <mutex>
#endif

// #define DEBUG

namespace RDKit {

static const int fingerprintSize = 512;
static const bool useTopologicalFingerprints = false;

// TODO scale variance by the number of separate attachments (as that variance
// will be counted for each attachment).

// Add fingerprint information to RGroupData
void addFingerprintToRGroupData(RGroupData *rgroupData) {
  if (rgroupData->fingerprint == nullptr) {
    RWMol mol(*rgroupData->combinedMol);
    for (auto atom : mol.atoms()) {
      // replace attachment atom by Boron
      // TODO- Handle multiple attachments differently?
      if (atom->getAtomicNum() == 0) {
        atom->setAtomicNum(5);
        if (atom->getIsotope() > 0) {
          atom->setIsotope(0);
        }
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

// Adds or subtracts a molecule match to the rgroup fingerprint bit counts
// vectors
void FingerprintVarianceScoreData::modifyVarianceData(
    int matchNumber, int permutationNumber,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels, bool add) {
  // For each label (group)
  const auto &match = matches.at(matchNumber).at(permutationNumber);
  for (int l : labels) {
    auto rg = match.rgroups.find(l);
    if (rg != match.rgroups.end()) {
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
  auto rgroupsMissing = match.numberMissingUserRGroups;
  if (add) {
    numberOfMissingUserRGroups += rgroupsMissing;
    numberOfMolecules++;
  } else {
    numberOfMissingUserRGroups -= rgroupsMissing;
    numberOfMolecules--;
  }
}

// Adds a molecule match to the rgroup fingerprint bit counts
// vectors
void FingerprintVarianceScoreData::addVarianceData(
    int matchNumber, int permutationNumber,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels) {
  modifyVarianceData(matchNumber, permutationNumber, matches, labels, true);
}

// Subtracts a molecule match from the rgroup fingerprint bit counts
// vectors
void FingerprintVarianceScoreData::removeVarianceData(
    int matchNumber, int permutationNumber,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels) {
  modifyVarianceData(matchNumber, permutationNumber, matches, labels, false);
}

// fingerprint variance score
// The arithmetic mean of the mean fingerprint bit variances for the
// fingerprints at each rgroup position.
double fingerprintVarianceScore(
    const std::vector<size_t> &permutation,
    const std::vector<std::vector<RGroupMatch>> &matches,
    const std::set<int> &labels,
    FingerprintVarianceScoreData *fingerprintVarianceScoreData) {
#ifdef DEBUG
  std::cerr << "---------------------------------------------------"
            << std::endl;
  std::cerr << "Fingerprint Scoring permutation "
            << " num matches: " << matches.size() << std::endl;
#endif
  CHECK_INVARIANT(permutation.size() <= matches.size(), "");
  FingerprintVarianceScoreData fingerprintVarianceScoreData2;
  if (!fingerprintVarianceScoreData) {
    fingerprintVarianceScoreData = &fingerprintVarianceScoreData2;
  }
  auto &labelsToVarianceData =
      fingerprintVarianceScoreData->labelsToVarianceData;

  // For each label (group)
  for (int l : labels) {
#ifdef DEBUG
    std::cerr << "Label: " << l << std::endl;
#endif

    std::shared_ptr<VarianceDataForLabel> variableDataForLabel;
    auto d = labelsToVarianceData.find(l);
    if (d == labelsToVarianceData.end()) {
      variableDataForLabel = std::make_shared<VarianceDataForLabel>(l);
      labelsToVarianceData.emplace(l, variableDataForLabel);
    } else {
      variableDataForLabel = d->second;
    }

    for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
      const auto &match = matches[m].at(permutation[m]);
      auto rg = match.rgroups.find(l);
      if (rg != match.rgroups.end()) {
        auto rgroupData = rg->second;
        variableDataForLabel->addRgroupData(rgroupData.get());
      }
    }
  }

  size_t numberMissingRGroups = 0;
  for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
    numberMissingRGroups +=
        matches[m].at(permutation[m]).numberMissingUserRGroups;
  }
  fingerprintVarianceScoreData->numberOfMissingUserRGroups +=
      numberMissingRGroups;
  fingerprintVarianceScoreData->numberOfMolecules += permutation.size();

  return fingerprintVarianceScoreData->fingerprintVarianceGroupScore();
}

// calculates fingerprint variance score from rgroup bit counts
double FingerprintVarianceScoreData::fingerprintVarianceGroupScore() {
  // arithmetic mean of scores for each label
#ifdef DEBUG
  std::cerr << "fingerprint variance score: ";
#endif
  auto sum = std::accumulate(
      labelsToVarianceData.cbegin(), labelsToVarianceData.cend(), 0.0,
      [](double sum,
         std::pair<int, std::shared_ptr<VarianceDataForLabel>> pair) {
        auto variance = pair.second->variance();
    // perhaps here the variance should be weighted by occupancy- so that
    // sparsely populated rgroups are penalized

    // e.g variance *= ((double) numberOfMolecules) /
    // ((double)pair.second->numberFingerprints);
#ifdef DEBUG
        std::cerr << variance << ',';
#endif
        return sum + variance;
      });

  // Heuristic correction of missing user r_groups - equivalent to a variance
  // penalty of 1 for each missing user R-group across the entire dataset
  CHECK_INVARIANT(numberOfMolecules > 0, "No compounds to be scored!");
  double rgroupPenalty =
      (double)numberOfMissingUserRGroups / (double)numberOfMolecules;
  // double the penalty to catch systems like
  // https://github.com/rdkit/rdkit/issues/3896
  auto score = sum + 2.0 * rgroupPenalty;
#ifdef DEBUG
  std::cerr << " sum " << sum << " rgroup penalty " << rgroupPenalty
            << " score " << score << std::endl;
#endif
  // want to minimize this score
  return -score;
}

VarianceDataForLabel::VarianceDataForLabel(const int &label,
                                           int numberFingerprints,
                                           std::vector<int> bitCounts)
    : label(label),
      numberFingerprints(numberFingerprints),
      bitCounts(std::move(bitCounts)) {}

VarianceDataForLabel::VarianceDataForLabel(const int &label) : label(label) {
  numberFingerprints = 0;
  bitCounts = std::vector<int>(fingerprintSize, 0.0);
}

#ifdef RDK_THREADSAFE_SSS
static std::mutex groupMutex;
#endif

// add an rgroup structure to a bit counts array
void VarianceDataForLabel::addRgroupData(RGroupData *rgroupData) {
  if (rgroupData->fingerprint == nullptr) {
#ifdef RDK_THREADSAFE_SSS
    const std::lock_guard<std::mutex> lock(groupMutex);
#endif
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
    if (bitCount == 0) {
      return sum;
    }
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

void FingerprintVarianceScoreData::clear() {
  numberOfMissingUserRGroups = 0;
  numberOfMolecules = 0;
  labelsToVarianceData.clear();
}

}  // namespace RDKit
