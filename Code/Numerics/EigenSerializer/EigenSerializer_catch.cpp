//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"

#include <RDGeneral/Invariant.h>
#include <Numerics/EigenSerializer/EigenSerializer.h>

TEST_CASE("Eigen Matrix Serialization Test") {
  SECTION("Test 1") {
    Eigen::ArrayXXd J(10, 5);
    J.setRandom();
    RDNumeric::EigenSerializer::serialize(J, "matrix.bin");

    Eigen::ArrayXXd JCopy;
    RDNumeric::EigenSerializer::deserialize(JCopy, "matrix.bin");

    CHECK(J.isApprox(JCopy));
  }
  SECTION("Test serializeAll and deserializeAll") {
    unsigned int numModels = 5;
    for (unsigned int i = 0; i < numModels; i++) {
      std::string fileName = "./model" + std::to_string(i);
      std::vector<std::string> atomTypes = {"H", "C", "N", "O"};
      std::vector<std::pair<
          std::string, std::vector<std::pair<std::string, Eigen::ArrayXXd>>>>
          weightsAndBiasesForEachAtomType;

      for (unsigned int j = 0; j < atomTypes.size(); j++) {
        unsigned int numLayers = 4;
        std::vector<std::pair<std::string, Eigen::ArrayXXd>> weights;
        for (unsigned int k = 0; k < numLayers; k++) {
          std::string weightType = "weight";
          weights.push_back(
              std::make_pair(weightType, Eigen::ArrayXXd::Random(10, 10)));
          std::string biasType = "bias";
          weights.push_back(
              std::make_pair(biasType, Eigen::ArrayXXd::Random(10, 1)));
        }
        weightsAndBiasesForEachAtomType.push_back(
            std::make_pair(atomTypes[j], weights));
      }
      RDNumeric::EigenSerializer::serializeAll(&weightsAndBiasesForEachAtomType,
                                               fileName);

      std::vector<std::pair<
          std::string, std::vector<std::pair<std::string, Eigen::ArrayXXd>>>>
          weightsAndBiasesForEachAtomTypeCopy;
      for (unsigned int j = 0; j < atomTypes.size(); j++) {
        std::vector<std::pair<std::string, Eigen::ArrayXXd>> weightsWithType;
        std::vector<Eigen::ArrayXXd> weights, biases;
        RDNumeric::EigenSerializer::deserializeAll(&weights, &biases, fileName,
                                                   atomTypes[j]);
        for (unsigned int k = 0; k < weights.size(); k++) {
          weightsWithType.push_back(std::make_pair("weight", weights[k]));
          weightsWithType.push_back(std::make_pair("bias", biases[k]));
        }
        weightsAndBiasesForEachAtomTypeCopy.push_back(
            std::make_pair(atomTypes[j], weightsWithType));
      }
      CHECK(weightsAndBiasesForEachAtomType.size() ==
            weightsAndBiasesForEachAtomTypeCopy.size());
      for (unsigned int i = 0; i < weightsAndBiasesForEachAtomType.size();
           i++) {
        CHECK(weightsAndBiasesForEachAtomType[i].first ==
              weightsAndBiasesForEachAtomTypeCopy[i].first);
        auto groundTruthWeights = weightsAndBiasesForEachAtomType[i].second;
        auto readWeights = weightsAndBiasesForEachAtomTypeCopy[i].second;
        CHECK(groundTruthWeights.size() == readWeights.size());
        for (unsigned int j = 0; j < groundTruthWeights.size(); j++) {
          CHECK(groundTruthWeights[j].first == readWeights[j].first);
          CHECK(groundTruthWeights[j].second.isApprox(readWeights[j].second));
        }
      }
    }
  }
}