//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "AtomicContrib.h"
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <Numerics/EigenSerializer/EigenSerializer.h>
#include <Eigen/Dense>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <fstream>
using namespace Eigen;

namespace ForceFields {
namespace ANI {
ANIAtomContrib::ANIAtomContrib(ForceField *owner, int atomType,
                               unsigned int atomIdx, VectorXi &speciesVec,
                               unsigned int numAtoms, unsigned int numLayers,
                               unsigned int ensembleSize,
                               std::string modelType) {
  PRECONDITION(owner, "Bad Owner")
  PRECONDITION(atomType == 0 || atomType == 1 || atomType == 2 || atomType == 3,
               "Atom Type not Supported");
  PRECONDITION(modelType == "ANI-1x" || modelType == "ANI-1ccx",
               "Model Not currently supported")
  PRECONDITION(ensembleSize > 0,
               "There must be at least 1 model in the ensemble");
  URANGE_CHECK(atomIdx, numAtoms);
  dp_forceField = owner;
  this->d_atomType = atomType;
  this->d_atomIdx = atomIdx;
  this->d_speciesVec = speciesVec;
  this->d_numAtoms = numAtoms;
  this->d_ensembleSize = ensembleSize;
  this->d_modelType = modelType;

  if (this->d_atomEncoding.find(this->d_atomType) !=
      this->d_atomEncoding.end()) {
    auto atomicSymbol = this->d_atomEncoding[this->d_atomType];
    for (unsigned int modelNum = 0; modelNum < ensembleSize; modelNum++) {
      std::vector<MatrixXd> currModelWeights;
      std::vector<MatrixXd> currModelBiases;
      Utils::loadFromBin(&currModelWeights, &currModelBiases, modelNum,
                         atomicSymbol, this->d_modelType);
      this->d_weights.push_back(currModelWeights);
      this->d_biases.push_back(currModelBiases);
    }
    Utils::loadSelfEnergy(&(this->d_selfEnergy), atomicSymbol,
                          this->d_modelType);
  } else {
    this->d_selfEnergy = 0;
  }

  // Different values for means of the gaussian symmetry functions
  std::string path = getenv("RDBASE");
  std::string paramFilePath =
      path + "/Data/ANIParams/" + modelType + "/AEVParams/";

  // Weights for the radial symmetry functions
  ArrayXd ShfR;
  RDNumeric::EigenSerializer::deserialize(ShfR, paramFilePath + "ShfR.bin");
  // Variance terms for the gaussian symmetry functions
  ArrayXd EtaR;
  RDNumeric::EigenSerializer::deserialize(EtaR, paramFilePath + "EtaR.bin");

  // Weights for the angular symmetry functions
  ArrayXd ShfZ;
  RDNumeric::EigenSerializer::deserialize(ShfZ, paramFilePath + "ShfZ.bin");
  ArrayXd ShfA;
  RDNumeric::EigenSerializer::deserialize(ShfA, paramFilePath + "ShfA.bin");
  // distance wise shifts in the distance term of the angular symmetry function

  ArrayXd zeta;
  RDNumeric::EigenSerializer::deserialize(zeta, paramFilePath + "zeta.bin");
  ArrayXd etaA;
  RDNumeric::EigenSerializer::deserialize(etaA, paramFilePath + "etaA.bin");

  this->d_aevParams.insert(std::make_pair("ShfR", ShfR));
  this->d_aevParams.insert(std::make_pair("EtaR", EtaR));
  this->d_aevParams.insert(std::make_pair("ShfZ", ShfZ));
  this->d_aevParams.insert(std::make_pair("ShfA", ShfA));
  this->d_aevParams.insert(std::make_pair("zeta", zeta));
  this->d_aevParams.insert(std::make_pair("etaA", etaA));
}

double ANIAtomContrib::forwardProp(ArrayXXd &aev) const {
  if (this->d_atomType == -1) {
    return 0;
  }

  if (aev.cols() != 1) {
    aev.transposeInPlace();
  }

  MatrixXd aevMat = aev.matrix();

  std::vector<double> energies;
  energies.reserve(this->d_weights.size());
  for (unsigned int modelNo = 0; modelNo < this->d_weights.size(); modelNo++) {
    auto temp = aevMat;
    for (unsigned int layer = 0; layer < this->d_weights[modelNo].size();
         layer++) {
      temp = ((this->d_weights[modelNo][layer] * temp) +
              this->d_biases[modelNo][layer])
                 .eval();
      if (layer < this->d_weights[modelNo].size() - 1) {
        Utils::CELU(temp, 0.1);
      }
    }
    energies.push_back(temp.coeff(0, 0));
  }
  return std::accumulate(energies.begin(), energies.end(), 0.0) /
         energies.size();
}

double ANIAtomContrib::getEnergy(double *pos) const {
  PRECONDITION(pos != nullptr, "Positions array is NULL");
  ArrayXXd aev;
  RDKit::Descriptors::ANI::AtomicEnvironmentVector(
      aev, pos, this->d_speciesVec, this->d_numAtoms, &(this->d_aevParams));
  ArrayXXd row = aev.row(this->d_atomIdx);
  return this->ANIAtomContrib::forwardProp(row) + this->d_selfEnergy;
}

double ANIAtomContrib::getEnergy(Eigen::ArrayXXd &aev) const {
  ArrayXXd row = aev.row(this->d_atomIdx);
  return this->ANIAtomContrib::forwardProp(row) + this->d_selfEnergy;
}

void ANIAtomContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(pos != nullptr, "Positions are null");
  PRECONDITION(grad != nullptr, "Gradient array is null");
  double displacement = 1e-5;

  // + - x movement
  pos[3 * this->d_atomIdx] += displacement;
  auto posXEnergy = this->dp_forceField->calcEnergy(pos);

  pos[3 * this->d_atomIdx] -= 2 * displacement;
  auto negXEnergy = this->dp_forceField->calcEnergy(pos);

  grad[3 * this->d_atomIdx] = (posXEnergy - negXEnergy) / (2 * displacement);
  pos[3 * this->d_atomIdx] += displacement;

  // + - Y movement
  pos[3 * this->d_atomIdx + 1] += displacement;
  auto posYEnergy = this->dp_forceField->calcEnergy(pos);

  pos[3 * this->d_atomIdx + 1] -= 2 * displacement;
  auto negYEnergy = this->dp_forceField->calcEnergy(pos);

  grad[3 * this->d_atomIdx + 1] =
      (posYEnergy - negYEnergy) / (2 * displacement);
  pos[3 * this->d_atomIdx + 1] += displacement;

  // + - Z movement
  pos[3 * this->d_atomIdx + 2] += displacement;
  auto posZEnergy = this->dp_forceField->calcEnergy(pos);

  pos[3 * this->d_atomIdx + 2] -= 2 * displacement;
  auto negZEnergy = this->dp_forceField->calcEnergy(pos);

  grad[3 * this->d_atomIdx + 2] =
      (posZEnergy - negZEnergy) / (2 * displacement);
  pos[3 * this->d_atomIdx + 2] += displacement;
}

namespace Utils {

void CELU(MatrixXd &input, double alpha) {
  input = input.unaryExpr([&](double val) {
    return std::max(0.0, val) +
           std::min(alpha * (std::exp(val / alpha) - 1), 0.0);
  });
}

std::vector<std::string> tokenize(const std::string &s) {
  boost::char_separator<char> sep(", \n\r\t");
  boost::tokenizer<boost::char_separator<char>> tok(s, sep);
  std::vector<std::string> tokens;
  std::copy(tok.begin(), tok.end(),
            std::back_inserter<std::vector<std::string>>(tokens));
  return tokens;
}

void loadFromBin(std::vector<MatrixXd> *weights, unsigned int model,
                 std::string weightType, unsigned int layer,
                 std::string atomType, std::string modelType) {
  PRECONDITION(weights != nullptr, "Vector of weights is NULL");
  std::string path = getenv("RDBASE");
  std::string paramFile = path + "/Data/ANIParams/" + modelType +
                          "/model" + std::to_string(model) + "/" + atomType +
                          "_" + std::to_string(layer) + "_" + weightType +
                          ".bin";
  MatrixXf weight;
  RDNumeric::EigenSerializer::deserialize(weight, paramFile);
  weights->push_back(weight.cast<double>());
}

void loadFromBin(std::vector<MatrixXd> *weights, std::vector<MatrixXd> *biases,
                 unsigned int model, std::string atomType,
                 std::string modelType) {
  PRECONDITION(weights != nullptr, "Weights array is null");
  PRECONDITION(biases != nullptr, "Biases array is null");
  std::string path = getenv("RDBASE");
  std::string paramFile = path + "/Data/ANIParams/" + modelType +
                          "/model" + std::to_string(model) + ".bin";
  std::vector<MatrixXf> floatWeights, floatBiases;
  RDNumeric::EigenSerializer::deserializeAll(&floatWeights, &floatBiases,
                                             paramFile, atomType);
  for (unsigned int i = 0; i < floatWeights.size(); i++) {
    weights->push_back(floatWeights[i].cast<double>());
    biases->push_back(floatBiases[i].cast<double>());
  }
}

void loadFromCSV(std::vector<MatrixXd> *weights, unsigned int model,
                 std::string weightType, unsigned int layer,
                 std::string atomType, std::string modelType) {
  PRECONDITION(weights != nullptr, "Weights array is null");
  std::string path = getenv("RDBASE");
  std::string paramFile = path + "/Data/ANIParams/" + modelType +
                          "/model" + std::to_string(model) + "/" + atomType +
                          "_" + std::to_string(layer) + "_" + weightType;

  std::ifstream instrmSF(paramFile.c_str());
  if (!instrmSF.good()) {
    throw ValueErrorException(paramFile + " Model File does not exist");
    return;
  }
  std::string line;
  std::vector<std::string> tokens;
  std::vector<std::vector<double>> weight;
  unsigned int cols = 1;
  while (!instrmSF.eof()) {
    std::getline(instrmSF, line);
    tokens = tokenize(line);
    std::vector<double> row;
    for (auto v : tokens) {
      std::istringstream os(v);
      double d;
      os >> d;
      row.push_back(d);
    }
    if (row.size() > 0) {
      cols = row.size();
      weight.push_back(row);
    }
  }

  MatrixXd param(weight.size(), cols);

  for (unsigned int i = 0; i < weight.size(); i++) {
    for (unsigned int j = 0; j < weight[i].size(); j++) {
      param(i, j) = weight[i][j];
    }
  }
  weights->push_back(param);
}

void loadSelfEnergy(double *energy, std::string atomType,
                    std::string modelType) {
  std::string path = getenv("RDBASE");
  std::string filePath =
      path + "/Data/ANIParams/" + modelType + "/selfEnergies";

  std::ifstream selfEnergyFile(filePath.c_str());
  if (!selfEnergyFile.good()) {
    throw ValueErrorException(filePath + " : File Does Not Exist");
    return;
  }
  std::string line;
  while (!selfEnergyFile.eof()) {
    std::getline(selfEnergyFile, line);
    boost::char_separator<char> sep(" ,=");
    boost::tokenizer<boost::char_separator<char>> tok(line, sep);
    std::vector<std::string> tokens;
    std::copy(tok.begin(), tok.end(),
              std::back_inserter<std::vector<std::string>>(tokens));

    if (tokens[0] == atomType) {
      std::istringstream os(tokens[2]);
      os >> *energy;
      break;
    }
  }
  selfEnergyFile.close();
}

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields