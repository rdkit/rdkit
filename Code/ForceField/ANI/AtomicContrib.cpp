#include "AtomicContrib.h"
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <eigen3/Eigen/Dense>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <fstream>
using namespace Eigen;

// #ifdef RDK_HAS_EIGEN3
namespace ForceFields {
namespace ANI {
ANIAtomContrib::ANIAtomContrib(ForceField *owner, int atomType,
                               unsigned int atomIdx, VectorXi speciesVec,
                               int numAtoms) {
  PRECONDITION(owner, "Bad Owner")
  PRECONDITION(atomType == -1 || atomType == 0 || atomType == 1 ||
                   atomType == 2 || atomType == 3,
               "Incorrect Atom Type");
  URANGE_CHECK(atomIdx, owner->positions().size());

  dp_forceField = owner;
  this->d_atomType = atomType;
  this->d_atomIdx = atomIdx;
  this->d_speciesVec = speciesVec;
  this->d_numAtoms = numAtoms;

  switch (d_atomType) {
    case 0:
      for (size_t i = 0; i < 4; i++) {
        this->d_weights.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_weights[i]), 0, "weight", i, 'H');
        this->d_biases.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_biases[i]), 0, "bias", i, 'H');
      }
      break;
    case 1:
      for (size_t i = 0; i < 4; i++) {
        this->d_weights.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_weights[i]), 0, "weight", i, 'C');
        this->d_biases.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_biases[i]), 0, "bias", i, 'C');
      }
      break;
    case 2:
      for (size_t i = 0; i < 4; i++) {
        this->d_weights.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_weights[i]), 0, "weight", i, 'N');
        this->d_biases.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_biases[i]), 0, "bias", i, 'N');
      }
      break;
    case 3:
      for (size_t i = 0; i < 4; i++) {
        this->d_weights.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_weights[i]), 0, "weight", i, 'O');
        this->d_biases.push_back(ArrayXXd::Zero(1, 1));
        Utils::loadFromCSV(&(this->d_biases[i]), 0, "bias", i, 'O');
      }
      break;
    case -1:
      break;

    default:
      throw ValueErrorException("Invalid values");
  }
}

double ANIAtomContrib::forwardProp(ArrayXXd aev) const {
  std::vector<ArrayXXd> layerOuts;
  if (aev.rows() != 1) {
    aev.transposeInPlace();
  }
  layerOuts.push_back(
      ((this->d_weights[0].matrix() * aev.matrix().transpose()).array() +
       this->d_biases[0])
          .transpose());
  Utils::CELU(layerOuts[0], 0.1);
  for (unsigned int layer = 1; layer < this->d_weights.size(); layer++) {
    layerOuts.push_back(((this->d_weights[layer].matrix() *
                          layerOuts[layer - 1].matrix().transpose())
                             .array() +
                         this->d_biases[layer])
                            .transpose());
    if (layer < this->d_weights.size() - 1) Utils::CELU(layerOuts[layer], 0.1);
  }
  auto size = layerOuts.size();
  return layerOuts[size - 1](0, 0);
}

double ANIAtomContrib::getEnergy(double *pos) const {
  auto aev = RDKit::Descriptors::ANI::AtomicEnvironmentVector(
      pos, this->d_speciesVec, this->d_numAtoms);
  return this->ANIAtomContrib::forwardProp(aev.row(this->d_atomIdx));
}

double ANIAtomContrib::getEnergy(Eigen::ArrayXXd aev) const {
  return this->ANIAtomContrib::forwardProp(aev.row(this->d_atomIdx));
}

void ANIAtomContrib::getGrad(double *pos, double *grad) const {}

namespace Utils {

double RELU(double val) { return std::max((double)0, val); }

double coeffMin(double val) { return std::min(double(0), val); }

void CELU(ArrayXXd &input, double alpha) {
  input = input.unaryExpr(&RELU) +
          (alpha * ((input / alpha).exp() - 1)).unaryExpr(&coeffMin);
}

std::vector<std::string> tokenize(const std::string &s) {
  boost::char_separator<char> sep(", \n\r\t");
  boost::tokenizer<boost::char_separator<char>> tok(s, sep);
  std::vector<std::string> tokens;
  std::copy(tok.begin(), tok.end(),
            std::back_inserter<std::vector<std::string>>(tokens));
  return tokens;
}

void loadFromCSV(ArrayXXd *param, unsigned int model, std::string type,
                 unsigned int layer, char atomType) {
  std::string path = getenv("RDBASE");
  std::string paramFile = path + "/Code/ForceField/ANI/Params/model" +
                          std::to_string(model) + "/" + atomType + "_" +
                          std::to_string(layer) + "_" + type;

  std::ifstream instrmSF(paramFile);
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

  param->resize(weight.size(), cols);

  for (size_t i = 0; i < weight.size(); i++) {
    for (size_t j = 0; j < weight[i].size(); j++) {
      (*param)(i, j) = weight[i][j];
    }
  }
}

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
   // #endif