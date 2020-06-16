#include "AtomicContrib.h"
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <eigen3/Eigen/Dense>
#include <GraphMol/Descriptors/AtomicEnvironmentVector.h>
#include <typeinfo>

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
      this->d_weights.push_back(ArrayXXd::Random(384, 160));
      this->d_weights.push_back(ArrayXXd::Random(160, 128));
      this->d_weights.push_back(ArrayXXd::Random(128, 96));
      this->d_weights.push_back(ArrayXXd::Random(96, 1));
      break;
    case 1:
      this->d_weights.push_back(ArrayXXd::Random(384, 144));
      this->d_weights.push_back(ArrayXXd::Random(144, 112));
      this->d_weights.push_back(ArrayXXd::Random(112, 96));
      this->d_weights.push_back(ArrayXXd::Random(96, 1));
      break;
    case 2:
      this->d_weights.push_back(ArrayXXd::Random(384, 128));
      this->d_weights.push_back(ArrayXXd::Random(128, 112));
      this->d_weights.push_back(ArrayXXd::Random(112, 96));
      this->d_weights.push_back(ArrayXXd::Random(96, 1));
      break;
    case 3:
      this->d_weights.push_back(ArrayXXd::Random(384, 144));
      this->d_weights.push_back(ArrayXXd::Random(144, 112));
      this->d_weights.push_back(ArrayXXd::Random(112, 96));
      this->d_weights.push_back(ArrayXXd::Random(96, 1));
      break;
    case -1:
      break;

    default:
      throw ValueErrorException("Invalid values");
  }
}

double ANIAtomContrib::forwardProp(ArrayXXd aev) const {
  auto vec = aev.matrix();
  std::vector<ArrayXXd> layerOuts;
  layerOuts.push_back((aev.matrix() * this->d_weights[0].matrix()).array());
  Utils::CELU(layerOuts[0], 0.1);
  for (unsigned int layer = 1; layer < this->d_weights.size(); layer++) {
    layerOuts.push_back(
        (layerOuts[layer - 1].matrix() * this->d_weights[layer].matrix())
            .array());
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

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
   // #endif