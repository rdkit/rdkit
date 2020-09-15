//
//  Copyright (C) 2020 Manan Goel, James Stevenson
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "AtomicContrib.h"

#include <ForceField/ForceField.h>
#include <GraphMol/PeriodicTable.h>
#include <Numerics/EigenSerializer/EigenSerializer.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

#include <Eigen/Dense>
#include <fstream>
using namespace Eigen;

namespace ForceFields {
namespace ANI {

// Model metaparameters
// Currently constants, not loaded from files
// Suitable for ANI-1x, ANI-1cc, ANI-2x,
//    and Schrodinger-ANI.
const double R_Rc = 5.2;
const double R_eta = 16.0;
const double A_Rc = 3.5;
const double A_eta = 8.0;
const double A_zeta = 32.0;
const int len_R_Rs = 16;
const double R_Rs[len_R_Rs] = {
    0.9,  1.16875, 1.4375, 1.70625, 1.975, 2.24375, 2.5125, 2.78125,
    3.05, 3.31875, 3.5875, 3.85625, 4.125, 4.39375, 4.6625, 4.93125};
const int len_A_thetas = 8;
const double A_thetas[len_A_thetas] = {0.19634954, 0.58904862, 0.9817477,
                                       1.3744468,  1.7671459,  2.1598449,
                                       2.552544,   2.9452431};
const int len_A_Rs = 4;
const double A_Rs[len_A_Rs] = {0.9, 1.55, 2.2, 2.85};
const double radial_prefactor = 0.25;
const double inner_product_prefactor = 0.95;

ANIAtomContrib::ANIAtomContrib(ForceField *owner, VectorXi &speciesVec,
                               std::string model) {
  PRECONDITION(owner, "This ForceField contrib needs an owner")
  dp_forceField = owner;
  // load element mapping and atomic self-energies
  this->paramDir = model;
  std::ifstream selfEnergyFile(paramDir + "/selfEnergies");
  if (!selfEnergyFile.is_open()) {  // try default location
    std::string rdbase = getenv("RDBASE");
    paramDir = rdbase + "/Code/ForceField/ANI/Params/" + model;
    selfEnergyFile.open(paramDir + "/selfEnergies");
  }
  PRECONDITION(selfEnergyFile.is_open(),
               "Missing ANI energy file " + paramDir + "/selfEnergies")

  std::map<int, std::string> elements_by_index;
  std::map<std::string, int> indices_by_element;
  // parse atomic energy lines of the form "H,0=-0.60095298000"
  std::string next;
  while (!selfEnergyFile.eof()) {
    std::getline(selfEnergyFile, next, ',');
    std::string element = next;
    std::getline(selfEnergyFile, next, '=');
    int index = atoi(next.c_str());
    std::getline(selfEnergyFile, next);
    double E = atof(next.c_str());
    elements_by_index[index] = element;
    indices_by_element[element] = index;
    this->d_selfEnergies.push_back(E);
  }
  selfEnergyFile.close();
  n_elements = elements_by_index.size();
  n_features =
      n_elements * ((n_elements - 1) / 2.0 + 1) * len_A_Rs * len_A_thetas +
      n_elements * len_R_Rs;

  // make d_atomTypes using indices_by_element and speciesVec
  auto table = RDKit::PeriodicTable::getTable();
  this->d_atomTypes.resize(speciesVec.size());
  for (long i = 0; i < speciesVec.size(); i++) {
    std::string element = table->getElementSymbol(speciesVec[i]);
    if (indices_by_element.find(element) == indices_by_element.end()) {
      throw ValueErrorException("Element not Supported: " + element);
      return;
    }
    this->d_atomTypes[i] = indices_by_element[element];
  }

  // load parameters by model number in ensemble, by element, and by layer
  int model_i = 0;
  while (true) {
    // get ensembleSize by counting files in param dir
    std::string paramFile =
        paramDir + "/model" + std::to_string(model_i) + ".bin";
    FILE *fp = fopen(paramFile.c_str(), "rb");  // check that file exists
    fclose(fp);
    if (!fp) {
      break;
    }
    this->d_weights.resize(model_i + 1);
    this->d_biases.resize(model_i + 1);
    this->d_weights[model_i].resize(n_elements);
    this->d_biases[model_i].resize(n_elements);
    for (int element_i = 0; element_i < n_elements; element_i++) {
      std::string element = elements_by_index[element_i];
      std::vector<MatrixXf> floatWeights, floatBiases;
      // note, deserializeAll is not efficient - should get all elements at once
      RDNumeric::EigenSerializer::deserializeAll(&floatWeights, &floatBiases,
                                                 paramFile, element);
      for (size_t i = 0; i < floatWeights.size(); i++) {
        this->d_weights[model_i][element_i].push_back(
            floatWeights[i].cast<double>());
        this->d_biases[model_i][element_i].push_back(
            floatBiases[i].cast<double>());
      }
    }
    model_i++;
  }
  this->d_ensembleSize = model_i;
  PRECONDITION(this->d_ensembleSize, "No ANI model files found in " + paramDir)

  /*
  TODO: load AEV metaparameters from files?
  Hyperparameters not yet included in these .bin files:
    const double radial_prefactor = 0.25;
    const double inner_product_prefactor = 0.95;
    const double R_Rc = 5.2;
    const double A_Rc = 3.5;

  std::string paramFilePath = paramDir + "/AEVParams/";
  ArrayXd ShfR;
  RDNumeric::EigenSerializer::deserialize(ShfR, paramFilePath + "ShfR.bin");
  ArrayXd EtaR;
  RDNumeric::EigenSerializer::deserialize(EtaR, paramFilePath + "EtaR.bin");
  ArrayXd ShfZ;
  RDNumeric::EigenSerializer::deserialize(ShfZ, paramFilePath + "ShfZ.bin");
  ArrayXd ShfA;
  RDNumeric::EigenSerializer::deserialize(ShfA, paramFilePath + "ShfA.bin");
  ArrayXd zeta;
  RDNumeric::EigenSerializer::deserialize(zeta, paramFilePath + "zeta.bin");
  ArrayXd etaA;
  RDNumeric::EigenSerializer::deserialize(etaA, paramFilePath + "etaA.bin");
  std::cout << "ShfR" << ShfR << "\n";  // R_Rs
  std::cout << "EtaR" << EtaR << "\n";  // R_eta
  std::cout << "ShfZ" << ShfZ << "\n";  // A_thetas
  std::cout << "ShfA" << ShfA << "\n";  // A_Rs
  std::cout << "zeta" << zeta << "\n";  // A_zeta
  std::cout << "etaA" << etaA << "\n";  // A_eta
  */
}

// Forward propagation for one atom and one model
// Puts all outputs into layer_values
double ANIAtomContrib::forwardPropOneAtom(MatrixXd &layer_values,
                                          const int modelNo,
                                          const int element_i) const {
  int n_layers = this->d_weights[0][0].size();  // note, assumes fixed n_layers
  for (int layer = 0; layer < n_layers; layer++) {
    // matrix-vector multiply: temp_out = biases + weights * temp_in
    int input_dim = this->d_weights[modelNo][element_i][layer].cols();
    int output_dim = this->d_weights[modelNo][element_i][layer].rows();
    for (int i = 0; i < output_dim; i++) {
      layer_values(layer * 2 + 1, i) =
          this->d_biases[modelNo][element_i][layer](i);
      for (int j = 0; j < input_dim; j++) {
        layer_values(layer * 2 + 1, i) +=
            layer_values.coeff(layer * 2, j) *
            this->d_weights[modelNo][element_i][layer].coeff(i, j);
      }
    }
    if (layer < n_layers - 1) {  // CELU activation, except on last layer
      for (int i = 0; i < output_dim; i++) {
        double x = layer_values(layer * 2 + 1, i);
        layer_values(layer * 2 + 2, i) = x > 0.0 ? x : 0.1 * (exp(x / 0.1) - 1);
      }
    }
  }
}

// Iterate across all atoms, running forward prop on each AEV row
// (note: might be faster to sort by element and run in blocks)
double ANIAtomContrib::forwardProp(MatrixXd &aev) const {
  double E_total = 0.0;
  int n_layers = this->d_weights[0][0].size();  // note, assumes fixed n_layers
  MatrixXd layer_values = MatrixXd::Zero(n_layers * 2, n_features);
  // Run forward and backprop for each atom and each model
  // Save layer values from forward prop for use in backprop
  // Output: mean across ensemble of dE_dAEV
  for (int atom_i = 0; atom_i < this->d_atomTypes.size(); atom_i++) {
    auto element_i = this->d_atomTypes[atom_i];
    for (int feature_i = 0; feature_i < n_features; feature_i++) {
      // layer_0 = features
      layer_values(0, feature_i) = aev(atom_i, feature_i);
    }
    for (size_t modelNo = 0; modelNo < this->d_weights.size(); modelNo++) {
      // forward prop for this atom and model
      this->forwardPropOneAtom(layer_values, modelNo, element_i);
      double E = layer_values(n_layers * 2 - 1, 0);  // last layer is E
      E_total += E;
    }
  }
  E_total /= this->d_weights.size();
  return E_total;
}

// This function takes xyz coords and atom types as inputs, and calculates
// the Atomic Environment Vectors or the gradient dE/dx.
// Calculating dE/dx requires the additional input dE_dAEV, which is
// operated on to output dE/dx using the chain rule (mathematically, we
// apply the Jacobian dAEV_i/dx_j to produce dE/dx).
void calcAEV(
    const double *pos,            // inputs: x y z. shape=(n_atoms,3).
    VectorXi const &atomTypes,    // type index for each atom, shape=(n_atoms)
    Eigen::MatrixXd &out_buffer,  // output: forces dE/dxyz shape=(n_atoms,3)
                                  // or AEVs, shape=(n_atoms,n_features)
    Eigen::MatrixXd const &dE_dAEV,  // Input dE/AEV_i, needed to make dE/dx.
                                     // shape=(n_atoms,n_features).
    const bool output_grad,  // whether to fill out_buffer with AEVs or dE/dx
    const int n_elements     // number of allowed elements (H, C, N, O, etc)
) {
  size_t numAtoms = atomTypes.size();
  // Note: identity of i0 matters: pair (i0,i1) != pair (i1,i0).
  for (size_t i0 = 0; i0 < numAtoms; i0++) {  // iterate over all atoms
    const double *xyz0 = &pos[i0 * 3];
    for (size_t i1 = 0; i1 < numAtoms; i1++) {  // iterate over all pairs
      if (i0 == i1) continue;                   // skip if same atom
      const double *xyz1 = &pos[i1 * 3];
      double x01[3] = {xyz1[0] - xyz0[0], xyz1[1] - xyz0[1], xyz1[2] - xyz0[2]};
      double r01 = sqrt(x01[0] * x01[0] + x01[1] * x01[1] + x01[2] * x01[2]);
      if (r01 > R_Rc) continue;  // skip pair if outside radial cutoff
      // radial part
      int element1_slot = atomTypes(i1) * len_R_Rs;
      double r01_scaled = r01 * M_PI / R_Rc;
      double fc = 0.5 * cos(r01_scaled) + 0.5;  // note, assumes r01 < R_Rc
      for (int R_i = 0; R_i < len_R_Rs; R_i++) {
        double R = R_Rs[R_i];
        double r_offset = r01 - R;
        double expR = exp(-R_eta * r_offset * r_offset);
        double radial_feature = radial_prefactor * expR * fc;
        int feature_i = element1_slot + R_i;
        if (!output_grad) {
          out_buffer(i0, feature_i) += radial_feature;
          continue;  // stop here, don't calculate gradient
        }
        // x, y, and z derivatives
        double fcd = 0.5 * M_PI * sin(r01_scaled) / (r01 * R_Rc);
        double r01_etad = 2 * R_eta * r_offset / r01;
        double expRd = r01_etad * expR;
        double d_radial_feature = radial_prefactor * (expRd * fc + expR * fcd);
        for (int nd = 0; nd < 3; nd++) {
          double dEdx = dE_dAEV(i0, feature_i) * x01[nd] * d_radial_feature;
          out_buffer(i0, nd) += dEdx;  // dE/dF_i0 * dF_i0/dx_i0
          out_buffer(i1, nd) -= dEdx;  // dE/dF_i0 * dF_i0/dx_i1
        }
      }
      if (r01 > A_Rc) continue;  // skip if outside angular cutoff
      // Calculate angular part
      // Iterate over triplets: we can skip duplicate (i1,i2)
      //    entries here because i1 and i2 are symmetric.
      for (size_t i2 = i1 + 1; i2 < numAtoms; i2++) {
        if (i0 == i2) continue;  // skip if same atom
        const double *xyz2 = &pos[i2 * 3];
        double x02[3] = {xyz2[0] - xyz0[0], xyz2[1] - xyz0[1],
                         xyz2[2] - xyz0[2]};
        double r02 = sqrt(x02[0] * x02[0] + x02[1] * x02[1] + x02[2] * x02[2]);
        if (r02 > A_Rc) continue;  // skip if outside angular cutoff
        // find cos theta from dot(A,B) / (|A|*|B|)
        double dot012 = x01[0] * x02[0] + x01[1] * x02[1] + x01[2] * x02[2];
        double cos_theta012 = dot012 / (r01 * r02) * inner_product_prefactor;
        double sin_theta012 = sqrt(1 - cos_theta012 * cos_theta012);
        double r01_scaled = r01 * (M_PI / A_Rc);
        double fc1 = 0.5 * cos(r01_scaled) + 0.5;
        double r02_scaled = r02 * (M_PI / A_Rc);
        double fc2 = 0.5 * cos(r02_scaled) + 0.5;
        double fc = fc1 * fc2;
        double mean_r012 = 0.5 * (r01 + r02);
        // sort, element1 <= element2
        int element1 = std::min(atomTypes(i1), atomTypes(i2));
        int element2 = std::max(atomTypes(i1), atomTypes(i2));
        int element_slot =
            (n_elements * (n_elements - 1) / 2 -
             (n_elements - element1) * (n_elements - element1 - 1) / 2 +
             element2) *
                len_A_Rs * len_A_thetas +
            n_elements * len_R_Rs;
        for (int A_i = 0; A_i < len_A_thetas; A_i++) {
          double cos_A_theta = cos(A_thetas[A_i]);
          double sin_A_theta = sin(A_thetas[A_i]);
          for (int R_i = 0; R_i < len_A_Rs; R_i++) {
            double A_R = A_Rs[R_i];
            double cos_theta012_minus_A =
                cos_theta012 * cos_A_theta + sin_theta012 * sin_A_theta;
            double zeta_arg = (1 + cos_theta012_minus_A) / 2;
            double theta_part = pow(zeta_arg, A_zeta);
            double r_eta_term = mean_r012 - A_R;
            double r_part = exp(-A_eta * r_eta_term * r_eta_term);
            double angular_feature = 2 * theta_part * r_part * fc;
            int feature_i = element_slot + R_i * len_A_thetas + A_i;
            if (!output_grad) {
              out_buffer(i0, feature_i) += angular_feature;
              continue;  // don't calculate gradient
            }
            // x, y, and z derivatives
            double zetad = A_zeta * pow(zeta_arg, A_zeta - 1);
            for (int nd = 0; nd < 3; nd++) {  // iterate over dimensions x, y, z
              // grad dF/dx for central atom, dF/dx0
              double r01d = -x01[nd] / r01;
              double r02d = -x02[nd] / r02;
              double dot012d = -x02[nd] - x01[nd];
              double cos_theta012d =
                  inner_product_prefactor *
                  (dot012d * r01 * r02 - dot012 * (r01d * r02 + r01 * r02d)) /
                  (r01 * r02 * (r01 * r02));
              // Note: don't need to check for sin_theta012 = 0.0,
              // because inner_product_prefactor is 0.95 and
              // thus cos_theta012 < 1.0
              double sin_theta012d =
                  -cos_theta012d * cos_theta012 / sin_theta012;
              double r01_scaledd = M_PI * r01d / A_Rc;
              double cos_theta012_minus_Ad =
                  cos_A_theta * cos_theta012d + sin_A_theta * sin_theta012d;
              double fc1d = -0.5 * r01_scaledd * sin(r01_scaled);
              double r02_scaledd = M_PI * r02d / A_Rc;
              double fc2d = -0.5 * r02_scaledd * sin(r02_scaled);
              double fcd = fc1d * fc2 + fc1 * fc2d;
              double theta_partd = cos_theta012_minus_Ad * 0.5 * zetad;
              double r_partd = -A_eta * (r01d + r02d) * r_eta_term * r_part;
              double dFdx0 =
                  2 * ((theta_partd * r_part + theta_part * r_partd) * fc +
                       theta_part * r_part * fcd);
              out_buffer(i0, nd) +=
                  dE_dAEV(i0, feature_i) * dFdx0;  // dE/dF_i0 * dF_i0/dx_i0

              // grad dF/dx terms for the other two atoms, dF/dx1 and dF/dx2
              r01d = x01[nd] / r01;
              dot012d = x02[nd];
              cos_theta012d = inner_product_prefactor *
                              (dot012d * r01 - dot012 * r01d) /
                              (r01 * r01 * r02);
              sin_theta012d = -cos_theta012 * cos_theta012d / sin_theta012;
              r01_scaledd = M_PI * r01d / A_Rc;
              fc1d = -0.5 * r01_scaledd * sin(r01_scaled);
              cos_theta012_minus_Ad =
                  cos_A_theta * cos_theta012d + sin_A_theta * sin_theta012d;
              fcd = fc2 * fc1d;
              theta_partd = cos_theta012_minus_Ad * 0.5 * zetad;
              r_partd = -A_eta * r01d * r_eta_term * r_part;
              double dFdx1 =
                  2 * ((theta_partd * r_part + theta_part * r_partd) * fc +
                       theta_part * r_part * fcd);
              out_buffer(i1, nd) +=
                  dE_dAEV(i0, feature_i) * dFdx1;  // dE/dF_i0 * dF_i0/dx_i1

              r02d = x02[nd] / r02;
              dot012d = x01[nd];
              cos_theta012d = inner_product_prefactor *
                              (dot012d * r02 - dot012 * r02d) /
                              (r01 * r02 * r02);
              sin_theta012d = -cos_theta012 * cos_theta012d / sin_theta012;
              r02_scaledd = M_PI * r02d / A_Rc;
              fc2d = -0.5 * r02_scaledd * sin(r02_scaled);
              cos_theta012_minus_Ad =
                  cos_A_theta * cos_theta012d + sin_A_theta * sin_theta012d;
              fcd = fc1 * fc2d;
              theta_partd = cos_theta012_minus_Ad * 0.5 * zetad;
              r_partd = -A_eta * r02d * r_eta_term * r_part;
              double dFdx2 =
                  2 * ((theta_partd * r_part + theta_part * r_partd) * fc +
                       theta_part * r_part * fcd);
              out_buffer(i2, nd) +=
                  dE_dAEV(i0, feature_i) * dFdx2;  // dE/dF_i0 * dF_i0/dx_i2
            }
          }
        }
      }
    }
  }
}

double ANIAtomContrib::getEnergy(double *pos) const {
  MatrixXd aev = MatrixXd::Zero(this->d_atomTypes.size(), n_features);
  MatrixXd dummy;
  calcAEV(pos, this->d_atomTypes, aev, dummy, false, n_elements);
  double E_pred = this->ANIAtomContrib::forwardProp(aev);
  for (int atom_i = 0; atom_i < this->d_atomTypes.size(); atom_i++) {
    E_pred += this->d_selfEnergies[this->d_atomTypes[atom_i]];
  }
  return E_pred;
}

void ANIAtomContrib::getGrad(double *pos, double *grad) const {
  MatrixXd aev = MatrixXd::Zero(this->d_atomTypes.size(), n_features);
  MatrixXd dummy;
  calcAEV(pos, this->d_atomTypes, aev, dummy, false, n_elements);
  // iterate across all atoms, running each AEV row individually
  // (note: might be faster to sort by element and run in blocks)
  // dE_dAEV size = (n_atoms, n_features), same size as AEV
  MatrixXd dE_dAEV = MatrixXd::Zero(this->d_atomTypes.size(), n_features);
  int n_layers = this->d_weights[0][0].size();  // note, assumes fixed n_layers
  MatrixXd layer_values = MatrixXd::Zero(n_layers * 2, n_features);
  VectorXd grad_in = VectorXd::Zero(n_features);
  VectorXd grad_out = VectorXd::Zero(n_features);
  // Run forward and backprop for each atom and each model
  // Save layer values from forward prop for use in backprop
  // Output: mean across ensemble of dE_dAEV
  for (int atom_i = 0; atom_i < this->d_atomTypes.size(); atom_i++) {
    auto element_i = this->d_atomTypes[atom_i];
    for (int feature_i = 0; feature_i < n_features; feature_i++) {
      // layer_0 = features
      layer_values(0, feature_i) = aev(atom_i, feature_i);
    }
    for (size_t modelNo = 0; modelNo < this->d_weights.size(); modelNo++) {
      // forward prop for this atom and model
      this->forwardPropOneAtom(layer_values, modelNo, element_i);
      // layer_values now contains all forward pass info for this atom and NN
      // note, E_atom is here if desired: E = layer_values(2*n_layers-1, 0)

      // backprop
      grad_in(0) = 1.0;  // start from dE/dE = 1.0, apply chain rule
      for (int layer = n_layers - 1; layer >= 0; layer--) {
        // matrix-vector multiply: grad_out = grad_in * weights^T
        int input_dim = this->d_weights[modelNo][element_i][layer].cols();
        int output_dim = this->d_weights[modelNo][element_i][layer].rows();
        for (int j = 0; j < input_dim; j++) {
          grad_out(j) = 0.0;  // zero out before sum below
          for (int i = 0; i < output_dim; i++) {
            grad_out(j) +=
                grad_in.coeff(i) *
                this->d_weights[modelNo][element_i][layer].coeff(i, j);
          }
        }
        if (layer > 0) {  // backprop CELU, except for top layer
          for (int j = 0; j < input_dim; j++) {
            double x = layer_values(2 * layer - 1, j);
            double celu_grad = x > 0.0 ? 1.0 : exp(x / 0.1);
            grad_in(j) = grad_out(j) * celu_grad;  // input for next pass
          }
        }
      }
      // assemble mean dE_dAEV
      for (int feature_i = 0; feature_i < n_features; feature_i++) {
        dE_dAEV(atom_i, feature_i) +=
            grad_out(feature_i) / this->d_weights.size();
      }
    }
  }

  // Pass gradient back through calcAEV() with output_grad=True
  // Multiplies dE_dAEV by sparse Jacobian dAEV_i/dx_j
  //   and returns result without putting the Jacobian in memory.
  MatrixXd grad_buffer = MatrixXd::Zero(this->d_atomTypes.size(), 3);  // output
  calcAEV(pos, this->d_atomTypes, grad_buffer, dE_dAEV, true, n_elements);
  // copy from buffer back into float* for output
  for (int atom_i = 0; atom_i < this->d_atomTypes.size(); atom_i++) {
    grad[3 * atom_i] = grad_buffer(atom_i, 0);
    grad[3 * atom_i + 1] = grad_buffer(atom_i, 1);
    grad[3 * atom_i + 2] = grad_buffer(atom_i, 2);
  }
}

}  // namespace ANI
}  // namespace ForceFields