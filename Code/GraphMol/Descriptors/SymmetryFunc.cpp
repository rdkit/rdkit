#include <GraphMol/RDKitBase.h>
#include <cmath>
#include "SymmetryFunc.h"
#include <eigen3/Eigen/Dense>

using namespace Eigen;

namespace RDKit {
namespace Descriptors {

ArrayXXd cosine_cutoff(ArrayXXd distances, double cutoff) {
  // Cosine cutoff function assuming all distances are less than the cutoff
  return 0.5 * (distances * (M_PI / cutoff)).cos() + 0.5;
}

VectorXd generateSpeciesVector(const ROMol &mol) {
  // Generate atom species vector as mentioned in torchani
  // H : 0
  // C : 1
  // N : 2
  // O : 3
  // All other atoms : -1
  int numAtoms = mol.getNumAtoms();
  VectorXd species(numAtoms);

  for (auto i = 0; i < numAtoms; i++) {
    auto atom = mol.getAtomWithIdx(i);
    if (atom->getAtomicNum() == 1) {
      species[i] = 0;
    }
    else if (atom->getAtomicNum() == 6) {
      species[i] = 1;
    }
    else if (atom->getAtomicNum() == 7) {
      species[i] = 2;
    }
    else if (atom->getAtomicNum() == 8) {
      species[i] = 3;
    }
    else {
      species[i] = -1;
    }
  }

  return species;
}

ArrayXd neighbor_pairs(ArrayXXd coordinates, VectorXd species, double cutoff, int numAtoms) {
  auto padding_mask = species.array() == -1;
  int cols = 0;
  MatrixXd upper_triag(2, numAtoms * (numAtoms - 1) / 2);
  for (auto i = 0; i < numAtoms; i++) {
    for (auto j = i + 1; j < numAtoms; j++) {
      upper_triag(0, cols) = i;
      upper_triag(1, cols) = j;
      cols++;
    }
  }


  ArrayXd upper_triag_flattened(numAtoms * (numAtoms - 1));
  int index = 0;
  for (auto i = 0; i < upper_triag.rows(); i++) {
    for (auto j = 0; j < upper_triag.cols(); j++) {
      upper_triag_flattened(index) = upper_triag(i, j);
      index++;
    }

  }

  ArrayXXd pair_coordinates(numAtoms * (numAtoms - 1), 3);
  for (auto i = 0; i < upper_triag_flattened.size(); i++) {
    pair_coordinates.row(i) = coordinates.row(upper_triag_flattened(i));
  }

  ArrayXXd distances(numAtoms * (numAtoms - 1) / 2, 3);
  int num_pairs = numAtoms * (numAtoms - 1) / 2;
  for (auto i = 0; i < numAtoms * (numAtoms - 1) / 2; i++) {
    distances.row(i) = pair_coordinates.row(i) - pair_coordinates.row(i + num_pairs);
  }

  auto distance_norm = distances.matrix().rowwise().norm().array();

  ArrayXXd pair_padded_mask(numAtoms * (numAtoms - 1), 1);
  for (auto i = 0; i < upper_triag_flattened.size(); i++) {
    pair_padded_mask(i, 0) = padding_mask(upper_triag_flattened(i), 0);
  }

  ArrayXd dist(distance_norm.size());

  for (auto i = 0; i < distance_norm.size(); i++) {
    if (pair_padded_mask(i, 0) == 1) {
      dist(i) = INFINITY;
    }
    else {
      dist(i) = distance_norm(i);
    }
  }

  auto x = dist.array() <= cutoff;

  std::vector <int> indices;
  std::vector <std::pair <int, int>> atom_index12;

  for (auto i = 0; i < x.size(); i++) {
    if (x(i) == 1) {
      atom_index12.push_back(std::make_pair(upper_triag(0, i), upper_triag(1, i)));
    }
  }

  ArrayXd atom_index12_array(atom_index12.size() * 2);
  index = 0;
  for (auto i : atom_index12) {
    atom_index12_array(index) = i.first;
    index++;
  }
  for (auto i : atom_index12) {
    atom_index12_array(index) = i.second;
    index++;
  }


  return atom_index12_array;

}

ArrayXXd radial_terms(double cutoff, ArrayXXd distances) {
  auto fc = cosine_cutoff(distances, cutoff);
  ArrayXd ShfR(16);
  ShfR << 0.9,    1.16875,
          1.4375, 1.70625,
          1.975,  2.24375,
          2.5125, 2.78125,
          3.05,   3.31875,
          3.5875, 3.85625, 
          4.1250, 4.39375,
          4.6625, 4.93125;
  double EtaR = 16;
  ArrayXXd ret(distances.rows(), 16);
  for (int i = 0; i < distances.rows(); i++) {
    auto intermediate =  0.25 * ((ShfR - distances(i)).pow(2) * EtaR * -1).exp() * fc(i);
    ret.row(i) = intermediate;
  }
  return ret;
}

ArrayXXd SymmetryFunc(const ROMol &mol, std::vector<std::vector<double>> &res, int confId) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers");

  int numAtoms = mol.getNumAtoms();
  
  const auto conf = mol.getConformer(confId);

  ArrayXXd coordinates(numAtoms, 3);
  VectorXd species = generateSpeciesVector(mol);
  for (auto i = 0; i < numAtoms; i++) {
    auto pos = conf.getAtomPos(i);
    coordinates.row(i) << pos.x, pos.y, pos.z;
  }

  auto atom_index12 = neighbor_pairs(coordinates, species, 5.2, numAtoms);
  ArrayXXd selected_coordinates(atom_index12.rows(), 3);
  for (auto i = 0; i < atom_index12.rows(); i++) {
    selected_coordinates.row(i) = coordinates.row(atom_index12(i));
  }

  int pairs = selected_coordinates.rows() / 2;
  ArrayXXd vec(pairs, 3);
  for (auto i = 0; i < pairs; i++) {
    vec.row(i) = selected_coordinates.row(i) - selected_coordinates.row(i + pairs);
  }
  
  auto distances = vec.matrix().rowwise().norm().array();
  ArrayXXd species12(2, pairs);
  ArrayXXd species12_flipped(2, pairs);
  ArrayXXd atom_index12_unflattened(2, pairs);
  for (int i = 0; i < pairs; i++) {
    species12(0, i) = species(atom_index12(i));
    species12(1, i) = species(atom_index12(i + pairs));
    species12_flipped(1, i) = species(atom_index12(i));
    species12_flipped(0, i) = species(atom_index12(i + pairs));
    atom_index12_unflattened(0, i) = atom_index12(i);
    atom_index12_unflattened(1, i) = atom_index12(i + pairs);
  }
  auto index12 =  atom_index12_unflattened * 4 + species12_flipped;

  auto radial_terms_ = radial_terms(5.2, distances);

  ArrayXXd radial_aev = ArrayXXd::Zero(4 * numAtoms, 16);
  for (int i = 0; i < index12.cols(); i++) {
    for (int j = 0; j < 4 * numAtoms; j++) {
      if (index12(0, i) == j) {
        radial_aev.row(j) += radial_terms_.row(i);
      }
    }
  }
  for (int i = 0; i < index12.cols(); i++) {
    for (int j = 0; j < 4 * numAtoms; j++) {
      if (index12(1, i) == j) {
        radial_aev.row(j) += radial_terms_.row(i);
      }
    }
  }
  ArrayXXd final_radial_aev(numAtoms, 64);
  int atom_idx = 0;
  for (int i = 0; i < radial_aev.rows(); i+=4) {
    final_radial_aev.row(atom_idx) << radial_aev.row(i), radial_aev.row(i + 1), radial_aev.row(i + 2), radial_aev.row(i + 3);
    atom_idx++;
  }

  return final_radial_aev;

}
}
}