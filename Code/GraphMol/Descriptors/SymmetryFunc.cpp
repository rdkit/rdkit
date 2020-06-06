//
//  Copyright (c) 2020, Manan Goel
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//

#include <GraphMol/RDKitBase.h>
#include <cmath>
#include "SymmetryFunc.h"
#include <Eigen/Dense>

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

ArrayXXd index_select(ArrayXXd vector1, ArrayXXd vector2, ArrayXd index, int dim) {
  for (auto i = 0; i < index.size(); i++) {
    if (dim == 0) {
      vector1.row(i) = vector2.row(index(i));
    }
    else if (dim == 1) {
      vector1.col(i) = vector2.col(index(i));
    }
  }
  return vector1;
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
  pair_coordinates = index_select(pair_coordinates, coordinates, upper_triag_flattened, 0);

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

ArrayXXd angular_terms(double cutoff, ArrayXXd vectors12) {


  ArrayXd ShfZ(8);
  ShfZ << 0.19634954, 0.58904862,
          0.98174770, 1.3744468,
          1.7671459 , 2.1598449,
          2.5525440 , 2.9452431;
  ArrayXd ShfA(4);
  ShfA << 0.90, 1.55, 2.20, 2.85;
  double zeta = 32;
  double etaA = 8;

  auto distances12 = vectors12.matrix().rowwise().norm().array();

  ArrayXXd cosine_angles(vectors12.rows() / 2, 1);
  for (int i = 0; i < vectors12.rows() / 2; i++) {
    cosine_angles(i, 0) = vectors12.matrix().row(i).dot(vectors12.matrix().row(i + vectors12.rows() / 2));
    cosine_angles(i, 0) = 0.95 * cosine_angles(i, 0) / (vectors12.matrix().row(i).norm() * vectors12.matrix().row(i + vectors12.rows() / 2).norm()); 
  }
  auto angles = cosine_angles.acos();
  auto fcj12 = cosine_cutoff(distances12, cutoff);

  ArrayXXd factor1(angles.rows(), 8);
  for (int i = 0; i < angles.rows(); i++) {
    factor1.row(i) =  (((-1 * (ShfZ.rowwise() - angles.row(i)).array()).cos().array() + 1).array() / 2).array().pow(zeta);
  }

  ArrayXXd distance12sum(distances12.rows() / 2, 1);


  for (int i = 0; i < distances12.rows() / 2; i++) {
    distance12sum(i, 0) = distances12(i, 0) + distances12(i + distances12.rows() / 2, 0);
  }

  ArrayXXd factor2(distance12sum.rows(), 4);
  for (int i = 0; i < distance12sum.rows(); i++) {
    factor2.row(i) = ((ShfA.rowwise() - (distance12sum.array() / 2).row(i)).pow(2).array() * -etaA).exp();
  }

  ArrayXXd fcj12prod(fcj12.rows() / 2, 1);
  for (int i = 0; i < fcj12.rows() / 2; i++) {
    fcj12prod(i, 0) = fcj12(i, 0) * fcj12(i + fcj12.rows() / 2, 0);
  }

  ArrayXXd ret(fcj12prod.rows(), 32);
  for (int i = 0; i < ret.rows(); i++) {
    int idx = 0;
    for (int j = 0; j < 8; j++) {
      for (int k = 0; k < 4; k++) {
        ret(i, idx) = 2 * factor1(i, j) * factor2(i, k) * fcj12prod(i, 0);
        idx++;
      }
    }
  }
  return ret;

}

std::vector <int> cumsum_from_zero(std::vector <int> count) {
  std::vector <int> cumsum{0};
  int index = 1;
  for (int i = 0; i < count.size() - 1; i++)  {
    cumsum.push_back(cumsum[index - 1] + count[i]);
    index++;
  }

  return cumsum;
}

std::pair <std::vector<int>, ArrayXXd> triple_by_molecules(ArrayXXd atom_index12_angular) {
  std::vector <std::pair<int, int>> atom_index12_angular_flattened;
  auto index = 0;
  for (int i = 0; i < atom_index12_angular.rows(); i++) {
    for (int j = 0; j < atom_index12_angular.cols(); j++) {
      atom_index12_angular_flattened.push_back(std::make_pair(atom_index12_angular(i, j), index));
      index++;
    }
  }
  std::sort(atom_index12_angular_flattened.begin(), atom_index12_angular_flattened.end());
  std::vector <int> rev_indices, sorted_ai;
  for (auto i : atom_index12_angular_flattened) {
    rev_indices.push_back(i.second);
    sorted_ai.push_back(i.first);
  }

  

  std::vector<int> unique_results(sorted_ai.size());
  
  auto ip = std::unique_copy(sorted_ai.begin(), sorted_ai.end(), unique_results.begin());
  unique_results.resize(std::distance(unique_results.begin(), ip));


  std::vector<int> counts(unique_results.size()), pair_sizes(unique_results.size());
  for (int i = 0; i < unique_results.size(); i++) {
    counts[i] = std::count(sorted_ai.begin(), sorted_ai.end(), unique_results[i]);
    pair_sizes[i] = counts[i] * (counts[i] - 1) / 2;
  }

  std::vector <int> pair_indices;

  for (int i = 0; i < pair_sizes.size(); i++) {
    int j = pair_sizes[i];
    while (j--) {
      pair_indices.push_back(i);
    }
  }

  std::vector <int> central_atom_index;
  for (int i = 0; i < pair_indices.size(); i++) {
    central_atom_index.push_back(unique_results[pair_indices[i]]);
  }

  int m;
  if (counts.size() > 0) {
    m = *std::max_element(counts.begin(), counts.end());
  }
  else {
    m = 0;
  }

  auto n = pair_sizes.size();
  ArrayXXd lower_triang(2, m * (m - 1) / 2);
  ArrayXXd intra_pair_indices(2 * n, m * (m - 1) / 2);
  index = 0;

  for (int i = 1; i < m; i++) {
    for (int j = 0; j < i; j++) {
      lower_triang.col(index) << i, j;
      index++;
    }
  }
  index = 0;
  for (int i = 0; i < lower_triang.rows(); i++) {
    int j = n;
    while (j--) {
      intra_pair_indices.row(index) << lower_triang.row(i);
      index++;
    }
  }

  ArrayXXd mask(1, pair_sizes.size() * lower_triang.cols());
  index = 0;
  for (int i = 0; i < pair_sizes.size(); i++) {
    for (int j = 0; j < lower_triang.cols(); j++) {
      mask(0, index) =  pair_sizes[i] > j;
      index++;
    }
  }

  ArrayXXd intra_pair_indices_flattened(2, n * m * (m - 1) / 2);
  if (intra_pair_indices.rows() != 0 && intra_pair_indices.cols() != 0) {
    index = 0;
    for (int i = 0; i < intra_pair_indices.rows() / 2; i++) {
      for (int j = 0; j < m * (m - 1) / 2; j++) {
        intra_pair_indices_flattened(0, index + j) = intra_pair_indices(i, j);
      }
      index += (m * (m - 1) / 2);
    }
    index = 0;

    for (int i = intra_pair_indices.rows() / 2; i < intra_pair_indices.rows(); i++) {
      for (int j = 0; j < m * (m - 1) / 2; j++) {
        intra_pair_indices_flattened(1, index + j) = intra_pair_indices(i, j);
      }
      index += (m * (m - 1) / 2);
    }
  }
  ArrayXXd sorted_local_index12(2, (mask.array() == 1).count());
  index = 0;
  for (int i = 0; i < mask.size(); i++) {
    if (mask(0, i) == 1) {
      sorted_local_index12.col(index) << intra_pair_indices_flattened.col(i);
      index++;
    }
  }

  auto cumsum_count = cumsum_from_zero(counts);

  VectorXd extra_local_index12 (pair_indices.size());
  index = 0;
  for (auto i : pair_indices) {
    extra_local_index12(index) = cumsum_count[i];
    index++;
  }

  sorted_local_index12 = (sorted_local_index12.matrix().rowwise() + extra_local_index12.transpose()).array();

  ArrayXXd local_index12(2, sorted_local_index12.cols());
  for (int j = 0; j < sorted_local_index12.cols(); j++) {
    local_index12(0, j) = rev_indices[sorted_local_index12(0, j)];
    local_index12(1, j) = rev_indices[sorted_local_index12(1, j)];
  }

  return std::make_pair(central_atom_index, local_index12);

}

ArrayXXd triu_index(int num_species) {
  std::vector <int> species1, species2, pair_index;

  for (int i = 0; i < num_species; i++) {
    for (int j = i; j < num_species; j++) {
      species1.push_back(i);
      species2.push_back(j);
    }
    pair_index.push_back(i);
  }
  ArrayXXd ret(num_species, num_species);
  int index1 = 0, index2 = pair_index.size() - 1;
  
  for (int i = 0; i < species1.size(); i++) {
    ret(species1[i], species2[i]) = index1;
    ret(species2[i], species1[i]) = index1;
    index1++;
  }

  return ret;
}

ArrayXXd index_add(ArrayXXd vector1, ArrayXXd vector2, ArrayXXd index, int multi, int numAtoms) {
  for (int idx_col = 0; idx_col < index.cols(); idx_col++) {
    for (int i = 0; i < index.rows(); i++) {
      for (int j = 0; j < multi * numAtoms; j++) {
        if (index(i, idx_col) == j) {
          vector1.row(j) += vector2.row(i);
        }
      }
    }
  }
  return vector1;
}

ArrayXXd SymmetryFunc(const ROMol &mol, int confId) {
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
  selected_coordinates = index_select(selected_coordinates, coordinates, atom_index12, 0);

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

  radial_aev = index_add(radial_aev, radial_terms_, index12.transpose(), 4, numAtoms);

  ArrayXXd final_radial_aev(numAtoms, 64);
  int atom_idx = 0;
  for (int i = 0; i < radial_aev.rows(); i+=4) {
    final_radial_aev.row(atom_idx) << radial_aev.row(i), radial_aev.row(i + 1), radial_aev.row(i + 2), radial_aev.row(i + 3);
    atom_idx++;
  }

  ArrayXd even_closer_indices((distances.array() <= 3.5).count());
  int idx = 0;
  for (int i = 0; i < distances.size(); i++) {
    if (distances(i) <= 3.5) {
      even_closer_indices(idx) = i;
      idx++;
    }
  }

  // Angular Terms

  ArrayXXd species12_angular(2, even_closer_indices.size());
  ArrayXXd atom_index12_angular(2, even_closer_indices.size());
  
  ArrayXXd vec_angular(even_closer_indices.size(), 3);

  species12_angular = index_select(species12_angular, species12, even_closer_indices, 1);
  atom_index12_angular = index_select(atom_index12_angular, atom_index12_unflattened, even_closer_indices, 1);
  vec_angular = index_select(vec_angular, vec, even_closer_indices, 0);


  auto n = even_closer_indices.size();
  auto ret = triple_by_molecules(atom_index12_angular);
  auto pair_index12 = ret.second;
  auto central_atom_index = ret.first;
  ArrayXXd sign12(2, pair_index12.cols());

  for (int i = 0; i < pair_index12.rows(); i++) {
    for (int j = 0; j < pair_index12.cols(); j++) {
      if (pair_index12(i, j) < n) {
        sign12(i, j) = 1;
      }
      else {
        sign12(i, j) = -1;
      }
    }
  }

  n = atom_index12_angular.cols();

  auto local_index = pair_index12.cast <int>();
  pair_index12 = (local_index.array() - (local_index.array() / n).array() * n).array().cast <double>();

  ArrayXd pair_index12_flattened(2 * pair_index12.cols());
  idx = 0;
  for (int i = 0; i < pair_index12.rows(); i++) {
    for (int j = 0; j < pair_index12.cols(); j++) {
      pair_index12_flattened(idx) = pair_index12(i, j);
      idx++;
    }
  }

  ArrayXXd vec_flattened(pair_index12_flattened.size(), 3);
  vec_flattened = index_select(vec_flattened, vec_angular, pair_index12_flattened, 0);

  ArrayXXd vec12(vec_flattened.rows(), 3);
  for (int i = 0; i < vec_flattened.rows() / 2; i++) {
    vec12.row(i) = vec_flattened.row(i) * sign12(0, i);
  }

  for (int i = vec_flattened.rows() / 2; i < vec_flattened.rows(); i++) {
    vec12.row(i) = vec_flattened.row(i) * sign12(1, i - vec_flattened.rows() / 2);
  }

  auto angular_terms_ = angular_terms(3.5, vec12);

  ArrayXXd central_atom_index_arr(central_atom_index.size(), 1);

  for (int i = 0; i < central_atom_index.size(); i++) {
    central_atom_index_arr.row(i) << central_atom_index[i];
  }

  ArrayXXd species12_small_1(2, pair_index12.cols());
  ArrayXXd species12_small_2(2, pair_index12.cols());

  for (int i = 0; i < pair_index12.rows(); i++) {
    for (int j = 0; j < pair_index12.cols(); j++) {
      species12_small_1(i, j) = species12_angular(0, pair_index12(i, j));
    }
  }

  for (int i = 0; i < pair_index12.rows(); i++) {
    for (int j = 0; j < pair_index12.cols(); j++) {
      species12_small_2(i, j) = species12_angular(1, pair_index12(i, j));
    }
  }

  ArrayXXd species12_(sign12.rows(), sign12.cols());

  for (int i = 0; i < sign12.rows(); i++) {
    for (int j = 0; j < sign12.cols(); j++) {
      if (sign12(i, j) == 1) {
        species12_(i, j) = species12_small_2(i, j);
      }
      else {
        species12_(i, j) = species12_small_1(i, j);
      }
    }
  }

  ArrayXXd index(species12_.cols(), 1);
  auto triu_indices = triu_index(4);

  for (int i = 0; i < species12_.cols(); i++) {
    index.row(i) = triu_indices(species12_(0, i), species12_(1, i));
  }

  index = index + (central_atom_index_arr.array() * 10).array();

  ArrayXXd angular_aev = ArrayXXd::Zero(10 * numAtoms, 32);
  angular_aev = index_add(angular_aev, angular_terms_, index, 10, numAtoms);

  ArrayXXd final_angular_aev(numAtoms, 320);
  atom_idx = 0;
  for (int i = 0; i < angular_aev.rows(); i+=10) {
    final_angular_aev.row(atom_idx) << angular_aev.row(i),     angular_aev.row(i + 1), 
                                       angular_aev.row(i + 2), angular_aev.row(i + 3);
                                       angular_aev.row(i + 4), angular_aev.row(i + 5);
                                       angular_aev.row(i + 6), angular_aev.row(i + 7);
                                       angular_aev.row(i + 8), angular_aev.row(i + 9);
    atom_idx++;
  }
  
  ArrayXXd final_aev(final_radial_aev.rows(), final_radial_aev.cols() + final_angular_aev.cols());
  final_aev << final_radial_aev,
               final_angular_aev;

  return final_aev;

}
}
}