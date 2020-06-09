//
//  Copyright (C) 2020 Manan Goel
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <cmath>
#include "AtomicEnvironmentVector.h"
#include <Eigen/Dense>

using namespace Eigen;

namespace RDKit {
namespace Descriptors {
namespace ANI {

//! Calculates the value a continuous smoothening function for a distance such
//! that values
// greater than the cutoff give 0
/*!
  \param distances A 2 dimensional array of pairwise distances
  \param cutoff    A double value signifying cutoff distance

  \return 2 dimensional array containing corresponding values computed by cutoff
  function
*/
ArrayXXd CosineCutoff(ArrayXXd distances, double cutoff) {
  // Cosine cutoff function assuming all distances are less than the cutoff
  PRECONDITION((distances <= cutoff).count() == distances.size(),
               "All distances must be less than the cutoff");
  return 0.5 * (distances * (M_PI / cutoff)).cos() + 0.5;
}

VectorXi GenerateSpeciesVector(const ROMol &mol) {
  // Generate atom species vector as mentioned in torchani
  // H : 0
  // C : 1
  // N : 2
  // O : 3
  // All other atoms : -1
  auto numAtoms = mol.getNumAtoms();
  VectorXi species(numAtoms);

  for (auto i = 0; i < numAtoms; i++) {
    auto atom = mol.getAtomWithIdx(i);

    switch (atom->getAtomicNum()) {
      case 1:
        species[i] = 0;
        break;
      case 6:
        species[i] = 1;
        break;
      case 7:
        species[i] = 2;
        break;
      case 8:
        species[i] = 3;
        break;
      default:
        species[i] = -1;
    }
  }

  return species;
}

//! Constructs a vector with values of another vector at specified indices along
//! given dimension
/*!
  \param vector1    Matrix in which values are to be stored
  \param vector2    Matrix from which values are to be taken
  \param index      Array which specifies indices of vector2
  \param dim        dimension along which indices are to be picked

  \return Matrix containing values at positions specified by index in vector2
*/
ArrayXXd IndexSelect(ArrayXXd vector1, ArrayXXd vector2, ArrayXi index,
                     unsigned int dim) {
  PRECONDITION(dim == 0 || dim == 1,
               "Only values 0 and 1 are accepted for dim");
  for (auto i = 0; i < index.size(); i++) {
    switch (dim) {
      case 0:
        vector1.row(i) = vector2.row(index(i));
        break;
      case 1:
        vector1.col(i) = vector2.col(index(i));
        break;
      default:
        throw ValueErrorException("Value of dim must be 0 or 1");
    }
  }
  return vector1;
}

//! Computes pairs of atoms that are neighbors bypassing duplication to make
//! calculation faster
/*!
  \param coordinates  A matrix of size atoms * 3 containing coordinates of each
  atom \param species      A vector of size atoms containing mapping from atom
  index to encoding \param cutoff       Maximum distance within which 2 atoms
  are considered to be neighbours

  \return 2 dimensional array with 2 rows with each column corresponding to a
  pair of atoms which are neighbours
*/
ArrayXi NeighborPairs(ArrayXXd coordinates, VectorXi species, double cutoff,
                      unsigned int numAtoms) {
  PRECONDITION(coordinates.rows() == numAtoms,
               "Number of coordinate vectors must be same as number of atoms");
  PRECONDITION(species.size() == numAtoms,
               "Size of species vector must be same as number of atoms");

  // Find atoms which are not H, C, N or O
  auto paddingMask = species.array() == -1;
  unsigned int cols = 0;
  auto numPairs = numAtoms * (numAtoms - 1) / 2;

  // Each column contains pair of indices of the upper
  // traingular matrix of size numAtoms
  MatrixXi upperTriag(2, numPairs);
  for (auto i = 0; i < numAtoms; i++) {
    for (auto j = i + 1; j < numAtoms; j++) {
      upperTriag.col(cols) << i, j;
      cols++;
    }
  }

  // Flattened previously mentioned matrix with indices
  // Example :  (0, 0, 1) -> (0, 0, 1, 1, 2, 2)
  //            (1, 2, 2)
  ArrayXi upperTriagFlattened(numPairs * 2);
  unsigned int index = 0;
  for (auto i = 0; i < upperTriag.rows(); i++) {
    for (auto j = 0; j < upperTriag.cols(); j++) {
      upperTriagFlattened(index) = upperTriag(i, j);
      index++;
    }
  }

  // Select indices in upperTriagFlattened from coordinates
  ArrayXXd pairCoordinates(numPairs * 2, 3);
  pairCoordinates =
      IndexSelect(pairCoordinates, coordinates, upperTriagFlattened, 0);

  // Find vector and distances between each pair of atoms
  ArrayXXd distances(numPairs, 3);
  for (auto i = 0; i < numPairs; i++) {
    distances.row(i) =
        pairCoordinates.row(i) - pairCoordinates.row(i + numPairs);
  }

  auto distanceNorm = distances.matrix().rowwise().norm().array();

  // Create Mask for atoms which are not H, C, N or O
  // In the mask 0 if atom is H, C, N or O and 1 otherwise
  ArrayXXi pairPaddedMask(numPairs * 2, 1);
  for (auto i = 0; i < upperTriagFlattened.size(); i++) {
    pairPaddedMask(i, 0) = paddingMask(upperTriagFlattened(i), 0);
  }

  ArrayXd dist(distanceNorm.size());

  // Wherever mask is 1, set distance as infinity
  for (auto i = 0; i < distanceNorm.size(); i++) {
    if (pairPaddedMask(i, 0) == 1) {
      dist(i) = INFINITY;
    } else {
      dist(i) = distanceNorm(i);
    }
  }

  auto x = dist.array() <= cutoff;

  // Store atom pairs distance between which is less than cutoff
  std::vector<int> indices;
  std::vector<std::pair<int, int>> atomIndex12;

  for (auto i = 0; i < x.size(); i++) {
    if (x(i) == 1) {
      atomIndex12.push_back(std::make_pair(upperTriag(0, i), upperTriag(1, i)));
    }
  }

  ArrayXi atomIndex12Array(atomIndex12.size() * 2);
  index = 0;
  for (auto i : atomIndex12) {
    atomIndex12Array(index) = i.first;
    index++;
  }
  for (auto i : atomIndex12) {
    atomIndex12Array(index) = i.second;
    index++;
  }

  return atomIndex12Array;
}

//! Computes radial terms of the torchANI style atom features
/*!
  \param cutoff     Maximum distance between 2 atoms to show if they contribute
  to each other's environments \param distances  Distances between pairs of
  atoms which are in each other's neighbourhoods \return Radial terms according
  to each pair of distances in the molecule calculated using hard coded
  parameters
*/
ArrayXXd RadialTerms(double cutoff, ArrayXXd distances) {
  // Find cutoff factor for each pair of atoms which lie in each other's
  // neghborhoods
  // All the constants were determined in torchANI and have been taken from
  // there
  auto fc = CosineCutoff(distances, cutoff);
  // Different values for means of the gaussian symmetry functions
  ArrayXd ShfR(16);
  ShfR << 0.9, 1.16875, 1.4375, 1.70625, 1.975, 2.24375, 2.5125, 2.78125, 3.05,
      3.31875, 3.5875, 3.85625, 4.1250, 4.39375, 4.6625, 4.93125;
  // Variance terms for the gaussian symmetry functions
  double EtaR = 16;

  // Size of ret is
  // distances.rows() * (number of means * number of variance terms)
  ArrayXXd ret(distances.rows(), 16);
  for (auto i = 0; i < distances.rows(); i++) {
    auto intermediate =
        0.25 * ((ShfR - distances(i)).pow(2) * EtaR * -1).exp() * fc(i);
    ret.row(i) = intermediate;
  }
  return ret;
}

//! Computes angular terms of the torchANI style atom features
/*!
  \param cutoff     Maximum distance between 2 atoms to show if they contribute
  to each other's environments \param vector12   Pairs of vectors generated by a
  triplet of 3 atoms which are in each other's neighbourhoods \return Angular
  terms according to each pair of distances in the molecule calculated using
  hard coded parameters
*/
ArrayXXd AngularTerms(double cutoff, ArrayXXd vectors12) {
  // All the constants were determined in torchANI and have been taken from
  // there
  // Angle wise shift in the trigonometric term of the angular symmetry function
  ArrayXd ShfZ(8);
  ShfZ << 0.19634954, 0.58904862, 0.98174770, 1.3744468, 1.7671459, 2.1598449,
      2.5525440, 2.9452431;
  ArrayXd ShfA(4);
  // distance wise shifts in the distance term of the angular symmetry function
  ShfA << 0.90, 1.55, 2.20, 2.85;
  double zeta = 32;
  double etaA = 8;

  auto distances12 = vectors12.matrix().rowwise().norm().array();

  // Each triplet gives one angle formed by 2 vectors. The angle is found by
  // taking dot product of the 2 vectors
  ArrayXXd cosineAngles(vectors12.rows() / 2, 1);
  for (auto i = 0; i < vectors12.rows() / 2; i++) {
    cosineAngles(i, 0) = vectors12.matrix().row(i).dot(
        vectors12.matrix().row(i + vectors12.rows() / 2));
    cosineAngles(i, 0) =
        0.95 * cosineAngles(i, 0) /
        (vectors12.matrix().row(i).norm() *
         vectors12.matrix().row(i + vectors12.rows() / 2).norm());
  }
  auto angles = cosineAngles.acos();
  auto fcj12 = CosineCutoff(distances12, cutoff);

  // Angle dependent factors for the angular symmetry functions
  ArrayXXd factor1(angles.rows(), 8);
  for (auto i = 0; i < angles.rows(); i++) {
    factor1.row(i) =
        (((-1 * (ShfZ.rowwise() - angles.row(i)).array()).cos().array() + 1)
             .array() /
         2)
            .array()
            .pow(zeta);
  }

  // Distance dependent factors for the angular symmetry functions from both
  // vectors given by the triplet
  ArrayXXd distance12sum(distances12.rows() / 2, 1);

  for (auto i = 0; i < distances12.rows() / 2; i++) {
    distance12sum(i, 0) =
        distances12(i, 0) + distances12(i + distances12.rows() / 2, 0);
  }

  ArrayXXd factor2(distance12sum.rows(), 4);
  for (auto i = 0; i < distance12sum.rows(); i++) {
    factor2.row(i) =
        ((ShfA.rowwise() - (distance12sum.array() / 2).row(i)).pow(2).array() *
         -etaA)
            .exp();
  }

  // cutoff terms from both vectors of the triplet
  ArrayXXd fcj12prod(fcj12.rows() / 2, 1);
  for (auto i = 0; i < fcj12.rows() / 2; i++) {
    fcj12prod(i, 0) = fcj12(i, 0) * fcj12(i + fcj12.rows() / 2, 0);
  }

  ArrayXXd ret(fcj12prod.rows(), 32);
  for (auto i = 0; i < ret.rows(); i++) {
    unsigned int idx = 0;
    for (auto j = 0; j < 8; j++) {
      for (auto k = 0; k < 4; k++) {
        ret(i, idx) = 2 * factor1(i, j) * factor2(i, k) * fcj12prod(i, 0);
        idx++;
      }
    }
  }
  return ret;
}

std::vector<int> cumsum(std::vector<int> count, bool fromZero) {
  std::vector<int> cumulativeSum;
  if (fromZero == true) {
    cumulativeSum.push_back(0);
  } else {
    cumulativeSum.push_back(count[0]);
  }
  unsigned int index = 1;
  for (size_t i = 0 + (unsigned int)(!fromZero); i < count.size() - 1; i++) {
    cumulativeSum.push_back(cumulativeSum[index - 1] + count[i]);
    index++;
  }

  return cumulativeSum;
}

//! Calculates triplets of atoms according to pairs of atoms close to each other
/*!
  \param atomIndex12Angular   Pairs of atoms close to each other according to
  defined cutoff

  \return pair of a vector containing indices of centrals atoms and Matrix
  containing pairs of their neighbours

  \verbatim
    Input: indices for pairs of atoms that are close to each other.
    each pair only appear once, i.e. only one of the pairs (1, 2) and
    (2, 1) exists.

    Output: indices for all central atoms and it pairs of neighbors. For
    example, if input has pair (0, 1), (0, 2), (0, 3), (0, 4), (1, 2),
    (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), then the output would have
    central atom 0, 1, 2, 3, 4 and for cental atom 0, its pairs of neighbors
    are (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
  \endverbatim
*/
std::pair<std::vector<int>, ArrayXXi> TripleByMolecules(
    ArrayXXi atomIndex12Angular) {
  // atomIndex12Angular is flattened and sorted and we keep track of initial
  // indices so that actual pairs can be calculated back later
  std::vector<std::pair<int, int>> atomIndex12AngularFlattened;
  auto index = 0;
  for (auto i = 0; i < atomIndex12Angular.rows(); i++) {
    for (auto j = 0; j < atomIndex12Angular.cols(); j++) {
      atomIndex12AngularFlattened.push_back(
          std::make_pair(atomIndex12Angular(i, j), index));
      index++;
    }
  }
  std::stable_sort(atomIndex12AngularFlattened.begin(),
                   atomIndex12AngularFlattened.end());
  std::vector<int> revIndices, sortedAi;
  for (auto i : atomIndex12AngularFlattened) {
    revIndices.push_back(i.second);
    sortedAi.push_back(i.first);
  }

  // Find unique keys
  std::vector<int> uniqueResults(sortedAi.size());

  auto ip =
      std::unique_copy(sortedAi.begin(), sortedAi.end(), uniqueResults.begin());
  uniqueResults.resize(std::distance(uniqueResults.begin(), ip));

  // Number of occurrences of each unique key
  std::vector<int> counts(uniqueResults.size()),
      pairSizes(uniqueResults.size());
  for (size_t i = 0; i < uniqueResults.size(); i++) {
    counts[i] = std::count(sortedAi.begin(), sortedAi.end(), uniqueResults[i]);
    pairSizes[i] = counts[i] * (counts[i] - 1) / 2;
  }
  // Compute vector of central atom indices
  std::vector<int> pairIndices;

  for (size_t i = 0; i < pairSizes.size(); i++) {
    auto j = pairSizes[i];
    while (j--) {
      pairIndices.push_back(i);
    }
  }

  std::vector<int> centralAtomIndex;
  for (size_t i = 0; i < pairIndices.size(); i++) {
    centralAtomIndex.push_back(uniqueResults[pairIndices[i]]);
  }

  // do local combinations within unique key, assuming sorted
  int m;
  if (counts.size() > 0) {
    m = *std::max_element(counts.begin(), counts.end());
  } else {
    m = 0;
  }

  auto n = pairSizes.size();
  ArrayXXi lowerTriang(2, m * (m - 1) / 2);
  ArrayXXi intraPairIndices(2 * n, m * (m - 1) / 2);
  index = 0;

  for (auto i = 1; i < m; i++) {
    for (auto j = 0; j < i; j++) {
      lowerTriang.col(index) << i, j;
      index++;
    }
  }
  index = 0;
  for (auto i = 0; i < lowerTriang.rows(); i++) {
    auto j = n;
    while (j--) {
      intraPairIndices.row(index) << lowerTriang.row(i);
      index++;
    }
  }

  ArrayXXi mask(1, pairSizes.size() * lowerTriang.cols());
  index = 0;
  for (size_t i = 0; i < pairSizes.size(); i++) {
    for (auto j = 0; j < lowerTriang.cols(); j++) {
      mask(0, index) = pairSizes[i] > j;
      index++;
    }
  }

  ArrayXXi intraPairIndicesFlattened(2, n * m * (m - 1) / 2);
  if (intraPairIndices.rows() != 0 && intraPairIndices.cols() != 0) {
    index = 0;
    for (auto i = 0; i < intraPairIndices.rows() / 2; i++) {
      for (auto j = 0; j < m * (m - 1) / 2; j++) {
        intraPairIndicesFlattened(0, index + j) = intraPairIndices(i, j);
      }
      index += (m * (m - 1) / 2);
    }
    index = 0;

    for (auto i = intraPairIndices.rows() / 2; i < intraPairIndices.rows();
         i++) {
      for (auto j = 0; j < m * (m - 1) / 2; j++) {
        intraPairIndicesFlattened(1, index + j) = intraPairIndices(i, j);
      }
      index += (m * (m - 1) / 2);
    }
  }
  ArrayXXi sortedLocalIndex12(2, (mask.array() == 1).count());
  index = 0;
  for (auto i = 0; i < mask.size(); i++) {
    if (mask(0, i) == 1) {
      sortedLocalIndex12.col(index) << intraPairIndicesFlattened.col(i);
      index++;
    }
  }

  auto cumsumCount = cumsum(counts, true);

  VectorXi extraLocalIndex12(pairIndices.size());
  index = 0;
  for (auto i : pairIndices) {
    extraLocalIndex12(index) = cumsumCount[i];
    index++;
  }

  sortedLocalIndex12 =
      (sortedLocalIndex12.matrix().rowwise() + extraLocalIndex12.transpose())
          .array();

  // unsort from last part
  ArrayXXi localIndex12(2, sortedLocalIndex12.cols());
  for (auto j = 0; j < sortedLocalIndex12.cols(); j++) {
    localIndex12(0, j) = revIndices[sortedLocalIndex12(0, j)];
    localIndex12(1, j) = revIndices[sortedLocalIndex12(1, j)];
  }

  return std::make_pair(centralAtomIndex, localIndex12);
}

ArrayXXi TriuIndex(unsigned int numSpecies) {
  std::vector<int> species1, species2, pairIndex;

  for (auto i = 0; i < numSpecies; i++) {
    for (auto j = i; j < numSpecies; j++) {
      species1.push_back(i);
      species2.push_back(j);
    }
    pairIndex.push_back(i);
  }
  ArrayXXi ret(numSpecies, numSpecies);
  unsigned int index1 = 0;

  for (size_t i = 0; i < species1.size(); i++) {
    ret(species1[i], species2[i]) = index1;
    ret(species2[i], species1[i]) = index1;
    index1++;
  }

  return ret;
}

//! Accumulate the elements of a Matrix into another Matrix by adding to the
//! indices in the order given in index.
/*!
  \param vector1    Matrix to which values are to be added
  \param vector2    Matrix from which values are to be added
  \param index      Indices in order of which values are added
  \param multi      Number of pairs to be considered
  \param numAtoms   Number of atoms in the molecules

  \return Matrix containing accumulated elements of vector2 into vector1
  according to order given in index

  \verbatim
    Index[i] == j, then the ith row of vector2 is added to the jth row of
  vector1
  \endverbatim
*/
ArrayXXd IndexAdd(ArrayXXd vector1, ArrayXXd vector2, ArrayXXi index,
                  unsigned int multi, unsigned int numAtoms) {
  PRECONDITION(vector1.rows() == multi * numAtoms, "Too few rows in vector1");
  PRECONDITION(vector2.rows() == index.rows(), "Too few rows in vector2");
  for (auto idxCol = 0; idxCol < index.cols(); idxCol++) {
    for (auto i = 0; i < index.rows(); i++) {
      for (auto j = 0; j < multi * numAtoms; j++) {
        if (index(i, idxCol) == j) {
          vector1.row(j) += vector2.row(i);
        }
      }
    }
  }
  return vector1;
}

ArrayXXd AtomicEnvironmentVector(const ROMol &mol, int confId) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers");

  auto numAtoms = mol.getNumAtoms();

  const auto conf = mol.getConformer(confId);

  ArrayXXd coordinates(numAtoms, 3);
  auto species = GenerateSpeciesVector(mol);
  for (auto i = 0; i < numAtoms; i++) {
    auto pos = conf.getAtomPos(i);
    coordinates.row(i) << pos.x, pos.y, pos.z;
  }
  // Fetch pairs of atoms which are neigbours which lie within the cutoff
  // distance 5.2 Angstroms. The constant was obtained by authors of torchANI
  auto atomIndex12 = NeighborPairs(coordinates, species, 5.2, numAtoms);
  ArrayXXd selectedCoordinates(atomIndex12.rows(), 3);
  selectedCoordinates =
      IndexSelect(selectedCoordinates, coordinates, atomIndex12, 0);

  // Vectors between pairs of atoms that lie in each other's neighborhoods
  unsigned int pairs = selectedCoordinates.rows() / 2;
  ArrayXXd vec(pairs, 3);
  for (auto i = 0; i < pairs; i++) {
    vec.row(i) =
        selectedCoordinates.row(i) - selectedCoordinates.row(i + pairs);
  }

  auto distances = vec.matrix().rowwise().norm().array();
  ArrayXXi species12(2, pairs);
  ArrayXXi species12Flipped(2, pairs);
  ArrayXXi atomIndex12Unflattened(2, pairs);
  for (auto i = 0; i < pairs; i++) {
    species12(0, i) = species(atomIndex12(i));
    species12(1, i) = species(atomIndex12(i + pairs));

    species12Flipped(1, i) = species(atomIndex12(i));
    species12Flipped(0, i) = species(atomIndex12(i + pairs));

    atomIndex12Unflattened(0, i) = atomIndex12(i);
    atomIndex12Unflattened(1, i) = atomIndex12(i + pairs);
  }

  // Obtain indices to use for constructing final radial aev
  // The constant 4 comes from the fact that each atom will have 4 kinds of
  // interactions. For Example, for H it is HC, HO, HN and HX where X is any
  // other atom
  auto index12 = atomIndex12Unflattened * 4 + species12Flipped;

  auto RadialTerms_ = RadialTerms(5.2, distances);

  ArrayXXd radialAEV = ArrayXXd::Zero(4 * numAtoms, 16);

  radialAEV =
      IndexAdd(radialAEV, RadialTerms_, index12.transpose(), 4, numAtoms);

  // Each atom finally has a total of 64 radial terms
  ArrayXXd finalRadialAEV(numAtoms, 64);
  unsigned int atomIdx = 0;
  for (auto i = 0; i < radialAEV.rows(); i += 4) {
    finalRadialAEV.row(atomIdx) << radialAEV.row(i), radialAEV.row(i + 1),
        radialAEV.row(i + 2), radialAEV.row(i + 3);
    atomIdx++;
  }

  // Distance cutoff for angular terms is smaller  than for radial terms, so we
  // construct a smaller neighbor list. The authors of torchANI found that the
  // cutoff 3.5 gave best results
  ArrayXi evenCloserIndices((distances.array() <= 3.5).count());
  unsigned int idx = 0;
  for (auto i = 0; i < distances.size(); i++) {
    if (distances(i) <= 3.5) {
      evenCloserIndices(idx) = i;
      idx++;
    }
  }

  // Angular Terms
  // Construct all the previously mentioned neighbor pair information for the
  // smaller neighbor list
  ArrayXXi species12Angular(2, evenCloserIndices.size());
  ArrayXXi atomIndex12Angular(2, evenCloserIndices.size());

  ArrayXXd vecAngular(evenCloserIndices.size(), 3);

  species12Angular =
      IndexSelect(species12Angular.matrix().cast<double>().array(),
                  species12.matrix().cast<double>().array(), evenCloserIndices,
                  1)
          .matrix()
          .cast<int>()
          .array();
  atomIndex12Angular =
      IndexSelect(atomIndex12Angular.matrix().cast<double>().array(),
                  atomIndex12Unflattened.matrix().cast<double>().array(),
                  evenCloserIndices, 1)
          .matrix()
          .cast<int>()
          .array();
  vecAngular = IndexSelect(vecAngular, vec, evenCloserIndices, 0);

  auto n = evenCloserIndices.size();

  // Find Triplets for which angular terms are to be found
  // TripleByMolecules returns an array of central atoms and the corresponding
  // pair of atoms in an STL pair
  auto tripletInfo = TripleByMolecules(atomIndex12Angular);
  auto pairIndex12 = tripletInfo.second;
  auto centralAtomIndex = tripletInfo.first;
  ArrayXXi sign12(2, pairIndex12.cols());

  // compute mapping between representation of central-other to pair
  for (auto i = 0; i < pairIndex12.rows(); i++) {
    for (auto j = 0; j < pairIndex12.cols(); j++) {
      if (pairIndex12(i, j) < n) {
        sign12(i, j) = 1;
      } else {
        sign12(i, j) = -1;
      }
    }
  }

  n = atomIndex12Angular.cols();

  // pairIndex12 = pairindex12 % n
  auto localIndex = pairIndex12.cast<int>();

  pairIndex12 =
      (localIndex.array() - (localIndex.array() / n).array() * n).array();

  ArrayXi pairIndex12Flattened(2 * pairIndex12.cols());
  idx = 0;
  for (auto i = 0; i < pairIndex12.rows(); i++) {
    for (auto j = 0; j < pairIndex12.cols(); j++) {
      pairIndex12Flattened(idx) = pairIndex12(i, j);
      idx++;
    }
  }

  ArrayXXd vecFlattened(pairIndex12Flattened.size(), 3);
  vecFlattened = IndexSelect(vecFlattened, vecAngular, pairIndex12Flattened, 0);

  ArrayXXd vec12(vecFlattened.rows(), 3);
  for (auto i = 0; i < vecFlattened.rows() / 2; i++) {
    vec12.row(i) = vecFlattened.row(i) * sign12(0, i);
  }

  for (auto i = vecFlattened.rows() / 2; i < vecFlattened.rows(); i++) {
    vec12.row(i) = vecFlattened.row(i) * sign12(1, i - vecFlattened.rows() / 2);
  }

  auto AngularTerms_ = AngularTerms(3.5, vec12);

  ArrayXXi centralAtomIndexArr(centralAtomIndex.size(), 1);

  for (size_t i = 0; i < centralAtomIndex.size(); i++) {
    centralAtomIndexArr.row(i) << centralAtomIndex[i];
  }

  ArrayXXi species12Small1(2, pairIndex12.cols());
  ArrayXXi species12Small2(2, pairIndex12.cols());

  for (auto i = 0; i < pairIndex12.rows(); i++) {
    for (auto j = 0; j < pairIndex12.cols(); j++) {
      species12Small1(i, j) = species12Angular(0, pairIndex12(i, j));
    }
  }

  for (auto i = 0; i < pairIndex12.rows(); i++) {
    for (auto j = 0; j < pairIndex12.cols(); j++) {
      species12Small2(i, j) = species12Angular(1, pairIndex12(i, j));
    }
  }

  ArrayXXi species12_(sign12.rows(), sign12.cols());

  for (auto i = 0; i < sign12.rows(); i++) {
    for (auto j = 0; j < sign12.cols(); j++) {
      if (sign12(i, j) == 1) {
        species12_(i, j) = species12Small2(i, j);
      } else {
        species12_(i, j) = species12Small1(i, j);
      }
    }
  }

  ArrayXXi index(species12_.cols(), 1);
  auto triuIndices = TriuIndex(4);

  for (auto i = 0; i < species12_.cols(); i++) {
    index.row(i) = triuIndices(species12_(0, i), species12_(1, i));
  }
  // The constant 10 comes from 10 pairs that can be formed
  index = index + (centralAtomIndexArr.array() * 10).array();

  ArrayXXd angularAEV = ArrayXXd::Zero(10 * numAtoms, 32);
  angularAEV = IndexAdd(angularAEV, AngularTerms_, index, 10, numAtoms);

  ArrayXXd finalAngularAEV(numAtoms, 320);
  atomIdx = 0;
  for (auto i = 0; i < angularAEV.rows(); i += 10) {
    finalAngularAEV.row(atomIdx) << angularAEV.row(i), angularAEV.row(i + 1),
        angularAEV.row(i + 2), angularAEV.row(i + 3);
    angularAEV.row(i + 4), angularAEV.row(i + 5);
    angularAEV.row(i + 6), angularAEV.row(i + 7);
    angularAEV.row(i + 8), angularAEV.row(i + 9);
    atomIdx++;
  }

  ArrayXXd finalAEV(finalRadialAEV.rows(),
                    finalRadialAEV.cols() + finalAngularAEV.cols());
  finalAEV << finalRadialAEV, finalAngularAEV;

  return finalAEV;
}
}  // namespace ANI
}  // namespace Descriptors
}  // namespace RDKit
