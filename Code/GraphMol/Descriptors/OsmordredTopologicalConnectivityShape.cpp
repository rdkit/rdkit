//  Copyright (c) 2025, Guillaume Godin Osmo Labs, PBC’s and others
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
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
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

#include "Osmordred.h"
#include "OsmordredHelpers.h"
#include <boost/functional/hash.hpp>  // For custom hashing of pairs
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <GraphMol/Descriptors/BCUT.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Atom.h>
#include <RDGeneral/export.h>
#include <RDGeneral/types.h>

#include <boost/graph/adjacency_list.hpp>

#include <set>
#include <cmath>  // For M_PI and pow
#include <tuple>
#include <map>
#include <string>
#include <utility>  // for std::pair
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <climits>
#include <queue>  // for fused rings
#include <stdexcept>
#include <iomanip>  // For std::fixed and std::setprecision
#include <sstream>  // For std::ostringstream
#include <iostream>
#include <cstring>  // For memcpy
#include <functional>
#include <numeric>
#include <stack>

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <GraphMol/RingInfo.h>
#include <RDGeneral/RDThreads.h>
#include <future>
#include <chrono>

// Define a custom hash function for std::pair<int, int>
namespace std {
template <>
struct hash<std::pair<int, int>> {
  size_t operator()(const std::pair<int, int> &p) const {
    // Combine the hash of the two elements of the pair
    return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
  }
};
}  // namespace std

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

std::vector<double> calcTopologicalChargeDescs(const ROMol &mol) {
  const int maxOrder = 10;
  std::vector<double> results(21, 0.0);

  // Calculate Distance and Adjacency Matrices
  Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(mol);
  Eigen::MatrixXd adjacencyMatrix = calculateAdjacencyMatrix(mol);

  // Step 1: Compute the Charge Term Matrix
  Eigen::MatrixXd CT =
      calculateChargeTermMatrix(adjacencyMatrix, distanceMatrix);

  // Step 3: Create a lower triangular distance matrix
  Eigen::MatrixXd lowerTriangularMask =
      Eigen::MatrixXd::Zero(distanceMatrix.rows(), distanceMatrix.cols());
  for (int i = 0; i < distanceMatrix.rows(); ++i) {
    for (int j = 0; j < i; ++j) {
      lowerTriangularMask(i, j) = 1.0;
    }
  }
  Eigen::MatrixXd D = distanceMatrix.cwiseProduct(lowerTriangularMask);
  D = D.unaryExpr([](double x) {
    return (x == 0) ? std::numeric_limits<double>::infinity() : x;
  });

  // Step 3: Calculate raw, mean, and global descriptors
  for (int order = 1; order <= maxOrder; ++order) {
    // Create mask for the current order
    Eigen::MatrixXd orderMask = (D.array() == order).cast<double>();

    // Apply mask to CT / D
    std::vector<double> filteredCT, filteredD;

    for (int i = 0; i < CT.rows(); ++i) {
      for (int j = 0; j < CT.cols(); ++j) {
        if (orderMask(i, j) > 0) {
          filteredCT.push_back(CT(i, j));
          filteredD.push_back(D(i, j));
        }
      }
    }

    // Raw descriptor: Absolute sum of filtered CT values
    double raw = 0.0;
    for (double val : filteredCT) {
      raw += std::abs(val);
    }
    results[order - 1] = raw;

    // Mean ie Normalize by frequencies
    if (!filteredCT.empty()) {
      // dict/Map to store frequencies of each distance value
      std::unordered_map<double, int> frequencies;
      for (double val : filteredD) {
        frequencies[val]++;
      }

      // get mean descriptor
      double mean = 0.0;
      for (size_t i = 0; i < filteredCT.size(); ++i) {
        mean += std::abs(filteredCT[i]) / frequencies[filteredD[i]];
      }
      results[10 + order - 1] = mean;
    }
  }

  // Global descriptor: full Sum of all mean values from 1 to maxOrder

  Eigen::MatrixXd orderMask = (D.array() <= maxOrder).cast<double>();

  // Apply mask to CT
  std::vector<double> gfilteredD, gfilteredCT;
  for (int i = 0; i < CT.rows(); ++i) {
    for (int j = 0; j < CT.cols(); ++j) {
      if (orderMask(i, j) > 0) {
        gfilteredCT.push_back(CT(i, j));
        gfilteredD.push_back(D(i, j));
      }
    }
  }

  // Mean descriptor: Normalize by frequencies
  if (!gfilteredCT.empty()) {
    // Map to store frequencies of each distance value
    std::unordered_map<double, int> gfrequencies;
    for (double val : gfilteredD) {
      gfrequencies[val]++;
    }
    double global = 0.0;
    // Calculate mean descriptor
    for (size_t i = 0; i < gfilteredCT.size(); ++i) {
      global += std::abs(gfilteredCT[i]) / gfrequencies[gfilteredD[i]];
    }
    results[20] = global;
  }

  return results;
}

/// not tested TopologocalIndex (4 outputs => need eccentricities computation)

// Function to calculate the radius (minimum of eccentricities)
double calculateRadius(const Eigen::VectorXd &eccentricities) {
  if (eccentricities.size() == 0) {
    throw std::invalid_argument("Eccentricity data is empty");
  }
  return eccentricities.minCoeff();  // Eigen function to get the minimum value
}

// Function to calculate the diameter (maximum of eccentricities)
double calculateDiameter(const Eigen::VectorXd &eccentricities) {
  if (eccentricities.size() == 0) {
    throw std::invalid_argument("Eccentricity data is empty");
  }
  return eccentricities.maxCoeff();  // Eigen function to get the maximum value
}

// Function to calculate the Topological Shape Index
double calculateTopologicalShapeIndex(double &r, const double &d) {
  if (r == 0.0) {
    return 0.;
    // throw std::runtime_error("Division by zero: radius is 0");
  }

  return (d - r) / r;
}

// Function to calculate the Petitjean Index
double calculatePetitjeanIndex(const double &r, const double &d) {
  if (d == 0.0) {
    return 0.;
    // throw std::runtime_error("Division by zero: diameter is 0");
  }

  return (d - r) / d;
}

std::vector<double> calcTopologicalIndex(const ROMol &mol) {
  Eigen::VectorXd Eccentricity = calculateEccentricity(mol);
  double Radius = calculateRadius(Eccentricity);
  double Diameter = calculateDiameter(Eccentricity);
  double TSI = calculateTopologicalShapeIndex(Radius, Diameter);
  double PJI = calculatePetitjeanIndex(Radius, Diameter);
  return {Diameter, Radius, TSI, PJI};
}


// first code slow V1
int calcPathsOfLengthN_(const ROMol &mol, int order) {
  // Extract and classify subgraphs for the current radius order
  auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
  int j = 0;
  for (const auto &[bonds, nodes, type] : classifiedPaths) {
    if (type == ChiType::Path) {
      j++;
    }
  }
  return j;
}

// second code  V2 faster than extractAndClassifyPaths
int calcPathsOfLengthN(const ROMol &mol, int order) {
  // Extract and classify subgraphs for the current radius order
  int j = 0;

  auto paths = findAllPathsOfLengthN(
      mol, order + 1, false, false, -1,
      false);  // Atoms indices we need +1 at it is linear path order bonds
               // equal order+1 atoms paths!!!!

  for (const auto &atomPath : paths) {
    std::unordered_set<int> visitedAtoms;  // Set to track visited atoms
    bool isDuplicate = false;

    for (size_t i = 0; i < atomPath.size(); ++i) {
      if (visitedAtoms.count(atomPath[i])) {
        isDuplicate = true;
        break;
      }
      visitedAtoms.insert(atomPath[i]);
    }

    if (!isDuplicate) {
      j++;  // Increment true path count and exclude Chain
    }
  }
  return j;
}

// thrid code no need for Iterator ... so very slightly faster then V2
int calcPathsOfLengthN__(const ROMol &mol, int order) {
  // Count for true paths
  int pathCount = 0;

  // Find linear atom paths (order + 1 atoms for a path of length `order`)
  auto atomPaths =
      findAllPathsOfLengthN(mol, order + 1, false, false, -1, false);

  for (const auto &atomPath : atomPaths) {
    // Check for duplicates in the atom path
    std::unordered_set<int> visitedAtoms(atomPath.begin(), atomPath.end());
    bool isChain = (visitedAtoms.size() < atomPath.size());

    if (!isChain) {
      pathCount++;  // Increment path count for true linear paths
    }
  }

  return pathCount;
}

// Function to convert bond paths into atom paths
std::vector<int> bondPathToAtomPath(const ROMol &mol,
                                    const std::vector<int> &bondPath) {
  std::unordered_set<int> visitedAtoms;
  std::vector<int> atomPath;

  for (int bondIdx : bondPath) {
    const auto *bond = mol.getBondWithIdx(bondIdx);
    int begin = bond->getBeginAtomIdx();
    int end = bond->getEndAtomIdx();

    if (visitedAtoms.insert(begin).second) atomPath.push_back(begin);
    if (visitedAtoms.insert(end).second) atomPath.push_back(end);
  }

  return atomPath;
}

struct pathHash {
  std::size_t operator()(const std::vector<int> &path) const {
    std::size_t seed = 0;
    for (int i : path) {
      seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

bool detectCycle(const std::vector<int> &atomPath) {
  std::unordered_set<int> visited;
  for (int atomIdx : atomPath) {
    if (visited.count(atomIdx)) {
      return true;  // Cycle detected
    }
    visited.insert(atomIdx);
  }
  return false;  // No cycle
}


std::vector<std::pair<std::vector<int>, ChiType>>
extractAndClassifyPathsAndSubgraphs(const ROMol &mol, unsigned int targetLength,
                                    bool useHs) {
  std::vector<std::pair<std::vector<int>, ChiType>> results;

  // Step 1: Classify linear atom paths
  auto atomPaths =
      findAllPathsOfLengthN(mol, targetLength + 1, useHs, false, -1, false);

  for (const auto &atomPath : atomPaths) {
    bool isChain = detectCycle(atomPath);  // Detect cycles in atom path
    ChiType type = isChain ? ChiType::Chain : ChiType::Path;

    results.emplace_back(atomPath, type);
  }

  // Step 2: Classify bond subgraphs
  auto bondSubgraphs = findAllSubgraphsOfLengthN(mol, targetLength, useHs, -1);

  for (const auto &bondPath : bondSubgraphs) {
    // Convert bond subgraph to atom path
    std::vector<int> atomPath = bondPathToAtomPath(mol, bondPath);

    // Check if the atom path was already classified as Path or Chain
    auto it =
        std::find_if(results.begin(), results.end(), [&](const auto &pair) {
          return pair.first == atomPath && (pair.second == ChiType::Path ||
                                            pair.second == ChiType::Chain);
        });

    if (it != results.end()) {
      continue;  // Skip already classified paths
    }

    // Classify the bond subgraph
    ChiType type = classifySubgraph(mol, bondPath);

    // Only add the subgraph if it is not a Path or Chain
    if (type != ChiType::Path && type != ChiType::Chain) {
      results.emplace_back(atomPath, type);
    }
  }

  return results;
}

// Calculate path counts and weighted product
std::pair<int, double> calculatePathCount(const ROMol &mol, int order) {
  int L = 0;           // Path count
  double piSum = 0.0;  // Weighted bond product sum

  // Get all paths of the given length
  auto paths = findAllPathsOfLengthN(
      mol, order + 1, false, false, -1,
      false);  // Atoms indices we need +1 at it is linear path order bonds
               // equal order+1 atoms paths!!!!

  for (const auto &atomPath : paths) {
    std::unordered_set<int> visitedAtoms;  // Set to track visited atoms
    bool isDuplicate = false;

    double bondProduct = 1.0;  // Initialize bond product

    for (size_t i = 0; i < atomPath.size(); ++i) {
      if (visitedAtoms.count(atomPath[i])) {
        isDuplicate = true;
        break;
      }
      visitedAtoms.insert(atomPath[i]);

      if (i > 0) {
        const auto *bond =
            mol.getBondBetweenAtoms(atomPath[i - 1], atomPath[i]);
        bondProduct *= bond->getBondTypeAsDouble();
      }
    }

    if (!isDuplicate) {
      L++;                   // Increment path count
      piSum += bondProduct;  // Add bond product to the sum
    }
  }

  return {L, piSum};
}

// Main function to calculate path descriptors
std::vector<double> calcPathCount(const ROMol &mol) {
  std::vector<double> results(21, 0.0);  // Output vector
  int totalMPC = mol.getNumAtoms();      // Initialize MPC1
  double cumulativePiSum =
      static_cast<double>(mol.getNumAtoms());  // Initialize piPC1

  for (int order = 1; order <= 10; ++order) {
    auto [L, piSum] = calculatePathCount(mol, order);

    if (order > 1) {
      results[order - 2] = L;  // MPC2-MPC10 in indices 0-8
    }
    results[10 + order - 1] =
        std::log(piSum + 1.0);  // piPC1-piPC10 in indices 10-19

    totalMPC += L;
    cumulativePiSum += piSum;
  }

  results[9] = totalMPC;                        // Total MPC in index 9
  results[20] = std::log(cumulativePiSum + 1);  // Total piPC in index 20

  return results;
}

// caution this is orignial kappa shape index so without the alpha term p 428 /
// 429!

double kappa1Helperfix(double P1, double A) {
  double denom = P1;
  double kappa = 0.0;
  if (denom) {
    kappa = (A) * (A - 1) * (A - 1) / (denom * denom);
  }
  return kappa;
}
double calcKappa1fix(const ROMol &mol) {
  double P1 = mol.getNumBonds();
  double A = mol.getNumHeavyAtoms();
  double kappa = kappa1Helperfix(P1, A);
  return kappa;
}

double kappa2Helperfix(double P2, double A) {
  double denom = (P2) * (P2);
  double kappa = 0.0;
  if (denom) {
    kappa = (A - 1) * (A - 2) * (A - 2) / denom;
  }
  return kappa;
}
double calcKappa2fix(const ROMol &mol) {
  // double P2  = findAllPathsOfLengthN(mol, 2).size();
  double P2 = static_cast<double>(calcPathsOfLengthN(mol, 2));

  double A = mol.getNumHeavyAtoms();
  double kappa = kappa2Helperfix(P2, A);
  return kappa;
}

double kappa3Helperfix(double P3, int A) {
  double denom = P3 * P3;
  double kappa = 0.0;
  if (denom) {
    if (A % 2) {
      kappa = (A - 1) * (A - 3) * (A - 3) / denom;
    } else {
      kappa = (A - 2) * (A - 2) * (A - 3) / denom;
    }
  }
  return kappa;
}

double calcKappa3fix(const ROMol &mol) {
  // double P3 = findAllPathsOfLengthN(mol, 3).size();
  double P3 = static_cast<double>(calcPathsOfLengthN(mol, 3));

  double A = mol.getNumHeavyAtoms();
  double kappa = kappa3Helperfix(P3, A);
  return kappa;
}

std::vector<double> calcKappaShapeIndex(const ROMol &mol) {
  std::vector<double> results(3, 0.0);
  results[0] = calcKappa1fix(mol);

  results[1] = calcKappa2fix(mol);

  results[2] = calcKappa3fix(mol);

  return results;
}

// it is not identical to "getRKHE"
void hkDeltas(const ROMol &mol, std::vector<double> &deltas, bool force) {
  PRECONDITION(deltas.size() >= mol.getNumAtoms(), "bad vector size");
  if (!force && mol.hasProp("_connectivityHKDeltas")) {
    mol.getProp("_connectivityHKDeltas", deltas);
    return;
  }

  // Compute the valence electrons for each atom
  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    const Atom *at = mol[*atBegin];
    double dv = getValenceElectrons(*at);  // Valence electrons

    // Replace previous delta calculation with valence electrons
    if (dv != 0.0) {
      deltas[at->getIdx()] = 1. / sqrt(dv);
    } else {
      deltas[at->getIdx()] = 0.0;
    }
    ++atBegin;
  }

  mol.setProp("_connectivityHKDeltas", deltas, true);
}

void nVals(const ROMol &mol, std::vector<double> &nVs, bool force) {
  PRECONDITION(nVs.size() >= mol.getNumAtoms(), "bad vector size");
  if (!force && mol.hasProp("_connectivityNVals")) {
    mol.getProp("_connectivityNVals", nVs);
    return;
  }

  // Compute the sigma electrons for each atom
  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    const Atom *at = mol[*atBegin];
    double sigma_electrons = getSigmaElectrons(*at);  // Sigma electrons

    if (sigma_electrons != 0.0) {
      nVs[at->getIdx()] = 1. / sqrt(sigma_electrons);
    } else {
      nVs[at->getIdx()] = 0.0;
    }
    ++atBegin;
  }

  mol.setProp("_connectivityNVals", nVs, true);
}

// optional alphaKappa version (for HeteroAtom differentiation)

double getAlpha(const Atom &atom, bool &found) {
  double res = 0.0;
  found = false;
  switch (atom.getAtomicNum()) {
    case 1:
      res = 0.0;
      found = true;
      break;
    case 6:
      switch (atom.getHybridization()) {
        case Atom::SP:
          res = -0.22;
          found = true;
          break;
        case Atom::SP2:
          res = -0.13;
          found = true;
          break;
        default:
          res = 0.00;
          found = true;
      };
      break;
    case 7:
      switch (atom.getHybridization()) {
        case Atom::SP:
          res = -0.29;
          found = true;
          break;
        case Atom::SP2:
          res = -0.20;
          found = true;
          break;
        default:
          res = -0.04;
          found = true;
          break;
      };
      break;
    case 8:
      switch (atom.getHybridization()) {
        case Atom::SP2:
          res = -0.20;
          found = true;
          break;
        default:
          res = -0.04;
          found = true;
          break;
      };
      break;
    case 9:
      switch (atom.getHybridization()) {
        default:
          res = -0.07;
          found = true;
          break;
      };
      break;
    case 15:
      switch (atom.getHybridization()) {
        case Atom::SP2:
          res = 0.30;
          found = true;
          break;
        default:
          res = 0.43;
          found = true;
          break;
      };
      break;
    case 16:
      switch (atom.getHybridization()) {
        case Atom::SP2:
          res = 0.22;
          found = true;
          break;
        default:
          res = 0.35;
          found = true;
          break;
      };
      break;
    case 17:
      switch (atom.getHybridization()) {
        default:
          res = 0.29;
          found = true;
          break;
      };
      break;
    case 35:
      switch (atom.getHybridization()) {
        default:
          res = 0.48;
          found = true;
          break;
      };
      break;
    case 53:
      switch (atom.getHybridization()) {
        default:
          res = 0.73;
          found = true;
          break;
      };
      break;
    default:
      break;
  }
  return res;
}

// old style from connectivity rdkit deps file
double calcHallKierAlpha(const ROMol &mol) {
  const PeriodicTable *tbl = PeriodicTable::getTable();
  double alphaSum = 0.0;
  double rC = tbl->getRb0(6);
  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    const Atom *at = mol[*atBegin];
    ++atBegin;
    unsigned int n = at->getAtomicNum();
    if (!n) continue;
    bool found;
    double alpha = getAlpha(*at, found);
    if (!found) {
      double rA = tbl->getRb0(n);
      alpha = rA / rC - 1.0;
    }
    alphaSum += alpha;
  }
  return alphaSum;
};

double alphakappa1Helper(double P1, double A, double alpha) {
  double denom = P1 + alpha;
  double kappa = 0.0;
  if (denom) {
    kappa = (A + alpha) * (A + alpha - 1) * (A + alpha - 1) / (denom * denom);
  }
  return kappa;
}
double alphakappa2Helper(double P2, double A, double alpha) {
  double denom = (P2 + alpha) * (P2 + alpha);
  double kappa = 0.0;
  if (denom) {
    kappa = (A + alpha - 1) * (A + alpha - 2) * (A + alpha - 2) / denom;
  }
  return kappa;
}

double alphakappa3Helperfix(double P3, int A, double alpha) {
  double denom = (P3 + alpha) * (P3 + alpha);
  double kappa = 0.0;
  if (denom) {
    if (A % 2) {
      kappa = (A + alpha - 1) * (A + alpha - 3) * (A + alpha - 3) / denom;
    } else {
      kappa = (A + alpha - 2) * (A + alpha - 2) * (A + alpha - 3) / denom;
    }
  }
  return kappa;
}

double calcalphaKappa1(const ROMol &mol) {
  double P1 = mol.getNumBonds();
  double A = mol.getNumHeavyAtoms();
  double alpha = calcHallKierAlpha(mol);
  double kappa = alphakappa1Helper(P1, A, alpha);
  return kappa;
}
double calcalphaKappa2(const ROMol &mol) {
  double P2 = static_cast<double>(calcPathsOfLengthN(mol, 2));
  double A = mol.getNumHeavyAtoms();
  double alpha = calcHallKierAlpha(mol);
  double kappa = alphakappa2Helper(P2, A, alpha);
  return kappa;
}
double calcalphaKappa3fix(const ROMol &mol) {
  double P3 = static_cast<double>(calcPathsOfLengthN(mol, 3));
  int A = mol.getNumHeavyAtoms();
  double alpha = calcHallKierAlpha(mol);
  double kappa = alphakappa3Helperfix(P3, A, alpha);
  return kappa;
}

double Flexibility(const ROMol &mol) {
  double AK1 = calcalphaKappa1(mol);
  double AK2 = calcalphaKappa2(mol);
  int numHeavyAtom = mol.getNumHeavyAtoms();
  return AK1 * AK2 / static_cast<double>(numHeavyAtom);
}

std::vector<double> calcFlexibility(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = Flexibility(mol);
  return res;
}

std::vector<double> calcAlphaKappaShapeIndex(const ROMol &mol) {
  std::vector<double> results(3, 0.0);
  results[0] = calcalphaKappa1(mol);

  results[1] = calcalphaKappa2(mol);

  results[2] = calcalphaKappa3fix(mol);

  return results;
}

// Compute for Path (Xp)
double calcXpDescriptor(const ROMol &mol, unsigned int order,
                        bool useValenceElectrons, bool force) {
  std::vector<double> electronVals(mol.getNumAtoms());
  if (useValenceElectrons) {
    hkDeltas(mol, electronVals, force);  // Valence Electrons (dv)
  } else {
    nVals(mol, electronVals, force);  // Sigma Electrons (d)
  }

  PATH_LIST ps = findAllPathsOfLengthN(mol, order, false);
  double result = 0.0;
  for (const auto &p : ps) {
    double accum = 1.0;
    for (int aidx : p) {
      accum *=
          electronVals[aidx];  // Multiply the electron values along the path
    }
    result += accum;
  }
  return result;
}

// Compute for Chain (Xch)
double calcXchDescriptor(const ROMol &mol, unsigned int order,
                         bool useValenceElectrons, bool force) {
  std::vector<double> electronVals(mol.getNumAtoms());
  if (useValenceElectrons) {
    hkDeltas(mol, electronVals, force);  // Valence Electrons (dv)
  } else {
    nVals(mol, electronVals, force);  // Sigma Electrons (d)
  }

  PATH_LIST ps = findAllPathsOfLengthN(
      mol, order,
      true);  // Chain specific path finding (true means chain only ???)
  double result = 0.0;
  for (const auto &p : ps) {
    double accum = 1.0;
    for (int aidx : p) {
      accum *= electronVals[aidx];
    }
    result += accum;
  }
  return result;
}

// Function to compute all descriptors for different ranges
double computeDescriptors(const ROMol &mol, bool force) {
  // Calculate for Xp (Path)
  for (unsigned int order = 0; order <= 7; ++order) {
    // Compute Xp-0d, Xp-1d, ..., Xp-7d
    double result_d = calcXpDescriptor(mol, order, false, force);
    std::cout << "Xp-" << order << "d: " << result_d << std::endl;

    // Compute Xp-0dv, Xp-1dv, ..., Xp-7dv
    double result_dv = calcXpDescriptor(mol, order, true, force);
    std::cout << "Xp-" << order << "dv: " << result_dv << std::endl;

    // Compute AXp-0d, AXp-1d, ..., AXp-7d
    // Assuming average means some form of averaging over the molecule
    double avg_result_d = result_d / mol.getNumAtoms();
    std::cout << "AXp-" << order << "d: " << avg_result_d << std::endl;

    // Compute AXp-0dv, AXp-1dv, ..., AXp-7dv
    double avg_result_dv = result_dv / mol.getNumAtoms();
    std::cout << "AXp-" << order << "dv: " << avg_result_dv << std::endl;
  }

  // Calculate for Xch (Chain)
  for (unsigned int order = 3; order <= 7; ++order) {
    double result_d = calcXchDescriptor(mol, order, false, force);
    std::cout << "Xch-" << order << "d: " << result_d << std::endl;

    double result_dv = calcXchDescriptor(mol, order, true, force);
    std::cout << "Xch-" << order << "dv: " << result_dv << std::endl;
  }

  return 1.0;
}

// test function for Chi
std::vector<int> calcChiPath(const ROMol &mol) {
  std::vector<int> res(3, 0);

  for (int i = 1; i <= 3; i++) {
    // Extract and classify subgraphs for the current radius
    auto classifiedPaths = extractAndClassifyPaths(mol, i, false);

    // std::vector<std::vector<int>> atomPaths;

    // Filter only the paths classified as ChiType::Path
    int j = 0;
    for (const auto &[bonds, nodes, type] : classifiedPaths) {
      if (type == ChiType::Path) {
        j++;
        // atomPaths.emplace_back(nodes.begin(), nodes.end());
      }
    }

    // Count the number of unique atom paths for this radius
    res[i - 1] = j;  // atomPaths.size();
  }

  return res;
}

// Function to compute Chi descriptor for a specific ChiType
double computeChiDescriptor(const ROMol &, const PATH_LIST &subgraphs,
                            const std::vector<double> &P,
                            bool averaged = false) {
  double chi_value = 0.0;

  for (const auto &nodes : subgraphs) {
    double c = 1.0;
    for (int node : nodes) {
      c *= P[node];  // Calculate the product of atomic properties
    }

    if (c <= 0) {
      std::cerr << "Error: Some properties are less than or equal to zero."
                << std::endl;
      return 0.0;  // Return 0 if any property is <= 0
    }

    chi_value += std::pow(c, -0.5);
  }

  // Optionally average the value by the number of subgraphs if `averaged` is
  // true
  if (averaged && !subgraphs.empty()) {
    chi_value /= subgraphs.size();
  }

  return chi_value;
}

//// DetourMatrix code here

// Function to compute longest path using DFS
void longestSimplePath(int u, const Eigen::MatrixXd &adjMatrix,
                       std::vector<double> &result, std::vector<bool> &visited,
                       double distance, int N) {
  visited[u] = true;
  result[u] = std::max(result[u], distance);

  for (int v = 0; v < N; ++v) {
    if (adjMatrix(u, v) > 0 && !visited[v]) {
      longestSimplePath(v, adjMatrix, result, visited,
                        distance + adjMatrix(u, v), N);
    }
  }

  visited[u] = false;
}

// Function to compute the detour matrix using an adjacency matrix
Eigen::MatrixXd computeDetourMatrix(int N, const Eigen::MatrixXd &adjMatrix) {
  Eigen::MatrixXd detourMatrix(N, N);
  detourMatrix.setZero();  // Initialize the matrix with zeros

  // Compute longest path from each node
  for (int i = 0; i < N; ++i) {
    std::vector<double> result(N, 0.0);
    std::vector<bool> visited(N, false);
    longestSimplePath(i, adjMatrix, result, visited, 0.0, N);

    for (int j = 0; j < N; ++j) {
      detourMatrix(i, j) = result[j];
    }
  }

  return detourMatrix;
}

// Compute Detour Index
int computeDetourIndex(const Eigen::MatrixXd &detourmatrix) {
  double sum = detourmatrix.sum();
  // int num_atoms = detourmatrix.rows();
  return static_cast<int>(0.5 * sum);
}

std::vector<double> calcDetourMatrixDescs(const ROMol &mol) {
  Eigen::MatrixXd AdjMat = calculateAdjacencyMatrix(mol);

  int numAtoms = mol.getNumAtoms();
  Eigen::Matrix detourMatrix = computeDetourMatrix(numAtoms, AdjMat);

  Eigen::VectorXd eigenvalues;

  Eigen::MatrixXd eigenvectors;

  compute_eigenvalues_and_eigenvectors(detourMatrix, eigenvalues, eigenvectors);

  // Compute descriptors
  double Sp_Abs = spAbs(eigenvalues);
  double Sp_Max = spMax(eigenvalues);
  double Sp_Diam = spDiam(eigenvalues);
  double Sp_Mean = spMean(eigenvalues);
  double Sp_AD = spAD(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE(eigenvalues);

  std::vector<std::pair<int, int>> bonds;

  // Iterate over all bonds in the molecule
  for (const auto &bond : mol.bonds()) {
    int beginAtomIdx = bond->getBeginAtomIdx();
    int endAtomIdx = bond->getEndAtomIdx();
    bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
  }

  // Calculate the descriptors
  double ve1 = VE1(detourMatrix, eigenvalues, eigenvectors);
  double ve2 = VE2(detourMatrix, numAtoms, eigenvalues, eigenvectors);
  double ve3 = VE3(detourMatrix, numAtoms, eigenvalues, eigenvectors);
  double vr1 = VR1(detourMatrix, bonds, eigenvalues, eigenvectors);
  double vr2 = VR2(detourMatrix, bonds, numAtoms, eigenvalues, eigenvectors);
  double vr3 = VR3(detourMatrix, bonds, numAtoms, eigenvalues, eigenvectors);
  double detourindex = static_cast<double>(computeDetourIndex(detourMatrix));

  double sm1 = SM1(detourMatrix);
  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3,    detourindex};
}

/// much faster lapack version

// Function to compute longest path using DFS
void longestSimplePathL(int u,
                        const std::vector<std::vector<double>> &adjMatrix,
                        std::vector<double> &result, std::vector<bool> &visited,
                        double distance, int N) {
  visited[u] = true;
  result[u] = std::max(result[u], distance);

  for (int v = 0; v < N; ++v) {
    if (adjMatrix[u][v] > 0 && !visited[v]) {
      longestSimplePathL(v, adjMatrix, result, visited,
                         distance + adjMatrix[u][v], N);
    }
  }

  visited[u] = false;
}

// Function to compute the detour matrix using an adjacency matrix
std::vector<std::vector<double>> computeDetourMatrixL(
    int N, const std::vector<std::vector<double>> &adjMatrix) {
  std::vector<std::vector<double>> detourMatrix(N, std::vector<double>(N, 0.0));

  // Compute longest path from each node
  for (int i = 0; i < N; ++i) {
    std::vector<double> result(N, 0.0);
    std::vector<bool> visited(N, false);
    longestSimplePathL(i, adjMatrix, result, visited, 0.0, N);

    for (int j = 0; j < N; ++j) {
      detourMatrix[i][j] = result[j];
    }
  }

  return detourMatrix;
}

// Compute Detour Index
int computeDetourIndexL(const std::vector<std::vector<double>> &detourMatrix) {
  double sum = 0.0;
  for (size_t i = 0; i < detourMatrix.size(); ++i) {
    for (size_t j = 0; j < detourMatrix[i].size(); ++j) {
      sum += detourMatrix[i][j];
    }
  }
  return static_cast<int>(0.5 * sum);
}

// Function to get the Adjacency Matrix using standard C++ structures
std::vector<std::vector<double>> calculateAdjacencyMatrixL(const ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();

  // Initialize a 2D vector with zeros
  std::vector<std::vector<double>> adjMatrix(nAtoms,
                                             std::vector<double>(nAtoms, 0.0));

  // Populate the adjacency matrix using RDKit's bond information
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j < nAtoms; ++j) {
      const Bond *bond = mol.getBondBetweenAtoms(i, j);
      adjMatrix[i][j] = (bond != nullptr) ? 1.0 : 0.0;
    }
  }

  return adjMatrix;
}

std::vector<double> calcDetourMatrixDescsL(const ROMol &mol) {
  auto adjMatrix = calculateAdjacencyMatrixL(
      mol);  // Convert adjacency matrix computation to LAPACK
  int numAtoms = mol.getNumAtoms();

  auto detourMatrix = computeDetourMatrixL(numAtoms, adjMatrix);

  // Convert detourMatrix to 1D array in column-major order for LAPACK
  std::vector<double> flatMatrix(numAtoms * numAtoms);
  for (int i = 0; i < numAtoms; ++i) {
    for (int j = 0; j < numAtoms; ++j) {
      flatMatrix[j * numAtoms + i] = detourMatrix[i][j];
    }
  }

  // Compute eigenvalues and eigenvectors using LAPACK
  std::vector<double> eigenvalues(numAtoms);
  std::vector<std::vector<double>> eigenvectors(numAtoms,
                                                std::vector<double>(numAtoms));

  compute_eigenvalues_and_eigenvectorsL(detourMatrix, eigenvalues,
                                        eigenvectors);

  // Compute descriptors
  double Sp_Abs = spAbsL(eigenvalues);
  double Sp_Max = spMaxL(eigenvalues);
  double Sp_Diam = spDiamL(eigenvalues);
  double Sp_Mean = spMeanL(eigenvalues);
  double Sp_AD = spADL(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE_stable(eigenvalues);

  std::vector<std::pair<int, int>> bonds;

  // Iterate over all bonds in the molecule
  for (const auto &bond : mol.bonds()) {
    int beginAtomIdx = bond->getBeginAtomIdx();
    int endAtomIdx = bond->getEndAtomIdx();
    bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
  }

  // Calculate descriptors
  double ve1 = VE1L(eigenvectors);
  double ve2 = VE2L(ve1, numAtoms);
  double ve3 = VE3L(ve1, numAtoms);
  double vr1 = VR1L(eigenvectors, bonds);
  double vr2 = VR2L(vr1, numAtoms);
  double vr3 = VR3L(vr1, numAtoms);
  double detourindex = static_cast<double>(computeDetourIndexL(detourMatrix));

  double sm1 = SM1L(detourMatrix);
  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3,    detourindex};
}

//// AdjacencyMatrix

std::vector<double> calcAdjMatrixDescs(const ROMol &mol) {
  Eigen::MatrixXd AdjMat = calculateAdjacencyMatrix(mol);

  int numAtoms = mol.getNumAtoms();

  Eigen::VectorXd eigenvalues;

  Eigen::MatrixXd eigenvectors;

  compute_eigenvalues_and_eigenvectors(AdjMat, eigenvalues, eigenvectors);

  std::vector<std::pair<int, int>> bonds;

  // Iterate over all bonds in the molecule
  for (const auto &bond : mol.bonds()) {
    int beginAtomIdx = bond->getBeginAtomIdx();
    int endAtomIdx = bond->getEndAtomIdx();
    bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
  }
  /*
      // Compute descriptors
      std::cout << "eigenvalues from Eigen: ";
      for (auto ei : eigenvalues) {
          std::cout << ei << " ";
      }
      std::cout << std::endl;
  */
  double Sp_Abs = spAbs(eigenvalues);
  double Sp_Max = spMax(eigenvalues);
  double Sp_Diam = spDiam(eigenvalues);
  double Sp_Mean =
      spMean(eigenvalues);  // tmp values not needed to export as result
  double Sp_AD = spAD(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE(
      eigenvalues);  // this one is not correct ??? do we need to add the bond
  double ve1 = VE1(AdjMat, eigenvalues, eigenvectors);
  double ve2 = VE2(AdjMat, numAtoms, eigenvalues, eigenvectors);
  double ve3 = VE3(AdjMat, numAtoms, eigenvalues, eigenvectors);
  double vr1 = VR1(AdjMat, bonds, eigenvalues, eigenvectors);
  double vr2 = VR2(AdjMat, bonds, numAtoms, eigenvalues, eigenvectors);
  double vr3 = VR3(AdjMat, bonds, numAtoms, eigenvalues, eigenvectors);
  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3};
}

std::vector<double> calcAdjMatrixDescsL(const ROMol &mol) {
  auto AdjMat =
      calculateAdjacencyMatrixL(mol);  // Use LAPACK-compatible adjacency matrix

  int numAtoms = mol.getNumAtoms();

  std::vector<double> eigenvalues;
  std::vector<std::vector<double>> eigenvectors;

  compute_eigenvalues_and_eigenvectorsL(
      AdjMat, eigenvalues, eigenvectors);  // LAPACK eigen computation

  std::vector<std::pair<int, int>> bonds;
  for (const auto &bond : mol.bonds()) {
    bonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
  }
  /*
          // Compute descriptors
          std::cout << "eigenvalues from Lapack: ";
          for (auto ei : eigenvalues) {
              std::cout << ei << " ";
          }
          std::cout << std::endl;
  */
  // Compute descriptors
  double Sp_Abs = spAbsL(eigenvalues);
  double Sp_Max = spMaxL(eigenvalues);
  double Sp_Diam = spDiamL(eigenvalues);
  double Sp_Mean = spMeanL(eigenvalues);
  double Sp_AD = spADL(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE_stable(eigenvalues);

  double ve1 = VE1L(eigenvectors);
  double ve2 = VE2L(ve1, numAtoms);
  double ve3 = VE3L(ve1, numAtoms);
  double vr1 = VR1L(eigenvectors, bonds);
  double vr2 = VR2L(vr1, numAtoms);
  double vr3 = VR3L(vr1, numAtoms);
  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3};
}

std::vector<double> calcDistMatrixDescs(const ROMol &mol) {
  Eigen::MatrixXd DMat = calculateDistanceMatrix(mol);

  int numAtoms = mol.getNumAtoms();

  Eigen::VectorXd eigenvalues;

  Eigen::MatrixXd eigenvectors;

  compute_eigenvalues_and_eigenvectors(DMat, eigenvalues, eigenvectors);

  std::vector<std::pair<int, int>> bonds;

  // Iterate over all bonds in the molecule
  for (const auto &bond : mol.bonds()) {
    int beginAtomIdx = bond->getBeginAtomIdx();
    int endAtomIdx = bond->getEndAtomIdx();
    bonds.push_back(std::make_pair(beginAtomIdx, endAtomIdx));
  }

  // Compute descriptors
  double Sp_Abs = spAbs(eigenvalues);
  double Sp_Max = spMax(eigenvalues);
  double Sp_Diam = spDiam(eigenvalues);
  double Sp_Mean =
      spMean(eigenvalues);  // tmp values not needed to export as result
  double Sp_AD = spAD(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE(
      eigenvalues);  // this one is not correct ??? do we need to add the bond
  double ve1 = VE1(DMat, eigenvalues, eigenvectors);
  double ve2 = VE2(DMat, numAtoms, eigenvalues, eigenvectors);
  double ve3 = VE3(DMat, numAtoms, eigenvalues, eigenvectors);
  double vr1 = VR1(DMat, bonds, eigenvalues, eigenvectors);
  double vr2 = VR2(DMat, bonds, numAtoms, eigenvalues, eigenvectors);
  double vr3 = VR3(DMat, bonds, numAtoms, eigenvalues, eigenvectors);

  // SM1(DMat) intentionally omitted for parity with L variant
  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3};
}

std::vector<std::vector<double>> calculateDistanceMatrixL(const ROMol &mol) {
  unsigned int nAtoms = mol.getNumAtoms();

  // Get the distance matrix using RDKit's MolOps::getDistanceMat
  double *distanceMat = MolOps::getDistanceMat(
      mol, false, false, false);  // No bond order, no weights, no hydrogens

  // Convert the raw pointer to a 2D vector
  std::vector<std::vector<double>> distMatrix(nAtoms,
                                              std::vector<double>(nAtoms, 0.0));
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j < nAtoms; ++j) {
      distMatrix[i][j] = distanceMat[i * nAtoms + j];
    }
  }

  return distMatrix;
}

std::vector<double> calcDistMatrixDescsL(const ROMol &mol) {
  auto DMat =
      calculateDistanceMatrixL(mol);  // Use LAPACK-compatible distance matrix

  int numAtoms = mol.getNumAtoms();

  std::vector<double> eigenvalues;
  std::vector<std::vector<double>> eigenvectors;

  compute_eigenvalues_and_eigenvectorsL(
      DMat, eigenvalues, eigenvectors);  // LAPACK eigen computation

  std::vector<std::pair<int, int>> bonds;
  for (const auto &bond : mol.bonds()) {
    bonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
  }

  // Compute descriptors
  double Sp_Abs = spAbsL(eigenvalues);
  double Sp_Max = spMaxL(eigenvalues);
  double Sp_Diam = spDiamL(eigenvalues);
  double Sp_Mean = spMeanL(eigenvalues);
  double Sp_AD = spADL(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE_stable(eigenvalues);
  double ve1 = VE1L(eigenvectors);
  double ve2 = VE2L(ve1, numAtoms);
  double ve3 = VE3L(ve1, numAtoms);
  double vr1 = VR1L(eigenvectors, bonds);
  double vr2 = VR2L(vr1, numAtoms);
  double vr3 = VR3L(vr1, numAtoms);
  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3};
}

// Molecular Distance Edge Calculation
std::vector<double> calcMolecularDistanceEdgeDescs(const ROMol &mol) {
  using namespace RDKit;

  // Initialize result container
  std::vector<double> results;

  // Supported atomic numbers and valence pairs
  std::vector<int> atomic_nums = {6, 8, 7};  // Carbon, Oxygen, Nitrogen

  std::map<int, std::vector<std::pair<int, int>>> valence_pairs = {
      {6,
       {{1, 1},
        {1, 2},
        {1, 3},
        {1, 4},
        {2, 2},
        {2, 3},
        {2, 4},
        {3, 3},
        {3, 4},
        {4, 4}}},
      {8, {{1, 1}, {1, 2}, {2, 2}}},
      {7, {{1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}}}};

  // Compute distance matrix
  const int num_atoms = mol.getNumAtoms();
  Eigen::MatrixXd distance_matrix = calculateDistanceMatrix(mol);
  std::vector<double> valences = calcValence(mol);

  // Iterate over atomic numbers
  for (const auto &atomic_num : atomic_nums) {
    const auto &pairs = valence_pairs[atomic_num];

    // Iterate over valence pairs
    for (const auto &[valence1, valence2] : pairs) {
      std::vector<double> Dv;

      // Filter distances based on valence and atomic number
      for (int i = 0; i < num_atoms; ++i) {
        const auto atom_i = mol.getAtomWithIdx(i);
        for (int j = i + 1; j < num_atoms; ++j) {
          const auto atom_j = mol.getAtomWithIdx(j);
          if ((valences[i] == valence1 && valences[j] == valence2) ||
              (valences[j] == valence1 && valences[i] == valence2)) {
            if (atom_j->getAtomicNum() == atomic_num &&
                atom_i->getAtomicNum() == atomic_num) {
              Dv.push_back(distance_matrix(i, j));
            }
          }
        }
      }

      // Compute descriptor if valid distances are found
      int n = Dv.size();
      if (n > 0) {
        double log_sum = 0.0;
        for (const auto &d : Dv) {
          log_sum += std::log(d);
        }
        double log_dx = log_sum / (2.0 * n);
        double dx = std::exp(log_dx);
        double descriptor_value = n / (dx * dx);

        results.push_back(descriptor_value);
      } else {
        results.push_back(0.0);  // Default value when no pairs are found
      }
    }
  }

  return results;
}

////////////////////////// ExtendedTopochemicalAtom  //////////////////////////

// ExtendedTopochemicalAtom  all codes
RWMol *modifyMolecule(const ROMol &mol, bool explicitHydrogens,
                      bool saturated) {
  RWMol *newMol = new RWMol();
  std::map<unsigned int, unsigned int> atomMap;

  // Add atoms to the new molecule
  for (const auto &atom : mol.atoms()) {
    if (atom->getAtomicNum() == 1) continue;  // Skip hydrogens

    Atom *newAtom = nullptr;

    if (saturated) {
      // Preserve atomic number and formal charge
      newAtom = new Atom(atom->getAtomicNum());
      newAtom->setFormalCharge(atom->getFormalCharge());
      newAtom->setIsAromatic(atom->getIsAromatic());  // Retain aromaticity
    } else {
      // Replace all non-hydrogen atoms with carbon
      newAtom = new Atom(6);          // Carbon
      newAtom->setIsAromatic(false);  // Reset aromaticity
    }

    unsigned int newIdx = newMol->addAtom(newAtom, false, true);
    atomMap[atom->getIdx()] = newIdx;
  }

  // Add bonds to the new molecule
  for (const auto &bond : mol.bonds()) {
    unsigned int startIdx = bond->getBeginAtomIdx();
    unsigned int endIdx = bond->getEndAtomIdx();

    if (atomMap.find(startIdx) != atomMap.end() &&
        atomMap.find(endIdx) != atomMap.end()) {
      Bond::BondType bondType;

      if (saturated) {
        // Retain original bond type if at least one atom is not carbon
        const auto *startAtom = mol.getAtomWithIdx(startIdx);
        const auto *endAtom = mol.getAtomWithIdx(endIdx);
        if (startAtom->getAtomicNum() != 6 || endAtom->getAtomicNum() != 6) {
          bondType = bond->getBondType();
        } else {
          bondType = Bond::SINGLE;  // Default to single bond
        }
      } else {
        bondType = Bond::SINGLE;  // Force single bonds for reference molecules
      }

      newMol->addBond(atomMap[startIdx], atomMap[endIdx], bondType);
    }
  }

  // Attempt sanitization
  unsigned int failedOps = 0;
  try {
    MolOps::sanitizeMol(*newMol, failedOps,
                        MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_PROPERTIES);
  } catch (const MolSanitizeException &e) {
    std::cerr << "Sanitization failed: " << e.what() << "\n";
  }

  if (failedOps > 0) {
    std::cerr
        << "Sanitization failed on some operations, continuing. Failed ops: "
        << failedOps << "\n";
  }

  // Add explicit hydrogens if requested
  if (explicitHydrogens) {
    MolOps::addHs(*newMol);
  }

  // Ensure the molecule is Kekulized
  try {
    MolOps::Kekulize(*newMol, false);
  } catch (const KekulizeException &e) {
    std::cerr << "Kekulization failed: " << e.what() << "\n";
  }

  return newMol;
}

// Function to check if any heteroatom has more than 4 neighbors
bool hasExcessiveNeighbors(const ROMol &mol) {
  for (auto atom : mol.atoms()) {
    int total_neighbors = atom->getDegree() + atom->getTotalNumHs();
    if (atom->getAtomicNum() != 6 && total_neighbors > 4) {
      std::cout
          << "Skipping molecule due to heteroatom with more than 4 neighbors (Atom Index: "
          << atom->getIdx() << " Atomic Num: " << atom->getAtomicNum()
          << " Neighbors: " << total_neighbors << ")\n";
      return true;  // A heteroatom with more than 4 neighbors found
    }
  }
  return false;  // Safe to proceed
}

// CAUTION : wrong for type 4 ie saturated True use modifyMolecule instead
RWMol *cloneAndModifyMolecule(const ROMol &originalMol, bool explicitHydrogens,
                              bool saturated) {
  try {
    // Check before proceeding
    if (hasExcessiveNeighbors(originalMol) && saturated) {
      return nullptr;
    }

    // Clone the molecule
    RWMol *clonedMol = new RWMol(originalMol);

    // Iterate over atoms and modify them
    for (auto atom : clonedMol->atoms()) {
      if (atom->getAtomicNum() == 1) {
        continue;
      }

      if (saturated) {
        // Retain the atomic number but adjust formal charges if necessary ...
        atom->setFormalCharge(atom->getFormalCharge());
      } else {
        // Change heavy atoms (non-hydrogen) to Carbon (atomic number 6)
        if (atom->getAtomicNum() != 6) {
          atom->setAtomicNum(6);
        }
      }
    }

    // Iterate over bonds and modify them
    for (auto bond : clonedMol->bonds()) {
      // Adjust bond types if not saturated
      if (!saturated) {
        bond->setBondType(Bond::SINGLE);
      }
    }

    // Perform sanitization with error handling
    unsigned int failedOps = 0;
    MolOps::sanitizeMol(*clonedMol, failedOps, MolOps::SANITIZE_ALL);

    if (failedOps > 0) {
      std::cerr
          << "Sanitization failed on properties, but continuing. Failed ops: "
          << failedOps << "\n";
    }

    if (explicitHydrogens) {
      MolOps::addHs(*clonedMol);
    }

    // Ensure the molecule is Kekulized
    MolOps::Kekulize(
        *clonedMol,
        false);  // same error as before !!! keep aromatic flags please!!!

    return clonedMol;  // Return the modified clone
  } catch (const MolSanitizeException &e) {
    std::cerr << "Sanitization failed: " << e.what() << "\n";
    return nullptr;  // Return nullptr if the modification fails
  } catch (const std::exception &e) {
    std::cerr << "Error cloning and modifying molecule: " << e.what() << "\n";
    return nullptr;
  }
}

// Helper function to get the core count
double getCoreCount(const Atom &atom) {
  int Z = atom.getAtomicNum();
  if (Z == 1) {
    return 0.0;
  }

  const PeriodicTable *tbl = PeriodicTable::getTable();
  double Zv = tbl->getNouterElecs(Z);
  int PN = GetPrincipalQuantumNumber(Z);

  return (Z - Zv) / (Zv * (PN - 1.));
}

// Helper function to calculate eta_epsilon
double getEtaEpsilon(const Atom &atom) {
  const PeriodicTable *tbl = PeriodicTable::getTable();
  double Zv = tbl->getNouterElecs(atom.getAtomicNum());
  return 0.3 * Zv - getCoreCount(atom);
}

// Helper function to calculate eta_beta_sigma
double getEtaBetaSigma(const Atom &atom) {
  double e = getEtaEpsilon(atom);
  double betaSigma = 0.0;

  for (const auto &neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
    if (neighbor->getAtomicNum() != 1) {
      double neighborEpsilon = getEtaEpsilon(*neighbor);
      betaSigma += (std::abs(neighborEpsilon - e) <= 0.3) ? 0.5 : 0.75;
    }
  }
  return betaSigma;
}

// Helper function to calculate non-sigma contribution
double getEtaNonSigmaContribute(const Bond &bond) {
  if (bond.getBondType() == Bond::BondType::SINGLE) {
    return 0.0;
  }
  double f = 1.0;
  if (bond.getBondType() == Bond::BondType::TRIPLE) {
    f = 2.0;
  }
  const Atom *a = bond.getBeginAtom();
  const Atom *b = bond.getEndAtom();
  double dEps = std::abs(getEtaEpsilon(*a) - getEtaEpsilon(*b));
  double y = 1.0;
  if (bond.getIsAromatic()) {
    y = 2.0;
  } else if (dEps > 0.3) {
    y = 1.5;
  }
  return y * f;
}

bool isAtomInRing(const Atom &atom) {
  const RingInfo *ringInfo = atom.getOwningMol().getRingInfo();
  if (ringInfo && ringInfo->isInitialized()) {
    return ringInfo->numAtomRings(atom.getIdx()) > 0;
  }
  return false;  // Not in a ring if no ring info is available
}

// Helper function to calculate eta_beta_delta
double getEtaBetaDelta(const Atom &atom) {
  const PeriodicTable *tbl = PeriodicTable::getTable();
  if (atom.getIsAromatic() || isAtomInRing(atom) ||
      (tbl->getNouterElecs(atom.getAtomicNum()) - atom.getTotalValence() <=
       0)) {
    return 0.0;
  }

  for (const auto &neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
    if (neighbor->getIsAromatic()) {
      return 0.5;
    }
  }
  return 0.0;
}

// Helper function to get the other atom in a bond
const Atom *getOtherAtom(const Bond &bond, const Atom &atom) {
  if (bond.getBeginAtom() == &atom) {
    return bond.getEndAtom();
  } else {
    return bond.getBeginAtom();
  }
}

// Helper function to calculate eta_beta_non_sigma
double getEtaBetaNonSigma(const Atom &atom) {
  double betaNonSigma = 0.0;

  for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
    const Atom *otherAtom = getOtherAtom(*bond, atom);
    if (otherAtom->getAtomicNum() != 1) {
      betaNonSigma += getEtaNonSigmaContribute(*bond);
    }
  }
  return betaNonSigma;
}

// Helper function to calculate eta_gamma
double getEtaGamma(const Atom &atom) {
  double beta =
      getEtaBetaSigma(atom) + getEtaBetaNonSigma(atom) + getEtaBetaDelta(atom);
  if (beta == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return getCoreCount(atom) / beta;
}

//  Etacorecount of reference can be merge with next function with an additional
//  parameter if needed...
double calculateEtaCoreCountRef(const ROMol &mol, bool averaged) {
  const ROMol *targetMol = &mol;  // Default to the input molecule
  // ROMol* molWithHs = nullptr;
  targetMol = cloneAndModifyMolecule(
      mol, false,
      false);  // this is false, false based on the python source code

  // Handle potential failure in cloning
  if (!targetMol) {
    // Suppress noisy warnings: return NaN silently
    delete targetMol;
    return std::numeric_limits<double>::quiet_NaN();  // Return NaN instead of
                                                      // crashing
  }

  double coreCount = 0.0;
  for (const auto &atom : targetMol->atoms()) {
    coreCount += getCoreCount(*atom);
  }
  if (averaged) {
    coreCount /= mol.getNumHeavyAtoms();
  }
  delete targetMol;
  return coreCount;
}

//  Etacorecount
std::vector<double> calculateEtaCoreCount(const ROMol &mol) {
  double coreCount = 0.0;
  for (const auto &atom : mol.atoms()) {
    coreCount += getCoreCount(*atom);
  }

  return {coreCount, coreCount / mol.getNumHeavyAtoms()};
}

std::vector<double> calculateEtaShapeIndex(const ROMol &mol, double alpha) {
  double shapeAlphaP = 0.0;
  double shapeAlphaY = 0.0;
  double shapeAlphaX = 0.0;

  for (const auto &atom : mol.atoms()) {
    switch (atom->getDegree()) {
      case 1:
        shapeAlphaP += getCoreCount(*atom);
        break;
      case 3:
        shapeAlphaY += getCoreCount(*atom);
        break;
      case 4:
        shapeAlphaX += getCoreCount(*atom);
        break;
      default:
        break;  // Ignore atoms with other degrees
    }
  }

  return {shapeAlphaP / alpha, shapeAlphaY / alpha, shapeAlphaX / alpha};
}

std::vector<double> computeEtaBetaDescriptors(const Atom &atom) {
  const PeriodicTable *tbl = PeriodicTable::getTable();

  double epsilon = getEtaEpsilon(atom);
  double betaSigma = 0.0;
  double betaNonSigma = 0.0;
  double betaDelta = 0.0;

  bool isAromatic = atom.getIsAromatic();
  bool inRing = isAtomInRing(atom);

  // Check if the atom satisfies delta condition
  if (!isAromatic && !inRing &&
      (tbl->getNouterElecs(atom.getAtomicNum()) - atom.getTotalValence() > 0)) {
    for (const auto &neighbor : atom.getOwningMol().atomNeighbors(&atom)) {
      if (neighbor->getIsAromatic()) {
        betaDelta = 0.5;
        break;
      }
    }
  }

  // Iterate through bonds
  for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
    const Atom *neighbor = getOtherAtom(*bond, atom);
    if (neighbor->getAtomicNum() != 1) {
      double neighborEpsilon = getEtaEpsilon(*neighbor);

      // Compute betaSigma (sigma contribution)
      betaSigma += (std::abs(neighborEpsilon - epsilon) <= 0.3) ? 0.5 : 0.75;

      // Compute betaNonSigma (non-sigma contribution)
      if (bond->getBondType() != Bond::SINGLE) {
        double bondFactor = (bond->getBondType() == Bond::TRIPLE) ? 2.0 : 1.0;
        double dEps = std::abs(epsilon - neighborEpsilon);
        double bondContribution =
            bond->getIsAromatic() ? 2.0 : (dEps > 0.3 ? 1.5 : 1.0);
        betaNonSigma += bondContribution * bondFactor;
      }
    }
  }

  return {betaSigma, betaNonSigma, betaDelta};
}

std::vector<double> calculateEtaVEMCount(const ROMol &mol) {
  double beta = 0.0;
  double beta_s = 0.0;
  double beta_ns = 0.0;
  double beta_ns_d = 0.0;

  int numAtoms = mol.getNumAtoms();

  for (const auto &atom : mol.atoms()) {
    std::vector<double> EBD = computeEtaBetaDescriptors(*atom);
    double betaDelta = EBD[2];
    double betaSigma = EBD[0] / 2.0;
    double betaNonSigma = EBD[1] / 2.0 + betaDelta;

    beta_s += betaSigma;
    beta_ns += betaNonSigma;
    beta_ns_d += betaDelta;
    beta += betaSigma + betaNonSigma;
  }

  // Calculate averaged values
  double avg_beta = beta / numAtoms;
  double avg_beta_s = beta_s / numAtoms;
  double avg_beta_ns = beta_ns / numAtoms;
  double avg_beta_ns_d = beta_ns_d / numAtoms;

  return {
      beta,          // ETA_beta
      avg_beta,      // AETA_beta
      beta_s,        // ETA_beta_s
      avg_beta_s,    // AETA_beta_s
      beta_ns,       // ETA_beta_ns
      avg_beta_ns,   // AETA_beta_ns
      beta_ns_d,     // ETA_beta_ns_d
      avg_beta_ns_d  // AETA_beta_ns_d
  };
}

// Function to Calculate ETA Composite Index
double calculateEtaCompositeIndex(const ROMol &mol, bool useReference,
                                  bool local, bool averaged) {
  // Fetch molecule (reference or input)

  const ROMol *targetMol = &mol;  // Default to the input molecule
  // ROMol* molWithHs = nullptr;

  if (useReference) {
    targetMol = cloneAndModifyMolecule(
        mol, false,
        false);  // this is false, false based on the python source code

    // Handle potential failure in cloning
    if (!targetMol) {
      // Suppress noisy warnings: return NaN silently
      delete targetMol;
      return std::numeric_limits<double>::quiet_NaN();  // Return NaN instead of
                                                        // crashing
    }
  }

  // Calculate distance matrix

  Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(*targetMol);
  int numAtoms = targetMol->getNumAtoms();

  // Define gamma values for each atom
  std::vector<double> gamma(numAtoms, 0.0);
  for (const auto &atom : targetMol->atoms()) {
    gamma[atom->getIdx()] = getEtaGamma(*atom);
  }

  // ETA calculation "triangle computation" ie  j = i + 1 trick to go faster
  double eta = 0.0;
  for (int i = 0; i < numAtoms; ++i) {
    for (int j = i + 1; j < numAtoms; ++j) {
      if (local && distanceMatrix(i, j) != 1.0) continue;
      if (!local && distanceMatrix(i, j) == 0.0) continue;
      eta += std::sqrt(gamma[i] * gamma[j] /
                       (distanceMatrix(i, j) * distanceMatrix(i, j)));
    }
  }

  // Averaged value if needed
  if (averaged) {
    eta /= numAtoms;
  }
  if (useReference) {
    delete targetMol;
  }
  return eta;
}

std::vector<double> calculateEtaCompositeIndices(const ROMol &mol) {
  std::vector<double> etaValues(8, 0.0);

  // Define options for different cases
  std::vector<std::tuple<bool, bool, bool>> options = {
      {false, false, false},  // ETA_eta
      {false, false, true},   // AETA_eta
      {false, true, false},   // ETA_eta_L
      {false, true, true},    // AETA_eta_L
      {true, false, false},   // ETA_eta_R
      {true, false, true},    // AETA_eta_R
      {true, true, false},    // ETA_eta_RL
      {true, true, true}      // AETA_eta_RL
  };

  for (size_t idx = 0; idx < options.size(); ++idx) {
    bool useReference = std::get<0>(options[idx]);
    bool local = std::get<1>(options[idx]);
    bool averaged = std::get<2>(options[idx]);

    // Fetch molecule (reference or input)
    const ROMol *targetMol = &mol;  // Default to the input molecule
    if (useReference) {
      targetMol = cloneAndModifyMolecule(
          mol, false, false);  // Clone with modifications (why not "True")

      // Handle potential failure in cloning
      if (!targetMol) {
        // Suppress noisy warnings: return NaN silently
        delete targetMol;
        return std::vector<double>(
            8, std::numeric_limits<double>::quiet_NaN());  // Return vector with
                                                           // NaN
      }
    }

    // Calculate distance matrix
    Eigen::MatrixXd distanceMatrix = calculateDistanceMatrix(*targetMol);
    int numAtoms = targetMol->getNumAtoms();

    // Define gamma values for each atom
    std::vector<double> gamma(numAtoms, 0.0);
    for (const auto &atom : targetMol->atoms()) {
      gamma[atom->getIdx()] = getEtaGamma(*atom);
    }

    // ETA calculation using the optimized "triangle computation"
    double eta = 0.0;
    for (int i = 0; i < numAtoms; ++i) {
      for (int j = i + 1; j < numAtoms; ++j) {
        if (local && distanceMatrix(i, j) != 1.0) continue;
        if (!local && distanceMatrix(i, j) == 0.0) continue;

        eta += std::sqrt(gamma[i] * gamma[j] /
                         (distanceMatrix(i, j) * distanceMatrix(i, j)));
      }
    }

    // Averaged value if required
    if (averaged) {
      eta /= numAtoms;
    }

    etaValues[idx] = eta;

    if (useReference) {
      delete targetMol;  // Clean up dynamically allocated molecule
    }
  }

  return etaValues;
}

std::vector<double> calculateEtaFunctionalityIndices(const ROMol &mol) {
  std::vector<double> etaFunctionalityValues(4, 0.0);

  // Define options for different cases
  std::vector<std::tuple<bool, bool>> options = {
      {false, false},  // ETA_eta_F
      {false, true},   // AETA_eta_F
      {true, false},   // ETA_eta_FL
      {true, true}     // AETA_eta_FL
  };

  for (size_t idx = 0; idx < options.size(); ++idx) {
    bool local = std::get<0>(options[idx]);
    bool averaged = std::get<1>(options[idx]);

    // Calculate eta without reference and with reference
    double eta = calculateEtaCompositeIndex(mol, false, local, false);
    double etaRef = calculateEtaCompositeIndex(mol, true, local, false);

    // Compute functionality index
    double etaF = etaRef - eta;

    // Apply averaging if needed
    if (averaged) {
      etaF /= mol.getNumAtoms();
    }

    etaFunctionalityValues[idx] = etaF;
  }

  return etaFunctionalityValues;
}

std::vector<double> calculateEtaBranchingIndices(const ROMol &mol) {
  std::vector<double> etaBranchingValues(4, 0.0);
  int atomCount = mol.getNumAtoms();

  if (atomCount <= 1) {
    return etaBranchingValues;  // Return zeros if the molecule has only one
                                // atom
  }

  // Calculate non-local branching term
  double eta_NL =
      (atomCount == 2) ? 1.0 : (std::sqrt(2.0) + 0.5 * (atomCount - 3));

  double eta_RL = calculateEtaCompositeIndex(mol, true, true, false);

  // Calculate ring count once
  double ringCount = Descriptors::calcNumRings(mol);

  // Define options for different cases
  std::vector<std::tuple<bool, bool>> options = {
      {false, false},  // ETA_eta_B
      {false, true},   // AETA_eta_B
      {true, false},   // ETA_eta_BR
      {true, true}     // AETA_eta_BR
  };

  for (size_t idx = 0; idx < options.size(); ++idx) {
    bool useRingCount = std::get<0>(options[idx]);
    bool averaged = std::get<1>(options[idx]);

    double etaBranch = eta_NL - eta_RL;
    if (useRingCount) {
      etaBranch += 0.086 * ringCount;
    }

    if (averaged) {
      etaBranch /= atomCount;
    }

    etaBranchingValues[idx] = etaBranch;
  }

  return etaBranchingValues;
}

std::vector<double> calculateEtaDeltaBetaAll(const ROMol &mol, double beta_ns,
                                             double beta_s) {
  std::vector<double> results(2, 0.0);

  // Calculate ETA_beta (delta_beta without averaging)
  double delta_beta = beta_ns - beta_s;
  results[0] = delta_beta;  // ETA_beta

  // Calculate AETA_beta (averaged delta_beta)
  if (mol.getNumAtoms() > 0) {
    results[1] = delta_beta / mol.getNumAtoms();  // AETA_beta
  } else {
    results[1] = 0.0;  // Handle division by zero case
  }

  return results;
}

// Function to calculate ETA Psi
double calculateEtaPsi(const ROMol &mol, double alpha, double epsilon) {
  if (epsilon == 0) {
    throw std::runtime_error("Epsilon cannot be zero");
  }
  return alpha / (mol.getNumAtoms() * epsilon);
}

std::vector<double> calculateEtaDeltaPsiAll(const ROMol &, double psi) {
  std::vector<double> results(2, 0.0);

  // Constants for reference values
  const double L = 0.714;
  const double R = psi;

  // Calculate ETA_dPsi_A (L - R)
  results[0] = std::max(L - R, 0.0);  // ETA_dPsi_A

  // Calculate ETA_dPsi_B (R - L)
  results[1] = std::max(R - L, 0.0);  // ETA_dPsi_B

  return results;
}

std::vector<double> calculateEtaDeltaAlpha(const ROMol &mol, double alpha,
                                           double alpha_R) {
  std::vector<double> etaDeltaAlphaValues(2, 0.0);
  int numAtoms = mol.getNumAtoms();

  if (numAtoms <= 0) {
    return etaDeltaAlphaValues;  // Return zeros if the molecule has no atoms
  }

  // Calculate dAlpha_A and dAlpha_B
  etaDeltaAlphaValues[0] =
      std::max((alpha - alpha_R) / numAtoms, 0.0);  // ETA_dAlpha_A
  etaDeltaAlphaValues[1] =
      std::max((alpha_R - alpha) / numAtoms, 0.0);  // ETA_dAlpha_B

  return etaDeltaAlphaValues;
}

std::vector<double> calculateEtaEpsilonAll(const ROMol &mol) {
  std::vector<double> etaEps(5, 0.0);

  // Step 1: Calculate ETA epsilon for type 2 (non-hydrogen atoms)
  double epsilon_sum_type2 = 0.0;
  for (const auto &atom : mol.atoms()) {
    epsilon_sum_type2 += getEtaEpsilon(*atom);
  }
  etaEps[1] = epsilon_sum_type2 / mol.getNumAtoms();  // ETA_epsilon_2

  // Step 2: Add hydrogens for types 1 and 5 calculations
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  // Step 3: Calculate ETA epsilon for type 1 and type 5 (hydrogen inclusion)
  double epsilon_sum_type1 = 0.0, epsilon_sum_type5 = 0.0;
  int count_type1 = 0, count_type5 = 0;
  // Debug counters
  int debug_num_h_in_hmol = 0;
  int debug_excluded_ch = 0;

  for (const auto &atom : hmol->atoms()) {
    // type 1
    double EtaEps = getEtaEpsilon(*atom);
    epsilon_sum_type1 += EtaEps;
    count_type1++;
    // type 5
    bool hasCarbonNeighbor = false;
    for (const auto &neighbor : hmol->atomNeighbors(atom)) {
      if (neighbor->getAtomicNum() == 6) {
        hasCarbonNeighbor = true;
        break;
      }
    }
    if (atom->getAtomicNum() == 1) {
      // count hydrogens in hmol
      debug_num_h_in_hmol++;
    }
    if (atom->getAtomicNum() != 1 || !hasCarbonNeighbor) {
      epsilon_sum_type5 += EtaEps;
      count_type5++;
    } else {
      // Hydrogen attached to carbon is excluded from epsilon_5
      if (atom->getAtomicNum() == 1 && hasCarbonNeighbor) {
        debug_excluded_ch++;
      }
    }
  }
  etaEps[0] = epsilon_sum_type1 / count_type1;  // ETA_epsilon_1
  etaEps[4] =
      count_type5 > 0 ? epsilon_sum_type5 / count_type5 : 0.0;  // ETA_epsilon_5

  // Optional debug logging (enabled when OSMO_ETA_DEBUG is set)
  const char *eta_dbg = std::getenv("OSMO_ETA_DEBUG");
  if (eta_dbg) {
    std::cerr << "[ETA_DEBUG] hmol_atoms=" << hmol->getNumAtoms()
              << " hydrogens_in_hmol=" << debug_num_h_in_hmol
              << " excluded_CH_for_eps5=" << debug_excluded_ch
              << " epsilon_sum_type1=" << epsilon_sum_type1
              << " epsilon_sum_type2=" << epsilon_sum_type2
              << " epsilon_sum_type5=" << epsilon_sum_type5
              << " count_type1=" << count_type1
              << " count_type5=" << count_type5 << std::endl;
  }

  // Step 4: Calculate ETA epsilon for types 3 and 4 (modified molecules)
  const ROMol *targetMol3 = modifyMolecule(mol, true, false);
  const ROMol *targetMol4 = modifyMolecule(mol, true, true);

  double epsilon_sum_type3 = 0.0;
  for (const auto &atom : targetMol3->atoms()) {
    epsilon_sum_type3 += getEtaEpsilon(*atom);
  }
  etaEps[2] = epsilon_sum_type3 / targetMol3->getNumAtoms();  // ETA_epsilon_3

  double epsilon_sum_type4 = 0.0;
  for (const auto &atom : targetMol4->atoms()) {
    epsilon_sum_type4 += getEtaEpsilon(*atom);
  }
  etaEps[3] = epsilon_sum_type4 / targetMol4->getNumAtoms();  // ETA_epsilon_4

  // Clean up dynamically allocated molecules
  delete targetMol3;
  delete targetMol4;

  return etaEps;
}

std::vector<double> calcExtendedTopochemicalAtom(const ROMol &mol) {
  std::vector<double> results;

  std::unique_ptr<RWMol> kekulizedMol(new RWMol(mol));
  MolOps::Kekulize(*kekulizedMol, false);

  //  Eta alpha  Descriptors correct  "EtaCoreCount" : alpha & shape can be
  //  group in one function
  std::vector<double> myalphas = calculateEtaCoreCount(*kekulizedMol);
  double alpha = myalphas[0];

  results.push_back(myalphas[0]);  // Eta alpha
  results.push_back(myalphas[1]);  // average Eta alpha
  // Shape Descriptors correct     "EtaShapeIndex",

  // working for all atoms types
  std::vector<double> ESP = calculateEtaShapeIndex(*kekulizedMol, myalphas[0]);

  results.push_back(ESP[0]);  // ETA_shape_p
  results.push_back(ESP[1]);  // ETA_shape_y
  results.push_back(ESP[2]);  // ETA_shape_x

  // Beta Descriptors ie EtaVEMCount correct : can be group in one function for
  // optimizatoin
  std::vector<double> EtaVEM = calculateEtaVEMCount(*kekulizedMol);

  results.push_back(EtaVEM[0]);  // ETA_beta
  results.push_back(EtaVEM[1]);  // AETA_beta
  double etabetas = EtaVEM[2];
  results.push_back(EtaVEM[2]);  // ETA_beta_s
  results.push_back(EtaVEM[3]);  // AETA_beta_s
  double etabetans = EtaVEM[4];
  results.push_back(EtaVEM[4]);  // ETA_beta_ns
  results.push_back(EtaVEM[5]);  // AETA_beta_ns
  results.push_back(EtaVEM[6]);  // ETA_beta_ns_d
  results.push_back(EtaVEM[7]);  // AETA_beta_ns_d

  // ETA Descriptors   ==  "EtaCompositeIndex"
  std::vector<double> ECI = calculateEtaCompositeIndices(*kekulizedMol);
  results.push_back(ECI[0]);  // ETA_eta
  results.push_back(ECI[1]);  // AETA_eta
  results.push_back(ECI[2]);  // ETA_eta_L
  results.push_back(ECI[3]);  // AETA_eta_L
  results.push_back(ECI[4]);  // ETA_eta_R
  results.push_back(ECI[5]);  // AETA_eta_R
  results.push_back(ECI[6]);  // ETA_eta_RL
  results.push_back(ECI[7]);  // AETA_eta_RL

  // Functionality and Branching EtaFunctionalityIndex not working for
  // heteroatom molecule...
  std::vector<double> EFI = calculateEtaFunctionalityIndices(*kekulizedMol);
  results.push_back(EFI[0]);  // ETA_eta_F
  results.push_back(EFI[1]);  // AETA_eta_F
  results.push_back(EFI[2]);  // ETA_eta_FL
  results.push_back(EFI[3]);  // AETA_eta_FL

  // EtaBranchingIndex  working
  std::vector<double> EBI = calculateEtaBranchingIndices(*kekulizedMol);
  results.push_back(EBI[0]);  // ETA_eta_B
  results.push_back(EBI[1]);  // AETA_eta_B
  results.push_back(EBI[2]);  // ETA_eta_BR
  results.push_back(EBI[3]);  // AETA_eta_BR

  //"EtaDeltaAlpha" :  dAlpha_A, dAlpha_B  correct
  double alpha_R = calculateEtaCoreCountRef(*kekulizedMol, false);
  std::vector<double> EDA =
      calculateEtaDeltaAlpha(*kekulizedMol, alpha, alpha_R);
  results.push_back(EDA[0]);
  results.push_back(EDA[1]);

  // "EtaEpsilon":
  std::vector<double> EtaEps = calculateEtaEpsilonAll(*kekulizedMol);
  results.push_back(EtaEps[0]);
  results.push_back(EtaEps[1]);
  results.push_back(EtaEps[2]);
  results.push_back(EtaEps[3]);
  results.push_back(EtaEps[4]);

  // "EtaDeltaEpsilon", A,B,C,D use EtaEps diff cases
  results.push_back(EtaEps[1 - 1] - EtaEps[3 - 1]);
  results.push_back(EtaEps[1 - 1] - EtaEps[4 - 1]);
  results.push_back(EtaEps[3 - 1] - EtaEps[4 - 1]);
  results.push_back(EtaEps[2 - 1] - EtaEps[5 - 1]);

  //"EtaDeltaBeta":
  std::vector<double> EDBA =
      calculateEtaDeltaBetaAll(*kekulizedMol, etabetans, etabetas);
  results.push_back(EDBA[0]);
  results.push_back(EDBA[1]);

  // EtaPsi:
  double EtaPsi = calculateEtaPsi(*kekulizedMol, alpha, EtaEps[1]);
  results.push_back(EtaPsi);  // ETA_psi_1

  // EtaDeltaPsi:
  std::vector<double> EDPA = calculateEtaDeltaPsiAll(*kekulizedMol, EtaPsi);
  results.push_back(EDPA[0]);
  results.push_back(EDPA[1]);

  return results;
}

      
//// new version faster but not correct!!!

std::vector<std::vector<int>> calculateTopologicalMatrix(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();
  std::vector<std::vector<int>> distanceMatrix(numAtoms,
                                               std::vector<int>(numAtoms, -1));

  for (int i = 0; i < numAtoms; ++i) {
    distanceMatrix[i][i] = 0;  // Distance to self is always zero
    for (int j = i + 1; j < numAtoms; ++j) {
      auto path = MolOps::getShortestPath(mol, i, j);
      distanceMatrix[i][j] = path.size();
      distanceMatrix[j][i] = distanceMatrix[i][j];  // Symmetric assignment
    }
  }
  return distanceMatrix;
}

// Function to calculate atomic properties
void calculateAtomicDescriptors(const ROMol &mol, std::vector<double> &alpha,
                                std::vector<double> &alphaR,
                                std::vector<double> &epsilon, double &Alpha_P,
                                double &Alpha_Y, double &Alpha_X) {
  int numAtoms = mol.getNumAtoms();
  alpha.resize(numAtoms, 0.0);
  alphaR.resize(numAtoms, 0.5);
  epsilon.resize(numAtoms, 0.3);
  const PeriodicTable *tbl = PeriodicTable::getTable();
  Alpha_P = 0.0;
  Alpha_Y = 0.0;
  Alpha_X = 0.0;

  for (int i = 0; i < numAtoms; ++i) {
    const auto atom = mol.getAtomWithIdx(i);
    int atomicNum = atom->getAtomicNum();

    if (atomicNum != 1) {
      int Z = atomicNum;
      int Zv = tbl->getNouterElecs(Z);
      int period = GetPrincipalQuantumNumber(Z);

      alpha[i] = ((Z - Zv) / static_cast<double>(Zv)) * (1.0 / (period - 1));
      epsilon[i] = -alpha[i] + 0.3 * Zv;

      int nonHNeighbors = 0;
      for (const auto &bond : mol.atomBonds(atom)) {
        const auto neighbor = bond->getOtherAtom(atom);
        if (neighbor->getAtomicNum() != 1) {
          nonHNeighbors++;
        }
      }

      if (nonHNeighbors == 1) {
        Alpha_P += alpha[i];
      } else if (nonHNeighbors == 3) {
        Alpha_Y += alpha[i];
      } else if (nonHNeighbors == 4) {
        Alpha_X += alpha[i];
      }
    }
  }
}

void computeVEMContributions(
    const ROMol &mol, const std::vector<double> &epsilon,
    std::vector<double> &betaS, std::vector<double> &betaNS,
    std::vector<double> &betaNSd, std::vector<double> &beta,
    std::vector<double> &gamma, std::vector<double> &gammaRef,
    const std::vector<double> &alpha, const std::vector<double> &alphaR) {
  int numAtoms = mol.getNumAtoms();
  betaS.resize(numAtoms, 0.0);
  betaNS.resize(numAtoms, 0.0);
  betaNSd.resize(numAtoms, 0.0);
  beta.resize(numAtoms, 0.0);
  gamma.resize(numAtoms, 0.0);
  gammaRef.resize(numAtoms, 0.0);

  const PeriodicTable *tbl = PeriodicTable::getTable();

  for (int i = 0; i < numAtoms; ++i) {
    const auto atom = mol.getAtomWithIdx(i);
    betaS[i] = 0.0;
    betaNS[i] = 0.0;
    betaNSd[i] = 0.0;
    double betaR = 0.0;
    bool aromatic = false;
    double nonconjugatedSum = 0.0;
    bool conjugated = false;
    bool hasConjugatedSingleBond = false;
    bool hasDoubleTripleBond = false;
    int conjugatedMultiplier = 1;
    bool delta = false;
    for (const auto &bond : mol.atomBonds(atom)) {
      const auto neighbor = bond->getOtherAtom(atom);
      if (neighbor->getAtomicNum() == 1) continue;
      betaR += 0.5;
      int j = neighbor->getIdx();
      double epsilonDiff = std::abs(epsilon[i] - epsilon[j]);
      if (epsilonDiff <= 0.3) {
        betaS[i] += 0.5 / 2;  // error in code have to divide by 2??? look like
                              // count twice ???
      } else {
        betaS[i] += 0.75 / 2;  // error in code have to divide by 2??? look like
                               // count twice ???
      }

      if (bond->getIsAromatic()) {
        aromatic = true;
      } else if (bond->getBondType() == Bond::SINGLE) {
        if (!conjugated) {
          if (neighbor->getAtomicNum() == 6) {  // Check if neighbor is carbon
            for (const auto &nbrBond : mol.atomBonds(neighbor)) {
              if (nbrBond->getBondType() == Bond::DOUBLE ||
                  nbrBond->getIsAromatic()) {
                hasConjugatedSingleBond = true;
                break;
              } else if (nbrBond->getBondType() == Bond::TRIPLE) {
                conjugatedMultiplier = 2;
                hasConjugatedSingleBond = true;
                break;
              }
            }
          } else if (tbl->getNouterElecs(neighbor->getAtomicNum()) -
                         neighbor->getFormalCharge() -
                         mol.getAtomDegree(neighbor) >=
                     2) {
            hasConjugatedSingleBond = true;
          }
          if (hasDoubleTripleBond && hasConjugatedSingleBond) conjugated = true;
        }
      } else if (bond->getBondType() == Bond::DOUBLE ||
                 bond->getBondType() == Bond::TRIPLE) {
        int multiplier = (bond->getBondType() == Bond::DOUBLE) ? 1 : 2;
        conjugatedMultiplier = multiplier;

        hasDoubleTripleBond = true;
        if (!conjugated) {
          if (hasConjugatedSingleBond) {
            conjugated = true;
          } else if (neighbor->getAtomicNum() == 6) {
            for (const auto &bond2 : mol.atomBonds(neighbor)) {
              if (bond2->getBondType() == Bond::SINGLE) {
                const auto atom3 = bond2->getOtherAtom(neighbor);
                if (atom3->getAtomicNum() == 6) {
                  for (const auto &bond3 : mol.atomBonds(atom3)) {
                    if (bond3->getBondType() == Bond::DOUBLE ||
                        bond3->getIsAromatic()) {
                      hasConjugatedSingleBond = true;
                      break;
                    } else if (bond3->getBondType() == Bond::TRIPLE) {
                      conjugatedMultiplier = 2;
                      hasConjugatedSingleBond = true;
                      break;
                    }
                  }
                } else if (tbl->getNouterElecs(atom3->getAtomicNum()) -
                               atom3->getFormalCharge() -
                               mol.getAtomDegree(atom3) >=
                           2) {
                  hasConjugatedSingleBond = true;
                }
                if (hasConjugatedSingleBond) {
                  conjugated = true;
                  break;
                }
              }
            }
          }
        }
        double g = (epsilonDiff <= 0.3) ? 1.0 : 1.5;
        nonconjugatedSum += g * multiplier;
      }
      if (tbl->getNouterElecs(atom->getAtomicNum()) - atom->getFormalCharge() -
                  mol.getAtomDegree(atom) >=
              2 &&
          !atom->getIsAromatic() && !isAtomInRing(*atom) &&
          neighbor->getIsAromatic() && isAtomInRing(*neighbor)) {
        delta = true;
      }
    }

    betaNS[i] += (aromatic ? 2.0 : 0.0) +
                 (conjugated ? 1.5 * conjugatedMultiplier : nonconjugatedSum) +
                 (delta ? 0.5 : 0.0);
    betaNSd[i] = (delta ? 0.5 : 0.0);
    beta[i] = betaS[i] + betaNS[i];
    // Compute gamma values
    gamma[i] = alpha[i] / beta[i];
    gammaRef[i] = alphaR[i] / betaR;
  }
}

// Main function to compute ETA descriptors
std::vector<double> calcETADescriptors(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::unique_ptr<RWMol> kekulizedMol(new RWMol(mol));
  MolOps::Kekulize(*kekulizedMol, false);

  // Step 1: Calculate atomic descriptors
  std::vector<double> alpha, alphaR, epsilon;
  double Alpha_P, Alpha_Y, Alpha_X;
  calculateAtomicDescriptors(*kekulizedMol, alpha, alphaR, epsilon, Alpha_P,
                             Alpha_Y, Alpha_X);

  // Step 2: Compute VEM contributions fixed by /2 the values why mistary ????
  std::vector<double> betaS, betaNS, betaNSd, beta, gamma, gammaRef;
  computeVEMContributions(*kekulizedMol, epsilon, betaS, betaNS, betaNSd, beta,
                          gamma, gammaRef, alpha, alphaR);

  // Step 3: Calculate topological matrix (distance matrix)
  auto pathLengths = calculateTopologicalMatrix(mol);

  // Step 4: Count rings
  int maxRings = MolOps::findSSSR(mol);

  // Step 5: Compute descriptor sums
  double alphaSum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
  double alphaRSum = std::accumulate(alphaR.begin(), alphaR.end(), 0.0);
  double epsilonSum = std::accumulate(epsilon.begin(), epsilon.end(), 0.0);
  double betaSSum = std::accumulate(betaS.begin(), betaS.end(), 0.0);
  double betaNSSum = std::accumulate(betaNS.begin(), betaNS.end(), 0.0);
  double betaNSdSum = std::accumulate(betaNSd.begin(), betaNSd.end(), 0.0);

  // Step 6: Compute various ETA values based on atomic and bond properties
  double etaSum = 0.0, etaRSum = 0.0, etaLSum = 0.0, etaRLSum = 0.0;

  for (int i = 0; i < numAtoms; ++i) {
    for (int j = i + 1; j < numAtoms; ++j) {
      if (pathLengths[i][j] > 0) {
        etaSum += std::sqrt(gamma[i] * gamma[j] /
                            (pathLengths[i][j] * pathLengths[i][j]));
        etaRSum += std::sqrt(gammaRef[i] * gammaRef[j] /
                             (pathLengths[i][j] * pathLengths[i][j]));

        if (pathLengths[i][j] == 1) {
          etaLSum += std::sqrt(gamma[i] * gamma[j]);
          etaRLSum += std::sqrt(gammaRef[i] * gammaRef[j]);
        }
      }
    }
  }

  // Step 7: Compute ETA indices

  int N = numAtoms;  // Total atoms count
  int Nv = numAtoms -
           std::count_if(mol.atoms().begin(), mol.atoms().end(),
                         [](const Atom *a) { return a->getAtomicNum() == 1; });
  int Nr =
      5 * N;  // Placeholder for number of reference hydrogens, adjust as needed
  int Nss = N;  // Total number of saturated atoms
  int Nxh = 0;  // Count heteroatom-bound hydrogens (modify logic if needed)
  double epsilonEHSum = epsilonSum;
  double epsilonRSum = epsilonSum;
  double epsilonSSSum = epsilonSum;
  double epsilonXHSum = 0.0;  // Placeholder for heteroatom-bonded H sum

  double PsiTemp = 0.71429;  // Fixed constant

  // EtaCoreCount
  double ETA_Alpha = alphaSum;
  double ETA_AlphaP = alphaSum / Nv;

  // EtaShapeIndex
  double ETA_Shape_P = Alpha_P / ETA_Alpha;
  double ETA_Shape_Y = Alpha_Y / ETA_Alpha;
  double ETA_Shape_X = Alpha_X / ETA_Alpha;

  // VEMCount
  double ETA_Beta = betaSSum + betaNSSum;
  double ETA_BetaP = ETA_Beta / Nv;
  double ETA_Beta_s = betaSSum;
  double ETA_BetaP_s = betaSSum / Nv;
  double ETA_Beta_ns = betaNSSum;
  double ETA_BetaP_ns = betaNSSum / Nv;
  double ETA_Beta_ns_d = betaNSdSum;
  double ETA_BetaP_ns_d = ETA_Beta_ns_d / Nv;

  //  "EtaCompositeIndex"
  double ETA_Eta = etaSum;
  double ETA_EtaP = etaSum / Nv;
  double ETA_Eta_L = etaLSum;
  double ETA_EtaP_L = etaLSum / Nv;
  double ETA_Eta_R = etaRSum;
  double ETA_EtaP_R = etaRSum / Nv;
  double ETA_Eta_R_L = etaRLSum;
  double ETA_EtaP_R_L = etaRLSum / Nv;

  // EtaFunctionalityIndex
  double ETA_Eta_F = etaRSum - etaSum;
  double ETA_EtaP_F = ETA_Eta_F / Nv;
  double ETA_Eta_F_L = etaRLSum - etaLSum;
  double ETA_EtaP_F_L = ETA_Eta_F_L / Nv;

  // EtaBranchingIndex
  double ETA_Eta_B = (Nv > 3 ? (1.414 + (Nv - 3) * 0.5) - etaRLSum : 0);
  double ETA_EtaP_B = ETA_Eta_B / Nv;
  double ETA_Eta_B_RC = ETA_Eta_B + 0.086 * maxRings;
  double ETA_EtaP_B_RC = ETA_Eta_B_RC / Nv;

  //"EtaDeltaAlpha" :  dAlpha_A, dAlpha_B
  double ETA_dAlpha_A = std::max((alphaSum - alphaRSum) / Nv, 0.0);
  double ETA_dAlpha_B = std::max((alphaRSum - alphaSum) / Nv, 0.0);

  // "EtaEpsilon":
  double ETA_Epsilon_1 = epsilonSum / N;
  double ETA_Epsilon_2 = epsilonEHSum / Nv;
  double ETA_Epsilon_3 = epsilonRSum / Nr;
  double ETA_Epsilon_4 = epsilonSSSum / Nss;
  double ETA_Epsilon_5 = (epsilonEHSum + epsilonXHSum) / (Nv + Nxh);

  // "EtaDeltaEpsilon", A,B,C,D use EtaEps diff cases
  double ETA_dEpsilon_A = ETA_Epsilon_1 - ETA_Epsilon_3;
  double ETA_dEpsilon_B = ETA_Epsilon_1 - ETA_Epsilon_4;
  double ETA_dEpsilon_C = ETA_Epsilon_3 - ETA_Epsilon_4;
  double ETA_dEpsilon_D = ETA_Epsilon_2 - ETA_Epsilon_5;

  //"EtaDeltaBeta" 2 Values
  double ETA_dBeta = betaNSSum - betaSSum;
  double ETA_dBetaP = ETA_dBeta / Nv;

  // EtaPsi
  double ETA_Psi_1 = alphaSum / epsilonEHSum;

  // EtaDeltaPsi
  double ETA_dPsi_A = std::max(PsiTemp - ETA_Psi_1, 0.0);
  double ETA_dPsi_B = std::max(ETA_Psi_1 - PsiTemp, 0.0);
  // Step 8: Return all calculated values
  return {ETA_Alpha,      ETA_AlphaP,     ETA_Shape_P,    ETA_Shape_Y,
          ETA_Shape_X,    ETA_Beta,       ETA_BetaP,      ETA_Beta_s,
          ETA_BetaP_s,    ETA_Beta_ns,    ETA_BetaP_ns,   ETA_Beta_ns_d,
          ETA_BetaP_ns_d, ETA_Eta,        ETA_EtaP,       ETA_Eta_L,
          ETA_EtaP_L,     ETA_Eta_R,      ETA_EtaP_R,     ETA_Eta_R_L,
          ETA_EtaP_R_L,   ETA_Eta_F,      ETA_EtaP_F,     ETA_Eta_F_L,
          ETA_EtaP_F_L,   ETA_Eta_B,      ETA_EtaP_B,     ETA_Eta_B_RC,
          ETA_EtaP_B_RC,  ETA_dAlpha_A,   ETA_dAlpha_B,   ETA_Epsilon_1,
          ETA_Epsilon_2,  ETA_Epsilon_3,  ETA_Epsilon_4,  ETA_Epsilon_5,
          ETA_dEpsilon_A, ETA_dEpsilon_B, ETA_dEpsilon_C, ETA_dEpsilon_D,
          ETA_dBeta,      ETA_dBetaP,     ETA_Psi_1,      ETA_dPsi_A,
          ETA_dPsi_B};
}

// start ringcount
std::vector<std::vector<int>> GetSSSR(const ROMol &mol) {
  ROMol mol_copy(mol);  // Create a non-const copy of the molecule
  std::vector<std::vector<int>> sssr;
  MolOps::symmetrizeSSSR(mol_copy, sssr);
  return sssr;
}

// Get all rings in the molecule
std::vector<std::set<int>> GetRings(const ROMol &,
                                    std::vector<std::vector<int>> sssr) {
  std::vector<std::set<int>> rings;
  for (const auto &ring : sssr) {
    rings.emplace_back(ring.begin(), ring.end());
  }
  return rings;
}

// Type aliases for readability
using Ring = std::set<int>;
using Rings = std::vector<Ring>;

// Helper function to find connected components in an undirected graph
std::vector<std::vector<int>> findConnectedComponents(
    const std::unordered_map<int, std::unordered_set<int>> &graph) {
  std::vector<std::vector<int>> components;
  std::unordered_set<int> visited;

  for (const auto &[node, _] : graph) {
    if (visited.count(node) == 0) {
      // Perform BFS to find all nodes in the current component
      std::vector<int> component;
      std::queue<int> toVisit;
      toVisit.push(node);
      visited.insert(node);

      while (!toVisit.empty()) {
        int current = toVisit.front();
        toVisit.pop();
        component.push_back(current);

        for (int neighbor : graph.at(current)) {
          if (visited.count(neighbor) == 0) {
            visited.insert(neighbor);
            toVisit.push(neighbor);
          }
        }
      }

      components.push_back(component);
    }
  }

  return components;
}

// Function to calculate fused rings
std::vector<std::vector<int>> calcFusedRings(const Rings &rings) {
  if (rings.size() < 2) {
    return {};
  }

  // Graph adjacency list representation
  std::unordered_map<int, std::unordered_set<int>> graph;

  size_t numRings = rings.size();
  for (size_t i = 0; i < numRings; ++i) {
    for (size_t j = i + 1; j < numRings; ++j) {
      std::vector<int> intersection;
      std::set_intersection(rings[i].begin(), rings[i].end(), rings[j].begin(),
                            rings[j].end(), std::back_inserter(intersection));

      if (intersection.size() >= 2) {
        graph[i].insert(j);
        graph[j].insert(i);
      }
    }
  }

  // Find connected components in the graph
  std::vector<std::vector<int>> components = findConnectedComponents(graph);

  // Map connected ring indices to their atom sets
  std::vector<std::vector<int>> fusedRings;
  for (const auto &component : components) {
    std::unordered_set<int> fusedRingAtoms;
    for (int ringIdx : component) {
      fusedRingAtoms.insert(rings[ringIdx].begin(), rings[ringIdx].end());
    }

    fusedRings.emplace_back(fusedRingAtoms.begin(), fusedRingAtoms.end());
  }

  return fusedRings;
}

// Calculate ring count descriptors
std::vector<int> calcRingDescriptors(const ROMol &mol) {
  std::vector<int> descriptors(138, 0);  // Placeholder for descriptor vector

  std::vector<std::vector<int>> sssr = GetSSSR(mol);

  // part one make rings inventory

  descriptors[0] = sssr.size();
  // Iterate over all rings and count descriptors based on size, aromaticity,
  // fused
  for (const auto &ring : sssr) {
    size_t ring_size = ring.size();

    bool is_aromatic = true;

    bool has_hetero = false;

    for (int atom_idx : ring) {
      const Atom *atom = mol.getAtomWithIdx(atom_idx);
      if (!atom->getIsAromatic()) {
        is_aromatic = false;
      }
      if (atom->getAtomicNum() != 6) {
        has_hetero = true;
      }
    }

    // Update ring counts for the respective features
    if (ring_size >= 3 && ring_size <= 12) {
      descriptors[ring_size - 3 + 1]++;  // n3Ring to n12Ring
      if (has_hetero) {
        descriptors[12]++;  // nHRing sum
        descriptors[12 + (ring_size - 3) + 1]++;
      }

      if (is_aromatic) {
        descriptors[24]++;  // naRing

        descriptors[24 + (ring_size - 3) + 1]++;  // n3aRing to n12aRing
      }

      if (is_aromatic && has_hetero) {
        descriptors[36]++;  // naHRing

        descriptors[36 + (ring_size - 3) + 1]++;  // n3aHRing to n12aHRing
      }

      if (!is_aromatic) {
        descriptors[48]++;                        // nARing
        descriptors[48 + (ring_size - 3) + 1]++;  // n3ARing to n12ARing
      }

      if (!is_aromatic && has_hetero) {
        descriptors[60]++;                        // nAHRing
        descriptors[60 + (ring_size - 3) + 1]++;  // n3AHRing to n12AHRing
      }
    }

    if (ring_size > 12) {
      // greater
      if (has_hetero) {
        descriptors[12]++;  // nHRing sum
        descriptors[23]++;  // nG12 HRing
      }
      if (is_aromatic) {
        descriptors[24]++;  // naRing  sum
        descriptors[35]++;  // nG12 aRing
      }

      if (is_aromatic && has_hetero) {
        descriptors[36]++;  // naHRing  sum
        descriptors[47]++;  // nG12 aHRing
      }

      if (!is_aromatic) {
        descriptors[48]++;  // nARing  sum
        descriptors[59]++;  // nG12 ARing
      }

      if (!is_aromatic && has_hetero) {
        descriptors[60]++;  // nAHRing  sum
        descriptors[71]++;  // nG12AHRing
      }
    }
  }

  // part two make fused ring inventory
  // maybe we don't need the rings step...
  auto rings = GetRings(mol, sssr);
  auto fusedRings = calcFusedRings(rings);

  descriptors[72] = fusedRings.size();

  for (const auto &fusedring : fusedRings) {
    size_t fused_ring_size = fusedring.size();

    bool is_aromatic = true;

    bool has_hetero = false;
    // we could also break it...
    for (int atom_idx : fusedring) {
      const Atom *atom = mol.getAtomWithIdx(atom_idx);
      if (!atom->getIsAromatic()) {
        is_aromatic = false;
      }
      if (atom->getAtomicNum() != 6) {
        has_hetero = true;
      }
    }

    // Update ring counts for the respective features
    if (fused_ring_size >= 4 && fused_ring_size <= 12) {
      descriptors[72 + fused_ring_size - 4 + 1]++;  // n4FRing to n12FRing
      if (has_hetero) {
        descriptors[83]++;                              // nFHRing sum
        descriptors[83 + (fused_ring_size - 4) + 1]++;  // n4FHRing to n12FHRing
      }
      if (is_aromatic) {
        descriptors[94]++;                              // nFHRing sum
        descriptors[94 + (fused_ring_size - 4) + 1]++;  // n4FaRing to n12FaRing
      }

      if (is_aromatic && has_hetero) {
        descriptors[105]++;  // naHRing

        descriptors[105 + (fused_ring_size - 4) +
                    1]++;  // n3aHRing to n12aHRing
      }

      if (!is_aromatic) {
        descriptors[116]++;                              // nARing
        descriptors[116 + (fused_ring_size - 4) + 1]++;  // n3ARing to n12ARing
      }

      if (!is_aromatic && has_hetero) {
        descriptors[127]++;  // nAHRing
        descriptors[127 + (fused_ring_size - 4) +
                    1]++;  // n3AHRing to n12AHRing
      }
    }
    if (fused_ring_size > 12) {
      // greater
      if (has_hetero) {
        descriptors[83]++;  // nFHRing sum
        descriptors[93]++;  // nG12FHRing
      }
      if (is_aromatic) {
        descriptors[94]++;   // nFaRing sum
        descriptors[104]++;  // nG12FaRing
      }

      if (is_aromatic && has_hetero) {
        descriptors[105]++;  // nFaHRing sum
        descriptors[115]++;  // nG12FaHRing
      }

      if (!is_aromatic) {
        descriptors[116]++;  // nFARing sum
        descriptors[126]++;  // nG12FARing
      }

      if (!is_aromatic && has_hetero) {
        descriptors[127]++;  // nFAHRing sum
        descriptors[137]++;  // nG12FAHRing
      }
    }
  }

  return descriptors;
}

static const std::vector<std::pair<std::string, std::string>> hsPatterns = {
    {"HsOH", "[OD1H]-*"},
    {"HdNH", "[ND1H]=*"},
    {"HsSH", "[SD1H]-*"},
    {"HsNH2", "[ND1H2]-*"},
    {"HssNH", "[ND2H](-*)-*"},
    {"HaaNH", "[nD2H](:*):*"},
    {"HsNH3p", "[ND1H3]-*"},
    {"HssNH2p", "[ND2H2](-*)-*"},
    {"HsssNHp", "[ND3H](-*)(-*)-*"},
    {"HtCH", "[CD1H]#*"},
    {"HdCH2", "[CD1H2]=*"},
    {"HdsCH", "[CD2H](=*)-*"},
    {"HaaCH", "[#6D2H](:*):*"},
    {"HCHnX", "[CX4;!H0]-[F,Cl,Br,I]"},
    {"HCsats", "[CX4;!H0]~[!a]"},
    {"HCsatu", "[CX4;!H0]-[*]:,=,#[*]"},
    {"HAvin", "[CX3H](=C)-[c]"},
    {"Hall", "[*;!H0]"},
    {"Hother", "[*;!H0]:,=,#[*]"},
    {"Hmisc", "[*;!H0]~[B,Si,P,Ge,As,Se,Sn,Pb]"}};

// Define the EState atom types and their SMARTS patterns
static const std::vector<std::pair<std::string, std::string>> esPatterns = {
    {"sLi", "[LiD1]-*"},
    {"ssBe", "[BeD2](-*)-*"},
    {"ssssBe", "[BeD4](-*)(-*)(-*)-*"},
    {"ssBH", "[BD2H](-*)-*"},
    {"sssB", "[BD3](-*)(-*)-*"},
    {"ssssB", "[BD4](-*)(-*)(-*)-*"},
    {"sCH3", "[CD1H3]-*"},
    {"dCH2", "[CD1H2]=*"},
    {"ssCH2", "[CD2H2](-*)-*"},
    {"tCH", "[CD1H]#*"},
    {"dsCH", "[CD2H](=*)-*"},
    {"aaCH", "[C,c;D2H](:*):*"},
    {"sssCH", "[CD3H](-*)(-*)-*"},
    {"ddC", "[CD2H0](=*)=*"},
    {"tsC", "[CD2H0](#*)-*"},
    {"dssC", "[CD3H0](=*)(-*)-*"},
    {"aasC", "[C,c;D3H0](:*)(:*)-*"},
    {"aaaC", "[C,c;D3H0](:*)(:*):*"},
    {"ssssC", "[CD4H0](-*)(-*)(-*)-*"},
    {"sNH3", "[ND1H3]-*"},
    {"sNH2", "[ND1H2]-*"},
    {"ssNH2", "[ND2H2](-*)-*"},
    {"dNH", "[ND1H]=*"},
    {"ssNH", "[ND2H](-*)-*"},
    {"aaNH", "[N,nD2H](:*):*"},
    {"tN", "[ND1H0]#*"},
    {"sssNH", "[ND3H](-*)(-*)-*"},
    {"dsN", "[ND2H0](=*)-*"},
    {"aaN", "[N,nD2H0](:*):*"},
    {"sssN", "[ND3H0](-*)(-*)-*"},
    {"ddsN", "[ND3H0](~[OD1H0])(~[OD1H0])-,:*"},
    {"aasN", "[N,nD3H0](:*)(:*)-,:*"},
    {"ssssN", "[ND4H0](-*)(-*)(-*)-*"},
    {"sOH", "[OD1H]-*"},
    {"dO", "[OD1H0]=*"},
    {"ssO", "[OD2H0](-*)-*"},
    {"aaO", "[O,oD2H0](:*):*"},
    {"sF", "[FD1]-*"},
    {"sSiH3", "[SiD1H3]-*"},
    {"ssSiH2", "[SiD2H2](-*)-*"},
    {"sssSiH", "[SiD3H1](-*)(-*)-*"},
    {"ssssSi", "[SiD4H0](-*)(-*)(-*)-*"},
    {"sPH2", "[PD1H2]-*"},
    {"ssPH", "[PD2H1](-*)-*"},
    {"sssP", "[PD3H0](-*)(-*)-*"},
    {"dsssP", "[PD4H0](=*)(-*)(-*)-*"},
    {"sssssP", "[PD5H0](-*)(-*)(-*)(-*)-*"},
    {"sSH", "[SD1H1]-*"},
    {"dS", "[SD1H0]=*"},
    {"ssS", "[SD2H0](-*)-*"},
    {"aaS", "[S,sD2H0](:*):*"},
    {"dssS", "[SD3H0](=*)(-*)-*"},
    {"ddssS", "[SD4H0](~[OD1H0])(~[OD1H0])(-*)-*"},
    {"sCl", "[ClD1]-*"},
    {"sGeH3", "[GeD1H3](-*)"},
    {"ssGeH2", "[GeD2H2](-*)-*"},
    {"sssGeH", "[GeD3H1](-*)(-*)-*"},
    {"ssssGe", "[GeD4H0](-*)(-*)(-*)-*"},
    {"sAsH2", "[AsD1H2]-*"},
    {"ssAsH", "[AsD2H1](-*)-*"},
    {"sssAs", "[AsD3H0](-*)(-*)-*"},
    {"sssdAs", "[AsD4H0](=*)(-*)(-*)-*"},
    {"sssssAs", "[AsD5H0](-*)(-*)(-*)(-*)-*"},
    {"sSeH", "[SeD1H1]-*"},
    {"dSe", "[SeD1H0]=*"},
    {"ssSe", "[SeD2H0](-*)-*"},
    {"aaSe", "[SeD2H0](:*):*"},
    {"dssSe", "[SeD3H0](=*)(-*)-*"},
    {"ddssSe", "[SeD4H0](=*)(=*)(-*)-*"},
    {"sBr", "[BrD1]-*"},
    {"sSnH3", "[SnD1H3]-*"},
    {"ssSnH2", "[SnD2H2](-*)-*"},
    {"sssSnH", "[SnD3H1](-*)(-*)-*"},
    {"ssssSn", "[SnD4H0](-*)(-*)(-*)-*"},
    {"sI", "[ID1]-*"},
    {"sPbH3", "[PbD1H3]-*"},
    {"ssPbH2", "[PbD2H2](-*)-*"},
    {"sssPbH", "[PbD3H1](-*)(-*)-*"},
    {"ssssPb", "[PbD4H0](-*)(-*)(-*)-*"}};

// Define the EState atom types and their SMARTS patterns
static const std::vector<std::pair<std::string, std::string>>
    esPatternsFromOEState = {
        {"sLi", "[LiD1]-*"},
        {"ssBe", "[BeD2](-*)-*"},
        {"ssssBe", "[BeD4](-*)(-*)(-*)-*"},
        {"ssBH", "[BD2H](-*)-*"},
        {"sssB", "[BD3](-*)(-*)-*"},
        {"ssssB", "[BD4](-*)(-*)(-*)-*"},
        {"sCH3", "[CD1H3]-*"},
        {"dCH2", "[CD1H2]=*"},
        {"ssCH2", "[CD2H2](-*)-*"},
        {"tCH", "[CD1H]#*"},
        {"dsCH", "[CD2H](=*)-*"},
        {"aaCH", "[C,c;D2H](:*):*"},
        {"sssCH", "[CD3H](-*)(-*)-*"},
        {"ddC", "[CD2H0](=*)=*"},
        {"tsC", "[CD2H0](#*)-*"},
        {"dssC", "[CD3H0](=*)(-*)-*"},
        {"aasC", "[C,c;D3H0](:*)(:*)-*"},
        {"aaaC", "[C,c;D3H0](:*)(:*):*"},
        {"ssssC", "[CD4H0](-*)(-*)(-*)-*"},
        {"sNH3", "[ND1H3]-*"},
        {"sNH2", "[ND1H2]-*"},
        {"ssNH2", "[ND2H2](-*)-*"},
        {"dNH", "[ND1H]=*"},
        {"ssNH", "[ND2H](-*)-*"},
        {"aaNH", "[N,nD2H](:*):*"},
        {"tN", "[ND1H0]#*"},
        {"sssNH", "[ND3H](-*)(-*)-*"},
        {"dsN", "[ND2H0](=*)-*"},
        {"aaN", "[N,nD2H0](:*):*"},
        {"sssN", "[ND3H0](-*)(-*)-*"},
        {"ddsN", "[ND3H0](~[OD1H0])(~[OD1H0])-,:*"},
        {"aasN", "[N,nD3H0](:*)(:*)-,:*"},
        {"ssssN", "[ND4H0](-*)(-*)(-*)-*"},
        {"sOH", "[OD1H]-*"},
        {"dO", "[OD1H0]=*"},
        {"ssO", "[OD2H0](-*)-*"},
        {"aaO", "[O,oD2H0](:*):*"},
        {"sF", "[FD1]-*"},
        {"sSiH3", "[SiD1H3]-*"},
        {"ssSiH2", "[SiD2H2](-*)-*"},
        {"sssSiH", "[SiD3H1](-*)(-*)-*"},
        {"ssssSi", "[SiD4H0](-*)(-*)(-*)-*"},
        {"sPH2", "[PD1H2]-*"},
        {"ssPH", "[PD2H1](-*)-*"},
        {"sssP", "[PD3H0](-*)(-*)-*"},
        {"dsssP", "[PD4H0](=*)(-*)(-*)-*"},
        {"sssssP", "[PD5H0](-*)(-*)(-*)(-*)-*"},
        {"sSH", "[SD1H1]-*"},
        {"dS", "[SD1H0]=*"},
        {"ssS", "[SD2H0](-*)-*"},
        {"aaS", "[S,sD2H0](:*):*"},
        {"dssS", "[SD3H0](=*)(-*)-*"},
        {"ddssS", "[SD4H0](~[OD1H0])(~[OD1H0])(-*)-*"},
        {"sCl", "[ClD1]-*"},
        {"sGeH3", "[GeD1H3](-*)"},
        {"ssGeH2", "[GeD2H2](-*)-*"},
        {"sssGeH", "[GeD3H1](-*)(-*)-*"},
        {"ssssGe", "[GeD4H0](-*)(-*)(-*)-*"},
        {"sAsH2", "[AsD1H2]-*"},
        {"ssAsH", "[AsD2H1](-*)-*"},
        {"sssAs", "[AsD3H0](-*)(-*)-*"},
        {"sssdAs", "[AsD4H0](=*)(-*)(-*)-*"},
        {"sssssAs", "[AsD5H0](-*)(-*)(-*)(-*)-*"},
        {"sSeH", "[SeD1H1]-*"},
        {"dSe", "[SeD1H0]=*"},
        {"ssSe", "[SeD2H0](-*)-*"},
        {"aaSe", "[SeD2H0](:*):*"},
        {"dssSe", "[SeD3H0](=*)(-*)-*"},
        {"ddssSe", "[SeD4H0](=*)(=*)(-*)-*"},
        {"sBr", "[BrD1]-*"},
        {"sSnH3", "[SnD1H3]-*"},
        {"ssSnH2", "[SnD2H2](-*)-*"},
        {"sssSnH", "[SnD3H1](-*)(-*)-*"},
        {"ssssSn", "[SnD4H0](-*)(-*)(-*)-*"},
        {"sI", "[ID1]-*"},
        {"sPbH3", "[PbD1H3]-*"},
        {"ssPbH2", "[PbD2H2](-*)-*"},
        {"sssPbH", "[PbD3H1](-*)(-*)-*"},
        {"ssssPb", "[PbD4H0](-*)(-*)(-*)-*"},
        {"sNH2(A)", "[#7;D1;X3;H2][CX4;A]"},
        {"sNH2(a)", "[#7;D1;X3;H2][c]"},
        {"sNH2(oth)", "[$([#7;D1;X3;H2][CX3])]"},
        {"ssNH(A)", "[$([#7;X3;D2;H]([CX4;A][CX4;A]))]"},
        {"ssNH(a)", "[$([#7;X3;D2;H]([c;a])-*)]"},
        {"ssNH(oth)", "[$([#7;X3;D2;H]([CX3])-*)]"},
        {"sssN(A)", "[$([#7;D3;X3;H0]([CX4;A])([CX4;A])[CX4;A])]"},
        {"sssN(a)", "[$([#7;D3;X3;H0]([c;a])(-*)-*)]"},
        {"sssN(oth)", "[$([#7D3;X3;H0]([CX3])(-*)-*)]"},
        {"ddsN(nitro)", "[$([#7X3](=O)=O),$([#7X3+](=O)[O-])][!#8]"},
        {"sOH(A)", "[$([#8;X2;D1;H][CX4;A])]"},
        {"sOH(a)", "[$([#8;X2;D1;H][c;a])]"},
        {"sOH(acid)", "[$([#8;X2;D1;H][CX3](=[OX1]))]"},
        {"sOH(zwit)", "[$([#8;X2;D1;H,OX1-][CX3](=[OX1])[NX3,NX4+])]"},
        {"ssO(ester)", "[$([#8;X2;D2;H0]([CX3]=[OX1H0])[#6])]"},
        {"dOfix", "[#8;D1;X1;H0]~*"},
        {"dO(keto)", "[$([#8;X1;D1;H0]=[#6X3]([#6])[#6])]"},
        {"dO(acid)", "[$([#8;X1;D1;H0]=[#6X3]([OX2H1]))]"},
        {"dO(ester)", "[$([#8;X1;D1;H0]=[#6X3]([OX2H0])[#6])]"},
        {"dO(amid)", "[$([#8;X1;D1;H0]=[#6X3]([#7])[#6])]"},
        {"dO(nitro)", "[$([#8;X1;D1;H0]~[#7X3]~[#8;X1;D1;H0])]"},
        {"dO(sulfo)",
         "[$([#8;X1]=[#16;X4]=[#8;X1]),$([#8;X1-][#16;X4+2][#8;X1-])]"}};

// Function to add SMARTS queries safely to the vector of pairs
}
}
}
