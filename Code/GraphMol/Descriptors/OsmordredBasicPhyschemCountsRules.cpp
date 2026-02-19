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
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
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
// v2.0: Filter function to check if a molecule is too large (will cause hangs)
// Returns true if molecule has >10 rings OR >200 heavy atoms
// This prevents Osmordred.Calculate from hanging on very complex molecules
bool isMoleculeTooLarge(const ROMol &mol) {
  const RingInfo *ri = mol.getRingInfo();
  int numRings = 0;
  if (ri && ri->isInitialized()) {
    numRings = ri->numRings();
  }

  int numHeavyAtoms = mol.getNumHeavyAtoms();

  // Filter: >10 rings OR >200 heavy atoms
  return (numRings > 10 || numHeavyAtoms > 200);
}

namespace {
void solveLinearSystem(const ROMol &mol, std::vector<double> &A,
                       std::vector<double> &B, int n, int nrhs, bool &success) {
  int lda = n;  // Leading dimension of A
  int ldb = n;  // Leading dimension of B
  int info;

  success = false;  // Initialize success flag

  // CRITICAL FIX v2.0: Save original RHS before any LAPACK calls modify B
  // LAPACK routines modify B in-place, even when they fail!
  std::vector<double> B_original = B;

  // First, try dposv (Cholesky factorization for positive definite matrices)
  std::vector<double> A_copy = A;  // Copy A because LAPACK modifies it
  info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, nrhs, A_copy.data(), lda,
                       B.data(), ldb);

  if (info == 0) {
    success = true;
    return;
  } else {
    // dposv failed; fall back to dgesv (LU factorization)
    // CRITICAL FIX v2.0: Restore original RHS before calling dgesv
    // dposv modified B even though it failed!
    B = B_original;

    std::vector<int> ipiv(n);  // Pivot array for dgesv
    A_copy = A;                // Reset A because it was modified by dposv
    info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, A_copy.data(), lda,
                         ipiv.data(), B.data(), ldb);

    if (info == 0) {
      success = true;
      return;
    } else {
      // dgesv failed (singular matrix); fall back to dgelss (pseudo-inverse via
      // SVD) CRITICAL FIX v2.0: Added dgelss fallback for singular matrices
      // This provides a minimum-norm least-squares solution when exact solution
      // doesn't exist
      B = B_original;  // Restore original RHS values

      std::vector<double> A_copy2 = A;  // Fresh copy for dgelss
      std::vector<double> B_copy = B;   // Copy B because dgelss modifies it

      // Allocate workspace for dgelss
      std::vector<double> s(n);  // Singular values
      int rank;                  // Rank of matrix
      double rcond = 1e-15;      // Condition number threshold

      // dgelss computes least-squares solution: min ||Ax - b||_2
      info = LAPACKE_dgelss(LAPACK_COL_MAJOR, n, n, nrhs, A_copy2.data(), lda,
                            B_copy.data(), ldb, s.data(), rcond, &rank);

      if (info == 0) {
        // dgelss succeeded - copy solution back to B
        B = B_copy;
        success = true;
        return;
      } else {
        // All solvers failed - this is a true error
        std::string outputSmiles = MolToSmiles(mol);
        std::cerr
            << "ERROR: All LAPACK solvers failed (dposv, dgesv, dgelss): info="
            << info << ", Smiles:" << outputSmiles << "\n";
      }
    }
  }
}
}  // namespace

// Function to count the number of endocyclic single bonds
int calcEndocyclicSingleBonds(const ROMol &mol) {
  const RingInfo *ri = mol.getRingInfo();
  if (!ri || !ri->isInitialized()) {
    return 0;  // No ring information available
  }

  std::unordered_set<int> bondIndices;

  // Collect all bond indices involved in rings
  for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    if (ri->numBondRings(i) > 0) {  // Check if the bond is part of a ring
      bondIndices.insert(i);
    }
  }

  int nbonds = 0;
  for (const auto &bondIdx : bondIndices) {
    const Bond *bond = mol.getBondWithIdx(bondIdx);
    if (bond->getBondType() == Bond::SINGLE) {
      nbonds++;
    }
  }

  return nbonds;
}

static const std::vector<std::string> acidicSMARTS = {
    "[O;H1]-[C,S,P]=O", "[*;-;!$(*~[*;+])]", "[NH](S(=O)=O)C(F)(F)F",
    "n1nnnc1"};

static const std::vector<std::shared_ptr<RWMol>> &GetAcidicSmarts() {
  static const std::vector<std::shared_ptr<RWMol>> compiledAcidicSMARTS = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smarts : acidicSMARTS) {
      auto mol = SmartsToMol(smarts);
      if (mol) {
        res.emplace_back(std::shared_ptr<RWMol>(mol));
      } else {
        std::cerr << "Invalid SMARTS: " << smarts << std::endl;
      }
    }
    return res;
  }();
  return compiledAcidicSMARTS;
}

static const std::vector<std::string> basicSMARTS = {"[NH2]-[CX4]",
                                                     "[NH](-[CX4])-[CX4]",
                                                     "N(-[CX4])(-[CX4])-[CX4]",
                                                     "[*;+;!$(*~[*;-])]",
                                                     "N=C-N",
                                                     "N-C=N"};

static const std::vector<std::shared_ptr<RWMol>> &GetBasicSmarts() {
  static const std::vector<std::shared_ptr<RWMol>> compiledBasicSMARTS = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smarts : basicSMARTS) {
      auto mol = SmartsToMol(smarts);
      if (mol) {
        res.emplace_back(std::shared_ptr<RWMol>(mol));
      } else {
        std::cerr << "Invalid SMARTS: " << smarts << std::endl;
      }
    }
    return res;
  }();
  return compiledBasicSMARTS;
}

// Function to count substructure matches
int countMatches(const ROMol &mol,
                 const std::vector<std::shared_ptr<RWMol>> &patterns) {
  int count = 0;
  for (const auto &pattern : patterns) {
    if (pattern) {
      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *pattern, matches);
      count += matches.size();
    }
  }
  return count;
}

// additional features

static const std::vector<std::string> alcoholsSMARTS = {
    "[#6;H2;!$(C=O)][OX2H]",
    "[#6;H1;!$(C=O)][OX2H]",
    "[#6;H0;!$(C=O)][OX2H]",
};

// Precompile SMARTS patterns for efficiency
static const std::vector<std::shared_ptr<RWMol>> &GetAlcoholSmarts() {
  static const std::vector<std::shared_ptr<RWMol>> compiledalcoholsSMARTS = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smarts : alcoholsSMARTS) {
      auto mol = SmartsToMol(smarts);
      if (mol) {
        res.emplace_back(std::shared_ptr<RWMol>(mol));
      } else {
        std::cerr << "Invalid SMARTS: " << smarts << std::endl;
      }
    }
    return res;
  }();
  return compiledalcoholsSMARTS;
}

// Function to count primary, secondary, and tertiary hydroxyl groups
std::vector<int> countHydroxylGroups(const ROMol &mol) {
  std::vector<int> results(3, 0);

  for (size_t i = 0; i < GetAlcoholSmarts().size(); ++i) {
    std::vector<MatchVectType> matches;
    SubstructMatch(mol, *GetAlcoholSmarts()[i], matches);
    results[i] = matches.size();
  }

  return results;
}

// Define static SMARTS patterns for bridged bonds, polyacids, and polyalcohols
static const std::vector<std::string> smartsPatterns = {
    "[*x3,*x4,*x5,*x6]",        // Bridged bonds
    "[CX3](=O)[OX1H0-,OX2H1]",  // Polyacid
    "[#6;!$(C=O)][OX2H]"        // Polyalcohol
};

// Precompile SMARTS patterns for efficiency
static const std::vector<std::shared_ptr<RWMol>> &GetSmarts() {
  static const std::vector<std::shared_ptr<RWMol>> compiledSMARTS = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smarts : smartsPatterns) {
      auto mol = SmartsToMol(smarts);
      if (mol) {
        res.emplace_back(std::shared_ptr<RWMol>(mol));
      } else {
        std::cerr << "Invalid SMARTS: " << smarts << std::endl;
      }
    }
    return res;
  }();
  return compiledSMARTS;
}

// Function to count the number of bridged bonds
int countBridgedBonds(const ROMol &mol) {
  int nBridgeheads = Descriptors::calcNumBridgeheadAtoms(mol);
  if (nBridgeheads > 0) {
    std::vector<MatchVectType> matches;
    int nbonds = 0;

    if (SubstructMatch(mol, *GetSmarts()[0],
                       matches)) {  // Use precompiled pattern for bridged bonds
      for (const auto &match : matches) {
        for (const auto &atom : match) {
          nbonds += std::distance(
              mol.getAtomNeighbors(mol.getAtomWithIdx(atom.second)).first,
              mol.getAtomNeighbors(mol.getAtomWithIdx(atom.second)).second);
        }
      }
    }
    return nbonds;
  }
  return 0;
}

// Function to check if a molecule is a polyacid
bool isPolyAcid(const ROMol &mol) {
  return SubstructMatch(mol, *GetSmarts()[1]).size() >
         1;  // Use precompiled pattern for polyacid
}

// Function to check if a molecule is a polyalcohol
bool isPolyAlcohol(const ROMol &mol) {
  return SubstructMatch(mol, *GetSmarts()[2]).size() >
         1;  // Use precompiled pattern for polyalcohol
}

std::vector<double> calcAddFeatures(const ROMol &mol) {
  std::vector<double> v(7, 0.);
  auto hydroxylCounts = countHydroxylGroups(mol);
  v[0] = hydroxylCounts[0];
  v[1] = hydroxylCounts[1];
  v[2] = hydroxylCounts[2];
  v[3] = static_cast<double>(countBridgedBonds(mol));
  v[4] = static_cast<double>(isPolyAcid(mol));
  v[5] = static_cast<double>(isPolyAlcohol(mol));
  v[6] = static_cast<double>(calcEndocyclicSingleBonds(mol));

  return v;
}

// Function to calculate the number of acidic groups in a molecule
int calcAcidicGroupCount(const ROMol &mol) {
  return countMatches(mol, GetAlcoholSmarts());
}

// Function to calculate the number of basic groups in a molecule
int calcBasicGroupCount(const ROMol &mol) {
  return countMatches(mol, GetBasicSmarts());
}

// Function to calculate both acidic and basic group counts
std::vector<int> calcAcidBase(const ROMol &mol) {
  return {calcAcidicGroupCount(mol), calcBasicGroupCount(mol)};
}

std::vector<double> calcABCIndex(const ROMol &mol) {
  std::vector<double> res(2, 0.);
  double ggAbcIndex = 0.0;
  double abcIndex = 0.0;

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  unsigned int numAtoms = mol.getNumAtoms();

  for (const auto &bond : mol.bonds()) {
    int u = bond->getBeginAtomIdx();
    int v = bond->getEndAtomIdx();
    auto atom1 = bond->getBeginAtom();
    auto atom2 = bond->getEndAtom();
    int nu = 0, nv = 0;

    double du = static_cast<double>(atom1->getDegree());
    double dv = static_cast<double>(atom2->getDegree());

    abcIndex += std::sqrt((du + dv - 2.0) / (du * dv));

    for (size_t i = 0; i < numAtoms; ++i) {
      if (distances[u * numAtoms + i] < distances[v * numAtoms + i]) nu++;
      if (distances[v * numAtoms + i] < distances[u * numAtoms + i]) nv++;
    }

    ggAbcIndex += std::sqrt((nu + nv - 2.0) / (nu * nv));
  }
  res[0] = abcIndex;
  res[1] = ggAbcIndex;

  return res;
}

// Function to calculate the number of aromatic atoms in a molecule
int countAromaticAtoms(const ROMol &mol) {
  int count = 0;
  for (const auto &atom : mol.atoms()) {
    if (atom->getIsAromatic()) {
      count++;
    }
  }
  return count;
}

// Function to calculate the number of aromatic bonds in a molecule
int countAromaticBonds(const ROMol &mol) {
  int count = 0;
  for (const auto &bond : mol.bonds()) {
    if (bond->getIsAromatic()) {
      count++;
    }
  }
  return count;
}

std::vector<int> calcAromatic(const ROMol &mol) {
  return {countAromaticAtoms(mol), countAromaticBonds(mol)};
}

static const std::unordered_map<std::string, int> elementMapAtomCounts = {
    {"H", 1}, {"B", 2}, {"C", 3},  {"N", 4},   {"O", 5}, {"S", 6},
    {"P", 7}, {"F", 8}, {"Cl", 9}, {"Br", 10}, {"I", 11}};

// Function to calculate the atom count descriptor
std::vector<int> calcAtomCounts(const ROMol &mol) {
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  // Initialize the counts for each atom type

  // nAtom,nHeavyAtom,nSpiro,nBridgehead,nHetero,nH,nB,nC,nN,nO,nS,nP,nF,nCl,nBr,nI,nX
  // version higher then v1 else
  // nAtom,nHeavyAtom,nSpiro,nBridgehead,nH,nB,nC,nN,nO,nS,nP,nF,nCl,nBr,nI,nX
  // for version 1

  // std::vector<int> counts(17, 0);

  int nAtoms = hmol->getNumAtoms();
  int nHeavy = Descriptors::calcNumHeavyAtoms(mol);
  int nSpiro = Descriptors::calcNumSpiroAtoms(mol);
  int nBrigde = Descriptors::calcNumBridgeheadAtoms(mol);

  // Halogen list (F, Cl, Br, I) not use in the logic ... also restrinction
  // std::vector<int> halogens = {9, 17, 35, 53};  // Atomic numbers of halogens
  // (F, Cl, Br, I)  remark: also 85, 117 in Mordred, but honesttly we have
  // rarely the case ...
  int nH = 0, nB = 0, nC = 0, nN = 0, nO = 0, nS = 0, nP = 0, nF = 0, nCl = 0,
      nBr = 0, nI = 0;

  // Iterate over all atoms in the molecule
  for (const auto &atom : hmol->atoms()) {
    std::string symbol = atom->getSymbol();

    auto it = elementMapAtomCounts.find(symbol);
    if (it != elementMapAtomCounts.end()) {
      switch (it->second) {
        case 1:
          nH++;
          break;
        case 2:
          nB++;
          break;
        case 3:
          nC++;
          break;
        case 4:
          nN++;
          break;
        case 5:
          nO++;
          break;
        case 6:
          nS++;
          break;
        case 7:
          nP++;
          break;
        case 8:
          nF++;
          break;
        case 9:
          nCl++;
          break;
        case 10:
          nBr++;
          break;
        case 11:
          nI++;
          break;
      }
    }
  }

  int nX = nF + nCl + nBr + nI;
  int nHetero = Descriptors::calcNumHeteroatoms(mol);
  return {nAtoms, nHeavy, nSpiro, nBrigde, nHetero, nH,  nB, nC, nN,
          nO,     nS,     nP,     nF,      nCl,     nBr, nI, nX};
}

// return vector sum over rows
std::vector<double> _VertexDegrees(const double *distances,
                                   const unsigned int numatoms) {
  std::vector<double> res(numatoms, 0.0);
  double sum;
  for (unsigned int i = 0; i < numatoms; ++i) {
    sum = 0.0;
    for (unsigned int j = 0; j < numatoms; ++j) {
      sum += distances[j * numatoms + i];
    }
    res[i] = sum;
  }
  return res;
}

double BalabanJ(const ROMol &mol) {
  double q = mol.getNumBonds();
  unsigned int n = mol.getNumAtoms();

  double *Topodistances =
      MolOps::getDistanceMat(mol, false, false, false,
                             "Balaban");  // we need bond order diagonal!  based
                                          // on the code mordred all false
  double *adjMat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");

  std::vector<double> s = _VertexDegrees(Topodistances, n);
  double mu = q - n + 1;

  double sum_ = 0.0;
  for (unsigned int i = 0; i < n; ++i) {
    double si = s[i];
    for (unsigned int j = i; j < n; ++j) {
      if (adjMat[i * n + j] == 1) {
        sum_ += 1.0 / sqrt(si * s[j]);
      }
    }
  }

  double J = (mu + 1 != 0) ? (q / (mu + 1)) * sum_ : 0.0;

  return J;
}

std::vector<double> calcBalabanJ(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = BalabanJ(mol);
  return res;
}

// bertyCT related functions (InfoGain can be found in rdkit ML/InfoGainFuncs
// part, but I have to change to input for matching python code)
template <typename... Args>
std::string makeKey(Args... args) {
  std::ostringstream oss;
  ((oss << args << "_"), ...);
  std::string key = oss.str();
  key.pop_back();  // Remove the trailing underscore
  return key;
}


// Function to assign symmetry classes to each atom based on the distance matrix
std::vector<int> assignSymmetryClasses(const ROMol &mol,
                                       const std::vector<std::vector<double>> &,
                                       int numAtoms, int cutoff) {
  std::vector<int> symList(numAtoms, 0);

  double *distances = MolOps::getDistanceMat(mol, true, false, true, "Balaban");
  std::vector<std::vector<double>> distMatrix(
      numAtoms, std::vector<double>(numAtoms, 0.0));

  // Fill the distance matrix
  for (int i = 0; i < numAtoms; ++i) {
    for (int j = i; j < numAtoms; ++j) {
      distMatrix[i][j] = distances[i * numAtoms + j];
      distMatrix[j][i] = distMatrix[i][j];
    }
  }

  // To store unique symmetry classes
  std::unordered_map<std::string, int> keysSeen;
  int currentClass = 1;

  // Assign symmetry classes based on distances
  for (int i = 0; i < numAtoms; ++i) {
    std::vector<double> tmpList(distMatrix[i].begin(), distMatrix[i].end());
    std::sort(tmpList.begin(), tmpList.end());

    // Take the first 'cutoff' distances to form the unique key
    std::string key = "";
    for (int j = 0; j < std::min(cutoff, static_cast<int>(tmpList.size()));
         ++j) {
      key += std::to_string(tmpList[j]) + ",";
    }

    if (keysSeen.find(key) == keysSeen.end()) {
      keysSeen[key] = currentClass++;
    }

    symList[i] = keysSeen[key];
  }

  return symList;
}

double lookUpBondOrder(
    int atom1Id, int atom2Id,
    const std::unordered_map<std::pair<int, int>, Bond::BondType> &bondDict) {
  std::pair<int, int> theKey = std::minmax(atom1Id, atom2Id);
  auto it = bondDict.find(theKey);

  if (it == bondDict.end()) {
    return 1.0;  // Default bond order if not found
  }

  Bond::BondType bondOrder = it->second;
  if (bondOrder == Bond::AROMATIC) {
    return 1.5;  // Aromatic bond order
  }
  return static_cast<double>(bondOrder);  // Convert bond type to numeric value
}

// Function to create bondDict, neighborList, and vdList efficiently
std::tuple<std::unordered_map<std::pair<int, int>, Bond::BondType>,
           std::vector<std::vector<int>>, std::vector<int>>
CreateBondDictEtc(const ROMol &mol, int numAtoms) {
  std::unordered_map<std::pair<int, int>, Bond::BondType> bondDict;
  std::vector<std::vector<int>> nList(
      numAtoms);                         // List of neighbors for each atom
  std::vector<int> vdList(numAtoms, 0);  // Valency list, initialized to 0

  // Iterate over bonds in the molecule
  for (const auto &bond : mol.bonds()) {
    int atom1 = bond->getBeginAtomIdx();
    int atom2 = bond->getEndAtomIdx();

    // Ensure atom1 < atom2 for consistent key order in bondDict
    if (atom1 > atom2) std::swap(atom1, atom2);

    // Add bond to bondDict (aromatic bonds are marked as AROMATIC)
    bondDict[std::make_pair(atom1, atom2)] = bond->getBondType();

    // Update neighbors for atom1
    if (std::find(nList[atom1].begin(), nList[atom1].end(), atom2) ==
        nList[atom1].end()) {
      nList[atom1].push_back(atom2);
    }

    // Update neighbors for atom2
    if (std::find(nList[atom2].begin(), nList[atom2].end(), atom1) ==
        nList[atom2].end()) {
      nList[atom2].push_back(atom1);
    }
  }

  // Calculate the valency (number of neighbors) for each atom
  for (int i = 0; i < numAtoms; ++i) {
    vdList[i] = nList[i].size();
  }
  return std::make_tuple(bondDict, nList, vdList);
}

// Function to calculate the entropies of connections and atom types
double CalculateEntropies(
    const std::unordered_map<std::string, double> &connectionDict,
    const std::unordered_map<int, double> &atomTypeDict, int numAtoms) {
  // Extract connection values into a list
  std::vector<double> connectionList;
  for (const auto &[key, value] : connectionDict) {
    connectionList.push_back(value);
  }

  if (connectionList.empty() || atomTypeDict.empty()) return 0.0;

  // Total connections
  double totConnections =
      std::accumulate(connectionList.begin(), connectionList.end(), 0.0);

  // Calculate connection entropy
  double connectionIE =
      totConnections *
      (InfoEntropy(connectionList) + std::log(totConnections) / std::log(2.0));

  // Extract atom type values into a list
  std::vector<double> atomTypeList;
  for (const auto &[key, value] : atomTypeDict) {
    atomTypeList.push_back(value);
  }

  // Calculate atom type entropy
  double atomTypeIE = numAtoms * InfoEntropy(atomTypeList);

  return atomTypeIE + connectionIE;
}

// Main BertzCT function (refactored)
double BertzCT(const ROMol &mol) {
  int cutoff = 100;
  std::unordered_map<int, double> atomTypeDict;  // Maps atom type to count
  std::unordered_map<std::string, double>
      connectionDict;  // Maps bond classes to count

  int numAtoms = mol.getNumAtoms();

  double *distances = MolOps::getDistanceMat(mol, true, false, true, "Balaban");
  std::vector<std::vector<double>> dMat(numAtoms,
                                        std::vector<double>(numAtoms, 0.0));

  // Fill the distance matrix
  for (int i = 0; i < numAtoms; ++i) {
    for (int j = i; j < numAtoms; ++j) {
      dMat[i][j] = distances[i * numAtoms + j];
      dMat[j][i] = dMat[i][j];
    }
  }

  if (numAtoms < 2) return 0.0;

  // Create bondDict, neighborList, and vdList
  auto [bondDict, neighborList, vdList] = CreateBondDictEtc(mol, numAtoms);
  // Assign symmetry classes
  auto symmetryClasses = assignSymmetryClasses(mol, dMat, numAtoms, cutoff);

  // Iterate over atoms to compute atomTypeDict and connectionDict
  for (int atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
    int hingeAtomNumber = mol.getAtomWithIdx(atomIdx)->getAtomicNum();
    atomTypeDict[hingeAtomNumber]++;

    int hingeAtomClass = symmetryClasses[atomIdx];
    int numNeighbors = vdList[atomIdx];

    for (int i = 0; i < numNeighbors; ++i) {
      int neighbor_iIdx = neighborList[atomIdx][i];
      int NiClass = symmetryClasses[neighbor_iIdx];
      double bond_i_order = lookUpBondOrder(atomIdx, neighbor_iIdx, bondDict);

      if (bond_i_order > 1 && neighbor_iIdx > atomIdx) {
        double numConnections = bond_i_order * (bond_i_order - 1) / 2;
        std::string connectionKey = makeKey(std::min(hingeAtomClass, NiClass),
                                            std::max(hingeAtomClass, NiClass));
        connectionDict[connectionKey] += numConnections;
      }

      for (int j = i + 1; j < numNeighbors; ++j) {
        int neighbor_jIdx = neighborList[atomIdx][j];
        int NjClass = symmetryClasses[neighbor_jIdx];
        double bond_j_order = lookUpBondOrder(atomIdx, neighbor_jIdx, bondDict);
        double numConnections = bond_i_order * bond_j_order;
        std::string connectionKey =
            makeKey(std::min(NiClass, NjClass), hingeAtomClass,
                    std::max(NiClass, NjClass));
        connectionDict[connectionKey] += numConnections;
      }
    }
  }

  if (connectionDict.empty()) return 0.0;
  if (atomTypeDict.empty()) return 0.0;

  // Calculate and return the final entropy-based complexity value
  return CalculateEntropies(connectionDict, atomTypeDict, numAtoms);
}

std::vector<double> calcBertzCT(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = BertzCT(mol);
  return res;
}

// bondCount

std::vector<int> calcBondCounts(const ROMol &mol) {
  // Vector to hold bond counts: [Any, Single, Double, Triple, Aromatic,
  // Multiple]
  std::vector<int> bondCounts(9, 0);

  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  bondCounts[0] = hmol->getNumBonds();

  // Loop through all bonds in the molecule
  for (const auto *bond : hmol->bonds()) {
    bool isAromatic = bond->getIsAromatic();
    auto bondType = bond->getBondType();

    // Count each bond type

    if (bondType == Bond::BondType::SINGLE) {
      ++bondCounts[2];  // Single bond
    }

    if (bondType == Bond::BondType::DOUBLE) {
      ++bondCounts[3];  // Double bond
    }

    if (bondType == Bond::BondType::TRIPLE) {
      ++bondCounts[4];  // Triple bond
    }

    if (isAromatic || bondType == Bond::BondType::AROMATIC) {
      ++bondCounts[5];  // Aromatic bond
    }

    if (isAromatic || bondType != Bond::BondType::SINGLE) {
      ++bondCounts[6];  // Multiple bond
    }

    if (bond->getBeginAtom()->getAtomicNum() > 1 &&
        bond->getEndAtom()->getAtomicNum() > 1) {
      ++bondCounts[1];  // Heavy bond
    }
  }

  RWMol *kekulizedMol = new RWMol(*hmol);
  MolOps::Kekulize(*kekulizedMol, true);

  for (const auto *bond : kekulizedMol->bonds()) {
    auto bondType = bond->getBondType();

    if (bondType == Bond::BondType::SINGLE) {
      ++bondCounts[7];  // Kekulize Single bond
    }

    if (bondType == Bond::BondType::DOUBLE) {
      ++bondCounts[8];  // Kekulize Double bond
    }
  }

  delete kekulizedMol;

  return bondCounts;
}

// CarbonTypes there is an issue in the code not sure why this is not the same
// as in python code!
// TODO: need debug
std::vector<double> calcCarbonTypes(const ROMol &mol) {
  int slots = 11;

  std::vector<double> carbonTypeCounts(
      slots, 0.);         // Store counts for C1SP1, C2SP1, C1SP2, ..., C4SP3
  unsigned int nSP3 = 0;  // Count of SP3 hybridized carbons
  unsigned int nSP2 = 0;  // Count of SP2 hybridized carbons

  // Make a copy of the molecule and Kekulize it
  // TO DO everywhere : must avoid segment fault so must replace all delete by
  // this method !!!
  std::unique_ptr<RWMol> kekulizedMol(new RWMol(mol));
  MolOps::Kekulize(*kekulizedMol, true);

  for (const Atom *atom : kekulizedMol->atoms()) {
    // Check if the atom is a carbon (atomic number 6)
    if (atom->getAtomicNum() == 6) {
      // Count only neighboring carbons not the real degree (aka ignore any
      // hetero atoms)
      int carbonNeighbors = 0;
      for (const auto &neighbor : kekulizedMol->atomNeighbors(atom)) {
        if (neighbor->getAtomicNum() == 6) {
          carbonNeighbors++;
        }
      }

      Atom::HybridizationType hybridization = atom->getHybridization();

      // Check SP1 hybridization
      if (hybridization == Atom::SP) {
        if (carbonNeighbors == 1) {
          carbonTypeCounts[0]++;  // C1SP1
        } else if (carbonNeighbors == 2) {
          carbonTypeCounts[1]++;  // C2SP1
        }
      }
      // Check SP2 hybridization
      else if (hybridization == Atom::SP2) {
        nSP2++;
        if (carbonNeighbors == 1) {
          carbonTypeCounts[2]++;  // C1SP2
        } else if (carbonNeighbors == 2) {
          carbonTypeCounts[3]++;  // C2SP2
        } else if (carbonNeighbors == 3) {
          carbonTypeCounts[4]++;  // C3SP2
        }
      }
      // Check SP3 hybridization
      else if (hybridization == Atom::SP3) {
        nSP3++;
        if (carbonNeighbors == 1) {
          carbonTypeCounts[5]++;  // C1SP3
        } else if (carbonNeighbors == 2) {
          carbonTypeCounts[6]++;  // C2SP3
        } else if (carbonNeighbors == 3) {
          carbonTypeCounts[7]++;  // C3SP3
        } else if (carbonNeighbors == 4) {
          carbonTypeCounts[8]++;  // C4SP3
        }
      }
    }
  }

  // Calculate Hybridization Ratio (HybRatio)
  if (nSP2 + nSP3 > 0) {
    carbonTypeCounts[9] = static_cast<double>(nSP3) / (nSP2 + nSP3);
  }

  double fcps3 = Descriptors::calcFractionCSP3(*kekulizedMol);
  carbonTypeCounts[10] = fcps3;
  return carbonTypeCounts;
}

// VertexAdjacencyInformation
double VertexAdjacencyInformation(const ROMol &mol) {
  int m = 0;

  // Count the number of heavy-heavy bonds
  for (const auto &bond : mol.bonds()) {
    const auto *beginAtom = bond->getBeginAtom();
    const auto *endAtom = bond->getEndAtom();

    if (beginAtom->getAtomicNum() != 1 && endAtom->getAtomicNum() != 1) {
      ++m;
    }
  }

  // Calculate the descriptor value
  return 1.0 + std::log2(static_cast<double>(m));
}

std::vector<double> calcVertexAdjacencyInformation(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = VertexAdjacencyInformation(mol);
  return res;
}

// WalkCount
// TODO: not correct output not sure why to be investigation equation ...
// Function to compute the adjacency matrix
Eigen::MatrixXd computeAdjacencyMatrix(const ROMol &mol,
                                       bool useBondOrder = false) {
  unsigned int nAtoms = mol.getNumAtoms();
  double *adjFlat =
      MolOps::getAdjacencyMatrix(mol, useBondOrder, false, false, "noBO");

  Eigen::MatrixXd adjMat(nAtoms, nAtoms);

  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = 0; j < nAtoms; ++j) {
      adjMat(i, j) = adjFlat[i * nAtoms + j];
    }
  }
  return adjMat;
}

void upperTriangularSelfProduct(const Eigen::MatrixXd &A,
                                Eigen::MatrixXd &result) {
  int n = A.rows();
  result = Eigen::MatrixXd::Zero(n, n);

  // Compute only upper triangle
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      for (int k = i; k <= j; ++k) {
        result(i, j) += A(i, k) * A(j, k);
      }
      if (i != j) {
        result(j, i) = result(i, j);  // Fill the symmetric part
      }
    }
  }
}

// Function to compute walk counts Need a Product Matrix order
std::vector<double> calcWalkCounts(const ROMol &mol) {
  const int maxOrder = 10;  // Maximum order of walks
  unsigned int nAtoms = mol.getNumAtoms();
  Eigen::MatrixXd adjMat = computeAdjacencyMatrix(mol);

  std::vector<double> results(21, 0.0);

  Eigen::MatrixXd powerMatrix = Eigen::MatrixXd::Identity(nAtoms, nAtoms);

  // we need to initialize at nAtoms both totals
  double totalMWC10 = nAtoms, totalSRW10 = nAtoms;

  // #pragma omp parallel for reduction(+:totalMWC10, totalSRW10)
  for (int order = 1; order <= maxOrder; ++order) {
    if (order == 1) {
      powerMatrix = adjMat;  // A^1
    } else {
      powerMatrix = powerMatrix * adjMat;  // A^order
    }

    // Compute MWC ie full
    double mwc =
        (order == 1) ? 0.5 * adjMat.sum() : std::log(powerMatrix.sum() + 1.0);
    results[order - 1] = mwc;  // Store MWC0n

    // Compute SRW ie trace
    double srw = std::log(powerMatrix.diagonal().sum() + 1.0);
    if (order > 1) {
      results[maxOrder + order - 1] = srw;  // Store SRW0n
    }

    // Accumulate totals for MWC10 and SRW10
    totalMWC10 += mwc;
    totalSRW10 += srw;
  }

  results[maxOrder] = totalMWC10;  // TMWC10
  results[20] = totalSRW10;        // TSRW10

  return results;
}

/*
    std::vector<double> calcWalkCountsBlas(const ROMol &mol) {
        const int maxOrder = 10;  // Maximum order of walks
        unsigned int nAtoms = mol.getNumAtoms();

        // Get adjacency matrix from RDKit
        std::vector<double> adjMat(nAtoms * nAtoms, 0.0);
        double* adjFlat = MolOps::getAdjacencyMatrix(mol, false, false, false,
   "noBO"); std::copy(adjFlat, adjFlat + nAtoms * nAtoms, adjMat.begin());

        std::vector<double> powerMatrix(adjMat); // Start with A^1
        std::vector<double> results(21, 0.0);

        // Initialize totals with the number of atoms
        double totalMWC10 = nAtoms, totalSRW10 = nAtoms;

        for (int order = 1; order <= maxOrder; ++order) {
            if (order > 1) {
                // Compute powerMatrix = powerMatrix * adjMat using BLAS (dgemm:
   C = alpha*A*B + beta*C) std::vector<double> tempMatrix(nAtoms * nAtoms, 0.0);
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            nAtoms, nAtoms, nAtoms,
                            1.0, powerMatrix.data(), nAtoms,
                                adjMat.data(), nAtoms,
                            0.0, tempMatrix.data(), nAtoms);
                powerMatrix = tempMatrix;
            }

            // Compute MWC (full matrix sum)
            double mwc = (order == 1) ? 0.5 * std::accumulate(adjMat.begin(),
   adjMat.end(), 0.0) : std::log(std::accumulate(powerMatrix.begin(),
   powerMatrix.end(), 0.0) + 1.0); results[order - 1] = mwc;

            // Compute SRW (sum of diagonal elements)
            double srw = 0.0;
            for (unsigned int i = 0; i < nAtoms; ++i) {
                srw += powerMatrix[i * nAtoms + i];
            }
            srw = std::log(srw + 1.0);

            if (order > 1) {
                results[maxOrder + order - 1] = srw;
            }

            // Accumulate totals
            totalMWC10 += mwc;
            totalSRW10 += srw;
        }

        results[maxOrder] = totalMWC10;  // Store TMWC10
        results[20] = totalSRW10;  // Store TSRW10

        return results;
    }
*/
// Weight
// trick is to add the Hs for the average not only heavy atoms!
// we need a function that can do this trick!!!
std::vector<double> calcWeight(const ROMol &mol) {
  std::vector<double> W(2, 0.);
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));
  W[0] = Descriptors::calcExactMW(mol);
  int fullatomsnumber = hmol->getNumAtoms();
  W[1] = W[0] / fullatomsnumber;
  return W;
}

// Wiener Index
// working!
std::vector<int> calcWienerIndex(const ROMol &mol) {
  std::vector<int> WI(2, 0);

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  unsigned int numAtoms = mol.getNumAtoms();

  for (size_t i = 0; i < numAtoms; ++i) {
    for (size_t j = 0; j < numAtoms; ++j) {
      if (distances[i * numAtoms + j] == 3) {
        WI[1] += 1.;
      }
      WI[0] += distances[i * numAtoms + j];
    }
  }
  // must return the half part (of the Distance matrix!) or only compute half
  // matrix... to speed up!!!!
  WI[0] = static_cast<int>(WI[0] * 0.5);
  WI[1] = static_cast<int>(WI[1] * 0.5);

  return WI;
}

// Perform eigen decomposition on a symmetric matrix
std::pair<std::vector<double>, std::vector<std::vector<double>>>
eigenDecompositionSymmetric(const std::vector<std::vector<double>> &matrix) {
  int n = matrix.size();
  assert(matrix.size() == matrix[0].size() && "Matrix must be square");

  // Convert std::vector<std::vector<double>> to a 1D array in column-major
  // order for LAPACK
  std::vector<double> matrixData(n * n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      matrixData[j * n + i] = matrix[i][j];  // Column-major order
    }
  }

  // Storage for eigenvalues
  std::vector<double> eigenValues(n);

  // Call LAPACK's dsyev to compute eigenvalues and eigenvectors
  int info;
  std::vector<double> work(1);
  int lwork = -1;  // Request optimal workspace size
  info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'U', n, matrixData.data(), n,
                            eigenValues.data(), work.data(), lwork);

  if (info != 0) {
    throw std::runtime_error(
        "Error querying optimal workspace size for LAPACK dsyev");
  }

  lwork = static_cast<int>(work[0]);
  work.resize(lwork);

  // Perform eigen decomposition
  info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'U', n, matrixData.data(), n,
                            eigenValues.data(), work.data(), lwork);
  if (info != 0) {
    throw std::runtime_error(
        "LAPACK dsyev failed to compute eigen decomposition");
  }

  // Convert the eigenvectors back to std::vector<std::vector<double>>
  std::vector<std::vector<double>> eigenVectors(n, std::vector<double>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      eigenVectors[i][j] = matrixData[j * n + i];  // Column-major to row-major
    }
  }

  return {eigenValues, eigenVectors};
}

// Zagreb
// TODO: debug YES
// this is correct but not matching the results myabe H's are added ???
std::vector<double> calcZagrebIndex(const ROMol &mol) {
  std::vector<double> zagrebIndex(4, 0.0);

  // Create a vector to store the degree of each atom
  std::vector<double> degrees = calcValence(mol);

  // Zagreb Index M1: sum of degrees raised to the power of 2 * lambda
  for (int degree : degrees) {
    zagrebIndex[0] += std::pow(degree, 2.);
    zagrebIndex[2] += std::pow(degree, -2.);
  }

  // Zagreb Index M2: sum of products of degrees of connected atoms raised to
  // the power of lambda of (-1,1)

  for (const Bond *bond : mol.bonds()) {
    zagrebIndex[1] += std::pow(
        degrees[bond->getBeginAtomIdx()] * degrees[bond->getEndAtomIdx()], 1.);
    zagrebIndex[3] += std::pow(
        degrees[bond->getBeginAtomIdx()] * degrees[bond->getEndAtomIdx()], -1.);
  }

  return zagrebIndex;
}

// Define Bondi radii for atomic contributions
const std::unordered_map<int, double> bondiRadii = {
    {1, 1.20},   // H
    {5, 2.13},   // B
    {6, 1.70},   // C
    {7, 1.55},   // N
    {8, 1.52},   // O
    {9, 1.47},   // F
    {14, 2.10},  // Si
    {15, 1.80},  // P
    {16, 1.80},  // S
    {17, 1.75},  // Cl
    {35, 1.85},  // Br
    {33, 1.85},  // As
    {34, 1.90},  // Se
    {53, 1.98}   // I missing in Mordred from paper/CDK implementation
};

// Precompute atomic contributions using 4/3 * pi * r^3
const std::unordered_map<int, double> atomContributions = []() {
  std::unordered_map<int, double> contribs;
  for (const auto &[atomicNum, radius] : bondiRadii) {
    contribs[atomicNum] = (4.0 / 3.0) * M_PI * std::pow(radius, 3);
  }
  return contribs;
}();

// VdwVolumeABC
// working "Need Hs explicit!"
double VdwVolumeABC(const ROMol &mol) {
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  // Nb is the number of bonds
  // NRa is the number of aromatic rings
  // NRA is the number of aliphatic rings
  int NRA = Descriptors::calcNumAliphaticRings(*hmol);
  int NRa = Descriptors::calcNumAromaticRings(*hmol);
  int Nb = hmol->getNumBonds();

  double ac = 0.0;

  // Calculate the sum of atomic contributions
  for (const auto &atom : hmol->atoms()) {
    int atomicNum = atom->getAtomicNum();
    auto it = atomContributions.find(atomicNum);
    if (it == atomContributions.end()) {
      throw std::runtime_error("Unknown atom type encountered in molecule.");
    }
    ac += it->second;
  }

  // Compute van der Waals volume
  return ac - 5.92 * Nb - 14.7 * NRa - 3.8 * NRA;
}

std::vector<double> calcVdwVolumeABC(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = VdwVolumeABC(mol);
  return res;
}

// TopoPSA (adding S & P atoms effect to rdkit version!)
// working for my molecules but need a S, P molecules to check ...
// Helper function to calculate hydrogen count
int hydrogenCount(const Atom *atom) {
  int nH = atom->getTotalNumHs();
  for (const auto &neighbor : atom->getOwningMol().atomNeighbors(atom)) {
    if (neighbor->getAtomicNum() == 1) {
      ++nH;  // Count the hydrogens attached to the atom
    }
  }
  return nH;
}

// Helper function to count bond types
std::unordered_map<Bond::BondType, int> bondTypeCount(const Atom *atom) {
  std::unordered_map<Bond::BondType, int> bondCounts;
  for (const auto &bond : atom->getOwningMol().atomBonds(atom)) {
    if (bond->getIsAromatic()) {
      bondCounts[Bond::BondType::AROMATIC] += 1;
    } else {
      bondCounts[bond->getBondType()] += 1;
    }
  }
  return bondCounts;
}

// Function to calculate the phosphorus contribution to TPSA
double getPhosphorusContribution(const Atom *atom) {
  int nH = hydrogenCount(atom);
  auto cnt = bondTypeCount(atom);

  if (atom->getFormalCharge() != 0 || atom->getIsAromatic()) {
    return 0.0;
  }

  if (nH == 1 &&
      cnt == std::unordered_map<Bond::BondType, int>{
                 {Bond::BondType::SINGLE, 3}, {Bond::BondType::DOUBLE, 1}}) {
    return 23.47;
  } else if (nH == 0) {
    if (cnt ==
        std::unordered_map<Bond::BondType, int>{{Bond::BondType::SINGLE, 3}}) {
      return 13.59;
    } else if (cnt ==
               std::unordered_map<Bond::BondType, int>{
                   {Bond::BondType::SINGLE, 1}, {Bond::BondType::DOUBLE, 1}}) {
      return 34.14;
    } else if (cnt ==
               std::unordered_map<Bond::BondType, int>{
                   {Bond::BondType::SINGLE, 3}, {Bond::BondType::DOUBLE, 1}}) {
      return 9.81;
    }
  }
  return 0.0;
}

// Function to calculate the sulfur contribution to TPSA
double getSulfurContribution(const Atom *atom) {
  int nH = hydrogenCount(atom);
  auto cnt = bondTypeCount(atom);

  if (atom->getFormalCharge() != 0) {
    return 0.0;
  }

  if (atom->getIsAromatic()) {
    if (nH == 0) {
      if (cnt == std::unordered_map<Bond::BondType, int>{
                     {Bond::BondType::AROMATIC, 2}}) {
        return 28.24;
      } else if (cnt == std::unordered_map<Bond::BondType, int>{
                            {Bond::BondType::AROMATIC, 2},
                            {Bond::BondType::DOUBLE, 1}}) {
        return 21.70;
      }
    }
  } else {
    if (nH == 1 && cnt == std::unordered_map<Bond::BondType, int>{
                              {Bond::BondType::SINGLE, 2}}) {
      return 38.80;
    } else if (nH == 0) {
      if (cnt == std::unordered_map<Bond::BondType, int>{
                     {Bond::BondType::SINGLE, 2}}) {
        return 25.30;
      } else if (cnt == std::unordered_map<Bond::BondType, int>{
                            {Bond::BondType::DOUBLE, 1}}) {
        return 32.09;
      } else if (cnt == std::unordered_map<Bond::BondType, int>{
                            {Bond::BondType::SINGLE, 2},
                            {Bond::BondType::DOUBLE, 1}}) {
        return 19.21;
      } else if (cnt == std::unordered_map<Bond::BondType, int>{
                            {Bond::BondType::SINGLE, 2},
                            {Bond::BondType::DOUBLE, 2}}) {
        return 8.38;
      }
    }
  }
  return 0.0;
}

// Main function to calculate the Topological Polar Surface Area (TPSA)
std::vector<double> calcTopoPSA(const ROMol &mol) {
  double tpsa = Descriptors::calcTPSA(mol);

  std::vector<double> res(2, 0.0);
  res[0] = tpsa;

  for (const auto &atom : mol.atoms()) {
    int atomicNum = atom->getAtomicNum();
    if (atomicNum == 15) {  // Phosphorus
      tpsa += getPhosphorusContribution(atom);
    } else if (atomicNum == 16) {  // Sulfur
      tpsa += getSulfurContribution(atom);
    }
  }
  res[1] = tpsa;

  return res;
}

// Main function to calculate the Topological Polar Surface Area (TPSA)
std::vector<double> calcSLogP(const ROMol &mol) {
  std::vector<double> res(2, 0.0);
  double logP;
  double MR;
  Descriptors::calcCrippenDescriptors(mol, logP, MR);

  res[0] = logP;
  res[1] = MR;

  return res;
}

// Hydrogen from Rdkit code
std::vector<double> calcHydrogenBond(const ROMol &mol) {
  std::vector<double> res(2, 0.);

  int nHBAcc = Descriptors::calcNumHBA(mol);

  int nHBDon = Descriptors::calcNumHBD(mol);
  res[0] = static_cast<double>(nHBAcc);
  res[1] = static_cast<double>(nHBDon);

  return res;
}

// MOE need to implement EState here ;-)

// Helper function to calculate VSA_EState
// Helper function for LabuteASA calculation
std::vector<double> VSAcontrib(const ROMol &mol, bool includeHs = true) {
  // Compute the Labute ASA contributions
  std::vector<double> atomContribs(mol.getNumAtoms(), 0.0);
  double totalASA = 0.0;
  Descriptors::getLabuteAtomContribs(mol, atomContribs, totalASA, includeHs);

  // Include the total hydrogen contribution as the first element
  std::vector<double> Vi(atomContribs.size() + 1, 0.0);

  Vi[0] = totalASA;

  std::copy(atomContribs.begin(), atomContribs.end(), Vi.begin() + 1);

  // Cache the result as a property in the molecule
  // mol.setProp("_LabuteContribs", Vi);

  return Vi;
}

// Function to calculate EState-based VSA contributions
std::vector<double> VSA_EState(const ROMol &mol,
                               const std::vector<double> &bins) {
  // Calculate EState indices
  std::vector<double> propContribs = calcEStateIndices(mol);

  // Calculate VSA contributions
  std::vector<double> volContribs = VSAcontrib(mol, true);

  // Initialize result vector
  std::vector<double> ans(bins.size() + 1, 0.0);

  // Assign contributions to bins
  for (size_t i = 0; i < propContribs.size(); ++i) {
    if (!std::isnan(propContribs[i])) {
      double volume = volContribs[i + 1];
      auto nbin =
          std::upper_bound(bins.begin(), bins.end(), volume) - bins.begin();
      ans[nbin] += propContribs[i];
    }
  }
  // Cache the result in the molecule
  // mol.setProp("_VSA_EState", ans);

  return ans;
}

// Function to calculate EState-based VSA contributions
std::vector<double> EState_VSA(const ROMol &mol,
                               const std::vector<double> &bins) {
  // Calculate EState indices
  std::vector<double> propContribs = calcEStateIndices(mol);

  // Calculate VSA contributions
  std::vector<double> volContribs = VSAcontrib(mol, true);

  // Initialize result vector
  std::vector<double> ans(bins.size() + 1, 0.0);

  // Assign contributions to bins
  for (size_t i = 0; i < propContribs.size(); ++i) {
    if (!std::isnan(propContribs[i])) {
      auto nbin = std::upper_bound(bins.begin(), bins.end(), propContribs[i]) -
                  bins.begin();
      ans[nbin] += volContribs[i + 1];
    }
  }
  // Cache the result in the molecule
  // mol.setProp("_EState_VSA", ans);

  return ans;
}

// this is fixed by p = distance + 1
std::vector<double> calcVSA_EState(const ROMol &mol) {
  const std::string customPropName =
      "EState";  // Assume EState values are stored in this property
  std::vector<double> vsaBins = {4.78, 5.00, 5.410, 5.740, 6.00,
                                 6.07, 6.45, 7.00,  11.0};
  std::vector<double> vsaEstate = VSA_EState(mol, vsaBins);
  return vsaEstate;
}

//
std::vector<double> calcEState_VSA(const ROMol &mol) {
  const std::string customPropName =
      "EState";  // Assume EState values are stored in this property
  std::vector<double> estateBins = {-0.390, 0.290, 0.717, 1.165, 1.540,
                                    1.807,  2.05,  4.69,  9.17,  15.0};

  // Calculate the EState contributions for each atom
  std::vector<double> estateIndices = calcEStateIndices(mol);
  std::vector<double> eStateVSA = EState_VSA(mol, estateBins);

  // Use calcCustomProp_VSA with EState bins
  // std::vector<double> eStateVSA = Descriptors::calcCustomProp_VSA(
  //    mol, customPropName, estateBins);

  return eStateVSA;
}

std::vector<double> calcMoeType(const ROMol &mol) {
  std::vector<double> res(54, 0.0);
  double LabuteASA = Descriptors::calcLabuteASA(mol);
  res[0] = LabuteASA;
  int p = 1;

  std::vector<double> PEOE_VSA = Descriptors::calcPEOE_VSA(mol);
  // std::cout << PEOE_VSA.size() << " ";
  for (size_t i = 0; i < PEOE_VSA.size(); ++i) {
    res[p + i] = PEOE_VSA[i];  // Copy PEOE_VSA to res[1:13]
  }
  p += PEOE_VSA.size() - 1;
  // std::cout << "- new p: " << p << " | ";

  std::vector<double> SMR_VSA = Descriptors::calcSMR_VSA(mol);
  // std::cout << SMR_VSA.size() << " ";

  for (size_t i = 0; i < SMR_VSA.size(); ++i) {
    res[p + i] = SMR_VSA[i];  // Copy SMR_VSA to res[14:23]
  }
  p += SMR_VSA.size() - 1;
  // std::cout << "- new p: " << p << " | ";

  std::vector<double> SlogP_VSA = Descriptors::calcSlogP_VSA(mol);
  // std::cout << SlogP_VSA.size() << " ";
  for (size_t i = 0; i < SlogP_VSA.size(); ++i) {
    res[p + i] = SlogP_VSA[i];  // Copy SlogP_VSA to res[23:34]
  }
  p += SlogP_VSA.size() - 1;
  // std::cout << "- new p: " <<  p << " | ";

  std::vector<double> EState_VSA = calcEState_VSA(mol);
  // std::cout << EState_VSA.size() << " ";
  for (size_t i = 0; i < EState_VSA.size(); ++i) {
    res[p + i] = EState_VSA[i];  // Copy EState_VSA
  }
  p += EState_VSA.size() - 1;
  // std::cout << "- new p: " << p << " | ";

  // EState (mordred) VSA_EState 1-9 & EState_VSA 1-10
  std::vector<double> VSA_EState = calcVSA_EState(mol);
  // std::cout << VSA_EState.size() << "\n";
  for (size_t i = 0; i < VSA_EState.size(); ++i) {
    res[p + i] = VSA_EState[i];  // Copy SlogP_VSA to res[24:34]
                                 // std::cout << " ( " << p+i;
  }
  // std::cout << ")\n";
  p += VSA_EState.size() - 1;
  // std::cout << "- new p: " << p << " end ";
  return res;
}
/// logS

// SMARTS patterns with their corresponding log contributions
static const std::vector<std::pair<std::string, double>> smartsLogs = {
    {"[NH0;X3;v3]", 0.71535}, {"[NH2;X3;v3]", 0.41056},
    {"[nH0;X3]", 0.82535},    {"[OH0;X2;v2]", 0.31464},
    {"[OH0;X1;v2]", 0.14787}, {"[OH1;X2;v2]", 0.62998},
    {"[CH2;!R]", -0.35634},   {"[CH3;!R]", -0.33888},
    {"[CH0;R]", -0.21912},    {"[CH2;R]", -0.23057},
    {"[ch0]", -0.37570},      {"[ch1]", -0.22435},
    {"F", -0.21728},          {"Cl", -0.49721},
    {"Br", -0.57982},         {"I", -0.51547},
};

// Precompile SMARTS patterns to avoid repeated parsing
static const std::vector<std::pair<std::shared_ptr<RWMol>, double>> &
GetSmartsLogs() {
  static const std::vector<std::pair<std::shared_ptr<RWMol>, double>>
      compiledSmartsLogs = [] {
        std::vector<std::pair<std::shared_ptr<RWMol>, double>> res;
        for (const auto &pair : smartsLogs) {
          auto mol = SmartsToMol(pair.first);
          if (mol) {
            res.emplace_back(std::shared_ptr<RWMol>(mol), pair.second);
          } else {
            std::cerr << "Invalid SMARTS: " << pair.first << std::endl;
          }
        }
        return res;
      }();
  return compiledSmartsLogs;
}
// Function to calculate LogS descriptor
double LogS(const ROMol &mol) {
  // Base formula contribution
  double molWeight = Descriptors::calcAMW(mol, false);  // Get molecular weight
  double logS = 0.89823 - 0.10369 * std::sqrt(molWeight);

  // Add contributions from precompiled SMARTS patterns
  for (const auto &pair : GetSmartsLogs()) {
    auto &smartsMol = pair.first;
    double logContribution = pair.second;

    if (!smartsMol) {
      continue;  // Skip invalid SMARTS
    }

    // Match SMARTS pattern
    std::vector<MatchVectType> matches;
    SubstructMatch(mol, *smartsMol, matches);

    // Add contributions for each match
    logS += matches.size() * logContribution;
  }

  return logS;
}

std::vector<double> calcLogS(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = LogS(mol);
  return res;
}

// Function to calculate Lipinski rule of 5
int calculateLipinski(double LogP, double MW, double HBDon, double HBAcc) {
  bool L = (HBDon <= 5 && HBAcc <= 10 && MW <= 500 && LogP <= 5);

  if (L) {
    return 1;
  } else {
    return 0;
  }
}

// Function to calculate Ghose filter
int calculateGhoseFilter(double MW, double LogP, double MR, int numAtoms) {
  bool G = (MW >= 160 && MW <= 480) && (numAtoms >= 20 && numAtoms <= 70) &&
           (LogP >= -0.4 && LogP <= 5.6) && (MR >= 40 && MR <= 130);

  if (G) {
    return 1;
  } else {
    return 0;
  }
}

// Main function to calculate Lipinski and Ghose filter
std::vector<int> calcLipinskiGhose(const ROMol &mol) {
  double MW = Descriptors::calcExactMW(mol);
  double HBDon = static_cast<double>(Descriptors::calcNumHBD(mol));
  double HBAcc = static_cast<double>(Descriptors::calcNumHBA(mol));
  double LogP;
  double MR;
  Descriptors::calcCrippenDescriptors(mol, LogP, MR);

  int lipinski = calculateLipinski(LogP, MW, HBDon, HBAcc);
  // must add Hs for Ghose

  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  int numAtoms = hmol->getNumAtoms();

  int ghoseFilter = calculateGhoseFilter(MW, LogP, MR, numAtoms);

  return {lipinski, ghoseFilter};
}


double McGowanVolume(const ROMol &mol) {
  // In Padel code this is /100 in order to match the Polarisability equation
  double res = 0.;
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  std::map<int, double> mgvmap = McGowanVolumAtomicMap();

  for (const auto &atom : hmol->atoms()) {
    int atomicNum = atom->getAtomicNum();
    if (mgvmap.find(atomicNum) != mgvmap.end()) {
      res += mgvmap.at(atomicNum);
    }
  }
  double finalres = res - hmol->getNumBonds() * 6.56;

  return finalres;
}

std::vector<double> calcMcGowanVolume(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = McGowanVolume(mol);
  return res;
}

// SMARTS patterns for fragments
static const std::vector<std::string> PolFrags = {
    "[CX4H3]",
    "[CX4H2]",
    "[CX4H1]",
    "F",
    "Cl",
    "Br",
    "I",
    "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "[N,n]",
    "[O,o]",
    "[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]",
    "[S,s]",
    "[P,p]"};

// Corresponding coefficients for the fragments
static const std::vector<double> coefPol = {
    10.152, 8.765, 5.702, 3.833,  16.557, 24.123, 38.506,
    10.488, 6.335, 4.307, 15.726, 22.366, 11.173};

// Precompile SMARTS patterns for efficiency

static const std::vector<std::shared_ptr<RWMol>> &GetCompiledPolFrags() {
  static const std::vector<std::shared_ptr<RWMol>> compiledPolFrags = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smarts : PolFrags) {
      auto mol = SmartsToMol(smarts);
      if (mol) {
        res.emplace_back(std::shared_ptr<RWMol>(mol));
      } else {
        std::cerr << "Invalid SMARTS: " << smarts << std::endl;
        res.emplace_back(nullptr);  // Placeholder for invalid SMARTS
      }
    }
    return res;
  }();

  return compiledPolFrags;
}

// Function to count hydrogen atoms in the molecule
int getNumHs(const ROMol &mol) {
  int nHs = 0;
  for (const auto atom : mol.atoms()) {
    nHs += atom->getTotalNumHs();
  }
  return nHs;
}

// Function to calculate polarity descriptor
double Polarity(const ROMol &mol) {
  double res = -1.529;  // Intercept value

  // Add contributions from precompiled SMARTS patterns
  for (size_t i = 0; i < GetCompiledPolFrags().size(); ++i) {
    auto &pattern = GetCompiledPolFrags()[i];
    if (!pattern) continue;  // Skip invalid patterns

    std::vector<MatchVectType> matches;
    SubstructMatch(mol, *pattern, matches, true);  // uniquify = true
    res += matches.size() * coefPol[i];
  }

  // Add hydrogen contribution
  res += 3.391 * static_cast<double>(getNumHs(mol));

  return res;
}

std::vector<double> calcPol(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = Polarity(mol);
  return res;
}

double MRvalue(const ROMol &mol) { return 4. / 3. * M_PI * Polarity(mol); }

std::vector<double> calcMR(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = MRvalue(mol);
  return res;
}

double ODT(const ROMol &) { return 1; }

std::vector<double> calcODT(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = ODT(mol);
  return res;
}

// Function to calculate the Schultz descriptor
double Schultz(const ROMol &mol) {
  // Get the number of atoms in the molecule
  int nAtoms = mol.getNumAtoms();
  if (nAtoms == 0) return 0.0;

  // Get the distance matrix
  const double *distMat = MolOps::getDistanceMat(mol, false, false, false);

  // Get the adjacency matrix
  const double *adjMat =
      MolOps::getAdjacencyMatrix(mol, false, false, false, "noBO");
  // Calculate the vertex degree for each atom
  std::vector<double> vertexDegree(nAtoms, 0.0);
  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      vertexDegree[i] += adjMat[i * nAtoms + j];
    }
  }

  // Compute the Schultz descriptor
  double schultz = 0.0;
  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      double distance = distMat[i * nAtoms + j];
      double adjacency = adjMat[i * nAtoms + j];
      schultz += (distance + adjacency) * vertexDegree[j];
    }
  }

  return schultz;
}
std::vector<double> calcSchultz(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = Schultz(mol);
  return res;
}


// Combined function for calculating both atomic and bond polarizability
std::vector<double> calcPolarizability(const ROMol &mol) {
  double atomicPol = 0.0;
  double bondPol = 0.0;
  const auto &polmap = Polarizability94AtomicMap();
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  for (const auto &atom : hmol->atoms()) {
    int atomicNum = atom->getAtomicNum();
    auto it = polmap.find(atomicNum);
    if (it != polmap.end()) {
      atomicPol += it->second;
    }
  }
  // Calculate bond polarizability
  for (const auto &bond : hmol->bonds()) {
    int atomicNum1 = bond->getBeginAtom()->getAtomicNum();
    int atomicNum2 = bond->getEndAtom()->getAtomicNum();

    auto it1 = polmap.find(atomicNum1);
    auto it2 = polmap.find(atomicNum2);

    if (it1 != polmap.end() && it2 != polmap.end()) {
      bondPol += std::abs(it1->second - it2->second);
    }
  }

  return {atomicPol, bondPol};
}

// Main Rotatabond
std::vector<double> calcRotatableBond(const ROMol &mol) {
  double rot = Descriptors::calcNumRotatableBonds(mol);

  std::vector<double> res(2, 0.0);
  res[0] = rot;

  int bondcountheavy = 0;
  for (const auto &bond : mol.bonds()) {
    if (bond->getBeginAtom()->getAtomicNum() != 1 &&
        bond->getEndAtom()->getAtomicNum() != 1) {
      bondcountheavy += 1;  // Heavy bond
    }
  }
  if (bondcountheavy > 0) {
    res[1] = static_cast<double>(static_cast<double>(rot) /
                                 static_cast<double>(bondcountheavy));
  }
  return res;
}

double FragmentComplexity(const ROMol &mol) {
  // Number of atoms (A)
  int A = mol.getNumAtoms();

  // Number of bonds (B)
  int B = mol.getNumBonds();

  // Number of heteroatoms (H) - atoms that are not Carbon
  int H = 0;
  for (const auto &atom : mol.atoms()) {
    if (atom->getAtomicNum() != 6) {  // Exclude Carbon atoms
      H++;
    }
  }

  // Calculate fragment complexity: |B^2 - A^2 + A| + H / 100
  double fragCpx = std::abs(std::pow(B, 2) - std::pow(A, 2) + A) + H / 100.0;

  return fragCpx;
}

std::vector<double> calcFragmentComplexity(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = FragmentComplexity(mol);
  return res;
}

// it is not correct!!!!
std::vector<double> calcEccentricConnectivityIndex(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  Eigen::VectorXd E = calculateEccentricity(mol);
  std::vector<double> V = calcValence(mol);

  if ((size_t)E.size() != V.size()) {
    throw std::invalid_argument(
        "Eccentricity and Valence vectors must have the same length.");
  }

  // Calculate element-wise product of E and V
  double productSum = 0.0;
  for (size_t i = 0; i < (size_t)E.size(); ++i) {
    productSum += static_cast<int>(E[i]) *
                  V[i];  // Cast E[i] to int and multiply with V[i]
  }
  res[0] = productSum;
  // Return the result as an integer
  return res;  // static_cast<int>(productSum);
}

// Function to find rings in a molecule using RDKit's ring detection
std::vector<std::vector<int>> findRings(const ROMol &mol) {
  std::vector<std::vector<int>> rings;
  RingInfo *ri = mol.getRingInfo();
  if (ri) {
    for (size_t i = 0; i < ri->numRings(); ++i) {
      std::vector<int> ring = ri->bondRings()[i];
      rings.push_back(ring);
    }
  }
  return rings;
}

// Constitutional
// code vs name
// Z    a.GetAtomicNum()
// m    mass[a.GetAtomicNum()]
// v    vdw_volume[a.GetAtomicNum()]
// se   sanderson[a.GetAtomicNum()]
// pe   pauling[a.GetAtomicNum()]
// are  allred_rocow[a.GetAtomicNum()]
// p    polarizability94[a.GetAtomicNum()] (last as default!!!)
// i    ionization_potentials[a.GetAtomicNum()]

std::vector<double> calcConstitutional(const ROMol &mol) {
  double SZ = 0.;
  double Sm = 0.;
  double Sv = 0.;
  double Sse = 0.;
  double Spe = 0.;
  double Sare = 0.;
  double Sp = 0.;
  double Si = 0.;
  double MZ, Mm, Mv, Mse, Mpe, Mare, Mp, Mi;
  const PeriodicTable *tbl = PeriodicTable::getTable();
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  double zcc = static_cast<double>(tbl->getAtomicNumber("C"));

  double mcc = static_cast<double>(tbl->getAtomicWeight("C"));

  const std::map<int, double> &vdwmap = VdWAtomicMap();
  double mvdwc = vdw_volume(vdwmap.at(6));

  const std::map<int, double> &semap = SandersonENAtomicMap();
  double msec = semap.at(6);

  const std::map<int, double> &pemap = PaulingENAtomicMap();
  double mpec = pemap.at(6);

  const std::map<int, double> &aremap = Allred_rocow_ENAtomicMap();
  double marec = aremap.at(6);

  const std::map<int, double> &pmap = Polarizability94AtomicMap();
  double mpc = pmap.at(6);

  const std::map<int, double> &imap = ionizationEnergyAtomicMap();
  double mic = imap.at(6);

  for (const auto &atom : hmol->atoms()) {
    std::string symbol = atom->getSymbol();
    int atzi = atom->getAtomicNum();
    SZ += static_cast<double>(atzi) / zcc;
    Sm += static_cast<double>(tbl->getAtomicWeight(symbol)) / mcc;
    // need to convert radius to volume!!!
    if (vdwmap.find(atzi) != vdwmap.end()) {
      Sv += vdw_volume(vdwmap.at(atzi)) / mvdwc;
    }
    if (semap.find(atzi) != semap.end()) {
      Sse += semap.at(atzi) / msec;
    }
    if (pemap.find(atzi) != pemap.end()) {
      Spe += pemap.at(atzi) / mpec;
    }
    if (aremap.find(atzi) != aremap.end()) {
      Sare += aremap.at(atzi) / marec;
    }
    if (pmap.find(atzi) != pmap.end()) {
      Sp += pmap.at(atzi) / mpc;
    }
    if (imap.find(atzi) != imap.end()) {
      Si += imap.at(atzi) / mic;
    }
  }

  double natoms = static_cast<double>(hmol->getNumAtoms());
  MZ = SZ / natoms;      // correct
  Mm = Sm / natoms;      // correct
  Mv = Sv / natoms;      // correct trick convert radius to volume!!!
  Mse = Sse / natoms;    // correct
  Mpe = Spe / natoms;    // correct
  Mare = Sare / natoms;  // correct
  Mp = Sp / natoms;      // correct
  Mi = Si / natoms;      // correct

  return {SZ, Sm, Sv, Sse, Spe, Sare, Sp, Si,
          MZ, Mm, Mv, Mse, Mpe, Mare, Mp, Mi};
}

////// Barysz Matrixes Eigen style

// Function to compute eigenvalues and eigenvectors
}
}
}
