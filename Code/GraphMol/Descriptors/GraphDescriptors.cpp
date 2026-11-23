//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include "GraphDescriptors.h"
#include <ML/InfoTheory/InfoGainFuncs.h>
#include "CharacteristicPolynomial.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace RDKit {
namespace Descriptors {

double calcChi0(const ROMol &mol) {
  double res = 0.0;
  for (const auto atom : mol.atoms()) {
    auto deg = atom->getDegree();
    if (deg) {
      res += std::sqrt(1. / static_cast<double>(deg));
    }
  }
  return res;
}

double calcChi1(const ROMol &mol) {
  double res = 0.0;
  for (const auto bond : mol.bonds()) {
    auto c1 =
        bond->getBeginAtom()->getDegree() * bond->getEndAtom()->getDegree();
    if (c1) {
      res += std::sqrt(1. / static_cast<double>(c1));
    }
  }
  return res;
}

double calcBalabanJ(const ROMol &mol, const double *distMat) {
  const unsigned int nAts = mol.getNumAtoms();
  if (!nAts) {
    return 0.0;
  }
  const double *dMat =
      distMat ? distMat
              : MolOps::getDistanceMat(mol, true, false, true, "Balaban");

  // column sums, accumulated row by row as the python code does
  std::vector<double> s(nAts, 0.0);
  for (unsigned int i = 0; i < nAts; ++i) {
    for (unsigned int j = 0; j < nAts; ++j) {
      s[j] += dMat[i * nAts + j];
    }
  }
  std::vector<std::vector<bool>> adj(nAts, std::vector<bool>(nAts, false));
  for (const auto bond : mol.bonds()) {
    auto bi = bond->getBeginAtomIdx();
    auto ei = bond->getEndAtomIdx();
    adj[bi][ei] = true;
    adj[ei][bi] = true;
  }

  const int q = mol.getNumBonds();
  const int mu = q - static_cast<int>(nAts) + 1;

  double sum = 0.0;
  for (unsigned int i = 0; i < nAts; ++i) {
    for (unsigned int j = i; j < nAts; ++j) {
      if (adj[i][j]) {
        sum += 1. / std::sqrt(s[i] * s[j]);
      }
    }
  }
  if (mu + 1 != 0) {
    return static_cast<double>(q) / static_cast<double>(mu + 1) * sum;
  }
  return 0.0;
}

double calcBertzCT(const ROMol &mol, unsigned int cutoff,
                   const double *distMat) {
  const unsigned int numAtoms = mol.getNumAtoms();
  if (numAtoms < 2) {
    return 0.0;
  }

  // neighbor lists (sorted) and bond orders
  std::vector<std::vector<unsigned int>> nList(numAtoms);
  std::map<std::pair<unsigned int, unsigned int>, double> bondOrders;
  for (const auto bond : mol.bonds()) {
    auto a1 = bond->getBeginAtomIdx();
    auto a2 = bond->getEndAtomIdx();
    if (a1 > a2) {
      std::swap(a1, a2);
    }
     double order;
    if (bond->getIsAromatic()) {
      order = 1.5;
    } else {
      order = bond->getBondTypeAsDouble();
    }
//double order = bond->getIsAromatic() ? 1.5 : bond->getBondTypeAsDouble();
    bondOrders[std::make_pair(a1, a2)] = order;
    if (std::find(nList[a1].begin(), nList[a1].end(), a2) == nList[a1].end()) {
      nList[a1].push_back(a2);
    }
    if (std::find(nList[a2].begin(), nList[a2].end(), a1) == nList[a2].end()) {
      nList[a2].push_back(a1);
    }
  }
  for (auto &nbrs : nList) {
    std::sort(nbrs.begin(), nbrs.end());
  }
  auto lookupBondOrder = [&bondOrders](unsigned int i, unsigned int j) {
    if (i > j) {
      std::swap(i, j);
    }
    return bondOrders.at(std::make_pair(i, j));
  };

  // symmetry classes from the bond-order weighted distance matrix
  const double *bdMat =
      distMat ? distMat
              : MolOps::getDistanceMat(mol, true, false, true, "Balaban");
  std::vector<std::vector<std::string>> keysSeen;
  std::vector<int> symClasses(numAtoms, 0);
  char buf[64];
  for (unsigned int i = 0; i < numAtoms; ++i) {
    std::vector<double> row(bdMat + i * numAtoms, bdMat + (i + 1) * numAtoms);
    std::sort(row.begin(), row.end());
    if (row.size() > cutoff) {
      row.resize(cutoff);
    }
    std::vector<std::string> key;
    key.reserve(row.size());
    for (auto v : row) {
      snprintf(buf, sizeof(buf), "%.4f", v);
      key.emplace_back(buf);
    }
    auto loc = std::find(keysSeen.begin(), keysSeen.end(), key);
    auto idx = static_cast<int>(loc - keysSeen.begin());
    if (loc == keysSeen.end()) {
      keysSeen.push_back(std::move(key));
    }
    symClasses[i] = idx + 1;
  }

  // dictionaries are kept in insertion order to match the python summation
  // order. Two-element keys are stored with a leading -1.
  using ConnKey = std::tuple<int, int, int>;
  std::map<ConnKey, size_t> connectionIdx;
  std::vector<double> connectionVals;
  auto addConnection = [&](const ConnKey &key, double val) {
    auto it = connectionIdx.find(key);
    if (it == connectionIdx.end()) {
      connectionIdx[key] = connectionVals.size();
      connectionVals.push_back(val);
    } else {
      connectionVals[it->second] += val;
    }
  };
  std::map<int, size_t> atomTypeIdx;
  std::vector<double> atomTypeVals;

  for (unsigned int atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
    int hingeAtomNumber = mol.getAtomWithIdx(atomIdx)->getAtomicNum();
    auto ait = atomTypeIdx.find(hingeAtomNumber);
    if (ait == atomTypeIdx.end()) {
      atomTypeIdx[hingeAtomNumber] = atomTypeVals.size();
      atomTypeVals.push_back(1.0);
    } else {
      atomTypeVals[ait->second] += 1.0;
    }

    int hingeAtomClass = symClasses[atomIdx];
    const auto &nbrs = nList[atomIdx];
    for (size_t i = 0; i < nbrs.size(); ++i) {
      auto nbrI = nbrs[i];
      int niClass = symClasses[nbrI];
      double boI = lookupBondOrder(atomIdx, nbrI);
      if (boI > 1 && nbrI > atomIdx) {
        double numConnections = boI * (boI - 1) / 2;
        addConnection(ConnKey(-1, std::min(hingeAtomClass, niClass),
                              std::max(hingeAtomClass, niClass)),
                      numConnections);
      }
      for (size_t j = i + 1; j < nbrs.size(); ++j) {
        auto nbrJ = nbrs[j];
        int njClass = symClasses[nbrJ];
        double boJ = lookupBondOrder(atomIdx, nbrJ);
        addConnection(ConnKey(std::min(niClass, njClass), hingeAtomClass,
                              std::max(niClass, njClass)),
                      boI * boJ);
      }
    }
  }
  if (connectionVals.empty()) {
    connectionVals.push_back(1.0);
  }

  double totConnections = 0.0;
  for (auto v : connectionVals) {
    totConnections += v;
  }
  double connectionIE =
      totConnections *
      (RDInfoTheory::InfoEntropy(connectionVals) + std::log(totConnections) / std::log(2.0));
  double atomTypeIE = numAtoms * RDInfoTheory::InfoEntropy(atomTypeVals);
  return atomTypeIE + connectionIE;
}

double calcIpc(const ROMol &mol, bool avg, const double *distMat) {
  const unsigned int nAts = mol.getNumAtoms();
  // The adjacency matrix: the bonds, or the entries of distMat equal to 1. Of
  // distMat only the lower triangle, diagonal included, is read, as the
  // symmetric eigensolver used before did.
  std::vector<std::pair<unsigned int, unsigned int>> edges;
  std::vector<bool> loops(nAts, false);
  if (distMat) {
    for (unsigned int i = 0; i < nAts; ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        if (distMat[i * nAts + j] == 1.0) {
          edges.emplace_back(i, j);
        }
      }
      loops[i] = distMat[i * nAts + i] == 1.0;
    }
  } else {
    for (const auto bond : mol.bonds()) {
      edges.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    }
  }

  // The coefficients are integers, computed exactly and rounded to double
  // once; see CharacteristicPolynomial.h.
  const auto cPoly =
      detail::characteristicPolynomialMagnitudes(nAts, edges, loops);
  if (!cPoly) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (avg) {
    return RDInfoTheory::InfoEntropy(*cPoly);
  }
  double sum = 0.0;
  for (auto v : *cPoly) {
    sum += v;
  }
  return sum * RDInfoTheory::InfoEntropy(*cPoly);
}

double calcAvgIpc(const ROMol &mol, const double *distMat) {
  return calcIpc(mol, true, distMat);
}

}  // namespace Descriptors
}  // namespace RDKit
