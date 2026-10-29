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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#ifdef RDK_HAS_EIGEN3
#include <Eigen/Eigenvalues>
#endif

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
      // the python code uses the integer value of the bond type enum
      order = static_cast<double>(bond->getBondType());
      order = bond->getBondTypeAsDouble();
      if (order != static_cast<double>(bond->getBondType()))
	std::cerr << bond->getBondTypeAsDouble() << " != " << static_cast<double>(bond->getBondType()) << std::endl;
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

namespace {
// eigenvalues of the symmetric matrix a (n x n, row major)
std::vector<double> symmetricEigenvalues(std::vector<double> a,
                                         unsigned int n) {
  if (n == 0) {
    return {};
  }
#ifdef RDK_HAS_EIGEN3
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                           Eigen::RowMajor>>
      m(a.data(), n, n);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      m, Eigen::EigenvaluesOnly);
  const auto &ev = solver.eigenvalues();
  return std::vector<double>(ev.data(), ev.data() + n);
#else
  // cyclic Jacobi rotations
  for (unsigned int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (unsigned int p = 0; p < n; ++p) {
      for (unsigned int q = p + 1; q < n; ++q) {
        off += a[p * n + q] * a[p * n + q];
      }
    }
    if (off == 0.0) {
      break;
    }
    for (unsigned int p = 0; p < n; ++p) {
      for (unsigned int q = p + 1; q < n; ++q) {
        double apq = a[p * n + q];
        if (apq == 0.0) {
          continue;
        }
        double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        double t = (theta >= 0 ? 1.0 : -1.0) /
                   (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        double c = 1.0 / std::sqrt(t * t + 1.0);
        double sn = t * c;
        for (unsigned int k = 0; k < n; ++k) {
          double akp = a[k * n + p];
          double akq = a[k * n + q];
          a[k * n + p] = c * akp - sn * akq;
          a[k * n + q] = sn * akp + c * akq;
        }
        for (unsigned int k = 0; k < n; ++k) {
          double apk = a[p * n + k];
          double aqk = a[q * n + k];
          a[p * n + k] = c * apk - sn * aqk;
          a[q * n + k] = sn * apk + c * aqk;
        }
        a[p * n + q] = a[q * n + p] = 0.0;
      }
    }
  }
  std::vector<double> res(n);
  for (unsigned int i = 0; i < n; ++i) {
    res[i] = a[i * n + i];
  }
  return res;
#endif
}
}  // namespace

double calcIpc(const ROMol &mol, bool avg, const double *distMat) {
  const unsigned int nAts = mol.getNumAtoms();
  // adjacency matrix from the topological distance matrix (distance == 1)
  std::vector<double> adj(nAts * nAts, 0.0);
  if (nAts) {
    const double *dMat =
        distMat ? distMat : MolOps::getDistanceMat(mol, false, false, true);
    for (unsigned int i = 0; i < nAts * nAts; ++i) {
      if (dMat[i] == 1.0) {
        adj[i] = 1.0;
      }
    }
  }

  // The characteristic polynomial is expanded from the eigenvalues, with the
  // roots multiplied in order of increasing magnitude, as in
  // rdkit.Chem.Graphs.CharacteristicPolynomial. The Le Verrier-Faddeev-Frame
  // recursion used before loses all precision beyond about 80 atoms.
  auto roots = symmetricEigenvalues(std::move(adj), nAts);
  std::sort(roots.begin(), roots.end(), [](double x, double y) {
    return std::fabs(x) < std::fabs(y);
  });
  std::vector<double> cPoly{1.0};
  cPoly.reserve(nAts + 1);
  for (auto r : roots) {
    cPoly.push_back(0.0);
    for (size_t i = cPoly.size() - 1; i > 0; --i) {
      cPoly[i] -= r * cPoly[i - 1];
    }
  }
  for (auto &v : cPoly) {
    v = std::fabs(v);
  }
  if (avg) {
    return RDInfoTheory::InfoEntropy(cPoly);
  }
  double sum = 0.0;
  for (auto v : cPoly) {
    sum += v;
  }
  return sum * RDInfoTheory::InfoEntropy(cPoly);
}

double calcAvgIpc(const ROMol &mol, const double *distMat) {
  return calcIpc(mol, true, distMat);
}

}  // namespace Descriptors
}  // namespace RDKit
