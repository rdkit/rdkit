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
#include <GraphMol/PeriodicTable.h>
#include "EState.h"
#include "MolSurf.h"

#include <algorithm>
#include <cmath>

namespace RDKit {
namespace Descriptors {

// These default bins were chosen using the PP3K solubility data set so that
// each bin has approximately the same number of atoms.
const std::vector<double> estateBins = {-0.390, 0.290, 0.717, 1.165, 1.540,
                                        1.807,  2.05,  4.69,  9.17,  15.0};
const std::vector<double> vsaEStateBins = {4.78, 5.00, 5.410, 5.740, 6.00,
                                           6.07, 6.45, 7.00,  11.0};

namespace {
unsigned int principalQuantumNumber(unsigned int atNum) {
  if (atNum <= 2) {
    return 1;
  } else if (atNum <= 10) {
    return 2;
  } else if (atNum <= 18) {
    return 3;
  } else if (atNum <= 36) {
    return 4;
  } else if (atNum <= 54) {
    return 5;
  } else if (atNum <= 86) {
    return 6;
  }
  return 7;
}

std::vector<double> labuteContribs(const ROMol &mol) {
  std::vector<double> Vi(mol.getNumAtoms(), 0.0);
  double hContrib = 0.0;
  getLabuteAtomContribs(mol, Vi, hContrib, true, true);
  return Vi;
}
}  // namespace

std::vector<double> calcEStateIndices(const ROMol &mol) {
  const unsigned int nAtoms = mol.getNumAtoms();
  std::vector<double> Is(nAtoms, 0.0);
  if (!nAtoms) {
    return Is;
  }
  const auto *tbl = PeriodicTable::getTable();
  for (const auto at : mol.atoms()) {
    auto d = at->getDegree();
    if (d > 0) {
      auto atNum = at->getAtomicNum();
      double dv = static_cast<int>(tbl->getNouterElecs(atNum)) -
                  static_cast<int>(at->getTotalNumHs());
      double N = principalQuantumNumber(atNum);
      Is[at->getIdx()] = (4. / (N * N) * dv + 1) / d;
    }
  }
  const double *dists = MolOps::getDistanceMat(mol, false, false, true);
  std::vector<double> accum(nAtoms, 0.0);
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      double p = dists[i * nAtoms + j] + 1;
      if (p < 1e6) {
        double tmp = (Is[i] - Is[j]) / (p * p);
        accum[i] += tmp;
        accum[j] -= tmp;
      }
    }
  }
  for (unsigned int i = 0; i < nAtoms; ++i) {
    accum[i] += Is[i];
  }
  return accum;
}

double calcMaxEStateIndex(const ROMol &mol) {
  auto vals = calcEStateIndices(mol);
  if (vals.empty()) {
    return 0.0;
  }
  return *std::max_element(vals.begin(), vals.end());
}

double calcMinEStateIndex(const ROMol &mol) {
  auto vals = calcEStateIndices(mol);
  if (vals.empty()) {
    return 0.0;
  }
  return *std::min_element(vals.begin(), vals.end());
}

double calcMaxAbsEStateIndex(const ROMol &mol) {
  auto vals = calcEStateIndices(mol);
  if (vals.empty()) {
    return 0.0;
  }
  double res = std::fabs(vals[0]);
  for (auto v : vals) {
    res = std::max(res, std::fabs(v));
  }
  return res;
}

double calcMinAbsEStateIndex(const ROMol &mol) {
  auto vals = calcEStateIndices(mol);
  if (vals.empty()) {
    return 0.0;
  }
  double res = std::fabs(vals[0]);
  for (auto v : vals) {
    res = std::min(res, std::fabs(v));
  }
  return res;
}

std::vector<double> calcEState_VSA(const ROMol &mol,
                                   const std::vector<double> *bins) {
  if (!bins) {
    bins = &estateBins;
  }
  auto propContribs = calcEStateIndices(mol);
  auto volContribs = labuteContribs(mol);
  std::vector<double> res(bins->size() + 1, 0.0);
  for (unsigned int i = 0; i < propContribs.size(); ++i) {
    auto nbin = std::upper_bound(bins->begin(), bins->end(), propContribs[i]) -
                bins->begin();
    res[nbin] += volContribs[i];
  }
  return res;
}

std::vector<double> calcVSA_EState(const ROMol &mol,
                                   const std::vector<double> *bins) {
  if (!bins) {
    bins = &vsaEStateBins;
  }
  auto propContribs = calcEStateIndices(mol);
  auto volContribs = labuteContribs(mol);
  std::vector<double> res(bins->size() + 1, 0.0);
  for (unsigned int i = 0; i < propContribs.size(); ++i) {
    auto nbin = std::upper_bound(bins->begin(), bins->end(), volContribs[i]) -
                bins->begin();
    res[nbin] += propContribs[i];
  }
  return res;
}

}  // namespace Descriptors
}  // namespace RDKit
