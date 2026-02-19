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

namespace RDKit {
namespace Descriptors {
namespace Osmordred {

// option using inspired CDS paper
// https://www.sciencedirect.com/org/science/article/pii/S2635098X24001426
//  full step matrix implements

std::vector<std::vector<float>> FastWeightedStepMatrix(
    int n_atom, std::vector<std::vector<float>> &Sa,
    std::unordered_map<int, std::set<int>> &Ladj,
    std::set<std::pair<int, int>> step_xy) {
  std::vector<std::vector<float>> SF = Sa;

  // Iterate step-wise from 1 to n_atom-2
  for (int m = 1; m < n_atom - 1; ++m) {
    if (step_xy.empty()) {
      break;
    }
    std::set<std::pair<int, int>> temp_step_xy{};
    for (const auto &item : step_xy) {
      int r = std::get<0>(item);
      for (int c : Ladj[std::get<1>(item)]) {
        // Update path weight if a shorter weighted path is found
        float newWeight =
            SF[r - 1][std::get<1>(item) - 1] + Sa[std::get<1>(item) - 1][c - 1];
        if ((SF[r - 1][c - 1] == 0 || newWeight < SF[r - 1][c - 1]) && r != c) {
          SF[r - 1][c - 1] = newWeight;
          temp_step_xy.insert(std::make_pair(r, c));
        }
      }
    }
    step_xy = temp_step_xy;
  }
  return SF;
}

std::vector<std::vector<float>> initializeWeightedAdjacency(
    const ROMol &mol, const std::vector<double> &atomicProps) {
  int numAtoms = mol.getNumAtoms();
  std::vector<std::vector<float>> Sa(numAtoms,
                                     std::vector<float>(numAtoms, 0.0));

  for (const auto &bond : mol.bonds()) {
    unsigned int i = bond->getBeginAtomIdx();
    unsigned int j = bond->getEndAtomIdx();
    double pi = bond->getBondTypeAsDouble();                  // Bond order
    Sa[i][j] = 1.0 / (atomicProps[i] * atomicProps[j] * pi);  // Weighted edge
    Sa[j][i] = Sa[i][j];
  }

  return Sa;
}

std::tuple<std::vector<std::vector<float>>,
           std::unordered_map<int, std::set<int>>,
           std::set<std::pair<int, int>>>
prepareAdjacencyData(const std::vector<std::vector<float>> &Sa) {
  int numAtoms = Sa.size();
  std::unordered_map<int, std::set<int>> Ladj;
  std::set<std::pair<int, int>> step_xy;

  for (int i = 0; i < numAtoms; ++i) {
    for (int j = 0; j < numAtoms; ++j) {
      if (Sa[i][j] > 0.0) {  // Edge exists
        Ladj[i + 1].insert(j + 1);
        step_xy.insert({i + 1, j + 1});
      }
    }
  }

  return {Sa, Ladj, step_xy};
}

Eigen::MatrixXd computeBaryszMatrix0(const ROMol &mol,
                                     const std::vector<double> &atomicProps) {
  const unsigned int numAtoms = mol.getNumAtoms();

  // Initialize Eigen matrix with a large value for non-edges
  double largeValue = 1e6;
  Eigen::MatrixXd baryszMatrix =
      Eigen::MatrixXd::Constant(numAtoms, numAtoms, largeValue);

  // Set diagonal elements to 0 initially already done in
  // floydWarshallRelaxation
  for (unsigned int i = 0; i < numAtoms; ++i) {
    baryszMatrix(i, i) = 0.0;
  }

  std::vector<double> reciprocalAtomicProps(numAtoms);
  for (unsigned int i = 0; i < numAtoms; ++i) {
    reciprocalAtomicProps[i] = 1.0 / atomicProps[i];
  }

  // Fill the matrix with edge weights from bonds
  for (const auto &bond : mol.bonds()) {
    unsigned int i = bond->getBeginAtomIdx();
    unsigned int j = bond->getEndAtomIdx();

    double pi =
        bond->getBondTypeAsDouble();  // Bond order it is already normalized
    // double w = (cw * cw) / (atomicProps[i] * atomicProps[j] * pi); // already
    // normalyzed so 1 instead of cw*cw
    double w = reciprocalAtomicProps[i] * reciprocalAtomicProps[j] / pi;
    baryszMatrix(i, j) = w;
    baryszMatrix(j, i) = w;
  }

  // Floyd-Warshall algorithm for shortest paths using Eigen
  Eigen::MatrixXd baryszMatrix_ = floydWarshall(baryszMatrix);
  // can we use this instead for speed
  // https://github.com/FangyouYan/Connectivity-Stepwise-Derivation-CSD-/blob/master/cpp_code/common_utils/MsUtils.cpp

  // Update diagonal with 1.0 - C / P[i] this is 1 / normP[i] in our case
  for (unsigned int i = 0; i < numAtoms; ++i) {
    baryszMatrix_(i, i) = 1.0 - 1 / atomicProps[i];
  }

  // Replace any remaining large values with 0 to ensure no invalid entries
  baryszMatrix_ =
      (baryszMatrix_.array() == largeValue).select(0.0, baryszMatrix_);

  return baryszMatrix_;
}

Eigen::MatrixXd computeBaryszMatrix1(const ROMol &mol,
                                     const std::vector<double> &atomicProps) {
  const unsigned int numAtoms = mol.getNumAtoms();

  auto Sa = initializeWeightedAdjacency(mol, atomicProps);

  // Prepare Ladj and step_xy structures for FastFullStepMatrix but this is the
  // same for all descriptors or ???

  auto [SaMatrix, Ladj, step_xy] = prepareAdjacencyData(Sa);

  // Compute shortest path matrix using FastFullStepMatrix
  auto shortestPathMatrix =
      FastWeightedStepMatrix(mol.getNumAtoms(), SaMatrix, Ladj, step_xy);

  // Convert shortestPathMatrix to Eigen matrix
  Eigen::MatrixXd baryszMatrix = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
  for (unsigned int i = 0; i < numAtoms; ++i) {
    for (unsigned int j = i + 1; j < numAtoms; ++j) {
      baryszMatrix(i, j) = shortestPathMatrix[i][j];
      baryszMatrix(j, i) = baryszMatrix(i, j);
    }
  }
  // Update diagonal with 1.0 - (1 / P[i])
  for (unsigned int i = 0; i < numAtoms; ++i) {
    baryszMatrix(i, i) = 1.0 - (1.0 / atomicProps[i]);
  }

  return baryszMatrix;
}

// Function to calculate the Barysz matrix
Eigen::MatrixXd computeBaryszMatrix2(const ROMol &mol,
                                     const std::vector<double> &w) {
  unsigned int natom = mol.getNumAtoms();
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(natom, natom);

  std::vector<double> reciprocalw(natom);
  for (unsigned int i = 0; i < natom; ++i) {
    reciprocalw[i] = 1.0 / w[i];
  }

  // Set diagonal entries
  for (unsigned int i = 0; i < natom; ++i) {
    matrix(i, i) = 1.0 - reciprocalw[i];
  }

  // Set off-diagonal entries
  for (unsigned int i = 0; i < natom; ++i) {
    for (unsigned int j = i + 1; j < natom; ++j) {
      std::list<int> path = MolOps::getShortestPath(mol, i, j);

      for (auto it = path.begin(); std::next(it) != path.end(); ++it) {
        int atomIdx1 = *it;
        int atomIdx2 = *std::next(it);

        // const Atom* atom1 = mol.getAtomWithIdx(atomIdx1);
        // const Atom* atom2 = mol.getAtomWithIdx(atomIdx2);
        const Bond *bond = mol.getBondBetweenAtoms(atomIdx1, atomIdx2);

        double weights = reciprocalw[atomIdx1] * reciprocalw[atomIdx2];
        if (bond) {
          if (bond->getIsAromatic()) {
            matrix(i, j) += weights / 1.5;
          } else {
            switch (bond->getBondType()) {
              case Bond::SINGLE:
                matrix(i, j) += weights;
                break;
              case Bond::DOUBLE:
                matrix(i, j) += 0.5 * weights;
                break;
              case Bond::TRIPLE:
                matrix(i, j) += weights / 3.0;
                break;
              default:
                break;
            }
          }
        }
      }
      matrix(j, i) = matrix(i, j);  // Ensure symmetry
    }
  }

  return matrix;
}

std::vector<double> MatrixDescs(const ROMol &mol, Eigen::MatrixXd Mat) {
  int numAtoms = mol.getNumAtoms();

  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors;

  compute_eigenvalues_and_eigenvectors(Mat, eigenvalues, eigenvectors);

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
  // look correct per say! ve1... vr3

  double sm1 = SM1(Mat);
  double ve1 = VE1(Mat, eigenvalues, eigenvectors);
  double ve2 = VE2(Mat, numAtoms, eigenvalues, eigenvectors);
  double ve3 = VE3(Mat, numAtoms, eigenvalues, eigenvectors);
  double vr1 = VR1(Mat, bonds, eigenvalues, eigenvectors);
  double vr2 = VR2(Mat, bonds, numAtoms, eigenvalues, eigenvectors);
  double vr3 = VR3(Mat, bonds, numAtoms, eigenvalues, eigenvectors);

  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3};  // 13 outputs
}

std::vector<double> calcBaryszMatrixDescs(const ROMol &mol) {
  int method = 1;
  const PeriodicTable *tbl = PeriodicTable::getTable();

  const std::map<int, double> &vdwmap = VdWAtomicMap();
  const std::map<int, double> &semap = SandersonENAtomicMap();
  const std::map<int, double> &pemap = PaulingENAtomicMap();
  const std::map<int, double> &aremap = Allred_rocow_ENAtomicMap();
  const std::map<int, double> &pmap = Polarizability94AtomicMap();
  const std::map<int, double> &imap = ionizationEnergyAtomicMap();

  double zcc = static_cast<double>(tbl->getAtomicNumber("C"));
  double mcc = static_cast<double>(tbl->getAtomicWeight("C"));
  double mvdwc = vdw_volume(vdwmap.at(6));
  double msec = semap.at(6);
  double mpec = pemap.at(6);
  double marec = aremap.at(6);
  double mpc = pmap.at(6);
  double mic = imap.at(6);

  // Prepare vectors to store the computed atomic properties
  std::vector<double> Zs_(mol.getNumAtoms(), 0.);
  std::vector<double> ms_(mol.getNumAtoms(), 0.);
  std::vector<double> vs_(mol.getNumAtoms(), 0.);
  std::vector<double> pes_(mol.getNumAtoms(), 0.);
  std::vector<double> ses_(mol.getNumAtoms(), 0.);
  std::vector<double> ares_(mol.getNumAtoms(), 0.);
  std::vector<double> ps_(mol.getNumAtoms(), 0.);
  std::vector<double> is_(mol.getNumAtoms(), 0.);

  for (const auto &atom : mol.atoms()) {
    int atzi = atom->getAtomicNum();
    // use autonormalized atomproperties so we simplify the code as it is
    // already (P / cw not just P)
    double zValue = static_cast<double>(atzi) / zcc;
    double mValue = static_cast<double>(tbl->getAtomicWeight(atzi)) / mcc;
    double vValue = (vdwmap.find(atzi) != vdwmap.end())
                        ? vdw_volume(vdwmap.at(atzi)) / mvdwc
                        : 0.0;
    double seValue =
        (semap.find(atzi) != semap.end()) ? semap.at(atzi) / msec : 0.0;
    double peValue =
        (pemap.find(atzi) != pemap.end()) ? pemap.at(atzi) / mpec : 0.0;
    double areValue =
        (aremap.find(atzi) != aremap.end()) ? aremap.at(atzi) / marec : 0.0;
    double pValue = (pmap.find(atzi) != pmap.end()) ? pmap.at(atzi) / mpc : 0.0;
    double iValue = (imap.find(atzi) != imap.end()) ? imap.at(atzi) / mic : 0.0;

    // You can use any of the computed values (zValue, mValue, etc.)
    // Here, choose the property you want to compute the Barysz matrix for:
    Zs_[atom->getIdx()] = zValue;      //  zValue
    ms_[atom->getIdx()] = mValue;      // mass
    vs_[atom->getIdx()] = vValue;      // volume
    ses_[atom->getIdx()] = seValue;    // sanderson
    pes_[atom->getIdx()] = peValue;    // pauling
    ares_[atom->getIdx()] = areValue;  // allred
    ps_[atom->getIdx()] = pValue;      // polarisability
    is_[atom->getIdx()] = iValue;      // ioni
  }

  Eigen::MatrixXd ZBaryszMat, mBaryszMat, vBaryszMat, seBaryszMat, peBaryszMat,
      areBaryszMat, pBaryszMat, iBaryszMat;
  // no need to pass the carbon property as it is autonormalized
  if (method == 0) {
    ZBaryszMat = computeBaryszMatrix0(mol, Zs_);
    mBaryszMat = computeBaryszMatrix0(mol, ms_);
    vBaryszMat = computeBaryszMatrix0(mol, vs_);
    seBaryszMat = computeBaryszMatrix0(mol, ses_);
    peBaryszMat = computeBaryszMatrix0(mol, pes_);
    areBaryszMat = computeBaryszMatrix0(mol, ares_);
    pBaryszMat = computeBaryszMatrix0(mol, ps_);
    iBaryszMat = computeBaryszMatrix0(mol, is_);

  }

  else if (method == 1) {
    ZBaryszMat = computeBaryszMatrix1(mol, Zs_);
    mBaryszMat = computeBaryszMatrix1(mol, ms_);
    vBaryszMat = computeBaryszMatrix1(mol, vs_);
    seBaryszMat = computeBaryszMatrix1(mol, ses_);
    peBaryszMat = computeBaryszMatrix1(mol, pes_);
    areBaryszMat = computeBaryszMatrix1(mol, ares_);
    pBaryszMat = computeBaryszMatrix1(mol, ps_);
    iBaryszMat = computeBaryszMatrix1(mol, is_);

  }

  else if (method == 2) {
    ZBaryszMat = computeBaryszMatrix2(mol, Zs_);
    mBaryszMat = computeBaryszMatrix2(mol, ms_);
    vBaryszMat = computeBaryszMatrix2(mol, vs_);
    seBaryszMat = computeBaryszMatrix2(mol, ses_);
    peBaryszMat = computeBaryszMatrix2(mol, pes_);
    areBaryszMat = computeBaryszMatrix2(mol, ares_);
    pBaryszMat = computeBaryszMatrix2(mol, ps_);
    iBaryszMat = computeBaryszMatrix2(mol, is_);
  }

  std::vector<double> zBaryszDesc = MatrixDescs(mol, ZBaryszMat);
  std::vector<double> mBaryszDesc = MatrixDescs(mol, mBaryszMat);
  std::vector<double> vBaryszDesc = MatrixDescs(mol, vBaryszMat);
  std::vector<double> seBaryszDesc = MatrixDescs(mol, seBaryszMat);
  std::vector<double> peBaryszDesc = MatrixDescs(mol, peBaryszMat);
  std::vector<double> areBaryszDesc = MatrixDescs(mol, areBaryszMat);
  std::vector<double> pBaryszDesc = MatrixDescs(mol, pBaryszMat);
  std::vector<double> iBaryszDesc = MatrixDescs(mol, iBaryszMat);

  // Concatenate all descriptors vector dimension is : 8*13
  std::vector<double> concatenatedDescriptors;
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 zBaryszDesc.begin(), zBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 mBaryszDesc.begin(), mBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 vBaryszDesc.begin(), vBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 seBaryszDesc.begin(), seBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 peBaryszDesc.begin(), peBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 areBaryszDesc.begin(), areBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 pBaryszDesc.begin(), pBaryszDesc.end());
  concatenatedDescriptors.insert(concatenatedDescriptors.end(),
                                 iBaryszDesc.begin(), iBaryszDesc.end());

  return concatenatedDescriptors;
}

std::vector<std::vector<double>> computeBaryszMatrix0L(
    const ROMol &mol, const std::vector<double> &atomicProps) {
  const unsigned int numAtoms = mol.getNumAtoms();
  double largeValue = 1e6;
  // int r = test();
  // Initialize adjacency matrix
  std::vector<std::vector<double>> baryszMatrix(
      numAtoms, std::vector<double>(numAtoms, largeValue));
  for (unsigned int i = 0; i < numAtoms; ++i) {
    baryszMatrix[i][i] = 0.0;
  }

  // Compute reciprocal atomic properties
  std::vector<double> reciprocalAtomicProps(numAtoms);
  for (unsigned int i = 0; i < numAtoms; ++i) {
    reciprocalAtomicProps[i] = 1.0 / atomicProps[i];
  }

  // Populate adjacency matrix
  for (const auto &bond : mol.bonds()) {
    unsigned int i = bond->getBeginAtomIdx();
    unsigned int j = bond->getEndAtomIdx();
    double pi = bond->getBondTypeAsDouble();
    double w = reciprocalAtomicProps[i] * reciprocalAtomicProps[j] / pi;
    baryszMatrix[i][j] = w;
    baryszMatrix[j][i] = w;
  }

  // Apply Floyd-Warshall
  baryszMatrix = floydWarshallL(baryszMatrix);

  // Update diagonal
  for (unsigned int i = 0; i < numAtoms; ++i) {
    baryszMatrix[i][i] = 1.0 - 1.0 / atomicProps[i];
  }

  // Replace large values with 0
  for (unsigned int i = 0; i < numAtoms; ++i) {
    for (unsigned int j = 0; j < numAtoms; ++j) {
      if (baryszMatrix[i][j] == largeValue) {
        baryszMatrix[i][j] = 0.0;
      }
    }
  }

  return baryszMatrix;
}

std::vector<double> MatrixDescsL(const ROMol &mol,
                                 std::vector<std::vector<double>> &matrix) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> eigenvalues;
  std::vector<std::vector<double>> eigenvectors;

  compute_eigenvalues_and_eigenvectorsL(matrix, eigenvalues, eigenvectors);

  std::vector<std::pair<int, int>> bonds;
  for (const auto &bond : mol.bonds()) {
    bonds.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
  }

  double Sp_Abs = spAbsL(eigenvalues);
  double Sp_Max = spMaxL(eigenvalues);
  double Sp_Diam = spDiamL(eigenvalues);
  double Sp_Mean = spMeanL(eigenvalues);
  double Sp_AD = spADL(eigenvalues, Sp_Mean);
  double Sp_MAD = Sp_AD / numAtoms;
  double Log_EE = logEE_stable(eigenvalues);
  double sm1 = SM1L(matrix);
  double ve1 = VE1L(eigenvectors);
  double ve2 = VE2L(ve1, numAtoms);
  double ve3 = VE3L(ve1, numAtoms);
  double vr1 = VR1L(eigenvectors, bonds);
  double vr2 = VR2L(vr1, numAtoms);
  double vr3 = VR3L(vr1, numAtoms);

  return {Sp_Abs, Sp_Max, Sp_Diam, Sp_AD, Sp_MAD, Log_EE, sm1,
          ve1,    ve2,    ve3,     vr1,   vr2,    vr3};
}

std::vector<double> calcBaryszMatrixDescsL(const ROMol &mol) {
  const PeriodicTable *tbl = PeriodicTable::getTable();

  std::map<int, double> vdwmap = VdWAtomicMap();
  std::map<int, double> semap = SandersonENAtomicMap();
  std::map<int, double> pemap = PaulingENAtomicMap();
  std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
  std::map<int, double> pmap = Polarizability94AtomicMap();
  std::map<int, double> imap = ionizationEnergyAtomicMap();

  double zcc = static_cast<double>(tbl->getAtomicNumber("C"));
  double mcc = static_cast<double>(tbl->getAtomicWeight("C"));
  double mvdwc = vdw_volume(vdwmap.at(6));
  double msec = semap.at(6);
  double mpec = pemap.at(6);
  double marec = aremap.at(6);
  double mpc = pmap.at(6);
  double mic = imap.at(6);

  std::vector<double> Zs_(mol.getNumAtoms(), 0.0);
  std::vector<double> ms_(mol.getNumAtoms(), 0.0);
  std::vector<double> vs_(mol.getNumAtoms(), 0.0);
  std::vector<double> ses_(mol.getNumAtoms(), 0.0);
  std::vector<double> pes_(mol.getNumAtoms(), 0.0);
  std::vector<double> ares_(mol.getNumAtoms(), 0.0);
  std::vector<double> ps_(mol.getNumAtoms(), 0.0);
  std::vector<double> is_(mol.getNumAtoms(), 0.0);

  for (const auto &atom : mol.atoms()) {
    int atzi = atom->getAtomicNum();
    Zs_[atom->getIdx()] = static_cast<double>(atzi) / zcc;
    ms_[atom->getIdx()] = static_cast<double>(tbl->getAtomicWeight(atzi)) / mcc;
    vs_[atom->getIdx()] = (vdwmap.find(atzi) != vdwmap.end())
                              ? vdw_volume(vdwmap.at(atzi)) / mvdwc
                              : 0.0;
    ses_[atom->getIdx()] =
        (semap.find(atzi) != semap.end()) ? semap.at(atzi) / msec : 0.0;
    pes_[atom->getIdx()] =
        (pemap.find(atzi) != pemap.end()) ? pemap.at(atzi) / mpec : 0.0;
    ares_[atom->getIdx()] =
        (aremap.find(atzi) != aremap.end()) ? aremap.at(atzi) / marec : 0.0;
    ps_[atom->getIdx()] =
        (pmap.find(atzi) != pmap.end()) ? pmap.at(atzi) / mpc : 0.0;
    is_[atom->getIdx()] =
        (imap.find(atzi) != imap.end()) ? imap.at(atzi) / mic : 0.0;
  }

  std::vector<std::vector<std::vector<double>>> matrices(8);
  matrices[0] = computeBaryszMatrix0L(mol, Zs_);
  matrices[1] = computeBaryszMatrix0L(mol, ms_);
  matrices[2] = computeBaryszMatrix0L(mol, vs_);
  matrices[3] = computeBaryszMatrix0L(mol, ses_);
  matrices[4] = computeBaryszMatrix0L(mol, pes_);
  matrices[5] = computeBaryszMatrix0L(mol, ares_);
  matrices[6] = computeBaryszMatrix0L(mol, ps_);
  matrices[7] = computeBaryszMatrix0L(mol, is_);

  std::vector<double> concatenatedDescriptors;
  for (auto &matrix : matrices) {
    auto desc = MatrixDescsL(mol, matrix);
    concatenatedDescriptors.insert(concatenatedDescriptors.end(), desc.begin(),
                                   desc.end());
  }

  return concatenatedDescriptors;
}

namespace {

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
}
/// TopologicalCharge
void addToQueries(
    std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> &queries,
    const std::string &key, std::shared_ptr<RWMol> mol) {
  for (auto &pair : queries) {
    if (pair.first == key) {  // Key already exists, update the value
      pair.second = mol;
      return;
    }
  }
  queries.emplace_back(key, mol);  // Key not found, add a new pair
}

// Precompile SMARTS patterns for efficiency
static const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> &
GetHsQueries() {
  static const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>>
      hsQueries = [] {
        std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> queries;
        for (const auto &entry : hsPatterns) {
          auto mol = SmartsToMol(entry.second);
          if (mol) {
            addToQueries(queries, entry.first, std::shared_ptr<RWMol>(mol));
          } else {
            BOOST_LOG(rdWarningLog) << "Invalid SMARTS: " << entry.second << std::endl;
          }
        }
        return queries;
      }();
  return hsQueries;
}

// Precompile SMARTS patterns for efficiency
static const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> &
GetesQueries() {
  static const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>>
      esQueries = [] {
        std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> queries;
        for (const auto &entry : esPatterns) {
          auto mol = SmartsToMol(entry.second);
          if (mol) {
            addToQueries(queries, entry.first, std::shared_ptr<RWMol>(mol));
          } else {
            BOOST_LOG(rdWarningLog) << "Invalid SMARTS: " << entry.second << std::endl;
          }
        }
        return queries;
      }();
  return esQueries;
}

// Precompile SMARTS patterns for efficiency
static const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> &
GetesExtQueries() {
  static const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>>
      esExtQueries = [] {
        std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> queries;
        for (const auto &entry : esPatternsFromOEState) {
          auto mol = SmartsToMol(entry.second);
          if (mol) {
            addToQueries(queries, entry.first, std::shared_ptr<RWMol>(mol));
          } else {
            BOOST_LOG(rdWarningLog) << "Invalid SMARTS: " << entry.second << std::endl;
          }
        }
        return queries;
      }();
  return esExtQueries;
}

// Function to switch between standard and extended SMARTS queries
const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>> &
getEStateQueries(bool extended) {
  return extended ? GetesExtQueries() : GetesQueries();
}

// Function to calculate EState fingerprints
std::vector<double> calcEStateDescs(const ROMol &mol, bool extended) {
  const auto &queries = getEStateQueries(extended);
  // const std::vector<std::pair<std::string, std::string>> esPat = extended ?
  // esPatternsFromOEState : esPatterns;
  size_t nPatts = queries.size();
  std::vector<int> counts(nPatts, 0);
  std::vector<double> sums(nPatts, 0.0);
  std::vector<double> maxValues(nPatts, std::numeric_limits<double>::lowest());
  std::vector<double> minValues(nPatts, std::numeric_limits<double>::max());
  // Parse SMARTS strings into RDKit molecule objects

  // Calculate EState indices for the molecule
  std::vector<double> esIndices = calcEStateIndices(mol);

  size_t i = 0;

  for (const auto &[name, pattern] : queries) {
    // Find all substructure matches
    std::vector<MatchVectType> matches;
    SubstructMatch(mol, *pattern, matches, true);

    // Update counts, sums, max, and min
    counts[i] = static_cast<int>(matches.size());
    for (const auto &match : matches) {
      int atomIdx = match[0].second;  // Atom index from the match
      double value = esIndices[atomIdx];
      sums[i] += value;
      maxValues[i] = std::max(maxValues[i], value);
      minValues[i] = std::min(minValues[i], value);
    }

    // Handle cases where there are no matches
    if (counts[i] == 0) {
      maxValues[i] = 0.0;
      minValues[i] = 0.0;
    }

    ++i;  // Increment the index
  }

  // Concatenate counts, sums, maxValues, and minValues into a single vector
  std::vector<double> results;
  results.reserve(4 * nPatts);
  results.insert(results.end(), counts.begin(), counts.end());  // Counts
  results.insert(results.end(), sums.begin(), sums.end());      // Sums
  results.insert(results.end(), maxValues.begin(),
                 maxValues.end());  // Max values
  results.insert(results.end(), minValues.begin(),
                 minValues.end());  // Min values

  return results;
}

std::vector<double> calcHBDHBAtDescs(const ROMol &mol,
                                     const std::vector<double> &esIndices) {
  // Indices for HBD, wHBD, HBA, and wHBA patterns
  // const std::vector<int> HBD = {10, 21, 23, 24, 25, 34, 48};
  // const std::vector<int> wHBD = {38, 54};
  // const std::vector<int> HBA = {21, 23, 24, 25, 29, 30, 34, 35, 36, 37, 38,
  // 50, 51, 54, 70}; const std::vector<int> wHBA = {18, 10, 11, 12, 14, 15, 16,
  // 17, 18};

  // Indices for HBD, wHBD, HBA, and wHBA patterns
  const std::vector<std::string> HBD = {"dNH",   "sNH2", "ssNH2", "aaNH",
                                        "sssNH", "ddsN", "aasN"};
  const std::vector<std::string> wHBD = {"sssN", "ssssN"};
  const std::vector<std::string> HBA = {
      "sOH",      "ssO",       "dO",       "aaO",       "ssO(ester)",
      "dO(acid)", "dO(ester)", "dO(amid)", "dO(nitro)", "dO(sulfo)"};
  const std::vector<std::string> wHBA = {"dO", "sOH", "aaO"};

  // Initialize accumulators for the groups
  int nHBd = 0, nwHBd = 0, nHBa = 0, nwHBa = 0;
  double SHBd = 0.0, SwHBd = 0.0, SHBa = 0.0, SwHBa = 0.0;

  // Function to find index of a SMARTS pattern in `esQueries`
  auto findIndex = [&](const std::string &smarts) -> int {
    auto esQueries = GetesQueries();
    for (size_t i = 0; i < esQueries.size(); ++i) {
      if (esQueries[i].first == smarts) return i;  // Found index
    }
    return -1;  // Not found
  };

  // Helper function to process matches
  auto processMatches = [&](const std::vector<std::string> &patterns,
                            int &count, double &sum) {
    for (const auto &smarts : patterns) {
      int idx = findIndex(smarts);
      if (idx == -1 || !GetesQueries()[idx].second)
        continue;  // Skip if not found or null pointer

      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *GetesQueries()[idx].second, matches, true);

      for (const auto &match : matches) {
        int atomIdx = match[0].second;
        double value = esIndices[atomIdx];
        count++;
        sum += value;
      }
    }
  };

  /*
  for (const auto& smarts : HBD) {
      auto it = esQueries.find(smarts);
      if (it == esQueries.end() || !it->second) continue;

      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *it->second, matches, true);

      for (const auto& match : matches) {
          int atomIdx = match[0].second;
          double value = esIndices[atomIdx];
          nHBd++;
          SHBd += value;
      }
  }

  for (const auto& smarts : wHBD) {
      auto it = esQueries.find(smarts);
      if (it == esQueries.end() || !it->second) continue;

      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *it->second, matches, true);

      for (const auto& match : matches) {
          int atomIdx = match[0].second;
          double value = esIndices[atomIdx];
          nwHBd++;
          SwHBd += value;
      }
  }

  for (const auto& smarts : HBA) {
      auto it = esQueries.find(smarts);
      if (it == esQueries.end() || !it->second) continue;

      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *it->second, matches, true);

      for (const auto& match : matches) {
          int atomIdx = match[0].second;
          double value = esIndices[atomIdx];
          nHBa++;
          SHBa += value;
      }
  }

  for (const auto& smarts : wHBA) {
      auto it = esQueries.find(smarts);
      if (it == esQueries.end() || !it->second) continue;

      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *it->second, matches, true);

      for (const auto& match : matches) {
          int atomIdx = match[0].second;
          double value = esIndices[atomIdx];
          nwHBa++;
          SwHBa += value;
      }
  }


  // Combine results into a single vector
  std::vector<double> results = {
      static_cast<double>(nHBd), SHBd,
      static_cast<double>(nwHBd), SwHBd,
      static_cast<double>(nHBa), SHBa,
      static_cast<double>(nwHBa), SwHBa
  };



  return results;
  */

  // Process each descriptor group
  processMatches(HBD, nHBd, SHBd);
  processMatches(wHBD, nwHBd, SwHBd);
  processMatches(HBA, nHBa, SHBa);
  processMatches(wHBA, nwHBa, SwHBa);

  // Combine results into a single vector
  return {static_cast<double>(nHBd), SHBd, static_cast<double>(nwHBd), SwHBd,
          static_cast<double>(nHBa), SHBa, static_cast<double>(nwHBa), SwHBa};
}

// Function to calculate HEState fingerprints + need to add the HBD, wHDBm HBA
// and wHBA patterns
std::vector<double> calcHEStateDescs(const ROMol &mol) {
  auto hsQueries = GetHsQueries();
  size_t nPatts = hsQueries.size();
  std::vector<int> counts(nPatts, 0);
  std::vector<double> sums(nPatts, 0.0);
  std::vector<double> maxValues(nPatts, 0.);
  std::vector<double> minValues(nPatts, 0.);
  // Parse SMARTS strings into RDKit molecule objects

  // Calculate EState indices for the molecule
  std::vector<double> hesIndices = CalcHEStateIndices(mol);

  int i = 0;
  for (const auto &[name, pattern] : hsQueries) {
    if (!pattern) continue;  // Skip invalid SMARTS patterns

    // Find all substructure matches
    std::vector<MatchVectType> matches;
    SubstructMatch(mol, *pattern, matches, true);

    // Update counts, sums, max, and min
    counts[i] = static_cast<int>(matches.size());
    for (const auto &match : matches) {
      int atomIdx = match[0].second;  // Atom index from the match
      double value = hesIndices[atomIdx];
      sums[i] += value;
      maxValues[i] = std::max(maxValues[i], value);
      minValues[i] = std::min(minValues[i], value);
    }

    // Handle cases where there are no matches
    if (counts[i] == 0) {
      maxValues[i] = 0.0;
      minValues[i] = 0.0;
    }
    i++;
  }

  // Concatenate counts, sums, maxValues, and minValues into a single vector
  std::vector<double> esIndices = CalcHEStateIndices(mol);
  std::vector<double> HDADescs = calcHBDHBAtDescs(mol, esIndices);

  std::vector<double> results;
  results.reserve(4 * nPatts + HDADescs.size());

  results.insert(results.end(), counts.begin(), counts.end());  // Counts
  results.insert(results.end(), sums.begin(), sums.end());      // Sums
  results.insert(results.end(), maxValues.begin(),
                 maxValues.end());  // Max values
  results.insert(results.end(), minValues.begin(),
                 minValues.end());  // Min values
  results.insert(results.end(), HDADescs.begin(),
                 HDADescs.end());  // HDADescs values

  return results;
}

// Chi computation Chain dv and d from 3 to 7 orders
std::vector<double> calcChichain(const ROMol &mol) {
  // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC
  // 1-10 and Total
  std::vector<double> results(10, 0.0);

  for (int order = 3; order <= 7; ++order) {
    auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
    double xd = 0.0, xdv = 0.0;
    for (const auto &[path, nodes, type] : classifiedPaths) {
      double cd = 1.0, cdv = 1.0;
      if (type == ChiType::Chain) {
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        xd += 1. / std::sqrt(cd);
        xdv += 1. / std::sqrt(cdv);
      }
    }
    results[order - 3] = xd;
    results[5 + order - 3] = xdv;
  }

  return results;
}

// Chi computation CLuster ie c for  dv and d from 3 to 6 orders
std::vector<double> calcChicluster(const ROMol &mol) {
  // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC
  // 1-10 and Total
  std::vector<double> results(8, 0.0);

  for (int order = 3; order <= 6; ++order) {
    auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
    double xd = 0.0, xdv = 0.0;
    for (const auto &[path, nodes, type] : classifiedPaths) {
      double cd = 1.0, cdv = 1.0;
      if (type == ChiType::Cluster) {
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        xd += 1. / std::sqrt(cd);
        xdv += 1. / std::sqrt(cdv);
      }
    }
    results[order - 3] = xd;
    results[4 + order - 3] = xdv;
  }

  return results;
}

// Chi computation Chain dv and d from 3 to 7 orders
std::vector<double> calcChipathcluster(const ROMol &mol) {
  // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC
  // 1-10 and Total
  std::vector<double> results(6, 0.0);

  for (int order = 4; order <= 6; ++order) {
    auto classifiedPaths = extractAndClassifyPaths(mol, order, false);
    double xd = 0.0, xdv = 0.0;
    for (const auto &[path, nodes, type] : classifiedPaths) {
      double cd = 1.0, cdv = 1.0;
      if (type == ChiType::PathCluster) {
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        xd += 1. / std::sqrt(cd);
        xdv += 1. / std::sqrt(cdv);
      }
    }
    results[order - 4] = xd;
    results[3 + order - 4] = xdv;
  }

  return results;
}

// Chi computation Chain dv and d from 3 to 7 orders
std::vector<double> calcChipath(const ROMol &mol) {
  // Output vector: First 10 elements are MPC 2-10 and Total, next 11 are piPC
  // 1-10 and Total
  std::vector<double> results(32, 0.0);

  double xd = 0.0, xdv = 0.0;

  for (const auto &at : mol.atoms()) {
    double cd = 1.0, cdv = 1.0;

    double d = getSigmaElectrons(*at);  // d
    cd *= d;
    double dv = getValenceElectrons(*at);  // dv
    cdv *= dv;

    xd += 1. / std::sqrt(cd);
    xdv += 1. / std::sqrt(cdv);
  }

  results[0] = xd;
  results[8] = xd / mol.getNumAtoms();
  results[16] = xdv;
  results[24] = xdv / mol.getNumAtoms();

  for (int order = 1; order <= 7; ++order) {
    auto classifiedPaths = extractAndClassifyPaths(mol, order, false);

    double xd = 0.0, xdv = 0.0;

    int nbnodes = 0;
    for (const auto &[bonds, nodes, type] : classifiedPaths) {
      double cd = 1.0, cdv = 1.0;
      if (type == ChiType::Path) {
        nbnodes += 1;
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        xd += 1. / std::sqrt(cd);
        xdv += 1. / std::sqrt(cdv);
      }
    }
    results[order] = xd;
    results[8 + order] = xd / nbnodes;
    results[16 + order] = xdv;
    results[24 + order] = xdv / nbnodes;
  }

  return results;
}

std::vector<double> calcAllChiDescriptors(const ROMol &mol) {
  std::vector<double> results(56,
                              0.0);  // Full results vector for all descriptors

  // Precompute sigma and valence electrons for all atoms
  std::vector<double> sigmaElectrons(mol.getNumAtoms(), 0.0);
  std::vector<double> valenceElectrons(mol.getNumAtoms(), 0.0);

  double path_0_xd = 0.0, path_0_xdv = 0.0;

  for (const auto &atom : mol.atoms()) {
    auto idx = atom->getIdx();
    sigmaElectrons[idx] = getSigmaElectrons(*atom);
    valenceElectrons[idx] = getValenceElectrons(*atom);
    path_0_xd += 1.0 / std::sqrt(sigmaElectrons[idx]);
    path_0_xdv += 1.0 / std::sqrt(valenceElectrons[idx]);
  }

  // Path-specific accumulators for totals and averages
  double path_xd_total = 0.0, path_xdv_total = 0.0;
  std::vector<double> path_xd(8, 0.0), path_xdv(8, 0.0);
  std::vector<int> path_node_counts(8, 0);

  // Store Path results for order 0
  results[24] = path_0_xd;                       // Total Xp-0d
  results[32] = path_0_xd / mol.getNumAtoms();   // Average Xp-0d
  results[40] = path_0_xdv;                      // Total Xp-0dv
  results[48] = path_0_xdv / mol.getNumAtoms();  // Average Xp-0dv

  for (int order = 1; order <= 7; ++order) {
    auto classifiedPaths = extractAndClassifyPaths(mol, order, false);

    double chain_xd = 0.0, chain_xdv = 0.0;
    double cluster_xd = 0.0, cluster_xdv = 0.0;
    double pathcluster_xd = 0.0, pathcluster_xdv = 0.0;

    for (const auto &[bonds, nodes, type] : classifiedPaths) {
      if (type == ChiType::Chain && order >= 3 && order <= 7) {
        double cd = 1.0, cdv = 1.0;
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        chain_xd += 1.0 / std::sqrt(cd);
        chain_xdv += 1.0 / std::sqrt(cdv);
      } else if (type == ChiType::Cluster && order >= 3 && order <= 6) {
        double cd = 1.0, cdv = 1.0;
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        cluster_xd += 1.0 / std::sqrt(cd);
        cluster_xdv += 1.0 / std::sqrt(cdv);
      } else if (type == ChiType::PathCluster && order >= 4 && order <= 6) {
        double cd = 1.0, cdv = 1.0;
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        pathcluster_xd += 1.0 / std::sqrt(cd);
        pathcluster_xdv += 1.0 / std::sqrt(cdv);
      } else if (type == ChiType::Path) {
        double cd = 1.0, cdv = 1.0;
        for (const auto &node : nodes) {
          const Atom *at = mol.getAtomWithIdx(node);
          double d = getSigmaElectrons(*at);  // d
          cd *= d;
          double dv = getValenceElectrons(*at);  // dv
          cdv *= dv;
        }
        path_node_counts[order] += 1;
        path_xd[order] += 1.0 / std::sqrt(cd);
        path_xdv[order] += 1.0 / std::sqrt(cdv);
      }
    }

    // Update total Path values
    if (order <= 7) {
      path_xd_total += path_xd[order];
      path_xdv_total += path_xdv[order];
    }

    // Store Chain results (order 3-7)
    if (order >= 3 && order <= 7) {
      results[order - 3] = chain_xd;
      results[5 + (order - 3)] = chain_xdv;
    }

    // Store Cluster results (order 3-6)
    if (order >= 3 && order <= 6) {
      results[10 + (order - 3)] = cluster_xd;
      results[14 + (order - 3)] = cluster_xdv;
    }

    // Store PathCluster results (order 4-6)
    if (order >= 4 && order <= 6) {
      results[18 + (order - 4)] = pathcluster_xd;
      results[21 + (order - 4)] = pathcluster_xdv;
    }
  }

  // Finalize Path results
  for (int order = 1; order <= 7; ++order) {
    results[24 + order] = path_xd[order];  // Total path xd
    results[32 + order] =
        path_xd[order] / path_node_counts[order];  // Average path xd
    results[40 + order] = path_xdv[order];         // Total path xdv
    results[48 + order] =
        path_xdv[order] / path_node_counts[order];  // Average path xdv
  }

  return results;
}

///////
// Define the graph as an adjacency list
using Graph = std::unordered_map<int, std::vector<std::pair<int, double>>>;

// Build the molecular graph
Graph buildGraph(const ROMol &mol) {
  Graph graph;
  for (const auto &bond : mol.bonds()) {
    int start = bond->getBeginAtomIdx();
    int end = bond->getEndAtomIdx();
    double weight = static_cast<double>(mol.getAtomWithIdx(start)->getDegree() *
                                        mol.getAtomWithIdx(end)->getDegree());

    graph[start].emplace_back(end, weight);
    graph[end].emplace_back(start, weight);
  }
  return graph;
}

// Recursive DFS for atomic ID computation
double computeAtomicId(const Graph &graph, int atomIdx, double epsilon,
                       double currentWeight, std::unordered_set<int> &visited,
                       double limit) {
  double id = 0.0;

  visited.insert(atomIdx);

  // Graphs can have single atoms
  auto res = graph.find(atomIdx);
  if ( res == graph.end() ) return id;
  
  for (const auto &[nextAtom, edgeWeight] : res->second) {
    if (visited.count(nextAtom)) continue;

    double combinedWeight = currentWeight * edgeWeight;

    id += 1.0 / std::sqrt(combinedWeight);  // Contribution to ID

    if (combinedWeight < limit) {
      id += computeAtomicId(graph, nextAtom, epsilon, combinedWeight, visited,
                            limit);
    }
  }

  visited.erase(atomIdx);  // Backtrack
  return id;
}

// Compute atomic IDs for all atoms in the molecule
std::vector<double> computeAtomicIds(const ROMol &mol, double epsilon) {
  int natoms = mol.getNumAtoms();
  std::vector<double> atomicIds(natoms, 0.0);
  if (natoms > 1) {
    Graph graph = buildGraph(mol);
    double limit = 1.0 / (epsilon * epsilon);

    for (int atomIdx = 0; atomIdx < natoms; ++atomIdx) {
      std::unordered_set<int> visited;
      double id = computeAtomicId(graph, atomIdx, epsilon, 1.0, visited, limit);
      atomicIds[atomIdx] = 1.0 + id / 2.0;  // Normalize
    }
  }
  return atomicIds;
}

std::vector<double> calcMolecularId(const ROMol &mol) {
  double epsilon = 1e-10;
  std::vector<double> atomicIds = computeAtomicIds(mol, epsilon);
  size_t numAtoms = mol.getNumAtoms();

  // 12 results:
  // MID (all), AMID (all averaged), MID_h, AMID_h, MID_C, AMID_C, MID_N,
  // AMID_N, MID_O, AMID_O, MID_X, AMID_X
  std::vector<double> results(12, 0.0);

  for (size_t i = 0; i < atomicIds.size(); ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    int atomicNum = atom->getAtomicNum();

    // All atoms
    results[0] += atomicIds[i];             // MID
    if (atomicNum > 1 && atomicNum != 6) {  // Heavy atoms
      results[2] += atomicIds[i];           // MID_h
    }

    // Specific element sums
    if (atomicNum == 6) results[4] += atomicIds[i];  // MID_C
    if (atomicNum == 7) results[6] += atomicIds[i];  // MID_N
    if (atomicNum == 8) results[8] += atomicIds[i];  // MID_O
    if (atomicNum == 9 || atomicNum == 17 || atomicNum == 35 ||
        atomicNum == 53 || atomicNum == 85) {
      // Halogens: F, Cl, Br, I, At
      results[10] += atomicIds[i];  // MID_X
    }
  }

  // Compute averages of each type
  results[1] = results[0] / numAtoms;    // AMID
  results[3] = results[2] / numAtoms;    // AMID_h
  results[5] = results[4] / numAtoms;    // AMID_C
  results[7] = results[6] / numAtoms;    // AMID_N
  results[9] = results[8] / numAtoms;    // AMID_O
  results[11] = results[10] / numAtoms;  // AMID_X

  return results;
}

// v2.0: Control function to check if Gasteiger parameters exist for all atoms
// BEFORE calling computeGasteigerCharges. This prevents crashes on unusual
// elements like Hg (mercury) or certain phosphorus environments.
//
// Returns:
//   - true: All atoms have parameters -> safe to call computeGasteigerCharges
//   - false: At least one atom lacks parameters -> return NaN instead
bool checkGasteigerParameters(const ROMol &mol) {
  PeriodicTable *table = PeriodicTable::getTable();
  const GasteigerParams *params = GasteigerParams::getParams();

  // Check each atom: if ONE atom doesn't have parameters, the ENTIRE
  // calculation is invalid
  for (auto &atom : mol.atoms()) {
    std::string elem = table->getElementSymbol(atom->getAtomicNum());
    std::string mode;

    // Reproduce EXACTLY the logic from computeGasteigerCharges to determine
    // mode
    switch (atom->getHybridization()) {
      case Atom::SP3:
        mode = "sp3";
        break;
      case Atom::SP2:
        mode = "sp2";
        break;
      case Atom::SP:
        mode = "sp";
        break;
      default:
        // For atoms with unspecified hybridization, determine mode based on
        // element
        if (atom->getAtomicNum() == 1) {
          mode = "*";  // Hydrogen always uses "*" mode
        } else if (atom->getAtomicNum() == 16) {
          // Sulfur: check oxygen neighbors
          ROMol::ADJ_ITER nbrIdx, endIdx;
          boost::tie(nbrIdx, endIdx) = mol.getAtomNeighbors(atom);
          int no = 0;
          while (nbrIdx != endIdx) {
            if (mol.getAtomWithIdx(*nbrIdx)->getAtomicNum() == 8) {
              no++;
            }
            nbrIdx++;
          }
          if (no == 2) {
            mode = "so2";
          } else if (no == 1) {
            mode = "so";
          } else {
            mode = "sp3";
          }
        } else {
          mode = "sp3";  // Default
        }
    }

    // Check if parameters exist for this element + mode
    try {
      params->getParams(elem, mode, true);
    } catch (...) {
      // Parameters don't exist (e.g., Hg, certain P environments)
      return false;
    }
  }

  return true;
}

std::vector<double> calcRNCG_RPCG(const ROMol &mol) {
  // subclass of CPSA using only 2D descriptors available in Mordred v1
  // v2.0: Check Gasteiger parameters first (on original mol, before adding H)
  if (!checkGasteigerParameters(mol)) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }

  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  double maxpos = 0;
  double maxneg = 0;
  double totalneg = 0;
  double totalpos = 0;

  // v2.0: Try-catch as safety net
  try {
    computeGasteigerCharges(*hmol, 12, true);
  } catch (...) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }

  for (auto &atom : hmol->atoms()) {
    double ch = atom->getProp<double>(common_properties::_GasteigerCharge);

    if (atom->hasProp(common_properties::_GasteigerHCharge)) {
      ch += atom->getProp<double>(common_properties::_GasteigerHCharge);
    }

    if (ch < 0) {
      totalneg += -ch;
      if (-ch > maxneg) {
        maxneg = -ch;
      }

    } else if (ch > 0) {
      totalpos += ch;
      if (ch > maxpos) {
        maxpos = ch;
      }
    }
  }
  if (totalneg == 0 || totalpos == 0) {
    return {0., 0.};
  }
  return {maxneg / totalneg, maxpos / totalpos};
}

// Function to compute BCUT descriptors for multiple properties
std::vector<double> calcBCUTs(const ROMol &mol) {
  // v2.0: Check Gasteiger parameters first
  bool gasteiger_ok = checkGasteigerParameters(mol);

  // Atomic properties to compute
  auto *tbl = PeriodicTable::getTable();
  std::map<int, double> vdwMap = VdWAtomicMap();
  std::map<int, double> sandersonENMap = SandersonENAtomicMap();
  std::map<int, double> paulingENMap = PaulingENAtomicMap();
  std::map<int, double> allredENMap = Allred_rocow_ENAtomicMap();
  std::map<int, double> polarizabilityMap = Polarizability94AtomicMap();
  std::map<int, double> ionizationMap = ionizationEnergyAtomicMap();

  size_t numAtoms = mol.getNumAtoms();
  std::vector<double> gasteigerCharges(numAtoms, 0.0);

  // v2.0: Only compute Gasteiger if parameters exist
  if (gasteiger_ok) {
    try {
      computeGasteigerCharges(mol, 12, true);
      for (size_t i = 0; i < numAtoms; ++i) {
        const auto *atom = mol.getAtomWithIdx(i);
        double ch = atom->getProp<double>(common_properties::_GasteigerCharge);
        if (atom->hasProp(common_properties::_GasteigerHCharge)) {
          ch += atom->getProp<double>(common_properties::_GasteigerHCharge);
        }
        gasteigerCharges[i] = ch;
      }
    } catch (...) {
      // Gasteiger failed - charges stay at 0.0
    }
  }
  // If gasteiger_ok is false, gasteigerCharges remains all zeros

  // Vector to store BCUT results
  std::vector<double> results;
  results.reserve(24);

  // List of atomic property vectors
  std::vector<std::vector<double>> atomicProperties(
      12, std::vector<double>(numAtoms, 0.0));

  // Populate atomic property vectors
  for (size_t i = 0; i < numAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    int atomNumber = atom->getAtomicNum();
    atomicProperties[0][i] = gasteigerCharges[i];  // Gasteiger charge (c)
    atomicProperties[1][i] =
        getValenceElectrons(*atom);  // Valence electrons (dv)
    atomicProperties[2][i] = getSigmaElectrons(*atom);  // Sigma electrons (d)
    atomicProperties[3][i] = getIntrinsicState(*atom);  // Intrinsic state (s)
    atomicProperties[4][i] =
        static_cast<double>(atomNumber);  // Atomic number (Z)
    atomicProperties[5][i] =
        tbl->getAtomicWeight(atomNumber);  // Atomic weight (m)
    atomicProperties[6][i] =
        vdw_volume(vdwMap[atomNumber]);  // Van der Waals volume need vdw_volume
                                         // to go from (r) to (v)
    atomicProperties[7][i] =
        sandersonENMap[atomNumber];  // Sanderson electronegativity (se)
    atomicProperties[8][i] =
        paulingENMap[atomNumber];  // Pauling electronegativity (pe)
    atomicProperties[9][i] =
        allredENMap[atomNumber];  // Allred-Rocow electronegativity (are)
    atomicProperties[10][i] =
        polarizabilityMap[atomNumber];  // Polarizability (p)
    atomicProperties[11][i] =
        ionizationMap[atomNumber];  // Ionization energy (i)
  }

  for (auto &result :
       BCUT2D(mol, atomicProperties, BCUTOptions::BURDEN_MATRIX)) {
    results.push_back(result.first);
    results.push_back(result.second);
  }

  return results;
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
// c    gasteiger charge

// autocorrelation

// Function to compute the ATS descriptors
std::vector<double> calcAutoCorrelationEigen(const ROMol &mol) {
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));

  double *dist = MolOps::getDistanceMat(*hmol, false);  // Topological matrix
  unsigned int numAtoms = hmol->getNumAtoms();

  // Lookup tables
  std::map<int, double> vdwmap = VdWAtomicMap();
  std::map<int, double> semap = SandersonENAtomicMap();
  std::map<int, double> pemap = PaulingENAtomicMap();
  std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
  std::map<int, double> pmap = Polarizability94AtomicMap();
  std::map<int, double> imap = ionizationEnergyAtomicMap();
  const auto *tbl = PeriodicTable::getTable();

  // Eigen vector for atomic properties
  Eigen::MatrixXd propertyMatrix(12, numAtoms);
  propertyMatrix.setZero();

  // Compute Gasteiger charges
  computeGasteigerCharges(*hmol, 12, true);
  std::vector<double> gasteigerCharges(numAtoms, 0.0);

  for (unsigned int i = 0; i < numAtoms; ++i) {
    const auto *atom = hmol->getAtomWithIdx(i);
    double charge = atom->getProp<double>(common_properties::_GasteigerCharge);
    if (atom->hasProp(common_properties::_GasteigerHCharge)) {
      charge += atom->getProp<double>(common_properties::_GasteigerHCharge);
    }
    gasteigerCharges[i] = charge;

    int atomNumber = atom->getAtomicNum();
    propertyMatrix(0, i) = getValenceElectrons(*atom);
    propertyMatrix(1, i) = getSigmaElectrons(*atom);
    propertyMatrix(2, i) = getIntrinsicState(*atom);
    propertyMatrix(3, i) = static_cast<double>(atomNumber);
    propertyMatrix(4, i) = tbl->getAtomicWeight(atomNumber);
    propertyMatrix(5, i) = vdw_volume(vdwmap[atomNumber]);
    propertyMatrix(6, i) = semap[atomNumber];
    propertyMatrix(7, i) = pemap[atomNumber];
    propertyMatrix(8, i) = aremap[atomNumber];
    propertyMatrix(9, i) = pmap[atomNumber];
    propertyMatrix(10, i) = imap[atomNumber];
    propertyMatrix(11, i) = gasteigerCharges[i];  // need to change order ...
  }

  // Initialize the topological distance matrix as an Eigen matrix
  Eigen::MatrixXd distanceMatrix(numAtoms, numAtoms);
  for (unsigned int i = 0; i < numAtoms; ++i) {
    for (unsigned int j = 0; j < numAtoms; ++j) {
      distanceMatrix(i, j) = dist[i * numAtoms + j];
    }
  }

  // Compute the ATS descriptors
  const int maxDistance = 9;
  std::vector<double> ATSDescriptors(11 * maxDistance,
                                     0.0);  // not charges descriptors
  std::vector<double> AATSDescriptors(11 * maxDistance,
                                      0.0);  // not charges descriptors
  std::vector<double> ATSCDescriptors(12 * maxDistance,
                                      0.0);  // we have now Charges too
  std::vector<double> AATSCDescriptors(12 * maxDistance,
                                       0.0);  // we have now Charges too
  std::vector<double> MATSDescriptors(
      12 * (maxDistance - 1),
      0.0);  // we have now Charges too but not zeros order
  std::vector<double> GATSDescriptors(
      12 * (maxDistance - 1),
      0.0);  // we have now Charges too but not zeros order

  // ATSC for centered ie avec - avec.mean()

  for (int k = 0; k < maxDistance; ++k) {
    Eigen::MatrixXd binaryMatrix =
        (distanceMatrix.array() == (k)).cast<double>();
    double gsum = binaryMatrix.array().sum();
    for (int t = 0; t < 12; ++t) {
      Eigen::VectorXd weights = propertyMatrix.row(t).transpose();
      Eigen::VectorXd weights_centered = weights.array() - weights.mean();
      Eigen::MatrixXd weights_col = weights.replicate(1, weights.size());
      Eigen::MatrixXd weights_row =
          weights.transpose().replicate(weights.size(), 1);
      Eigen::MatrixXd diffSquared =
          (weights_row.array() - weights_col.array()).square();

      double gatsdenominator = weights_centered.array().square().sum() /
                               (weights_centered.size() - 1);
      if (k > 0) {
        if (t < 11) {
          ATSDescriptors[t * maxDistance + k] =
              0.5 * (weights.transpose() * binaryMatrix * weights).value();
          AATSDescriptors[t * maxDistance + k] =
              ATSDescriptors[t * maxDistance + k] / (0.5 * gsum);
          ATSCDescriptors[(t + 1) * maxDistance + k] =
              0.5 *
              (weights_centered.transpose() * binaryMatrix * weights_centered)
                  .value();
          AATSCDescriptors[(t + 1) * maxDistance + k] =
              ATSCDescriptors[(t + 1) * maxDistance + k] / (0.5 * gsum);
          MATSDescriptors[(t + 1) * (maxDistance - 1) + k - 1] =
              weights.size() * AATSCDescriptors[(t + 1) * maxDistance + k] /
              weights_centered.array().square().sum();

          // Compute n
          GATSDescriptors[(t + 1) * (maxDistance - 1) + k - 1] =
              ((binaryMatrix.array() * diffSquared.array()).sum() /
               (2.0 * gsum)) /
              gatsdenominator;

        } else {
          ATSCDescriptors[k] = 0.5 * (weights_centered.transpose() *
                                      binaryMatrix * weights_centered)
                                         .value();
          AATSCDescriptors[k] = ATSCDescriptors[k] / (0.5 * gsum);
          MATSDescriptors[k - 1] = weights.size() * AATSCDescriptors[k] /
                                   weights_centered.array().square().sum();
          GATSDescriptors[k - 1] =
              ((binaryMatrix.array() * diffSquared.array()).sum() /
               (2.0 * gsum)) /
              gatsdenominator;
        }
      } else {
        if (t < 11) {
          ATSDescriptors[t * maxDistance + k] = weights.array().square().sum();
          AATSDescriptors[t * maxDistance + k] =
              ATSDescriptors[t * maxDistance + k] / gsum;
          ATSCDescriptors[(t + 1) * maxDistance + k] =
              weights_centered.array().square().sum();
          AATSCDescriptors[(t + 1) * maxDistance + k] =
              ATSCDescriptors[(t + 1) * maxDistance + k] / gsum;
        } else {
          ATSCDescriptors[k] = weights_centered.array().square().sum();
          AATSCDescriptors[k] = ATSCDescriptors[k] / gsum;
        }
      }
    }
  }

  std::vector<double> autocorrDescriptors;
  autocorrDescriptors.reserve(606);  // Pre-allocate memory for 606 elements

  // Append ATSDescriptors
  for (const auto &val : ATSDescriptors) {
    autocorrDescriptors.push_back(val);
  }

  // Append AATSDescriptors
  for (const auto &val : AATSDescriptors) {
    autocorrDescriptors.push_back(val);
  }

  // Append ATSCDescriptors
  for (const auto &val : ATSCDescriptors) {
    autocorrDescriptors.push_back(val);
  }

  // Append AATSCDescriptors
  for (const auto &val : AATSCDescriptors) {
    autocorrDescriptors.push_back(val);
  }

  // Append MATSDescriptors
  for (const auto &val : MATSDescriptors) {
    autocorrDescriptors.push_back(val);
  }

  // Append GATSDescriptors
  for (const auto &val : GATSDescriptors) {
    autocorrDescriptors.push_back(val);
  }

  return autocorrDescriptors;
}

// Function to compute the ATS descriptors without Eigen
std::vector<double> calcAutoCorrelation(const ROMol &mol) {
  // v2.0: Check Gasteiger parameters first (on original mol, before adding H)
  bool gasteiger_ok = checkGasteigerParameters(mol);

  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));
  double *dist = MolOps::getDistanceMat(*hmol, false);  // Topological matrix
  const unsigned int numAtoms = hmol->getNumAtoms();
  const unsigned int numProperties = 12;
  // Lookup tables
  std::map<int, double> vdwmap = VdWAtomicMap();
  std::map<int, double> semap = SandersonENAtomicMap();
  std::map<int, double> pemap = PaulingENAtomicMap();
  std::map<int, double> aremap = Allred_rocow_ENAtomicMap();
  std::map<int, double> pmap = Polarizability94AtomicMap();
  std::map<int, double> imap = ionizationEnergyAtomicMap();
  const auto *tbl = PeriodicTable::getTable();
  // Property matrix as a vector of vectors
  std::vector<std::vector<double>> propertyMatrix(
      numProperties, std::vector<double>(numAtoms, 0.0));  // correct

  // v2.0: Compute Gasteiger charges only if parameters exist
  std::vector<double> gasteigerCharges(numAtoms, 0.0);
  if (gasteiger_ok) {
    try {
      computeGasteigerCharges(*hmol, 12, true);
      for (unsigned int i = 0; i < numAtoms; ++i) {
        const auto *atom = hmol->getAtomWithIdx(i);
        double charge =
            atom->getProp<double>(common_properties::_GasteigerCharge);
        if (atom->hasProp(common_properties::_GasteigerHCharge)) {
          charge += atom->getProp<double>(common_properties::_GasteigerHCharge);
        }
        gasteigerCharges[i] = charge;
      }
    } catch (...) {
      // Gasteiger failed - charges stay at 0.0
    }
  }

  for (unsigned int i = 0; i < numAtoms; ++i) {
    const auto *atom = hmol->getAtomWithIdx(i);
    int atomNumber = atom->getAtomicNum();
    propertyMatrix[0][i] = gasteigerCharges[i];  // not for ATS & AATS output
    propertyMatrix[1][i] = getValenceElectrons(*atom);
    propertyMatrix[2][i] = getSigmaElectrons(*atom);
    propertyMatrix[3][i] = getIntrinsicState(*atom);
    propertyMatrix[4][i] = static_cast<double>(atomNumber);
    propertyMatrix[5][i] = tbl->getAtomicWeight(atomNumber);
    propertyMatrix[6][i] = vdw_volume(vdwmap[atomNumber]);
    propertyMatrix[7][i] = semap[atomNumber];
    propertyMatrix[8][i] = pemap[atomNumber];
    propertyMatrix[9][i] = aremap[atomNumber];
    propertyMatrix[10][i] = pmap[atomNumber];
    propertyMatrix[11][i] = imap[atomNumber];
  }

  // Initialize the topological symetric distance matrix without diagonal
  std::vector<std::vector<double>> distanceMatrix(
      numAtoms, std::vector<double>(numAtoms, 0.0));
  for (unsigned int i = 0; i < numAtoms; ++i) {
    for (unsigned int j = i + 1; j < numAtoms; ++j) {
      distanceMatrix[i][j] = dist[i * numAtoms + j];
    }
  }

  // Compute the ATS descriptors
  const int maxDistance = 8;

  std::vector<std::vector<double>> ATS_(
      maxDistance + 1, std::vector<double>(numProperties - 1, 0.0));
  std::vector<std::vector<double>> AATS_(
      maxDistance + 1, std::vector<double>(numProperties - 1, 0.0));
  std::vector<std::vector<double>> ATSC(
      maxDistance + 1, std::vector<double>(numProperties, 0.0));
  std::vector<std::vector<double>> AATSC(
      maxDistance + 1, std::vector<double>(numProperties, 0.0));
  std::vector<std::vector<double>> MATS(
      maxDistance + 1,
      std::vector<double>(numProperties,
                          0.0));  // we skip the first one at the end!
  std::vector<std::vector<double>> GATS(
      maxDistance + 1,
      std::vector<double>(numProperties,
                          0.0));  // we skip the first one at the end!

  // Centered property values
  std::vector<std::vector<double>> centeredProperties(
      numProperties, std::vector<double>(numAtoms, 0.0));
  for (unsigned int t = 0; t < numProperties; ++t) {
    double sum = 0.0;
    for (unsigned int i = 0; i < numAtoms; ++i) {
      sum += propertyMatrix[t][i];
    }
    double mean = sum / numAtoms;
    for (unsigned int i = 0; i < numAtoms; ++i) {
      centeredProperties[t][i] = propertyMatrix[t][i] - mean;
    }
  }

  // Lag 0: self-correlations
  for (unsigned int t = 0; t < numProperties; ++t) {
    for (unsigned int i = 0; i < numAtoms; ++i) {
      if (t > 0) {
        // we skip Charges
        ATS_[0][t - 1] += propertyMatrix[t][i] * propertyMatrix[t][i];
      }
      ATSC[0][t] += centeredProperties[t][i] * centeredProperties[t][i];
    }
    if (t > 0) {
      // we skip Charges
      AATS_[0][t - 1] = ATS_[0][t - 1] / numAtoms;
    }
    AATSC[0][t] = ATSC[0][t] / numAtoms;
    MATS[0][t] = AATSC[0][t];
    GATS[0][t] = ATSC[0][t] / (numAtoms - 1);
  }

  // Lags 1 to maxLag: pairwise correlations
  for (int k = 1; k <= maxDistance; ++k) {
    int maxkVertexPairs = 0;
    for (unsigned int i = 0; i < numAtoms; ++i) {
      for (unsigned int j = i + 1; j < numAtoms; ++j) {
        if (distanceMatrix[i][j] == k) {
          ++maxkVertexPairs;
          for (unsigned int t = 0; t < numProperties; ++t) {
            double diff = propertyMatrix[t][i] - propertyMatrix[t][j];
            if (t > 0) {
              ATS_[k][t - 1] += propertyMatrix[t][i] * propertyMatrix[t][j];
            }
            ATSC[k][t] += centeredProperties[t][i] * centeredProperties[t][j];
            GATS[k][t] += diff * diff;
          }
        }
      }
    }

    if (maxkVertexPairs > 0) {
      for (unsigned int t = 0; t < numProperties; ++t) {
        if (t > 0) {
          AATS_[k][t - 1] = ATS_[k][t - 1] / maxkVertexPairs;
        }
        AATSC[k][t] = ATSC[k][t] / maxkVertexPairs;
        if (MATS[0][t] > 0.0) {
          MATS[k][t] = AATSC[k][t] / MATS[0][t];
        }
        if (GATS[0][t] > 0.0) {
          GATS[k][t] /= (2 * maxkVertexPairs * GATS[0][t]);
        }
      }
    }
  }
  std::vector<double> descriptors;

  // Flatten the descriptors into a single vector
  for (const auto &vec : {ATS_, AATS_}) {
    for (unsigned int t = 0; t < numProperties - 1; ++t) {
      for (unsigned int k = 0; k <= maxDistance; ++k) {
        descriptors.push_back(vec[k][t]);
      }
    }
  }

  // Flatten the descriptors into a single vector
  for (const auto &vec : {ATSC, AATSC}) {
    for (unsigned int t = 0; t < numProperties; ++t) {
      for (unsigned int k = 0; k <= maxDistance; ++k) {
        descriptors.push_back(vec[k][t]);
      }
    }
  }

  // Flatten the descriptors into a single vector
  for (const auto &vec : {MATS, GATS}) {
    for (unsigned int t = 0; t < numProperties; ++t) {
      for (unsigned int k = 1; k <= maxDistance; ++k) {
        descriptors.push_back(vec[k][t]);
      }
    }
  }

  return descriptors;
}

std::unordered_set<int> findLinkersWithBFS(
    const ROMol &mol, const std::unordered_set<int> &ringAtoms) {
  const RingInfo *ringInfo = mol.getRingInfo();
  std::unordered_set<int> nonRingAtoms;
  std::unordered_set<int> linkers;

  // Collect all non-ring atoms
  int numAtoms = mol.getNumAtoms();
  for (int i = 0; i < numAtoms; ++i) {
    if (ringAtoms.find(i) == ringAtoms.end()) {
      nonRingAtoms.insert(i);
    }
  }

  // Perform BFS for each pair of fused rings
  std::vector<std::vector<int>> atomRings = ringInfo->atomRings();
  for (unsigned int i = 0; i < atomRings.size(); ++i) {
    for (unsigned int j = i + 1; j < atomRings.size(); ++j) {
      // Identify start atoms (non-ring neighbors of the first ring)
      std::queue<std::pair<int, std::vector<int>>>
          q;  // Queue of (current atom, path)
      std::unordered_set<int> visited;

      for (int ringAtom : atomRings[i]) {
        const Atom *atom = mol.getAtomWithIdx(ringAtom);
        for (const Atom *neighbor : mol.atomNeighbors(atom)) {
          int neighborIdx = neighbor->getIdx();
          if (nonRingAtoms.find(neighborIdx) != nonRingAtoms.end()) {
            q.push({neighborIdx, {neighborIdx}});
            visited.insert(neighborIdx);
          }
        }
      }

      // BFS traversal
      while (!q.empty()) {
        auto [current, path] = q.front();
        q.pop();

        const Atom *currentAtom = mol.getAtomWithIdx(current);
        for (const Atom *neighbor : mol.atomNeighbors(currentAtom)) {
          int neighborIdx = neighbor->getIdx();

          if (visited.find(neighborIdx) != visited.end()) continue;

          // Check if we've reached the second ring
          if (ringAtoms.find(neighborIdx) != ringAtoms.end() &&
              std::find(atomRings[j].begin(), atomRings[j].end(),
                        neighborIdx) != atomRings[j].end()) {
            // Valid path found
            linkers.insert(path.begin(), path.end());
            break;
          }

          // Continue BFS through non-ring atoms
          if (nonRingAtoms.find(neighborIdx) != nonRingAtoms.end()) {
            std::vector<int> newPath = path;
            newPath.push_back(neighborIdx);
            q.push({neighborIdx, newPath});
            visited.insert(neighborIdx);
          }
        }
      }
    }
  }

  return linkers;
}

// Function to calculate the FMF ratio
double Framework(const ROMol &mol) {
  const RingInfo *ringInfo = mol.getRingInfo();
  std::unordered_set<int> ringAtoms;

  // Collect all atoms that are part of rings
  std::vector<std::vector<int>> atomRings = ringInfo->atomRings();
  for (const auto &ring : atomRings) {
    for (int atom : ring) {
      ringAtoms.insert(atom);
    }
  }

  // Find linkers using BFS on non-ring atoms
  std::unordered_set<int> linkers = findLinkersWithBFS(mol, ringAtoms);

  // Total number of atoms (including hydrogens)
  std::unique_ptr<ROMol> hmol(MolOps::addHs(mol));
  int totalAtoms = hmol->getNumAtoms();

  // Number of framework atoms: linkers + ring atoms
  int frameworkAtoms = linkers.size() + ringAtoms.size();

  // Calculate FMF
  double FMF = static_cast<double>(frameworkAtoms) / totalAtoms;

  return FMF;
}

std::vector<double> calcFramework(const ROMol &mol) {
  std::vector<double> res(1, 0.);
  res[0] = Framework(mol);
  return res;
}

// BRStates: Tetko version only organis !

const std::vector<std::string> &getOrganicBondKeys() {
  static const std::vector<std::string> organicbondkeys = {
      "7-S-7",   "9-S-7",   "11-S-7",  "13-S-7",  "15-S-7",  "16-S-7",
      "17-S-7",  "19-S-7",  "20-S-7",  "21-S-7",  "22-S-7",  "24-S-7",
      "27-S-7",  "28-S-7",  "30-S-7",  "31-S-7",  "33-S-7",  "34-S-7",
      "36-S-7",  "38-S-7",  "48-S-7",  "50-S-7",  "52-S-7",  "53-S-7",
      "54-S-7",  "70-S-7",  "75-S-7",  "39-S-7",  "8-D-8",   "11-D-8",
      "14-D-8",  "16-D-8",  "23-D-8",  "28-D-8",  "35-D-8",  "49-D-8",
      "52-D-8",  "9-S-9",   "11-S-9",  "13-S-9",  "15-S-9",  "16-S-9",
      "17-S-9",  "19-S-9",  "20-S-9",  "21-S-9",  "22-S-9",  "24-S-9",
      "27-S-9",  "28-S-9",  "30-S-9",  "31-S-9",  "33-S-9",  "34-S-9",
      "36-S-9",  "38-S-9",  "48-S-9",  "50-S-9",  "52-S-9",  "53-S-9",
      "54-S-9",  "70-S-9",  "75-S-9",  "39-S-9",  "10-T-10", "15-T-10",
      "26-T-10", "11-D-11", "11-S-11", "13-S-11", "14-D-11", "15-S-11",
      "16-D-11", "16-S-11", "17-S-11", "19-S-11", "20-S-11", "21-S-11",
      "22-S-11", "23-D-11", "24-S-11", "27-S-11", "28-D-11", "28-S-11",
      "30-S-11", "31-S-11", "33-S-11", "34-S-11", "35-D-11", "36-S-11",
      "38-S-11", "48-S-11", "49-D-11", "50-S-11", "52-D-11", "52-S-11",
      "54-S-11", "70-S-11", "75-S-11", "39-S-11", "12-A-12", "17-A-12",
      "18-A-12", "25-A-12", "29-A-12", "32-A-12", "37-A-12", "13-S-13",
      "15-S-13", "16-S-13", "17-S-13", "19-S-13", "20-S-13", "21-S-13",
      "22-S-13", "24-S-13", "27-S-13", "28-S-13", "30-S-13", "31-S-13",
      "33-S-13", "34-S-13", "36-S-13", "38-S-13", "48-S-13", "50-S-13",
      "52-S-13", "53-S-13", "54-S-13", "70-S-13", "75-S-13", "39-S-13",
      "14-D-14", "16-D-14", "23-D-14", "28-D-14", "35-D-14", "49-D-14",
      "52-D-14", "53-D-14", "15-T-15", "15-S-15", "17-S-15", "19-S-15",
      "20-S-15", "21-S-15", "22-S-15", "24-S-15", "26-T-15", "27-S-15",
      "28-S-15", "30-S-15", "31-S-15", "33-S-15", "34-S-15", "36-S-15",
      "38-S-15", "48-S-15", "50-S-15", "52-S-15", "53-S-15", "54-S-15",
      "70-S-15", "75-S-15", "39-S-15", "16-D-16", "16-S-16", "17-S-16",
      "19-S-16", "20-S-16", "21-S-16", "22-S-16", "23-D-16", "24-S-16",
      "27-S-16", "28-D-16", "28-S-16", "30-S-16", "31-D-16", "34-S-16",
      "35-D-16", "36-S-16", "38-D-16", "48-S-16", "49-D-16", "50-S-16",
      "52-D-16", "52-S-16", "53-D-16", "53-S-16", "54-S-16", "70-S-16",
      "75-S-16", "39-S-16", "17-A-17", "17-S-17", "18-A-17", "19-S-17",
      "20-S-17", "21-S-17", "22-S-17", "24-S-17", "28-S-17", "30-S-17",
      "31-S-17", "33-S-17", "34-S-17", "36-S-17", "37-A-17", "38-S-17",
      "48-S-17", "50-S-17", "51-A-17", "52-S-17", "53-S-17", "54-S-17",
      "70-S-17", "75-S-17", "39-S-17", "18-A-18", "25-A-18", "29-A-18",
      "32-A-18", "37-A-18", "51-A-18", "19-S-19", "20-S-19", "21-S-19",
      "22-S-19", "24-S-19", "27-S-19", "28-S-19", "30-S-19", "31-S-19",
      "33-S-19", "34-S-19", "36-S-19", "38-S-19", "48-S-19", "50-S-19",
      "52-S-19", "53-S-19", "54-S-19", "70-S-19", "75-S-19", "39-S-19",
      "24-S-20", "28-S-20", "30-S-20", "50-S-20", "52-S-20", "53-S-20",
      "21-S-21", "24-S-21", "28-S-21", "52-S-21", "53-S-21", "24-S-22",
      "28-S-22", "30-S-22", "50-S-22", "52-S-22", "53-S-22", "28-D-23",
      "28-S-23", "52-D-23", "52-S-23", "53-D-23", "53-D-23", "24-S-24",
      "28-S-24", "52-S-24", "53-S-24", "25-A-25", "29-S-25", "37-A-25",
      "51-A-25", "28-S-27", "30-S-27", "50-S-27", "52-S-27", "53-S-27",
      "28-D-28", "28-S-28", "30-S-28", "31-S-28", "35-D-28", "36-S-28",
      "38-S-28", "49-D-28", "52-D-28", "52-S-28", "53-D-28", "53-S-28",
      "29-A-29", "37-A-29", "51-A-29", "30-S-30", "31-S-30", "34-S-30",
      "36-S-30", "37-S-30", "52-S-30", "53-S-30", "31-S-31", "35-D-31",
      "35-S-32", "50-S-33", "52-S-33", "53-S-33", "34-S-34", "36-S-34",
      "52-S-34", "53-S-34", "52-D-35", "53-D-35", "39-S-36", "50-S-36",
      "52-S-36", "53-S-36", "54-S-36", "70-S-36", "75-S-36", "37-A-37",
      "51-A-37", "38-S-38", "39-S-38", "50-S-38", "52-S-38", "53-S-38",
      "54-S-38", "70-S-38", "75-S-38", "39-S-39", "48-S-39", "50-S-39",
      "52-S-39", "53-S-39", "54-S-39", "70-S-39", "75-S-39", "48-S-48",
      "50-S-48", "52-S-48", "53-S-48", "52-D-49", "53-D-49", "50-S-50",
      "52-S-50", "53-S-50", "54-S-50", "70-S-50", "75-S-50", "52-S-52",
      "53-S-52", "54-S-52", "70-S-52", "75-S-52", "54-S-54", "70-S-54",
      "75-S-54", "70-S-70", "75-S-70", "75-S-75"};
  return organicbondkeys;
}

struct BondEStateResult {
  std::vector<std::string> SumKeys;
  std::vector<double> BEStotal;
  std::vector<double> SumBES;
  std::vector<double> nBES;
  std::vector<double> minBES;
  std::vector<double> maxBES;
};

// retreive the positional of the
std::vector<int> NamePosES(
    const ROMol &mol,
    const std::vector<std::pair<std::string, std::shared_ptr<RWMol>>>
        &queries) {
  size_t nAtoms = mol.getNumAtoms();
  std::vector<int> pos(nAtoms, 0);  // Initialize positions with 0
  for (unsigned int idx = 0; idx < queries.size(); idx++) {
    auto entry = queries[idx];
    if (!entry.second) continue;  // Skip invalid SMARTS patterns

    std::vector<MatchVectType> qmatches;
    if (SubstructMatch(mol, *entry.second, qmatches, true)) {
      for (unsigned int i = 0; i < qmatches.size(); ++i) {
        int atomIdx = qmatches[i][0].second;
        // std::cout << entry.first << ":" << atomIdx << ": " << idx+1 << "\n";
        pos[atomIdx] = idx + 1;
      }
    }
  }

  /*std::cout << "C++ nPatts: " << queries.size() << std::endl;
  std::cout << "pos : ";
  for (int i=0; i< pos.size();i++) {
      std::cout << pos[i] << " ";
  }
  std::cout << "\n";
  */
  return pos;
}

// Tetko version
BondEStateResult getBEStateFeatures(const ROMol &mol, bool extended) {
  size_t nBonds = mol.getNumBonds();
  size_t nAtoms = mol.getNumAtoms();

  // Precompute atomic EState indices
  std::vector<double> Is = calcIStateIndices(mol);

  /*std::cout << "IStatesIndices : ";
  for (int i=0; i< Is.size();i++) {
      std::cout << Is[i] << " ";
  }
  std::cout << "\n";
  */
  // Get the distance matrix using RDKit's MolOps::getDistanceMat
  double *dists = MolOps::getDistanceMat(
      mol, false, false, false);  // no bond order, no weights, no hydrogens

  /*
  for (int i=0; i<nAtoms; i++){
      for (int j=0; j<nAtoms; j++){
          std::cout << dists[j * nAtoms + i]+1. << ", ";
      }
      std::cout << "\n";
  }
  std::cout << "\n";
  */

  // Bond-specific indices
  std::vector<double> Iij(nBonds, 0.0);
  std::vector<std::string> BEScode(nBonds);

  const auto &queries = extended ? GetesExtQueries() : GetesQueries();
  const std::vector<int> pos = NamePosES(mol, queries);

  // Compute bond contributions
  for (size_t t = 0; t < nBonds; ++t) {
    const Bond *bt = mol.getBondWithIdx(t);
    int i = bt->getBeginAtomIdx();
    int j = bt->getEndAtomIdx();

    // Bond type and code
    std::string btype;
    switch (bt->getBondType()) {
      case Bond::SINGLE:
        btype = "S";
        break;
      case Bond::DOUBLE:
        btype = "D";
        break;
      case Bond::TRIPLE:
        btype = "T";
        break;
      default:
        btype = "A";
        break;
    }

    // Calculate Iij
    Iij[t] = 0.5 * (Is[i] + Is[j]);

    // Generate bond code
    int posi = pos[i];
    int posj = pos[j];
    if (posi >= posj) {
      BEScode[t] =
          std::to_string(posi) + "-" + btype + "-" + std::to_string(posj);
    } else {
      BEScode[t] =
          std::to_string(posj) + "-" + btype + "-" + std::to_string(posi);
    }
  }

  // Initialize BES
  std::vector<double> BES(nBonds, 0.0);

  // Compute edge EState contributions
  for (size_t i = 0; i < nBonds; ++i) {
    const Bond *bt1 = mol.getBondWithIdx(i);
    int i1 = bt1->getBeginAtomIdx();
    int i2 = bt1->getEndAtomIdx();

    for (size_t j = 0; j < i; ++j) {
      const Bond *bt2 = mol.getBondWithIdx(j);
      int j1 = bt2->getBeginAtomIdx();
      int j2 = bt2->getEndAtomIdx();
      double d11 = dists[i1 * nAtoms + j1];
      double d12 = dists[i1 * nAtoms + j2];
      double d21 = dists[i2 * nAtoms + j1];
      double d22 = dists[i2 * nAtoms + j2];

      // Topological distance adding one or ???
      double p = 1 + (d11 + d12 + d21 + d22) / 4;
      double dI = (Iij[i] - Iij[j]) / (p * p);
      BES[i] += dI;
      BES[j] -= dI;
    }
  }

  // Total EState contributions
  std::vector<double> BEStotal(nBonds);
  for (size_t i = 0; i < nBonds; ++i) {
    BEStotal[i] = BES[i] + Iij[i];
  }
  /*
  std::cout << "Iij : ";
  for (int i=0; i< Iij.size();i++) {
      std::cout << Iij[i] << " ";
  }
  std::cout << "\n";


  std::cout << "BES : ";
  for (int i=0; i< BES.size();i++) {
      std::cout << BES[i] << " ";
  }
  std::cout << "\n";


  std::cout << "BEStotal : ";
  for (int i=0; i< BEStotal.size();i++) {
      std::cout << BEStotal[i] << " ";
  }
  std::cout << "\n";
  */

  // Summ contributions by bond code
  std::unordered_map<std::string, std::vector<double>> codeMap;
  for (size_t i = 0; i < nBonds; ++i) {
    codeMap[BEScode[i]].push_back(BEStotal[i]);
  }

  // Aggregated results
  std::vector<std::string> SumKeys;
  std::vector<double> SumBES, minBES, maxBES, nBES;

  for (const auto &[key, values] : codeMap) {
    SumKeys.push_back("S" + key);
    SumBES.push_back(std::accumulate(values.begin(), values.end(), 0.0));
    nBES.push_back(static_cast<double>(values.size()));
    minBES.push_back(*std::min_element(values.begin(), values.end()));
    maxBES.push_back(*std::max_element(values.begin(), values.end()));
  }
  /*
  std::sort(SumKeys.begin(), SumKeys.end());
  for (const auto &k : SumKeys) {
      std::cout << k << " ";
  }
  std::cout << "\n";
  */
  return {SumKeys, BEStotal, SumBES, nBES, minBES, maxBES};
}

// Function to compute Bond E-State fingerprints
std::vector<double> calcBEStateDescs(const ROMol &mol) {
  // Call the function to calculate Bond E-State descriptors using the extended
  // patterns (aka true)
  auto [SumKeys_i, BEStotal_i, SumBES_i, nBES_i, minBES_i, maxBES_i] =
      getBEStateFeatures(mol, true);

  const auto &orgbondkeys = getOrganicBondKeys();

  // Initialize results with size equal to organicbondkeys plus one (for
  // unmatched patterns)
  size_t nKeys = orgbondkeys.size();

  std::vector<double> sumsBES(nKeys + 1, 0.);
  std::vector<double> nBES(nKeys + 1, 0.);
  std::vector<double> minBES(nKeys + 1, 99999);
  std::vector<double> maxBES(nKeys + 1, 0.);

  // Process each descriptor
  for (size_t i = 0; i < SumBES_i.size(); ++i) {
    // std::cout << "Key: " << SumKeys_i[i] << " BEStotal: " << BEStotal_i[i] <<
    // " SumBES:" << SumBES_i[i] << " nBES:" << nBES_i[i] << " minBES:" <<
    // minBES_i[i] << " maxBES:" << maxBES_i[i] << "\n";

    const auto &pattern = SumKeys_i[i];
    const auto &descriptor =
        pattern.substr(1);  // Extract the key after "S" only one to remove!

    // Check if the descriptor exists in organicbondkeys
    auto it = std::find(orgbondkeys.begin(), orgbondkeys.end(), descriptor);
    if (it != orgbondkeys.end()) {
      // Get the position index
      size_t posidx = std::distance(orgbondkeys.begin(), it);

      // std::cout << "found at position: " << posidx <<"\n";

      if (posidx < sumsBES.size()) {
        // Update the corresponding values
        sumsBES[posidx] += SumBES_i[i];
        nBES[posidx] += nBES_i[i];
        minBES[posidx] = std::min(minBES[posidx], minBES_i[i]);
        maxBES[posidx] = std::max(maxBES[posidx], maxBES_i[i]);
      } else {
        BOOST_LOG(rdWarningLog) << "Out-of-bounds access detected for posIdx: " << posidx
                  << " : " << sumsBES.size() << "\n";
      }

    } else {
      // Update the last index (unmatched patterns)
      size_t unmatchedIdx = nKeys;
      sumsBES[unmatchedIdx] += SumBES_i[i];
      nBES[unmatchedIdx] += nBES_i[i];
      minBES[unmatchedIdx] = std::min(minBES[unmatchedIdx], minBES_i[i]);
      maxBES[unmatchedIdx] = std::max(maxBES[unmatchedIdx], maxBES_i[i]);
    }
  }

  for (unsigned int i = 0; i < nKeys + 1; i++) {
    if (minBES[i] == 99999) {
      minBES[i] = 0.;
    }
  }

  // Concatenate all vectors into a single vector<double>
  std::vector<double> concatenatedResult;
  concatenatedResult.reserve(sumsBES.size() + nBES.size() + minBES.size() +
                             maxBES.size());

  concatenatedResult.insert(concatenatedResult.end(), sumsBES.begin(),
                            sumsBES.end());
  concatenatedResult.insert(concatenatedResult.end(), nBES.begin(), nBES.end());
  concatenatedResult.insert(concatenatedResult.end(), minBES.begin(),
                            minBES.end());
  concatenatedResult.insert(concatenatedResult.end(), maxBES.begin(),
                            maxBES.end());

  return concatenatedResult;
}

static const std::vector<std::string> AFragments = {
    "[C][OX2H]",
    "[c][OX2H]",
    "[C][NX3;H2]",
    "[c][NX3;H2;!$(NC=O)]",
    "[C][NX3;H1;!R][C]",
    "[C][NX3;H1;R][C]",
    "[c][NX3;H1;!$(NC=O)][C]",
    "[c][nX3;H1][c]",
    "[CX3](=O)[OX1H0-,OX2H1]",
    "[CX3](=[OX1])[NX3;H2]",
    "[CX3](=[OX1])[NX3;H1][C]",
    "[CX3](=[OX1])[NX3;H1][c]",
    "[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]",
    "[NX3;H1]C(=[OX1])[NX3;H1]",
    "[NX3;H0]C(=[OX1])[NX3;H1]",
    "[NX3;H1]C(=[OX1])O",
    "[NX3;H1]C(=N)[NX3;H0]",
    "[C]#[CH]",
    "P[OH,O-]",
    "[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
    "[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
    "[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]",
    "[nX3;H1]:n",
    "[nX3;H1]:c:n",
    "[OX2;H1]CC[O,N]",
    "[OX2;H1]C[C,N]=[O,S]",
    "[OX2;H1]c1ccccc1[O,NX3]",
    "[OX2;H1]c1ccccc1C=[O,S]",
    "[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "[NH,NH2,NH3+]CC[O,N]",
    "[NH,NH2,NH3+]c1ccccc1[O,N]",
    "[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]",
    "[OX2H]c1ccccc1[Cl,Br,I]",
    "[OX1]=[C,c]~[C,c]C[OH]",
    "[OH]c1cccc2cccnc12",
    "[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
    "[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
    "[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
    "[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
    "[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
    "[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
    "[OH]c1c([CX4])cccc1[CX4]",
    "[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]",
    "[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1",
    "[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1",
    "[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1",
    "[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1",
    "[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1"};

static const std::vector<std::string> BSELFragments = {
    "[CX4H3]",
    "[CX4H2]",
    "[CX4H1]",
    "[CX4H0]",
    "*=[CX3H2]",
    "[$(*=[CX3H1]),$([cX3H1](a)a)]",
    "[$(*=[CX3H0]),$([cX3H0](a)(a)A)]",
    "c(a)(a)a",
    "*#C",
    "[C][NX3;H2]",
    "[c][NX3;H2]",
    "[C][NX3;H1][C]",
    "[c][NX3;H1]",
    "[c][nX3;H1][c]",
    "[C][NX3;H0](C)[C]",
    "[c][NX3;H0](C)[C]",
    "[c][nX3;H0][c]",
    "*=[Nv3;!R]",
    "*=[Nv3;R]",
    "[nX2H0,nX3H1+](a)a",
    "N#C[A;!#1]",
    "N#C[a;!#1]",
    "[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]",
    "[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]",
    "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]",
    "[OH]",
    "[OX2;H0;!R]",
    "[OX2;H0;R]",
    "[oX2](a)a",
    "*=O",
    "[SX2](*)*",
    "[sX2](a)a",
    "*=[SX1]",
    "*=[SX3]",
    "[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]",
    "[S,s]",
    "[P,p]",
    "FA",
    "Fa",
    "Cl",
    "Br",
    "I",
    "[CX3;!R](=[OX1])[OX2H0]",
    "[CX3;R](=[OX1])[OX2H0;R]",
    "P(=[OX1])(O)(O)O",
    "[CX3](=[OX1])([OX2H0])[OX2H0]",
    "[CX3](=O)[OX1H0-,OX2H1]",
    "nC=[OX1]",
    "[N;!R]C=[OX1]",
    "[N;R][C;R]=[OX1]",
    "[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]",
    "NC(=[OX1])N",
    "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]",
    "[CX3](=[OX1])[NX3][CX3](=[OX1])",
    "C1(=[OX1])C=CC(=[OX1])C=C1",
    "[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
    "*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1",
    "[OX2H]CC[O,N]",
    "[OX2H]C[C,N]=[O,S]",
    "[OX2H]c1ccccc1[O,Nv3]",
    "[OX2H]c1ccccc1C=[O,S]",
    "[OX2H]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "[NH,NH2,NH3+]CC[O,N]",
    "[NH,NH2,NH3+]c1ccccc1[O,N]",
    "[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]",
    "[OX2H]c1ccccc1[Cl,Br,I]",
    "[CX4]([OH])[CX4][OH]",
    "n:n",
    "o:n",
    "n:c:n",
    "o:c:n",
    "n:c:c:n",
    "[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]",
    "[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]",
    "[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]",
    "P(=[OX1])N",
    "Nc:n",
    "[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]",
    "[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
    "[OX2]-c:c-[OX2]"};

// Precompile SMARTS patterns for efficiency
static const std::vector<std::shared_ptr<RWMol>> &GetQueriesA() {
  static const std::vector<std::shared_ptr<RWMol>> queriesA = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smi : AFragments) {
      auto mol = SmartsToMol(smi);
      if (mol) {
        res.emplace_back(std::move(mol));
      } else {
        BOOST_LOG(rdWarningLog) << "Invalid SMARTS: " << smi << std::endl;
      }
    }
    return res;
  }();
  return queriesA;
}

static const std::vector<std::shared_ptr<RWMol>> &GetQueriesB() {
  static const std::vector<std::shared_ptr<RWMol>> queriesB = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smi : BSELFragments) {
      auto mol = SmartsToMol(smi);
      if (mol) {
        res.emplace_back(std::move(mol));
      } else {
        BOOST_LOG(rdWarningLog) << "Invalid SMARTS: " << smi << std::endl;
      }
    }
    return res;
  }();
  return queriesB;
}

static const std::vector<double> coefAFragments = {
    0.345,  0.543,  0.177,  0.247,  0.087,  0.321,  0.194,  0.371,  0.243,
    0.275,  0.281,  -0.091, 0.356,  -0.165, -0.119, -0.105, 0.170,  0.082,
    0.493,  0.019,  0.050,  -0.362, 0.118,  0.1,    0.051,  0.194,  0.042,
    -0.089, -0.161, -0.251, -0.418, -0.45,  -0.155, 0.,     -0.093, -0.11,
    -0.601, -0.475, 0.119,  0.176,  0.08,   0.084,  0.085,  0.055,  -0.162,
    -0.181, 0.195,  -0.203, 0.096,  0.185,  0.203,  0.003};

static const std::vector<double> coefBHFragments = {
    0.007,  0.,     0.011,  0.037,  0.019,  0.011,  0.,     0.019,  0.028,
    0.481,  0.275,  0.541,  0.415,  0.316,  0.653,  0.321,  0.392,  0.200,
    0.596,  0.321,  0.242,  0.103,  -0.476, -0.525, -0.204, 0.307,  0.211,
    0.331,  0.047,  0.334,  0.168,  0.043,  0.071,  0.448,  -0.188, 0.,
    1.183,  -0.036, 0.,     0.,     -0.011, 0.,     -0.206, -0.214, -0.394,
    -0.267, -0.308, -0.095, -0.287, -0.231, -0.446, -0.076, -0.252, -0.148,
    -0.051, -0.014, 0.013,  0.267,  0.,     -0.068, -0.079, -0.387, -0.126,
    0.,     -0.059, -0.045, -0.130, 0.,     -0.132, -0.157, -0.098, -0.170,
    -0.089, 0.031,  -0.035, -0.023, -0.668, -0.042, 0.131,  -0.408, -0.216,
    0.071};

static const std::vector<double> coefBOFragments = {
    0.000,  0.000,  0.020,  0.047,  0.024,  0.012,  0.000,  0.018,  0.032,
    0.486,  0.326,  0.543,  0.426,  0.267,  0.655,  0.338,  0.338,  0.202,
    0.589,  0.300,  0.245,  0.093,  -0.595, -0.533, -0.202, 0.311,  0.226,
    0.330,  0.060,  0.339,  0.175,  0.083,  0.069,  0.319,  -0.190, 0.000,
    1.189,  -0.033, 0.000,  0.000,  0.000,  0.000,  -0.223, -0.169, -0.408,
    -0.298, -0.312, -0.038, -0.292, -0.242, -0.443, -0.054, -0.251, -0.149,
    -0.050, -0.016, 0.010,  0.218,  0.000,  -0.090, -0.122, -0.403, -0.120,
    0.000,  -0.027, -0.069, -0.130, -0.018, -0.094, -0.141, -0.113, -0.184,
    -0.073, 0.025,  -0.033, -0.025, -0.668, -0.057, 0.129,  -0.405, -0.218,
    0.064};

static const std::vector<double> coefSFragments = {
    -0.075, 0.,     0.036,  0.071,  -0.085, 0.050,  0.101,  0.121,  0.034,
    0.175,  0.383,  0.265,  0.311,  0.221,  0.323,  0.295,  0.265,  0.125,
    0.254,  0.223,  0.694,  0.390,  0.,     -0.231, -0.476, 0.247,  0.185,
    0.185,  0.,     0.370,  0.189,  0.,     0.618,  1.065,  -0.505, 0.643,
    0.703,  -0.042, 0.,     0.082,  0.161,  0.198,  -0.225, 0.360,  -0.240,
    -0.190, -0.412, -0.076, 0.175,  -0.1,   -0.569, -0.553, -0.588, -0.510,
    -0.411, -0.050, 0.000,  1.029,  -0.067, -0.095, -0.237, -0.344, -0.276,
    -0.102, 0.,     -0.140, -0.120, 0.052,  0.024,  0.047,  -0.040, 0.087,
    -0.051, -0.043, -0.038, 0.,     -0.452, 0.098,  0.,     0.434,  0.380,
    0.277};

static const std::vector<double> coefEFragments = {
    -0.104, 0.,     0.089,  0.187, -0.045, 0.068,  0.18,   0.3,   0.04,
    0.085,  0.163,  0.138,  0.192, -0.03,  0.22,   0.346,  0.083, 0.117,
    0.121,  0.046,  0.,     0.,    0.2,    0.21,   0.,     0.061, 0.014,
    0.013,  -0.125, -0.041, 0.33,  0.116,  0.364,  0.413,  0.,    0.465,
    0.295,  -0.18,  -0.23,  0.023, 0.196,  0.533,  -0.113, 0.,    -0.1,
    0.,     -0.192, 0.221,  0.,    0.061,  -0.111, -0.11,  0.,    0.,
    0.,     -0.017, 0.012,  0.285, 0.029,  0.,     -0.069, 0.,    0.,
    0.,     0.,     0.,     -0.1,  -0.043, 0.092,  -0.113, 0.,    0.052,
    0.,     0.,     0.,     0.,    -0.08,  0.185,  0.,     0.,    0.,
    0.248};

static const std::vector<double> coefLFragments = {
    0.321,  0.499,  0.449,  0.443, 0.244,  0.469,  0.624,  0.744,  0.332,
    0.781,  0.949,  0.568,  0.912, 1.25,   0.4,    0.869,  0.794,  -0.235,
    -0.24,  0.574,  0.757,  0.732, 0.278,  0.347,  0.,     0.672,  0.36,
    0.359,  0.057,  0.495,  1.258, 0.848,  0.954,  2.196,  0.,     0.554,
    2.051,  -0.143, -0.147, 0.669, 1.097,  1.590,  -0.39,  0.406,  -0.483,
    0.,     -0.369, 0.,     0.603, 0.583,  0.,     0.,     0.,     0.,
    0.,     -0.111, 0.054,  0.488, -0.072, -0.337, 0.,     -0.303, -0.364,
    0.062,  0.,     0.169,  -0.4,  0.1,    -0.179, 0.,     0.042,  0.209,
    -0.058, -0.081, -0.026, 0.,    0.,     0.149,  -0.145, 0.,     0.,
    0.13};

std::vector<double> calcAbrahams(const ROMol &mol) {
  std::vector<double> retval(6, 0.0);

  try {
    // Calculate A descriptor
    auto queriesA = GetQueriesA();
    for (size_t i = 0; i < queriesA.size(); ++i) {
      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *queriesA[i], matches, true);  // uniquify = true
      retval[0] += matches.size() * coefAFragments[i];
    }

    // Calculate BSEL descriptors
    int sulphurCount = 0;
    auto queriesB = GetQueriesB();
    for (size_t i = 0; i < queriesB.size(); ++i) {
      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *queriesB[i], matches, true);  // uniquify = true

      int uniqueMatches = matches.size();
      if (30 <= i && i <= 34) {
        sulphurCount += uniqueMatches;
      } else if (i == 35) {
        uniqueMatches -= sulphurCount;
      }

      retval[1] += uniqueMatches * coefBHFragments[i];  // BH
      retval[2] += uniqueMatches * coefBOFragments[i];  // BO
      retval[3] += uniqueMatches * coefSFragments[i];   // S
      retval[4] += uniqueMatches * coefEFragments[i];   // E
      retval[5] += uniqueMatches * coefLFragments[i];   // L
    }

    // Add intercepts
    if (coefAFragments.size() > queriesA.size()) {
      retval[0] += coefAFragments.back();  // A intercept
    }
    if (coefBHFragments.size() > queriesB.size()) {
      retval[1] += coefBHFragments.back();  // BH intercept
      retval[2] += coefBOFragments.back();  // BO intercept
      retval[3] += coefSFragments.back();   // S intercept
      retval[4] += coefEFragments.back();   // E intercept
      retval[5] += coefLFragments.back();   // L intercept
    }

  } catch (const std::exception &e) {
    BOOST_LOG(rdWarningLog) << "Error in SMARTS matching: " << e.what() << std::endl;
    throw std::runtime_error("Error in SMARTSQueryTool");
  }

  return retval;
}

//// IC based on the Roy, Basak, Harriss, Magnuson paper Neighorhoo complexities
///and symmetry of chemical graphs and their biological applications
// in this algorithm we use the cluster list to iterate per radius not the whole
// atoms accross clusters. so we don't need to compare all the keys over all
// atoms but only in a cluster the logic is to split the cluster until we get a
// singleton cluster or we rich radius 5 We don't need to generate a special key
// for a radius because we compute all the radius in one iterative step. this
// make the code simpler and sorting and comparison more intuitive.

// Function to generate a key for an atom's environment based on the paper we
// don't need the Neighbor degree for a delta radius key cluster separation. but
// we may need the "truncated leaf" how to encode them ? A truncated leaf
// dictionnary ? for me the logic would be: truncated leaf key*100 + radius of
// the truncated leaf to be included in the remaining extension radius...
// TODO : check if "-2" case is properly used in cluster because by definition
// an empty key is a key too to discrimitate!

int generateKey(int rootNum, int rootDeg, int bondOrder, int neighNum) {
  // return (rootNum * 10 + rootDeg) * 1000 + bondOrder * 100 + neighNum;
  return (rootNum * 10 + rootDeg) * 1000 + bondOrder * 100 + neighNum;
}

int getbondtypeint(const Bond::BondType &bd) {
  if (bd == Bond::BondType::AROMATIC) {
    return 4;
  } else {
    return static_cast<int>(bd);
  }
}

// Function to get bond type as double
double getbondtypeindouble(const Bond::BondType &bd) {
  switch (bd) {
    case Bond::BondType::SINGLE:
      return 1.0;
    case Bond::BondType::DOUBLE:
      return 2.0;
    case Bond::BondType::TRIPLE:
      return 3.0;
    case Bond::BondType::QUADRUPLE:
      return 4.0;
    case Bond::BondType::AROMATIC:
      return 1.5;
    default:
      return static_cast<double>(bd);
  }
}

// Initialize adjacency and shortest-path matrices
std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
initializeMatrixAndSP(int nAtoms, int maxradius) {
  std::vector<std::vector<int>> M(
      nAtoms, std::vector<int>(nAtoms, -1));  // not sure of N+1 here ???
  std::vector<std::vector<int>> SP(nAtoms, std::vector<int>(maxradius + 1, -1));

  for (int i = 0; i < nAtoms; ++i) {
    M[i][0] = i;
    SP[i][0] = 0;
  }

  return {M, SP};
}

// Find the last occupied position in a matrix row
int findLastOccupied(const std::vector<std::vector<int>> &M, int atomIdx) {
  auto it = std::find(M[atomIdx].begin(), M[atomIdx].end(), -1);
  return it == M[atomIdx].begin() ? 0
                                  : std::distance(M[atomIdx].begin(), it) - 1;
}

// Update clusters by grouping atoms with identical keys
std::vector<std::vector<int>> updateClustersWithKeys(
    const std::vector<std::pair<int, std::vector<int>>> &keys) {
  std::vector<std::pair<int, std::vector<int>>> sortedKeys = keys;
  std::sort(sortedKeys.begin(), sortedKeys.end(),
            [](auto &a, auto &b) { return a.second < b.second; });

  std::vector<std::vector<int>> newClusters;
  std::vector<int> currentCluster;
  std::vector<int> currentKey;

  for (auto &[atomIdx, key] : sortedKeys) {
    if (currentCluster.empty() || key != currentKey) {
      if (!currentCluster.empty()) {
        newClusters.push_back(currentCluster);
      }
      currentCluster = {atomIdx};
      currentKey = key;
    } else {
      currentCluster.push_back(atomIdx);
    }
  }
  if (!currentCluster.empty()) newClusters.push_back(currentCluster);
  return newClusters;
}

// Main pipeline
std::map<int, std::vector<std::vector<int>>> computePipeline(
    RWMol &mol, int maxRadius, bool addDeadKeys = false) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  if (nAtoms == 0) {
    BOOST_LOG(rdWarningLog) << "Error: Molecule has no atoms." << std::endl;
    return {};
  }

  std::string smi = MolToSmiles(mol);

  bool debug = false;  // (smi=="FP(F)F" || smi=="BrBr");

  if (debug) {
    BOOST_LOG(rdWarningLog) << "Debugging enabled for molecule: " << smi
              << "n & NumAtoms: " << nAtoms << std::endl;
  }

  auto [M, SP] = initializeMatrixAndSP(nAtoms, maxRadius);

  if (debug) {
    BOOST_LOG(rdDebugLog) << "Initial M matrix:\n";
    for (const auto &row : M) {
      for (int val : row) {
        BOOST_LOG(rdWarningLog) << val << " ";
      }
      BOOST_LOG(rdWarningLog) << "\n";
    }

    BOOST_LOG(rdDebugLog) << "Initial SP matrix:\n";
    for (const auto &row : SP) {
      for (int val : row) {
        BOOST_LOG(rdDebugLog) << val << " ";
      }
      BOOST_LOG(rdDebugLog) << "\n";
    }
  }

  std::map<int, std::vector<std::vector<int>>> CN;  // Combined CN and AN
  std::map<int, std::vector<int>> clustersByAN;

  // Radius 0: Group atoms by atomic number
  for (auto &atom : mol.atoms()) {
    clustersByAN[atom->getAtomicNum()].push_back(atom->getIdx());
  }

  std::vector<std::vector<int>> clusters;
  for (auto &[_, atoms] : clustersByAN) {
    clusters.push_back(atoms);
  }

  // Initialize CN[0]
  CN[0].resize(2);  // Two vectors: sizes and last atomic values
  for (auto &cluster : clusters) {
    if (cluster.empty()) continue;  // Skip empty clusters
    CN[0][0].push_back(static_cast<int>(cluster.size()));  // Cluster size
    int atomIdx = cluster.back();
    if (atomIdx < 0 || atomIdx >= nAtoms) {
      BOOST_LOG(rdWarningLog) << "Error: Invalid atom index: " << atomIdx << std::endl;
      continue;
    }
    CN[0][1].push_back(
        mol.getAtomWithIdx(atomIdx)->getAtomicNum());  // Last atomic number
  }

  if (debug) {
    BOOST_LOG(rdDebugLog) << "\nM matrix before radius 1:\n";
    for (const auto &row : M) {
      for (int val : row) {
        BOOST_LOG(rdDebugLog) << val << " ";
      }
      BOOST_LOG(rdDebugLog) << "\n";
    }

    BOOST_LOG(rdDebugLog) << "\nSP matrix before radius 1:\n";
    for (const auto &row : SP) {
      for (int val : row) {
        BOOST_LOG(rdDebugLog) << val << " ";
      }
      BOOST_LOG(rdDebugLog) << "\n";
    }
  }

  for (int r = 1; r <= maxRadius; ++r) {
    bool stopExpansion = true;

    std::vector<std::vector<int>> newClusters;

    for (auto &cluster : clusters) {
      if (cluster.size() == 1) {
        newClusters.push_back(cluster);
        continue;
      }

      std::vector<std::pair<int, std::vector<int>>> clusterKeys;

      for (int atomIdx : cluster) {
        int start =
            SP[atomIdx]
              [r - 1];  // get the start neighbor to visite at this radius
        int stop = findLastOccupied(M, atomIdx);  // get the end

        if (start == -2) {
          continue;  // Skip further processing but allow iteration
        }

        std::vector<int> eqKeys;
        std::vector<int> neighbors;

        if (debug) {
          BOOST_LOG(rdDebugLog) << "Radius " << r << ", Atom " << atomIdx
                    << " .Symbol: " << mol.getAtomWithIdx(atomIdx)->getSymbol()
                    << ", Start " << start << ", Stop " << stop << std::endl;
        }

        for (int pos = start; pos <= stop; ++pos) {
          int rootIdx = M[atomIdx][pos];
          if (rootIdx < 0 || rootIdx >= nAtoms) {
            std::cout << "Atom index out of boundaries:" << rootIdx
                      << "smile:" << smi << "\n";
            continue;
          }
          const Atom *rootAtom = mol.getAtomWithIdx(rootIdx);
          int rootNum = rootAtom->getAtomicNum();
          int rootDeg = rootAtom->getDegree();

          for (const auto &nb : mol.atomNeighbors(rootAtom)) {
            int nbIdx = nb->getIdx();
            if (std::find(M[atomIdx].begin(), M[atomIdx].begin() + stop + 1,
                          nbIdx) != M[atomIdx].begin() + stop + 1) {
              continue;  // Already visited no effect of adding truncated branch
                         // this was already seen in previous radius!
            }

            neighbors.push_back(nbIdx);
            const Bond *bond = mol.getBondBetweenAtoms(rootIdx, nbIdx);
            int bondOrder = getbondtypeint(
                bond->getBondType());  // don't need kekulize like in Mordred
            int neighNum = mol.getAtomWithIdx(nbIdx)->getAtomicNum();
            // the logic is to look at the dgree of the source atom not the
            // destination as we are at a delta order comparison (relation not
            // absolute detection)
            eqKeys.push_back(
                generateKey(rootNum, rootDeg, bondOrder, neighNum));
          }
        }

        // option one for the deadkeys
        if (neighbors.empty()) {
          if (SP[atomIdx][r] == -1) {
            SP[atomIdx][r] = -2;  // Mark as exhausted
          }
          // Add dead key if enabled
          if (addDeadKeys) {
            eqKeys.push_back(-2);
          }
        }

        std::sort(eqKeys.begin(), eqKeys.end());
        clusterKeys.emplace_back(atomIdx, eqKeys);

        // Appends only validated visited new neighbors, not all neighbors
        // (i.e., exclude backward, leaf, or cycle paths)
        if (!neighbors.empty()) {
          stopExpansion =
              false;  // Continue expansion if new neighbors are found
          for (int i = 0; i < static_cast<int>(neighbors.size()); ++i) {
            int freeSlot = findLastOccupied(M, atomIdx) + 1;
            if (i == 0) SP[atomIdx][r] = (freeSlot < nAtoms) ? freeSlot : -2;
            if (freeSlot < nAtoms) {
              M[atomIdx][freeSlot] = neighbors[i];
            }
          }
        }
        // option two for the deadkeys
        // else {
        //    if (SP[atomIdx][r] == -1) {
        //        SP[atomIdx][r] = -2;  // Mark as exhausted, no further
        //        expansion
        //    }
        //}
      }

      auto subClusters = updateClustersWithKeys(clusterKeys);
      newClusters.insert(newClusters.end(), subClusters.begin(),
                         subClusters.end());
    }

    clusters = newClusters;

    if (debug) {
      BOOST_LOG(rdDebugLog) << "M Matrix after radius " << r << ":\n";
      for (const auto &row : M) {
        for (int val : row) {
          BOOST_LOG(rdDebugLog) << val << " ";
        }
        BOOST_LOG(rdDebugLog) << std::endl;
      }

      BOOST_LOG(rdDebugLog) << "SP Matrix after radius " << r << ":\n";
      for (const auto &row : SP) {
        for (int val : row) {
          BOOST_LOG(rdDebugLog) << val << " ";
        }
        BOOST_LOG(rdDebugLog) << std::endl;
      }
    }

    CN[r].resize(2);  // Two vectors: sizes and last atomic values
    for (auto &cluster : clusters) {
      if (cluster.empty()) continue;  // Skip empty clusters
      int atomIdx = cluster.back();   // Atom index

      if (atomIdx < 0 || atomIdx >= nAtoms) {
        std::cout << "Atom index out of boundaries:" << atomIdx << "\n";
        continue;
      }

      int atomicNum =
          mol.getAtomWithIdx(atomIdx)->getAtomicNum();  // Last atomic mass

      CN[r][0].push_back(
          static_cast<int>(cluster.size()));  // store Cluster size
      CN[r][1].push_back(atomicNum);          // store atomic mass
    }

    if (stopExpansion || std::all_of(clusters.begin(), clusters.end(),
                                     [](auto &c) { return c.size() == 1; })) {
      if (debug) {
        BOOST_LOG(rdDebugLog) << "Stopping expansion at radius " << r << " - Reason: "
                  << (stopExpansion ? "No new neighbors"
                                    : "All clusters are singletons")
                  << std::endl;
      }

      break;
    }
  }

  // Finalize CN by padding remaining radii

  if (!CN.empty()) {
    auto lastRadius = CN.rbegin()->first;
    for (int rr = lastRadius + 1; rr <= maxRadius; ++rr) {
      CN[rr] = CN[lastRadius];
    }
  }

  if (debug) {
    BOOST_LOG(rdDebugLog) << "Final CN values for " << smi << ":\n";
    for (const auto &[r, values] : CN) {
      BOOST_LOG(rdDebugLog) << "Radius " << r << ": ";
      for (const auto &val : values[0]) {
        BOOST_LOG(rdDebugLog) << val << " ";
      }
      BOOST_LOG(rdDebugLog) << std::endl;
    }
  }

  // clean up the memory!
  // M.clear();
  // SP.clear();
  // clustersByAN.clear();
  // clusters.clear();

  return CN;
}

std::vector<double> ShannonEntropies(
    std::map<int, std::vector<std::vector<int>>> CN, int maxradius,
    double log2nA, double log2nB, int nAtoms) {
  const auto *tbl = PeriodicTable::getTable();

  int singleoutputsize = maxradius + 1;

  std::vector<double> icvalues(7 * singleoutputsize, 0.);

  for (auto &[r, data] : CN) {
    icvalues[r] = InfoEntropy(data[0]);                     // IC
    icvalues[r + singleoutputsize] = icvalues[r] * nAtoms;  // TIC
    if (log2nA > 0) {
      icvalues[r + 2 * singleoutputsize] = icvalues[r] / log2nA;  // SIC
    }

    if (log2nB > 0) {
      icvalues[r + 3 * singleoutputsize] = icvalues[r] / log2nB;  // BIC
    }

    icvalues[r + 4 * singleoutputsize] = log2nA - icvalues[r];  // CIC

    std::vector<double> w(data[1].size(), 0.);

    for (unsigned int j = 0; j < data[1].size(); j++) {
      w[j] = tbl->getAtomicWeight(data[1][j]);  // catch Atomic mass
    }

    icvalues[r + 5 * singleoutputsize] = WeightedInfoEntropy(
        data[0], w);  // MIC  use the atomic mass weighted Shannon Entropy
    icvalues[r + 6 * singleoutputsize] = WeightedCrossInfoEntropy(
        data[0],
        data[1]);  // ZMIC  use the atomic mass weighted cross Shannon Entropy
    w.clear();
  }

  return icvalues;
}

std::vector<double> calcInformationContent(const ROMol &mol, int maxradius) {
  std::unique_ptr<RWMol> hmol(new RWMol(mol));
  MolOps::addHs(*hmol);

  int nAtoms = hmol->getNumAtoms();

  if (nAtoms == 0) {
    BOOST_LOG(rdWarningLog) << "Error: Molecule has no atoms after adding hydrogens."
              << std::endl;
    return {};
  }

  double nBonds = 0.;
  for (auto &bond : hmol->bonds()) {
    nBonds += getbondtypeindouble(bond->getBondType());
  }

  double log2nA = std::log(static_cast<double>(nAtoms)) / std::log(2);
  double log2nB = std::log(static_cast<double>(nBonds)) / std::log(2);

  auto CN = computePipeline(*hmol, maxradius);

  if (CN.empty()) {
    BOOST_LOG(rdWarningLog) << "Error: ComputePipeline returned empty CN." << std::endl;
    return {};
  }

  return ShannonEntropies(CN, maxradius, log2nA, log2nB, nAtoms);
}

std::vector<double> calcInformationContent_(const ROMol &mol) {
  int maxradius = 5;
  // Dynamically allocate RWMol using new
  RWMol *hmol = new RWMol(mol);

  try {
    // Add hydrogens
    MolOps::addHs(*hmol);

    int nAtoms = hmol->getNumAtoms();
    if (nAtoms == 0) {
      BOOST_LOG(rdWarningLog) << "Error: Molecule has no atoms after adding hydrogens."
                << std::endl;
      delete hmol;  // Clean up memory
      return {};
    }

    double nBonds = 0.0;
    for (auto &bond : hmol->bonds()) {
      nBonds += getbondtypeindouble(bond->getBondType());
    }

    double log2nA = std::log(static_cast<double>(nAtoms)) / std::log(2);
    double log2nB =
        (nBonds > 0) ? std::log(static_cast<double>(nBonds)) / std::log(2) : 0;

    auto CN = computePipeline(*hmol, maxradius);
    if (CN.empty()) {
      BOOST_LOG(rdWarningLog) << "Error: ComputePipeline returned empty CN." << std::endl;
      delete hmol;  // Clean up memory
      return {};
    }

    // Calculate Shannon Entropies
    std::vector<double> icvalues =
        ShannonEntropies(CN, maxradius, log2nA, log2nB, nAtoms);

    delete hmol;  // Clean up memory
    return icvalues;

  } catch (const std::exception &e) {
    BOOST_LOG(rdWarningLog) << "Error: " << e.what() << std::endl;
    delete hmol;  // Clean up memory
    return {};
  }
}

// triplet example AZ

std::vector<double> TIn(const ROMol &mol, const std::vector<double> b) {
  double ti1 = std::accumulate(b.begin(), b.end(), 0.0);
  double ti2 =
      std::accumulate(b.begin(), b.end(), 0.0,
                      [](double acc, double yi) { return acc + yi * yi; });
  double ti3 = std::accumulate(
      b.begin(), b.end(), 0.0,
      [](double acc, double yi) { return acc + std::sqrt(yi); });

  double ti4 = 0;
  for (unsigned int i = 0; i < b.size(); ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      const auto *bond = mol.getBondBetweenAtoms(i, j);
      if (bond != nullptr) {
        ti4 += std::pow(b[i] * b[j], -0.5);
      }
    }
  }

  double ti5 = b.size() * std::pow(std::accumulate(b.begin(), b.end(), 1.0,
                                                   std::multiplies<double>()),
                                   1.0 / b.size());

  return {ti1, ti2, ti3, ti4, ti5};
}

// triplet example AZ

std::vector<double> TIn(const ROMol &mol, const std::vector<double> b,
                        int nhrs = 1) {
  int n = mol.getNumAtoms();
  std::vector res(5 * nhrs, 0.0);
  for (int i = 0; i < nhrs; i++) {
    std::vector<double> bi(b.begin() + i * n, b.begin() + (i + 1) * n);

    double ti1 = std::accumulate(bi.begin(), bi.end(), 0.0);
    double ti2 =
        std::accumulate(bi.begin(), bi.end(), 0.0,
                        [](double acc, double yi) { return acc + yi * yi; });
    double ti3 = std::accumulate(
        bi.begin(), bi.end(), 0.0,
        [](double acc, double yi) { return acc + std::sqrt(yi); });

    double ti4 = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        const auto *bond = mol.getBondBetweenAtoms(i, j);
        if (bond != nullptr) {
          ti4 += std::pow(bi[i] * bi[j], -0.5);
        }
      }
    }

    double ti5 = n * std::pow(std::accumulate(bi.begin(), bi.end(), 1.0,
                                              std::multiplies<double>()),
                              1.0 / n);
    res[i * 5 + 0] = ti1;
    res[i * 5 + 1] = ti2;
    res[i * 5 + 2] = ti3;
    res[i * 5 + 3] = ti4;
    res[i * 5 + 4] = ti5;
  }
  return res;
}

// triplet AN*x = V :  S,V,Z,I,N
std::vector<double> calcANMat(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A
  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // Modify diagonal to include atomic numbers

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] =
        static_cast<double>(nAtoms);  // Atomic number on diagonal  == Z
  }
  int n = nAtoms, nrhs = 5;

  std::vector<double> V(nAtoms * nrhs, 0.0);  // B matrix for 5 right-hand sides

  for (int j = 0; j < nrhs; ++j) {      // Iterate over columns
    for (int i = 0; i < nAtoms; ++i) {  // Iterate over rows
      switch (j) {
        case 0:
          V[i + j * nAtoms] = DistSum[i];
          break;  // S
        case 1:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getDegree());
          break;  // V
        case 2:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum());
          break;  // Z
        case 3:
          V[i + j * nAtoms] = static_cast<double>(1);
          break;  // I
        case 4:
          V[i + j * nAtoms] = static_cast<double>(nAtoms);
          break;  // N
      }
    }
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);
  bool success = false;

  std::vector<double> B(V);

  solveLinearSystem(mol, A_flat, B, n, nrhs, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, B, nrhs);
}

// triplet AZ*x = V : V,S,N
std::vector<double> calcAZMat(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A
  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // Modify diagonal to include atomic numbers

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(
        atom->getAtomicNum());  // Atomic number on diagonal  == Z
  }
  int n = nAtoms, nrhs = 3;

  std::vector<double> V(nAtoms * nrhs, 0.0);  // B matrix for 5 right-hand sides

  for (int j = 0; j < nrhs; ++j) {      // Iterate over columns
    for (int i = 0; i < nAtoms; ++i) {  // Iterate over rows
      switch (j) {
        case 0:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getDegree());
          break;  // V
        case 1:
          V[i + j * nAtoms] = DistSum[i];
          break;  // S
        case 2:
          V[i + j * nAtoms] = static_cast<double>(nAtoms);
          break;  // N
      }
    }
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);
  bool success = false;

  std::vector<double> B(V);

  solveLinearSystem(mol, A_flat, B, n, nrhs, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, B, nrhs);
}

// triplet AS*x = V: N,V,Z,I
std::vector<double> calcASMat(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A
  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // Modify diagonal to include atomic numbers

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);  // Atomic number on diagonal  == Z
  }
  int n = nAtoms, nrhs = 4;

  std::vector<double> V(nAtoms * nrhs, 0.0);  // B matrix for 5 right-hand sides

  for (int j = 0; j < nrhs; ++j) {      // Iterate over columns
    for (int i = 0; i < nAtoms; ++i) {  // Iterate over rows
      switch (j) {
        case 0:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getDegree());
          break;  // V
        case 1:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum());
          break;  // Z
        case 2:
          V[i + j * nAtoms] = static_cast<double>(1);
          break;  // I
        case 3:
          V[i + j * nAtoms] = static_cast<double>(nAtoms);
          break;  // N
      }
    }
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);
  bool success = false;

  std::vector<double> B(V);

  solveLinearSystem(mol, A_flat, B, n, nrhs, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, B, nrhs);
}

// triplet DS*x = V: V,I,N,Z
std::vector<double> calcDSMat(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // Modify diagonal to include atomic numbers

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);  // Atomic number on diagonal  == Z
  }
  int n = nAtoms, nrhs = 4;

  std::vector<double> V(nAtoms * nrhs, 0.0);  // B matrix for 5 right-hand sides

  for (int j = 0; j < nrhs; ++j) {      // Iterate over columns
    for (int i = 0; i < nAtoms; ++i) {  // Iterate over rows
      switch (j) {
        case 0:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getDegree());
          break;  // V
        case 1:
          V[i + j * nAtoms] = static_cast<double>(1);
          break;  // I
        case 2:
          V[i + j * nAtoms] = static_cast<double>(nAtoms);
          break;  // N
        case 3:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum());
          break;  // Z
      }
    }
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);
  bool success = false;

  std::vector<double> B(V);

  solveLinearSystem(mol, A_flat, B, n, nrhs, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, B, nrhs);
}

// triplet DN2*x = V: S,I,N,Z
std::vector<double> calcDN2Mat(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // Modify diagonal to include atomic numbers

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] =
        nAtoms * nAtoms;  // Atomic number on diagonal  == Z
  }
  int n = nAtoms, nrhs = 4;

  std::vector<double> V(nAtoms * nrhs, 0.0);  // B matrix for 5 right-hand sides

  for (int j = 0; j < nrhs; ++j) {      // Iterate over columns
    for (int i = 0; i < nAtoms; ++i) {  // Iterate over rows
      switch (j) {
        case 0:
          V[i + j * nAtoms] = DistSum[i];
          break;  // S
        case 1:
          V[i + j * nAtoms] = static_cast<double>(1);
          break;  // I
        case 2:
          V[i + j * nAtoms] = static_cast<double>(nAtoms);
          break;  // N
        case 3:
          V[i + j * nAtoms] =
              static_cast<double>(mol.getAtomWithIdx(i)->getAtomicNum());
          break;  // Z
      }
    }
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);
  bool success = false;

  std::vector<double> B(V);

  solveLinearSystem(mol, A_flat, B, n, nrhs, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, B, nrhs);
}

// triplet AZ*x = V

std::vector<double> calcAZV(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A
  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  // Modify diagonal to include atomic numbers

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(
        atom->getAtomicNum());  // Atomic number on diagonal  == Z
  }

  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    V[i] = static_cast<double>(atom->getDegree());  // Degree of the atom  == V
  }

  int n = nAtoms, nrhs = 1;

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, n, nrhs, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcASV(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<int> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<int>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);            // Distance Sum == S
    V[i] = static_cast<double>(atom->getDegree());  // Degree of the atom  == V
  }

  int n = nAtoms;

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, n, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDSV(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<int> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<int>(distances[i * nAtoms + j]);  // "S"
    }
  }
  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    DistMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);            // Distance Sum == S
    V[i] = static_cast<double>(atom->getDegree());  // Degree of the atom  == V
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcAZS(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(atom->getAtomicNum());  // Z
    // V[i] = static_cast<double>(atom->getDegree());   // Degree of the atom ==
    // V
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(DistSum);  // S

  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcASZ(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]);  // S
    V[i] = static_cast<double>(atom->getAtomicNum());             // Z
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // Z
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDN2S(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // D
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }
  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] = nAtoms * nAtoms;  // N2
    // V[i] = static_cast<double>(atom->getDegree());   // Degree of the atom ==
    // V
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(DistSum);  // S
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDN2I(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // D
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] = nAtoms * nAtoms;  // N2
    V[i] = static_cast<double>(1);                 // I
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);  // I
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcASI(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]);  // S
    V[i] = static_cast<double>(1);                                // I
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // Z
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDSI(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<int> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<int>(distances[i * nAtoms + j]);  // "S"
    }
  }
  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);  // Distance Sum == S
    V[i] = static_cast<double>(1);        // I
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcASN(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] = static_cast<double>(DistSum[i]);  // S
    V[i] = static_cast<double>(nAtoms);                           // N
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // Z
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDSN(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<int> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<int>(distances[i * nAtoms + j]);  // "S"
    }
  }
  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);  // Distance Sum == S
    V[i] = static_cast<double>(nAtoms);   // N
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDN2N(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // D
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    DistMatVec[i * nAtoms + i] = nAtoms * nAtoms;  // N2
    V[i] = static_cast<double>(nAtoms);            // N
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);  // I
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcANS(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"

  std::vector<double> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<double>(distances[i * nAtoms + j]);  // "S"
    }
  }

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms);  // N
    V[i] = static_cast<double>(DistSum[i]);                   // S
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // Z
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcANV(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms);  // N
    V[i] = static_cast<double>(atom->getDegree());            //  V
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // V

  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcAZN(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A
  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(atom->getAtomicNum());  // Z
    V[i] = static_cast<double>(nAtoms);  //  N
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // N
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcANZ(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms);  // N
    V[i] = static_cast<double>(atom->getAtomicNum());         // Z
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // Z
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcANI(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms);  // N
    V[i] = static_cast<double>(1);                            // I
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // I
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDSZ(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<int> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<int>(distances[i * nAtoms + j]);  // "S"
    }
  }
  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    DistMatVec[i * nAtoms + i] =
        static_cast<double>(DistSum[i]);               // Distance Sum == S
    V[i] = static_cast<double>(atom->getAtomicNum());  // Z
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcANN(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> adjMatVec(nAtoms * nAtoms, 0.0);  // A

  double *Mat = MolOps::getAdjacencyMatrix(mol, false, false, false, "NoBO");
  std::copy(Mat, Mat + (nAtoms * nAtoms), adjMatVec.begin());

  // std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());
  //  Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    adjMatVec[i * nAtoms + i] = static_cast<double>(nAtoms);  // N
    V[i] = static_cast<double>(nAtoms);                       // N
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(adjMatVec);

  std::vector<double> b(V);  // I
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

std::vector<double> calcDN2Z(const ROMol &mol) {
  int nAtoms = rdcast<int>(mol.getNumAtoms());

  // Use RDKit's built-in function to get the adjacency matrix
  std::vector<double> DistMatVec(nAtoms * nAtoms, 0.0);  // "D"
  double *distances = MolOps::getDistanceMat(
      mol, false, false, false);  // no need for "Bond order"
  std::copy(distances, distances + (nAtoms * nAtoms), DistMatVec.begin());

  std::vector<int> DistSum(nAtoms, 0.);

  for (int i = 0; i < nAtoms; ++i) {
    for (int j = 0; j < nAtoms; ++j) {
      DistSum[i] += static_cast<int>(distances[i * nAtoms + j]);  // "S"
    }
  }
  // Modify diagonal to include atomic numbers
  std::vector<double> V(nAtoms, 0.0);

  for (int i = 0; i < nAtoms; ++i) {
    const auto *atom = mol.getAtomWithIdx(i);
    DistMatVec[i * nAtoms + i] =
        static_cast<double>(nAtoms * nAtoms);          // Distance Sum == S
    V[i] = static_cast<double>(atom->getAtomicNum());  // Z
  }

  // Copy the adjacency matrix into a working buffer
  std::vector<double> A_flat(DistMatVec);

  std::vector<double> b(V);
  bool success = false;

  solveLinearSystem(mol, A_flat, b, nAtoms, 1, success);

  if (!success) {
    return {0., 0., 0., 0., 0.};
  }

  return TIn(mol, b, 1);
}

static const std::vector<std::string> frags = {
    "[CX4H3]",
    "[CX4H2]",
    "[CX4H1]",
    "[CX4H0]",
    "*=[CX3H2]",
    "[$(*=[CX3H1]),$([cX3H1](a)a)]",
    "[$(*=[CX3H0]),$([cX3H0](a)(a)A)]",
    "c(a)(a)a",
    "*#C",
    "[C][NX3;H2]",
    "[c][NX3;H2]",
    "[C][NX3;H1][C]",
    "[c][NX3;H1]",
    "[c][nX3;H1][c]",
    "[C][NX3;H0](C)[C]",
    "[c][NX3;H0](C)[C]",
    "[c][nX3;H0][c]",
    "*=[Nv3;!R]",
    "*=[Nv3;R]",
    "[nX2H0,nX3H1+](a)a",
    "N#C[A;!#1]",
    "N#C[a;!#1]",
    "[$([A;!#1][NX3](=O)=O),$([A;!#1][NX3+](=O)[O-])]",
    "[$([a;!#1][NX3](=O)=O),$([a;!#1][NX3+](=O)[O-])]",
    "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]",
    "[OH]",
    "[OX2;H0;!R]",
    "[OX2;H0;R]",
    "[oX2](a)a",
    "*=O",
    "[SX2](*)*",
    "[sX2](a)a",
    "*=[SX1]",
    "*=[SX3]",
    "[$([#16X4](=[OX1])(=[OX1])([!#8])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([!#8])[OX2H0])]",
    "[S,s]",
    "[P,p]",
    "FA",
    "Fa",
    "Cl",
    "Br",
    "I",
    "[CX3;!R](=[OX1])[OX2H0]",
    "[CX3;R](=[OX1])[OX2H0;R]",
    "P(=[OX1])(O)(O)O",
    "[CX3](=[OX1])([OX2H0])[OX2H0]",
    "[CX3](=O)[OX1H0-,OX2H1]",
    "nC=[OX1]",
    "[N;!R]C=[OX1]",
    "[N;R][C;R]=[OX1]",
    "[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]",
    "NC(=[OX1])N",
    "[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]",
    "[CX3](=[OX1])[NX3][CX3](=[OX1])",
    "C1(=[OX1])C=CC(=[OX1])C=C1",
    "[$([CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
    "*1~*2~*(~*3~*(~*~*~*~*3)~*1)~*~*~*1~*2~*~*~*1",
    "[OX2;H1]CC[O,N]",
    "[OX2;H1]C[C,N]=[O,S]",
    "[OX2;H1]c1ccccc1[O,NX3]",
    "[OX2;H1]c1ccccc1C=[O,S]",
    "[OX2;H1]c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "[NH,NH2,NH3+]CC[O,N]",
    "[NH,NH2,NH3+]c1ccccc1[O,N]",
    "[NH,NH2,NH3+]c1ccccc1[C,N]=[O,S]",
    "[OX2H]c1ccccc1[Cl,Br,I]",
    "[CX4]([OH])[CX4][OH]",
    "n:n",
    "o:n",
    "n:c:n",
    "o:c:n",
    "n:c:c:n",
    "[F,Cl,Br,I,N,O,S]-c:c-[F,Cl,Br,I,N,O,S]",
    "[F,Cl,Br,I,N,O,S]-c:c:c-[F,Cl,Br,I,N,O,S]",
    "[F,Cl,Br,I,N,O,S]-c:c:c:c-[F,Cl,Br,I,N,O,S]",
    "P(=[OX1])N",
    "Nc:n",
    "[$(cC[OH]);!$(c[CX3](=O)[OX1H0-,OX2H1])]",
    "[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]",
    "[OX2]-c:c-[OX2]",
    "[C][OX2H]",
    "[c][OX2H]",
    "[C][NX3;H1;!R][C]",
    "[C][NX3;H1;R][C]",
    "[c][NX3;H1;!$(NC=O)][C]",
    "[CX3](=[OX1])[NX3H2]",
    "[CX3](=[OX1])[NX3;H1][C]",
    "[CX3](=[OX1])[NX3;H1][c]",
    "[$([SX4](=[OX1])(=[OX1])([!O])[NH,NH2,NH3+]),$([SX4+2]([OX1-])([OX1-])([!O])[NH,NH2,NH3+])]",
    "[NX3;H1]C(=[OX1])[NX3;H1]",
    "[NX3;H0]C(=[OX1])[NX3;H1]",
    "[NX3;H1]C(=[OX1])O",
    "[NX3;H1]C(=N)[NX3;H0]",
    "[C]#[CH]",
    "P[OH,O-]",
    "[CH][F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
    "[CH]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)]",
    "[CX4]([CX3](=O)[OX1H0-,OX2H1])[CX4][CX3](=O)[OX1H0-,OX2H1]"
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX3](=O)[OX1H0-,OX2H1]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[OH]",
    "[CX4]([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])[CX4][OH]",
    "[nX3;H1]:n",
    "[nX3;H1]:c:n",
    "[OX1]=[C,c]~[C,c]C[OH]",
    "[OH]c1cccc2cccnc12",
    "[OH]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
    "[OH]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
    "[NH,NH2,NH3+]c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
    "[NH,NH2,NH3+]c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
    "[CX3](=O)([OX1H0-,OX2H1])c1cc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])ccc1",
    "[CX3](=O)([OX1H0-,OX2H1])c1ccc([F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cc1",
    "[OH]c1c([CX4])cccc1[CX4]",
    "[NH,NH2,NH3+]c1c([CX4])cccc1[CX4]",
    "[OH]c1c(C[F,Cl,Br,I,$([NX3](=O)=O),$([NX3+](=O)[O-]),$(C#N),$([CX4](F)(F)F)])cccc1",
    "[OH]c1cc([CX3](=O)[OX1H0-,OX2H1])ccc1",
    "[OH]c1ccc([CX3](=O)[OX1H0-,OX2H1])cc1",
    "[OH]c1cc([$([CH](=O)),$(C(=O)C)])ccc1",
    "[OH]c1ccc([$([CH](=O)),$(C(=O)C)])cc1",
    "[OX2H0+0]-[cX3H0;$(*-A)]:[cX3H0;$(*-a)]-[cX3H0;$(*-a)]:[cX3H1]:[cX3H1]",
    "[CX4H3]-[nX3H0+0;$(*-A)]",
    "[FX1H0]-[CX4H1](-[nX3H0+0;$(*-A)]1:[cX3H0;$(*-A)](-[CX4H3]):[nX2H0+0;!$(*-a);!$(*~A)]:[nX3H0+0;$(*-a)](:[cX3H0;$(*=A)]:1=[OX1H0+0])-[cX3H0;$(*-a)]1:[cX3H1]:[cX3H0;$(*-A)](:[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]:1-[ClX1H0])-[NX3H1+0]-[SX4H0](=[OX1H0+0])(=[OX1H0+0])-[CX4H3])-[FX1H0]",
    "[cX3H1]:[cX3H1]:[cX3H1]:[cX3H0;$(*-A)]-[CX3H1]=[CX3H1]-[CX2H0]#[NX1H0+0]",
    "[OX2H1+0]-[cX3H0;$(*-A)](:[cX3H1]):[cX3H1]",
    "[CX4H3]-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]",
    "[CX4H3]-[cX3H0;$(*-A)]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H0;!$(*-a);!$(*~A)]",
    "[CX4H3]-[OX2H0+0]-[cX3H0;$(*-A)]",
    "[CX4H2]-[SX2H0]",
    "[CX4H2]-[CX4H3]",
    "[cX3H1]:[cX3H1]:[cX3H0;!$(*-a);!$(*~A)](:[cX3H0;!$(*-a);!$(*~A)]):[cX3H0;!$(*-a);!$(*~A)]",
    "[cX3H0;!$(*-a);!$(*~A)]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H1]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H0;!$(*-a);!$(*~A)]:[cX3H1]:[cX3H1]",
    "[nX2H0+0;!$(*-a);!$(*~A)]:[nX3H0+0;$(*-A)]",
    "[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[nX2H0+0;!$(*-a);!$(*~A)]",
    "[CX4H3]-[SX2H0]",
    "[CX4H3]-[NX3H0+0]-[CX4H3]",
    "[CX4H2]-[CX4H2]-[cX3H0;$(*-A)]",
    "[cX3H1]:[cX3H0;$(*-A)]-[CX2H0]#[NX1H0+0]",
    "[OX2H0+0]-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H1]:[cX3H1]",
    "[CX4H1]-[CX4H2]",
    "[CX4H1]-[cX3H0;$(*-A)]",
    "[cX3H1]:[cX3H0;$(*-A)]:[nX2H0+0;!$(*-a);!$(*~A)]",
    "[cX3H1]:[cX3H1]:[cX3H0;$(*-A)]-[CX3H0]=[OX1H0+0]",
    "[NX3H1+0]-[CX3H0]=[OX1H0+0]",
    "[FX1H0]-[CX4H0](-[FX1H0])-[FX1H0]",
    "[cX3H1]:[cX3H0;$(*-a)](:[cX3H1])-[cX3H0;$(*-a)](:[cX3H1]):[cX3H1]",
    "[NX3H1+0]-[cX3H0;$(*-A)](:[cX3H1]):[cX3H1]",
    "[cX3H1]:[cX3H1]:[cX3H0;$(*-A)]-[NX3H0+0](=[OX1H0+0])=[OX1H0+0]",
    "[CX4H3]-[cX3H0;$(*-A)]:[cX3H1]:[cX3H1]:[cX3H1]",
    "[CX4H3]-[CX4H0]-[CX4H3]",
    "[CX4H1]-[CX4H1]",
    "[CX4H2]-[CX3H0](=[OX1H0+0])-[NX3H0+0]",
    "[NX3H0+0]-[cX3H0;$(*-A)](:[cX3H0;$(*-A)]):[cX3H0;$(*-A)]",
    "[OX2H1+0]-[cX3H0;$(*-A)](:[cX3H1]):[cX3H0;$(*-A)]-[ClX1H0]",
    "[CX4H2]-[OX2H0+0]-[CX3H0]=[OX1H0+0]",
    "[CX4H3]-[OX2H0+0]-[PX4H0](=[SX1H0])-[OX2H0+0]-[CX4H3]",
    "[OX2H0+0]-[CX3H0](=[OX1H0+0])-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]",
    "[NX3H0+0]-[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]",
    "[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]",
    "[CX4H2]-[cX3H0;$(*-A)]:[cX3H1]:[cX3H0;$(*-A)]",
    "[CX4H0]-[CX4H0]-[ClX1H0]",
    "[cX3H0;$(*-A)]1:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:[cX3H0;$(*-A)]:1",
    "[A!#1x0+0]#[A!#1x0+0]",
    "[SX2H0]",
    "[ax3+0;$(*-[A!#1])]",
    "[NX3H0+0]",
    "[CX3H1]=[CX3H2]",
    "[CX3H1]=[OX1H0+0]",
    "[cX3H0;$(*-A)]",
    "[ax3+0;$(*-a)]",
    "[nX3H1+0]",
    "[OX2H0+0]",
    "[A!#1x0+0]",
    "[#8]",
    "[cX3H0;!$(*-a);!$(*~A)]",
    "[#7]",
    "[#6]",
    "[SX2H1]",
    "[CX3](=O)[OX2H1]",
    "[$([CX3H][#6]),$([CX3H2])]=[OX1]",
    "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]",
    "[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]",
    "[CX4](F)(F)F",
    "[NX3]=[CX3]",
    "[NX3][CX3]=[NX3]",
    "[NX1]#[CX2]",
    "[CX3]=[OX1]",
    "[#6][CX3](=O)[#6]",
    "[CX3H1](=O)[#6]",
    "[NX3][CX3](=[OX1])[#6]",
    "[NX3][CX3](=[OX1])[#5]",
    "[NX3][CX3](=[OX1])[OX2H0]",
    "[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]",
    "[#6][CX3](=[OX1])[OX2H0][#6]",
    "[CX3](=[OX1])[OX1-]",
    "[CX3](=O)[OX2H1]",
    "[OX1]=[CX3]([OX2])[OX2]",
    "[CX3]=[SX1]",
    "[NX3][NX3]",
    "[NX2]=N",
    "[NX2]=[OX1]",
    "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
    "[OX1]=[NX2][OX2]",
    "[OX2,OX1-][OX2,OX1-]",
    "[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]",
    "[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]",
    "[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]",
    "[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]",
    "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]",
    "[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]",
    "[#16X2H0][#16X2H0]",
    "[PX5](=[OX1])([OX1-])[OX1-]",
    "[PX5](=[OX1])([OX2H])[OX2H]",
    "[PX6](=[OX1])([OX1-])([OX1-])[OX1-]"};

static const std::vector<std::shared_ptr<RWMol>> GetQueriesFrags() {
  static const std::vector<std::shared_ptr<RWMol>> queriesFrags = [] {
    std::vector<std::shared_ptr<RWMol>> res;
    for (const auto &smi : frags) {
      auto mol = SmartsToMol(smi);
      if (mol) {
        res.emplace_back(std::move(mol));
      } else {
        BOOST_LOG(rdWarningLog) << "Invalid SMARTS: " << smi << std::endl;
      }
    }
    return res;
  }();
  return queriesFrags;
}

std::vector<double> calcFrags(const ROMol &mol) {
  auto queriesFrags = GetQueriesFrags();
  std::vector<double> retval(queriesFrags.size(), 0.0);

  try {
    // Calculate A descriptor
    for (size_t i = 0; i < queriesFrags.size(); ++i) {
      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *queriesFrags[i], matches, true);  // uniquify = true
      retval[i] = matches.size();
    }

  } catch (const std::exception &e) {
    BOOST_LOG(rdWarningLog) << "Error in SMARTS matching: " << e.what() << std::endl;
    throw std::runtime_error("Error in SMARTSQueryTool");
  }

  return retval;
}

}  // namespace Osmordred
}  // namespace Descriptors
}  // namespace RDKit
