//
//  Copyright (c) 2016, Guillaume GODIN
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
//     * Neither the name of Institue of Cancer Research.
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
// for build & set RDBASE! => export RDBASE=/Users/GVALMTGG/Github/rdkit_mine/

#include <GraphMol/RDKitBase.h>

#include "WHIM.h"
#include "MolData3Ddescriptors.h"

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace Eigen;

namespace RDKit {
namespace Descriptors {
namespace {

MolData3Ddescriptors moldata3D;

double roundn(double in, int factor) {
  return std::round(in * pow(10., factor)) / pow(10., factor);
}

MatrixXd GetCenterMatrix(MatrixXd &Mat) {
  VectorXd v = Mat.colwise().mean();
  MatrixXd X = Mat.rowwise() - v.transpose();
  return X;
}

MatrixXd GetCovMatrix(MatrixXd &X, MatrixXd &Weight, double weight) {
  return X.transpose() * Weight * X / weight;
}

JacobiSVD<MatrixXd> *getSVD(MatrixXd &Mat) {
  auto *svd = new JacobiSVD<MatrixXd>(Mat, ComputeThinU | ComputeThinV);
  return svd;
}

std::vector<double> getWhimD(std::vector<double> weightvector,
                             MatrixXd MatOrigin, int numAtoms, double th) {
  double *weightarray = &weightvector[0];

  Map<VectorXd> Weight(weightarray, numAtoms);
  // std::cerr << "Weight:\n" << Weight << "\n";

  MatrixXd WeightMat = Weight.asDiagonal();

  double weight = WeightMat.diagonal().sum();
  // fix issue if the sum is close to zeros
  // only for the charges cases normally
  if (fabs(weight) < 1e-4) {
    weight = 1.0;
    // std::cerr << "fix weight sum:\n";
  }

  MatrixXd Xmean = GetCenterMatrix(MatOrigin);
  // std::cerr << "Xmean:\n" << Xmean << "\n";

  MatrixXd covmat = GetCovMatrix(Xmean, WeightMat, weight);

  JacobiSVD<MatrixXd> *svd = getSVD(covmat);

  std::vector<double> w(18);
  // prepare data for Whim parameter computation

  const double *SingVal = svd->singularValues().data();
  MatrixXd Scores = Xmean * svd->matrixV();  //  V is similar

  // compute parameters
  w[0] = SingVal[0];
  w[1] = SingVal[1];
  w[2] = SingVal[2];
  w[3] = SingVal[0] + SingVal[1] + SingVal[2];  // T
  w[4] = SingVal[0] * SingVal[1] + SingVal[0] * SingVal[2] +
         SingVal[1] * SingVal[2];                              // A
  w[5] = w[3] + w[4] + SingVal[0] * SingVal[1] * SingVal[2];   // V
  w[6] = SingVal[0] / (SingVal[0] + SingVal[1] + SingVal[2]);  // P1
  w[7] = SingVal[1] / (SingVal[0] + SingVal[1] + SingVal[2]);  // p2
  w[8] = SingVal[2] / (SingVal[0] + SingVal[1] + SingVal[2]);  // P3

  double res = 0.0;
  for (int i = 0; i < 3; i++) {
    res += std::fabs(w[i] / w[3] - 1.0 / 3.0);
  }

  w[9] = 3.0 / 4.0 * res;  // K

  // center original matrix
  VectorXd v1 = Scores.col(0);
  VectorXd v2 = Scores.col(1);
  VectorXd v3 = Scores.col(2);

  //  inverse of the kurtosis
  if (v1.array().pow(4).sum() > 0) {
    w[10] = numAtoms * pow(w[0], 2) / v1.array().pow(4).sum();  // E1
  } else {
    w[10] = 0.0;
  }

  if (v2.array().pow(4).sum() > 0) {
    w[11] = numAtoms * pow(w[1], 2) / v2.array().pow(4).sum();  // E2
  } else {
    w[11] = 0.0;
  }

  if (v3.array().pow(4).sum() > 0) {
    w[12] = numAtoms * pow(w[2], 2) / v3.array().pow(4).sum();  // E3
  } else {
    w[12] = 0.0;
  }

  w[13] = (w[10] + w[11] + w[12]) / 3.0;  // mean total density of the atoms
                                          // called D is used on Dragon 6 not
                                          // just the sum!

  // check if the molecule is fully symmetrical "like CH4" using Canonical Rank
  // Index and/or Sphericity !

  double gamma[3];  // Gamma values
  auto nAT = (double)numAtoms;

  // check if two atoms are symmetric versus the new axis ie newx,newy,newz a
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < numAtoms; j++) {
      Scores(j, i) = roundn(Scores(j, i),
                            3);  // round the matrix! same as eigen tolerance !
    }
  }

  //  we should take into account atoms that are in the axis too!!! which is not
  //  trivial
  for (int i = 0; i < 3; i++) {
    std::vector<double> Symmetric(2 * numAtoms, 0.0);
    double ns = 0.0;
    double na = 0.0;
    for (int j = 0; j < numAtoms; j++) {
      bool amatch = false;
      for (int k = 0; k < numAtoms; k++) {
        if (j == k) {
          continue;
        }
        // those that are close opposite & not close to the axis!
        if (std::fabs(Scores(j, i) + Scores(k, i)) <= th) {
          // check only once the symmetric none null we need to add +2!
          // (reduce the loop duration)
          ns += 1;
          amatch = true;
          Symmetric[j] = 1.0;
          Symmetric[j + numAtoms] = 2.0;
          Symmetric[k] = 1.0;
          Symmetric[k + numAtoms] = 2.0;
          break;
        }
      }
      if (!amatch) {
        na += 1;
        Symmetric[j] = 0.0;
        Symmetric[j + numAtoms] = std::fabs(Scores(j, i));
      }
    }
    // take into account the atoms close to the axis
    for (int aj = 0; aj < numAtoms; aj++) {
      if (Symmetric[aj + numAtoms] < th && Symmetric[aj] < 1.0) {
        ns += 1;
        na -= 1;
      }
    }
    gamma[i] = 0.0;
    double gammainv = 1.0;
    if (ns == 0) {
      gammainv = 1.0 - (na / nAT) * log(1.0 / nAT) / log(2.);
    }
    if (ns > 0) {
      gammainv = 1.0 - ((ns / nAT) * log(ns / nAT) / log(2.) +
                        (na / nAT) * log(1.0 / nAT) / log(2.));
    }
    gamma[i] = 1.0 / gammainv;
  }
  w[14] = gamma[0];  // G1
  w[15] = gamma[1];  // G2
  w[16] = gamma[2];  // G3
  w[17] = pow(gamma[0] * gamma[1] * gamma[2], 1.0 / 3.0);
  delete svd;

  return w;
}

void GetWHIMs(const Conformer &conf, std::vector<double> &result,
              double *Vpoints, double th) {
  std::vector<double> wu(18);
  std::vector<double> wm(18);
  std::vector<double> wv(18);
  std::vector<double> we(18);
  std::vector<double> wp(18);
  std::vector<double> wi(18);
  std::vector<double> ws(18);

  int numAtoms = conf.getNumAtoms();
  Map<MatrixXd> matorigin(Vpoints, 3, numAtoms);
  MatrixXd MatOrigin = matorigin.transpose();
  std::vector<double> weightvector;

  // intermediate 18 values stored in this order per weighted vector :
  // "L1","L2","L3","T","A","V","P1","P2","P3","K","E1","E2","E3","D","G1","G2","G3","G"
  weightvector = moldata3D.GetUn(numAtoms);
  wu = getWhimD(weightvector, MatOrigin, numAtoms, th);

  weightvector = moldata3D.GetRelativeMW(conf.getOwningMol());
  wm = getWhimD(weightvector, MatOrigin, numAtoms, th);

  weightvector = moldata3D.GetRelativeVdW(conf.getOwningMol());
  wv = getWhimD(weightvector, MatOrigin, numAtoms, th);

  weightvector = moldata3D.GetRelativeENeg(conf.getOwningMol());
  we = getWhimD(weightvector, MatOrigin, numAtoms, th);

  weightvector = moldata3D.GetRelativePol(conf.getOwningMol());
  wp = getWhimD(weightvector, MatOrigin, numAtoms, th);

  weightvector = moldata3D.GetRelativeIonPol(conf.getOwningMol());
  wi = getWhimD(weightvector, MatOrigin, numAtoms, th);

  weightvector = moldata3D.GetIState(conf.getOwningMol());
  ws = getWhimD(weightvector, MatOrigin, numAtoms, th);

  result.clear();
  result.resize(126);

  for (int i = 0; i < 18; i++) {
    result[i + 18 * 0] = wu[i];
    result[i + 18 * 1] = wm[i];
    result[i + 18 * 2] = wv[i];
    result[i + 18 * 3] = we[i];
    result[i + 18 * 4] = wp[i];
    result[i + 18 * 5] = wi[i];
    result[i + 18 * 6] = ws[i];
  }
  wu.clear();
  wm.clear();
  wv.clear();
  we.clear();
  wp.clear();
  wi.clear();
  ws.clear();
}

void GetWHIMsCustom(const Conformer &conf, std::vector<double> &result,
                    double *Vpoints, double th,
                    const std::string &customAtomPropName) {
  std::vector<double> wc(18);

  result.clear();
  result.resize(18);

  int numAtoms = conf.getNumAtoms();
  Map<MatrixXd> matorigin(Vpoints, 3, numAtoms);
  MatrixXd MatOrigin = matorigin.transpose();

  // intermediate 18 values stored in this order per weighted vector :
  // "L1","L2","L3","T","A","V","P1","P2","P3","K","E1","E2","E3","D","G1","G2","G3","G"
  std::vector<double> weightvector =
      moldata3D.GetCustomAtomProp(conf.getOwningMol(), customAtomPropName);

  wc = getWhimD(weightvector, MatOrigin, numAtoms, th);

  for (int i = 0; i < 18; i++) {
    result[i] = wc[i];
  }
  wc.clear();
}

void getWHIM(const ROMol &mol, std::vector<double> &res, int confId,
             double th) {
  int numAtoms = mol.getNumAtoms();
  const Conformer &conf = mol.getConformer(confId);
  auto *Vpoints = new double[3 * numAtoms];

  for (int i = 0; i < numAtoms; ++i) {
    Vpoints[3 * i] = conf.getAtomPos(i).x;
    Vpoints[3 * i + 1] = conf.getAtomPos(i).y;
    Vpoints[3 * i + 2] = conf.getAtomPos(i).z;
  }

  std::vector<double> w(126);
  GetWHIMs(conf, w, Vpoints, th);
  delete[] Vpoints;

  // Dragon extract only this list in this order : L1 L2 L3 P1 P2 G1 G2 G3 E1 E2
  // E3
  int map1[11] = {0, 1, 2, 6, 7, 14, 15, 16, 10, 11, 12};

  for (int k = 0; k < 7; k++) {
    for (int i = 0; i < 11; i++) {
      res[i + 11 * k] = roundn(w[map1[i] + 18 * k], 3);
    }
  }

  for (int i = 0; i < 2; i++) {
    res[i + 13 * 7] = roundn(w[17 + 18 * i], 3);  // 92  93 for Gu  Gm
  }

  for (int i = 0; i < 7; i++) {
    res[i + 11 * 7] =
        roundn(w[3 + 18 * i],
               3);  // 78  79  80  81  82  83  84  for Tu  Tm  Tv  Te  Tp  Ti Ts
    res[i + 12 * 7] =
        roundn(w[4 + 18 * i],
               3);  // 85  86  87  88  89  90  91  for Tu  Am  Av  Ae  Ap  Ai As
    res[i + 13 * 7 + 2] =
        roundn(w[9 + 18 * i],
               3);  // 94  95  96  97  98  99  100 for Ku  Km  Kv  Ke  Kp  Ki Ks
    res[i + 14 * 7 + 2] =
        roundn(w[13 + 18 * i],
               3);  // 101 102 103 104 105 106 107 for Du  Dm  Dv  De  Dp  Di Ds
    res[i + 15 * 7 + 2] =
        roundn(w[5 + 18 * i],
               3);  // 108 109 110 111 112 113 114 for Vu  Vm  Vv  Ve  Vp  Vi Vs
  }
}

void getWHIMone(const ROMol &mol, std::vector<double> &res, int confId,
                double th, const std::string &customAtomPropName) {
  int numAtoms = mol.getNumAtoms();
  const Conformer &conf = mol.getConformer(confId);
  auto *Vpoints = new double[3 * numAtoms];

  for (int i = 0; i < numAtoms; ++i) {
    Vpoints[3 * i] = conf.getAtomPos(i).x;
    Vpoints[3 * i + 1] = conf.getAtomPos(i).y;
    Vpoints[3 * i + 2] = conf.getAtomPos(i).z;
  }

  std::vector<double> w(18);
  GetWHIMsCustom(conf, w, Vpoints, th, customAtomPropName);
  delete[] Vpoints;

  // Dragon extract only this list in this order : L1 L2 L3 P1 P2 G1 G2 G3 E1 E2
  // E3
  // "L1","L2","L3","T","A","V","P1","P2","P3","K","E1","E2","E3","D","G1","G2","G3","G"
  int map1[17] = {0, 1, 2, 6, 7, 14, 15, 16, 10, 11, 12, 3, 4, 17, 9, 13, 5};

  for (int i = 0; i < 17; i++) {
    res[i] = roundn(w[map1[i]], 3);
  }
}

}  // end of anonymous namespace

void WHIM(const ROMol &mol, std::vector<double> &res, int confId, double th,
          const std::string &customAtomPropName) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  // Dragon final list is: L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u
  // Tu   Au    Gu   Ku    Du   Vu
  if (customAtomPropName != "") {
    res.clear();
    res.resize(17);
    getWHIMone(mol, res, confId, th, customAtomPropName);
  } else {
    res.clear();
    res.resize(114);
    getWHIM(mol, res, confId, th);
  }
}
}  // namespace Descriptors
}  // namespace RDKit
