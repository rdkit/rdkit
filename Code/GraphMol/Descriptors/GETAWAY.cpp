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
// Adding GETAWAY descriptors by Guillaume Godin
// for build & set RDBASE! => export RDBASE=/Users/mbp/Github/rdkit_mine/

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "GETAWAY.h"
#include "PBF.h"
#include "MolData3Ddescriptors.h"

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <Numerics/EigenSolvers/PowerEigenSolver.h>

#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <deque>
#include <Eigen/Core>
#include <Eigen/QR>

using namespace Eigen;

namespace RDKit {

namespace Descriptors {

namespace {

MolData3Ddescriptors moldata3D;

double roundn(double in, int factor) {
  return std::round(in * pow(10., factor)) / pow(10., factor);
}

double* retreiveMat(MatrixXd matrix) {
  double* arrayd = matrix.data();
  return arrayd;
}

double* retreiveVect(VectorXd matrix) {
  double* arrayd = matrix.data();
  return arrayd;
}

VectorXd getEigenVect(std::vector<double> v) {
  double* varray_ptr = &v[0];
  Map<VectorXd> V(varray_ptr, v.size());
  return V;
}

double round_to_n_digits_(double x, int n) {
  char buff[32];
  sprintf(buff, "%.*g", n, x);  // use g to have significant digits!
  return atof(buff);
}

double round_to_n_digits(double x, int n) {
  double scale = pow(10.0, ceil(log10(fabs(x))) + n);
  return std::round(x * scale) / scale;
}

bool IsClose(double a, double b, int n) {
  return (fabs(a - b) <= pow(0.1, n) * 1.1);
}

int countZeros(std::string ta) {
  int nbzero = 0;
  for (char i : ta) {
    if (i != '0') {
      break;
    }
    nbzero = nbzero + 1;
  }

  return nbzero;
}

bool IsClose2(double a, double b, unsigned int) {
  bool isclose = false;

  std::string sa, sb;
  std::string ta, tb;
  std::stringstream outa, outb;
  outa << a;
  sa = outa.str();
  outb << b;
  sb = outb.str();
  ta = sa.substr(sa.find('.') + 1);
  tb = sb.substr(sb.find('.') + 1);

  // same number of digits!
  if (countZeros(ta) == countZeros(tb)) {
    if (ta.length() > tb.length()) {
      tb += '0';
    }
    if (ta.length() < tb.length()) {
      ta += '0';
    }
  }
  // std::cout << tb << ":"<<  atoi(tb.c_str()) << ":" << countZeros(ta) <<
  // "\n";
  // std::cout << ta << ":"<<  atoi(ta.c_str()) << ":" << countZeros(ta) <<
  // "\n";

  if (ta.length() == tb.length()) {
    // last digit +/-1 deviation only!)
    if (abs(atoi(ta.c_str()) - atoi(tb.c_str())) < 2) {
      // std::cout << a << "," << b << " ta:"  << ta << "," << "tb:" << tb<<
      // "\n";
      isclose = true;
    }
  }
  return isclose;
}

bool IsClose3(double a, double b, unsigned int) {
  bool isclose = false;
  std::string sa, sb;
  std::string ta, tb;
  std::stringstream outa, outb;
  outa << a;
  sa = outa.str();
  outb << b;
  sb = outb.str();
  ta = sa.substr(sa.find('.') + 1);
  tb = sb.substr(sb.find('.') + 1);
  // same number of digits!
  if (ta.length() == tb.length()) {
    // last digit +/-1 deviation only!)
    if (abs(atoi(ta.c_str()) - atoi(tb.c_str())) < 2) {
      // std::cout << a << "," << b << " ta:"  << ta << "," << "tb:" << tb<<
      // "\n";
      isclose = true;
    }
  }
  return isclose;
}

// need to clean that code to have always the same output which is not the case
// the classes used the last 2 digit +/- 1 to cluster in the same group (with 3
// digit precision)
// 0.012 and 0.013 are in the same group while 0.014 are not!
std::vector<double> clusterArray(std::vector<double> data, double precision) {
  std::vector<double> Store;

  // sort the input data
  std::sort(data.begin(), data.end());

  // find the difference between each number and its predecessor
  std::vector<double> diffs;

  std::adjacent_difference(data.begin(), data.end(), std::back_inserter(diffs));

  // convert differences into percentage changes
  // std::transform(diffs.begin(), diffs.end(), data.begin(), diffs.begin(),
  //    std::divides<double>());

  int j = 0;
  int count = 0;
  for (unsigned int i = 0; i < data.size(); i++) {
    // std::cout << diffs[i] << ",";
    count++;
    // if a difference exceeds 0.01 <=> 1%, start a new group: if transform is
    // used!
    // use diff not ratio (with 0.003 precision look like it's what it's used in
    // Dragon 6!)
    if (diffs[i] > pow(0.1, precision)) {
      Store.push_back(count);
      count = 0;
      j++;
    }
  }

  return Store;
}

// need to clean that code to have always the same output which is not the case
// the classes used the last 2 digit +/- 1 to cluster in the same group (with 3
// digit precision)
// 0.012 and 0.013 are in the same group while 0.014 are not!
// alternative clustering method not working as in Dragon
std::vector<double> clusterArray2(std::vector<double> data,
                                  unsigned int precision) {
  std::vector<double> Store;

  // sort the input data descend order!
  std::sort(data.begin(), data.end(), std::greater<double>());
  std::deque<double> B(data.begin(), data.end());

  int count = 0;
  while (!B.empty()) {
    double dk = B.front();
    B.pop_front();
    count++;
    // case where the last value is not grouped it goes to the single last group
    if (B.empty()) {
      Store.push_back(count);
    }

    for (unsigned int i = 0; i < B.size(); i++) {
      if (IsClose2(dk, B[i], precision)) {
        count++;
      } else {
        Store.push_back(count);
        for (unsigned int j = 0; j < i; j++) {
          B.pop_front();
        }
        count = 0;
        break;
      }

      if (i == B.size() - 1) {
        Store.push_back(count);
        B.pop_front();
        count = 0;
      }
    }
  }
  return Store;
}

std::vector<double> GetGeodesicMatrix(const double* dist, int lag,
                                      int numAtoms) {
  PRECONDITION(dist != nullptr, "bad array");

  int sizeArray = numAtoms * numAtoms;
  std::vector<double> Geodesic;
  Geodesic.reserve(sizeArray);
  std::transform(dist, dist + sizeArray, Geodesic.begin(),
                 [lag](double dist) { return int(dist == lag); });

  return Geodesic;
}

JacobiSVD<MatrixXd> getSVD(const MatrixXd& A) {
  JacobiSVD<MatrixXd> mysvd(A, ComputeThinU | ComputeThinV);
  return mysvd;
}

MatrixXd GetPinv(const MatrixXd& A) {
  JacobiSVD<MatrixXd> svd = getSVD(A);
  double pinvtoler = 1.e-3;  // choose your tolerance wisely!
  VectorXd vs = svd.singularValues();
  VectorXd vsinv = svd.singularValues();

  for (unsigned int i = 0; i < A.cols(); ++i) {
    if (vs(i) > pinvtoler) {
      vsinv(i) = 1.0 / vs(i);
    } else {
      vsinv(i) = 0.0;
    }
  }

  MatrixXd S = vsinv.asDiagonal();
  MatrixXd Ap = svd.matrixV() * S * svd.matrixU().transpose();
  return Ap;
}

MatrixXd GetCenterMatrix(MatrixXd Mat) {
  VectorXd v = Mat.colwise().mean();
  MatrixXd X = Mat.rowwise() - v.transpose();
  return X;
}

MatrixXd GetHmatrix(MatrixXd X) {
  MatrixXd Weighted = X.transpose() * X;
  return X * GetPinv(Weighted) * X.transpose();
}

MatrixXd GetRmatrix(MatrixXd H, MatrixXd DM, int numAtoms) {
  MatrixXd R = MatrixXd::Zero(numAtoms, numAtoms);
  for (int i = 0; i < numAtoms - 1; i++) {
    for (int j = i + 1; j < numAtoms; j++) {
      R(i, j) = sqrt(H(i, i) * H(j, j)) / DM(i, j);
      R(j, i) = R(i, j);
    }
  }

  return R;
}

std::vector<int> GetHeavyList(const ROMol& mol) {
  int numAtoms = mol.getNumAtoms();
  std::vector<int> HeavyList;
  HeavyList.reserve(numAtoms);
  for (int i = 0; i < numAtoms; ++i) {
    const RDKit::Atom* atom = mol.getAtomWithIdx(i);
    HeavyList.push_back(int(atom->getAtomicNum() > 1));
  }
  return HeavyList;
}

double* AppendDouble(double* w, const double* Append, int length, int pos) {
  PRECONDITION(w != nullptr, "bad array");
  PRECONDITION(Append != nullptr, "bad array");

  for (int i = pos; i < pos + length; i++) {
    w[i] = Append[i - pos];
  }

  return w;
}

double getRCON(MatrixXd R, MatrixXd Adj, int numAtoms) {
  // similar implementation of J. Chem. Inf. Comput. Sci. 2004, 44, 200-209
  // equation 1 or 2 page 201
  // we use instead of atomic absolute values the atomic relative ones as in
  // Dragon
  double RCON = 0.0;
  VectorXd VSR = R.rowwise().sum();
  for (int i = 0; i < numAtoms - 1; ++i) {
    for (int j = i + 1; j < numAtoms; ++j) {
      if (Adj(i, j) > 0) {
        RCON +=
            sqrt(VSR(i) * VSR(j));  // the sqrt is in the sum not over the sum!
      }
    }
  }
  return RCON;
}

double getHATS(double W1, double W2, double H1, double H2) {
  return W1 * H1 * W2 * H2;
}

double getH(double W1, double W2, double H) { return W1 * H * W2; }

double getMax(const double* Rk) {
  PRECONDITION(Rk != nullptr, "bad rK");
  double RTp = 0;
  for (int j = 0; j < 8; j++) {
    if (Rk[j] > RTp) {
      RTp = Rk[j];
    }
  }
  return RTp;
}

void getGETAWAYDescCustom(MatrixXd H, MatrixXd R, MatrixXd Adj, int numAtoms,
                          std::vector<int> Heavylist, const ROMol& mol,
                          std::vector<double>& res, unsigned int precision,
                          const std::string& customAtomPropName) {
  // prepare data for Getaway parameter computation
  // compute parameters

  VectorXd Lev = H.diagonal();
  std::vector<double> heavyLev;

  for (int i = 0; i < numAtoms; i++) {
    if (Heavylist[i] == 1) {
      heavyLev.push_back(round_to_n_digits_(Lev(i), precision));
    }
  }

  std::vector<double> Clus = clusterArray2(heavyLev, precision);

  double numHeavy = heavyLev.size();
  double ITH0 = numHeavy * log(numHeavy) / log(2.0);
  double ITH = ITH0;
  for (double Clu : Clus) {
    ITH -= Clu * log(Clu) / log(2.0);
  }
  res[0] = roundn(ITH, 3);  // issue sometime with this due to cluster
  double ISH = ITH / ITH0;
  res[1] = roundn(ISH, 3);  // issue sometime with this due to cluster

  // use the PBF to determine 2D vs 3D (with Threshold)
  // determine if this is a plane molecule
  double pbf = RDKit::Descriptors::PBF(mol);
  double D;
  if (pbf < 1.e-5) {
    D = 2.0;
  } else {
    D = 3.0;
  }

  // what about linear molecule ie (D=1) ?
  double HIC = 0.0;
  for (int i = 0; i < numAtoms; i++) {
    HIC -= H(i, i) / D * log(H(i, i) / D) / log(2);
  }
  res[2] = roundn(HIC, 3);

  double HGM = 1.0;
  for (int i = 0; i < numAtoms; i++) {
    HGM = HGM * H(i, i);
  }
  HGM = 100.0 * pow(HGM, 1.0 / numAtoms);
  res[3] = roundn(HGM, 3);

  double RARS = R.rowwise().sum().sum() / numAtoms;

  JacobiSVD<MatrixXd> mysvd = getSVD(R);

  VectorXd EIG = mysvd.singularValues();

  double rcon = getRCON(R, std::move(Adj), numAtoms);

  std::vector<double> customAtomArray =
      moldata3D.GetCustomAtomProp(mol, customAtomPropName);
  VectorXd Wc = getEigenVect(customAtomArray);

  MatrixXd Bi;
  MatrixXd RBw;
  double HATSc, HATSct, H0ct, R0ct, H0c, R0c, Rkmaxc, tmpc;
  double HATSk[9];
  double Hk[9];
  double Rk[8];
  double Rp[8];

  double* dist =
      MolOps::getDistanceMat(mol, false);  // need to be be set to false to have
  // topological distance not weighted!

  Map<MatrixXd> D2(dist, numAtoms, numAtoms);

  double Dmax = D2.colwise().maxCoeff().maxCoeff();

  HATSct = 0.0;
  H0ct = 0.0;
  R0ct = 0.0;

  // need to loop other all the D values so <= not <!
  for (int i = 0; i <= Dmax; i++) {
    if (i == 0) {
      Bi = H.diagonal().asDiagonal();
    }

    HATSc = 0.0;
    H0c = 0.0;

    if (i == 0) {  // use Bi
      for (int j = 0; j < numAtoms; j++) {
        for (int k = j; k < numAtoms; k++) {
          if (Bi(j, k) > 0) {
            HATSc += getHATS((double)Wc(j), (double)Wc(j), (double)H(j, j),
                             (double)H(j, j));

            if (H(j, k) > 0) {
              H0c += getH((double)Wc(j), (double)Wc(k), (double)H(j, k));
            }
          }
        }
      }
    }
    auto Bimat = GetGeodesicMatrix(dist, i, numAtoms);
    Map<MatrixXd> Bj(Bimat.data(), numAtoms, numAtoms);
    if (i > 0 && i < 9) {  // use Bj
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            HATSc += getHATS((double)Wc(j), (double)Wc(k), (double)H(j, j),
                             (double)H(k, k));

            if (H(j, k) > 0) {
              H0c += getH((double)Wc(j), (double)Wc(k), (double)H(j, k));
            }
          }
        }
      }
    }

    if (i >= 9) {  // Totals missing part
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            HATSct += 2 * getHATS((double)Wc(j), (double)Wc(k), (double)H(j, j),
                                  (double)H(k, k));

            if (H(j, k) > 0) {
              H0ct += 2 * getH((double)Wc(j), (double)Wc(k), (double)H(j, k));
            }
          }
        }
      }
    }

    if (i < 9) {
      HATSk[i] = HATSc;
      Hk[i] = H0c;
    }

    R0c = 0.0;
    Rkmaxc = 0.0;

    // from 1 to 8
    if (i > 0 && i < 9) {
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            tmpc = getH((double)Wc(j), (double)Wc(k),
                        (double)R(j, k));  // Use same function but on all R not
            // "H>0" like in the previous loop &
            // i>0!

            R0c += tmpc;

            if (tmpc > Rkmaxc) {
              Rkmaxc = tmpc;
            }
          }
        }
      }

      Rk[i - 1] = R0c;
      Rp[i - 1] = Rkmaxc;
    }

    if (i >= 9) {
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            R0ct += 2 * getH((double)Wc(j), (double)Wc(k),
                             (double)R(j, k));  // Use same function but on all
            // R not "H>0" like in the
            // previous loop & i>0!
          }
        }
      }
    }
  }

  // can be column vs row selected that can explain the issue!
  double HATSTc = HATSk[0] + HATSct;
  for (int ii = 1; ii < 9; ii++) {
    HATSTc += 2 * HATSk[ii];
  }

  double HTc = Hk[0] + H0ct;
  for (int ii = 1; ii < 9; ii++) {
    HTc += 2.0 * Hk[ii];
  }

  double RTc = 0.0;
  for (double rk_ii : Rk) {
    RTc += 2.0 * rk_ii;
  }

  double RTMc = getMax(Rp);

  // create the output vector...
  for (int i = 0; i < 9; i++) {
    res[i + 4] = roundn(Hk[i], 3);
    res[i + 14] = roundn(HATSk[i], 3);
  }

  res[13] = roundn(HTc, 3);
  res[23] = roundn(HATSTc, 3);

  res[24] = roundn(rcon, 3);
  res[25] = roundn(RARS, 3);
  res[26] = roundn(EIG(0), 3);

  for (int i = 0; i < 8; i++) {
    res[i + 27] = roundn(Rk[i], 3);
    res[i + 36] = roundn(Rp[i], 3);
  }

  res[35] = roundn(RTc + R0ct, 3);
  res[44] = roundn(RTMc, 3);
}

void getGETAWAYDesc(MatrixXd H, MatrixXd R, MatrixXd Adj, int numAtoms,
                    std::vector<int> Heavylist, const ROMol& mol,
                    std::vector<double>& res, unsigned int precision) {
  // prepare data for Getaway parameter computation
  // compute parameters

  VectorXd Lev = H.diagonal();
  std::vector<double> heavyLev;

  for (int i = 0; i < numAtoms; i++) {
    if (Heavylist[i] == 1) {
      heavyLev.push_back(round_to_n_digits_(Lev(i), precision));
    }
  }

  std::vector<double> Clus = clusterArray2(heavyLev, precision);

  double numHeavy = heavyLev.size();
  double ITH0 = numHeavy * log(numHeavy) / log(2.0);
  double ITH = ITH0;
  for (double Clu : Clus) {
    ITH -= Clu * log(Clu) / log(2.0);
  }
  res[0] = roundn(ITH, 3);  // issue sometime with this due to cluster
  double ISH = ITH / ITH0;
  res[1] = roundn(ISH, 3);  // issue sometime with this due to cluster

  // use the PBF to determine 2D vs 3D (with Threshold)
  // determine if this is a plane molecule
  double pbf = RDKit::Descriptors::PBF(mol);
  double D;
  if (pbf < 1.e-5) {
    D = 2.0;
  } else {
    D = 3.0;
  }

  // what about linear molecule ie (D=1) ?
  double HIC = 0.0;
  for (int i = 0; i < numAtoms; i++) {
    HIC -= H(i, i) / D * log(H(i, i) / D) / log(2.);
  }
  res[2] = roundn(HIC, 3);

  double HGM = 1.0;
  for (int i = 0; i < numAtoms; i++) {
    HGM = HGM * H(i, i);
  }
  HGM = 100.0 * pow(HGM, 1.0 / numAtoms);
  res[3] = roundn(HGM, 3);

  double RARS = R.rowwise().sum().sum() / numAtoms;

  JacobiSVD<MatrixXd> mysvd = getSVD(R);

  VectorXd EIG = mysvd.singularValues();

  double rcon = getRCON(R, std::move(Adj), numAtoms);

  std::vector<double> wp = moldata3D.GetRelativePol(mol);

  VectorXd Wp = getEigenVect(wp);

  std::vector<double> wm = moldata3D.GetRelativeMW(mol);

  VectorXd Wm = getEigenVect(wm);

  std::vector<double> wi = moldata3D.GetRelativeIonPol(mol);

  VectorXd Wi = getEigenVect(wi);

  std::vector<double> wv = moldata3D.GetRelativeVdW(mol);

  VectorXd Wv = getEigenVect(wv);

  std::vector<double> we = moldata3D.GetRelativeENeg(mol);

  VectorXd We = getEigenVect(we);

  std::vector<double> wu = moldata3D.GetUn(numAtoms);

  VectorXd Wu = getEigenVect(wu);

  std::vector<double> ws = moldata3D.GetIState(mol);

  VectorXd Ws = getEigenVect(ws);

  MatrixXd Bi;
  MatrixXd RBw;
  double HATSu, HATSm, HATSv, HATSe, HATSp, HATSi, HATSs;
  double HATSut, HATSmt, HATSvt, HATSet, HATSpt, HATSit, HATSst;
  double H0ut, H0mt, H0vt, H0et, H0pt, H0it, H0st;
  double R0ut, R0mt, R0vt, R0et, R0pt, R0it, R0st;
  double H0u, H0m, H0v, H0e, H0p, H0i, H0s;
  double R0u, R0m, R0v, R0e, R0p, R0i, R0s;
  double Rkmaxu, Rkmaxm, Rkmaxv, Rkmaxe, Rkmaxp, Rkmaxi, Rkmaxs;
  double tmpu, tmpm, tmpv, tmpe, tmpp, tmpi, tmps;
  double HATSk[7][9];
  double Hk[7][9];
  double Rk[7][8];
  double Rp[7][8];

  double* dist =
      MolOps::getDistanceMat(mol, false);  // need to be be set to false to have
                                           // topological distance not weighted!

  Map<MatrixXd> D2(dist, numAtoms, numAtoms);

  double Dmax = D2.colwise().maxCoeff().maxCoeff();

  HATSut = 0.0;
  HATSmt = 0.0;
  HATSvt = 0.0;
  HATSet = 0.0;
  HATSpt = 0.0;
  HATSit = 0.0;
  HATSst = 0.0;

  H0ut = 0.0;
  H0mt = 0.0;
  H0vt = 0.0;
  H0et = 0.0;
  H0pt = 0.0;
  H0it = 0.0;
  H0st = 0.0;

  R0ut = 0.0;
  R0mt = 0.0;
  R0vt = 0.0;
  R0et = 0.0;
  R0pt = 0.0;
  R0it = 0.0;
  R0st = 0.0;

  // need to loop other all the D values so <= not <!
  for (int i = 0; i <= Dmax; i++) {
    if (i == 0) {
      Bi = H.diagonal().asDiagonal();
    }

    HATSu = 0.0;
    HATSm = 0.0;
    HATSv = 0.0;
    HATSe = 0.0;
    HATSp = 0.0;
    HATSi = 0.0;
    HATSs = 0.0;

    H0u = 0.0;
    H0m = 0.0;
    H0v = 0.0;
    H0e = 0.0;
    H0p = 0.0;
    H0i = 0.0;
    H0s = 0.0;

    if (i == 0) {  // use Bi
      for (int j = 0; j < numAtoms; j++) {
        for (int k = j; k < numAtoms; k++) {
          if (Bi(j, k) > 0) {
            HATSu += getHATS((double)Wu(j), (double)Wu(j), (double)H(j, j),
                             (double)H(j, j));
            HATSm += getHATS((double)Wm(j), (double)Wm(j), (double)H(j, j),
                             (double)H(j, j));
            HATSv += getHATS((double)Wv(j), (double)Wv(j), (double)H(j, j),
                             (double)H(j, j));
            HATSe += getHATS((double)We(j), (double)We(j), (double)H(j, j),
                             (double)H(j, j));
            HATSp += getHATS((double)Wp(j), (double)Wp(j), (double)H(j, j),
                             (double)H(j, j));
            HATSi += getHATS((double)Wi(j), (double)Wi(j), (double)H(j, j),
                             (double)H(j, j));
            HATSs += getHATS((double)Ws(j), (double)Ws(j), (double)H(j, j),
                             (double)H(j, j));

            if (H(j, k) > 0) {
              H0u += getH((double)Wu(j), (double)Wu(k), (double)H(j, k));
              H0m += getH((double)Wm(j), (double)Wm(k), (double)H(j, k));
              H0v += getH((double)Wv(j), (double)Wv(k), (double)H(j, k));
              H0e += getH((double)We(j), (double)We(k), (double)H(j, k));
              H0p += getH((double)Wp(j), (double)Wp(k), (double)H(j, k));
              H0i += getH((double)Wi(j), (double)Wi(k), (double)H(j, k));
              H0s += getH((double)Ws(j), (double)Ws(k), (double)H(j, k));
            }
          }
        }
      }
    }
    auto Bimat = GetGeodesicMatrix(dist, i, numAtoms);
    Map<MatrixXd> Bj(Bimat.data(), numAtoms, numAtoms);
    if (i > 0 && i < 9) {  // use Bj
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            HATSu += getHATS((double)Wu(j), (double)Wu(k), (double)H(j, j),
                             (double)H(k, k));
            HATSm += getHATS((double)Wm(j), (double)Wm(k), (double)H(j, j),
                             (double)H(k, k));
            HATSv += getHATS((double)Wv(j), (double)Wv(k), (double)H(j, j),
                             (double)H(k, k));
            HATSe += getHATS((double)We(j), (double)We(k), (double)H(j, j),
                             (double)H(k, k));
            HATSp += getHATS((double)Wp(j), (double)Wp(k), (double)H(j, j),
                             (double)H(k, k));
            HATSi += getHATS((double)Wi(j), (double)Wi(k), (double)H(j, j),
                             (double)H(k, k));
            HATSs += getHATS((double)Ws(j), (double)Ws(k), (double)H(j, j),
                             (double)H(k, k));

            if (H(j, k) > 0) {
              H0u += getH((double)Wu(j), (double)Wu(k), (double)H(j, k));
              H0m += getH((double)Wm(j), (double)Wm(k), (double)H(j, k));
              H0v += getH((double)Wv(j), (double)Wv(k), (double)H(j, k));
              H0e += getH((double)We(j), (double)We(k), (double)H(j, k));
              H0p += getH((double)Wp(j), (double)Wp(k), (double)H(j, k));
              H0i += getH((double)Wi(j), (double)Wi(k), (double)H(j, k));
              H0s += getH((double)Ws(j), (double)Ws(k), (double)H(j, k));
            }
          }
        }
      }
    }

    if (i >= 9) {  // Totals missing part
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            HATSut += 2 * getHATS((double)Wu(j), (double)Wu(k), (double)H(j, j),
                                  (double)H(k, k));
            HATSmt += 2 * getHATS((double)Wm(j), (double)Wm(k), (double)H(j, j),
                                  (double)H(k, k));
            HATSvt += 2 * getHATS((double)Wv(j), (double)Wv(k), (double)H(j, j),
                                  (double)H(k, k));
            HATSet += 2 * getHATS((double)We(j), (double)We(k), (double)H(j, j),
                                  (double)H(k, k));
            HATSpt += 2 * getHATS((double)Wp(j), (double)Wp(k), (double)H(j, j),
                                  (double)H(k, k));
            HATSit += 2 * getHATS((double)Wi(j), (double)Wi(k), (double)H(j, j),
                                  (double)H(k, k));
            HATSst += 2 * getHATS((double)Ws(j), (double)Ws(k), (double)H(j, j),
                                  (double)H(k, k));

            if (H(j, k) > 0) {
              H0ut += 2 * getH((double)Wu(j), (double)Wu(k), (double)H(j, k));
              H0mt += 2 * getH((double)Wm(j), (double)Wm(k), (double)H(j, k));
              H0vt += 2 * getH((double)Wv(j), (double)Wv(k), (double)H(j, k));
              H0et += 2 * getH((double)We(j), (double)We(k), (double)H(j, k));
              H0pt += 2 * getH((double)Wp(j), (double)Wp(k), (double)H(j, k));
              H0it += 2 * getH((double)Wi(j), (double)Wi(k), (double)H(j, k));
              H0st += 2 * getH((double)Ws(j), (double)Ws(k), (double)H(j, k));
            }
          }
        }
      }
    }

    if (i < 9) {
      HATSk[0][i] = HATSu;
      HATSk[1][i] = HATSm;
      HATSk[2][i] = HATSv;
      HATSk[3][i] = HATSe;
      HATSk[4][i] = HATSp;
      HATSk[5][i] = HATSi;
      HATSk[6][i] = HATSs;

      Hk[0][i] = H0u;
      Hk[1][i] = H0m;
      Hk[2][i] = H0v;
      Hk[3][i] = H0e;
      Hk[4][i] = H0p;
      Hk[5][i] = H0i;
      Hk[6][i] = H0s;
    }

    R0u = 0.0;
    R0m = 0.0;
    R0v = 0.0;
    R0e = 0.0;
    R0p = 0.0;
    R0i = 0.0;
    R0s = 0.0;

    Rkmaxu = 0;
    Rkmaxm = 0;
    Rkmaxv = 0;
    Rkmaxe = 0;
    Rkmaxp = 0;
    Rkmaxi = 0;
    Rkmaxs = 0;

    // from 1 to 8
    if (i > 0 && i < 9) {
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            tmpu = getH((double)Wu(j), (double)Wu(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            tmpm = getH((double)Wm(j), (double)Wm(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            tmpv = getH((double)Wv(j), (double)Wv(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            tmpe = getH((double)We(j), (double)We(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            tmpp = getH((double)Wp(j), (double)Wp(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            tmpi = getH((double)Wi(j), (double)Wi(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            tmps = getH((double)Ws(j), (double)Ws(k),
                        (double)R(j, k));  // Use same function but on all R not
                                           // "H>0" like in the previous loop &
                                           // i>0!
            R0u += tmpu;
            R0m += tmpm;
            R0v += tmpv;
            R0e += tmpe;
            R0p += tmpp;
            R0i += tmpi;
            R0s += tmps;
            if (tmpu > Rkmaxu) {
              Rkmaxu = tmpu;
            }
            if (tmpm > Rkmaxm) {
              Rkmaxm = tmpm;
            }
            if (tmpv > Rkmaxv) {
              Rkmaxv = tmpv;
            }
            if (tmpe > Rkmaxe) {
              Rkmaxe = tmpe;
            }
            if (tmpp > Rkmaxp) {
              Rkmaxp = tmpp;
            }
            if (tmpi > Rkmaxi) {
              Rkmaxi = tmpi;
            }
            if (tmps > Rkmaxs) {
              Rkmaxs = tmps;
            }
          }
        }
      }

      Rk[0][i - 1] = R0u;
      Rk[1][i - 1] = R0m;
      Rk[2][i - 1] = R0v;
      Rk[3][i - 1] = R0e;
      Rk[4][i - 1] = R0p;
      Rk[5][i - 1] = R0i;
      Rk[6][i - 1] = R0s;

      Rp[0][i - 1] = Rkmaxu;
      Rp[1][i - 1] = Rkmaxm;
      Rp[2][i - 1] = Rkmaxv;
      Rp[3][i - 1] = Rkmaxe;
      Rp[4][i - 1] = Rkmaxp;
      Rp[5][i - 1] = Rkmaxi;
      Rp[6][i - 1] = Rkmaxs;
    }

    if (i >= 9) {
      for (int j = 0; j < numAtoms - 1; j++) {
        for (int k = j + 1; k < numAtoms; k++) {
          if (Bj(j, k) == 1) {
            R0ut += 2 * getH((double)Wu(j), (double)Wu(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
            R0mt += 2 * getH((double)Wm(j), (double)Wm(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
            R0vt += 2 * getH((double)Wv(j), (double)Wv(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
            R0et += 2 * getH((double)We(j), (double)We(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
            R0pt += 2 * getH((double)Wp(j), (double)Wp(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
            R0it += 2 * getH((double)Wi(j), (double)Wi(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
            R0st += 2 * getH((double)Ws(j), (double)Ws(k),
                             (double)R(j, k));  // Use same function but on all
                                                // R not "H>0" like in the
                                                // previous loop & i>0!
          }
        }
      }
    }
  }

  // can be column vs row selected that can explain the issue!
  double HATSTu = HATSk[0][0] + HATSut;
  double HATSTm = HATSk[1][0] + HATSmt;
  double HATSTv = HATSk[2][0] + HATSvt;
  double HATSTe = HATSk[3][0] + HATSet;
  double HATSTp = HATSk[4][0] + HATSpt;
  double HATSTi = HATSk[5][0] + HATSit;
  double HATSTs = HATSk[6][0] + HATSst;

  for (int ii = 1; ii < 9; ii++) {
    HATSTu += 2 * HATSk[0][ii];
    HATSTm += 2 * HATSk[1][ii];
    HATSTv += 2 * HATSk[2][ii];
    HATSTe += 2 * HATSk[3][ii];
    HATSTp += 2 * HATSk[4][ii];
    HATSTi += 2 * HATSk[5][ii];
    HATSTs += 2 * HATSk[6][ii];
  }

  double HTu = Hk[0][0] + H0ut;
  double HTm = Hk[1][0] + H0mt;
  double HTv = Hk[2][0] + H0vt;
  double HTe = Hk[3][0] + H0et;
  double HTp = Hk[4][0] + H0pt;
  double HTi = Hk[5][0] + H0it;
  double HTs = Hk[6][0] + H0st;

  for (int ii = 1; ii < 9; ii++) {
    HTu += 2.0 * Hk[0][ii];
    HTm += 2.0 * Hk[1][ii];
    HTv += 2.0 * Hk[2][ii];
    HTe += 2.0 * Hk[3][ii];
    HTp += 2.0 * Hk[4][ii];
    HTi += 2.0 * Hk[5][ii];
    HTs += 2.0 * Hk[6][ii];
  }

  // 2*(Rk[1]+Rk[2]+Rk[3]+Rk[4]+Rk[5]+Rk[6]+Rk[7]+Rk[8]) // this is not true
  // this is on all the element

  double RTu = 0.0;
  double RTm = 0.0;
  double RTv = 0.0;
  double RTe = 0.0;
  double RTp = 0.0;
  double RTi = 0.0;
  double RTs = 0.0;

  for (int ii = 0; ii < 8; ii++) {
    RTu += 2.0 * Rk[0][ii];
    RTm += 2.0 * Rk[1][ii];
    RTv += 2.0 * Rk[2][ii];
    RTe += 2.0 * Rk[3][ii];
    RTp += 2.0 * Rk[4][ii];
    RTi += 2.0 * Rk[5][ii];
    RTs += 2.0 * Rk[6][ii];
  }

  double RTMu = getMax(Rp[0]);
  double RTMm = getMax(Rp[1]);
  double RTMv = getMax(Rp[2]);
  double RTMe = getMax(Rp[3]);
  double RTMp = getMax(Rp[4]);
  double RTMi = getMax(Rp[5]);
  double RTMs = getMax(Rp[6]);

  // create the output vector...
  for (int i = 0; i < 9; i++) {
    res[i + 4] = roundn(Hk[0][i], 3);
    res[i + 14] = roundn(HATSk[0][i], 3);
    res[i + 24] = roundn(Hk[1][i], 3);
    res[i + 34] = roundn(HATSk[1][i], 3);
    res[i + 44] = roundn(Hk[2][i], 3);
    res[i + 54] = roundn(HATSk[2][i], 3);
    res[i + 64] = roundn(Hk[3][i], 3);
    res[i + 74] = roundn(HATSk[3][i], 3);
    res[i + 84] = roundn(Hk[4][i], 3);
    res[i + 94] = roundn(HATSk[4][i], 3);
    res[i + 104] = roundn(Hk[5][i], 3);
    res[i + 114] = roundn(HATSk[5][i], 3);
    res[i + 124] = roundn(Hk[6][i], 3);
    res[i + 134] = roundn(HATSk[6][i], 3);
  }

  res[13] = roundn(HTu, 3);
  res[23] = roundn(HATSTu, 3);
  res[33] = roundn(HTm, 3);
  res[43] = roundn(HATSTm, 3);
  res[53] = roundn(HTv, 3);
  res[63] = roundn(HATSTv, 3);
  res[73] = roundn(HTe, 3);
  res[83] = roundn(HATSTe, 3);
  res[93] = roundn(HTp, 3);
  res[103] = roundn(HATSTp, 3);
  res[113] = roundn(HTi, 3);
  res[123] = roundn(HATSTi, 3);
  res[133] = roundn(HTs, 3);
  res[143] = roundn(HATSTs, 3);

  res[144] = roundn(rcon, 3);
  res[145] = roundn(RARS, 3);
  res[146] = roundn(EIG(0), 3);

  for (int i = 0; i < 8; i++) {
    res[i + 147] = roundn(Rk[0][i], 3);
    res[i + 156] = roundn(Rp[0][i], 3);
    res[i + 165] = roundn(Rk[1][i], 3);
    res[i + 174] = roundn(Rp[1][i], 3);
    res[i + 183] = roundn(Rk[2][i], 3);
    res[i + 192] = roundn(Rp[2][i], 3);
    res[i + 201] = roundn(Rk[3][i], 3);
    res[i + 210] = roundn(Rp[3][i], 3);
    res[i + 219] = roundn(Rk[4][i], 3);
    res[i + 228] = roundn(Rp[4][i], 3);
    res[i + 237] = roundn(Rk[5][i], 3);
    res[i + 246] = roundn(Rp[5][i], 3);
    res[i + 255] = roundn(Rk[6][i], 3);
    res[i + 264] = roundn(Rp[6][i], 3);
  }

  res[155] = roundn(RTu + R0ut, 3);
  res[164] = roundn(RTMu, 3);
  res[173] = roundn(RTm + R0mt, 3);
  res[182] = roundn(RTMm, 3);
  res[191] = roundn(RTv + R0vt, 3);
  res[200] = roundn(RTMv, 3);
  res[209] = roundn(RTe + R0et, 3);
  res[218] = roundn(RTMe, 3);
  res[227] = roundn(RTp + R0pt, 3);
  res[236] = roundn(RTMp, 3);
  res[245] = roundn(RTi + R0it, 3);
  res[254] = roundn(RTMi, 3);
  res[263] = roundn(RTs + R0st, 3);
  res[272] = roundn(RTMs, 3);
}

/*
std::vector<std::string>
GETAWAYNAMES={"ITH","ISH","HIC","HGM","H0u","H1u","H2u","H3u","H4u","H5u","H6u","H7u","H8u","HTu",
"HATS0u","HATS1u","HATS2u","HATS3u","HATS4u","HATS5u","HATS6u","HATS7u","HATS8u","HATSu","H0m","H1m","H2m","H3m","H4m","H5m",
"H6m","H7m","H8m","HTm","HATS0m","HATS1m","HATS2m","HATS3m","HATS4m","HATS5m","HATS6m","HATS7m","HATS8m","HATSm","H0v","H1v",
"H2v","H3v","H4v","H5v","H6v","H7v","H8v","HTv","HATS0v","HATS1v","HATS2v","HATS3v","HATS4v","HATS5v","HATS6v","HATS7v","HATS8v",
"HATSv","H0e","H1e","H2e","H3e","H4e","H5e","H6e","H7e","H8e","HTe","HATS0e","HATS1e","HATS2e","HATS3e","HATS4e","HATS5e","HATS6e",
"HATS7e","HATS8e","HATSe","H0p","H1p","H2p","H3p","H4p","H5p","H6p","H7p","H8p","HTp","HATS0p","HATS1p","HATS2p","HATS3p","HATS4p",
"HATS5p","HATS6p","HATS7p","HATS8p","HATSp","H0i","H1i","H2i","H3i","H4i","H5i","H6i","H7i","H8i","HTi","HATS0i","HATS1i","HATS2i",
"HATS3i","HATS4i","HATS5i","HATS6i","HATS7i","HATS8i","HATSi","H0s","H1s","H2s","H3s","H4s","H5s","H6s","H7s","H8s","HTs","HATS0s",
"HATS1s","HATS2s","HATS3s","HATS4s","HATS5s","HATS6s","HATS7s","HATS8s","HATSs","RCON","RARS","REIG","R1u","R2u","R3u","R4u","R5u",
"R6u","R7u","R8u","RTu","R1u+","R2u+","R3u+","R4u+","R5u+","R6u+","R7u+","R8u+","RTu+","R1m","R2m","R3m","R4m","R5m","R6m","R7m",
"R8m","RTm","R1m+","R2m+","R3m+","R4m+","R5m+","R6m+","R7m+","R8m+","RTm+","R1v","R2v","R3v","R4v","R5v","R6v","R7v","R8v","RTv",
"R1v+","R2v+","R3v+","R4v+","R5v+","R6v+","R7v+","R8v+","RTv+","R1e","R2e","R3e","R4e","R5e","R6e","R7e","R8e","RTe","R1e+","R2e+",
"R3e+","R4e+","R5e+","R6e+","R7e+","R8e+","RTe+","R1p","R2p","R3p","R4p","R5p","R6p","R7p","R8p","RTp","R1p+","R2p+","R3p+","R4p+",
"R5p+","R6p+","R7p+","R8p+","RTp+","R1i","R2i","R3i","R4i","R5i","R6i","R7i","R8i","RTi","R1i+","R2i+","R3i+","R4i+","R5i+","R6i+",
"R7i+","R8i+","RTi+","R1s","R2s","R3s","R4s","R5s","R6s","R7s","R8s","RTs","R1s+","R2s+","R3s+","R4s+","R5s+","R6s+","R7s+","R8s+","RTs+"};
 */

void GetGETAWAYone(double* dist3D, double* AdjMat, std::vector<double> Vpoints,
                   const ROMol& mol, const Conformer& conf,
                   std::vector<int> Heavylist, std::vector<double>& res,
                   unsigned int precision,
                   const std::string& customAtomPropName) {
  PRECONDITION(dist3D != nullptr, "no distance matrix");
  PRECONDITION(AdjMat != nullptr, "no adjacency matrix");
  int numAtoms = conf.getNumAtoms();

  Map<MatrixXd> ADJ(AdjMat, numAtoms, numAtoms);

  Map<MatrixXd> DM(dist3D, numAtoms, numAtoms);

  double* vpoints = &Vpoints[0];

  Map<MatrixXd> matorigin(vpoints, 3, numAtoms);

  MatrixXd MatOrigin = matorigin.transpose();

  MatrixXd Xmean = GetCenterMatrix(MatOrigin);

  MatrixXd H = GetHmatrix(Xmean);

  MatrixXd R = GetRmatrix(H, DM, numAtoms);

  getGETAWAYDescCustom(H, R, ADJ, numAtoms, std::move(Heavylist), mol, res,
                       precision, customAtomPropName);
}

void GetGETAWAY(double* dist3D, double* AdjMat, std::vector<double> Vpoints,
                const ROMol& mol, const Conformer& conf,
                std::vector<int> Heavylist, std::vector<double>& res,
                unsigned int precision) {
  PRECONDITION(dist3D != nullptr, "no distance matrix");
  PRECONDITION(AdjMat != nullptr, "no adjacency matrix");

  int numAtoms = conf.getNumAtoms();

  Map<MatrixXd> ADJ(AdjMat, numAtoms, numAtoms);

  Map<MatrixXd> DM(dist3D, numAtoms, numAtoms);

  double* vpoints = &Vpoints[0];

  Map<MatrixXd> matorigin(vpoints, 3, numAtoms);

  MatrixXd MatOrigin = matorigin.transpose();

  MatrixXd Xmean = GetCenterMatrix(MatOrigin);

  MatrixXd H = GetHmatrix(Xmean);

  MatrixXd R = GetRmatrix(H, DM, numAtoms);

  getGETAWAYDesc(H, R, ADJ, numAtoms, std::move(Heavylist), mol, res,
                 precision);
}

}  // end of anonymous namespace

void GETAWAY(const ROMol& mol, std::vector<double>& res, int confId,
             unsigned int precision, const std::string& customAtomPropName) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")

  int numAtoms = mol.getNumAtoms();

  const Conformer& conf = mol.getConformer(confId);

  std::vector<double> Vpoints(3 * numAtoms);

  for (int i = 0; i < numAtoms; ++i) {
    Vpoints[3 * i] = conf.getAtomPos(i).x;
    Vpoints[3 * i + 1] = conf.getAtomPos(i).y;
    Vpoints[3 * i + 2] = conf.getAtomPos(i).z;
  }

  std::vector<int> Heavylist = GetHeavyList(mol);
  // should be the same as the List size upper!
  // int nHeavyAt= mol.getNumHeavyAtoms();

  double* dist3D = MolOps::get3DDistanceMat(mol, confId);

  double* AdjMat = MolOps::getAdjacencyMatrix(
      mol, false, 0, false,
      nullptr);  // false to have only the 1,0 matrix unweighted

  if (!customAtomPropName.empty()) {
    res.clear();
    res.resize(45);
    GetGETAWAYone(dist3D, AdjMat, Vpoints, mol, conf, Heavylist, res, precision,
                  customAtomPropName);
  } else {
    res.clear();
    res.resize(273);

    GetGETAWAY(dist3D, AdjMat, Vpoints, mol, conf, Heavylist, res, precision);
  }
}

}  // namespace Descriptors
}  // namespace RDKit
