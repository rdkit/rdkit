//
//  Copyright (c) 2016, Guillaume GODIN.
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
// Guillaume GODIN access the AutoCorrelation 3D descriptors names in Dragon TDB

#include <GraphMol/RDKitBase.h>

#include "AUTOCORR3D.h"
#include "MolData3Ddescriptors.h"

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>

using namespace Eigen;
namespace RDKit {
namespace Descriptors {

namespace {

MolData3Ddescriptors moldata3D;

VectorXd getEigenVect(std::vector<double> v) {
  double* varray_ptr = &v[0];
  Map<VectorXd> V(varray_ptr, v.size());
  return V;
}

double* GetGeodesicMatrix(double* dist, int lag, int numAtoms) {
  int sizeArray = numAtoms * numAtoms;
  auto* Geodesic = new double[sizeArray];
  for (int i = 0; i < sizeArray; i++) {
    if (dist[i] == lag) {
      Geodesic[i] = 1.0;
    } else {
      Geodesic[i] = 0.0;
    }
  }
  return Geodesic;
}

// matrix prod that mimic the v2 version
// the code is in respect to Dragon 6 descriptors.
// replace the number of "Bicount" vertex per lag by a (numAtoms * (numAtoms -
// 1))!
// provided by Kobe team!
void get3DautocorrelationDesc(double* dist3D, double* topologicaldistance,
                              int numAtoms, const ROMol& mol,
                              std::vector<double>& res) {
  Map<MatrixXd> dm(dist3D, numAtoms, numAtoms);
  Map<MatrixXd> di(topologicaldistance, numAtoms, numAtoms);

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

  std::vector<double> wr = moldata3D.GetRelativeRcov(mol);
  VectorXd Wr = getEigenVect(wr);

  MatrixXd Bi;
  MatrixXd tmp;
  double TDBmat[8][10];
  double dtmp;

  for (int i = 0; i < 10; i++) {
    double* Bimat = GetGeodesicMatrix(topologicaldistance, i + 1, numAtoms);
    Map<MatrixXd> Bi(Bimat, numAtoms, numAtoms);
    MatrixXd RBi = Bi.cwiseProduct(dm);
    // double Bicount = (double)Bi.sum();

    tmp = Wu.transpose() * RBi * Wu;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[0][i] = dtmp;

    tmp = Wm.transpose() * RBi * Wm;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[1][i] = dtmp;

    tmp = Wv.transpose() * RBi * Wv;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[2][i] = dtmp;

    tmp = We.transpose() * RBi * We;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[3][i] = dtmp;

    tmp = Wp.transpose() * RBi * Wp;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[4][i] = dtmp;

    tmp = Wi.transpose() * RBi * Wi;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[5][i] = dtmp;

    tmp = Ws.transpose() * RBi * Ws;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[6][i] = dtmp;

    tmp = Wr.transpose() * RBi * Wr;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[7][i] = dtmp;
    delete[] Bimat;
  }

  // update the Output vector!
  for (unsigned int j = 0; j < 8; ++j) {
    for (unsigned int i = 0; i < 10; ++i) {
      res[j * 10 + i] =
          std::round(1000 * TDBmat[j][i] / (numAtoms * (numAtoms - 1))) / 1000;
    }
  }
}

void get3DautocorrelationDescCustom(double* dist3D, double* topologicaldistance,
                                    int numAtoms, const ROMol& mol,
                                    std::vector<double>& res,
                                    const std::string& customAtomPropName) {
  Map<MatrixXd> dm(dist3D, numAtoms, numAtoms);
  Map<MatrixXd> di(topologicaldistance, numAtoms, numAtoms);

  std::vector<double> customAtomArray =
      moldata3D.GetCustomAtomProp(mol, customAtomPropName);
  VectorXd Wc = getEigenVect(customAtomArray);

  MatrixXd Bi;
  MatrixXd tmp;
  double TDBmat[10];
  double dtmp;

  for (int i = 0; i < 10; i++) {
    double* Bimat = GetGeodesicMatrix(topologicaldistance, i + 1, numAtoms);
    Map<MatrixXd> Bi(Bimat, numAtoms, numAtoms);
    MatrixXd RBi = Bi.cwiseProduct(dm);

    tmp = Wc.transpose() * RBi * Wc;
    dtmp = (double)tmp(0);
    if (std::isnan(dtmp)) {
      dtmp = 0.0;
    }
    TDBmat[i] = dtmp;
    delete[] Bimat;
  }
  // update the Output vector!
  for (unsigned int i = 0; i < 10; ++i) {
    res[i] = std::round(1000 * TDBmat[i] / (numAtoms * (numAtoms - 1))) / 1000;
  }
}

void Get3Dauto(double* dist3D, double* topologicaldistance, int numAtoms,
               const ROMol& mol, std::vector<double>& res) {
  // AUTOCORRNAMES={"TDB01u","TDB02u","TDB03u","TDB04u","TDB05u","TDB06u","TDB07u","TDB08u","TDB09u","TDB10u","TDB01m","TDB02m","TDB03m","TDB04m","TDB05m","TDB06m","TDB07m","TDB08m","TDB09m","TDB10m","TDB01v","TDB02v","TDB03v","TDB04v","TDB05v","TDB06v","TDB07v","TDB08v","TDB09v","TDB10v","TDB01e","TDB02e","TDB03e","TDB04e","TDB05e","TDB06e","TDB07e","TDB08e","TDB09e","TDB10e","TDB01p","TDB02p","TDB03p","TDB04p","TDB05p","TDB06p","TDB07p","TDB08p","TDB09p","TDB10p","TDB01i","TDB02i","TDB03i","TDB04i","TDB05i","TDB06i","TDB07i","TDB08i","TDB09i","TDB10i","TDB01s","TDB02s","TDB03s","TDB04s","TDB05s","TDB06s","TDB07s","TDB08s","TDB09s","TDB10s","TDB01r","TDB02r","TDB03r","TDB04r","TDB05r","TDB06r","TDB07r","TDB08r","TDB09r","TDB10r"};
  get3DautocorrelationDesc(dist3D, topologicaldistance, numAtoms, mol, res);
}

void Get3Dautoone(double* dist3D, double* topologicaldistance, int numAtoms,
                  const ROMol& mol, std::vector<double>& res,
                  const std::string& customAtomPropName) {
  get3DautocorrelationDescCustom(dist3D, topologicaldistance, numAtoms, mol,
                                 res, customAtomPropName);
}

}  // end of anonymous namespace

void AUTOCORR3D(const ROMol& mol, std::vector<double>& res, int confId,
                const std::string& customAtomPropName) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  double* topologicaldistance =
      MolOps::getDistanceMat(mol, false);  // topological matrix
  double* dist3D =
      MolOps::get3DDistanceMat(mol, confId, false, true);  // 3D distance matrix
  if (customAtomPropName != "") {
    res.clear();
    res.resize(10);
    Get3Dautoone(dist3D, topologicaldistance, numAtoms, mol, res,
                 customAtomPropName);
  } else {
    res.clear();
    res.resize(80);
    Get3Dauto(dist3D, topologicaldistance, numAtoms, mol, res);
  }
}
}  // namespace Descriptors
}  // namespace RDKit
