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
// Guillaume GODIN access the AutoCorrelation 2D descriptors names in Dragon TDB

#include <GraphMol/RDKitBase.h>

#include "AUTOCORR2D.h"
#include "MolData3Ddescriptors.h"
#include <cmath>
#include <iostream>

namespace RDKit {
namespace Descriptors {

namespace {

MolData3Ddescriptors moldata3D;

// this is the Broto-Moreau 2D descriptors (centered or not)
void get2DautocorrelationDesc(const double* dist, unsigned int numAtoms,
                              const ROMol& mol, std::vector<double>& res) {
  std::vector<double> wp = moldata3D.GetRelativePol(mol);
  std::vector<double> wm = moldata3D.GetRelativeMW(mol);
  std::vector<double> wv = moldata3D.GetRelativeVdW(mol);
  std::vector<double> wi = moldata3D.GetRelativeIonPol(mol);
  std::vector<double> we = moldata3D.GetRelativeENeg(mol);
  std::vector<double> ws = moldata3D.GetIState(mol);
  std::vector<double> w(6 * numAtoms, 0.0);
  std::vector<double> wmean(6, 0.0);

  for (unsigned int i = 0; i < numAtoms; i++) {
    w[0 * numAtoms + i] = wm[i];
    w[1 * numAtoms + i] = wv[i];
    w[2 * numAtoms + i] = we[i];
    w[3 * numAtoms + i] = wp[i];
    w[4 * numAtoms + i] = wi[i];
    w[5 * numAtoms + i] = ws[i];
  }

  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int t = 0; t < 6; ++t) {
      wmean[t] += w[t * numAtoms + i] / (double)numAtoms;
    }
  }

  std::vector<double> squaresumdiff(6, 0.0);
  for (unsigned int i = 0; i < numAtoms; i++) {
    for (unsigned int t = 0; t < 6; ++t) {
      squaresumdiff[t] +=
          (w[t * numAtoms + i] - wmean[t]) * (w[t * numAtoms + i] - wmean[t]);
    }
  }

  std::vector<double> TDBmat(48, 0.0);
  std::vector<double> TDBmatC(48, 0.0);
  std::vector<double> TDBmatM(48, 0.0);
  std::vector<double> TDBmatG(48, 0.0);

  for (unsigned int k = 0; k < 8; k++) {
    int maxkVertexPairs = 0;
    for (unsigned int i = 0; i < numAtoms; ++i) {
      for (unsigned int j = i + 1; j < numAtoms; ++j) {
        if (dist[j * numAtoms + i] == k + 1) {
          for (unsigned int t = 0; t < 6; ++t) {
            TDBmatM[t * 8 + k] += (w[t * numAtoms + i] - wmean[t]) *
                                  (w[t * numAtoms + j] - wmean[t]);  // ATSC

            TDBmatG[t * 8 + k] += (w[t * numAtoms + i] - w[t * numAtoms + j]) *
                                  (w[t * numAtoms + i] - w[t * numAtoms + j]);

            TDBmat[t * 8 + k] += w[t * numAtoms + i] * w[t * numAtoms + j];
            TDBmatC[t * 8 + k] += fabs(w[t * numAtoms + i] - wmean[t]) *
                                  fabs(w[t * numAtoms + j] - wmean[t]);
          }
          maxkVertexPairs += 1;
        }
      }
    }

    for (unsigned int t = 0; t < 6; ++t) {
      if (maxkVertexPairs > 0 && squaresumdiff[t] > 0) {
        TDBmat[t * 8 + k] = log1p(TDBmat[t * 8 + k]);
        TDBmatG[t * 8 + k] = TDBmatG[t * 8 + k] / squaresumdiff[t] /
                             maxkVertexPairs * (numAtoms - 1) / 2.0;
        TDBmatM[t * 8 + k] =
            TDBmatM[t * 8 + k] / squaresumdiff[t] / maxkVertexPairs * numAtoms;
      } else {
        TDBmat[t * 8 + k] = 0.0;
        TDBmatC[t * 8 + k] = 0.0;
        TDBmatM[t * 8 + k] = 0.0;
        TDBmatG[t * 8 + k] = 0.0;
      }
    }
  }

  // update the Output vector!
  for (unsigned int t = 0; t < 6; ++t) {
    for (unsigned int k = 0; k < 8; ++k) {
      res[t * 8 + k] = std::round(1000 * TDBmat[k + t * 8]) / 1000;
      res[t * 8 + k + 48] = std::round(1000 * TDBmatC[k + t * 8]) / 1000;
      res[t * 8 + k + 96] = std::round(1000 * TDBmatM[k + t * 8]) / 1000;
      res[t * 8 + k + 144] = std::round(1000 * TDBmatG[k + t * 8]) / 1000;
    }
  }

  TDBmat.clear();
  TDBmatC.clear();
  TDBmatM.clear();
  TDBmatG.clear();

  wp.clear();
  wm.clear();
  wv.clear();
  wi.clear();
  we.clear();
  ws.clear();
  w.clear();
  wmean.clear();
}

// this is the Broto-Moreau 2D descriptors (centered or not)
void get2DautocorrelationDescCustom(const double* dist, unsigned int numAtoms,
                                    const ROMol& mol, std::vector<double>& res,
                                    const std::string& customAtomPropName) {
  std::vector<double> wc = moldata3D.GetCustomAtomProp(mol, customAtomPropName);
  std::vector<double> w(numAtoms, 0.0);
  double wmean = 0.0;

  for (unsigned int i = 0; i < numAtoms; i++) {
    wmean += wc[+i] / (double)numAtoms;
  }

  double squaresumdiff = 0.0;
  for (unsigned int i = 0; i < numAtoms; i++) {
    squaresumdiff += (wc[i] - wmean) * (wc[i] - wmean);
  }

  std::vector<double> TDBmat(8, 0.0);
  std::vector<double> TDBmatC(8, 0.0);
  std::vector<double> TDBmatM(8, 0.0);
  std::vector<double> TDBmatG(8, 0.0);

  for (unsigned int k = 0; k < 8; k++) {
    int maxkVertexPairs = 0;
    for (unsigned int i = 0; i < numAtoms; ++i) {
      for (unsigned int j = i + 1; j < numAtoms; ++j) {
        if (dist[j * numAtoms + i] == k + 1) {
          TDBmatM[k] += (wc[i] - wmean) * (wc[j] - wmean);
          TDBmatG[k] += (wc[i] - wc[j]) * (wc[i] - w[j]);
          TDBmat[k] += wc[i] * wc[j];
          TDBmatC[k] += fabs(wc[i] - wmean) * fabs(wc[j] - wmean);
          maxkVertexPairs += 1;
        }
      }
    }

    for (unsigned int t = 0; t < 1; ++t) {
      if (maxkVertexPairs > 0) {
        TDBmat[k] = log1p(TDBmat[k]);
        TDBmatG[k] =
            TDBmatG[k] / squaresumdiff / maxkVertexPairs * (numAtoms - 1) / 2.0;
        TDBmatM[k] = TDBmatM[k] / squaresumdiff / maxkVertexPairs * numAtoms;

      } else {
        TDBmat[k] = 0.0;
        TDBmatC[k] = 0.0;
        TDBmatM[k] = 0.0;
        TDBmatG[k] = 0.0;
      }
    }
  }

  // update the Output vector!
  for (unsigned int k = 0; k < 8; ++k) {
    res[k] = std::round(1000 * TDBmat[k]) / 1000;
    res[k + 8] = std::round(1000 * TDBmatC[k]) / 1000;
    res[k + 16] = std::round(1000 * TDBmatM[k]) / 1000;
    res[k + 24] = std::round(1000 * TDBmatG[k]) / 1000;
  }

  TDBmat.clear();
  TDBmatC.clear();
  TDBmatM.clear();
  TDBmatG.clear();

  wc.clear();
}

void Get2Dauto(const double* dist, unsigned int numAtoms, const ROMol& mol,
               std::vector<double>& res) {
  get2DautocorrelationDesc(dist, numAtoms, mol, res);
}
void Get2Dautoone(const double* dist, unsigned int numAtoms, const ROMol& mol,
                  std::vector<double>& res,
                  const std::string& customAtomPropName) {
  get2DautocorrelationDescCustom(dist, numAtoms, mol, res, customAtomPropName);
}
}  // end of anonymous namespace

void AUTOCORR2D(const ROMol& mol, std::vector<double>& result,
                const std::string& customAtomPropName) {
  unsigned int numAtoms = mol.getNumAtoms();
  double* dist = MolOps::getDistanceMat(mol, false);  // topological matrix
  if (!customAtomPropName.empty()) {
    result.clear();
    result.resize(32);
    Get2Dautoone(dist, numAtoms, mol, result, customAtomPropName);
  } else {
    result.clear();
    result.resize(192);
    Get2Dauto(dist, numAtoms, mol, result);
  }
}
}  // namespace Descriptors
}  // namespace RDKit
