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
#include <math.h>
#include <iostream>

namespace RDKit {
namespace Descriptors {

namespace {

MolData3Ddescriptors moldata3D;

// this is the Broto-Moreau 2D descriptors (centered or not)
void get2DautocorrelationDesc(double* dist, int numAtoms, const ROMol& mol,
                              std::vector<double>& res) {
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

  /*
                              for (unsigned int i = 0; i < numAtoms; i++) {
                                for (unsigned int t = 0; t < 6; ++t) {
                                  wmean[t] += w[t][i] / (double) numAtoms;
                                }
                              }

                              std::vector<double> squaresumdiff(6,0.0);
                              for (unsigned int i = 0; i < numAtoms; i++) {
                                for (unsigned int t = 0; t < 6; ++t) {
                                  squaresumdiff[t]+=(w[t][i]-wmean[t])*(w[t][i]-wmean[t]);
                                }
                              }
  */

  std::vector<double> TDBmat(48, 0.0);
  // double TDBmatM[6][8] = { {0.0} };
  // double TDBmatG[6][8] = { {0.0} };
  // double TDBmatC[6][8] = { {0.0} };

  for (unsigned int k = 0; k < 8; k++) {
    int maxkVertexPairs = 0;
    for (unsigned int i = 0; i < numAtoms - 1; ++i) {
      for (unsigned int j = i + 1; j < numAtoms; ++j) {
        if (dist[j * numAtoms + i] == k + 1) {
          for (unsigned int t = 0; t < 6; ++t) {
            //  TDBmatM[t][k] += (w[t][i]-wmean[t]) * (w[t][j]-wmean[t]);
            //  TDBmatG[t][k] += (w[t][i] - w[t][j]) * (w[t][i] - w[t][j]);
            TDBmat[t * 8 + k] += w[t * numAtoms + i] * w[t * numAtoms + j];
            //  TDBmatC[t][k] += abs(w[t][i]-wmean[t]) * abs(w[t][j]-wmean[t]);
          }
          maxkVertexPairs += 1;
        }
      }
    }

    for (unsigned int t = 0; t < 6; ++t) {
      if (maxkVertexPairs > 0) {
        TDBmat[t * 8 + k] = log(TDBmat[t * 8 + k] + 1);
        // TDBmatG[t][k] =
        // TDBmatG[t][k]/squaresumdiff[t]/maxkVertexPairs/2/(numAtoms-1);
        // TDBmatM[t][k] =
        // TDBmatM[t][k]/squaresumdiff[t]/maxkVertexPairs/(numAtoms);
        // TDBmatC[t][k] = TDBmatC[t][k];
      } else {
        TDBmat[t * 8 + k] = 0.0;
        //  TDBmatC[t][k] =  0.0;
        //  TDBmatM[t][k] =  0.0;
        //  TDBmatG[t][k] =  0.0;
      }
    }
  }

  // update the Output vector!
  for (unsigned int j = 0; j < 6; ++j) {
    for (unsigned int i = 0; i < 8; ++i) {
      res[j * 8 + i] = round(1000 * TDBmat[i + j * 8]) / 1000;
    }
  }

  /*
       for (unsigned int j = 0; j < 6; ++j) {
       for (unsigned int i = 0; i < 8; ++i) {
           res[j * 8 + i+48] = TDBmatC[j][i];
     }
   }


   for (unsigned int j = 0; j < 6; ++j) {
       for (unsigned int i = 0; i < 8; ++i) {
         res[j * 8 + i+96] = TDBmatM[j][i];
       }
   }

   for (unsigned int j = 0; j < 6; ++j) {
       for (unsigned int i = 0; i < 8; ++i) {
         res[j * 8 + i+144] = TDBmatG[j][i];
       }
   }

TDBmat.clear();
wp.clear();
wm.clear();
wv.clear();
wi.clear();
we.clear();
ws.clear();
w.clear();
wmean.clear();

*/
}

void Get2Dauto(double* dist, int numAtoms, const ROMol& mol,
               std::vector<double>& res) {
  get2DautocorrelationDesc(dist, numAtoms, mol, res);
}

}  // end of anonymous namespace

void AUTOCORR2D(const ROMol& mol, std::vector<double>& result) {
  int numAtoms = mol.getNumAtoms();
  double* dist = MolOps::getDistanceMat(mol, false);  // topological matrix

  result.clear();
  result.resize(48);

  Get2Dauto(dist, numAtoms, mol, result);
}
}  // end of Descriptors namespace
}  // end of RDKit namespace
