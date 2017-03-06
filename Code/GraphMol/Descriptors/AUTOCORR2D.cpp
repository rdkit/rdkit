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
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "AUTOCORR2D.h"
#include "MolData3Ddescriptors.h"

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <Numerics/EigenSolvers/PowerEigenSolver.h>

#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <iostream>


using namespace Eigen;
namespace RDKit {
    namespace Descriptors {

        namespace {

            MolData3Ddescriptors moldata3D;



          // this is the Broto-Moreau 2D descriptors (centered or not)
          void get2DautocorrelationDesc(double* dist, int numAtoms,
                           const ROMol& mol, std::vector<double>& res) {


                            std::vector<double> wp = moldata3D.GetRelativePol(mol);
                            std::vector<double> wm = moldata3D.GetRelativeMW(mol);
                            std::vector<double> wv = moldata3D.GetRelativeVdW(mol);
                            std::vector<double> wi = moldata3D.GetRelativeIonPol(mol);
                            std::vector<double> we = moldata3D.GetRelativeENeg(mol);
                            std::vector<double> ws = moldata3D.GetIState(mol);
                            double w[6][numAtoms];
                            std::vector<double> wmean(6,0);


                            for (unsigned int i = 0; i < numAtoms; i++) {
                                w[0][i] = wm[i];
                                w[1][i] = wv[i];
                                w[2][i] = we[i];
                                w[3][i] = wp[i];
                                w[4][i] = wi[i];
                                w[5][i] = ws[i];
                            }

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

                            double TDBmat[6][8] = { {0.0} };
                            double TDBmatM[6][8] = { {0.0} };
                            double TDBmatG[6][8] = { {0.0} };
                            double TDBmatC[6][8] = { {0.0} };

                            for (unsigned int k = 0; k < 8; k++) {
                                int maxkVertexPairs = 0;
                                for (unsigned int i = 0; i < numAtoms - 1 ; ++i)
                                {
                                    for (unsigned int j = i + 1; j < numAtoms; ++j)
                                    {
                                        if (dist[j * numAtoms + i]==k + 1)
                                        {
                                            for (unsigned int t = 0; t < 6; ++t)
                                            {
                                                TDBmatM[t][k] += (w[t][i]-wmean[t]) * (w[t][j]-wmean[t]);
                                                TDBmatG[t][k] += (w[t][i] - w[t][j]) * (w[t][i] - w[t][j]);
                                                TDBmat[t][k] += w[t][i] * w[t][j];
                                                TDBmatC[t][k] += abs(w[t][i]-wmean[t]) * abs(w[t][j]-wmean[t]);
                                            }
                                            maxkVertexPairs += 1;
                                        }
                                    }
                                }

                                for (unsigned int t = 0; t < 6; ++t)
                                {
                                    if (maxkVertexPairs>0) {
                                        TDBmat[t][k] = log(TDBmat[t][k]+1) ;
                                        TDBmatG[t][k] = TDBmatG[t][k]/squaresumdiff[t]/maxkVertexPairs/2/(numAtoms-1);
                                        TDBmatM[t][k] = TDBmatM[t][k]/squaresumdiff[t]/maxkVertexPairs/(numAtoms);
                                        TDBmatC[t][k] = TDBmatC[t][k];
                                    }
                                    else {
                                        TDBmat[t][k] =  0.0;
                                        TDBmatC[t][k] =  0.0;
                                        TDBmatM[t][k] =  0.0;
                                        TDBmatG[t][k] =  0.0;
                                    }
                                }
                            }

                            for (unsigned int j = 0; j < 6; ++j) {
                                for (unsigned int i = 0; i < 8; ++i) {
                                  std::cout <<   TDBmatC[j][i]  << ",";
                                  std::cout <<   TDBmatM[j][i]  << ",";
                                  std::cout <<   TDBmatG[j][i]  << ",";
                                  std::cout <<   TDBmat[j][i]  << ",";

                                }
                            }



                            // update the Output vector!
                            for (unsigned int j = 0; j < 6; ++j) {
                                for (unsigned int i = 0; i < 8; ++i) {
                                    res[j * 8 + i] = TDBmat[j][i];


                                }
                            }

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

                        }




            void Get2Dauto(double* topologicaldistance, int numAtoms, const ROMol& mol,
                           std::vector<double>& res) {
                // AUTOCORRNAMES={"TDB01u","TDB02u","TDB03u","TDB04u","TDB05u","TDB06u","TDB07u","TDB08u","TDB09u","TDB10u","TDB01m","TDB02m","TDB03m","TDB04m","TDB05m","TDB06m","TDB07m","TDB08m","TDB09m","TDB10m","TDB01v","TDB02v","TDB03v","TDB04v","TDB05v","TDB06v","TDB07v","TDB08v","TDB09v","TDB10v","TDB01e","TDB02e","TDB03e","TDB04e","TDB05e","TDB06e","TDB07e","TDB08e","TDB09e","TDB10e","TDB01p","TDB02p","TDB03p","TDB04p","TDB05p","TDB06p","TDB07p","TDB08p","TDB09p","TDB10p","TDB01i","TDB02i","TDB03i","TDB04i","TDB05i","TDB06i","TDB07i","TDB08i","TDB09i","TDB10i","TDB01s","TDB02s","TDB03s","TDB04s","TDB05s","TDB06s","TDB07s","TDB08s","TDB09s","TDB10s","TDB01r","TDB02r","TDB03r","TDB04r","TDB05r","TDB06r","TDB07r","TDB08r","TDB09r","TDB10r"};

                get2DautocorrelationDesc(topologicaldistance, numAtoms, mol, res);
            }

        }  // end of anonymous namespace

        void AUTOCORR2D(const ROMol& mol, std::vector<double>& res, int confId) {
            int numAtoms = mol.getNumAtoms();
            const Conformer& conf = mol.getConformer(confId);
            double* topologicaldistance = MolOps::getDistanceMat(mol, false);  // topological matrix

            res.clear();
            res.resize(192);
            //Get3Dauto(dist3D, topologicaldistance, numAtoms, mol, res);
            Get2Dauto( topologicaldistance, numAtoms, mol, res);

        }
    }  // end of Descriptors namespace
}  // end of RDKit namespace
