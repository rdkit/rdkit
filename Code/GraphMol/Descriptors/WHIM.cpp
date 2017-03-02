//
//  Copyright (c) 2016, Guillaume GODIN
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are
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
//       products derived from this software without specific prior written permission.
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
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "WHIM.h"
#include "PBF.h"
#include "MolData3Ddescriptors.h"
#include <GraphMol/new_canon.h>

#include <Numerics/EigenSolvers/PowerEigenSolver.h>

#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace Eigen;

#define EIGEN_TOLERANCE 1.0e-2

namespace RDKit {
    namespace Descriptors{
        namespace {

          MolData3Ddescriptors moldata3D;

          double roundn(double in,int factor) {
            return round(in*pow(10,factor))/pow(10,factor);
          }

          double* retreiveMat(Eigen::MatrixXd matrix) {
             double* arrayd = matrix.data();
             return arrayd;
          }

          double* retreiveVect(Eigen::VectorXd matrix) {
             double* arrayd = matrix.data();
             return arrayd;
          }

          std::vector<double> CopyVect(VectorXd v1){
            std::vector<double> v2;
            v2.resize(v1.size());
            VectorXd::Map(&v2[0], v1.size()) = v1;
            return v2;
          }



          MatrixXd GetCenterMatrix(MatrixXd Mat){
              VectorXd v = Mat.colwise().mean();
              MatrixXd X=Mat.rowwise() - v.transpose();
              return X;
          }

          MatrixXd GetCovMatrix(MatrixXd X, MatrixXd Weigth, double weigth){
              return X.transpose() * Weigth * X / weigth;
          }

          JacobiSVD<MatrixXd> getSVD(MatrixXd Mat) {
              JacobiSVD<MatrixXd> svd(Mat,  ComputeThinU | ComputeThinV);
              return svd;
          }

          // this method use the new Canonical Ranking Atoms to gene we don't need it there but can be useful in another package 
          std::vector<unsigned int> CanonicalRankAtoms(const ROMol &mol,
                                                       bool breakTies = true,
                                                       bool includeChirality = true,
                                                       bool includeIsotopes = true) {
              std::vector<unsigned int> ranks(mol.getNumAtoms());
              Canon::rankMolAtoms(mol, ranks, breakTies, includeChirality, includeIsotopes);
              return ranks;
          }

          // center Score to the closest Atom per axis
          std::vector<double> centeringVector(std::vector<double> V){
            double D;
            double AD=1000000.0;
            for (int i = 0 ; i < V.size() ; i++){
                if (std::abs(V[i])<AD) {
                  D=V[i];
                  AD=abs(V[i]);
                }
            }
            for (int i=0;i<V.size();i++){
                  V[i]=V[i]-D;
                }
            return V;
          }

          std::vector<double> getWhimD(std::vector<double> weigthvector, MatrixXd MatOrigin, int numAtoms, double th, bool printscore) {

              double* weigtharray = &weigthvector[0];

              Map<VectorXd> Weigth(weigtharray,numAtoms);

              MatrixXd WeigthMat = Weigth.asDiagonal();

              double weigth=WeigthMat.diagonal().sum();

              MatrixXd Xmean = GetCenterMatrix(MatOrigin);

              MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

              JacobiSVD<MatrixXd> svd = getSVD(covmat);

              std::vector<double> w(18);
              // prepare data for Whim parameter computation

              double *SingVal = retreiveVect(svd.singularValues());
              MatrixXd Scores= Xmean*svd.matrixV(); //  V is similar

              // compute parameters
              w[0] = SingVal[0];
              w[1] = SingVal[1];
              w[2] = SingVal[2];
              w[3] = SingVal[0] + SingVal[1] + SingVal[2];  // T
              w[4] = SingVal[0] * SingVal[1] + SingVal[0] * SingVal[2] + SingVal[1] * SingVal[2]; // A
              w[5] = w[3] + w[4] + SingVal[0] * SingVal[1] * SingVal[2]; // V
              w[6] = SingVal[0] / (SingVal[0] + SingVal[1] + SingVal[2]); // P1
              w[7] = SingVal[1] / (SingVal[0] + SingVal[1] + SingVal[2]); // p2
              w[8] = SingVal[2] / (SingVal[0] + SingVal[1] + SingVal[2]); // P3

              double res=0.0;
              for (int i = 0 ; i < 3 ; i++)
              {
                res += std::abs( w[i] / w[3] - 1.0 / 3.0);
              }

              w[9]= 3.0 / 4.0 * res; // K

              // center original matrix

              VectorXd v1 = Scores.col(0);
              VectorXd v2 = Scores.col(1);
              VectorXd v3 = Scores.col(2);

             //  inverse of the kurtosis
              if (v1.array().pow(4).sum()>0) {
                  w[10] = numAtoms * pow( w[0] , 2) / v1.array().pow(4).sum(); // E1
              }
              else {
                w[10] =0.0;
              }

              if (v2.array().pow(4).sum()>0) {
                w[11] = numAtoms * pow( w[1] , 2) / v2.array().pow(4).sum(); // E2
              }
              else {
                w[11] =0.0;
              }

              if (v3.array().pow(4).sum()>0) {
                w[12] = numAtoms * pow( w[2] , 2) / v3.array().pow(4).sum(); // E3
              }
              else {
                w[12] =0.0;
              }

              w[13] = (w[10] + w[11] + w[12] ) / 3.0; // mean total density of the atoms called D is used on Dragon 6 not just the sum!

              // check if the molecule is fully symmetrical "like CH4" using Canonical Rank Index and/or Sphericity !


              double gamma[3]; // Gamma values
              double nAT = (double) numAtoms;

              // check if two atoms are symetric versus the new axis ie newx,newy,newz a
              for (int i = 0 ; i < 3 ; i++) {
                for (int j = 0 ; j < numAtoms ; j++) {
                  Scores(j,i) = roundn( Scores(j,i) , 2); // round the matrix! same as eigen tolerance !
                }
              }
              /*
              std::vector<double> VS;
              // need to make a copy of the vector using "CopyVect" to std:vector center and reinject to vector

              for (int i = 0 ; i < 3 ; i++) {
                  VS = centeringVector(CopyVect(Scores.col(i)));
                  for (int j = 0 ; j < numAtoms ; j++){
                    Scores(j,i) = VS[j];
                  }
              }*/
              for (int i = 0 ; i < 3 ; i++) {
                double ns=0.0;
                for (int j = 0 ; j < numAtoms-1 ; j++) {
                  for (int k = j+1 ; k < numAtoms ; k++){
                    //if (j==k) continue;
                    if (std::abs(Scores(j,i) + Scores(k,i))< th and (std::abs(Scores(j,i))>=th or std::abs(Scores(k,i))>=th )) {
                      // those that are close opposite & not close to the axis!
                        ns+=2; // check only once the symetric none null we need to add +2! (reduce the loop duration)
                      break;
                    }
                  }
                }

                // for all atoms close to the axis we need to add +1!
                for (int j = 0 ; j < numAtoms ; j++) {
                  if (std::abs(Scores(j,i))<th) {
                      ns++; // atom close to the the axis are symetric!
                  }
                }

                double na=0.0;
                na = (double) nAT - ns;

                gamma[i] =0.0;
                if (ns>0) {
                     gamma[i] = (ns / nAT) * log(ns / nAT) / log(2) + (na / nAT) * log(1.0 / nAT) / log(2); // log2 base used
                     gamma[i] = 1.0 / (1.0 - gamma[i]);
                }
                // trick to have the WHIM value always set ns to one if there is no symetry on the axis!
                if (ns==0) {
                     ns=1;
                     na=nAT-ns;
                     gamma[i] = (ns / nAT) * log(ns / nAT) / log(2) + (na / nAT) * log(1.0 / nAT) / log(2); // log2 base used
                     gamma[i] = 1.0 / (1.0 - gamma[i]);
                }

              }

              // case of complete symetry of two Components than there are set to 1!
              /*if (SingVal[0]==SingVal[1]) {
                gamma[0]=1;
                gamma[1]=1;
              }
              if (SingVal[1]==SingVal[2]) {
                gamma[1]=1;
                gamma[2]=1;
              }
              */
              w[14] = gamma[0]; // G1
              w[15] = gamma[1]; // G2
              w[16] = gamma[2]; // G3
              w[17] = pow( gamma[0] * gamma[1] * gamma[2] , 1.0 / 3.0);

              return w;
          }

          std::vector<double>  GetWHIMs(const Conformer &conf, double Vpoints[], double th){
              std::vector<double> w(18);
              std::vector<double> res(126);
              //ROMol& mol = conf.getOwningMol();
              //MolOps::removeHs(mol, false, false); // special case
              //int numAtoms = mol.getNumAtoms();
              int numAtoms = conf.getNumAtoms();
              Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);
              MatrixXd MatOrigin=matorigin.transpose();
              std::vector<double> weigthvector;

              // 18 values stored in this order : "L1u","L2u","L3u","Tu","Au","Vu","P1u","P2u","P3u","Ku","E1u","E2u","E3u","Du","G1u","G2u","G3u","Gu"
              weigthvector = moldata3D.GetUn(numAtoms);
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*0] = w[i];
              }
              w.clear();
              w.resize(18);

              weigthvector = moldata3D.GetRelativeMW(conf.getOwningMol());
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*1] = w[i];
              }
              w.clear();
              w.resize(18);

              weigthvector = moldata3D.GetRelativeVdW(conf.getOwningMol());
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*2] = w[i];
              }
              w.clear();
              w.resize(18);

              weigthvector = moldata3D.GetRelativeENeg(conf.getOwningMol());
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*3] = w[i];
              }
              w.clear();
              w.resize(18);

              weigthvector = moldata3D.GetRelativePol(conf.getOwningMol());
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*4] = w[i];
              }
              w.clear();
              w.resize(18);

              weigthvector = moldata3D.GetRelativeIonPol(conf.getOwningMol());
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*5] = w[i];
              }

              w.clear();
              w.resize(18);
              weigthvector = moldata3D.GetIState(conf.getOwningMol()); // caution not only neighours hum not sure on this based on the paper! Also this should be only on Heavy atoms
              w= getWhimD(weigthvector, MatOrigin, numAtoms, th, false);
              for (int i = 0 ; i < 18 ; i++) {
                res[i+18*6] = w[i];
              }
              w.clear();
              return res;
          }
        } //end of anonymous namespace


        std::vector<double> WHIM(const ROMol& mol, int confId, double th){
            PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
            int numAtoms = mol.getNumAtoms();
            //double pbf=RDKit::Descriptors::PBF(mol);
            const Conformer &conf = mol.getConformer(confId);
            double Vpoints[3*numAtoms];

            for(int i=0; i<numAtoms; ++i){
               Vpoints[3*i]   =conf.getAtomPos(i).x;
               Vpoints[3*i+1] =conf.getAtomPos(i).y;
               Vpoints[3*i+2] =conf.getAtomPos(i).z;
            }

            std::vector<double> res(114);
            std::vector<double>  w(126);
            w= GetWHIMs(conf, Vpoints, th);

            // takes only L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u
            int map1[11] = {0,1,2,6,7,14,15,16,10,11,12};

            for (int k = 0 ; k < 7 ; k++) {
              for (int i = 0 ; i < 11 ; i++) {
                  res[i + 11 * k ] = roundn( w[ map1[i] + 18 * k ] , 3 );
              }
            }

            for (int i = 0 ; i < 2; i++) {
              res[i + 13*7] = roundn(w[17 + 18* i], 3);  //92  93 //Gu  Gm
            }

            for (int i = 0 ; i < 7; i++) {
              res[i + 11*7] = roundn(w[3 + 18* i], 3);     //78  79  80  81  82  83  84 //Tu  Tm  Tv  Te  Tp  Ti  Ts
              res[i + 12*7] = roundn(w[4 + 18* i], 3);     //85  86  87  88  89  90  91 //Au  Am  Av  Ae  Ap  Ai  As
              res[i + 13*7+2] = roundn(w[9 + 18* i], 3);   //94  95  96  97  98  99  100 //Ku  Km  Kv  Ke  Kp  Ki  Ks
              res[i + 14*7+2] = roundn(w[13 + 18* i], 3);  //101 102 103 104 105 106 107 //Du  Dm  Dv  De  Dp  Di  Ds
              res[i + 15*7+2] = roundn(w[5 + 18* i], 3);   //108 109 110 111 112 113 114 //Vu  Vm  Vv  Ve  Vp  Vi  Vs
            }
          // res output are : L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u L1m L2m L3m P1m P2m G1m G2m G3m E1m E2m E3m L1v L2v L3v P1v P2v G1v G2v G3v E1v E2v E3v L1e L2e L3e P1e P2e G1e G2e G3e E1e E2e E3e L1p L2p L3p P1p P2p G1p G2p G3p E1p E2p E3p L1i L2i L3i P1i P2i G1i G2i G3i E1i E2i E3i L1s L2s L3s P1s P2s G1s G2s G3s E1s E2s E3s Tu  Tm  Tv  Te  Tp  Ti  Ts  Au  Am  Av  Ae  Ap  Ai  As  Gu  Gm  Ku  Km  Kv  Ke  Kp  Ki  Ks  Du  Dm  Dv  De  Dp  Di  Ds  Vu  Vm  Vv  Ve  Vp  Vi  Vs
            return res;
          }
    } // end of Descriptors namespace
} // end of RDKit namespace