//
//  Copyright (c) 2012, Institue of Cancer Research.
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
// For more information on the Plane of Best Fit please see http://pubs.acs.org/doi/abs/10.1021/ci300293f
//
//  If this code has been useful to you, please include the reference
//  in any work which has made use of it:

//  Plane of Best Fit: A Novel Method to Characterize the Three-Dimensionality of Molecules, Nicholas C. Firth, Nathan Brown, and Julian Blagg, Journal of Chemical Information and Modeling 2012 52 (10), 2516-2525

//
//
// Created by Nicholas Firth, November 2011
// Modified by Greg Landrum for inclusion in the RDKit distribution November 2012
// Further modified by Greg Landrum for inclusion in the RDKit core September 2016
// Adding RBF descriptors to 3D descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "WHIM.h"

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"

#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>


using namespace Eigen;


double Pol2[]={0.67,0,24.3,5.60,3.03,1.76,1.10,0.80,0.56,0,23.6,10.6,6.80,5.38,3.63,2.90,2.18,0,43.4,22.8,0,0,0,11.60,9.40,8.40,7.50,6.80,6.10,7.10,8.12,6.07,4.31,3.73,3.05,0,47.3,27.6,0,0,0,12.80,0,0,0,0,7.20,7.20,10.20,7.70,6.60,5.50,5.35,0,0,0,0,0,0,0,0,0,0,23.50,0,0,0,0,0,0,0,0,0,0,0,0,0,6.50,5.80,5.70,7.60,6.80,7.40};
double ElectroNeg2[]={2.59, 0, 0.89, 1.81, 2.28, 2.75, 3.19, 3.65, 4.0, 0, 0.56, 1.32, 1.71, 2.14, 2.52, 2.96, 3.48, 0, 0.45, 0.95, 0, 0, 0, 1.66, 2.2, 2.2, 2.56, 1.94, 1.95, 2.23, 2.42, 2.62, 2.82, 3.01, 3.22, 0, 0.31, 0.72, 0, 0, 0, 1.15, 0, 0, 0, 0, 1.83, 1.98, 2.14, 2.3, 2.46, 2.62, 2.78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.28, 2.65, 2.2, 2.25, 2.29, 2.34};
double VdW2[]={6.71, 0, 25.25, 0.0, 17.88, 22.45, 15.6, 11.49, 9.2, 0, 49.0, 21.69, 36.51, 31.98, 26.52, 24.43, 22.45, 0, 87.11, 0.0, 0, 0, 0, 44.6, 43.4, 41.05, 35.04, 17.16, 11.49, 11.25, 27.39, 28.73, 26.52, 28.73, 31.06, 0, 0.0, 0.0, 0, 0, 0, 33.51, 0, 0, 0, 0, 21.31, 16.52, 30.11, 45.83, 38.79, 36.62, 38.79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 72.78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 22.45, 19.16, 15.6, 31.54, 34.53, 38.79};

double relativeMw[]={0.084,0,0,0,0.900,1.000,1.166,1.332,1.582,0,0,0, 2.246, 2.339, 2.579,2.670, 2.952,0,0,0,0,0,0,0,0, 4.650, 4.907, 4.887, 5.291, 5.445,0,0,0,0, 6.653,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 9.884,0,0, 10.56};
double relativeVdW[]={0.299,0,0,0,0.796,1.000,0.695,0.512,0.410,0,0,0, 1.626, 1.424, 1.181,1.088, 1.035,0,0,0,0,0,0,0,0, 1.829, 1.561, 0.764, 0.512, 1.708,0,0,0,0, 1.384,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 2.042,0,0, 1.728};
double relativeNeg[]={0.944,0,0,0,0.828,1.000,1.163,1.331,1.457,0,0,0, 0.624, 0.779, 0.916,1.077, 1.265,0,0,0,0,0,0,0,0, 0.728, 0.728, 0.728, 0.740, 0.810,0,0,0,0, 1.172,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.837,0,0, 1.012};
double relativePol[]={0.379,0,0,0,1722,1.000,0.625,0.456,0.316,0,0,0, 3.864, 3.057, 2.063,1.648, 1.239,0,0,0,0,0,0,0,0, 4.773, 4.261, 3.864, 3.466, 4.034,0,0,0,0, 1.733,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 4.375,0,0, 3.040};



namespace RDKit {
namespace Descriptors{

namespace {

double roundn(double in,int factor) {

  return round(in*pow(10,factor))/pow(10,factor);
}


std::vector<double> GetCharges(const ROMol& mol){

  std::vector<double> charges(mol.getNumAtoms(), 0);
  // use 12 iterations... can be more
  computeGasteigerCharges(mol, charges, 12, true);
  return charges;
}


std::vector<double> GetRelativeMW(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativeMw[mol.getAtomWithIdx(i)->getAtomicNum()-1];
  }


  return pol;
}


std::vector<double> GetRelativeENeg(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> neg(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    neg[i]=relativeNeg[mol.getAtomWithIdx(i)->getAtomicNum()-1];
}

  return neg;
}


std::vector<double> GetRelativeVdW(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> vdw(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    vdw[i]=relativeVdW[mol.getAtomWithIdx(i)->getAtomicNum()-1];

  }

  return vdw;
}

std::vector<double> GetRelativePol(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativePol[mol.getAtomWithIdx(i)->getAtomicNum()-1];

  }

  return pol;
}


double* retreiveMat(Eigen::MatrixXd matrix) {
   double* arrayd = matrix.data();
   return arrayd;
}

double* retreiveVect(Eigen::VectorXd matrix) {
   double* arrayd = matrix.data();
   return arrayd;

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





double getKu(double* w, int numAtoms) {

  double res=0.0;
  double sums=0.0;
  for (int i=0;i<numAtoms;i++) 
  {
      sums+=w[i];
  }

  for (int i=0;i<numAtoms;i++) 
  {
    res+=std::abs(i/sums-1.0/3.0);
  }

  double Ku = 3.0/4.0*res;
  return round(1000*Ku)/1000;
}

double* getWhimDesc(JacobiSVD<MatrixXd> svd, MatrixXd Xmean, int numAtoms, double th) {

    double *w = new double[18];
    // prepare data for Whim parameter computation
    
    double *SingVal = retreiveVect(svd.singularValues());
    MatrixXd Scores= Xmean*svd.matrixV(); //  V is similar

    // compute parameters
    w[0]=SingVal[0];
    w[1]=SingVal[1];
    w[2]=SingVal[2];
    w[3]=SingVal[0]+SingVal[1]+SingVal[2];  // T
    w[4]=SingVal[0]*SingVal[1]+SingVal[0]*SingVal[2]+SingVal[1]*SingVal[2]; // A
    w[5]=w[3]+w[4]+SingVal[0]*SingVal[1]*SingVal[2]; // V
    w[6]=SingVal[0]/(SingVal[0]+SingVal[1]+SingVal[2]); // P1
    w[7]=SingVal[1]/(SingVal[0]+SingVal[1]+SingVal[2]); // p2
    w[8]=SingVal[2]/(SingVal[0]+SingVal[1]+SingVal[2]); // P3

    double res=0.0;
    for (int i=0;i<3;i++) 
    {
      res+=std::abs(w[i]/w[3]-1.0/3.0);
    }

    w[9]=3.0/4.0*res; // K

    // center original matrix

    VectorXd v1= Scores.col(0);
    VectorXd v2= Scores.col(1);
    VectorXd v3= Scores.col(2);


    w[10] = numAtoms*pow(w[0],2)/ v1.array().pow(4).sum(); // E1
    w[11] = numAtoms*pow(w[1],2)/ v2.array().pow(4).sum(); // E2
    w[12] = numAtoms*pow(w[2],2)/ v3.array().pow(4).sum(); // E3
    w[13] = w[10]+w[11]+w[12]; //D


    double gamma[3]; // Gamma values
    double nAT = (double) numAtoms;
    // check if two atoms are symetric versus the new axis ie newx,newy,newz a
    for (int i=0;i<3; i++) {
      double ns=0.0;
      for (int j=0;j< numAtoms-1;j++) {
        for (int k=j+1;k< numAtoms;k++){
          if (std::abs(roundn(Scores(j,i),3)+ roundn(Scores(k,i),3))<2*th and std::abs(roundn(Scores(k,i),3))>2.0*th  and std::abs(roundn(Scores(j,i),3))>2.0*th) {  // those that are close opposite but close to the axis!
            ns+=2;
            break;
          }
        }
        if (std::abs(roundn(Scores(j,i),3))<th) {
        ns++; // atom close to the the axis are symetric! 
        std::cout << "close to axis,";
        } 
      }

      double na=0.0;
      na = nAT-ns;
      std::cout << "ns:" << ns << ",";
      gamma[i] =0.0;
      if (ns>0) {
           gamma[i] = ((ns / nAT) * log(ns / nAT) / log(2.0) + (na / nAT) * log(1.0 / nAT) / log(2.0));
           gamma[i] = 1.0 / (1.0 - gamma[i]);
      }
    }
    std::cout  << "\n";
    w[14]=gamma[0]; // G1
    w[15]=gamma[1]; // G2
    w[16]=gamma[2]; // G3
    w[17]=pow(gamma[0]*gamma[1]*gamma[2],1.0/3.0);

    return w;

}


int GetWHIMU(const Conformer &conf, double Vpoints[], double th){
    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Weigth;

    Weigth.setIdentity(numAtoms,numAtoms);

    double weigth=Weigth.diagonal().sum();

    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,Weigth,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    double *w= getWhimDesc(svd, Xmean, numAtoms, th);

    std::vector<std::string> descnames={"L1u","L2u","L3u","Tu","Au","Vu","P1u","P2u","P3u","Ku","E1u","E2u","E3u","Du","G1u","G2u","G3u","Gu"};

    for (int i=0;i<18;i++) {

      std::cout << descnames[i] << ":"<< roundn(w[i],3) << ",";
     }
    std::cout << "\n";



    return 1;
  }



int GetWHIMMass(const Conformer &conf, double Vpoints[], double th){
    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = GetRelativeMW(conf.getOwningMol());

    double* weigtharray = &weigthvector[0];

    Map<VectorXd> Weigth(weigtharray,numAtoms);

    MatrixXd WeigthMat = Weigth.asDiagonal();


    double weigth=WeigthMat.diagonal().sum();


    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    double *w= getWhimDesc(svd, Xmean, numAtoms, th);

    std::vector<std::string> descnames={"L1m","L2m","L3m","Tm","Am","Vm","P1m","P2m","P3m","Km","E1m","E2m","E3m","Dm","G1m","G2m","G3m","Gm"};

    for (int i=0;i<18;i++) {

      std::cout << descnames[i] << ":"<< roundn(w[i],3) << ",";
     }
    std::cout << "\n";


    return 1;
  }

int GetWHIMvdw(const Conformer &conf, double Vpoints[], double th){
    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = GetRelativeVdW(conf.getOwningMol());

    double* weigtharray = &weigthvector[0];

    Map<VectorXd> Weigth(weigtharray,numAtoms);

    MatrixXd WeigthMat = Weigth.asDiagonal();


    double weigth=WeigthMat.diagonal().sum();


    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    double *w= getWhimDesc(svd, Xmean, numAtoms, th);

     std::vector<std::string> descnames={"L1v","L2v","L3v","Tv","Av","Vv","P1v","P2v","P3v","Kv","E1v","E2v","E3v","Dv","G1v","G2v","G3v","Gv"};

    for (int i=0;i<18;i++) {

      std::cout << descnames[i] << ":"<< roundn(w[i],3) << ",";
     }
    std::cout << "\n";


    return 1;
  }


int GetWHIMneg(const Conformer &conf, double Vpoints[], double th){
    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = GetRelativeENeg(conf.getOwningMol());

    double* weigtharray = &weigthvector[0];

    Map<VectorXd> Weigth(weigtharray,numAtoms);

    MatrixXd WeigthMat = Weigth.asDiagonal();


    double weigth=WeigthMat.diagonal().sum();


    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    double *w= getWhimDesc(svd, Xmean, numAtoms, th);

     std::vector<std::string> descnames={"L1e","L2e","L3e","Te","Ae","Ve","P1e","P2e","P3e","Ke","E1e","E2e","E3e","De","G1e","G2e","G3e","Ge"};

    for (int i=0;i<18;i++) {

      std::cout << descnames[i] << ":"<< roundn(w[i],3) << ",";
     }
    std::cout << "\n";


    return 1;
  }




int GetWHIMpol(const Conformer &conf, double Vpoints[], double th){
    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = GetRelativePol(conf.getOwningMol());

    double* weigtharray = &weigthvector[0];

    Map<VectorXd> Weigth(weigtharray,numAtoms);

    MatrixXd WeigthMat = Weigth.asDiagonal();


    double weigth=WeigthMat.diagonal().sum();


    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    double *w= getWhimDesc(svd, Xmean, numAtoms, th);

     std::vector<std::string> descnames={"L1p","L2p","L3p","Tp","Ap","Vp","P1p","P2p","P3p","Kp","E1p","E2p","E3p","Dp","G1p","G2p","G3p","Gp"};

    for (int i=0;i<18;i++) {

      std::cout << descnames[i] << ":"<< roundn(w[i],3) << ",";
     }
    std::cout << "\n";


    return 1;
  }




} //end of anonymous namespace






double WHIM(const ROMol& mol,int confId, double th){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double Vpoints[3*numAtoms];

  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  double vxyz= GetWHIMU(conf, Vpoints, th);
  double vxyz1= GetWHIMMass(conf, Vpoints, th);
  double vxyz2= GetWHIMvdw(conf, Vpoints, th);
  double vxyz3= GetWHIMneg(conf, Vpoints, th);
  double vxyz4= GetWHIMpol(conf, Vpoints, th);

  return vxyz;
}

} // end of Descriptors namespace
} // end of RDKit namespace
