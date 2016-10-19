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


namespace RDKit {
namespace Descriptors{

namespace {



void computeCovarianceTerms(const Conformer &conf,
  const RDGeom::Point3D &center,
  double &xx, double &xy, double &xz, double &yy, double &yz, double &zz,
  bool normalize, bool ignoreHs,
  const std::vector<double> *weights ){
    PRECONDITION(!weights || weights->size()>=conf.getNumAtoms(), "bad weights vector");

    xx = xy = xz = yy = yz = zz = 0.0;
    const ROMol &mol = conf.getOwningMol();
    double wSum = 0.0;
    for (ROMol::ConstAtomIterator cai = mol.beginAtoms();
        cai != mol.endAtoms(); cai++) {
      if (((*cai)->getAtomicNum() == 1) && (ignoreHs)) {
        continue;
      }
      RDGeom::Point3D loc = conf.getAtomPos((*cai)->getIdx());
      loc -= center;
      if(weights) {
        double w = (*weights)[(*cai)->getIdx()];
        wSum += w;
        loc *= w;
      } else {
        wSum+=1.0;
      }
      xx += loc.x * loc.x;
      xy += loc.x * loc.y;
      xz += loc.x * loc.z;
      yy += loc.y * loc.y;
      yz += loc.y * loc.z;
      zz += loc.z * loc.z;

    }
    if (normalize) {
      xx /= wSum;
      xy /= wSum;
      xz /= wSum;
      yy /= wSum;
      yz /= wSum;
      zz /= wSum;
    }
}

std::vector<double> getG(int n){
  std::vector<double> res;
  for (int i=0;i<n;i++) {
    res[i] = 1+i*n/2;
  }
  return res;
}


std::vector<double> GetCharges(const ROMol& mol){

  std::vector<double> charges(mol.getNumAtoms(), 0);
  // use 12 iterations... can be more
  computeGasteigerCharges(mol, charges, 12, true);
  return charges;
}


std::vector<double> GetRelativePol(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=Pol2[mol.getAtomWithIdx(i)->getAtomicNum()]/Pol2[6];
  }

  return pol;
}


std::vector<double> GetRelativeElectroNeg(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> REN(numAtoms, 0);
  for( int i=0; i<numAtoms; ++i){

    REN[i]=ElectroNeg2[mol.getAtomWithIdx(i)->getAtomicNum()]/ElectroNeg2[6];
  }

  return REN;
}


std::vector<double> GetRelativeVdW(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> vdw(numAtoms, 0);
  for( int i=0; i<numAtoms; ++i){

    vdw[i]=VdW2[mol.getAtomWithIdx(i)->getAtomicNum()]/VdW2[6];
  }

  return vdw;
}

std::vector<double> GetAbsPol(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();
  std::vector<double> pol(numAtoms, 0);
  for( int i=0; i<numAtoms; ++i){

      pol[i]=Pol2[mol.getAtomWithIdx(i)->getAtomicNum()];
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


/*
double* GetCenterMatrix(double* Mat,int numAtoms){

    Map<MatrixXd> mat(Mat,numAtoms,3);  // row , col

    v = Xorigin.colwise().mean();

    MatrixXd X=Xorigin.rowwise() - v.transpose();

    return retreiveMat(X);
}
*/

JacobiSVD<MatrixXd> getSVDEig(double Mat[][3], int numAtoms) {





    Map<MatrixXd> mat(*Mat, numAtoms,3);
    VectorXd v = mat.colwise().mean();

    // center matrix
    MatrixXd X=mat.rowwise() - v.transpose();

    MatrixXd Weigth;

    Weigth.setIdentity(numAtoms,numAtoms);

    double weigth=Weigth.diagonal().sum();
    // svd of the centered matrix product by W:
    JacobiSVD<MatrixXd> svd(X.transpose() * Weigth * X / weigth ,  ComputeThinU | ComputeThinV);

    std::cout << "Eig done\n";

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



int GetWHIMU(const Conformer &conf, double Vpoints[]){
    int numAtoms = conf.getNumAtoms();

    
    bool weights=false;
    bool ignoreHs=false;

  RDGeom::Point3D origin(0,0,0);
  double wSum=0.0;
  for(unsigned int i=0;i<conf.getNumAtoms();++i){
    double w=1.0;
    wSum+=w;
    origin+=conf.getAtomPos(i)*w;
  }
  // std::cerr<<"  origin: "<<origin<<" "<<wSum<<std::endl;
  origin /= wSum;

  double sumXX,sumXY,sumXZ,sumYY,sumYZ,sumZZ;
  computeCovarianceTerms(conf,origin,sumXX,sumXY,sumXZ,sumYY,sumYZ,sumZZ,true,false,NULL);

  Eigen::Matrix3d mat;
  mat << sumXX, sumXY, sumXZ,
         sumXY, sumYY, sumYZ,
         sumXZ, sumYZ, sumZZ;

    std::cout << mat << "\n";

    JacobiSVD<MatrixXd> svd(mat,  ComputeThinU | ComputeThinV);
    std::cout << "-------------------" << "\n";
    std::cout << "U: "<< svd.matrixU() << "\n";
    std::cout << "V: "<< svd.matrixV() << "\n";
    std::cout << "S: "<< svd.singularValues() << "\n";
    std::cout << "-------------------" << "\n";

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    std::cout << "X original" << matorigin.transpose() << "\n";

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Weigth;

    Weigth.setIdentity(numAtoms,numAtoms);

    double weigth=Weigth.diagonal().sum();

    JacobiSVD<MatrixXd> svd2(MatOrigin.transpose() * Weigth * MatOrigin / weigth ,  ComputeThinU | ComputeThinV);
    std::cout << "-------------------" << "\n";
    std::cout << "Cov:" << MatOrigin.transpose() * Weigth * MatOrigin / weigth << "\n";
    std::cout << "U2: "<< svd2.matrixU() << "\n";
    std::cout << "V2: "<< svd2.matrixV() << "\n";
    std::cout << "S2: "<< svd2.singularValues() << "\n";
    std::cout << "-------------------" << "\n";
  
    double *SingVal = retreiveVect(svd.singularValues());

    double w[18];
    w[0]=round(SingVal[0]*1000.0)/1000.0;
    w[1]=round(SingVal[1]*1000.0)/1000.0;
    w[2]=round(SingVal[2]*1000.0)/1000.0;

    w[3]=round((SingVal[0]+SingVal[1]+SingVal[2])*1000.0)/1000.0;  // T

    w[4]=round((SingVal[0]*SingVal[1]+SingVal[0]*SingVal[2]+SingVal[1]*SingVal[2])*1000.0)/1000.0; // A

    w[5]=round((w[3]+w[4]+SingVal[0]*SingVal[1]*SingVal[2])*1000.0)/1000.0; // V

    w[6]=SingVal[0]/(SingVal[0]+SingVal[1]+SingVal[2]); // P1

    w[7]=SingVal[1]/(SingVal[0]+SingVal[1]+SingVal[2]); // p2

    w[8]=SingVal[2]/(SingVal[0]+SingVal[1]+SingVal[2]); // P3

    // get K 
    double res=0.0;
    for (int i=0;i<3;i++) 
    {
      res+=std::abs(w[i]/w[3]-1.0/3.0);
    }

    double Ku = 3.0/4.0*res;

    w[9]=round(1000.0*Ku)/1000.0; // K


    // center original matrix
    VectorXd v = MatOrigin.colwise().mean();
    MatrixXd Xmean= MatOrigin.rowwise() - v.transpose();

    MatrixXd Scores= Xmean*svd.matrixV(); // 

    VectorXd v1= Scores.col(0);
    VectorXd v2= Scores.col(1);
    VectorXd v3= Scores.col(2);

    w[10] = numAtoms*pow(w[0],2)/ v1.array().pow(4).sum(); // E1
    w[11] = numAtoms*pow(w[1],2)/ v2.array().pow(4).sum(); // E2
    w[12] = numAtoms*pow(w[2],2)/ v3.array().pow(4).sum(); // E3
    w[13] = w[10]+w[11]+w[12]; //D


     std::cout << "Scores" << Scores << "\n";

    double gamma[3]; // G

    for (int i=0;i<3; i++) {
      double ns=0.0;
      double na=0.0;

      for (int j=0;j< numAtoms;j++) {
        bool Found = false;
        for (int k=0;k< numAtoms;k++){
          if (k==j) continue;
          if (Scores(j,i) == -Scores(k,i)) {
            ns++;
            Found=true;
            break;
          }
        }
        if (!Found) na++;
      }

        std::cout << "Na:" << na << "\n";
        std::cout << "Ns:" << ns << "\n";

        gamma[i] = -1.0* ((ns / numAtoms) * log(ns / numAtoms) / log(2.0) + (na / numAtoms) * log(1.0 / numAtoms) / log(2.0));
        gamma[i] = 1.0 / (1.0 + gamma[i]);

    }

    w[14]=gamma[0]; // G1
    w[15]=gamma[1]; // G2
    w[16]=gamma[2]; // G3
    w[17]=pow(gamma[0]*gamma[1]*gamma[2],1.0/3.0);

    std::vector<std::string> descnames={"L1u","L2u","L3u","Tu","Au","Vu","P1u","P2u","P3u","Ku","E1u","E2u","E3u","Du","G1u","G2u","G3u","Gu"};


    for (int i=0;i<18;i++) {

      std::cout << descnames[i] << ":"<< w[i] << ",";
     }
    std::cout << "\n";

    return 1;
  }

} //end of anonymous namespace


double WHIM(const ROMol& mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double Vpoints[3*numAtoms];

  for(unsigned int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  double vxyz= GetWHIMU(conf, Vpoints);

  return vxyz;
}

} // end of Descriptors namespace
} // end of RDKit namespace
