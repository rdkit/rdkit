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
// Adding WHIM descriptors to 3D descriptors by Guillaume Godin
// for build & set RDBASE! => export RDBASE=/Users/mbp/Github/rdkit_mine/

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "WHIM.h"

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
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

double relativeMw[]={0.084,0,0,0,0.900,1.000,1.166,1.332,1.582,0,0,0, 2.246, 2.339, 2.579,2.670, 2.952,0,0,0,0,0,0,0,0, 4.650, 4.907, 4.887, 5.291, 5.445,0,0,0,0, 6.653,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 9.884,0,0, 10.56};
double relativeVdW[]={0.299,0,0,0,0.796,1.000,0.695,0.512,0.410,0,0,0, 1.626, 1.424, 1.181,1.088, 1.035,0,0,0,0,0,0,0,0, 1.829, 1.561, 0.764, 0.512, 1.708,0,0,0,0, 1.384,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 2.042,0,0, 1.728};
double relativeNeg[]={0.944,0,0,0,0.828,1.000,1.163,1.331,1.457,0,0,0, 0.624, 0.779, 0.916,1.077, 1.265,0,0,0,0,0,0,0,0, 0.728, 0.728, 0.728, 0.740, 0.810,0,0,0,0, 1.172,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.837,0,0, 1.012};
double relativePol[]={0.379,0,0,0,1722,1.000,0.625,0.456,0.316,0,0,0, 3.864, 3.057, 2.063,1.648, 1.239,0,0,0,0,0,0,0,0, 4.773, 4.261, 3.864, 3.466, 4.034,0,0,0,0, 1.733,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 4.375,0,0, 3.040};
double relativeIonPol[]={1.208,0,0.479,0.828,0.737,1,1.291,1.209,1.547,0,0.456,0.679,0.532,0.724,0.931,0.92,1.152,0,0.386,0.543,0,0,0,0.601,0.66,0.702,0.7,0.679,0.686,0.834,0.533,0.702,0.872,0.866,1.049,0,0.371,0.506,0,0,0,0.63,0,0,0,0,0.673,0.799,0.514,0.652,0.767,0.8,0.928,0,0,0,0,0,0,0,0,0,0,0.546,0,0,0,0,0,0,0,0,0,0,0,0,0,0.799,0.819,0.927,0.542,0.659,0.647};


namespace RDKit {
namespace Descriptors{

namespace {

double roundn(double in,int factor) {

  return round(in*pow(10,factor))/pow(10,factor);
}

int GetPrincipalQuantumNumber(int AtomicNum) {
  if (AtomicNum<=2) return 1;
  else if (AtomicNum<=10) return 2;
  else if (AtomicNum<=18) return 3;
  else if (AtomicNum<=36) return 4;
  else if (AtomicNum<=54) return 5;
  else if (AtomicNum<=86) return 6;
  else return 7;
}

// adaptation from EState.py 
// we need the Is value only there
std::vector<double> GetIState(const ROMol &mol){
  int numAtoms = mol.getNumAtoms();
  std::vector<double> Is;

  for (int i = 0; i < numAtoms; ++i) {
    const RDKit::Atom * atom= mol.getAtomWithIdx(i);
    int atNum=atom->getAtomicNum();
    int d = atom->getDegree();
    if (d>0 and atNum>1) {
      int h = atom->getTotalNumHs();
      int dv = RDKit::PeriodicTable::getTable()->getNouterElecs(atNum)-h;
      int N = GetPrincipalQuantumNumber(atNum);
      Is.push_back(round(1000*(4.0/(N*N)*dv+1.0)/d)/1000);  // WHIM-P5.pdf paper 1997  => +7 & NoHydrogens is used!
    }
    else Is.push_back(0.0);
 }

  double tmp,p;
  double *dist = MolOps::getDistanceMat(mol,false,false); // topological distance not 3D distance!!!
  double accum[numAtoms];
  for (int i=0;i<numAtoms;i++) {
    for (int j=i+1;j<numAtoms;j++) {
       p = dist[i * numAtoms + j]+1;

      if (p < 1e6) {
        tmp = (Is[i] - Is[j]) / (p * p);
        accum[i] += tmp;
        accum[j] -= tmp;
      }
    }
  }

 for (int i=0;i<numAtoms;i++) {
    Is[i]+=accum[i];
 }


 return Is;
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

std::vector<double> GetRelativeIonPol(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativeIonPol[mol.getAtomWithIdx(i)->getAtomicNum()-1];

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

double* getWhimDesc(JacobiSVD<MatrixXd> svd, MatrixXd Xmean, int numAtoms, double th,bool printscore) {

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

if (printscore and th ==0.03) {
      std::cout << v1.transpose() << "\n";
      std::cout << v2.transpose() << "\n";
      std::cout << v3.transpose() << "\n";
}




    w[10] = numAtoms*pow(w[0],2)/ v1.array().pow(4).sum(); // E1
    w[11] = numAtoms*pow(w[1],2)/ v2.array().pow(4).sum(); // E2
    w[12] = numAtoms*pow(w[2],2)/ v3.array().pow(4).sum(); // E3
    w[13] = w[10]+w[11]+w[12]; //D


    double gamma[3]; // Gamma values
    double nAT = (double) numAtoms;
    // check if two atoms are symetric versus the new axis ie newx,newy,newz a
    for (int i=0;i<3; i++) {
      for (int j =0;j<numAtoms;j++) {
        Scores(j,i)=roundn(Scores(j,i),2); // round the matrix! same as eigen tolerance !
      }
    }

    //std::cout << Scores.transpose() <<"\n";

    for (int i=0;i<3; i++) {
      double ns=0.0;
      for (int j=0;j< numAtoms;j++) {
        for (int k=0;k< numAtoms;k++){

          if (j==k) continue;

          if (std::abs(Scores(j,i) + Scores(k,i))< th and (std::abs(Scores(j,i))>0 or std::abs(Scores(k,i))>0 )) {  // those that are close opposite & not close to the axis!
            ns++; // check only once the symetric none null we need to add +2!
            break;
          }
        }
      }

      // for all atoms close to the axis we need to add +1!
      for (int j=0;j< numAtoms;j++) {
        if (Scores(j,i)==0) {
        ns++; // atom close to the the axis are symetric! 
        } 
      }

      double na=0.0;
      na = nAT-ns;
    //  std::cout << "ns:" << ns << ",";
      gamma[i] =0.0;
      if (ns>0) {
           gamma[i] = ((ns / nAT) * log(ns / nAT) / log(2.0) + (na / nAT) * log(1.0 / nAT) / log(2.0));
           gamma[i] = 1.0 / (1.0 - gamma[i]);
      }
    }
  //  std::cout  << "\n";
    w[14]=gamma[0]; // G1
    w[15]=gamma[1]; // G2
    w[16]=gamma[2]; // G3
    w[17]=pow(gamma[0]*gamma[1]*gamma[2],1.0/3.0);

    return w;

}


double* GetWHIMU(const Conformer &conf, double Vpoints[], double th){

    double *w = new double[18];

    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Weigth;

    Weigth.setIdentity(numAtoms,numAtoms);

    double weigth=Weigth.diagonal().sum();

    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,Weigth,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);


    /* // compare both methods results identicals!
        Matrix3d covmat3 =covmat;
        std::cout << svd.matrixV() << "\n";
        std::cout << svd.singularValues() << "\n";


        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covmat3);
        if(eigensolver.info()!=Eigen::Success){
          BOOST_LOG(rdErrorLog)<<"eigenvalue calculation did not converge"<<std::endl;
          return 0;
        }


        Eigen::Matrix3d axes = eigensolver.eigenvectors();
        Eigen::Vector3d moments = eigensolver.eigenvalues();

        std::cout << axes << "\n";
        std::cout << moments << "\n";
    */


    w= getWhimDesc(svd, Xmean, numAtoms, th, true);

  


    return w;
  }



double*  GetWHIMMass(const Conformer &conf, double Vpoints[], double th){
      double *w = new double[18];

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

    w= getWhimDesc(svd, Xmean, numAtoms, th,false);



    return w;
  }

double*  GetWHIMvdw(const Conformer &conf, double Vpoints[], double th){
    double *w = new double[18];

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

    w= getWhimDesc(svd, Xmean, numAtoms, th,false);



    return w;
  }


double*  GetWHIMneg(const Conformer &conf, double Vpoints[], double th){
    double *w = new double[18];

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

    w= getWhimDesc(svd, Xmean, numAtoms, th, false);



    return w;
  }




double*  GetWHIMpol(const Conformer &conf, double Vpoints[], double th){
    double *w = new double[18];

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

    w= getWhimDesc(svd, Xmean, numAtoms, th, false);



    return w;
  }



double*  GetWHIMIonPol(const Conformer &conf, double Vpoints[], double th){
      double *w = new double[18];

    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = GetRelativeIonPol(conf.getOwningMol());

    double* weigtharray = &weigthvector[0];

    Map<VectorXd> Weigth(weigtharray,numAtoms);

    MatrixXd WeigthMat = Weigth.asDiagonal();


    double weigth=WeigthMat.diagonal().sum();


    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    w= getWhimDesc(svd, Xmean, numAtoms, th,false);



    return w;
  }


double*  GetWHIMIState(const Conformer &conf, double Vpoints[], double th){
    double *w = new double[18];

    ROMol& mol = conf.getOwningMol();

    MolOps::removeHs(mol, false, false);

    int numAtoms = mol.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = GetIState(mol);

    double* weigtharray = &weigthvector[0];

    Map<VectorXd> Weigth(weigtharray,numAtoms);

    MatrixXd WeigthMat = Weigth.asDiagonal();


    double weigth=WeigthMat.diagonal().sum();


    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,WeigthMat,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);

    w= getWhimDesc(svd, Xmean, numAtoms, th,false);

    return w;
  }



} //end of anonymous namespace



double* WHIM(const ROMol& mol,int confId, double th){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double Vpoints[3*numAtoms];

  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  double *res= new double[114];

  double* wu= GetWHIMU(conf, Vpoints, th);
  double* wm= GetWHIMMass(conf, Vpoints, th);
  double* wv= GetWHIMvdw(conf, Vpoints, th);
  double* we= GetWHIMneg(conf, Vpoints, th);
  double* wp= GetWHIMpol(conf, Vpoints, th);
  double* wi= GetWHIMIonPol(conf, Vpoints, th);
  double* ws= GetWHIMIState(conf, Vpoints, th);
 
  // takes only L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u
  int map1[11] = {0,1,2,6,7,14,15,16,10,11,12};
  // std::vector<std::string> descnames={"L1u","L2u","L3u","Tu","Au","Vu","P1u","P2u","P3u","Ku","E1u","E2u","E3u","Du","G1u","G2u","G3u","Gu"};
    for (int i=0;i<11;i++) {
        res[i]=roundn(wu[map1[i]],3);

     }

   // descnames={"L1m","L2m","L3m","Tm","Am","Vm","P1m","P2m","P3m","Km","E1m","E2m","E3m","Dm","G1m","G2m","G3m","Gm"};

    for (int i=0;i<11;i++) {
        res[i+11]=roundn(wm[map1[i]],3);

     }

  //   descnames={"L1v","L2v","L3v","Tv","Av","Vv","P1v","P2v","P3v","Kv","E1v","E2v","E3v","Dv","G1v","G2v","G3v","Gv"};

    for (int i=0;i<11;i++) {
        res[i+11*2] = roundn(wv[map1[i]],3);
     }

   // descnames={"L1e","L2e","L3e","Te","Ae","Ve","P1e","P2e","P3e","Ke","E1e","E2e","E3e","De","G1e","G2e","G3e","Ge"};

    for (int i=0;i<11;i++) {
        res[i+11*3] = roundn(we[map1[i]],3);
     }


   // descnames={"L1p","L2p","L3p","Tp","Ap","Vp","P1p","P2p","P3p","Kp","E1p","E2p","E3p","Dp","G1p","G2p","G3p","Gp"};

    for (int i=0;i<11;i++) {
        res[i+11*4] = roundn(wp[map1[i]],3);

     }
    

    for (int i=0;i<11;i++) {
        res[i+11*5] = roundn(wi[map1[i]],3);

     }


    for (int i=0;i<11;i++) {
        res[i+11*6] = roundn(ws[map1[i]],3);

     }

       //78  79  80  81  82  83  84
       //Tu  Tm  Tv  Te  Tp  Ti  Ts
      res[77]=roundn(wu[3],3);
      res[78]=roundn(wm[3],3);
      res[79]=roundn(wv[3],3);
      res[80]=roundn(we[3],3);
      res[81]=roundn(wp[3],3);
      res[82]=roundn(wi[3],3);
      res[83]=roundn(ws[3],3);

      //85  86  87  88  89  90  91
      //Au  Am  Av  Ae  Ap  Ai  As
      res[84]=roundn(wu[4],3);
      res[85]=roundn(wm[4],3);
      res[86]=roundn(wv[4],3);
      res[87]=roundn(we[4],3);
      res[88]=roundn(wp[4],3);
      res[89]=roundn(wi[4],3);
      res[90]=roundn(ws[4],3);

      //92  93
      //Gu  Gm
      res[91]=roundn(wu[17],3);
      res[92]=roundn(wm[17],3);

      //94  95  96  97  98  99  100
      //Ku  Km  Kv  Ke  Kp  Ki  Ks
      res[93]=roundn(wu[9],3);
      res[94]=roundn(wm[9],3);
      res[95]=roundn(wv[9],3);
      res[96]=roundn(we[9],3);
      res[97]=roundn(wp[9],3);
      res[98]=roundn(wi[9],3);
      res[99]=roundn(ws[9],3);

      //101 102 103 104 105 106 107
       //Du  Dm  Dv  De  Dp  Di  Ds
      res[100]=roundn(wu[13],3);
      res[101]=roundn(wm[13],3);
      res[102]=roundn(wv[13],3);
      res[103]=roundn(we[13],3);
      res[104]=roundn(wp[13],3);
      res[105]=roundn(wi[13],3);
      res[106]=roundn(ws[13],3);


      //108 109 110 111 112 113 114
       //Vu  Vm  Vv  Ve  Vp  Vi  Vs
      res[107]=roundn(wu[5],3);
      res[108]=roundn(wm[5],3);
      res[109]=roundn(wv[5],3);
      res[110]=roundn(we[5],3);
      res[111]=roundn(wp[5],3);
      res[112]=roundn(wi[5],3);
      res[113]=roundn(wi[5],3);



  return res;
}

} // end of Descriptors namespace
} // end of RDKit namespace
