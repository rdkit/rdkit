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
// Adding WHIM descriptors to 3D descriptors by Guillaume Godin
// for build & set RDBASE! => export RDBASE=/Users/mbp/Github/rdkit_mine/

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "WHIM.h"
#include "MolData3Ddescriptors.h"

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
      gamma[i] =0.0;
      if (ns>0) {
           gamma[i] = ((ns / nAT) * log(ns / nAT) / log(2.0) + (na / nAT) * log(1.0 / nAT) / log(2.0));
           gamma[i] = 1.0 / (1.0 - gamma[i]);
      }
    }
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

    std::vector<double> weigthvector = moldata3D.GetRelativeMW(conf.getOwningMol());

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

    std::vector<double> weigthvector = moldata3D.GetRelativeVdW(conf.getOwningMol());

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

    std::vector<double> weigthvector = moldata3D.GetRelativeENeg(conf.getOwningMol());

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

    std::vector<double> weigthvector = moldata3D.GetRelativePol(conf.getOwningMol());

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

    std::vector<double> weigthvector = moldata3D.GetRelativeIonPol(conf.getOwningMol());

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

    std::vector<double> weigthvector = moldata3D.GetIState(mol);

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




std::vector<double> WHIM(const ROMol& mol,int confId, double th){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double Vpoints[3*numAtoms];

  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  std::vector<double> res;

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
        res.push_back(roundn(wu[map1[i]],3));

     }

   // descnames={"L1m","L2m","L3m","Tm","Am","Vm","P1m","P2m","P3m","Km","E1m","E2m","E3m","Dm","G1m","G2m","G3m","Gm"};

    for (int i=0;i<11;i++) {
        res.push_back(roundn(wm[map1[i]],3));

     }

  //   descnames={"L1v","L2v","L3v","Tv","Av","Vv","P1v","P2v","P3v","Kv","E1v","E2v","E3v","Dv","G1v","G2v","G3v","Gv"};

    for (int i=0;i<11;i++) {
        res.push_back(roundn(wv[map1[i]],3));
     }

   // descnames={"L1e","L2e","L3e","Te","Ae","Ve","P1e","P2e","P3e","Ke","E1e","E2e","E3e","De","G1e","G2e","G3e","Ge"};

    for (int i=0;i<11;i++) {
        res.push_back( roundn(we[map1[i]],3));
     }


   // descnames={"L1p","L2p","L3p","Tp","Ap","Vp","P1p","P2p","P3p","Kp","E1p","E2p","E3p","Dp","G1p","G2p","G3p","Gp"};

    for (int i=0;i<11;i++) {
        res.push_back( roundn(wp[map1[i]],3));

     }
    

    for (int i=0;i<11;i++) {
       res.push_back( roundn(wi[map1[i]],3));

     }


    for (int i=0;i<11;i++) {
        res.push_back(roundn(ws[map1[i]],3));

     }

       //78  79  80  81  82  83  84
       //Tu  Tm  Tv  Te  Tp  Ti  Ts
      res.push_back(roundn(wu[3],3));
      res.push_back(roundn(wm[3],3));
      res.push_back(roundn(wv[3],3));
      res.push_back(roundn(we[3],3));
      res.push_back(roundn(wp[3],3));
      res.push_back(roundn(wi[3],3));
      res.push_back(roundn(ws[3],3));

      //85  86  87  88  89  90  91
      //Au  Am  Av  Ae  Ap  Ai  As
      res.push_back(roundn(wu[4],3));
      res.push_back(roundn(wm[4],3));
      res.push_back(roundn(wv[4],3));
      res.push_back(roundn(we[4],3));
      res.push_back(roundn(wp[4],3));
      res.push_back(roundn(wi[4],3));
      res.push_back(roundn(ws[4],3));

      //92  93
      //Gu  Gm
      res.push_back(roundn(wu[17],3));
      res.push_back(roundn(wm[17],3));

      //94  95  96  97  98  99  100
      //Ku  Km  Kv  Ke  Kp  Ki  Ks
      res.push_back(roundn(wu[9],3));
      res.push_back(roundn(wm[9],3));
      res.push_back(roundn(wv[9],3));
      res.push_back(roundn(we[9],3));
      res.push_back(roundn(wp[9],3));
      res.push_back(roundn(wi[9],3));
      res.push_back(roundn(ws[9],3));

      //101 102 103 104 105 106 107
       //Du  Dm  Dv  De  Dp  Di  Ds
      res.push_back(roundn(wu[13],3));
      res.push_back(roundn(wm[13],3));
      res.push_back(roundn(wv[13],3));
      res.push_back(roundn(we[13],3));
      res.push_back(roundn(wp[13],3));
      res.push_back(roundn(wi[13],3));
      res.push_back(roundn(ws[13],3));


      //108 109 110 111 112 113 114
       //Vu  Vm  Vv  Ve  Vp  Vi  Vs
      res.push_back(roundn(wu[5],3));
      res.push_back(roundn(wm[5],3));
      res.push_back(roundn(wv[5],3));
      res.push_back(roundn(we[5],3));
      res.push_back(roundn(wp[5],3));
      res.push_back(roundn(wi[5],3));
      res.push_back(roundn(wi[5],3));



  return res;
}

} // end of Descriptors namespace
} // end of RDKit namespace
