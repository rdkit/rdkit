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

#include "GETAWAY.h"

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
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>

using namespace Eigen;


double relativeMw1[]={0.084,0,0,0,0.900,1.000,1.166,1.332,1.582,0,0,0, 2.246, 2.339, 2.579,2.670, 2.952,0,0,0,0,0,0,0,0, 4.650, 4.907, 4.887, 5.291, 5.445,0,0,0,0, 6.653,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 9.884,0,0, 10.56};
double relativeVdW1[]={0.299,0,0,0,0.796,1.000,0.695,0.512,0.410,0,0,0, 1.626, 1.424, 1.181,1.088, 1.035,0,0,0,0,0,0,0,0, 1.829, 1.561, 0.764, 0.512, 1.708,0,0,0,0, 1.384,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 2.042,0,0, 1.728};
double relativeNeg1[]={0.944,0,0,0,0.828,1.000,1.163,1.331,1.457,0,0,0, 0.624, 0.779, 0.916,1.077, 1.265,0,0,0,0,0,0,0,0, 0.728, 0.728, 0.728, 0.740, 0.810,0,0,0,0, 1.172,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.837,0,0, 1.012};
double relativePol1[]={0.379,0,0,0,1722,1.000,0.625,0.456,0.316,0,0,0, 3.864, 3.057, 2.063,1.648, 1.239,0,0,0,0,0,0,0,0, 4.773, 4.261, 3.864, 3.466, 4.034,0,0,0,0, 1.733,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 4.375,0,0, 3.040};
double ionpol[]={1.208,0,0.479,0.828,0.737,1,1.291,1.209,1.547,0,0.456,0.679,0.532,0.724,0.931,0.92,1.152,0,0.386,0.543,0,0,0,0.601,0.66,0.702,0.7,0.679,0.686,0.834,0.533,0.702,0.872,0.866,1.049,0,0.371,0.506,0,0,0,0.63,0,0,0,0,0.673,0.799,0.514,0.652,0.767,0.8,0.928,0,0,0,0,0,0,0,0,0,0,0.546,0,0,0,0,0,0,0,0,0,0,0,0,0,0.799,0.819,0.927,0.542,0.659,0.647};
double rcov[]={0.37,0,1.34,0.90,0.82,0.77,0.73,0.71,0,1.54,1.30,1.18,1.11,1.06,1.02,0.99,0,1.96,1.74,0,0,0,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.19,1.16,1.14,0,2.11,1.92,0,0,0,1.45,0,0,0,0,1.53,1.48,1.44,1.41,1.38,1.35,1.33,0,0,0,0,0,0,0,0,0,0,1.79,0,0,0,0,0,0,0,0,0,0,0,0,0,1.28,1.44,1.49,1.48,1.47, 1.46};


namespace RDKit {
namespace Descriptors{

namespace {

double* retreiveMat(MatrixXd matrix) {
   double* arrayd = matrix.data();
   return arrayd;
}

double* retreiveVect(VectorXd matrix) {
   double* arrayd = matrix.data();
   return arrayd;

}

VectorXd getEigenVect(std::vector<double> v){
    
    double* varray_ptr = &v[0];
    Map<VectorXd> V(varray_ptr,v.size());
    return V;
}


/*
VectorXd clusterArray(VectorXd Data, double diff) {

    std::vector<double> data(Data.data(), Data.data()+Data.size());

    // sort the input data
    std::sort(data.begin(), data.end());

    // find the difference between each number and its predecessor
    std::vector<double> diffs;
    std::adjacent_difference(data.begin(), data.end(), std::back_inserter(diffs));

    // convert differences to percentage changes
    std::transform(diffs.begin(), diffs.end(), data.begin(), diffs.begin(),
        std::divides<double>());

    // print out the results
    for (int i = 0; i < data.size(); i++) {

        // if a difference exceeds 40%, start a new group:
        if (diffs[i] > diff)  // diff=0.4 <=> 40%
            std::cout << "\n";

        // print out an item:
        std::cout << data[i] << "\t";
    }
    Map<VectorXd> Res(&data, Data.size());

    return Res;
}
*/


MatrixXd GetPinv(MatrixXd A){
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);

    double  pinvtoler=1.e-6; // choose your tolerance wisely!
    VectorXd vs=svd.singularValues();
    VectorXd vsinv=svd.singularValues();

    for ( long i=0; i<A.cols(); ++i) {
        if ( vs(i) > pinvtoler )
           vsinv(i)=1.0/vs(i);
       else vsinv(i)=0;
    }

    MatrixXd S =  vsinv.asDiagonal();
    MatrixXd Ap = svd.matrixV() * S * svd.matrixU().transpose();
    return Ap;
}


MatrixXd GetCenterMatrix(MatrixXd Mat){

    VectorXd v = Mat.colwise().mean();

    MatrixXd X=Mat.rowwise() - v.transpose();

    return X;
}


MatrixXd GetHmatrix(MatrixXd X){

    MatrixXd  Weigthed = X.transpose()*X;
    return X * GetPinv(Weigthed) * X.transpose();
}


MatrixXd GetRmatrix(MatrixXd H, MatrixXd DM, int numAtoms){

    MatrixXd R = MatrixXd::Zero(numAtoms,numAtoms);
    for (int i=0;i<numAtoms-1;i++) {
      for (int j=i+1;j<numAtoms;j++){
          R(i,j)=sqrt(H(i,i)*H(j,j)) / DM(i,j) ;
          R(j,i)=R(i,j);
      }
    }

    return R;
}

/* for 3Dautocorrelation staff */
double* GetGeodesicMatrix(double* dist, int lag,int numAtoms){
    int sizeArray=numAtoms*numAtoms;
    double *Geodesic = new double[sizeArray];
    for (int i=0; i<sizeArray;i++) {
      if (dist[i]==lag) Geodesic[i]=1;
      else  Geodesic[i]=0;
    }

    return Geodesic;
}

/* for 3Dautocorrelation staff */

std::vector<double> GetRelativeMW(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativeMw1[mol.getAtomWithIdx(i)->getAtomicNum()-1];
  }


  return pol;
}

std::vector<double> GetRelativeRcov(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> wroc(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    wroc[i]=rcov[mol.getAtomWithIdx(i)->getAtomicNum()-1]/rcov[5];
  }


  return wroc;
}



std::vector<double> GetRelativeENeg(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> neg(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    neg[i]=relativeNeg1[mol.getAtomWithIdx(i)->getAtomicNum()-1];
}

  return neg;
}



std::vector<double> GetUn(int numAtoms){

  std::vector<double> u(numAtoms, 1.0);

  return u;
}


std::vector<double> GetRelativeVdW(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> vdw(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    vdw[i]=relativeVdW1[mol.getAtomWithIdx(i)->getAtomicNum()-1];

  }

  return vdw;
}

std::vector<double> GetRelativePol(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=relativePol1[mol.getAtomWithIdx(i)->getAtomicNum()-1];

  }

  return pol;
}

std::vector<double> GetRelativeIonPol(const ROMol& mol){
   int numAtoms= mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0.0);
  for( int i=0; i<numAtoms; ++i){

    pol[i]=ionpol[mol.getAtomWithIdx(i)->getAtomicNum()-1]/ionpol[5];

  }

  return pol;
}



/*
double GetRCON(MatrixXd H, MatrixXd DM, int numAtoms){

    for (int i=0;i<numAtoms-1;i++) {
      for (int j=i+1;j<numAtoms;j++){
          R(i,j)=sqrt(H(i,i)*H(j,j)) / DM(i,j) ;

      }
    }

    return R;
}
*/



double getTDB(int k,int numAtoms,std::vector<double>  w, double* distTopo, double* dist3D){
// similar implementation of J. Chem. Inf. Comput. Sci. 2004, 44, 200-209 equation 1 or 2 page 201
// we use instead of absolute values the relative ones as in Dragon
        double tbd=0.0;

        int numtopoatoms = 0;
            for (int i=0; i<numAtoms-1; ++i)
            {
                for (int j=i+1; j<numAtoms; ++j)
                {
                    if (distTopo[i*numAtoms+j]==k)
                    {
                        tbd += w[i] * dist3D[i*numAtoms+j] * w[j];
                        numtopoatoms++;
                    }
                }
            }

            if (numtopoatoms>0)
                tbd = tbd / numtopoatoms;

  return tbd;
}







JacobiSVD<MatrixXd> getSVD(MatrixXd Mat) {

    JacobiSVD<MatrixXd> svd(Mat,  ComputeThinU | ComputeThinV);

    return svd;
}



double* getGetawayDesc(MatrixXd H, MatrixXd R, int numAtoms) {

    double *w = new double[1];
    // prepare data for Whim parameter computation
    // compute parameters
    
std::cout << H.diagonal() << "\n";

    w[0]=1.0;
    double HGM=1.0;
    for (int i=0;i<numAtoms;i++) {
      std::cout << H(i,i) << ",";
      HGM=HGM*H(i,i);
    }

    HGM=100.0*pow(HGM,1.0/numAtoms);

    double HIC=0.0;
   for (int i=0;i<numAtoms;i++) {
      HIC-=H(i,i)/2.0*log(H(i,i)/2.0)/log(2);
    }    
 
    double RARS=R.rowwise().sum().sum()/numAtoms;
 
    JacobiSVD<MatrixXd> svd = getSVD(R);

    VectorXd EIG = svd.singularValues();


    std::cout << "HGM:"<< HGM << "| HIC:"<< HIC <<"| RARS:"<< RARS << "| REIG:"<< EIG(0) << "\n";

    return w;



}


double* GetGETAWAY(const Conformer &conf, double Vpoints[], MatrixXd DM, MatrixXd AM){

    double *w = new double[18];

    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd H = GetHmatrix(Xmean);

    MatrixXd R = GetRmatrix(H, DM, numAtoms);
    std::cout << "H:" << H << "\n";
    std::cout << "R" << R << "\n";
    std::cout <<  "\n";

    w= getGetawayDesc(H,R, numAtoms);

    return w;
}

} //end of anonymous namespace



double* GETAWAY(const ROMol& mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double Vpoints[3*numAtoms];

  double *dist3D = MolOps::get3DDistanceMat(mol, confId);
  double *dist = MolOps::getDistanceMat(mol, false); // need to be be set to false to have topological distance not weigthed!
  
  double *AdjMat = MolOps::getAdjacencyMatrix(mol,false,0,false,0); // false to have only the 1,0 matrix unweighted

  Map<MatrixXd> am(AdjMat, numAtoms,numAtoms);

  Map<MatrixXd> dm(dist3D, numAtoms,numAtoms);

  Map<MatrixXd> dmtopo(dist, numAtoms,numAtoms);



  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  double *res= new double[1];


   std::vector<double>  wp = GetRelativePol(mol);

   VectorXd Wp = getEigenVect(wp);

   std::vector<double>  wm = GetRelativeMW(mol);

   VectorXd Wm = getEigenVect(wm);

   std::vector<double>  wi = GetRelativeIonPol(mol);

   VectorXd Wi = getEigenVect(wi);

   std::vector<double>  wv = GetRelativeVdW(mol);

   VectorXd Wv = getEigenVect(wv);

   std::vector<double>  we = GetRelativeENeg(mol);

   VectorXd We = getEigenVect(we);

   std::vector<double>  wu = GetUn(numAtoms);

   VectorXd Wu = getEigenVect(wu);

   std::vector<double>  wr = GetRelativeRcov(mol);

   VectorXd Wr = getEigenVect(wr);


   MatrixXd Bi;
   VectorXd tmp;
  for (int i=1;i<11;i++){
    double * Bimat = GetGeodesicMatrix(dist, i, numAtoms);
    Map<MatrixXd> Bi(Bimat, numAtoms,numAtoms);
    MatrixXd RBi=Bi.cwiseProduct(dm);
    double Bicount=(double) Bi.sum();


    tmp=Wu.transpose() * RBi * Wu / Bicount;

    std::cout << "TDBu" << i << ":"<<  tmp << "\n";

    tmp=Wm.transpose() * RBi * Wm / Bicount;

    std::cout << "TDBm" << i << ":"<<  tmp << "\n";

    tmp=Wv.transpose() * RBi * Wv / Bicount;

    std::cout << "TDBv" << i << ":"<<  tmp << "\n";
    
    tmp=We.transpose() * RBi * We / Bicount;

    std::cout << "TDBe" << i << ":"<<  tmp << "\n";

    tmp=Wp.transpose() * RBi * Wp / Bicount;

    std::cout << "TDBp" << i << ":"<<  tmp << "\n";

    tmp=Wi.transpose() * RBi * Wi / Bicount;

    std::cout << "TDBi" << i << ":"<<  tmp << "\n";

    tmp=Wr.transpose() * RBi * Wr / Bicount;

    std::cout << "TDBr" << i << ":"<<  tmp << "\n";


    //double tbd2= getTDB(i,numAtoms, wp,  dist, dist3D);
    //std::cout << "loopTDBp" << i << ":"<< tbd2<< "\n";

    



  }


  double* wus= GetGETAWAY(conf, Vpoints, dm,am);

  return res;
}

} // end of Descriptors namespace
} // end of RDKit namespace
