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
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>

using namespace Eigen;


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
  
  double *AdjMat = MolOps::getAdjacencyMatrix(mol,false,0,false,0); // false to have only the 1,0 matrix unweighted

  Map<MatrixXd> am(AdjMat, numAtoms,numAtoms);

  Map<MatrixXd> dm(dist3D, numAtoms,numAtoms);

  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  double *res= new double[1];

  double* wu= GetGETAWAY(conf, Vpoints, dm,am);

  return res;
}

} // end of Descriptors namespace
} // end of RDKit namespace
