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

#include "auto3Dcorr.h"

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


double Pol3[] = {
    0.67,  0,    24.3, 5.60,  3.03, 1.76,  1.10, 0.80, 0.56, 0,    23.6, 10.6,
    6.80,  5.38, 3.63, 2.90,  2.18, 0,     43.4, 22.8, 0,    0,    0,    11.60,
    9.40,  8.40, 7.50, 6.80,  6.10, 7.10,  8.12, 6.07, 4.31, 3.77, 3.05, 0,
    47.3,  27.6, 0,    0,     0,    12.80, 0,    0,    0,    0,    7.20, 7.20,
    10.20, 7.70, 6.60, 5.50,  5.35, 0,     0,    0,    0,    0,    0,    0,
    0,     0,    0,    23.50, 0,    0,     0,    0,    0,    0,    0,    0,
    0,     0,    0,    0,     0,    6.50,  5.80, 5.70, 7.60, 6.80, 7.40};
double ElectroNeg3[] = {
    2.59, 0,    0.89, 1.81, 2.28, 2.75, 3.19, 3.65, 4.0,  0,    0.56, 1.32,
    1.71, 2.14, 2.52, 2.96, 3.48, 0,    0.45, 0.95, 0,    0,    0,    1.66,
    2.2,  2.2,  2.56, 1.94, 1.98, 2.23, 2.42, 2.62, 2.82, 3.01, 3.22, 0,
    0.31, 0.72, 0,    0,    0,    1.15, 0,    0,    0,    0,    1.83, 1.98,
    2.14, 2.3,  2.46, 2.62, 2.78, 0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    2.0,  0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    2.28, 2.54, 2.2,  2.25, 2.29, 2.34};

double VdW3[] = {
    6.71,  0,     25.25, 0.0,   17.88, 22.45, 15.6,  11.49, 9.20,   0,     49.0,
    21.69, 36.51, 31.98, 26.52, 24.43, 22.45, 0,     87.11, 0.0,   0,     0,
    0,     44.6,  43.4,  41.05, 35.04, 17.16, 11.49, 11.25, 27.39, 28.73, 26.52,
    28.73, 31.06, 0,     0.0,   0.0,   0,     0,     0,     33.51, 0,     0,
    0,     0,     21.31, 16.52, 30.11, 45.83, 38.79, 36.62, 32.52, 0,     0,
    0,     0,     0,     0,     0,     0,     0,     0,     72.78, 0,     0,
    0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    22.45, 19.16, 15.6,  31.54, 34.53, 38.79};

double IonPolarizability3[] = {
  13.598, 0, 5.392, 9.323, 8.298, 11.26, 14.534, 13.618, 17.423, 0, 5.139, 7.646,
  5.986, 8.152, 10.487, 10.360, 12.968, 0, 4.341, 6.113, 0, 0, 0, 6.767, 7.434,
  7.902, 7.881, 7.640, 7.723, 9.394, 5.999, 7.900, 9.815, 9.752, 11.814, 0, 4.177,
  5.695, 0, 0, 0, 7.092, 0, 0, 0, 0, 7.576, 8.994, 5.786, 7.344, 8.640, 9.010,
  10.451, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 9.000, 9.226, 10.438, 6.108, 7.417, 7.289};



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

VectorXd getEigenVect(std::vector<double> v){
    
    double* varray_ptr = &v[0];
    Map<VectorXd> V(varray_ptr,v.size());
    return V;
}



std::vector<double> GetRelativePol(const ROMol &mol) {
  int numAtoms = mol.getNumAtoms();

  std::vector<double> pol(numAtoms, 0);
  for (int i = 0; i < numAtoms; ++i) {
    pol[i] = Pol3[mol.getAtomWithIdx(i)->getAtomicNum()-1] / Pol3[5];
  }

  return pol;
}


MatrixXd GetGeodesicMatrix(double* dist, int lag,int numAtoms){

    double *Geodesic = new double(numAtoms*numAtoms);
    for (int i=0; i<dist.size();i++) {
      if (dist[i]==lag) Geodesic[i]=1;
    }
    Map<MatrixXd> R(Geodesic, numAtoms,numAtoms);

    return R;
}


double* get3DautocorrelationDesc(MatrixXd DM3D,MatrixXd DM, int numAtoms) {

  double *w = new double[1];
    
   std::vector<double>  wp = GetRelativePol(mol);

   VectorXd W = getEigenVect(wp);

   double* TDB = new double(10);

  for (int i=1;i<11;i++){
    MatrixXd Bi= GetGeodesicMatrix(dist, i, numAtoms);
    TDB[i]= W.transpose() * Bi * W;
    std::cout << "TDB" << i << ":"<< TDB[i]*0.5 << "\n";
  }

  return w;

}


double* Get3Dauto(const Conformer &conf, MatrixXd DM3D, MatrixXd DM){

    double *w = new double[18];

    double *w1= get3DautocorrelationDesc(DM3D, DM, numAtoms);

    return w;
}

} //end of anonymous namespace



double* AUTO3DCORR(const ROMol& mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double *dist3D = MolOps::get3DDistanceMat(mol, confId);
  
  double *dist = MolOps::getDistanceMat(mol, false);

  Map<MatrixXd> dm(dist3D, numAtoms,numAtoms);

  double* wpol= Get3Dauto(conf, dm,dist);

  return wpol;
}

} // end of Descriptors namespace
} // end of RDKit namespace
