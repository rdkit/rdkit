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
// Guillaume GODIN access the AutoCorrelation 3D descriptors names in Dragon TDB


#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "AUTOCORR3D.h"
#include "MolData3Ddescriptors.h"

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

MolData3Ddescriptors moldata3D;

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


double* AppendDouble(double *w, double* Append, int length, int pos){
    for (int i=pos;i<pos+length;i++) {
          w[i]=Append[i-pos];
        }

    return w;
}


double* GetGeodesicMatrix(double* dist, int lag,int numAtoms){
    int sizeArray=numAtoms*numAtoms;
    double *Geodesic = new double[sizeArray];
    for (int i=0; i<sizeArray;i++) {
      if (dist[i]==lag) Geodesic[i]=1;
      else  Geodesic[i]=0;
    }

    return Geodesic;
}


double* get3DautocorrelationDesc(double* dist3D, double* dist, int numAtoms, const ROMol& mol) {

   Map<MatrixXd> dm(dist3D, numAtoms,numAtoms);

   double *w= new double[80];

   std::vector<double> wp= moldata3D.GetRelativePol(mol);

   VectorXd Wp = getEigenVect(wp);

   std::vector<double> wm= moldata3D.GetRelativeMW(mol);

   VectorXd Wm = getEigenVect(wm);

   std::vector<double> wi= moldata3D.GetRelativeIonPol(mol);

   VectorXd Wi = getEigenVect(wi);

   std::vector<double> wv= moldata3D.GetRelativeVdW(mol);

   VectorXd Wv = getEigenVect(wv);

   std::vector<double> we= moldata3D.GetRelativeENeg(mol);

   VectorXd We = getEigenVect(we);

   std::vector<double>  wu = moldata3D.GetUn(numAtoms);

   VectorXd Wu = getEigenVect(wu);

   std::vector<double>  ws =  moldata3D.GetIState(mol);

   VectorXd Ws = getEigenVect(ws);

   std::vector<double> wr= moldata3D.GetRelativeRcov(mol);

   VectorXd Wr = getEigenVect(wr);

 
   MatrixXd Bi;
   MatrixXd tmp;
   double TDBmat[8][10];
   double dtmp;

  for (int i=0;i<10;i++){
    double * Bimat = GetGeodesicMatrix(dist, i+1, numAtoms);
    Map<MatrixXd> Bi(Bimat, numAtoms,numAtoms);
    MatrixXd RBi=Bi.cwiseProduct(dm);
    double Bicount=(double) Bi.sum();

    tmp=Wu.transpose() * RBi * Wu / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[0][i]=dtmp;

    tmp=Wm.transpose() * RBi * Wm / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[1][i]=dtmp;

    tmp=Wv.transpose() * RBi * Wv / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[2][i]=dtmp;

    tmp=We.transpose() * RBi * We / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[3][i]=dtmp;

    tmp=Wp.transpose() * RBi * Wp / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[4][i]=dtmp;

    tmp=Wi.transpose() * RBi * Wi / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[5][i]=dtmp;

    tmp=Ws.transpose() * RBi * Ws / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[6][i]=dtmp;

    tmp=Wr.transpose() * RBi * Wr / Bicount;
    dtmp=(double)tmp(0);
    if  (std::isnan(dtmp)) dtmp=0.0;
    TDBmat[7][i]=dtmp;

  }

  // create the Output vector!
    w= AppendDouble(w, TDBmat[0], 10, 0);
    w= AppendDouble(w, TDBmat[1], 10, 10);
    w= AppendDouble(w, TDBmat[2], 10, 20);
    w= AppendDouble(w, TDBmat[3], 10, 30);
    w= AppendDouble(w, TDBmat[4], 10, 40);
    w= AppendDouble(w, TDBmat[5], 10, 50);
    w= AppendDouble(w, TDBmat[6], 10, 60);
    w= AppendDouble(w, TDBmat[7], 10, 70);

    return w;

}

  double* Get3Dauto(double* dist3D,double* dist, int numAtoms, const ROMol& mol){

      std::vector<std::string> AUTOCORRNAMES={"TDB01u","TDB02u","TDB03u","TDB04u","TDB05u","TDB06u","TDB07u","TDB08u","TDB09u","TDB10u","TDB01m","TDB02m","TDB03m","TDB04m","TDB05m","TDB06m","TDB07m","TDB08m","TDB09m","TDB10m","TDB01v","TDB02v","TDB03v","TDB04v","TDB05v","TDB06v","TDB07v","TDB08v","TDB09v","TDB10v","TDB01e","TDB02e","TDB03e","TDB04e","TDB05e","TDB06e","TDB07e","TDB08e","TDB09e","TDB10e","TDB01p","TDB02p","TDB03p","TDB04p","TDB05p","TDB06p","TDB07p","TDB08p","TDB09p","TDB10p","TDB01i","TDB02i","TDB03i","TDB04i","TDB05i","TDB06i","TDB07i","TDB08i","TDB09i","TDB10i","TDB01s","TDB02s","TDB03s","TDB04s","TDB05s","TDB06s","TDB07s","TDB08s","TDB09s","TDB10s","TDB01r","TDB02r","TDB03r","TDB04r","TDB05r","TDB06r","TDB07r","TDB08r","TDB09r","TDB10r"};
      double *res= new double[80];
      res= get3DautocorrelationDesc(dist3D,dist, numAtoms, mol);
      return res;
  }

} //end of anonymous namespace


std::vector<double> AUTOCORR3D(const ROMol& mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double *dist = MolOps::getDistanceMat(mol, false);
  double *dist3D = MolOps::get3DDistanceMat(mol, confId);
  double *wpol= Get3Dauto(dist3D, dist, numAtoms, mol);

  std::vector<double> dataVec;
  for (int i=0;i<80;i++){
    dataVec.push_back(wpol[i]);
  }


  return dataVec;
}

} // end of Descriptors namespace
} // end of RDKit namespace