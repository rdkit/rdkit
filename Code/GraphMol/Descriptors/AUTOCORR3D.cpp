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

   double *w= new double[1];

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

    tmp=Ws.transpose() * RBi * Ws / Bicount;

    std::cout << "TDBs" << i << ":"<<  tmp << "\n";

    tmp=Wr.transpose() * RBi * Wr / Bicount;

    std::cout << "TDBr" << i << ":"<<  tmp << "\n";

    //double tbd2= getTDB(i,numAtoms, wp,  dist, dist3D);
    //std::cout << "loopTDBp" << i << ":"<< tbd2<< "\n";
  }

    return w;

  }


  double* Get3Dauto(double* dist3D,double* dist, int numAtoms, const ROMol& mol){


      double *w1= get3DautocorrelationDesc(dist3D,dist, numAtoms, mol);

      return w1;
  }

} //end of anonymous namespace


double* AUTOCORR3D(const ROMol& mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double *dist = MolOps::getDistanceMat(mol, false);
  double *dist3D = MolOps::get3DDistanceMat(mol, confId);

  double* wpol= Get3Dauto(dist3D,dist, numAtoms, mol);

  return wpol;
}

} // end of Descriptors namespace
} // end of RDKit namespace