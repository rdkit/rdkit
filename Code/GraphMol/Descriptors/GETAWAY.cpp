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
#include "Data3Ddescriptors.h"
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

/*
Data3Ddescriptors data3D;

double* relativeMw1=data3D.getMW();
double* relativeVdW1=data3D.getVDW();
double* relativeNeg1=data3D.getNEG();
double* relativePol1=data3D.getPOL();
double* ionpol=data3D.getIonPOL();
double* rcov=data3D.getRCOV();
*/

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
    int degree = atom->getDegree();
    if (degree>0 and atNum>1) {
      int h = atom->getTotalNumHs();
      int Zv = RDKit::PeriodicTable::getTable()->getNouterElecs(atNum);
      double dv =(double) Zv-h;
      dv = dv / (double) (atNum-Zv-1);
      int N = GetPrincipalQuantumNumber(atNum);
      Is.push_back(round(1000*(4.0/(N*N)*dv+1.0)/degree)/1000);  // WHIM-P5.pdf paper 1997  => +7 & NoHydrogens is used!
    }
    else Is.push_back(0);
 }


 return Is;
}



std::vector<double>  clusterArray(std::vector<double> data) {
    std::vector<double> Store;

    // sort the input data
    std::sort(data.begin(), data.end());

    // find the difference between each number and its predecessor
    std::vector<double> diffs;

    std::adjacent_difference(data.begin(), data.end(), std::back_inserter(diffs));

    // convert differences to percentage changes
    std::transform(diffs.begin(), diffs.end(), data.begin(), diffs.begin(),
        std::divides<double>());

    // print out the results
    int j=0;
    int count=0;
    for (int i = 0; i < data.size(); i++) {
        count++;
        // if a difference exceeds 40%, start a new group:
        if (diffs[i] > 0.001)  {// diff=0.4 <=> 40%
            Store.push_back(count);
            //std::cout << "\n" << Store[j];
            count=0;
            j++;
            //std::cout << "\n";
        }
        // print out an item:
        
       //std::cout << data[i] << "\t";
    }
    
    //std::cout << "\n";


    return Store;
}


double GetLevClassNumber(VectorXd LevHeavy, int numHeavyAtoms, std::vector<int> HeavyList){

  int ClassNum=numHeavyAtoms;
  for (int i=0;i<numHeavyAtoms-1;i++){
    for (int j=i+1;j<numHeavyAtoms;j++){
        if (std::abs(LevHeavy[i]-LevHeavy[j])<0.01) {
          ClassNum--;
          break;  // found one pair so move to next i!
        }
    }
  }
  return ClassNum;
}



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


MatrixXd GetRowwiseProdMatVect(MatrixXd A, VectorXd V, int numAtoms){

    std::cout << "A:" << A << "\n";
    std::cout << "V:" << V << "\n";

    MatrixXd R;

    R = A.array().colwise() * V.array();
    //R= A.cwiseProduct( V.replicate( 1, A.rows()).transpose());

    std::cout << "A.*V:" << R << "\n";

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




std::vector<double> GetUn(int numAtoms){

  std::vector<double> u(numAtoms, 1.0);

  return u;
}





int*  GetHeavyList(const ROMol& mol){
  int numAtoms= mol.getNumAtoms();

  int* HeavyList = new int[numAtoms];

  for (int i = 0; i < numAtoms; ++i) {
    const RDKit::Atom * atom= mol.getAtomWithIdx(i);
    int atNum=atom->getAtomicNum();
    if ( atNum>1)  HeavyList[i]=1;
    else  HeavyList[i]=0;
       
  }
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
// we use instead of atomic absolute values the atomic relative ones as in Dragon
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



double getRCON(MatrixXd R, MatrixXd Adj,int numAtoms){
// similar implementation of J. Chem. Inf. Comput. Sci. 2004, 44, 200-209 equation 1 or 2 page 201
// we use instead of atomic absolute values the atomic relative ones as in Dragon
        double RCON=0.0;
        VectorXd VSR = R.rowwise().sum();
        std::cout << "VSR:" << VSR << "\n";
        for (int i=0; i<numAtoms; ++i)
            {
                for (int j=0; j<numAtoms; ++j)
                {
                    if (Adj(i,j)>0)
                    {
                        RCON += VSR(i)*VSR(j);
                    }
                }
            }
  return sqrt(RCON);
}





JacobiSVD<MatrixXd> getSVD(MatrixXd Mat) {

    JacobiSVD<MatrixXd> svd(Mat,  ComputeThinU | ComputeThinV);

    return svd;
}



double* getGetawayDesc(MatrixXd H, MatrixXd R, MatrixXd Adj, int numAtoms,   int* Heavylist,const ROMol& mol) {

    double *w = new double[1];
    // prepare data for Whim parameter computation
    // compute parameters

    VectorXd Lev=H.diagonal();
    std::vector<double> heavyLev;
    for (int i=0;i<numAtoms;i++){
      if (Heavylist[i]>0)
        heavyLev.push_back(Lev(i));
    }

    std::vector<double> Clus= clusterArray(heavyLev);
    double numHeavy=heavyLev.size();

    std::cout << "HeavyAtoms:" << numHeavy << "\n";

    double ITH0 = numHeavy*log(numHeavy)/log(2);
    double ITH=ITH0;
    for (int j=0;j<Clus.size();j++){
      ITH -= Clus[j]*log(Clus[j])/log(2);
    }

    double ISH=ITH/ITH0;

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

    double rcon= getRCON(R,  Adj, numAtoms);

    std::cout <<  "ISH:"<< ISH << "| ITH:"<< ITH <<  " |HGM:"<< HGM << "| HIC:"<< HIC <<"| RARS:"<< RARS << "| REIG:"<< EIG(0) << "| RCON:"<< rcon <<"\n";


   std::vector<double> wp= moldata3D.GetRelativePol(mol);

   VectorXd Wp = getEigenVect(wp);


   std::vector<double> wm= moldata3D.GetRelativeMW(mol);

   //std::vector<double>  wm = GetRelativeMW(mol);



   VectorXd Wm = getEigenVect(wm);

   std::vector<double> wi= moldata3D.GetRelativeIonPol(mol);
   //std::vector<double>  wi = GetRelativeIonPol(mol);

   VectorXd Wi = getEigenVect(wi);

   std::vector<double> wv= moldata3D.GetRelativeVdW(mol);
   //std::vector<double>  wv = GetRelativeVdW(mol);

   VectorXd Wv = getEigenVect(wv);
   
   std::vector<double> we= moldata3D.GetRelativeENeg(mol);

   //std::vector<double>  we = GetRelativeENeg(mol);

   VectorXd We = getEigenVect(we);

   std::vector<double>  wu = GetUn(numAtoms);

   VectorXd Wu = getEigenVect(wu);

   std::vector<double>  ws = GetIState(mol);

   VectorXd Ws = getEigenVect(ws);

   std::vector<double> wr= moldata3D.GetRelativeRcov(mol);

   //std::vector<double>  wr = GetRelativeRcov(mol);

   VectorXd Wr = getEigenVect(wr);


// HATS definition is  (hi*wi) * (hj*wj) using topological distance
// H (0,k) definition is hij*wj*wi using topological distance
// HT = H0+ 2*sum(Hk)
// R like H
// RT like HT

   MatrixXd Bi;
   MatrixXd tmp;
   MatrixXd RBw ;
  double HATS;
  double * HATSk= new double[9];
  double H0;
  double * Hk= new double[9];
  

   double *dist = MolOps::getDistanceMat(mol, false); // need to be be set to false to have topological distance not weigthed!



  for (int i=0;i<9;i++){
    if (i==0) {
      Bi = H.diagonal().asDiagonal();
    }
    
      double* Bimat = GetGeodesicMatrix(dist, i, numAtoms);
      Map<MatrixXd> Bj(Bimat, numAtoms,numAtoms);



    MatrixXd R;
  HATS =0.0;
  H0=0.0;
if (i==0) {
    for (int j=0;j<numAtoms;j++){
      for (int k=j;k<numAtoms;k++){
        if (Bi(j,k)>0){
              HATS+=(Wu(j)*Wu(j)*H(j,j)*H(j,j));
              if (H(j,k)>0)
              H0+=Wu(j)*Wu(k)*H(j,k);

        }
      }
    }
  }


if (i>0) {
    for (int j=0;j<numAtoms;j++){
      for (int k=j;k<numAtoms;k++){
        if (Bj(j,k)==1){
            HATS+=Wu(j)*H(j,j)*Wu(k)*H(k,k);
            if (H(j,k)>0)
              H0+=Wu(j)*Wu(k)*H(j,k);

        }
      }
    }
  }
  HATSk[i]=HATS;
  Hk[i]=H0;


    std::cout << "HATSu" << i << ":"<<  HATS << "\n";
    std::cout << "Hu" << i << ":"<<  H0 << "\n";

/*    RBw = GetRowwiseProdMatVect(Bi,  Wm, numAtoms);
    tmp=RBw.transpose() * RBw;

    std::cout << "HATSm" << i << ":"<<  tmp << "\n";

    RBw = GetRowwiseProdMatVect(Bi,  Wv, numAtoms);
    tmp=RBw.transpose() * RBw;

    std::cout << "HATSv" << i << ":"<<  tmp << "\n";

    RBw = GetRowwiseProdMatVect(Bi,  We, numAtoms);
    tmp=RBw.transpose() * RBw;

    std::cout << "HATSe" << i << ":"<<  tmp << "\n";

    RBw = GetRowwiseProdMatVect(Bi,  We, numAtoms);
    tmp=RBw.transpose() * RBw;

    std::cout << "HATSp" << i << ":"<<  tmp << "\n";

    RBw = GetRowwiseProdMatVect(Bi,  Wi, numAtoms);
    tmp=RBw.transpose() * RBw;

    std::cout << "HATSi" << i << ":"<<  tmp << "\n";


    RBw = GetRowwiseProdMatVect(Bi,  Ws, numAtoms);
    tmp=RBw.transpose() * RBw;

    std::cout << "HATSs" << i << ":"<<  tmp << "\n";

    RBw = GetRowwiseProdMatVect(Bi,  Wr, numAtoms);
    tmp=RBw.transpose() * RBw;


    std::cout << "HATSr" << i << ":"<<  tmp << "\n";
    */
  }
    //double tbd2= getTDB(i,numAtoms, wp,  dist, dist3D);
    //std::cout << "loopTDBp" << i << ":"<< tbd2<< "\n";
  

  double HATST = HATSk[0]+2*(HATSk[1]+HATSk[2]+HATSk[3]+HATSk[4]+HATSk[5]+HATSk[6]+HATSk[7]+HATSk[8]);
    std::cout << "HATStotla:"<<  HATST << "\n";

  double HT = Hk[0]+2*(Hk[1]+Hk[2]+Hk[3]+Hk[4]+Hk[5]+Hk[6]+Hk[7]+Hk[8]);
    std::cout << "Htotla:"<<  HT << "\n";

  return w;



}


double* GetGETAWAY(const Conformer &conf, double Vpoints[], MatrixXd DM, MatrixXd ADJ,  int* Heavylist){

    double *w = new double[18];

    int numAtoms = conf.getNumAtoms();
    const ROMol mol= conf.getOwningMol();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd H = GetHmatrix(Xmean);

    MatrixXd R = GetRmatrix(H, DM, numAtoms);
    std::cout << "H:" << H << "\n";
    std::cout << "R" << R << "\n";
    std::cout <<  "\n";

    w= getGetawayDesc(H,R, ADJ, numAtoms, Heavylist, mol);

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

  Map<MatrixXd> adj(AdjMat, numAtoms,numAtoms);

  std::cout << adj << "\n";
  
  int* Heavylist= GetHeavyList(mol);
  int nHeavyAt= mol.getNumHeavyAtoms(); // should be the same as the List size upper!

  
  Map<MatrixXd> dm(dist3D, numAtoms,numAtoms);

  Map<MatrixXd> dmtopo(dist, numAtoms,numAtoms);



  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  double *res= new double[1];

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

   std::vector<double>  wu = GetUn(numAtoms);

   VectorXd Wu = getEigenVect(wu);

   std::vector<double>  ws = GetIState(mol);

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

// HATS like  0 = leverage >0  => autocorr
// 
 

  

  double* wus= GetGETAWAY(conf, Vpoints, dm,adj, Heavylist);

  return res;
}

} // end of Descriptors namespace
} // end of RDKit namespace
