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

#include "GETAWAY.h"
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


// code from web 
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

    int j=0;
    int count=0;
    for (unsigned int i = 0; i < data.size(); i++) {
        count++;
        // if a difference exceeds 1%, start a new group:
        if (diffs[i] > 0.001)  {// diff=0.01 <=> 1%
            Store.push_back(count);
            count=0;
            j++;
        }
    }

    return Store;
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


MatrixXd GetPinv(MatrixXd A){
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    double  pinvtoler=1.e-2; // choose your tolerance wisely!
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


double* AppendDouble(double *w, double* Append, int length, int pos){
    for (int i=pos;i<pos+length;i++) {
          w[i]=Append[i-pos];
        }

    return w;
}


double getRCON(MatrixXd R, MatrixXd Adj,int numAtoms){
// similar implementation of J. Chem. Inf. Comput. Sci. 2004, 44, 200-209 equation 1 or 2 page 201
// we use instead of atomic absolute values the atomic relative ones as in Dragon
        double RCON=0.0;
        VectorXd VSR = R.rowwise().sum();
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


double getHATS(double W1, double W2, double H1, double H2){
    return W1*H1*W2*H2;

}

double getH(double W1, double W2, double H){
    return W1*H*W2;

}

double getHtotal(double * Hk){

  return Hk[0]+2*(Hk[1]+Hk[2]+Hk[3]+Hk[4]+Hk[5]+Hk[6]+Hk[7]+Hk[8]);

}

double getRtotal(double * Rk){

  return 2*(Rk[0]+Rk[1]+Rk[2]+Rk[3]+Rk[4]+Rk[5]+Rk[6]+Rk[7]);

}


double getMax(double * Rk){
 double RTp=0;
 for (int j=0;j<8;j++){
  if (Rk[j]>RTp) RTp=Rk[j];
 }
 return RTp;

}



JacobiSVD<MatrixXd> getSVD(MatrixXd Mat) {

    JacobiSVD<MatrixXd> svd(Mat,  ComputeThinU | ComputeThinV);

    return svd;
}



std::vector<double> getGetawayDesc(MatrixXd H, MatrixXd R, MatrixXd Adj, int numAtoms,   int* Heavylist,const ROMol& mol) {

    //double *w = new double[273];

    std::vector<double> res;
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


    double ITH0 = numHeavy*log(numHeavy)/log(2);
    double ITH=ITH0;
    for (unsigned int j=0;j<Clus.size();j++){
      ITH -= Clus[j]*log(Clus[j])/log(2);
    }
    res.push_back(ITH);
    double ISH=ITH/ITH0;
    res.push_back(ISH);

    double HIC=0.0;
    for (int i=0;i<numAtoms;i++) {
      HIC-=H(i,i)/2.0*log(H(i,i)/2.0)/log(2);
    }
    res.push_back(HIC);

    double HGM=1.0;
    for (int i=0;i<numAtoms;i++) {
      HGM=HGM*H(i,i);
    }
    HGM=100.0*pow(HGM,1.0/numAtoms);
    res.push_back(HGM);


    double RARS=R.rowwise().sum().sum()/numAtoms;

    JacobiSVD<MatrixXd> svd = getSVD(R);

    VectorXd EIG = svd.singularValues();

    double rcon= getRCON(R,  Adj, numAtoms);

 

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



  MatrixXd Bi;
  MatrixXd RBw;
  double HATSu,HATSm,HATSv,HATSe,HATSp,HATSi,HATSs;
  double H0u,H0m,H0v,H0e,H0p,H0i,H0s;
  double R0u,R0m,R0v,R0e,R0p,R0i,R0s;
  double Rkmaxu,Rkmaxm,Rkmaxv,Rkmaxe,Rkmaxp,Rkmaxi,Rkmaxs;
  double tmpu,tmpm,tmpv,tmpe,tmpp,tmpi,tmps;
  double HATSk[7][9];
  double Hk[7][9];
  double Rk[7][8];
  double Rp[7][8];
  

   double *dist = MolOps::getDistanceMat(mol, false); // need to be be set to false to have topological distance not weigthed!

  for (int i=0;i<9;i++){
    if (i==0) {
      Bi = H.diagonal().asDiagonal();
    }
      double* Bimat = GetGeodesicMatrix(dist, i, numAtoms);
      Map<MatrixXd> Bj(Bimat, numAtoms,numAtoms);



      HATSu =0.0;
      HATSm =0.0;
      HATSv =0.0;
      HATSe =0.0;
      HATSp =0.0;
      HATSi =0.0;
      HATSs =0.0;
      
      H0u=0.0;
      H0m=0.0;
      H0v=0.0;
      H0e=0.0;
      H0p=0.0;
      H0i=0.0;
      H0s=0.0;

      if (i==0) {
          for (int j=0;j<numAtoms;j++){
            for (int k=j;k<numAtoms;k++){
              if (Bi(j,k)>0){
                    HATSu+=getHATS((double)Wu(j), (double)Wu(j), (double)H(j,j), (double)H(j,j));
                    HATSm+=getHATS((double)Wm(j), (double)Wm(j), (double)H(j,j), (double)H(j,j));
                    HATSv+=getHATS((double)Wv(j), (double)Wv(j), (double)H(j,j), (double)H(j,j));
                    HATSe+=getHATS((double)We(j), (double)We(j), (double)H(j,j), (double)H(j,j));
                    HATSp+=getHATS((double)Wp(j), (double)Wp(j), (double)H(j,j), (double)H(j,j));
                    HATSi+=getHATS((double)Wi(j), (double)Wi(j), (double)H(j,j), (double)H(j,j));
                    HATSs+=getHATS((double)Ws(j), (double)Ws(j), (double)H(j,j), (double)H(j,j));

                    if (H(j,k)>0) { 
                        H0u+=getH((double)Wu(j), (double)Wu(k), (double)H(j,k));
                        H0m+=getH((double)Wm(j), (double)Wm(k), (double)H(j,k));
                        H0v+=getH((double)Wv(j), (double)Wv(k), (double)H(j,k));
                        H0e+=getH((double)We(j), (double)We(k), (double)H(j,k));
                        H0p+=getH((double)Wp(j), (double)Wp(k), (double)H(j,k));
                        H0i+=getH((double)Wi(j), (double)Wi(k), (double)H(j,k));
                        H0s+=getH((double)Ws(j), (double)Ws(k), (double)H(j,k));
                    }
              }
            }
          }
        }


      if (i>0) {
          for (int j=0;j<numAtoms;j++){
            for (int k=j;k<numAtoms;k++){
              if (Bj(j,k)==1){
                    HATSu+=getHATS((double)Wu(j), (double)Wu(k), (double)H(j,j), (double)H(k,k));
                    HATSm+=getHATS((double)Wm(j), (double)Wm(k), (double)H(j,j), (double)H(k,k));
                    HATSv+=getHATS((double)Wv(j), (double)Wv(k), (double)H(j,j), (double)H(k,k));
                    HATSe+=getHATS((double)We(j), (double)We(k), (double)H(j,j), (double)H(k,k));
                    HATSp+=getHATS((double)Wp(j), (double)Wp(k), (double)H(j,j), (double)H(k,k));
                    HATSi+=getHATS((double)Wi(j), (double)Wi(k), (double)H(j,j), (double)H(k,k));
                    HATSs+=getHATS((double)Ws(j), (double)Ws(k), (double)H(j,j), (double)H(k,k));



                  if (H(j,k)>0) {
                        H0u+=getH((double)Wu(j), (double)Wu(k), (double)H(j,k));
                        H0m+=getH((double)Wm(j), (double)Wm(k), (double)H(j,k));
                        H0v+=getH((double)Wv(j), (double)Wv(k), (double)H(j,k));
                        H0e+=getH((double)We(j), (double)We(k), (double)H(j,k));
                        H0p+=getH((double)Wp(j), (double)Wp(k), (double)H(j,k));
                        H0i+=getH((double)Wi(j), (double)Wi(k), (double)H(j,k));
                        H0s+=getH((double)Ws(j), (double)Ws(k), (double)H(j,k));
                  }
              }
            }
          }
        }


        HATSk[0][i]=HATSu;
        HATSk[1][i]=HATSm;
        HATSk[2][i]=HATSv;
        HATSk[3][i]=HATSe;
        HATSk[4][i]=HATSp;
        HATSk[5][i]=HATSi;
        HATSk[6][i]=HATSs;

        Hk[0][i]=H0u;
        Hk[1][i]=H0m;
        Hk[2][i]=H0v;
        Hk[3][i]=H0e;
        Hk[4][i]=H0p;
        Hk[5][i]=H0i;
        Hk[6][i]=H0s;

        R0u=0.0;
        R0m=0.0;
        R0v=0.0;
        R0e=0.0;
        R0p=0.0;
        R0i=0.0;
        R0s=0.0;

        Rkmaxu=0;
        Rkmaxm=0;
        Rkmaxv=0;
        Rkmaxe=0;
        Rkmaxp=0;
        Rkmaxi=0;
        Rkmaxs=0;

       if (i>0) {
       for (int j=0;j<numAtoms-1;j++){
            for (int k=j+1;k<numAtoms;k++){
              if (Bj(j,k)==1){
                  tmpu = getH((double)Wu(j), (double)Wu(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  tmpm = getH((double)Wm(j), (double)Wm(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  tmpv = getH((double)Wv(j), (double)Wv(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  tmpe = getH((double)We(j), (double)We(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  tmpp = getH((double)Wp(j), (double)Wp(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  tmpi = getH((double)Wi(j), (double)Wi(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  tmps = getH((double)Ws(j), (double)Ws(k), (double)R(j,k)); // Use same function but on all R not "H>0" like in the previous loop & i>0!
                  R0u+=tmpu;
                  R0m+=tmpm;
                  R0v+=tmpv;
                  R0e+=tmpe;
                  R0p+=tmpp;
                  R0i+=tmpi;
                  R0s+=tmps;
                  if (tmpu>Rkmaxu) {
                    Rkmaxu=tmpu;
                  }
                  if (tmpm>Rkmaxm) {
                    Rkmaxm=tmpm;
                  }
                  if (tmpv>Rkmaxv) {
                    Rkmaxv=tmpv;
                  }
                  if (tmpe>Rkmaxe) {
                    Rkmaxe=tmpe;
                  }
                  if (tmpp>Rkmaxp) {
                    Rkmaxp=tmpp;
                  }
                  if (tmpi>Rkmaxi) {
                    Rkmaxi=tmpi;
                  }
                  if (tmps>Rkmaxs) {
                    Rkmaxs=tmps;
                  }



              }
            }
          }
          Rk[0][i-1]=R0u;
          Rk[1][i-1]=R0m;
          Rk[2][i-1]=R0v;
          Rk[3][i-1]=R0e;
          Rk[4][i-1]=R0p;
          Rk[5][i-1]=R0i;
          Rk[6][i-1]=R0s;

          Rp[0][i-1]=Rkmaxu;
          Rp[1][i-1]=Rkmaxm;
          Rp[2][i-1]=Rkmaxv;
          Rp[3][i-1]=Rkmaxe;
          Rp[4][i-1]=Rkmaxp;
          Rp[5][i-1]=Rkmaxi;
          Rp[6][i-1]=Rkmaxs;

        }
  }

  //HATSk[0]+2*(HATSk[1]+HATSk[2]+HATSk[3]+HATSk[4]+HATSk[5]+HATSk[6]+HATSk[7]+HATSk[8]);
  double HATSTu =   getHtotal(HATSk[0]);
  double HATSTm =   getHtotal(HATSk[1]);
  double HATSTv =   getHtotal(HATSk[2]);
  double HATSTe =   getHtotal(HATSk[3]);
  double HATSTp =   getHtotal(HATSk[4]);
  double HATSTi =   getHtotal(HATSk[5]);
  double HATSTs =   getHtotal(HATSk[6]);


  // Hk[0]+2*(Hk[1]+Hk[2]+Hk[3]+Hk[4]+Hk[5]+Hk[6]+Hk[7]+Hk[8]);
  double HTu = getHtotal(Hk[0]);
  double HTm = getHtotal(Hk[1]);
  double HTv = getHtotal(Hk[2]);
  double HTe = getHtotal(Hk[3]);
  double HTp = getHtotal(Hk[4]);
  double HTi = getHtotal(Hk[5]);
  double HTs = getHtotal(Hk[6]);

  //
  //2*(Rk[1]+Rk[2]+Rk[3]+Rk[4]+Rk[5]+Rk[6]+Rk[7]+Rk[8]);
  double RTu = getRtotal(Rk[0]);
  double RTm = getRtotal(Rk[1]);
  double RTv = getRtotal(Rk[2]);
  double RTe = getRtotal(Rk[3]);
  double RTp = getRtotal(Rk[4]);
  double RTi = getRtotal(Rk[5]);
  double RTs = getRtotal(Rk[6]);


 double RTMu=getMax(Rp[0]);
 double RTMm=getMax(Rp[1]);
 double RTMv=getMax(Rp[2]);
 double RTMe=getMax(Rp[3]);
 double RTMp=getMax(Rp[4]);
 double RTMi=getMax(Rp[5]);
 double RTMs=getMax(Rp[6]);

// create the output vector...
for (int i=0;i<9;i++){
  res.push_back(Hk[0][i]);
 }
res.push_back(HTu);

for (int i=0;i<9;i++){
  res.push_back(HATSk[0][i]);
 }
res.push_back(HATSTu);

for (int i=0;i<9;i++){
  res.push_back(Hk[1][i]);
 }
res.push_back(HTm);

for (int i=0;i<9;i++){
  res.push_back(HATSk[1][i]);
 }
res.push_back(HATSTm);

for (int i=0;i<9;i++){
  res.push_back(Hk[2][i]);
 }
res.push_back(HTv);

for (int i=0;i<9;i++){
  res.push_back(HATSk[2][i]);
 }
res.push_back(HATSTv);

for (int i=0;i<9;i++){
  res.push_back(Hk[3][i]);
 }
res.push_back(HTe);

for (int i=0;i<9;i++){
  res.push_back(HATSk[3][i]);
 }
res.push_back(HATSTe);

for (int i=0;i<9;i++){
  res.push_back(Hk[4][i]);
 }
res.push_back(HTp);

for (int i=0;i<9;i++){
  res.push_back(HATSk[4][i]);
 }
res.push_back(HATSTp);

for (int i=0;i<9;i++){
  res.push_back(Hk[5][i]);
 }
res.push_back(HTi);

for (int i=0;i<9;i++){
  res.push_back(HATSk[5][i]);
 }
res.push_back(HATSTi);

for (int i=0;i<9;i++){
  res.push_back(Hk[6][i]);
 }
res.push_back(HTs);

for (int i=0;i<9;i++){
  res.push_back(HATSk[6][i]);
 }
res.push_back(HATSTs);

res.push_back(rcon); // this is not the same as in Dragon
res.push_back(RARS);
res.push_back(EIG(0));

 for (int i=0;i<8;i++){
  res.push_back(Rk[0][i]);
 }
res.push_back(RTu);

for (int i=0;i<8;i++){
  res.push_back(Rp[0][i]);
 }
res.push_back(RTMu);

 for (int i=0;i<8;i++){
  res.push_back(Rk[1][i]);
 }
res.push_back(RTm);

for (int i=0;i<8;i++){
  res.push_back(Rp[1][i]);
 }
res.push_back(RTMm);

 for (int i=0;i<8;i++){
  res.push_back(Rk[2][i]);
 }
res.push_back(RTv);

for (int i=0;i<8;i++){
  res.push_back(Rp[2][i]);
 }
res.push_back(RTMv);

 for (int i=0;i<8;i++){
  res.push_back(Rk[3][i]);
 }
res.push_back(RTe);

for (int i=0;i<8;i++){
  res.push_back(Rp[3][i]);
 }
res.push_back(RTMe);

 for (int i=0;i<8;i++){
  res.push_back(Rk[4][i]);
 }
res.push_back(RTp);

for (int i=0;i<8;i++){
  res.push_back(Rp[4][i]);
 }
res.push_back(RTMp);

 for (int i=0;i<8;i++){
  res.push_back(Rk[5][i]);
 }
res.push_back(RTi);

for (int i=0;i<8;i++){
  res.push_back(Rp[5][i]);
 }
res.push_back(RTMi);

// wrong vs Dragon
 for (int i=0;i<8;i++){
  res.push_back(Rk[6][i]);
 }
 // wrong vs Dragon
res.push_back(RTs);
// wrong vs Dragon
for (int i=0;i<8;i++){
  res.push_back(Rp[6][i]);
 }
 // wrong vs Dragon
res.push_back(RTMs);



  return res;
}


std::vector<double> GetGETAWAY(const Conformer &conf, double Vpoints[], MatrixXd DM, MatrixXd ADJ,  int* Heavylist){

    std::vector<std::string> GETAWAYNAMES={"ITH","ISH","HIC","HGM","H0u","H1u","H2u","H3u","H4u","H5u","H6u","H7u","H8u","HTu","HATS0u","HATS1u","HATS2u","HATS3u","HATS4u","HATS5u","HATS6u","HATS7u","HATS8u","HATSu","H0m","H1m","H2m","H3m","H4m","H5m","H6m","H7m","H8m","HTm","HATS0m","HATS1m","HATS2m","HATS3m","HATS4m","HATS5m","HATS6m","HATS7m","HATS8m","HATSm","H0v","H1v","H2v","H3v","H4v","H5v","H6v","H7v","H8v","HTv","HATS0v","HATS1v","HATS2v","HATS3v","HATS4v","HATS5v","HATS6v","HATS7v","HATS8v","HATSv","H0e","H1e","H2e","H3e","H4e","H5e","H6e","H7e","H8e","HTe","HATS0e","HATS1e","HATS2e","HATS3e","HATS4e","HATS5e","HATS6e","HATS7e","HATS8e","HATSe","H0p","H1p","H2p","H3p","H4p","H5p","H6p","H7p","H8p","HTp","HATS0p","HATS1p","HATS2p","HATS3p","HATS4p","HATS5p","HATS6p","HATS7p","HATS8p","HATSp","H0i","H1i","H2i","H3i","H4i","H5i","H6i","H7i","H8i","HTi","HATS0i","HATS1i","HATS2i","HATS3i","HATS4i","HATS5i","HATS6i","HATS7i","HATS8i","HATSi","H0s","H1s","H2s","H3s","H4s","H5s","H6s","H7s","H8s","HTs","HATS0s","HATS1s","HATS2s","HATS3s","HATS4s","HATS5s","HATS6s","HATS7s","HATS8s","HATSs","RCON","RARS","REIG","R1u","R2u","R3u","R4u","R5u","R6u","R7u","R8u","RTu","R1u+","R2u+","R3u+","R4u+","R5u+","R6u+","R7u+","R8u+","RTu+","R1m","R2m","R3m","R4m","R5m","R6m","R7m","R8m","RTm","R1m+","R2m+","R3m+","R4m+","R5m+","R6m+","R7m+","R8m+","RTm+","R1v","R2v","R3v","R4v","R5v","R6v","R7v","R8v","RTv","R1v+","R2v+","R3v+","R4v+","R5v+","R6v+","R7v+","R8v+","RTv+","R1e","R2e","R3e","R4e","R5e","R6e","R7e","R8e","RTe","R1e+","R2e+","R3e+","R4e+","R5e+","R6e+","R7e+","R8e+","RTe+","R1p","R2p","R3p","R4p","R5p","R6p","R7p","R8p","RTp","R1p+","R2p+","R3p+","R4p+","R5p+","R6p+","R7p+","R8p+","RTp+","R1i","R2i","R3i","R4i","R5i","R6i","R7i","R8i","RTi","R1i+","R2i+","R3i+","R4i+","R5i+","R6i+","R7i+","R8i+","RTi+","R1s","R2s","R3s","R4s","R5s","R6s","R7s","R8s","RTs","R1s+","R2s+","R3s+","R4s+","R5s+","R6s+","R7s+","R8s+","RTs+"};

    std::vector<double> w;

    int numAtoms = conf.getNumAtoms();
    const ROMol mol= conf.getOwningMol();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd H = GetHmatrix(Xmean);

    MatrixXd R = GetRmatrix(H, DM, numAtoms);

    w= getGetawayDesc(H,R, ADJ, numAtoms, Heavylist, mol);

    return w;
}


} //end of anonymous namespace


std::vector<double> GETAWAY(const ROMol& mol,int confId){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  std::vector<double> res;

  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

  double Vpoints[3*numAtoms];

  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  int* Heavylist= GetHeavyList(mol); //int nHeavyAt= mol.getNumHeavyAtoms(); // should be the same as the List size upper!

  double *dist3D = MolOps::get3DDistanceMat(mol, confId);

  double *dist = MolOps::getDistanceMat(mol, false); // need to be be set to false to have topological distance not weigthed!

  double *AdjMat = MolOps::getAdjacencyMatrix(mol,false,0,false,0); // false to have only the 1,0 matrix unweighted

  Map<MatrixXd> adj(AdjMat, numAtoms,numAtoms);

  Map<MatrixXd> dm(dist3D, numAtoms,numAtoms);

  Map<MatrixXd> dmtopo(dist, numAtoms,numAtoms);

  res = GetGETAWAY(conf, Vpoints, dm, adj, Heavylist);

  return res;
}


} // end of Descriptors namespace
} // end of RDKit namespace
