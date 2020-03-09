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
#include <GraphMol/new_canon.h>

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

std::vector<double> CopyVect(VectorXd v1){
  std::vector<double> v2;
  v2.resize(v1.size());
  VectorXd::Map(&v2[0], v1.size()) = v1;
  return v2;
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




// this method use the new Canonical Ranking Atoms to gene we don't need it there but can be useful in another package 
std::vector<unsigned int> CanonicalRankAtoms(const ROMol &mol,
                                             bool breakTies = true,
                                             bool includeChirality = true,
                                             bool includeIsotopes = true) {

  
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  Canon::rankMolAtoms(mol, ranks, breakTies, includeChirality, includeIsotopes);
  return ranks;
}


// center Score to the closest Atom per axis
std::vector<double> centeringVector(std::vector<double> V){
double D;
double AD=1000000.0;
for (int i=0;i<V.size();i++){
    if (std::abs(V[i])<AD) {
      D=V[i];
      AD=abs(V[i]);
    }
}

for (int i=0;i<V.size();i++){
      V[i]=V[i]-D;
    }


return V;

}


double* getWhimDesc(JacobiSVD<MatrixXd> svd, MatrixXd Xmean, int numAtoms, double th, bool printscore) {


  // takes only L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u

    double *w = new double[17];
    // prepare data for Whim parameter computation
    
    double *SingVal = retreiveVect(svd.singularValues());
    MatrixXd Scores= Xmean*svd.matrixV(); //  V is similar

    // Part I : directionnal descriptors
    // Singular values
    // directional WHIM size descriptors (or d-WSIZ indices)
    w[0]=SingVal[0]; // lambda
    w[1]=SingVal[1];
    w[2]=SingVal[2];

    // directional WHIM shape descriptors (or d-WSHA indices)
    double totalSingVal = SingVal[0]+SingVal[1]+SingVal[2];
    w[3]=SingVal[0]/totalSingVal; // P1 (nu)
    w[4]=SingVal[1]/totalSingVal; // p2
    
    // Start Gamma
    // directional WHIM symmetry descriptors (or d-WSYM indices)

    double gamma[3]; // Gamma values
    double nAT = (double) numAtoms;
    // check if two atoms are symetric versus the new axis ie newx,newy,newz a
    for (int i=0;i<3; i++) {
      for (int j =0;j<numAtoms;j++) {
        Scores(j,i)=roundn(Scores(j,i),3); // round the matrix! same as eigen tolerance !
      }
    }


    // Copy vectors to std:vector center & reset Score vector
    std::vector<double> VS;

    for (int i=0;i<3; i++) {
        VS = centeringVector(CopyVect(Scores.col(i)));
        for (int j=0;j<numAtoms;j++){
          Scores(j,i)=VS[j];
        }
    }


    // Get symetry based on centered Scores
    for (int i=0;i<3; i++) {
      double ns=0.0;
      for (int j=0;j< numAtoms;j++) {
        for (int k=0;k< numAtoms;k++){

          if (j==k) continue;

          if (std::abs(Scores(j,i) + Scores(k,i))< 2*th and (std::abs(Scores(j,i))>th or std::abs(Scores(k,i))>th )) {  // those that are close opposite & not close to the axis!
            ns++; // check only once the symetric none null we need to add +2!
            break;
          }
        }
      }

      // for all atoms close to the axis we need to add +1!
      for (int j=0;j< numAtoms;j++) {
        if (std::abs(Scores(j,i))<=th) {
        ns++; // atom close to the the axis are symetric! 
        } 
      }

      double na=0.0;
      na = nAT-ns;
      gamma[i] =0.0;
      if (ns>0) {
           gamma[i] = (ns / nAT) * log(ns / nAT) / log(2.0) + (na / nAT) * log(1.0 / nAT) / log(2.0);
           gamma[i] = 1.0 / (1.0 - gamma[i]);
      }
    }


    // check if the molecule is fully symmetrical "like CH4" using Canonical Rank Index and/or Sphericity !
    // case of complete symetry of two Components than there are set to 1!
    if (SingVal[0]==SingVal[1]) {
      gamma[0]=1;
      gamma[1]=1;
    }

    if (SingVal[1]==SingVal[2]) {
      gamma[1]=1;
      gamma[2]=1;
    }

    w[5]=gamma[0]; // G1 (gamma)
    w[6]=gamma[1]; // G2
    w[7]=gamma[2]; // G3
    // end gamma
    
    // directional WHIM density descriptors (or d-WDEN indices, or WHIM emptiness)
    // the inverse of the kurtosis  Eta
    MatrixXd Scores2= Xmean*svd.matrixV(); //  V is similar

    VectorXd v1= Scores2.col(0);
    VectorXd v2= Scores2.col(1);
    VectorXd v3= Scores2.col(2);
    w[8] =  numAtoms*pow(w[0],2)/ v1.array().pow(4).sum(); // E1 (eta)
    w[9] =  numAtoms*pow(w[1],2)/ v2.array().pow(4).sum(); // E2
    w[10] = numAtoms*pow(w[2],2)/ v3.array().pow(4).sum(); // E3

   
    // Part II : global WHIM T,A,G,K,D,V
    // linear contributions to the total molecular size
    w[11]=SingVal[0]+SingVal[1]+SingVal[2];  // T
    // quadratic contributions to the total molecular size
    w[12]=SingVal[0]*SingVal[1]+SingVal[0]*SingVal[2]+SingVal[1]*SingVal[2]; // A
    // WHIM symmetry
    w[13]=pow(gamma[0]*gamma[1]*gamma[2],1.0/3.0); // G
    
    // WHIM shape (or WSHA index)
    double res=0.0;
    for (int i=0;i<3;i++) 
    {
      res+=std::abs(w[i]/w[11]-1.0/3.0);
    }
    w[14]=res*3.0/4.0; // K
    // WHIM density or  WSYM index
    w[15] = (w[8]+w[9]+w[10])/3; // mean total density of the atoms called D of the E1, E2,E3 in dragon it's /3!
    // V is the complete expansion
    w[16]=w[11]+w[12]+SingVal[0]*SingVal[1]*SingVal[2]; // V

    return w;

}


double* GetWHIMU(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){

    double *w = new double[17];

    int numAtoms = conf.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);


    //std::cout << "X:\n";

    //std::cout << matorigin <<"\n";

    MatrixXd MatOrigin=matorigin.transpose();

    MatrixXd Weigth;

    Weigth.setIdentity(numAtoms,numAtoms);

    double weigth=Weigth.diagonal().sum();

    MatrixXd Xmean = GetCenterMatrix(MatOrigin);

    MatrixXd covmat = GetCovMatrix(Xmean,Weigth,weigth);

    JacobiSVD<MatrixXd> svd = getSVD(covmat);


    //std::cout << "covmat\n";  
    //std::cout << covmat <<"\n";


    //std::cout << "Xmean\n";  
    //std::cout << Xmean <<"\n";

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



double*  GetWHIMMass(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){
    double *w = new double[17];

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

double*  GetWHIMvdw(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){
    double *w = new double[17];

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


double*  GetWHIMneg(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){
    double *w = new double[17];

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




double*  GetWHIMpol(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){
    double *w = new double[17];

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



double*  GetWHIMIonPol(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){
      double *w = new double[17];

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


double*  GetWHIMIState(const Conformer &conf, double Vpoints[], double th, std::vector<unsigned int> CRA){
    double *w = new double[17];

    ROMol& mol = conf.getOwningMol();

    MolOps::removeHs(mol, false, false);

    int numAtoms = mol.getNumAtoms();

    Map<MatrixXd> matorigin(Vpoints, 3,numAtoms);

    MatrixXd MatOrigin=matorigin.transpose();

    std::vector<double> weigthvector = moldata3D.GetEState2(mol); // caution not only neighours hum not sure on this based on the paper! Also this should be only on Heavy atoms

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




std::vector<double> WHIM(const ROMol& mol, int confId, double th){
  PRECONDITION(mol.getNumConformers()>=1,"molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  const Conformer &conf = mol.getConformer(confId);

//std::cout <<"CanonicalRankAtoms:\n";

std::vector<unsigned int> CRA = CanonicalRankAtoms(mol,false,true, true); // BreakTies to false
//for (int i=0;i<numAtoms;i++){
//  std::cout << CRA[i] << ",";

//}
//std::cout <<"\n";


// uv not yet created
// std::map<unsigned int, int> freq_map;
// for (auto const & x : CRA)
//     ++freq_map[x];
// std::vector<unsigned int> uv;
// std::vector<int> freq_uv;
// for (auto const & p : freq_map)
// {
//     uv.push_back(p.first);
//     freq_uv.push_back(p.second);
//     std::cout << p.first << "," << p.second << "|";
// }
// std::cout <<"\n";






  double Vpoints[3*numAtoms];

  for(int i=0; i<numAtoms; ++i){
     Vpoints[3*i]   =conf.getAtomPos(i).x;
     Vpoints[3*i+1] =conf.getAtomPos(i).y;
     Vpoints[3*i+2] =conf.getAtomPos(i).z;
  }

  std::vector<double> res;

  double* wu= GetWHIMU(conf, Vpoints, th, CRA);
  double* wm= GetWHIMMass(conf, Vpoints, th, CRA);
  double* wv= GetWHIMvdw(conf, Vpoints, th, CRA);
  double* we= GetWHIMneg(conf, Vpoints, th, CRA);
  double* wp= GetWHIMpol(conf, Vpoints, th, CRA);
  double* wi= GetWHIMIonPol(conf, Vpoints, th, CRA);
  double* ws= GetWHIMIState(conf, Vpoints, th, CRA);


  // takes only L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u

  // std::vector<std::string> descnames={"L1u","L2u","L3u","Tu","Au","Vu","P1u","P2u","P3u","Ku","E1u","E2u","E3u","Du","G1u","G2u","G3u","Gu"};
    for (int i=0;i<11;i++) {
        res.push_back(roundn(wu[i],3));

     }

   // descnames={"L1m","L2m","L3m","Tm","Am","Vm","P1m","P2m","P3m","Km","E1m","E2m","E3m","Dm","G1m","G2m","G3m","Gm"};

    for (int i=0;i<11;i++) {
        res.push_back(roundn(wm[i],3));

     }

  //   descnames={"L1v","L2v","L3v","Tv","Av","Vv","P1v","P2v","P3v","Kv","E1v","E2v","E3v","Dv","G1v","G2v","G3v","Gv"};

    for (int i=0;i<11;i++) {
        res.push_back(roundn(wv[i],3));
     }

   // descnames={"L1e","L2e","L3e","Te","Ae","Ve","P1e","P2e","P3e","Ke","E1e","E2e","E3e","De","G1e","G2e","G3e","Ge"};

    for (int i=0;i<11;i++) {
        res.push_back( roundn(we[i],3));
     }


   // descnames={"L1p","L2p","L3p","Tp","Ap","Vp","P1p","P2p","P3p","Kp","E1p","E2p","E3p","Dp","G1p","G2p","G3p","Gp"};

    for (int i=0;i<11;i++) {
        res.push_back( roundn(wp[i],3));

     }
    
   // descnames={"L1i","L2i","L3i","Ti","Ai","Vi","P1i","P2i","P3i","Ki","E1i","E2i","E3i","Di","G1i","G2i","G3i","Gi"};

    for (int i=0;i<11;i++) {
       res.push_back( roundn(wi[i],3));

     }

   //descnames={"L1s","L2s","L3s","Ts","As","Vs","P1s","P2s","P3s","Ks","E1s","E2s","E3s","Ds","G1s","G2s","G3s","Gs"};

    for (int i=0;i<11;i++) {
        res.push_back(roundn(ws[i],3));

     }

       //78  79  80  81  82  83  84
       //Tu  Tm  Tv  Te  Tp  Ti  Ts
      res.push_back(roundn(wu[11],3));
      res.push_back(roundn(wm[11],3));
      res.push_back(roundn(wv[11],3));
      res.push_back(roundn(we[11],3));
      res.push_back(roundn(wp[11],3));
      res.push_back(roundn(wi[11],3));
      res.push_back(roundn(ws[11],3));

      //85  86  87  88  89  90  91
      //Au  Am  Av  Ae  Ap  Ai  As
      res.push_back(roundn(wu[12],3));
      res.push_back(roundn(wm[12],3));
      res.push_back(roundn(wv[12],3));
      res.push_back(roundn(we[12],3));
      res.push_back(roundn(wp[12],3));
      res.push_back(roundn(wi[12],3));
      res.push_back(roundn(ws[12],3));

      //92  93
      //Gu  Gm
      res.push_back(roundn(wu[13],3));
      res.push_back(roundn(wm[13],3));

      //94  95  96  97  98  99  100
      //Ku  Km  Kv  Ke  Kp  Ki  Ks
      res.push_back(roundn(wu[14],3));
      res.push_back(roundn(wm[14],3));
      res.push_back(roundn(wv[14],3));
      res.push_back(roundn(we[14],3));
      res.push_back(roundn(wp[14],3));
      res.push_back(roundn(wi[14],3));
      res.push_back(roundn(ws[14],3));

      //101 102 103 104 105 106 107
       //Du  Dm  Dv  De  Dp  Di  Ds
      res.push_back(roundn(wu[15],3));
      res.push_back(roundn(wm[15],3));
      res.push_back(roundn(wv[15],3));
      res.push_back(roundn(we[15],3));
      res.push_back(roundn(wp[15],3));
      res.push_back(roundn(wi[15],3));
      res.push_back(roundn(ws[15],3));


      //108 109 110 111 112 113 114
       //Vu  Vm  Vv  Ve  Vp  Vi  Vs
      res.push_back(roundn(wu[16],3));
      res.push_back(roundn(wm[16],3));
      res.push_back(roundn(wv[16],3));
      res.push_back(roundn(we[16],3));
      res.push_back(roundn(wp[16],3));
      res.push_back(roundn(wi[16],3));
      res.push_back(roundn(ws[16],3));

//L1u L2u L3u P1u P2u G1u G2u G3u E1u E2u E3u L1m L2m L3m P1m P2m G1m G2m G3m E1m E2m E3m L1v L2v L3v P1v P2v G1v G2v G3v E1v E2v E3v L1e L2e L3e P1e P2e G1e G2e G3e E1e E2e E3e L1p L2p L3p P1p P2p G1p G2p G3p E1p E2p E3p L1i L2i L3i P1i P2i G1i G2i G3i E1i E2i E3i L1s L2s L3s P1s P2s G1s G2s G3s E1s E2s E3s Tu  Tm  Tv  Te  Tp  Ti  Ts  Au  Am  Av  Ae  Ap  Ai  As  Gu  Gm  Ku  Km  Kv  Ke  Kp  Ki  Ks  Du  Dm  Dv  De  Dp  Di  Ds  Vu  Vm  Vv  Ve  Vp  Vi  Vs

  return res;
}

} // end of Descriptors namespace
} // end of RDKit namespace
