// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ForceField.h"
#include "Contrib.h"

#include <RDGeneral/Invariant.h>
#include <Numerics/Optimizer/BFGSOpt.h>

namespace ForceFieldsHelper {
  class calcEnergy {
  public:
    calcEnergy(ForceFields::ForceField *ffHolder) :
      mp_ffHolder(ffHolder) {} ;
    double operator() (double *pos) const {
      return mp_ffHolder->calcEnergy(pos);
    }
  private:
    ForceFields::ForceField *mp_ffHolder;    
  };

  class calcGradient {
  public:
    calcGradient(ForceFields::ForceField *ffHolder) :
      mp_ffHolder(ffHolder) {} ;
    double operator() (double *pos,double *grad) const {
      double res = 1.0;
      // the contribs to the gradient function use +=, so we need
      // to zero the grad out before moving on:
      for(unsigned int i=0;i<mp_ffHolder->numPoints()*mp_ffHolder->dimension();i++){
        grad[i] = 0.0;
      }
      mp_ffHolder->calcGrad(pos,grad);

#if 1
      // FIX: this hack reduces the gradients so that the
      // minimizer is more efficient.
      double maxGrad=-1e8;
      double gradScale=0.1;
      for(unsigned int i=0;i<mp_ffHolder->numPoints()*mp_ffHolder->dimension();i++){
        grad[i] *= gradScale;
        if(grad[i]>maxGrad) maxGrad=grad[i];
      }
      // this is a continuation of the same hack to avoid
      // some potential numeric instabilities:
      if(maxGrad>10.0){
        while(maxGrad*gradScale>10.0){
          gradScale*=.5;
        }
        for(unsigned int i=0;i<mp_ffHolder->numPoints()*mp_ffHolder->dimension();i++){
          grad[i] *= gradScale;
        }
      }
      res=gradScale;
#endif
      return res;

    }
  private:
    ForceFields::ForceField *mp_ffHolder;    
  };
}

namespace ForceFields {
  ForceField::~ForceField(){
    d_numPoints=0;
    d_positions.clear();
    d_contribs.clear();
    delete [] dp_distMat;
    dp_distMat=0;
  }

  double ForceField::distance(unsigned int i,unsigned int j,double *pos) {
    PRECONDITION(df_init,"not initialized");
    RANGE_CHECK(0,i,d_numPoints-1);
    RANGE_CHECK(0,j,d_numPoints-1);
    if(j<i){
      int tmp=j;
      j = i;
      i = tmp;
    }
    unsigned int idx=i+j*(j+1)/2;
    CHECK_INVARIANT(idx<d_matSize,"Bad index");
    double &res=dp_distMat[idx];
    if(res<0.0){
      // we need to calculate this distance:
      if(!pos){
        res = 0.0;
        for (unsigned int idx = 0; idx < d_dimension; ++idx) {
          double tmp = (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
          res += tmp*tmp;
        }
      } else {
        res = 0.0;
        for(unsigned int idx=0;idx<d_dimension;idx++){
          double tmp=pos[d_dimension*i+idx]-pos[d_dimension*j+idx];
          res += tmp*tmp;
        }
      }
      res = sqrt(res);
      dp_distMat[idx]=res;
    }
    return res;
  }

  double ForceField::distance(unsigned int i,unsigned int j,double *pos) const {
    PRECONDITION(df_init,"not initialized");
    RANGE_CHECK(0,i,d_numPoints-1);
    RANGE_CHECK(0,j,d_numPoints-1);
    if(j<i){
      int tmp=j;
      j = i;
      i = tmp;
    }
    double res;
    if(!pos){
      res = 0.0;
      for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        double tmp = (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        res += tmp*tmp;
      }
    } else {
      res = 0.0;
      for(unsigned int idx=0;idx<d_dimension;idx++){
        double tmp=pos[d_dimension*i+idx]-pos[d_dimension*j+idx];
        res += tmp*tmp;
      }
    }
    res = sqrt(res);
    return res;
  }

  void ForceField::initialize(){
    // clean up if we have used this already:
    df_init=false;
    delete [] dp_distMat;
    dp_distMat=0;
    
    d_numPoints = d_positions.size();
    d_matSize=d_numPoints*(d_numPoints+1)/2;
    dp_distMat = new double[d_matSize];
    this->initDistanceMatrix();
    df_init=true;
  }

  int ForceField::minimize(unsigned int maxIts,double forceTol,double energyTol){
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(static_cast<unsigned int>(d_numPoints)==d_positions.size(),"size mismatch");
    if(d_contribs.empty()) return 0;

    unsigned int numIters=0;
    unsigned int dim=this->d_numPoints*d_dimension;
    double finalForce;
    double *points=new double[dim];

    this->scatter(points);
    ForceFieldsHelper::calcEnergy eCalc(this);
    ForceFieldsHelper::calcGradient gCalc(this);
    
    int res = BFGSOpt::minimize(dim,points,forceTol,numIters,finalForce,
                                eCalc,gCalc,
                                energyTol,maxIts);
    this->gather(points);

    delete [] points;
    return res;
  }

  double ForceField::calcEnergy() const{
    PRECONDITION(df_init,"not initialized");
    double res = 0.0;
    if(d_contribs.empty()) return res;

    unsigned int N = d_positions.size();
    double *pos = new double[d_dimension*N];
    this->scatter(pos);
    // now loop over the contribs
    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
        contrib != d_contribs.end();contrib++){
      res += (*contrib)->getEnergy(pos);
    }
    delete [] pos;
    return res;    
  }

  double ForceField::calcEnergy(double *pos){
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");
    double res = 0.0;

    this->initDistanceMatrix();
    if(d_contribs.empty()) return res;

    // now loop over the contribs
    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
        contrib != d_contribs.end();contrib++){
      double E=(*contrib)->getEnergy(pos);
      res += E;
    }
    return res;
  }

  void ForceField::calcGrad(double *grad) const {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(grad,"bad gradient vector");
    if(d_contribs.empty()) return;

    unsigned int N = d_positions.size();
    double *pos = new double[d_dimension*N];
    this->scatter(pos);
    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
        contrib != d_contribs.end();contrib++){
      (*contrib)->getGrad(pos,grad);
    }
    // zero out gradient values for any fixed points:
    for(INT_VECT::const_iterator it=d_fixedPoints.begin();
        it!=d_fixedPoints.end();it++){
      CHECK_INVARIANT(static_cast<unsigned int>(*it)<d_numPoints,"bad fixed point index");
      unsigned int idx=d_dimension*(*it);
      for (unsigned int di = 0; di < this->dimension(); ++di) {
        grad[idx+di] = 0.0;
      }
    }
    delete [] pos;
  }
  void ForceField::calcGrad(double *pos,double *grad) {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");
    PRECONDITION(grad,"bad gradient vector");
    if(d_contribs.empty()) return;

    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
        contrib != d_contribs.end();contrib++){
      (*contrib)->getGrad(pos,grad);
    }

    for(INT_VECT::const_iterator it=d_fixedPoints.begin();
        it!=d_fixedPoints.end();it++){
      CHECK_INVARIANT(static_cast<unsigned int>(*it)<d_numPoints,"bad fixed point index");
      unsigned int idx=d_dimension*(*it);
      for (unsigned int di = 0; di < this->dimension(); ++di) {
        grad[idx+di] = 0.0;
      }
    }
  }

  void ForceField::scatter(double *pos) const {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");

    unsigned int tab=0;
    for(unsigned int i=0;i<d_positions.size();i++){
      for (unsigned int di=0; di < this->dimension(); ++di){
        pos[tab+di] = (*d_positions[i])[di]; //->x;
      }
      tab+=this->dimension();
    }
    POSTCONDITION(tab==this->dimension()*d_positions.size(),"bad index");
  }

  void ForceField::gather(double *pos) {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");
    
    unsigned int tab=0;
    for(unsigned int i=0;i<d_positions.size();i++){
      for (unsigned int di=0; di < this->dimension(); ++di){
        (*d_positions[i])[di] = pos[tab+di];
      }
      tab+=this->dimension();
    }
  }


  void ForceField::initDistanceMatrix(){
    PRECONDITION(d_numPoints,"no points");
    PRECONDITION(dp_distMat,"no distance matrix");
    PRECONDITION(static_cast<unsigned int>(d_numPoints*(d_numPoints+1)/2)<=d_matSize,"matrix size mismatch");
    for(unsigned int i=0;i<d_numPoints*(d_numPoints+1)/2;i++){
      dp_distMat[i]=-1.0;
    }
  }
  

}
