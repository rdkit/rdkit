// $Id: ForceField.cpp 4954 2006-02-17 23:33:34Z glandrum $
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "ForceField.h"
#include "Contrib.h"

#include <RDGeneral/Invariant.h>
#include <Numerics/Optimizer/BFGSOpt.h>

namespace ForceFieldsHelper {
  ForceFields::ForceField *_ffHolder;
  double calcEnergy(double *pos){
    return _ffHolder->calcEnergy(pos);
  };
  void calcGrad(double *pos,double *grad){
    // the contribs to the gradient function use +=, so we need
    // to zero the grad out before moving on:
    for(int i=0;i<_ffHolder->numPoints()*3;i++){
      grad[i] = 0.0;
    }
    _ffHolder->calcGrad(pos,grad);

    double maxGrad=-1e8;
    // FIX: this hack reduces the gradients so that the
    // minimizer is more efficient.
    double gradScale=0.1;
    for(int i=0;i<_ffHolder->numPoints()*3;i++){
      grad[i] *= gradScale;
      if(grad[i]>maxGrad) maxGrad=grad[i];
    }
    // this is a continuation of the same hack to avoid
    // some potential numeric instabilities:
    if(maxGrad>10.0){
      while(maxGrad>10.0){
	maxGrad*=gradScale;
      }
      for(int i=0;i<_ffHolder->numPoints()*3;i++){
	grad[i] *= gradScale;
      }
    }
  }
}

namespace ForceFields {
  ForceField::~ForceField(){
    //std::cerr << " *** destroy force field " << std::endl;
    d_numPoints=0;
    d_positions.clear();
    d_forces.clear();
    d_contribs.clear();
    if(dp_distMat) delete [] dp_distMat;
    dp_distMat=0;
    //std::cerr << " *** destroy force field DONE" << std::endl;
  }

  double ForceField::distance(int i,int j,double *pos) {
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
    //std::cerr << "        d: " << i << "," << j <<"=" << res << std::endl;
    if(res<0.0){
      // we need to calculate this distance:
      if(!pos){
	RDGeom::Point3D v = *(this->positions()[i])-*(this->positions()[j]);
	res = v.length();
      } else {
	res = 0.0;
	for(int idx=0;idx<3;idx++){
	  double tmp=pos[3*i+idx]-pos[3*j+idx];
	  res += tmp*tmp;
	}
	res = sqrt(res);
      }
    }
    //std::cerr << "        d: " << i << "," << j <<"=" << res << std::endl;
    return res;
  }

  double ForceField::distance(int i,int j,double *pos) const {
    PRECONDITION(df_init,"not initialized");
    RANGE_CHECK(0,i,d_numPoints-1);
    RANGE_CHECK(0,j,d_numPoints-1);
    if(j<i){
      int tmp=j;
      j = i;
      i = tmp;
    }
    double res; //=dp_distMat[i+j*(j+1)/2];
    //if(res<0.0){
    if(!pos){
      RDGeom::Point3D v = *(this->positions()[i])-*(this->positions()[j]);
      res = v.length();
    } else {
      res = 0.0;
      for(int idx=0;idx<3;idx++){
        double tmp=pos[3*i+idx]-pos[3*j+idx];
        res += tmp*tmp;
      }
      res = sqrt(res);
    }
    //}
    return res;
  }

  void ForceField::initialize(){
    //PRECONDITION(d_positions.size()==d_forces.size(),"bad position or force vector");

    // clean up if we have used this already:
    df_init=false;
    if(dp_distMat) delete [] dp_distMat;
    dp_distMat=0;
    
    d_numPoints = d_positions.size();
    d_matSize=d_numPoints*(d_numPoints+1)/2;
    dp_distMat = new double[d_matSize];
    this->initDistanceMatrix();
    df_init=true;
  }

  int ForceField::minimize(int maxIts,double forceTol,double energyTol){
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(static_cast<unsigned int>(d_numPoints)==d_positions.size(),"size mismatch");
    int numIters=0;
    int dim=this->d_numPoints*3;
    double finalForce;
    double *points=new double[dim];
    // FIX: nothing is ever being done with these forces, and they cannot
    // currently be updated.
    double *forces=0;
    //forces = new double[dim];

    this->scatter(points,forces);
    ForceFieldsHelper::_ffHolder=this;
    int res = BFGSOpt::minimize(dim,points,forceTol,numIters,finalForce,
			    ForceFieldsHelper::calcEnergy,
			    ForceFieldsHelper::calcGrad,energyTol,maxIts);
    this->gather(points,forces);

    delete [] points;
    //delete [] forces;
    return res;
  }

  double ForceField::calcEnergy() const{
    PRECONDITION(df_init,"not initialized");
    double res = 0.0;

    //this->initDistanceMatrix();
    int N = d_positions.size();
    double *pos = new double[3*N];
    this->scatter(pos,0);
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

#if 0
    std::cout << "************ Es" << std::endl;
    for(unsigned int i=0;i<d_positions.size();++i){
      std::cout << "\t" << i << " " << pos[3*i] << " " << pos[3*i+1] << " " << pos[3*i+2] << std::endl;
    }
#endif
    // now loop over the contribs
    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
	contrib != d_contribs.end();contrib++){
      double E=(*contrib)->getEnergy(pos);
      //std::cout << "\t" << count++ << " " << E << std::endl;
      res += E;
    }
    //std::cout << "  E: " << res << std::endl;
    return res;
  }

  void ForceField::calcGrad(double *grad) const {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(grad,"bad gradient vector");
    //this->initDistanceMatrix();
    int N = d_positions.size();
    double *pos = new double[3*N];
    this->scatter(pos,0);
    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
	contrib != d_contribs.end();contrib++){
      (*contrib)->getGrad(pos,grad);
    }
    // zero out gradient values for any fixed points:
    for(INT_VECT::const_iterator it=d_fixedPoints.begin();
	it!=d_fixedPoints.end();it++){
      CHECK_INVARIANT(*it<d_numPoints,"bad fixed point index");
      unsigned int idx=3*(*it);
      //std::cerr << "&&&& set: " << *it << " -> " << idx << std::endl;
      grad[idx] = 0.0;
      grad[idx+1] = 0.0;
      grad[idx+2] = 0.0;
    }
    delete [] pos;
  }
  void ForceField::calcGrad(double *pos,double *grad) {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");
    PRECONDITION(grad,"bad gradient vector");

    for(ContribPtrVect::const_iterator contrib=d_contribs.begin();
	contrib != d_contribs.end();contrib++){
      (*contrib)->getGrad(pos,grad);
#if 0
      for(unsigned int i=0;i<d_positions.size();++i){
	if(grad[3*i]!=grad[3*i]){
	  std::cout << "NaN found" << std::endl;
	  ContribPtr contribP=*contrib;
	  std::cout << "\t\t" << typeid(*(contribP.get())).name() << std::endl;
	}
      }
#endif

    }


#if 0
    std::cout << "GRAD:" << std::endl;
    for(unsigned int i=0;i<d_positions.size();++i){
      std::cout << "\t" << i << " " << grad[3*i] << " " << grad[3*i+1] << " " << grad[3*i+2] << " NAN? " << (grad[3*i]!=grad[3*i]) << std::endl;
    }
#endif
    for(INT_VECT::const_iterator it=d_fixedPoints.begin();
	it!=d_fixedPoints.end();it++){
      CHECK_INVARIANT(*it<d_numPoints,"bad fixed point index");
      unsigned int idx=3*(*it);
      //std::cerr << "&&&& set: " << *it << " -> " << idx << std::endl;
      grad[idx] = 0.0;
      grad[idx+1] = 0.0;
      grad[idx+2] = 0.0;
    }

  }

  void ForceField::scatter(double *pos,double *grad) const {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");
    //PRECONDITION(grad,"bad gradient vector");

    unsigned int tab=0;
    for(unsigned int i=0;i<d_positions.size();i++){
      pos[tab] = d_positions[i]->x;
      pos[tab+1] = d_positions[i]->y;
      pos[tab+2] = d_positions[i]->z;
#if 0
      grad[tab] = d_forces[i]->x;
      grad[tab+1] = d_forces[i]->y;
      grad[tab+2] = d_forces[i]->z;
#endif
      tab+=3;
    }
    POSTCONDITION(tab==3*d_positions.size(),"bad index");
  }

  void ForceField::gather(double *pos,double *grad) {
    PRECONDITION(df_init,"not initialized");
    PRECONDITION(pos,"bad position vector");
    //PRECONDITION(grad,"bad gradient vector");

    unsigned int tab=0;
    for(unsigned int i=0;i<d_positions.size();i++){
      d_positions[i]->x = pos[tab];
      d_positions[i]->y = pos[tab+1];
      d_positions[i]->z = pos[tab+2];
#if 0
      d_forces[i]->x = grad[tab];
      d_forces[i]->y = grad[tab+1];
      d_forces[i]->z = grad[tab+2];
#endif
      tab+=3;
    }
  }


  void ForceField::initDistanceMatrix(){
    PRECONDITION(d_numPoints,"no points");
    PRECONDITION(dp_distMat,"no distance matrix");
    PRECONDITION(static_cast<unsigned int>(d_numPoints*(d_numPoints+1)/2)<=d_matSize,"matrix size mismatch");
    for(int i=0;i<d_numPoints*(d_numPoints+1)/2;i++){
      dp_distMat[i]=-1.0;
    }
  }
  

}
