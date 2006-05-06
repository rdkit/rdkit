// $Id$
//
// Copyright (C)  2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <math.h>
#include <RDGeneral/Invariant.h>

#include "BFGSOpt.h"

namespace BFGSOpt {




#define CLEANUP() { delete [] grad; delete [] dGrad; delete [] hessDGrad;\
 delete [] newPos; delete [] xi; delete [] invHessian; }

  // ------------------------------------------------------------
  // an adaptation of the function dfpmin (function minimization using the BFGS
  // algorithm) from Chapter 10 of Numerical Recipes in C
  //
  // NOTE: the funcTol argument is *not* used.
  // ------------------------------------------------------------
  int minimize(const int dim,double *pos,
		double gradTol,int &numIters,
		double &funcVal,
		double (*func)(double *),
		void (*gradFunc)(double *,double*),
		double funcTol,int maxIts){
    PRECONDITION(dim>0,"bad dimension");
    PRECONDITION(pos,"bad input array");
    PRECONDITION(gradTol>0,"bad tolerance");
    PRECONDITION(func,"bad function");
    PRECONDITION(gradFunc,"bad function");

    double sum,maxStep,fp;

    double *grad,*dGrad,*hessDGrad;
    double *newPos,*xi;
    double *invHessian;

    grad = new double[dim];
    dGrad = new double[dim];
    hessDGrad = new double[dim];
    newPos = new double[dim];
    xi = new double[dim];
    invHessian = new double[dim*dim];

    // evaluate the function and gradient in our current position:
    fp=func(pos);
    gradFunc(pos,grad);

    sum = 0.0;
    for(int i=0;i<dim;i++){
      int itab=i*dim;
      // initialize the inverse hessian to be identity:
      for(int j=0;j<dim;j++) invHessian[itab+j]=0.0;
      invHessian[itab+i]=1.0;
      // the first line dir is -grad:
      xi[i] = -grad[i];
      sum += pos[i]*pos[i];
    }
    // pick a max step size:
    maxStep = MAXSTEP * maxVal(sqrt(sum),static_cast<double>(dim));


    for(int iter=1;iter<=maxIts;iter++){
      numIters=iter;
      int status;

      // do the line search:
      linearSearch(dim,pos,fp,grad,xi,newPos,funcVal,func,maxStep,status);
      CHECK_INVARIANT(status>=0,"bad direction in linearSearch");

      // save the function value for the next search:
      fp = funcVal;

      // set the direction of this line and save the gradient:
      double test=0.0;
      for(int i=0;i<dim;i++){
	xi[i] = newPos[i]-pos[i];
	pos[i] = newPos[i];
	double temp=fabs(xi[i])/maxVal(fabs(pos[i]),1.0);
	if(temp>test) test=temp;
	dGrad[i] = grad[i];
      }
      if(test<TOLX) {
	CLEANUP();
	return 0;
      }

      // update the gradient:
      gradFunc(pos,grad);
	
      // is the gradient converged?
      test=0.0;
      double term=maxVal(funcVal,1.0);
      for(int i=0;i<dim;i++){
	double temp=fabs(grad[i])*maxVal(fabs(pos[i]),1.0)/term;
	if(temp>test) test=temp;
      }

      //std::cout << "------->>>>>>> Iter: " << iter << " test: " << test << std::endl;
      //if(test<gradTol && deltaFunc<funcTol){
      if(test<gradTol){
	CLEANUP();
	return 0;
      }

      // figure out how much the gradient changed:
      for(int i=0;i<dim;i++){
	dGrad[i] = grad[i]-dGrad[i];
      }
    
      // compute hessian*dGrad:
      double fac=0,fae=0,sumDGrad=0,sumXi=0;
      for(int i=0;i<dim;i++){
	int itab=i*dim;
	hessDGrad[i] = 0.0;
	for(int j=0;j<dim;j++){
	  hessDGrad[i] += invHessian[itab+j]*dGrad[j];
	}

	fac += dGrad[i]*xi[i];
	fae += dGrad[i]*hessDGrad[i];
	sumDGrad += dGrad[i]*dGrad[i];
	sumXi += xi[i]*xi[i];
      }
      if(fac > sqrt(EPS*sumDGrad*sumXi)){
	fac = 1.0/fac;
	double fad = 1.0/fae;
	for(int i=0;i<dim;i++){
	  dGrad[i] = fac*xi[i] - fad*hessDGrad[i];
	}
	for(int i=0;i<dim;i++){
	  int itab=i*dim;
	  for(int j=i;j<dim;j++){
	    invHessian[itab+j] += fac*xi[i]*xi[j] -
	      fad*hessDGrad[i]*hessDGrad[j] +
	      fae*dGrad[i]*dGrad[j];
	    invHessian[j*dim+i] = invHessian[itab+j];
	  }
	}
      }
      // generate the next direction to move:
      for(int i=0;i<dim;i++){
	int itab=i*dim;
	xi[i] = 0.0;
	for(int j=0;j<dim;j++){
	  xi[i] -= invHessian[itab+j]*grad[j];
	}
      }
    }
    CLEANUP();
    return 1;
  }
}
