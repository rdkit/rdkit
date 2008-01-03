// $Id$
//
// Copyright (C)  2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <math.h>
#include <RDGeneral/Invariant.h>

#include "BFGSOpt.h"

namespace BFGSOpt {

  // ------------------------------------------------------------
  // an adaptation of the function lnsearch from Chapter 9 of
  // Numerical Recipes in C
  //   possible return values:
  //    -1: bad direction
  //     0: converged
  //     1: step got too small, probably converged
  // ------------------------------------------------------------
  void linearSearch(unsigned int dim,double *oldPt,double oldVal,
		    double *grad,double *dir,double *newPt,
		    double &newVal,
		    double (*func)(double *),
		    double maxStep,int &resCode){
    PRECONDITION(oldPt,"bad input array");
    PRECONDITION(grad,"bad input array");
    PRECONDITION(dir,"bad input array");
    PRECONDITION(newPt,"bad input array");
    PRECONDITION(func,"bad function");

    double sum=0.0,slope=0.0,test=0.0,lambda=0.0;
    double lambda2=0.0,lambdaMin=0.0,tmpLambda=0.0,val2=0.0;

    resCode=0;

    // get the length of the direction vector:
    sum=0.0;
    for(unsigned int i=0;i<dim;i++)
      sum +=dir[i]*dir[i];
    sum=sqrt(sum);

    // rescale if we're trying to move too far:
    if(sum>maxStep){
      for(unsigned int i=0;i<dim;i++)
	dir[i] *= maxStep/sum;
    }
      
    // make sure our direction has at least some component along
    // -grad
    slope=0.0;
    for(unsigned int i=0;i<dim;i++){
      slope += dir[i]*grad[i];
    }
    if(slope>=0.0){
      resCode=-1;
      return;
    }

    test=0.0;
    for(unsigned int i=0;i<dim;i++){
      double temp=fabs(dir[i])/maxVal(fabs(oldPt[i]),1.0);
      if(temp>test) test=temp;
    }

    lambdaMin = MOVETOL/test;
    lambda = 1.0;
    while(1){
      for(unsigned int i=0;i<dim;i++){
	newPt[i]=oldPt[i]+lambda*dir[i];
      }
      newVal = func(newPt);
      //std::cout << "\t" << lambda << " " << lambdaMin << std::endl;
      if(lambda<lambdaMin){
	// the position change is too small... set the resCode and return
	for(unsigned int i=0;i<dim;i++){
	  newPt[i]=oldPt[i];
        }
	resCode=1;
	return;
      } else if( newVal-oldVal <= FUNCTOL*lambda*slope ){
	// we're converged on the function:
	return;
      } else {
	// if we made it this far, we need to backtrack:
	if(lambda==1.0){
	  // it's the first step:
	  tmpLambda = -slope / (2.0*(newVal-oldVal-slope));
	} else {
	  double rhs1 = newVal-oldVal-lambda*slope;
	  double rhs2 = val2-oldVal-lambda2*slope;
	  double a = (rhs1/(lambda*lambda) - rhs2/(lambda2*lambda2))/
	    (lambda-lambda2);
	  double b = (-lambda2*rhs1/(lambda*lambda)+lambda*rhs2/(lambda2*lambda2))/
	    (lambda-lambda2);
	  if( a==0.0 ){
	    tmpLambda = -slope/(2.0*b);
	  } else {
	    double disc=b*b-3*a*slope;
	    if(disc<0.0){
	      tmpLambda = 0.5*lambda;
	    } else if(b<=0.0) {
	      tmpLambda = (-b+sqrt(disc))/(3.0*a);
	    } else {
	      tmpLambda = -slope/(b+sqrt(disc));
	    }
	  }
	  if( tmpLambda > 0.5*lambda ){
	    tmpLambda = 0.5*lambda;
	  }
	}
      }
      lambda2 = lambda;
      val2 = newVal;
      lambda = maxVal(tmpLambda,0.1*lambda);
    }
  }

}
