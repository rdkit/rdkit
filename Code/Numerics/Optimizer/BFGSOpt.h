//
// Copyright (C)  2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

namespace BFGSOpt {
  const double FUNCTOL=1e-4;  //!< Default tolerance for function convergence in the minimizer
  const double MOVETOL=1e-7;  //!< Default tolerance for x changes in the minimizer
  const int    MAXITS=200;    //!< Default maximum number of iterations
  const double EPS=3e-8;      //!< Default gradient tolerance in the minimizer
  const double TOLX=4.*EPS;   //!< Default direction vector tolerance in the minimizer
  const double MAXSTEP=100.0; //!< Default maximim step size in the minimizer

  //! Do a Quasi-Newton minimization along a line.  
  /*!
    See Numerical Recipes in C, Section 9.7 for a description of the algorithm.
  
     \param dim     the dimensionality of the space.
     \param oldPt   the current position, as an array.
     \param oldVal  the current function value.
     \param grad    the value of the function gradient at oldPt
     \param dir     the minimization direction
     \param newPt   used to return the final position
     \param newVal  used to return the final function value
     \param func    the function to minimize
     \param maxStep the maximum allowable step size
     \param resCode used to return the results of the search.
    
     Possible values for resCode are on return are:
      -  0: success
      -  1: the stepsize got too small.  This probably indicates success.
      - -1: the direction is bad (orthogonal to the gradient)
  */
  void linearSearch(unsigned int dim,double *oldPt,double oldVal,
                    double *grad,double *dir,double *newPt,
                    double &newVal,
                    double (*func)(double *),
                    double maxStep,int &resCode);

  //! Do a BFGS minimization of a function.
  /*!
     See Numerical Recipes in C, Section 10.7 for a description of the algorithm.
    
     \param dim     the dimensionality of the space.
     \param pos   the starting position, as an array.
     \param gradTol tolerance for gradient convergence
     \param numIters used to return the number of iterations required
     \param funcVal  used to return the final function value
     \param func    the function to minimize
     \param gradFunc  calculates the gradient of func
     \param funcTol tolerance for changes in the function value for convergence.
     \param maxIts   maximum number of iterations allowed
    
     \return a flag indicating success (or type of failure). Possible values are:
      -  0: success
      -  1: too many iterations were required
  */
  int minimize(unsigned int dim,double *pos,
               double gradTol,
               unsigned int &numIters,
               double &funcVal,
               double (*func)(double *),
               void (*gradFunc)(double *,double*),
               double funcTol=TOLX,
               unsigned int maxIts=MAXITS);
}
