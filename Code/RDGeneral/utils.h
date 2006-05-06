//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
//
#ifndef __RD_UTILS_H__
#define __RD_UTILS_H__

#include "types.h"

#include <boost/random.hpp>

namespace RDKit{
  const int NUM_PRIMES_AVAIL = 1000; //!< the number of primes available and stored
  extern int firstThousandPrimes[NUM_PRIMES_AVAIL];
  
  const int FILE_MAXLINE=256;  //!< an assumed maximum length for lines read from files
  
  //! \brief compute the product of the set of primes corresponding to the
  //!        values in an INT_VECT 
  double computeIntVectPrimesProduct(const INT_VECT &ring);

  //! floating point comparison with a tolerance
  bool feq(double v1,double v2,double tol=1e-4);
  
  typedef boost::minstd_rand rng_type;
  typedef boost::uniform_int<> uniform_int;
  typedef boost::uniform_real<> uniform_double;
  typedef boost::variate_generator<rng_type &,uniform_int> int_source_type;
  typedef boost::variate_generator<rng_type &,uniform_double> double_source_type;

  //! Optionally seed and return a reference to the global (Boost) random generator
  rng_type &getRandomGenerator(int seed=-1);

  //! Return a random double value between 0.0 and 1.0
  //! Optionally seed the random number generator 
  double getRandomVal(int seed = -1);
  
}



#endif
