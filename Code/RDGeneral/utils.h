//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#ifndef __RD_UTILS_H__
#define __RD_UTILS_H__

#include "types.h"
#include <RDGeneral/Invariant.h>

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

  //! return a reference to the global (Boost) random source
  double_source_type &getDoubleRandomSource();
  

  template <class T>
  unsigned int countSwapsToInterconvert(const T &ref,T probe){
    PRECONDITION(ref.size()==probe.size(),"size mismatch");
    typename T::const_iterator refIt=ref.begin();
    typename T::iterator probeIt=probe.begin();
    typename T::iterator probeIt2;

    unsigned int nSwaps = 0;
    while(refIt!=ref.end()){
      if((*probeIt)!=(*refIt)){
        bool foundIt=false;
        probeIt2=probeIt;
        while((*probeIt2)!=(*refIt) && probeIt2!=probe.end()){
          ++probeIt2;
        }
        if(probeIt2 != probe.end()){
          foundIt=true;
        }
        CHECK_INVARIANT(foundIt,"could not find probe element");

        std::swap(*probeIt,*probeIt2);
        nSwaps++;
      }
      ++probeIt;
      ++refIt;
    }
    return nSwaps;
  }
}

// contribution from dkoes
template<unsigned n>
inline double int_pow(double x) {
  double half =  int_pow<n/2>(x);
  if(n%2 == 0) //even
    return half*half;
  else
    return half*half*x;
}

template<>
inline double int_pow<0>(double) {
  return 1;
}

template<>
inline double int_pow<1>(double x) {
  return x; // this does a series of muls
}


#endif
