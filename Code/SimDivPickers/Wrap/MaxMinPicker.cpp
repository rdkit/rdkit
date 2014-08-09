// $Id$
//
//  Copyright (C) 2003-2008  Greg Landrum and Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY

#define PY_ARRAY_UNIQUE_SYMBOL rdpicker_array_API
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <boost/python/numeric.hpp>
#include "numpy/oldnumeric.h"
#include <map>

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>
#include <SimDivPickers/HierarchicalClusterPicker.h>
#include <iostream>

namespace python = boost::python;
namespace RDPickers {
  
  // REVIEW: the poolSize can be pulled from the numeric array
  RDKit::INT_VECT MaxMinPicks(MaxMinPicker *picker, 
                              python::object distMat,
                              int poolSize, 
                              int pickSize,
                              python::object firstPicks,
                              int seed) {
    if(pickSize>=poolSize){
      throw ValueErrorException("pickSize must be less than poolSize");
    }

    if (!PyArray_Check(distMat.ptr())){
      throw ValueErrorException("distance mat argument must be a numpy matrix");
    }

    PyArrayObject *copy;
    copy = (PyArrayObject *)PyArray_ContiguousFromObject(distMat.ptr(), 
							 PyArray_DOUBLE, 1,1);
    double *dMat = (double *)copy->data;
    
    RDKit::INT_VECT firstPickVect;
    for(unsigned int i=0;i<python::extract<unsigned int>(firstPicks.attr("__len__")());++i){
      firstPickVect.push_back(python::extract<int>(firstPicks[i]));
    }
    RDKit::INT_VECT res=picker->pick(dMat, poolSize, pickSize,firstPickVect,seed);
    Py_DECREF(copy);
    return res;
  }
                        
  class pyobjFunctor {
  public:
    pyobjFunctor(python::object obj,bool useCache) : dp_obj(obj), dp_cache(NULL) {
      if(useCache) dp_cache= new std::map<std::pair<unsigned int,unsigned int>,double>();
    }
    ~pyobjFunctor() {
      delete dp_cache;
    }
    double operator()(unsigned int i,unsigned int j) {
      double res;
      std::pair<unsigned int ,unsigned int> idxPair(i,j);
      if(dp_cache && dp_cache->count(idxPair)>0){
        res = (*dp_cache)[idxPair];
      } else {
        res=python::extract<double>(dp_obj(i,j));
        if(dp_cache)
          (*dp_cache)[idxPair]=res;
      }
      return res;
    }
  private:
    python::object dp_obj;
    std::map<std::pair<unsigned int,unsigned int>,double> *dp_cache;
  };

  RDKit::INT_VECT LazyMaxMinPicks(MaxMinPicker *picker, 
                                  python::object distFunc,
                                  int poolSize, 
                                  int pickSize,
                                  python::object firstPicks,
                                  int seed,
                                  bool useCache) {
    RDKit::INT_VECT firstPickVect;
    for(unsigned int i=0;i<python::extract<unsigned int>(firstPicks.attr("__len__")());++i){
      firstPickVect.push_back(python::extract<int>(firstPicks[i]));
    }
    RDKit::INT_VECT res;
    pyobjFunctor functor(distFunc,useCache);
    res=picker->lazyPick(functor, poolSize, pickSize,firstPickVect,seed);
    return res;
  }

  // NOTE: TANIMOTO and DICE provably return the same results for the diversity picking
  //    this is still here just in case we ever later want to support other methods.
  typedef enum {
    TANIMOTO=1,
    DICE
  } DistanceMethod;

  template <typename BV>
  class pyBVFunctor {
  public:
    pyBVFunctor(const std::vector<const BV *> &obj,DistanceMethod method,bool useCache) : d_obj(obj), d_method(method), dp_cache(NULL) {
      if(useCache) dp_cache = new std::map<std::pair<unsigned int,unsigned int>,double>();
    }
    ~pyBVFunctor() {
      delete dp_cache;
    }
    double operator()(unsigned int i,unsigned int j) {
      double res=0.0;
      std::pair<unsigned int ,unsigned int> idxPair(i,j);
      if(dp_cache && dp_cache->count(idxPair)>0){
        res = (*dp_cache)[idxPair];
      } else {
        switch(d_method){
        case TANIMOTO:
          res = 1.-TanimotoSimilarity(*d_obj[i],*d_obj[j]);
          break;
        case DICE:
          res = 1.-DiceSimilarity(*d_obj[i],*d_obj[j]);
          break;
        default:
          throw_value_error("unsupported similarity value");
        }
        if(dp_cache){
          (*dp_cache)[idxPair]=res;
        }
      }
      return res;
    }
  private:
    const std::vector<const BV *> &d_obj;
    DistanceMethod d_method;
    std::map<std::pair<unsigned int,unsigned int>,double> *dp_cache;
  };
  
  RDKit::INT_VECT LazyVectorMaxMinPicks(MaxMinPicker *picker, 
                                        python::object objs,
                                        int poolSize, 
                                        int pickSize,
                                        python::object firstPicks,
                                        int seed,
                                        bool useCache
                                        ) {
    std::vector<const ExplicitBitVect *> bvs(poolSize);
    for(unsigned int i=0;i<poolSize;++i){
      bvs[i]=python::extract<const ExplicitBitVect *>(objs[i]);
    }
    pyBVFunctor<ExplicitBitVect> functor(bvs,TANIMOTO,useCache);
    RDKit::INT_VECT firstPickVect;
    for(unsigned int i=0;i<python::extract<unsigned int>(firstPicks.attr("__len__")());++i){
      firstPickVect.push_back(python::extract<int>(firstPicks[i]));
    }
    RDKit::INT_VECT res=picker->lazyPick(functor, poolSize, pickSize,firstPickVect,seed);
    return res;
  }
} // end of namespace RDPickers

  struct MaxMin_wrap {
    static void wrap() {
      python::class_<RDPickers::MaxMinPicker>("MaxMinPicker", 
                                   "A class for diversity picking of items using the MaxMin Algorithm\n")
        .def("Pick", RDPickers::MaxMinPicks,
             (python::arg("self"),python::arg("distMat"),python::arg("poolSize"),
              python::arg("pickSize"),python::arg("firstPicks")=python::tuple(),
              python::arg("seed")=-1),
             "Pick a subset of items from a pool of items using the MaxMin Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 \n\n"
             "ARGUMENTS:\n"
             "  - distMat: 1D distance matrix (only the lower triangle elements)\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - firstPicks: (optional) the first items to be picked (seeds the list)\n"
             "  - seed: (optional) seed for the random number generator\n"
             )

        .def("LazyPick", RDPickers::LazyMaxMinPicks,
             (python::arg("self"),python::arg("distFunc"),python::arg("poolSize"),
              python::arg("pickSize"),python::arg("firstPicks")=python::tuple(),
              python::arg("seed")=-1,python::arg("useCache")=true),
             "Pick a subset of items from a pool of items using the MaxMin Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 \n"
             "ARGUMENTS:\n\n"
             "  - distFunc: a function that should take two indices and return the\n"
             "              distance between those two points.\n"
             "              NOTE: the implementation caches distance values, so the\n"
             "              client code does not need to do so; indeed, it should not.\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - firstPicks: (optional) the first items to be picked (seeds the list)\n"
             "  - seed: (optional) seed for the random number generator\n"
             "  - useCache: (optional) toggles use of a cache for the distance calculation\n"
             "              This trades memory usage for speed.\n"
             )
        .def("LazyBitVectorPick", RDPickers::LazyVectorMaxMinPicks,
             (python::arg("self"),python::arg("objects"),python::arg("poolSize"),
              python::arg("pickSize"),python::arg("firstPicks")=python::tuple(),
              python::arg("seed")=-1,
              python::arg("useCache")=true),
             "Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm\n"
             "Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 \n"
             "ARGUMENTS:\n\n"
             "  - vectors: a sequence of the bit vectors that should be picked from.\n"
             "  - poolSize: number of items in the pool\n"
             "  - pickSize: number of items to pick from the pool\n"
             "  - firstPicks: (optional) the first items to be picked (seeds the list)\n"
             "  - seed: (optional) seed for the random number generator\n"
             "  - useCache: (optional) toggles use of a cache for the distance calculation\n"
             "              This trades memory usage for speed.\n"

             )
        ;
    };
  };


void wrap_maxminpick() {
  MaxMin_wrap::wrap();
}

  
