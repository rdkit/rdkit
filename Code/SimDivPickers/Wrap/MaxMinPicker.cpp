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


#include <SimDivPickers/DistPicker.h>
#include <SimDivPickers/MaxMinPicker.h>
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
    pyobjFunctor(python::object obj) : dp_obj(obj) {}
    double operator()(unsigned int i,unsigned int j) {
      double res;
      std::pair<unsigned int ,unsigned int> idxPair(i,j);
      if(this->d_cache.count(idxPair)>0){
        res = this->d_cache[idxPair];
      } else {
        res=python::extract<double>(dp_obj(i,j));
        this->d_cache[idxPair]=res;
      }
      return res;
    }
  private:
    python::object dp_obj;
    std::map<std::pair<unsigned int,unsigned int>,double> d_cache;
  };

  RDKit::INT_VECT LazyMaxMinPicks(MaxMinPicker *picker, 
                                  python::object distFunc,
                                  int poolSize, 
                                  int pickSize,
                                  python::object firstPicks,
                                  int seed) {
    pyobjFunctor functor(distFunc);
    RDKit::INT_VECT firstPickVect;
    for(unsigned int i=0;i<python::extract<unsigned int>(firstPicks.attr("__len__")());++i){
      firstPickVect.push_back(python::extract<int>(firstPicks[i]));
    }
    RDKit::INT_VECT res=picker->lazyPick(functor, poolSize, pickSize,firstPickVect,seed);
    return res;
  }

  struct MaxMin_wrap {
    static void wrap() {
      python::class_<MaxMinPicker>("MaxMinPicker", 
                                   "A class for diversity picking of items using the MaxMin Algorithm\n")
        .def("Pick", MaxMinPicks,
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
             "  - seed: (optional) seed for the random number genrator\n"
             )

        .def("LazyPick", LazyMaxMinPicks,
             (python::arg("self"),python::arg("distFunc"),python::arg("poolSize"),
              python::arg("pickSize"),python::arg("firstPicks")=python::tuple(),
              python::arg("seed")=-1),
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
             "  - seed: (optional) seed for the random number genrator\n"
             )
        ;
    };
  };
}

void wrap_maxminpick() {
  RDPickers::MaxMin_wrap::wrap();
}

  
