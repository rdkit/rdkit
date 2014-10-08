// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmetric_array_API
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include "numpy/arrayobject.h"

#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

#include <RDGeneral/types.h>

#include <DataManip/MetricMatrixCalc/MetricMatrixCalc.h>
#include <DataManip/MetricMatrixCalc/MetricFuncs.h>
#include <DataStructs/BitVects.h>
#include <string>

using namespace RDDataManip;

void wrap_MMcalc();

namespace python = boost::python;
namespace RDDataManip {
  
  PyObject *getEuclideanDistMat(python::object descripMat) {
    // Bit of a pain involved here, we accept three types of PyObjects here
    // 1. A Numeric Array 
    //     - first find what 'type' of entry we have (float, double and int is all we recognize for now)
    //     - then point to contiguous piece of memory from the array that contains the data with a type*
    //     - then make a new type** pointer so that double index into this contiguous memory will work
    //       and then pass it along to the distance calculator
    // 2. A list of Numeric Vector (or 1D arrays) 
    //     - in this case wrap descripMat with a PySequenceHolder<type*> where type is the 
    //       type of entry in vector (accepted types are int, double and float
    //     - Then pass the PySequenceHolder to the metrci calculator
    // 3. A list (or tuple) of lists (or tuple)
    //     - In this case other than wrapping descripMat with a PySequenceHolder
    //       each of the indivual list in there are also wrapped by a PySequenceHolder
    //     - so the distance calculator is passed in a "PySequenceHolder<PySequenceHolder<double>>"
    //     - FIX: not that we always convert entry values to double here, even if we passed
    //       in a list of list of ints (or floats). Given that lists can be heterogeneous, I do not
    //       know how to ask a list what type of entries if contains.
    //
    //  OK my brain is going to explode now
    
    // first deal with situation where we have an Numeric Array
    PyObject *descMatObj = descripMat.ptr();
    PyArrayObject *distRes;
    if (PyArray_Check(descMatObj)) {
      // get the dimensions of the array
      int nrows = ((PyArrayObject *)descMatObj)->dimensions[0];
      int ncols = ((PyArrayObject *)descMatObj)->dimensions[1];
      int i;
      CHECK_INVARIANT((nrows > 0) && (ncols > 0), "");

      npy_intp dMatLen = nrows*(nrows-1)/2;
      
      // now that we have the dimensions declare the distance matrix which is always a 
      // 1D double array
      distRes = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
      
      // grab a pointer to the data in the array so that we can directly put values in there
      // and avoid copying : 
      double *dMat = (double *)distRes->data;

      PyArrayObject *copy;
      copy = (PyArrayObject *)PyArray_ContiguousFromObject(descMatObj, 
							   ((PyArrayObject *)descMatObj)->descr->type_num,
							   2,2);
      // if we have double array
      if (((PyArrayObject *)descMatObj)->descr->type_num == NPY_DOUBLE) {
        double *desc = (double *)copy->data;
        
	// REVIEW: create an adaptor object to hold a double * and support
	//  operator[]() so that we don't have to do this stuff:

        // here is the 2D array trick this so that when the distance calaculator
        // asks for desc2D[i] we basically get the ith row as double*
        double **desc2D = new double*[nrows];
        for (i = 0; i < nrows; i++) {
          desc2D[i] = desc;
          desc += ncols;
        }
        MetricMatrixCalc<double**, double*> mmCalc;
        mmCalc.setMetricFunc(&EuclideanDistanceMetric<double *, double *>);
        mmCalc.calcMetricMatrix(desc2D, nrows, ncols, dMat);
        
        delete [] desc2D;
        // we got the distance matrix we are happy so return
        return PyArray_Return(distRes);
      }
      
      // if we have a float array
      else if (((PyArrayObject *)descMatObj)->descr->type_num == NPY_FLOAT) {
        float* desc = (float *)copy->data;
        float **desc2D = new float*[nrows];
        for (i = 0; i < nrows; i++) {
          desc2D[i] = desc;
          desc += ncols;
        }
        MetricMatrixCalc<float**, float*> mmCalc;
        mmCalc.setMetricFunc(&EuclideanDistanceMetric<float *, float*>);
        mmCalc.calcMetricMatrix(desc2D, nrows, ncols, dMat);
        delete [] desc2D;
        return PyArray_Return(distRes);
      }

      // if we have an interger array
      else if (((PyArrayObject *)descMatObj)->descr->type_num == NPY_INT) {
        int *desc = (int *)copy->data;
        int **desc2D = new int*[nrows];
        for (i = 0; i < nrows; i++) {
          desc2D[i] = desc;
          desc += ncols;
        }
        MetricMatrixCalc<int**, int*> mmCalc;
        mmCalc.setMetricFunc(&EuclideanDistanceMetric<int *, int*>);
        mmCalc.calcMetricMatrix(desc2D, nrows, ncols, dMat);
        delete [] desc2D;
        return PyArray_Return(distRes);
      }
      else {
        // unreconiged type for the matrix, throw up
        throw_value_error("The array has to be of type int, float, or double for GetEuclideanDistMat");
      }
    } // done with an array input
    else {
      // REVIEW: removed a ton of code here

      // we have probably have a list or a tuple
      
      unsigned int ncols = 0;
      unsigned int nrows = python::extract<unsigned int>(descripMat.attr("__len__")());
      CHECK_INVARIANT(nrows > 0, "Empty list passed in");

      npy_intp dMatLen = nrows*(nrows-1)/2;
      distRes = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
      double *dMat = (double *)distRes->data;
      
      // assume that we a have a list of list of values (that can be extracted to double)
      std::vector<PySequenceHolder<double> > dData;
      dData.reserve(nrows);
      for (unsigned int i = 0; i < nrows; i++) {
	//PySequenceHolder<double> row(seq[i]);
        PySequenceHolder<double> row(descripMat[i]);
        if(i==0){
	  ncols = row.size();
	} else if( row.size() != ncols ){
	  throw_value_error("All subsequences must be the same length");
	}
	dData.push_back(row);
      }
      
      MetricMatrixCalc< std::vector<PySequenceHolder<double> >, PySequenceHolder<double> > mmCalc;
      mmCalc.setMetricFunc(&EuclideanDistanceMetric< PySequenceHolder<double>, PySequenceHolder<double> >);
      mmCalc.calcMetricMatrix(dData, nrows, ncols, dMat);
    }
    return PyArray_Return(distRes);
  }
  
  PyObject *getTanimotoDistMat(python::object bitVectList) {
    // we will assume here that we have a either a list of ExplicitBitVectors or
    // SparseBitVects
    int nrows = python::extract<int>(bitVectList.attr("__len__")());
    CHECK_INVARIANT(nrows > 1, "");

    // First check what type of vector we have
    python::object v1 = bitVectList[0];
    python::extract<ExplicitBitVect> ebvWorks(v1);
    python::extract<SparseBitVect> sbvWorks(v1);
    if(!ebvWorks.check() && !sbvWorks.check()){
      throw_value_error("GetTanimotoDistMat can only take a sequence of ExplicitBitVects or SparseBitvects");
    }

    npy_intp dMatLen = nrows*(nrows-1)/2;
    PyArrayObject *simRes = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
    double *sMat = (double *)simRes->data;
    
    if (ebvWorks.check()) {
      PySequenceHolder<ExplicitBitVect> dData(bitVectList);
      MetricMatrixCalc<PySequenceHolder<ExplicitBitVect>, ExplicitBitVect> mmCalc;
      mmCalc.setMetricFunc(&TanimotoDistanceMetric<ExplicitBitVect, ExplicitBitVect>);
      mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
    }
    else if (sbvWorks.check()) {
      PySequenceHolder<SparseBitVect> dData(bitVectList);
      MetricMatrixCalc<PySequenceHolder<SparseBitVect>, SparseBitVect> mmCalc;
      mmCalc.setMetricFunc(&TanimotoDistanceMetric<SparseBitVect, SparseBitVect>);
      mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
    }
    return PyArray_Return(simRes);
  }

  PyObject *getTanimotoSimMat(python::object bitVectList) {
    // we will assume here that we have a either a list of ExplicitBitVectors or
    // SparseBitVects
    int nrows = python::extract<int>(bitVectList.attr("__len__")());
    CHECK_INVARIANT(nrows > 1, "");

    // First check what type of vector we have
    python::object v1 = bitVectList[0];
    python::extract<ExplicitBitVect> ebvWorks(v1);
    python::extract<SparseBitVect> sbvWorks(v1);
    if(!ebvWorks.check() && !sbvWorks.check()){
      throw_value_error("GetTanimotoDistMat can only take a sequence of ExplicitBitVects or SparseBitvects");
    }

    npy_intp dMatLen = nrows*(nrows-1)/2;
    PyArrayObject *simRes = (PyArrayObject *)PyArray_SimpleNew(1, &dMatLen, NPY_DOUBLE);
    double *sMat = (double *)simRes->data;
    
    if (ebvWorks.check()) {
      PySequenceHolder<ExplicitBitVect> dData(bitVectList);
      MetricMatrixCalc<PySequenceHolder<ExplicitBitVect>, ExplicitBitVect> mmCalc;
      mmCalc.setMetricFunc(&TanimotoSimilarityMetric<ExplicitBitVect, ExplicitBitVect>);
      mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
    }
    else if (sbvWorks.check()) {
      PySequenceHolder<SparseBitVect> dData(bitVectList);
      MetricMatrixCalc<PySequenceHolder<SparseBitVect>, SparseBitVect> mmCalc;
      mmCalc.setMetricFunc(&TanimotoSimilarityMetric<SparseBitVect, SparseBitVect>);
      mmCalc.calcMetricMatrix(dData, nrows, 0, sMat);
    }
    return PyArray_Return(simRes);
  }
}
                                            
BOOST_PYTHON_MODULE(rdMetricMatrixCalc) 
{
  python::scope().attr("__doc__") =
    "Module containing the calculator for metric matrix calculation, \n"
    "e.g. similarity and distance matrices"
    ;

  rdkit_import_array();
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  
  std::string docString;
  docString = "Compute the distance matrix from a descriptor matrix using the Euclidean distance metric\n\n\
  ARGUMENTS: \n\
\n\
    descripMat - A python object of any one of the following types \n\
                   1. A numeric array of dimensions n by m where n is the number of items in the data set \n\
                       and m is the number of descriptors \n\
                   2. A list of Numeric Vectors (or 1D arrays), each entry in the list corresponds \n\
                       to descriptor vector for one item \n\
                   3. A list (or tuple) of lists (or tuples) of values, where the values can be extracted to \n\
                       double. \n\n\
  RETURNS: \n\
    A numeric one-dimensional array containing the lower triangle elements of the symmetric distance matrix\n\n";
  python::def("GetEuclideanDistMat", RDDataManip::getEuclideanDistMat, 
              docString.c_str());

  docString = "Compute the distance matrix from a list of BitVects using the Tanimoto distance metric\n\n\
  ARGUMENTS: \n\
\n\
    bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors, \n\
                  needs to be expanded to support a list of SparseBitVects\n\n\
  RETURNS: \n\
    A numeric 1 dimensional array containing the lower triangle elements of the\n\
    symmetric distance matrix\n\n";
  python::def("GetTanimotoDistMat", RDDataManip::getTanimotoDistMat,
              docString.c_str());
  
  docString = "Compute the similarity matrix from a list of BitVects \n\n\
  ARGUMENTS: \n\
\n\
    bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors, \n\
                  needs to be expanded to support a list of SparseBitVects\n\n\
  RETURNS: \n\
    A numeric 1 dimensional array containing the lower triangle elements of the symmetric similarity matrix\n\n";
  python::def("GetTanimotoSimMat", RDDataManip::getTanimotoSimMat,
              docString.c_str());
}
