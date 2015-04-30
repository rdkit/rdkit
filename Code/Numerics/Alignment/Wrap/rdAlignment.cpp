// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdalignment_array_API
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include "numpy/arrayobject.h"
#include <RDBoost/import_array.h>

#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/Wrap.h>
#include <Geometry/point.h>

#include <Numerics/Alignment/AlignPoints.h>

namespace python = boost::python;

namespace RDNumeric {
  namespace Alignments {
    void GetPointsFromPythonSequence(python::object &points,
                                     RDGeom::Point3DConstPtrVect &pts){
      PyObject *pyObj=points.ptr();
      unsigned int nrows, ncols;
      double *data;
      if (PyArray_Check(pyObj)) {
        // get the dimensions of the array
        PyArrayObject *ptsMat = reinterpret_cast<PyArrayObject *>(pyObj);
        nrows = ptsMat->dimensions[0];
        ncols = ptsMat->dimensions[1];
        
        if (ncols != 3) 
          throw_value_error("Wrong dimension for the points array");
        
        data = reinterpret_cast<double *>(ptsMat->data);

        for (unsigned int i = 0; i < nrows; i++) {
          RDGeom::Point3D *rpt = new RDGeom::Point3D(data[i*3],
                                                     data[i*3+1],
                                                     data[i*3+2]);
          pts.push_back(rpt);
        }
      } else if(PySequence_Check(pyObj)) {
        nrows = PySequence_Size(pyObj);
        if (nrows <= 0) 
          throw_value_error("Empty sequence passed in");
        python::extract<RDGeom::Point3D> ptOk(points[0]);
        if(!ptOk.check()){
          for (unsigned int i = 0; i < nrows; i++) {
            PySequenceHolder<double> row(points[i]);
            if (row.size() != 3) 
              throw_value_error("Wrong number of entries in the list of lists");
            RDGeom::Point3D *rpt = new RDGeom::Point3D(row[0],row[1],row[2]);
            pts.push_back(rpt);
          }
        }else{
            for (unsigned int i = 0; i < nrows; i++) {
              python::extract<RDGeom::Point3D> pt(points[i]);
              if(pt.check()){
                RDGeom::Point3D *rpt = new RDGeom::Point3D(pt);
                pts.push_back(rpt);
              } else {
                throw_value_error("non-Point3D found in sequence of points");
              }
            }

        }
      } else {
        throw_value_error("non-sequence argument provided");
      }
    }

    PyObject *AlignPointPairs(python::object refPoints, python::object probePoints, 
                              python::object weights=python::list(), bool reflect=false, 
                              unsigned int maxIterations=50) {
      // The reference and probe points can be specifed in two formats
      // 1. A Numeric array of dimensions (N,3) where N is the number of points
      // 2. A list (or tuple) or lists (or tuples)
      // 
      // The similar thing applies to weights
      // 1. Can be a numeric vector of size N
      // 2. A list of doubles of size N
      
      // first deal with situation where we have Numerics arrays
      RDGeom::Point3DConstPtrVect refPts, probePts;

      GetPointsFromPythonSequence(refPoints,refPts);
      GetPointsFromPythonSequence(probePoints,probePts);
      
      unsigned int npt = refPts.size();
      if (npt != probePts.size()) 
        throw_value_error("Mis-match in number of points");
      
      PyObject *weightsObj = weights.ptr();
      RDNumeric::DoubleVector *wtsVec;
      wtsVec = 0;
      double *data;
      if (PyArray_Check(weightsObj)) {
        PyArrayObject *wtsMat = reinterpret_cast<PyArrayObject *>(weightsObj);
        unsigned int nwts = wtsMat->dimensions[0];
        if (nwts != npt) 
          throw_value_error("Number of weights supplied do not match the number of points");
        wtsVec = new RDNumeric::DoubleVector(nwts);
        data = reinterpret_cast<double *>(wtsMat->data);
        for (unsigned int i = 0; i < nwts; i++) {
          wtsVec->setVal(i, data[i]);
        }
      } else {
        PySequenceHolder<double> wts(weights);
        unsigned int nwts = wts.size();
        
        if (nwts > 0) {
          if (nwts != npt) 
            throw_value_error("Number of weights supplied do not match the number of points");
          wtsVec = new RDNumeric::DoubleVector(nwts);
          for (unsigned int i = 0; i < npt; i++) {
            wtsVec->setVal(i, wts[i]);
          }
        }
      }

      RDGeom::Transform3D trans;
      double ssd = AlignPoints(refPts, probePts, trans, wtsVec, reflect, maxIterations);
      
      npy_intp dims[2];
      dims[0] = 4;
      dims[1] = 4;
      PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
      double *resData=reinterpret_cast<double *>(res->data);
      const double *tdata = trans.getData();
      for(unsigned int i=0; i < trans.numRows(); ++i){
        unsigned int itab = i*4;
        for (unsigned int j = 0; j < trans.numRows(); ++j) {
          resData[itab + j] = tdata[itab+j];
        }
      }
      
      delete wtsVec;

      for (unsigned int i = 0; i < npt; ++i) {
        delete probePts[i];
        delete refPts[i];
      }
      PyObject *resTup = PyTuple_New(2);
      PyObject *ssdItem = PyFloat_FromDouble(ssd);
      PyTuple_SetItem(resTup,0,ssdItem);
      PyTuple_SetItem(resTup,1,PyArray_Return(res));
      return resTup;
    }
  }
}

BOOST_PYTHON_MODULE(rdAlignment) {

  rdkit_import_array();
  python::scope().attr("__doc__") =
    "Module containing functions to align pairs of points in 3D"
    ;

  std::string docString = "Compute the optimal alignment (minimum RMSD) between two set of points \n\n\
 \n\
 ARGUMENTS:\n\n\
    - refPoints : reference points sepcified as a N by 3 Numeric array or \n\
                  sequence of 3-sequences or sequence of Point3Ds \n\
    - probePoints : probe points to align to reference points - same format \n\
                  restrictions as reference points apply here \n\
    - weights : optional numeric vector or list of weights to associate to each pair of points\n\
    - reflect : reflect the probe points before attempting alignment\n\
    - maxIteration : maximum number of iterations to try to minimize RMSD \n\
                  \n\
 RETURNS:\n\n\
    a 2-tuple:\n\
      - SSD value for the alignment\n\
      - the 4x4 transform matrix, as a Numeric array\n\
\n";
  python::def("GetAlignmentTransform", RDNumeric::Alignments::AlignPointPairs,
              (python::arg("refPoints"), python::arg("probePoints"),
               python::arg("weights")=python::list(),
               python::arg("reflect")=false, python::arg("maxIterations")=50),
              docString.c_str());
}

