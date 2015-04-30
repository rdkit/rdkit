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

#include <boost/python.hpp>

#define PY_ARRAY_UNIQUE_SYMBOL DistGeom_array_API
#include "numpy/arrayobject.h"
#include <RDBoost/Wrap.h>
#include <RDBoost/pyint_api.h>
#include <RDBoost/import_array.h>

#include <Geometry/point.h>
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/ChiralSet.h>

#include <ForceField/ForceField.h>

#include <vector>
#include <map>

namespace python = boost::python;

namespace RDKit {
  bool doTriangleSmoothing(python::object boundsMatArg,double tol){
    PyObject *boundsMatObj = boundsMatArg.ptr();
    if(!PyArray_Check(boundsMatObj))
      throw_value_error("Argument isn't an array");
    
    PyArrayObject *boundsMat=reinterpret_cast<PyArrayObject *>(boundsMatObj);
    // get the dimensions of the array
    int nrows = boundsMat->dimensions[0];
    int ncols = boundsMat->dimensions[1];
    if(nrows!=ncols)
      throw_value_error("The array has to be square");
    if(nrows<=0)
      throw_value_error("The array has to have a nonzero size");
    if (boundsMat->descr->type_num != PyArray_DOUBLE)
      throw_value_error("Only double arrays are currently supported");

    unsigned int dSize = nrows*nrows;
    double *cData = new double[dSize];
    double *inData = reinterpret_cast<double *>(boundsMat->data);
    memcpy(static_cast<void *>(cData), 
           static_cast<const void *>(inData),
           dSize*sizeof(double));
    DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
    DistGeom::BoundsMatrix bm(nrows,sdata);

    bool res=DistGeom::triangleSmoothBounds(&bm,tol);
    memcpy(static_cast<void *>(inData), 
           static_cast<const void *>(cData),
           dSize*sizeof(double));
    return res;
  }

  PyObject *embedBoundsMatrix(python::object boundsMatArg,int maxIters=10,
                              bool randomizeOnFailure=false,int numZeroFail=2,
                              python::list weights=python::list(),
                              int randomSeed=-1){
    PyObject *boundsMatObj = boundsMatArg.ptr();
    if(!PyArray_Check(boundsMatObj))
      throw_value_error("Argument isn't an array");
    
    PyArrayObject *boundsMat=reinterpret_cast<PyArrayObject *>(boundsMatObj);
    // get the dimensions of the array
    unsigned int nrows = boundsMat->dimensions[0];
    unsigned int ncols = boundsMat->dimensions[1];
    if(nrows!=ncols)
      throw_value_error("The array has to be square");
    if(nrows<=0)
      throw_value_error("The array has to have a nonzero size");
    if (boundsMat->descr->type_num != PyArray_DOUBLE)
      throw_value_error("Only double arrays are currently supported");

    unsigned int dSize = nrows*nrows;
    double *cData = new double[dSize];
    double *inData = reinterpret_cast<double *>(boundsMat->data);
    memcpy(static_cast<void *>(cData), 
           static_cast<const void *>(inData),
           dSize*sizeof(double));

    DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
    DistGeom::BoundsMatrix bm(nrows,sdata);

    RDGeom::Point3D *positions=new RDGeom::Point3D[nrows];
    std::vector<RDGeom::Point *> posPtrs;
    for (unsigned int i = 0; i < nrows; i++) {
      posPtrs.push_back(&positions[i]);
    }

    RDNumeric::DoubleSymmMatrix distMat(nrows, 0.0);

    // ---- ---- ---- ---- ---- ---- ---- ---- ---- 
    // start the embedding:
    bool gotCoords=false;
    for(int iter=0;iter<maxIters && !gotCoords;iter++){
      // pick a random distance matrix
      DistGeom::pickRandomDistMat(bm,distMat,randomSeed);

      // and embed it:
      gotCoords=DistGeom::computeInitialCoords(distMat,posPtrs,randomizeOnFailure,
                                               numZeroFail,randomSeed);

      // update the seed:
      if(randomSeed>=0) randomSeed+=iter*999;
    }

    if(gotCoords){
      std::map<std::pair<int,int>,double> weightMap;
      unsigned int nElems=PySequence_Size(weights.ptr());
      for(unsigned int entryIdx=0;entryIdx<nElems;entryIdx++){
        PyObject *entry=PySequence_GetItem(weights.ptr(),entryIdx);
        if(!PySequence_Check(entry) || PySequence_Size(entry)!=3){
          throw_value_error("weights argument must be a sequence of 3-sequences");
        }
        int idx1=PyInt_AsLong(PySequence_GetItem(entry,0));
        int idx2=PyInt_AsLong(PySequence_GetItem(entry,1));
        double w=PyFloat_AsDouble(PySequence_GetItem(entry,2));
        weightMap[std::make_pair(idx1,idx2)]=w;
      }
      DistGeom::VECT_CHIRALSET csets;
      ForceFields::ForceField *field = DistGeom::constructForceField(bm,posPtrs,csets,0.0, 0.0,
                                                                     &weightMap);
      CHECK_INVARIANT(field,"could not build dgeom force field");
      field->initialize();
      if(field->calcEnergy()>1e-5){
        int needMore=1;
        while(needMore){
          needMore=field->minimize();
        }
      }
      delete field;
    } else {
      throw_value_error("could not embed matrix");
    }

    // ---- ---- ---- ---- ---- ---- ---- ---- ---- 
    // construct the results matrix:
    npy_intp dims[2];
    dims[0] = nrows;
    dims[1] = 3;
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    double *resData=reinterpret_cast<double *>(res->data);
    for(unsigned int i=0;i<nrows;i++){
      unsigned int iTab=i*3;
      for (unsigned int j = 0; j < 3; ++j) {
        resData[iTab + j]=positions[i][j]; //.x;
      }
    }
    delete [] positions;

    return PyArray_Return(res);
  }
}

BOOST_PYTHON_MODULE(DistGeom) {
  python::scope().attr("__doc__") =
    "Module containing functions for basic distance geometry operations"
    ;

  rdkit_import_array();
  python::register_exception_translator<ValueErrorException>(&translate_value_error);

  std::string docString;

  docString = "Do triangle smoothing on a bounds matrix\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mat: a square Numeric array of doubles containing the bounds matrix, this matrix\n\
           *is* modified by the smoothing\n\
 \n\
 RETURNS:\n\n\
    a boolean indicating whether or not the smoothing worked.\n\
\n";
  python::def("DoTriangleSmoothing", RDKit::doTriangleSmoothing,
              (python::arg("boundsMatrix"),python::arg("tol")=0.),
              docString.c_str());


  docString = "Embed a bounds matrix and return the coordinates\n\n\
 \n\
 ARGUMENTS:\n\n\
    - boundsMatrix: a square Numeric array of doubles containing the bounds matrix, this matrix\n\
           should already be smoothed\n\
    - maxIters: (optional) the maximum number of random distance matrices to try\n\
    - randomizeOnFailure: (optional) toggles using random coords if a matrix fails to embed\n\
    - numZeroFail: (optional) sets the number of zero eigenvalues to be considered a failure\n\
    - weights: (optional) a sequence of 3 sequences (i,j,weight) indicating elements of \n\
       the bounds matrix whose weights should be adjusted\n\
    - randomSeed: (optional) sets the random number seed used for embedding\n\
 \n\
 RETURNS:\n\n\
    a Numeric array of doubles with the coordinates\n\
\n";
  python::def("EmbedBoundsMatrix", RDKit::embedBoundsMatrix,
              (python::arg("boundsMatrix"),python::arg("maxIters")=10,
               python::arg("randomizeOnFailure")=false,
               python::arg("numZeroFail")=2,
               python::arg("weights")=python::list(),
               python::arg("randomSeed")=-1),
              docString.c_str());


}
