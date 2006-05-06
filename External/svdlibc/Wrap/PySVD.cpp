//-----------------------------------------------
//
// Copyright (c) 2004 Rational Discovery LLC
//   All Rights Reserved
//
//-----------------------------------------------
// """ Code to wrap svdlibc for Python
// """

#ifdef WIN32
#pragma warning (disable: 4786) // warning: long & complicated stl warning
#pragma warning (disable: 4305) // warning: truncation from 'const double' to 'const float'
#pragma warning (disable: 4042) // warning: which has bad storage class
#pragma warning (disable: 4275)
#pragma warning (disable: 4800)
#endif

#include <boost/python.hpp>
#include <RDGeneral/types.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <RDBoost/PySequenceHolder.h>
#define PY_ARRAY_UNIQUE_SYMBOL pySVD_array_API
#include "Numeric/arrayobject.h"
extern "C" {
#include <svdlibc/svdlib.h>
};

namespace python = boost::python;

static long imin(long a, long b) {return (a < b) ? a : b;}

PyObject *GetSingularValues(SMat sparseMat,int iterations,int nVects,
			    double las2end[2],double kappa,bool returnU){
  int i,j;
  PyObject *res=0;
  SVDVerbosity = 0;
  //SVDRec svdRes=svdLAS2(sparseMat,nVects,iterations,las2end,kappa);
  SVDRec svdRes=svdLAS2(sparseMat,iterations,nVects,las2end,kappa);
  if(svdRes){
    res = PyTuple_New(3);
    int dims[2];
    PyArrayObject *tmp;
    double *d;

    dims[0] = svdRes->Vt->rows;
    dims[1] = svdRes->Vt->cols;
    tmp = (PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE);
    d =(double *)tmp->data;
    for(i=0;i<dims[0];i++){
      for(j=0;j<dims[1];j++){
	d[i*dims[1]+j] = svdRes->Vt->value[i][j];
      }
    }
    PyTuple_SetItem(res,0,(PyObject *)tmp);

    tmp = (PyArrayObject *)PyArray_FromDims(1,dims,PyArray_DOUBLE);
    memcpy(tmp->data,(void *)svdRes->S,dims[0]*sizeof(double));
    PyTuple_SetItem(res,1,(PyObject *)tmp);

    if(returnU){
      dims[0] = svdRes->Ut->rows;
      dims[1] = svdRes->Ut->cols;
      //std::cout.precision(3);
      tmp = (PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE);
      d =(double *)tmp->data;
      for(i=0;i<dims[0];i++){
	for(j=0;j<dims[1];j++){
	  d[i*dims[1]+j] = svdRes->Ut->value[i][j];
	  //std::cout << svdRes->Ut->value[i][j] << " ";
	}
	//std::cout << std::endl;
      }
      PyTuple_SetItem(res,2,(PyObject *)tmp);
    } else {
      PyTuple_SetItem(res,2,Py_None);
    }
  } else {
    res = PyTuple_New(3);
    PyTuple_SetItem(res,0,Py_None);
    PyTuple_SetItem(res,1,Py_None);
    PyTuple_SetItem(res,2,Py_None);
  }
  return res;
}


PyObject *SparseSVD(python::object indices,python::object data,
		    int nRows,int nCols,int nVects=-1,
		    bool returnU=false) {
  // -----------------
  // defaults:
  int minV = imin(nRows,nCols);
  if(nVects<0) nVects=minV;
  if(nVects>minV)
    throw_value_error("Number of values requested exceeds number available");
  double las2end[2] = {-1.0e-30, 1.0e-30};
  double kappa = 1e-6;
  //int iterations = imin(imin(nRows,nCols),nVects*4);
  int iterations = imin(nRows,nCols);

  PySequenceHolder<double> vals(data);
  int totSize=vals.size();
  int i,j;

  // NOTE the transpose here:
  SMat sparseMat = svdNewSMat(nCols,nRows,totSize);
  int nSeen = 0;
  for(i=0;i<nRows;i++) {
    python::object tpl=indices[i];
    PySequenceHolder<int> vect(tpl);
    sparseMat->pointr[i]= nSeen;
    for(j=0;j<vect.size();j++) {
      sparseMat->rowind[nSeen]=vect[j];
      sparseMat->value[nSeen]=vals[nSeen];
      nSeen++;
    }
  }
  sparseMat->pointr[nRows]= totSize;
  
  PyObject *res=GetSingularValues(sparseMat,iterations,nVects,las2end,kappa,
				  returnU);
  // FIX: is this required?
  Py_XINCREF(res);
  return res;
};


PyObject *SparseSVDBitMatrix(python::object indices,
			     int nRows,int nCols,int nVects=-1,
			     bool returnU=false) {
  // -----------------
  // defaults:
  int minV = imin(nRows,nCols);
  if(nVects<0) nVects=minV;
  if(nVects>minV)
    throw_value_error("Number of values requested exceeds number available");
  double las2end[2] = {-1.0e-30, 1.0e-30};
  double kappa = 1e-6;
  int iterations = imin(nRows,nCols);

  int i,j;
  int totSize=0;
  for(i=0;i<nRows;i++){
    totSize += python::extract<int>(indices[i].attr("__len__")());
  }

  // NOTE the transpose here:
  SMat sparseMat = svdNewSMat(nCols,nRows,totSize);
  int nSeen = 0;
  for(i=0;i<nRows;i++) {
    python::object tpl=indices[i];
    PySequenceHolder<int> vect(tpl);
    sparseMat->pointr[i]= nSeen;
    for(j=0;j<vect.size();j++) {
      sparseMat->rowind[nSeen]=vect[j];
      sparseMat->value[nSeen]=1;
      nSeen++;
    }
  }
  sparseMat->pointr[nRows]= totSize;
  
  PyObject *res=GetSingularValues(sparseMat,iterations,nVects,las2end,kappa,
				  returnU);
  // FIX: is this required?
  Py_XINCREF(res);
  return res;
};


BOOST_PYTHON_FUNCTION_OVERLOADS(sparseSVD_ol,SparseSVD, 4,6)
BOOST_PYTHON_FUNCTION_OVERLOADS(sparseSVDBits_ol,SparseSVDBitMatrix, 3,5)

  std::string docString;
BOOST_PYTHON_MODULE(cSVD)
{
  import_array();


  python::scope().attr("__doc__") =
    "Module containing fast sparse SVD functionality."
    "\n"
    ;

  docString="SVD (or partial SVD) on a sparse matrix\n\
\n\
  ARGUMENTS:\n\
\n\
    - indices: sequence of nRows sequences with indices of nonzero column indices.\n\
\n\
    - data: sequence of doubles with matrix values\n\
\n\
    - nRows: number of rows in the matrix\n\
\n\
    - nCols: number of columns in the matrix\n\
\n\
    - nVects: (optional) number of singular values to return\n\
      Defaults to nRow.\n\
\n\
    - returnU: (optional) if nonzero, the matrix U^T (see module documentation for\n\
      notation) is returned, otherwise U^T is None.  **NOTE:** this matrix can be\n\
      very large.\n\
      Defaults to zero.\n\
\n\
  RETURNS: a 3-tuple of Numeric Arrays:\n\
\n\
     1) the V matrix [nVects x nRows] \n\
\n\
     2) the S vector [nVects] \n\
\n\
     3) the U^T matrix [nVects x nCols] (or None)\n\
\n\
\n";
  python::def("SparseSVD",SparseSVD,
	      sparseSVD_ol(python::args("indices","values","nRows","nCols","nVects",
					"returnU"),
			   docString.c_str()));
  docString="SVD (or partial SVD) on a sparse matrix of bit vectors\n\
\n\
  ARGUMENTS:\n\
\n\
    - indices: sequence of nRows sequences with indices of nonzero column indices.\n\
\n\
    - nRows: number of rows in the matrix\n\
\n\
    - nCols: number of columns in the matrix\n\
\n\
    - nVects: (optional) number of singular values to return\n\
      Defaults to nRow.\n\
\n\
    - returnU: (optional) if nonzero, the matrix U^T (see module documentation for\n\
      notation) is returned, otherwise U^T is None.  **NOTE:** this matrix can be\n\
      very large.\n\
      Defaults to zero.\n\
\n\
  RETURNS: a 3-tuple of Numeric Arrays:\n\
\n\
     1) the V matrix [nVects x nRows] \n\
\n\
     2) the S vector [nVects] \n\
\n\
     3) the U^T matrix [nVects x nCols] (or None)\n\
\n\
\n";
  python::def("SparseSVDBitMatrix",SparseSVDBitMatrix,
	      sparseSVDBits_ol(python::args("indices","nRows","nCols","nVects",
					    "returnU"),
			       docString.c_str()));



}
