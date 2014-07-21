// $Id$
//
//  Copyright (C) 2002-2010 Greg Landrum and Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_Clustering

#include <RDBoost/Wrap.h>
#include <boost/cstdint.hpp>

namespace python = boost::python;

#include <numpy/arrayobject.h>

#include <RDBoost/import_array.h>

typedef double real;

extern "C"
void distdriver_(boost::int64_t *n,boost::int64_t *len,
		 real *dists,
		 boost::int64_t *toggle,
		 boost::int64_t *ia,boost::int64_t *ib,real *crit);

//
// Rather than deal with any nonsense like trying to get
// the distance matrix built properly on the f2c side of things
// (thus drowning in the waves of f2c hate), we'll generate
// the distance matrix on our own here and then call distdriver_
//
void clusterit(real *dataP,boost::int64_t n,boost::int64_t m,boost::int64_t iopt,
	       boost::int64_t *ia,boost::int64_t *ib,real *crit){
  real *dists;
  boost::int64_t len;
  boost::int64_t pos = 0;
  boost::int64_t i,j,k,iTab,jTab;
  double tmp;
  len = (n*(n-1))/2;
  dists = (real *)calloc(len,sizeof(real));
  for(i=1;i<n;i++){
    iTab = i*m;
    for(j=0;j<i;j++){
      jTab = j*m;
      for(k=0;k<m;k++){
	tmp = dataP[iTab+k]-dataP[jTab+k];
	dists[pos] += tmp*tmp;
      }
      pos++;
    }
  }
  distdriver_(&n,&len,dists,&iopt,ia,ib,crit);
  free(dists);
};

static PyObject *
Clustering_MurtaghCluster(python::object data, int nPts, int sz, int option)
{
  PyArrayObject *dataContig;
  boost::int64_t *ia,*ib;
  real *crit;
  PyObject *res;
  PyObject *tmp;
  npy_intp dims[2];

  if (PyArray_Check(data.ptr())) {
    dataContig 
      = reinterpret_cast<PyArrayObject *>(PyArray_ContiguousFromObject(data.ptr(),PyArray_DOUBLE,2,2));
  }
  else {
    throw_value_error("PyArray_Type expected as input");
  }

  ia = (boost::int64_t *)calloc(nPts,sizeof(boost::int64_t));
  ib = (boost::int64_t *)calloc(nPts,sizeof(boost::int64_t));
  crit = (real *)calloc(nPts,sizeof(real));

  clusterit((real *)dataContig->data,nPts,sz,option,ia,ib,crit);

  dims[0] = nPts;
  res = PyTuple_New(3);

  //  NOTE: these operations maintain pointers to the respective arrays,
  //  that's why it's ok that we do not free them in this function,
  //  Python will take care of it for us.
  //
  tmp = PyArray_SimpleNewFromData(1,dims,NPY_LONG,(void *)ia);
  PyTuple_SetItem(res,0,(PyObject *)tmp);

  tmp = PyArray_SimpleNewFromData(1,dims,NPY_LONG,(void *)ib);
  PyTuple_SetItem(res,1,(PyObject *)tmp);

  tmp = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,(void *)crit);
  PyTuple_SetItem(res,2,(PyObject *)tmp);

  return res;
};



void distclusterit(real *dists,boost::int64_t n,boost::int64_t iopt,
		   boost::int64_t *ia,boost::int64_t *ib,real *crit){
  boost::int64_t len;

  len = (n*(n-1))/2;
  distdriver_(&n,&len,dists,&iopt,ia,ib,crit);
};


static PyObject *
Clustering_MurtaghDistCluster(python::object data, int nPts, int option)
{
  PyArrayObject *dataContig;
  boost::int64_t *ia,*ib;
  real *crit;
  PyObject *res=PyTuple_New(3);
  PyObject *tmp;
  npy_intp dims[] = {1};

  if (PyArray_Check(data.ptr())) {
    dataContig 
      = reinterpret_cast<PyArrayObject *>(PyArray_ContiguousFromObject(data.ptr(),PyArray_DOUBLE,1,1));
  }
  else {
    throw_value_error("PyArray_Type expected as input");
  }

  ia = (boost::int64_t *)calloc(nPts,sizeof(boost::int64_t));
  ib = (boost::int64_t *)calloc(nPts,sizeof(boost::int64_t));
  crit = (real *)calloc(nPts,sizeof(real));
  distclusterit((real *)dataContig->data,nPts,option,ia,ib,crit);

  dims[0] = nPts;

  //
  //  NOTE: these operations maintain pointers to the respective arrays,
  //  that's why it's ok that we do not free them in this function,
  //  Python will take care of it for us.
  //
  tmp = PyArray_SimpleNewFromData(1,dims,NPY_LONG,(void *)ia);
  PyTuple_SetItem(res,0,tmp);

  tmp = PyArray_SimpleNewFromData(1,dims,NPY_LONG,(void *)ib);
  PyTuple_SetItem(res,1,tmp);

  tmp = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,(void *)crit);
  PyTuple_SetItem(res,2,tmp);
  
  return res;
};


BOOST_PYTHON_MODULE(Clustering) {

  rdkit_import_array();

  python::def("MurtaghCluster", Clustering_MurtaghCluster,
	      ( python::arg("data"), python::arg("nPts"), 
		python::arg("sz"), python::arg("option") ),
	      "TODO: provide docstring");
  python::def("MurtaghDistCluster", Clustering_MurtaghDistCluster,
	      ( python::arg("data"), python::arg("nPts"), 
		python::arg("option") ),
	      "TODO: provide docstring");
}

