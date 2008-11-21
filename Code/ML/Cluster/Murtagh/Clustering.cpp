// $Id$
//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//     All Rights Reserved
//
#ifdef WIN32
#define CALGORITHMS_EXPORTS
#endif
#define PYTH_FILE_WITH_INIT
#include "Clustering.h"
#include <numpy/arrayobject.h>

#ifdef WIN32
BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
    }
    return TRUE;
}
#endif

typedef double real;



extern "C"
void distdriver_(long int *n,long int *len,
		 real *dists,
		 long int *toggle,
		 long int *ia,long int *ib,real *crit);

//
// Rather than deal with any nonsense like trying to get
// the distance matrix built properly on the f2c side of things
// (thus drowning in the waves of f2c hate), we'll generate
// the distance matrix on our own here and then call distdriver_
//
void clusterit(real *dataP,long int n,long int m,long int iopt,
	       long int *ia,long int *ib,real *crit){
  real *dists;
  long int len;
  long int pos = 0;
  long int i,j,k,iTab,jTab;
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
Clustering_MurtaghCluster(PyObject *self,PyObject *args)
{
  PyArrayObject *dataContig;
  PyObject *data;
  int nPts,sz,option;
  long int *ia,*ib;
  real *crit;
  PyObject *res;
  PyObject *tmp;
  npy_intp dims[2];

  if(!PyArg_ParseTuple(args,"Oiii",&data,&nPts,&sz,&option))
    return NULL;
  dataContig = (PyArrayObject *)PyArray_ContiguousFromObject(data,PyArray_DOUBLE,2,2);

  ia = (long int *)calloc(nPts,sizeof(long int));
  ib = (long int *)calloc(nPts,sizeof(long int));
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

  Py_DECREF(dataContig);
  return res;
};



void distclusterit(real *dists,long int n,long int iopt,
		   long int *ia,long int *ib,real *crit){
  long int len;

  len = (n*(n-1))/2;
  distdriver_(&n,&len,dists,&iopt,ia,ib,crit);
};


static PyObject *
Clustering_MurtaghDistCluster(PyObject *self,PyObject *args)
{
  PyArrayObject *dataContig;
  PyObject *data;
  long int nPts,option;
  long int *ia,*ib;
  real *crit;
  PyObject *res=PyTuple_New(3);
  PyObject *tmp;
  npy_intp dims[] = {1};

  if(!PyArg_ParseTuple(args,"Oii",&data,&nPts,&option))
    return NULL;

  dataContig = (PyArrayObject *)PyArray_ContiguousFromObject(data,PyArray_DOUBLE,1,1);
  ia = (long int *)calloc(nPts,sizeof(long int));
  ib = (long int *)calloc(nPts,sizeof(long int));
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
  
  Py_DECREF(dataContig);
  return res;
};

static PyMethodDef locMethods[] = {
  {"MurtaghCluster",Clustering_MurtaghCluster,METH_VARARGS},
  {"MurtaghDistCluster",Clustering_MurtaghDistCluster,METH_VARARGS},
  {NULL,NULL}
};

CALGORITHMS_API void initClustering()
{
  (void) Py_InitModule("Clustering",locMethods);
#ifndef NO_IMPORT_ARRAY
  import_array();
#endif
}

