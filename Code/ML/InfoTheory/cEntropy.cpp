// $Id$
//
//  Copyright (C) 2001-2008 greg Landrum and Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef WIN32
#define CENTROPY_EXPORTS
#endif
#include "cEntropy.h"
#include <numpy/oldnumeric.h>
#include <math.h>

#if defined(WIN32) && defined(CENTROPY_EXPORTS)
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

/************
   calculates the informational entropy of the values in an array

 **Arguments**
    
    - tPtr: pointer to a long int array containing the data

    - dim: long int containing the length of the _tPtr_ array.

  **Returns**

    a double


************/
template<class T>
double
InfoEntropy(T *tPtr,long int dim)
{
  int i;
  T nInstances = 0;
  double accum=0.0,d;

  for(i=0;i<dim;i++){
    nInstances += tPtr[i];
  }

  if(nInstances != 0){
    for(i=0;i<dim;i++){
      d = (double)tPtr[i]/nInstances;
      if(d != 0){
	accum += -d*log(d);
      }
    }
  }
  return accum/log(2.0);
}


/************
   calculates the informational entropy of the values in an Numeric array

   **Arguments**

     - resultsArray: a Numeric Array object
  
      For example, if a function has 3 possible results, and the
       variable in question hits them 5, 6 and 1 times each,
       resultsArray would be [5,6,1]

   **Returns**

     - a Python float object

   **Notes**

     - this is a dropin replacement for _PyInfoEntropy()_ in entropy.py

************/
static PyObject *
cEntropy_InfoEntropy(PyObject *self, PyObject *args)
{
  PyArrayObject *resultsContig;
  PyObject *resultsArray;
  double res;
  

  if (!PyArg_ParseTuple(args, "O!",&PyArray_Type, &resultsArray))
    return NULL;
  if(((PyArrayObject *)resultsArray)->descr->type_num == PyArray_DOUBLE ||
     ((PyArrayObject *)resultsArray)->descr->type_num == PyArray_FLOAT){
    resultsContig = (PyArrayObject *)PyArray_ContiguousFromObject(resultsArray,PyArray_DOUBLE,1,1);
    res = InfoEntropy((double *)(resultsContig->data),
		      (long int)(resultsContig->dimensions[0]));
  } else {
    resultsContig = (PyArrayObject *)PyArray_ContiguousFromObject(resultsArray,PyArray_LONG,1,1);
    res = InfoEntropy((long int *)(resultsContig->data),
		      (long int)(resultsContig->dimensions[0]));
  }

  Py_DECREF(resultsContig);
  return Py_BuildValue("d",res);
}

CENTROPY_API double
InfoGain(long int *dMat,long int dim1,long int dim2)
{
  int i,j;
  long int *variableRes, *overallRes;
  double gain,term2;
  int tSum;

  variableRes = (long int *)calloc(dim1,sizeof(long int));
  // do the row sums
  for(i=0;i<dim1;i++){
    int idx1 = i*dim2;
    for(j=0;j<dim2;j++){
      variableRes[i] += dMat[idx1+j];
    }
  }

  overallRes = (long int *)calloc(dim2,sizeof(long int));
  // do the col sums
  for(i=0;i<dim2;i++){
    for(j=0;j<dim1;j++){
      overallRes[i] += dMat[j*dim2+i];
    }
  }

  term2 = 0.0;
  for(i=0;i<dim1;i++){
    long int *tPtr;
    tPtr = dMat + i*dim2;
    term2 += variableRes[i] * InfoEntropy(tPtr,dim2);
  }
  tSum = 0;
  for(i=0;i<dim2;i++){
    tSum += overallRes[i];
  }

  if(tSum != 0){
    term2 /= tSum;
    gain = InfoEntropy(overallRes,dim2) - term2;
  }
  else{
    gain = 0.0;
  }

  free(overallRes);
  free(variableRes);

  return gain;
}

/************
   calculates the information gain for a variable

   **Arguments**

     - varMat: a Numeric Array object

       varMat is a Numeric array with the number of possible occurances
         of each result for reach possible value of the given variable.

       So, for a variable which adopts 4 possible values and a result which
         has 3 possible values, varMat would be 4x3
  
   **Returns**

     - a Python float object

   **Notes**

     - this is a dropin replacement for _PyInfoGain()_ in entropy.py

************/
static PyObject *
cEntropy_InfoGain(PyObject *self, PyObject *args)
{
  PyArrayObject *varMatContig;
  PyObject *varMat;
  long int dim1,dim2;
  long int *dMat;
  double gain;
  

  // FIX: this crashes if we pass in anything other than an int array
  if (!PyArg_ParseTuple(args, "O!",&PyArray_Type, &varMat))
    return NULL;
  varMatContig = (PyArrayObject *)PyArray_ContiguousFromObject(varMat,
  							       PyArray_LONG,2,2);
  dMat = (long int *)varMatContig->data;
  dim1 = varMatContig->dimensions[0];
  dim2 = varMatContig->dimensions[1];

  gain = InfoGain(dMat,dim1,dim2);

  Py_DECREF(varMatContig);
  return Py_BuildValue("d",gain);
}




//  ------------- Initialization foo  --------------------

static PyMethodDef cEntropyMethods[] = {
  {"InfoEntropy",cEntropy_InfoEntropy,METH_VARARGS},
  {"InfoGain",cEntropy_InfoGain,METH_VARARGS},
  {NULL,NULL}
};

CENTROPY_API void initcEntropy()
{
  (void) Py_InitModule("cEntropy",cEntropyMethods);
  import_array();
}



