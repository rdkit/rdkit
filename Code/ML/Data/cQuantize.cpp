// $Id$
//
// Copyright 2003-2008 Rational Discovery LLC and Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef WIN32
#define CQUANTIZE_EXPORTS
#endif
#include "cQuantize.h"
#include <cstring>
#include <numpy/oldnumeric.h>
#include <ML/InfoTheory/InfoGainFuncs.h>
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

/***********************************************

   constructs a variable table for the data passed in
   The table for a given variable records the number of times each possible value
    of that variable appears for each possible result of the function.

  **Arguments**  

   - vals: pointer to double, contains the values of the variable,
     should be sorted

   - nVals: int, the length of _vals_

   - cuts: pointer to int, the indices of the quantization bounds

   - nCuts: int, the length of _cuts_

   - starts: pointer to int, the potential starting points for quantization bounds

   - nStarts: int, the length of _starts_

   - results: poitner to int, the result codes

   - nPossibleRes: int, the number of possible result codes
 

  **Returns**

    _varTable_ (a pointer to int), which is also modified in place.

  **Notes:**

    - _varTable_ is modified in place

    - the _results_ array is assumed to be _nVals_ long

 ***********************************************/
long int *
GenVarTable(double *vals,int nVals,long int *cuts,int nCuts,long int *starts,
	    long int *results,int nPossibleRes,long int *varTable)
{
  int nBins = nCuts + 1;
  int idx,i,iTab;

  memset(varTable,0,nBins*nPossibleRes*sizeof(long int));
  idx = 0;
  for(i=0;i<nCuts;i++){
    int cut = cuts[i];
    iTab = i*nPossibleRes;
    while(idx<starts[cut]){
      varTable[iTab+results[idx]] += 1;
      idx++;
    }
  }
  iTab = nCuts*nPossibleRes;
  while(idx<nVals){
    varTable[iTab+results[idx]] += 1;
    idx++;
  }
  return varTable;  
}

/***********************************************

 This actually does the recursion required by *cQuantize_RecurseOnBounds()*,
  we do things this way to avoid having to convert things back and forth
  from Python objects

  **Arguments**  

   - vals: pointer to double, contains the values of the variable,
     should be sorted

   - nVals: int, the length of _vals_

   - cuts: pointer to int, the indices of the quantization bounds

   - nCuts: int, the length of _cuts_

   - which: int, the quant bound being modified here

   - starts: pointer to int, the potential starting points for quantization bounds

   - nStarts: int, the length of _starts_

   - results: poitner to int, the result codes

   - nPossibleRes: int, the number of possible result codes
 

  **Returns**

    a double, the expected information gain for the best bounds found
      (which are found in _cuts_ )

  **Notes:**

    - _cuts_ is modified in place

    - the _results_ array is assumed to be _nVals_ long

 ***********************************************/
double
RecurseHelper(double *vals,int nVals,long int *cuts,int nCuts,int which,
	      long int *starts,int nStarts,long int *results,int nPossibleRes)
{
  double maxGain=-1e6,gainHere;
  long int *bestCuts,*tCuts;
  long int *varTable=0;
  int highestCutHere = nStarts - nCuts + which;
  int i,nBounds=nCuts;
  
  varTable = (long int *)calloc((nCuts+1)*nPossibleRes,sizeof(long int));
  bestCuts = (long int *)calloc(nCuts,sizeof(long int));
  tCuts = (long int *)calloc(nCuts,sizeof(long int));
  GenVarTable(vals,nVals,cuts,nCuts,starts,results,nPossibleRes,varTable);
  while(cuts[which] <= highestCutHere){
    gainHere = RDInfoTheory::InfoEntropyGain(varTable,nCuts+1,nPossibleRes);
    if(gainHere > maxGain){
      maxGain = gainHere;
      memcpy(bestCuts,cuts,nCuts*sizeof(long int));
    }

    // recurse on the next vars if needed
    if(which < nBounds-1){
      memcpy(tCuts,cuts,nCuts*sizeof(long int));
      gainHere = RecurseHelper(vals,nVals,tCuts,nCuts,which+1,starts,nStarts,
			       results,nPossibleRes);
      if(gainHere > maxGain){
        maxGain = gainHere;
	memcpy(bestCuts,tCuts,nCuts*sizeof(long int));
      }
    }

    // update this cut
    int oldCut = cuts[which];
    cuts[which] += 1;
    int top,bot;
    bot = starts[oldCut];
    if(oldCut+1 < nStarts)
      top = starts[oldCut+1];
    else
      top = starts[nStarts-1];
    for(i=bot;i<top;i++) {
      int v=results[i];
      varTable[which*nPossibleRes+v] += 1;
      varTable[(which+1)*nPossibleRes+v] -= 1;
    }
    for(i=which+1;i<nBounds;i++){
      if(cuts[i] == cuts[i-1]) cuts[i] += 1;
    }
  }
  memcpy(cuts,bestCuts,nCuts*sizeof(long int));
  free(tCuts);
  free(bestCuts);
  free(varTable);
  return maxGain;
}


/***********************************************
 
   Recursively finds the best quantization boundaries

   **Arguments**

     - vals: a 1D Numeric array with the values of the variables,
       this should be sorted

     - cuts: a list with the indices of the quantization bounds
       (indices are into _starts_ )

     - which: an integer indicating which bound is being adjusted here
       (and index into _cuts_ )

     - starts: a list of potential starting points for quantization bounds

     - results: a 1D Numeric array of integer result codes

     - nPossibleRes: an integer with the number of possible result codes

   **Returns**

     - a 2-tuple containing:

       1) the best information gain found so far

       2) a list of the quantization bound indices ( _cuts_ for the best case)
   
   **Notes**

    - this is not even remotely efficient, which is why a C replacement
      was written

    - this is a drop-in replacement for *ML.Data.Quantize._PyRecurseBounds*
						
 ***********************************************/
static PyObject *
cQuantize_RecurseOnBounds(PyObject *self, PyObject *args)
{
  PyObject *vals,*results;
  PyObject *pyCuts,*pyStarts;
  int which,nPossibleRes;

  PyArrayObject *contigVals,*contigResults;
  long int *cuts,*starts;
  PyObject *res,*cutObj;
  double gain;
  int i,nCuts,nStarts;

  if (!PyArg_ParseTuple(args, "OOiOOi",&vals,
			&pyCuts,&which,&pyStarts,
			&results,&nPossibleRes))
    return NULL;

  /*
    -------

    Setup code

    -------
  */
  contigVals = (PyArrayObject *)PyArray_ContiguousFromObject(vals,PyArray_DOUBLE,1,1);
  if(!contigVals){
    PyErr_SetString(PyExc_ValueError,
                    "could not convert value argument");
    return 0;
  }
  contigResults = (PyArrayObject *)PyArray_ContiguousFromObject(results,PyArray_LONG,1,1);
  if(!contigResults){
    PyErr_SetString(PyExc_ValueError,
                    "could not convert results argument");
    return 0;
  }
  
  nCuts = PyList_Size(pyCuts);
  cuts = (long int *)calloc(nCuts,sizeof(long int));
  for(i=0;i<nCuts;i++){
    cuts[i] = PyInt_AsLong(PyList_GetItem(pyCuts,i));
  }
  nStarts = PyList_Size(pyStarts);
  starts = (long int *)calloc(nStarts,sizeof(long int));
  for(i=0;i<nStarts;i++){
    starts[i] = PyInt_AsLong(PyList_GetItem(pyStarts,i));
  }

  // do the real work
  gain = RecurseHelper((double *)contigVals->data,contigVals->dimensions[0],
		       cuts,nCuts,which,starts,nStarts,
		       (long int *)contigResults->data,nPossibleRes);
		       
  /*
    -------

    Construct the return value

    -------
  */
  res = PyTuple_New(2);
  PyTuple_SetItem(res,0,PyFloat_FromDouble(gain));
  
  cutObj = PyList_New(nCuts);
  for(i=0;i<nCuts;i++){
    PyList_SetItem(cutObj,i,PyInt_FromLong((long int)cuts[i]));
  }
  PyTuple_SetItem(res,1,cutObj);
  free(cuts);
  free(starts);
  Py_DECREF(contigVals);
  Py_DECREF(contigResults);
  return res;
    
}

static PyObject *
cQuantize_FindStartPoints(PyObject *self, PyObject *args)
{
  PyObject *values,*results;
  PyArrayObject *contigVals,*contigResults;
  int nData;

  if (!PyArg_ParseTuple(args, "OOi",&values,&results,&nData)){
    return NULL;
  }
  PyObject *startPts = PyList_New(0);

  if(nData<2){
    return startPts;
  }

  contigVals = (PyArrayObject *)PyArray_ContiguousFromObject(values,PyArray_DOUBLE,1,1);
  if(!contigVals){
    PyErr_SetString(PyExc_ValueError,
                    "could not convert value argument");
    return 0;
  }
  double *vals=(double *)contigVals->data;

  contigResults = (PyArrayObject *)PyArray_ContiguousFromObject(results,PyArray_LONG,1,1);
  if(!contigResults){
    PyErr_SetString(PyExc_ValueError,
                    "could not convert results argument");
    return 0;
  }
  long *res=(long *)contigResults->data;

  bool firstBlock=true;
  long lastBlockAct=-2,blockAct=res[0];
  int lastDiv=-1;
  double tol=1e-8;

  int i=1;
  while(i<nData){
    while(i<nData && vals[i]-vals[i-1]<=tol){
      if(res[i]!=blockAct){
        blockAct=-1;
      }
      ++i;
    }
    if(firstBlock){
      firstBlock=false;
      lastBlockAct=blockAct;
      lastDiv=i;
    } else {
      if(blockAct==-1 || lastBlockAct==-1 || blockAct!=lastBlockAct){
        PyObject *pyint=PyInt_FromLong(lastDiv);
        PyList_Append(startPts,pyint);
        Py_DECREF(pyint);
        lastDiv=i;
        lastBlockAct=blockAct;
      } else {
        lastDiv=i;
      }
    }
    if(i<nData) blockAct=res[i];
    ++i; 
  }

  // catch the case that the last point also sets a bin:
  if( blockAct != lastBlockAct ){
    PyObject *pyint=PyInt_FromLong(lastDiv);
    PyList_Append(startPts,pyint);
    Py_DECREF(pyint);
  }
  Py_DECREF(contigVals);
  Py_DECREF(contigResults);
  return startPts;
}
static PyMethodDef cQuantizeMethods[] = {
  {"_RecurseOnBounds",cQuantize_RecurseOnBounds,METH_VARARGS},
  {"_FindStartPoints",cQuantize_FindStartPoints,METH_VARARGS},
  {NULL,NULL}
};

CQUANTIZE_API void initcQuantize()
{
  (void) Py_InitModule("cQuantize",cQuantizeMethods);
  import_array();
}



