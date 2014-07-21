// $Id$
//
// Copyright 2003-2008 Rational Discovery LLC and Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>

#include <RDBoost/Wrap.h>
#include <numpy/oldnumeric.h>
#include <RDBoost/import_array.h>

namespace python = boost::python;

#include <ML/InfoTheory/InfoGainFuncs.h>

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
static python::tuple
cQuantize_RecurseOnBounds(python::object vals, python::list pyCuts, int which, 
			  python::list pyStarts, python::object results, int nPossibleRes)
{
  PyArrayObject *contigVals,*contigResults;
  long int *cuts,*starts;

  /*
    -------

    Setup code

    -------
  */
  contigVals 
    = reinterpret_cast<PyArrayObject *>(PyArray_ContiguousFromObject(vals.ptr(),PyArray_DOUBLE,1,1));
  if(!contigVals){
    throw_value_error("could not convert value argument");
  }

  contigResults 
    = reinterpret_cast<PyArrayObject *>(PyArray_ContiguousFromObject(results.ptr(),PyArray_LONG,1,1));
  if(!contigResults){
    throw_value_error("could not convert results argument");
  }

  python::ssize_t nCuts = python::len(pyCuts);
  cuts = (long int *)calloc(nCuts,sizeof(long int));
  for (python::ssize_t i=0; i<nCuts; i++) {
    python::object elem = pyCuts[i];
    cuts[i] = python::extract<long int>(elem);
  }

  python::ssize_t nStarts = python::len(pyStarts);
  starts = (long int *)calloc(nStarts,sizeof(long int));
  for (python::ssize_t i=0; i<nStarts; i++){
    python::object elem = pyStarts[i];
    starts[i] = python::extract<long int>(elem);
  }

  // do the real work
  double gain 
    = RecurseHelper((double *)contigVals->data,contigVals->dimensions[0],
		    cuts,nCuts,which,starts,nStarts,
		    (long int *)contigResults->data,nPossibleRes);
		       
  /*
    -------

    Construct the return value

    -------
  */
  python::list cutObj;
  for (python::ssize_t i=0; i<nCuts; i++) {
    cutObj.append(cuts[i]);
  }
  free(cuts);
  free(starts);
  return python::make_tuple(gain, cutObj); 
}

static python::list
cQuantize_FindStartPoints(python::object values, python::object results, 
			  int nData)
{
  python::list startPts;

  if(nData<2){
    return startPts;
  }

  PyArrayObject *contigVals 
    = reinterpret_cast<PyArrayObject *>(PyArray_ContiguousFromObject(values.ptr(),PyArray_DOUBLE,1,1));
  if(!contigVals){
    throw_value_error("could not convert value argument");
  }

  double *vals=(double *)contigVals->data;

  PyArrayObject *contigResults 
    = reinterpret_cast<PyArrayObject *>(PyArray_ContiguousFromObject(results.ptr(),PyArray_LONG,1,1));
  if(!contigResults){
    throw_value_error("could not convert results argument");
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
	startPts.append(lastDiv);
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
    startPts.append(lastDiv);
  }

  return startPts;
}

BOOST_PYTHON_MODULE(cQuantize) {

  rdkit_import_array();

  python::def("_RecurseOnBounds", cQuantize_RecurseOnBounds,
	      ( python::arg("vals"), python::arg("pyCuts"), 
		python::arg("which"), python::arg("pyStarts"), 
		python::arg("results"), python::arg("nPossibleRes") ),
	      "TODO: provide docstring");
  python::def("_FindStartPoints", cQuantize_FindStartPoints,
	      ( python::arg("values"), python::arg("results"), 
		python::arg("nData") ),
	      "TODO: provide docstring");
}



