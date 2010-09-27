// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>
#include <RDBoost/Wrap.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>

namespace python = boost::python;

#include "fastentropy.h"

template <typename BV> 
PyObject *runCMIM(python::list examples,unsigned int nToTake){
  //std::cerr << "select " << std::endl;std::cerr.flush();
  unsigned int nExamples=python::extract<int>(examples.attr("__len__")());
  //std::cerr << "NEX " << nExamples << std::endl;std::cerr.flush();
  python::list example=python::extract<python::list>(examples[0]);
  const BV &tmpBV=python::extract<BV>(example[1]);
  unsigned int exSize=python::extract<int>(example.attr("__len__")());
  //std::cerr << "EXSIZE " << exSize << std::endl;std::cerr.flush();
  unsigned int nBits=tmpBV.getNumBits();
  //std::cerr << "nBits " << nBits << std::endl;std::cerr.flush();

  unsigned int sz=(nExamples+31)/32;
  uint32_t *y = new uint32_t[sz];
  memset(y,0,sizeof(uint32_t)*sz);
  uint32_t *raw = new uint32_t[sz*nBits];
  memset(raw,0,sizeof(uint32_t)*sz*nBits);
  uint32_t **x = new uint32_t *[nBits];
  for(unsigned int i=0;i<nBits;i++){
    x[i] = &raw[sz*i];
  }
  
  for(unsigned int i=0;i<nExamples;i++){
    example=python::extract<python::list>(examples[i]);
    fe_set_bit(i,y,python::extract<uint32_t>(example[exSize-1])>0);

    const BV &bv=python::extract<BV>(example[1]);
    //std::cerr << BitVectToText(bv) << std::endl;
    for(unsigned int j=0;j<nBits;j++){
      fe_set_bit(i,x[j],bv.getBit(j));
    }
  }

  int *sels=new int[nToTake];
  fe_selection_cmim(nExamples,nBits,x,y,nToTake,sels);
  delete [] y;y=0;
  delete [] x;x=0;
  delete [] raw;raw=0;

  PyObject *res = PyTuple_New(nToTake);
  for(unsigned int i=0;i<nToTake;i++){
    PyTuple_SetItem(res,i,PyInt_FromLong(sels[i]));
  }
  delete [] sels;sels=0;
  return res;

}

PyObject * selectCMIM(python::list &examples,unsigned int nToTake){
  PyObject *res=0;

  python::list example=python::extract<python::list>(examples[0]);
  
  python::extract<ExplicitBitVect> conv(example[1]);
  if(conv.check()){
    res = runCMIM<ExplicitBitVect>(examples,nToTake);
  } else {
    res = runCMIM<SparseBitVect>(examples,nToTake);
  }
  return res;
}



BOOST_PYTHON_MODULE(rdFeatSelect)
{
  fe_init_tables();
  python::scope().attr("__doc__") =
    "Module containing functions for feature selection"
    ;
  
  std::string docString="";

  python::def("selectCMIM", selectCMIM,
              docString.c_str());
}

