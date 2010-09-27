// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and  Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//


#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#define PY_ARRAY_UNIQUE_SYMBOL rdinfotheory_array_API
#include "numpy/arrayobject.h"

#include <RDBoost/Wrap.h>
#include <RDBoost/PySequenceHolder.h>
#include <ML/InfoTheory/CorrMatGenerator.h>
#include <RDGeneral/types.h>

namespace python = boost::python;

namespace RDInfoTheory {
  
  PyObject *getCorrMatrix(BitCorrMatGenerator *cmGen) {
    double *dres = cmGen->getCorrMat();
    unsigned int nb = cmGen->getCorrBitList().size();
    npy_intp dim = nb*(nb-1)/2;
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(1,&dim,NPY_DOUBLE);
    memcpy(static_cast<void *>(res->data),
           static_cast<void *>(dres), dim*sizeof(double));
    return PyArray_Return(res);
  }

  void setBitList(BitCorrMatGenerator *cmGen, python::object bitList) {
    PySequenceHolder<int> blist(bitList);
    unsigned int nb = blist.size();
    RDKit::INT_VECT res;
    res.reserve(nb);
    for (unsigned int i = 0; i < nb; i++) {
      res.push_back(blist[i]);
    }
    cmGen->setBitIdList(res);
  }

  void CollectVotes(BitCorrMatGenerator *cmGen, python::object bitVect) {
    python::extract<ExplicitBitVect> ebvWorks(bitVect);
    python::extract<SparseBitVect> sbvWorks(bitVect);
    if (ebvWorks.check()) {
      ExplicitBitVect ev = python::extract<ExplicitBitVect>(bitVect);
      cmGen->collectVotes(ev);
    }
    else if (sbvWorks.check()) {
      SparseBitVect sv = python::extract<SparseBitVect>(bitVect);
      cmGen->collectVotes(sv);
    }
    else {
      throw_value_error("CollectVote can only take ExplicitBitVects or SparseBitVects");
    }
  }

  struct corrmat_wrap {
    static void wrap() {
      std::string docString = "A class to generate a pariwise correlation matrix between a list of bits\n"
        "The mode of operation for this class is something like this\n"
        "   >>> cmg = BitCorrMatGenerator() \n"
        "   >>> cmg.SetBitList(blist) \n"
        "   >>> for fp in fpList:  \n"
        "   >>>    cmg.CollectVotes(fp)  \n"
        "   >>> corrMat = cmg.GetCorrMatrix() \n"
        "    \n"
        "   The resulting correlation matrix is a one dimensional nummeric array containing the \n"
        "   lower triangle elements\n";
      python::class_<BitCorrMatGenerator>("BitCorrMatGenerator",
                                          docString.c_str())
        .def("SetBitList", setBitList,
             "Set the list of bits that need to be correllated\n\n"
             " This may for example be ther top ranking ensemble bits\n\n"
             "ARGUMENTS:\n\n"
             "  - bitList : an integer list of bit IDs\n")
        .def("CollectVotes", CollectVotes,
             "For each pair of on bits (bi, bj) in fp increase the correlation count for the pair by 1\n\n"
             "ARGUMENTS:\n\n"
             "  - fp : a bit vector to collect the fingerprints from\n")
        .def("GetCorrMatrix", getCorrMatrix,
             "Get the correlation matrix following the collection of votes from a bunch of fingerprints\n")
        ;
    };
  };
}

void wrap_corrmatgen() {
  RDInfoTheory::corrmat_wrap::wrap();
}


