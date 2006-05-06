//
//  Copyright (C) 2003  Rational Discovery LLC
//   All Rights Reserved
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#define PY_ARRAY_UNIQUE_SYMBOL rdinfotheory_array_API
#include "Numeric/arrayobject.h"

#include <RDBoost/Wrap.h>

#include <ML/InfoTheory/InfoBitRanker.h>
#include <DataStructs/BitVects.h>
#include <RDBoost/PySequenceHolder.h>

namespace python = boost::python;

namespace RDInfoTheory {
  
  PyObject *getTopNbits(InfoBitRanker *ranker, int num){// int ignoreNoClass=-1) {
    double *dres = ranker->getTopN(num);
    int dims[2];
    dims[0] = num;
    dims[1] = ranker->getNumClasses() + 2;
    PyArrayObject *res = (PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE);
    memcpy(static_cast<void *>(res->data),
           static_cast<void *>(dres), dims[0]*dims[1]*sizeof(double));
    return (PyObject *) res;
  }
  //BOOST_PYTHON_FUNCTION_OVERLOADS(getTopNbits_ol, getTopNbits, 2, 3)

  void AccumulateVotes(InfoBitRanker *ranker, python::object bitVect, int label) {
    python::extract<ExplicitBitVect> ebvWorks(bitVect);
    python::extract<SparseBitVect> sbvWorks(bitVect);
    if (ebvWorks.check()) {
      ExplicitBitVect ev = python::extract<ExplicitBitVect>(bitVect);
      ranker->accumulateVotes(ev, label);
    }
    else if (sbvWorks.check()) {
      SparseBitVect sv = python::extract<SparseBitVect>(bitVect);
      ranker->accumulateVotes(sv, label);
    }
    else {
      throw_value_error("Accumulate Vote can only take a explicitBitVects or SparseBitvects");
    }
  }
  
  void SetBiasList(InfoBitRanker *ranker, python::object classList) {
    RDKit::INT_VECT cList;
    PySequenceHolder<int> bList(classList);
    cList.reserve(bList.size());
    for (unsigned int i = 0; i < bList.size(); i++) {
      cList.push_back(bList[i]);
    }
    ranker->setBiasList(cList);
  }

  void SetMaskBits(InfoBitRanker *ranker, python::object maskBits) {
    RDKit::INT_VECT cList;
    PySequenceHolder<int> bList(maskBits);
    cList.reserve(bList.size());
    for (unsigned int i = 0; i < bList.size(); i++) {
      cList.push_back(bList[i]);
    }
    ranker->setMaskBits(cList);
  }

  void tester(InfoBitRanker *ranker, python::object bitVect) {
    python::extract<SparseBitVect> sbvWorks(bitVect);
    if (sbvWorks.check()){
      SparseBitVect sv = python::extract<SparseBitVect>(bitVect);
      std::cout << "Num of on bits: " << sv.GetNumOnBits() << "\n";
    }
  }

  struct ranker_wrap {
    static void wrap() {
      std::string docString = "A class to rank the bits from a series of labelled fingerprints\n";
      python::class_<InfoBitRanker>("InfoBitRanker",
                                    docString.c_str(),
                                    python::init<int, int>(python::args("nBits", "nClasses")))
        .def(python::init<int, int, InfoBitRanker::InfoType>
                                    (python::args("nBits", "nClasses", "infoType")))
        .def("AccumulateVotes", AccumulateVotes,
             "Accumulate the votes for all the bits turned on in a bit vector\n\n"
             "ARGUMENTS:\n\n"
             "  - bv : bit vector either ExplicitBitVect or SparseBitVect operator\n"
             "  - label : the class label for the bit vector. It is assumed that 0 <= class < nClasses \n")
        .def ("SetBiasList", SetBiasList,
              "Set the classes to which the entropy calculation should be biased\n\n"
              "This list contains a set of class ids used when in the BIASENTROPY mode of ranking bits. \n"
              "In this mode, a bit must be correllated higher with one of the biased classes than all the \n"
              "other classes. For example, in a two class problem with actives and inactives, the fraction of \n"
              "actives that hit the bit has to be greater than the fraction of inactives that hit the bit\n\n"
              "ARGUMENTS: \n\n"
              "  - classList : list of class ids that we want a bias towards\n")
        .def ("SetMaskBits", SetMaskBits,
              "Set the mask bits for the calculation\n\n"
              "ARGUMENTS: \n\n"
              "  - maskBits : list of mask bits to use\n")
        .def("GetTopN", getTopNbits,
             "Returns the top n bits ranked by the information metric\n"
             "This is actually the function where most of the work of ranking is happening\n\n"
             "ARGUMENTS:\n\n"
             "  - num : the number of top ranked bits that are required\n")
        .def("WriteTopBitsToFile", &InfoBitRanker::writeTopBitsToFile,
             "Write the bits that have been ranked to a file")
        .def("Tester", tester)
        ;
      
      python::enum_<InfoBitRanker::InfoType>("InfoType")
        .value("ENTROPY", InfoBitRanker::ENTROPY)
        .value("BIASENTROPY", InfoBitRanker::BIASENTROPY)
        .value("CHISQUARE", InfoBitRanker::CHISQUARE)
        .value("BIASCHISQUARE", InfoBitRanker::BIASCHISQUARE)
        ;
    };
  };
}

void wrap_ranker() {
  RDInfoTheory::ranker_wrap::wrap();
}

