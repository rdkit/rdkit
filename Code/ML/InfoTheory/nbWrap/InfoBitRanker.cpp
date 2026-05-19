//
//  Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL rdinfotheory_array_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <iostream>
#include <nanobind/nanobind.h>
#include <ML/InfoTheory/InfoBitRanker.h>
#include <DataStructs/BitVects.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDInfoTheory {

nb::object getTopNbits(InfoBitRanker *ranker, int num) {
  double *dres = ranker->getTopN(num);
  npy_intp dims[2];
  dims[0] = num;
  dims[1] = ranker->getNumClasses() + 2;
  auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  memcpy(static_cast<void *>(PyArray_DATA(res)), static_cast<void *>(dres),
         dims[0] * dims[1] * sizeof(double));
  return nb::steal<nb::object>(PyArray_Return(res));
}

void AccumulateVotes(InfoBitRanker *ranker, nb::object bitVect, int label) {
  if (nb::isinstance<ExplicitBitVect>(bitVect)) {
    ranker->accumulateVotes(nb::cast<ExplicitBitVect &>(bitVect), label);
  } else if (nb::isinstance<SparseBitVect>(bitVect)) {
    ranker->accumulateVotes(nb::cast<SparseBitVect &>(bitVect), label);
  } else {
    throw nb::value_error(
        "Accumulate Vote can only take a explicitBitVects or SparseBitvects");
  }
}

void SetBiasList(InfoBitRanker *ranker, nb::iterable classList) {
  RDKit::INT_VECT cList;
  for (auto item : classList) {
    cList.push_back(nb::cast<int>(item));
  }
  ranker->setBiasList(cList);
}

void SetMaskBits(InfoBitRanker *ranker, nb::iterable maskBits) {
  RDKit::INT_VECT cList;
  for (auto item : maskBits) {
    cList.push_back(nb::cast<int>(item));
  }
  ranker->setMaskBits(cList);
}

void tester(InfoBitRanker *, nb::object bitVect) {
  if (nb::isinstance<SparseBitVect>(bitVect)) {
    std::cout << "Num of on bits: "
              << nb::cast<SparseBitVect &>(bitVect).getNumOnBits() << "\n";
  }
}

}  // namespace RDInfoTheory

void wrap_ranker(nb::module_ &m) {
  nb::class_<RDInfoTheory::InfoBitRanker>(m, "InfoBitRanker",
      R"DOC(A class to rank the bits from a series of labelled fingerprints
A simple demonstration may help clarify what this class does.
Here's a small set of vectors:

>>> for i,bv in enumerate(bvs): print(bv.ToBitString(),acts[i])
...
0001 0
0101 0
0010 1
1110 1

Default ranker, using infogain:

>>> ranker = InfoBitRanker(4,2)
>>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
...
>>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
...
3 1.000 2 0
2 1.000 0 2
0 0.311 0 1

Using the biased infogain:

>>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASENTROPY)
>>> ranker.SetBiasList((1,))
>>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
...
>>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
...
2 1.000 0 2
0 0.311 0 1
1 0.000 1 1

A chi squared ranker is also available:

>>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.CHISQUARE)
>>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
...
>>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
...
3 4.000 2 0
2 4.000 0 2
0 1.333 0 1

As is a biased chi squared:

>>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASCHISQUARE)
>>> ranker.SetBiasList((1,))
>>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
...
>>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
...
2 4.000 0 2
0 1.333 0 1
1 0.000 1 1
)DOC")
      .def(nb::init<unsigned int, unsigned int>(), "nBits"_a, "nClasses"_a)
      .def(nb::init<unsigned int, unsigned int, RDInfoTheory::InfoBitRanker::InfoType>(),
           "nBits"_a, "nClasses"_a, "infoType"_a)
      .def("AccumulateVotes", RDInfoTheory::AccumulateVotes, "bitVect"_a,
           "label"_a,
           R"DOC(Accumulate the votes for all the bits turned on in a bit vector

ARGUMENTS:

  - bv : bit vector either ExplicitBitVect or SparseBitVect operator
  - label : the class label for the bit vector. It is assumed that 0 <= class < nClasses
)DOC")
      .def("SetBiasList", RDInfoTheory::SetBiasList, "classList"_a,
           R"DOC(Set the classes to which the entropy calculation should be biased

This list contains a set of class ids used when in the BIASENTROPY mode of ranking bits.
In this mode, a bit must be correlated higher with one of the biased classes than all the
other classes. For example, in a two class problem with actives and inactives, the fraction of
actives that hit the bit has to be greater than the fraction of inactives that hit the bit

ARGUMENTS:

  - classList : list of class ids that we want a bias towards
)DOC")
      .def("SetMaskBits", RDInfoTheory::SetMaskBits, "maskBits"_a,
           R"DOC(Set the mask bits for the calculation

ARGUMENTS:

  - maskBits : list of mask bits to use
)DOC")
      .def("GetTopN", RDInfoTheory::getTopNbits, "num"_a,
           R"DOC(Returns the top n bits ranked by the information metric
This is actually the function where most of the work of ranking is happening

ARGUMENTS:

  - num : the number of top ranked bits that are required
)DOC")
      .def("WriteTopBitsToFile", &RDInfoTheory::InfoBitRanker::writeTopBitsToFile,
           "fileName"_a,
           "Write the bits that have been ranked to a file")
      .def("Tester", RDInfoTheory::tester, "bitVect"_a);

  nb::enum_<RDInfoTheory::InfoBitRanker::InfoType>(m, "InfoType")
      .value("ENTROPY", RDInfoTheory::InfoBitRanker::ENTROPY)
      .value("BIASENTROPY", RDInfoTheory::InfoBitRanker::BIASENTROPY)
      .value("CHISQUARE", RDInfoTheory::InfoBitRanker::CHISQUARE)
      .value("BIASCHISQUARE", RDInfoTheory::InfoBitRanker::BIASCHISQUARE)
      .export_values();
}
