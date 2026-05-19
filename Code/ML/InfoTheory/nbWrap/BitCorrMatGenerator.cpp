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

#include <nanobind/nanobind.h>
#include <ML/InfoTheory/CorrMatGenerator.h>
#include <RDGeneral/types.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDInfoTheory {

nb::object getCorrMatrix(BitCorrMatGenerator *cmGen) {
  double *dres = cmGen->getCorrMat();
  unsigned int nb_size = cmGen->getCorrBitList().size();
  npy_intp dim = nb_size * (nb_size - 1) / 2;
  auto *res = (PyArrayObject *)PyArray_SimpleNew(1, &dim, NPY_DOUBLE);
  memcpy(static_cast<void *>(PyArray_DATA(res)), static_cast<void *>(dres),
         dim * sizeof(double));
  return nb::steal<nb::object>(PyArray_Return(res));
}

void setBitList(BitCorrMatGenerator *cmGen, nb::iterable bitList) {
  RDKit::INT_VECT res;
  for (auto item : bitList) {
    res.push_back(nb::cast<int>(item));
  }
  cmGen->setBitIdList(res);
}

void CollectVotes(BitCorrMatGenerator *cmGen, nb::object bitVect) {
  if (nb::isinstance<ExplicitBitVect>(bitVect)) {
    cmGen->collectVotes(nb::cast<ExplicitBitVect &>(bitVect));
  } else if (nb::isinstance<SparseBitVect>(bitVect)) {
    cmGen->collectVotes(nb::cast<SparseBitVect &>(bitVect));
  } else {
    throw nb::value_error(
        "CollectVote can only take ExplicitBitVects or SparseBitVects");
  }
}

}  // namespace RDInfoTheory

void wrap_corrmatgen(nb::module_ &m) {
  nb::class_<RDInfoTheory::BitCorrMatGenerator>(m, "BitCorrMatGenerator",
      R"DOC(A class to generate a pairwise correlation matrix between a list of bits
The mode of operation for this class is something like this

   >>> cmg = BitCorrMatGenerator()
   >>> cmg.SetBitList(blist)
   >>> for fp in fpList:
   >>>    cmg.CollectVotes(fp)
   >>> corrMat = cmg.GetCorrMatrix()

   The resulting correlation matrix is a one dimensional nummeric array containing the
   lower triangle elements
)DOC")
      .def(nb::init<>())
      .def("SetBitList", RDInfoTheory::setBitList, "bitList"_a,
           R"DOC(Set the list of bits that need to be correllated

This may for example be their top ranking ensemble bits

ARGUMENTS:

  - bitList : an integer list of bit IDs
)DOC")
      .def("CollectVotes", RDInfoTheory::CollectVotes, "bitVect"_a,
           R"DOC(For each pair of on bits (bi, bj) in fp increase the correlation count for the pair by 1

ARGUMENTS:

  - fp : a bit vector to collect the fingerprints from
)DOC")
      .def("GetCorrMatrix", RDInfoTheory::getCorrMatrix,
           R"DOC(Get the correlation matrix following the collection of votes from a bunch of fingerprints
)DOC");
}
