// $Id$
//
//  Copyright (C) 2007-2008 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <RDGeneral/BoostStartInclude.h>
#include <cstdint>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <DataStructs/SparseIntVect.h>

using namespace RDKit;

namespace {
template <typename IndexType>
python::object SIVToBinaryText(const SparseIntVect<IndexType> &siv) {
  std::string res = siv.toString();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}
}  // namespace

template <typename IndexType>
struct siv_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const SparseIntVect<IndexType> &self) {
    return python::make_tuple(SIVToBinaryText(self));
  };
};

namespace {
template <typename IndexType>
void pyUpdateFromSequence(SparseIntVect<IndexType> &vect, python::object &seq) {
  PySequenceHolder<IndexType> seqL(seq);
  for (unsigned int i = 0; i < seqL.size(); ++i) {
    IndexType idx = seqL[i];
    vect.setVal(idx, vect[idx] + 1);
  }
}
template <typename IndexType>
python::dict pyGetNonzeroElements(SparseIntVect<IndexType> &vect) {
  python::dict res;
  auto iter = vect.getNonzeroElements().begin();
  while (iter != vect.getNonzeroElements().end()) {
    res[iter->first] = iter->second;
    ++iter;
  }
  return res;
}

template <typename IndexType>
python::list pyToList(SparseIntVect<IndexType> &vect) {
  python::list res;
  res.append(0);
  res *= vect.getLength();
  for (auto iter : vect.getNonzeroElements()) {
    res[iter.first] = iter.second;
  }
  return res;
}

template <typename T>
python::list BulkDice(const T &siv1, python::list sivs, bool returnDistance) {
  python::list res;
  unsigned int nsivs = python::extract<unsigned int>(sivs.attr("__len__")());
  for (unsigned int i = 0; i < nsivs; ++i) {
    double simVal;
    const T *siv2 = python::extract<T *>(sivs[i])();
    simVal = DiceSimilarity(siv1, *siv2, returnDistance);
    res.append(simVal);
  }
  return res;
}
template <typename T>
python::list BulkTanimoto(const T &siv1, python::list sivs,
                          bool returnDistance) {
  python::list res;
  unsigned int nsivs = python::extract<unsigned int>(sivs.attr("__len__")());
  for (unsigned int i = 0; i < nsivs; ++i) {
    double simVal;
    const T *siv2 = python::extract<T *>(sivs[i])();
    simVal = TanimotoSimilarity(siv1, *siv2, returnDistance);
    res.append(simVal);
  }
  return res;
}

template <typename T>
python::list BulkTversky(const T &siv1, python::list sivs, double a, double b,
                         bool returnDistance) {
  python::list res;
  unsigned int nsivs = python::extract<unsigned int>(sivs.attr("__len__")());
  for (unsigned int i = 0; i < nsivs; ++i) {
    double simVal;
    const T *siv2 = python::extract<T *>(sivs[i])();
    simVal = TverskySimilarity(siv1, *siv2, a, b, returnDistance);
    res.append(simVal);
  }
  return res;
}
}  // namespace

std::string sparseIntVectDoc =
    "A container class for storing integer\n\
values within a particular range.\n\
\n\
The length of the vector is set at construction time.\n\
\n\
As you would expect, _SparseIntVects_ support a set of binary operations\n\
so you can do things like:\n\
  Arithmetic:\n\
  siv1 += siv2\n\
  siv3 = siv1 + siv2\n\
  siv1 -= siv3\n\
  siv3 = siv1 - siv2\n\
  \"Fuzzy\" binary operations:\n\
  siv3 = siv1 & siv2  the result contains the smallest value in each entry\n\
  siv3 = siv1 | siv2  the result contains the largest value in each entry\n\
\n\
Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])\n\
\n";

struct sparseIntVec_wrapper {
  template <typename IndexType>
  static void wrapOne(const char *className) {
    python::class_<SparseIntVect<IndexType>,
                   boost::shared_ptr<SparseIntVect<IndexType>>>(
        className, sparseIntVectDoc.c_str(),
        python::init<IndexType>("Constructor"))
        .def(python::init<std::string>(python::args("self", "pkl")))
        // Note: we cannot support __len__ because, at least at the moment
        // (BPL v1.34.1), it must return an int.
        .def("__setitem__", &SparseIntVect<IndexType>::setVal,
             python::args("self"), "Set the value at a specified location")
        .def("__getitem__", &SparseIntVect<IndexType>::getVal,
             python::args("self"),
             "Get the value at a specified location")
        .def(python::self & python::self)
        .def(python::self | python::self)
        .def(python::self - python::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(python::self -=
             python::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def(python::self + python::self)
        .def(python::self += python::self)
        .def(python::self == python::self)
        .def(python::self != python::self)
        //.def(python::self - int())
        .def(python::self -= int())
        //.def(python::self + int())
        .def(python::self += int())
        //.def(python::self / int())
        .def(python::self /= int())
        //.def(python::self * int())
        .def(python::self *= int())
        .def("GetTotalVal", &SparseIntVect<IndexType>::getTotalVal,
             ((python::args("self"), python::args("useAbs") = false)),
             "Get the sum of the values in the vector, basically L1 norm")
        .def("GetLength", &SparseIntVect<IndexType>::getLength,
             python::args("self"), "Returns the length of the vector")
        .def("ToBinary", &SIVToBinaryText<IndexType>, python::args("self"),
             "returns a binary (pickle) representation of the vector")
        .def("UpdateFromSequence", &pyUpdateFromSequence<IndexType>,
             python::args("self", "seq"),
             "update the vector based on the values in the list or tuple")
        .def("GetNonzeroElements", &pyGetNonzeroElements<IndexType>,
             python::args("self"),
             "returns a dictionary of the nonzero elements")
        .def("ToList", pyToList<IndexType>, python::args("self"),
             "Return the SparseIntVect as a python list")
        .def_pickle(siv_pickle_suite<IndexType>());

    python::def(
        "DiceSimilarity", &DiceSimilarity<IndexType>,
        (python::args("siv1"), python::args("siv2"),
         python::args("returnDistance") = false, python::args("bounds") = 0.0),
        "return the Dice similarity between two vectors");
    python::def("BulkDiceSimilarity", &BulkDice<SparseIntVect<IndexType>>,
                (python::args("v1"), python::args("v2"),
                 python::args("returnDistance") = false),
                "return the Dice similarities between one vector and a "
                "sequence of others");
    python::def(
        "TanimotoSimilarity", &TanimotoSimilarity<IndexType>,
        (python::args("siv1"), python::args("siv2"),
         python::args("returnDistance") = false, python::args("bounds") = 0.0),
        "return the Tanimoto similarity between two vectors");
    python::def("BulkTanimotoSimilarity",
                &BulkTanimoto<SparseIntVect<IndexType>>,
                (python::args("v1"), python::args("v2"),
                 python::args("returnDistance") = false),
                "return the Tanimoto similarities between one vector and a "
                "sequence of others");
    python::def("TverskySimilarity", &TverskySimilarity<IndexType>,
                (python::args("siv1"), python::args("siv2"), python::args("a"),
                 python::args("b"), python::args("returnDistance") = false,
                 python::args("bounds") = 0.0),
                "return the Tversky similarity between two vectors");
    python::def("BulkTverskySimilarity", &BulkTversky<SparseIntVect<IndexType>>,
                (python::args("v1"), python::args("v2"), python::args("a"),
                 python::args("b"), python::args("returnDistance") = false),
                "return the Tversky similarities between one vector and a "
                "sequence of others");
  }

  static void wrap() {
    wrapOne<std::int32_t>("IntSparseIntVect");
    wrapOne<boost::int64_t>("LongSparseIntVect");
    wrapOne<std::uint32_t>("UIntSparseIntVect");
    wrapOne<boost::uint64_t>("ULongSparseIntVect");
  }
};

void wrap_sparseIntVect() { sparseIntVec_wrapper::wrap(); }
