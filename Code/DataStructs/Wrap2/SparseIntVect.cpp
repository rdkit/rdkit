//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstdint>
#include <RDGeneral/Invariant.h>
#include <DataStructs/SparseIntVect.h>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

namespace {
template <typename IndexType>
nb::bytes SIVToBinaryText(const SparseIntVect<IndexType> &siv) {
  std::string res = siv.toString();
  return nb::bytes(res.c_str(), res.length());
}
}  // namespace

namespace {
template <typename IndexType>
void pyUpdateFromSequence(SparseIntVect<IndexType> &vect,
                          const nb::iterable &seq) {
  for (auto item : seq) {
    IndexType idx = nb::cast<IndexType>(item);
    vect.setVal(idx, vect[idx] + 1);
  }
}

template <typename IndexType>
nb::dict pyGetNonzeroElements(SparseIntVect<IndexType> &vect) {
  nb::dict res;
  auto iter = vect.getNonzeroElements().begin();
  while (iter != vect.getNonzeroElements().end()) {
    res[nb::cast(iter->first)] = nb::cast(iter->second);
    ++iter;
  }
  return res;
}

template <typename IndexType>
nb::list pyToList(SparseIntVect<IndexType> &vect) {
  nb::list res;
  for (IndexType i = 0; i < vect.getLength(); ++i) {
    res.append(0);
  }
  for (auto iter : vect.getNonzeroElements()) {
    res[static_cast<size_t>(iter.first)] = iter.second;
  }
  return res;
}

template <typename T>
nb::list BulkDice(const T &siv1, const nb::iterable &sivs,
                  bool returnDistance) {
  nb::list res;
  for (auto siv : sivs) {
    const auto &siv2 = nb::cast<const T &>(siv);
    auto simVal = DiceSimilarity(siv1, siv2, returnDistance);
    res.append(simVal);
  }
  return res;
}

template <typename T>
nb::list BulkTanimoto(const T &siv1, const nb::iterable &sivs,
                      bool returnDistance) {
  nb::list res;
  for (auto siv : sivs) {
    const auto &siv2 = nb::cast<const T &>(siv);
    auto simVal = TanimotoSimilarity(siv1, siv2, returnDistance);
    res.append(simVal);
  }
  return res;
}

template <typename T>
nb::list BulkTversky(const T &siv1, const nb::iterable &sivs, double a,
                     double b, bool returnDistance) {
  nb::list res;
  for (auto siv : sivs) {
    const auto &siv2 = nb::cast<const T &>(siv);
    auto simVal = TverskySimilarity(siv1, siv2, a, b, returnDistance);
    res.append(simVal);
  }
  return res;
}
}  // namespace

std::string sparseIntVectDoc =
    R"DOC(A container class for storing integer
values within a particular range.

The length of the vector is set at construction time.

As you would expect, _SparseIntVects_ support a set of binary operations
so you can do things like:
  Arithmetic:
  siv1 += siv2
  siv3 = siv1 + siv2
  siv1 -= siv3
  siv3 = siv1 - siv2
  "Fuzzy" binary operations:
  siv3 = siv1 & siv2  the result contains the smallest value in each entry
  siv3 = siv1 | siv2  the result contains the largest value in each entry

Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])

)DOC";

struct sparseIntVec_wrapper {
  template <typename IndexType>
  static void wrapOne(nb::module_ &m, const char *className) {
    nb::class_<SparseIntVect<IndexType>>(m, className, sparseIntVectDoc.c_str())
        .def(nb::init<IndexType>(), R"DOC(Constructor)DOC")
        .def(
            "__init__",
            [](SparseIntVect<IndexType> *t, nb::bytes b) {
              new (t) SparseIntVect<IndexType>(
                  std::string(static_cast<const char *>(b.data()),
                              static_cast<size_t>(b.size())));
            },
            "pkl"_a)
        // Note: we cannot support __len__ because, at least at the moment
        // (BPL v1.34.1), it must return an int.
        .def("__setitem__", &SparseIntVect<IndexType>::setVal,
             R"DOC(Set the value at a specified location)DOC")
        .def("__getitem__", &SparseIntVect<IndexType>::getVal,
             R"DOC(Get the value at a specified location)DOC")
        .def(nb::self & nb::self)
        .def(nb::self | nb::self)
        .def(nb::self - nb::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(nb::self -= nb::self)  // clang warns incorrectly on this construct
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def(nb::self + nb::self)
        .def(nb::self += nb::self)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self -= int())
        .def(nb::self += int())
        .def(nb::self /= int())
        .def(nb::self *= int())
        .def(
            "GetTotalVal", &SparseIntVect<IndexType>::getTotalVal,
            "useAbs"_a = false,
            R"DOC(Get the sum of the values in the vector, basically L1 norm)DOC")
        .def("GetLength", &SparseIntVect<IndexType>::getLength,
             R"DOC(Returns the length of the vector)DOC")
        .def("ToBinary", &SIVToBinaryText<IndexType>,
             R"DOC(returns a binary (pickle) representation of the vector)DOC")
        .def(
            "UpdateFromSequence", &pyUpdateFromSequence<IndexType>, "seq"_a,
            R"DOC(update the vector based on the values in the list or tuple)DOC")
        .def("GetNonzeroElements", &pyGetNonzeroElements<IndexType>,
             R"DOC(returns a dictionary of the nonzero elements)DOC")
        .def("ToList", pyToList<IndexType>,
             R"DOC(Return the SparseIntVect as a python list)DOC")
        .def("__getstate__",
             [](const SparseIntVect<IndexType> &siv) {
               const auto pkl = siv.toString();
               return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
             })
        .def("__setstate__",
             [](SparseIntVect<IndexType> &siv,
                const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&siv) SparseIntVect<IndexType>(pkl);
             })
        .def("__setstate__",
             [](SparseIntVect<IndexType> &siv,
                const std::tuple<std::string> &state) {
               new (&siv) SparseIntVect<IndexType>(std::get<0>(state));
             })
        .doc() = sparseIntVectDoc.c_str();

    m.def("DiceSimilarity", &DiceSimilarity<IndexType>, "siv1"_a, "siv2"_a,
          "returnDistance"_a = false, "bounds"_a = 0.0,
          R"DOC(return the Dice similarity between two vectors)DOC");
    m.def(
        "BulkDiceSimilarity", &BulkDice<SparseIntVect<IndexType>>, "v1"_a,
        "v2"_a, "returnDistance"_a = false,
        R"DOC(return the Dice similarities between one vector and a sequence of others)DOC");
    m.def("TanimotoSimilarity", &TanimotoSimilarity<IndexType>, "siv1"_a,
          "siv2"_a, "returnDistance"_a = false, "bounds"_a = 0.0,
          R"DOC(return the Tanimoto similarity between two vectors)DOC");
    m.def(
        "BulkTanimotoSimilarity", &BulkTanimoto<SparseIntVect<IndexType>>,
        "v1"_a, "v2"_a, "returnDistance"_a = false,
        R"DOC(return the Tanimoto similarities between one vector and a sequence of others)DOC");
    m.def("TverskySimilarity", &TverskySimilarity<IndexType>, "siv1"_a,
          "siv2"_a, "a"_a, "b"_a, "returnDistance"_a = false, "bounds"_a = 0.0,
          R"DOC(return the Tversky similarity between two vectors)DOC");
    m.def(
        "BulkTverskySimilarity", &BulkTversky<SparseIntVect<IndexType>>, "v1"_a,
        "v2"_a, "a"_a, "b"_a, "returnDistance"_a = false,
        R"DOC(return the Tversky similarities between one vector and a sequence of others)DOC");
  }

  static void wrap(nb::module_ &m) {
    wrapOne<std::int32_t>(m, "IntSparseIntVect");
    wrapOne<std::int64_t>(m, "LongSparseIntVect");
    wrapOne<std::uint32_t>(m, "UIntSparseIntVect");
    wrapOne<std::uint64_t>(m, "ULongSparseIntVect");
  }
};

void wrap_sparseIntVect(nb::module_ &m) { sparseIntVec_wrapper::wrap(m); }
