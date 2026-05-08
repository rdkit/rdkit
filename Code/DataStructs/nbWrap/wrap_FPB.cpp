//
//  Copyright (C) 2026 greg Landrum and other RDKit contributors
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <DataStructs/FPBReader.h>
#include <RDGeneral/Invariant.h>
#include <DataStructs/MultiFPBReader.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {
const std::uint8_t *bytesToFP(const nb::bytes &bytes) {
  return reinterpret_cast<const std::uint8_t *>(bytes.data());
}

nb::tuple toTuple(const nb::list &l) {
  return nb::steal<nb::tuple>(PySequence_Tuple(l.ptr()));
}

nb::tuple taniNbrHelper(const FPBReader *self, const nb::bytes &bytes,
                        double threshold) {
  const auto *bv = bytesToFP(bytes);
  std::vector<std::pair<double, unsigned int>> nbrs =
      self->getTanimotoNeighbors(bv, threshold);
  nb::list result;
  for (auto &nbr : nbrs) {
    result.append(nb::make_tuple(nbr.first, nbr.second));
  }
  return toTuple(result);
}

nb::tuple tverskyNbrHelper(const FPBReader *self, const nb::bytes &bytes,
                           double ca, double cb, double threshold) {
  const auto *bv = bytesToFP(bytes);
  std::vector<std::pair<double, unsigned int>> nbrs =
      self->getTverskyNeighbors(bv, ca, cb, threshold);
  nb::list result;
  for (auto &nbr : nbrs) {
    result.append(nb::make_tuple(nbr.first, nbr.second));
  }
  return toTuple(result);
}

nb::tuple containingNbrHelper(const FPBReader *self, const nb::bytes &bytes) {
  const auto *bv = bytesToFP(bytes);
  std::vector<unsigned int> nbrs = self->getContainingNeighbors(bv);
  nb::list result;
  for (auto &nbr : nbrs) {
    result.append(nbr);
  }
  return toTuple(result);
}

nb::tuple multiTaniNbrHelper(const MultiFPBReader *self, const nb::bytes &bytes,
                             double threshold, unsigned int numThreads) {
  const auto *bv = bytesToFP(bytes);
  std::vector<MultiFPBReader::ResultTuple> nbrs =
      self->getTanimotoNeighbors(bv, threshold, numThreads);
  nb::list result;
  for (auto &nbr : nbrs) {
    result.append(
        nb::make_tuple(std::get<0>(nbr), std::get<1>(nbr), std::get<2>(nbr)));
  }
  return toTuple(result);
}

nb::tuple multiTverskyNbrHelper(const MultiFPBReader *self,
                                const nb::bytes &bytes, double ca, double cb,
                                double threshold, unsigned int numThreads) {
  const auto *bv = bytesToFP(bytes);
  std::vector<MultiFPBReader::ResultTuple> nbrs =
      self->getTverskyNeighbors(bv, ca, cb, threshold, numThreads);
  nb::list result;
  for (auto &nbr : nbrs) {
    result.append(
        nb::make_tuple(std::get<0>(nbr), std::get<1>(nbr), std::get<2>(nbr)));
  }
  return toTuple(result);
}

nb::tuple multiContainingNbrHelper(const MultiFPBReader *self,
                                   const nb::bytes &bytes,
                                   unsigned int numThreads) {
  const auto *bv = bytesToFP(bytes);
  std::vector<std::pair<unsigned int, unsigned int>> nbrs =
      self->getContainingNeighbors(bv, numThreads);
  nb::list result;
  for (auto &nbr : nbrs) {
    result.append(nb::make_tuple(nbr.first, nbr.second));
  }
  return toTuple(result);
}

nb::bytes getBytesHelper(const FPBReader *self, unsigned int which) {
  boost::shared_array<std::uint8_t> bv = self->getBytes(which);
  return nb::bytes(reinterpret_cast<const char *>(bv.get()), self->nBits() / 8);
}

double getTaniHelper(const FPBReader *self, unsigned int which,
                     const nb::bytes &bytes) {
  const auto *bv = bytesToFP(bytes);
  return self->getTanimoto(which, bv);
}

ExplicitBitVect getFPHelper(const FPBReader *self, unsigned int which) {
  auto fp = self->getFP(which);
  if (!fp) {
    throw nb::value_error("null fingerprint pointer");
  }
  return *fp;
}

nb::tuple getItemHelper(const FPBReader *self, unsigned int which) {
  std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> v = (*self)[which];
  if (!v.first) {
    throw nb::value_error("null fingerprint pointer");
  }
  return nb::make_tuple(*v.first, v.second);
}

double getTverskyHelper(const FPBReader *self, unsigned int which,
                        const nb::bytes &bytes, double ca, double cb) {
  const auto *bv = bytesToFP(bytes);
  return self->getTversky(which, bv, ca, cb);
}
}  // namespace

struct FPB_wrapper {
  static void wrap(nb::module_ &m) {
    std::string FPBReaderClassDoc =
        R"DOC(A class for reading and searching FPB files from Andrew Dalke's chemfp.
    Note that this functionality is still experimental and the API may
    change in future releases.
)DOC";
    nb::class_<FPBReader>(m, "FPBReader", FPBReaderClassDoc.c_str())
        .def(nb::init<const std::string &, bool>(), "filename"_a,
             "lazy"_a = false, R"DOC(docstring)DOC")
        .def("Init", &FPBReader::init,
             R"DOC(Read the fingerprints from the file. This can take a while.
)DOC")
        .def("__len__", &FPBReader::length)
        .def("__getitem__", &getItemHelper, "which"_a)
        .def("GetNumBits", &FPBReader::nBits,
             R"DOC(returns the number of bits in a fingerprint)DOC")
        .def("GetFP", &getFPHelper, "idx"_a,
             R"DOC(returns a particular fingerprint as an ExplicitBitVect)DOC")
        .def("GetBytes", &getBytesHelper, "which"_a,
             R"DOC(returns a particular fingerprint as bytes)DOC")
        .def("GetId", &FPBReader::getId, "idx"_a,
             R"DOC(returns the id of a particular fingerprint)DOC")
        .def(
            "GetTanimoto", &getTaniHelper, "which"_a, "bytes"_a,
            R"DOC(return the tanimoto similarity of a particular fingerprint to the bytes provided)DOC")
        .def(
            "GetTanimotoNeighbors", &taniNbrHelper, "bv"_a, "threshold"_a = 0.7,
            R"DOC(returns tanimoto similarities to and indices of all neighbors above the specified threshold)DOC")
        .def(
            "GetTversky", &getTverskyHelper, "which"_a, "bytes"_a, "ca"_a,
            "cb"_a,
            R"DOC(return the Tverksy similarity of a particular fingerprint to the bytes provided)DOC")
        .def(
            "GetTverskyNeighbors", &tverskyNbrHelper, "bv"_a, "ca"_a, "cb"_a,
            "threshold"_a = 0.7,
            R"DOC(returns Tversky similarities to and indices of all neighbors above the specified threshold)DOC")
        .def(
            "GetContainingNeighbors", &containingNbrHelper, "bv"_a,
            R"DOC(returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set))DOC");

    std::string MultiFPBReaderClassDoc =
        R"DOC(A class for reading and searching multiple FPB files from Andrew Dalke's chemfp.
    Note that this functionality is still experimental and the API may
    change in future releases.
)DOC";

    nb::class_<MultiFPBReader>(m, "MultiFPBReader",
                               MultiFPBReaderClassDoc.c_str())
        .def(nb::init<bool>(), "initOnSearch"_a = false, R"DOC(docstring)DOC")
        .def("Init", &MultiFPBReader::init,
             R"DOC(Call Init() on each of our children. This can take a while.
)DOC")
        .def("__len__", &MultiFPBReader::length)
        .def("GetNumBits", &MultiFPBReader::nBits,
             R"DOC(returns the number of bits in a fingerprint)DOC")
        .def("AddReader", &MultiFPBReader::addReader, nb::keep_alive<1, 2>(),
             "rdr"_a, R"DOC(adds an FPBReader to our set of readers)DOC")
        .def("GetReader", &MultiFPBReader::getReader, "which"_a,
             nb::rv_policy::reference_internal,
             R"DOC(returns one of our readers)DOC")
        .def(
            "GetTanimotoNeighbors", &multiTaniNbrHelper, "bv"_a,
            "threshold"_a = 0.7, "numThreads"_a = 1,
            R"DOC(returns tanimoto similarities to and indices of all neighbors above the specified threshold)DOC")
        .def(
            "GetTverskyNeighbors", &multiTverskyNbrHelper, "bv"_a, "ca"_a,
            "cb"_a, "threshold"_a = 0.7, "numThreads"_a = 1,
            R"DOC(returns Tversky similarities to and indices of all neighbors above the specified threshold)DOC")
        .def(
            "GetContainingNeighbors", &multiContainingNbrHelper, "bv"_a,
            "numThreads"_a = 1,
            R"DOC(returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set))DOC");
  }
};

void wrap_FPB(nb::module_ &m) { FPB_wrapper::wrap(m); }
