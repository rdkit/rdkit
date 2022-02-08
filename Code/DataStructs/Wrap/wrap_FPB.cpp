//
//  Copyright (C) 2016 greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>
#include <DataStructs/FPBReader.h>
#include <DataStructs/MultiFPBReader.h>
#include <RDBoost/PySequenceHolder.h>
#include "wrap_helpers.h"

namespace python = boost::python;
using namespace RDKit;

namespace {
python::tuple taniNbrHelper(const FPBReader *self, const std::string &bytes,
                            double threshold) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  std::vector<std::pair<double, unsigned int>> nbrs =
      self->getTanimotoNeighbors(bv, threshold);
  python::list result;
  for (auto &nbr : nbrs) {
    result.append(python::make_tuple(nbr.first, nbr.second));
  }
  return python::tuple(result);
}
python::tuple tverskyNbrHelper(const FPBReader *self, const std::string &bytes,
                               double ca, double cb, double threshold) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  std::vector<std::pair<double, unsigned int>> nbrs =
      self->getTverskyNeighbors(bv, ca, cb, threshold);
  python::list result;
  for (auto &nbr : nbrs) {
    result.append(python::make_tuple(nbr.first, nbr.second));
  }
  return python::tuple(result);
}
python::tuple containingNbrHelper(const FPBReader *self,
                                  const std::string &bytes) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  std::vector<unsigned int> nbrs = self->getContainingNeighbors(bv);
  python::list result;
  for (auto &nbr : nbrs) {
    result.append(nbr);
  }
  return python::tuple(result);
}

python::tuple multiTaniNbrHelper(const MultiFPBReader *self,
                                 const std::string &bytes, double threshold,
                                 unsigned int numThreads) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  std::vector<MultiFPBReader::ResultTuple> nbrs =
      self->getTanimotoNeighbors(bv, threshold, numThreads);
  python::list result;
  for (auto &nbr : nbrs) {
    result.append(python::make_tuple(nbr.get<0>(), nbr.get<1>(), nbr.get<2>()));
  }
  return python::tuple(result);
}
python::tuple multiTverskyNbrHelper(const MultiFPBReader *self,
                                    const std::string &bytes, double ca,
                                    double cb, double threshold,
                                    unsigned int numThreads) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  std::vector<MultiFPBReader::ResultTuple> nbrs =
      self->getTverskyNeighbors(bv, ca, cb, threshold, numThreads);
  python::list result;
  for (auto &nbr : nbrs) {
    result.append(python::make_tuple(nbr.get<0>(), nbr.get<1>(), nbr.get<2>()));
  }
  return python::tuple(result);
}
python::tuple multiContainingNbrHelper(const MultiFPBReader *self,
                                       const std::string &bytes,
                                       unsigned int numThreads) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  std::vector<std::pair<unsigned int, unsigned int>> nbrs =
      self->getContainingNeighbors(bv, numThreads);
  python::list result;
  for (auto &nbr : nbrs) {
    result.append(python::make_tuple(nbr.first, nbr.second));
  }
  return python::tuple(result);
}

python::object getBytesHelper(const FPBReader *self, unsigned int which) {
  boost::shared_array<std::uint8_t> bv = self->getBytes(which);
  python::object retval =
      python::object(python::handle<>(PyBytes_FromStringAndSize(
          reinterpret_cast<const char *>(bv.get()), self->nBits() / 8)));
  return retval;
}

double getTaniHelper(const FPBReader *self, unsigned int which,
                     const std::string &bytes) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  return self->getTanimoto(which, bv);
}
python::tuple getItemHelper(const FPBReader *self, unsigned int which) {
  std::pair<boost::shared_ptr<ExplicitBitVect>, std::string> v = (*self)[which];
  return python::make_tuple(v.first, v.second);
}
double getTverskyHelper(const FPBReader *self, unsigned int which,
                        const std::string &bytes, double ca, double cb) {
  const auto *bv = reinterpret_cast<const std::uint8_t *>(bytes.c_str());
  return self->getTversky(which, bv, ca, cb);
}
}  // namespace

struct FPB_wrapper {
  static void wrap() {
    std::string FPBReaderClassDoc =
        "A class for reading and searching FPB files from Andrew Dalke's chemfp.\n\
    Note that this functionality is still experimental and the API may\n\
    change in future releases.\n";
    python::class_<FPBReader, boost::noncopyable>(
        "FPBReader", FPBReaderClassDoc.c_str(),
        python::init<std::string, python::optional<bool>>(
            (python::arg("filename"), python::arg("lazy") = false),
            "docstring"))
        .def("Init", &FPBReader::init,
             "Read the fingerprints from the file. This can take a while.\n")
        .def("__len__", &FPBReader::length)
        .def("__getitem__", &getItemHelper)
        .def("GetNumBits", &FPBReader::nBits,
             "returns the number of bits in a fingerprint")
        .def("GetFP", &FPBReader::getFP,
             "returns a particular fingerprint as an ExplicitBitVect")
        .def("GetBytes", &getBytesHelper,
             "returns a particular fingerprint as bytes")
        .def("GetId", &FPBReader::getId,
             "returns the id of a particular fingerprint")
        .def("GetTanimoto", &getTaniHelper,
             "return the tanimoto similarity of a particular fingerprint to "
             "the bytes provided")
        .def("GetTanimotoNeighbors", &taniNbrHelper,
             (python::arg("bv"), python::arg("threshold") = 0.7),
             "returns tanimoto similarities to and indices of all neighbors "
             "above the specified threshold")
        .def("GetTversky", &getTverskyHelper,
             "return the Tverksy similarity of a particular fingerprint to "
             "the bytes provided")
        .def("GetTverskyNeighbors", &tverskyNbrHelper,
             (python::arg("bv"), python::arg("ca"), python::arg("cb"),
              python::arg("threshold") = 0.7),
             "returns Tversky similarities to and indices of all neighbors "
             "above the specified threshold")
        .def(
            "GetContainingNeighbors", &containingNbrHelper, (python::arg("bv")),
            "returns indices of neighbors that contain this fingerprint (where "
            "all bits from this fingerprint are also set)");

    std::string MultiFPBReaderClassDoc =
        "A class for reading and searching multiple FPB files from Andrew Dalke's chemfp.\n\
    Note that this functionality is still experimental and the API may\n\
    change in future releases.\n";

    python::class_<MultiFPBReader, boost::noncopyable>(
        "MultiFPBReader", MultiFPBReaderClassDoc.c_str(),
        python::init<python::optional<bool>>(
            (python::arg("initOnSearch") = false), "docstring"))
        .def("Init", &MultiFPBReader::init,
             "Call Init() on each of our children. This can take a while.\n")
        .def("__len__", &MultiFPBReader::length)
        .def("GetNumBits", &MultiFPBReader::nBits,
             "returns the number of bits in a fingerprint")
        .def("AddReader", &MultiFPBReader::addReader,
             python::with_custodian_and_ward<1, 2>(),
             "adds an FPBReader to our set of readers")
        .def("GetReader", &MultiFPBReader::getReader,
             python::return_value_policy<python::reference_existing_object>(),
             "returns one of our readers")
        .def("GetTanimotoNeighbors", &multiTaniNbrHelper,
             (python::arg("bv"), python::arg("threshold") = 0.7,
              python::arg("numThreads") = 1),
             "returns tanimoto similarities to and indices of all neighbors "
             "above the specified threshold")
        .def("GetTverskyNeighbors", &multiTverskyNbrHelper,
             (python::arg("bv"), python::arg("ca"), python::arg("cb"),
              python::arg("threshold") = 0.7, python::arg("numThreads") = 1),
             "returns Tversky similarities to and indices of all neighbors "
             "above the specified threshold")
        .def(
            "GetContainingNeighbors", &multiContainingNbrHelper,
            (python::arg("bv"), python::arg("numThreads") = 1),
            "returns indices of neighbors that contain this fingerprint (where "
            "all bits from this fingerprint are also set)");
  }
};

void wrap_FPB() { FPB_wrapper::wrap(); }
