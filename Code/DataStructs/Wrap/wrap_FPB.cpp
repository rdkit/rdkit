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
#include <RDBoost/PySequenceHolder.h>
#include "wrap_helpers.h"

namespace python = boost::python;
using namespace RDKit;

namespace {
python::tuple taniNbrHelper(const FPBReader *self, const std::string &bytes,
                            double threshold, unsigned int topN) {
  const boost::uint8_t *bv = (const boost::uint8_t *)bytes.c_str();
  std::vector<std::pair<double, unsigned int> > nbrs =
      self->getTanimotoNeighbors(bv, threshold, topN);
  python::list result;
  for (unsigned int i = 0; i < nbrs.size(); ++i) {
    result.append(python::make_tuple(nbrs[i].first, nbrs[i].second));
  }
  return python::tuple(result);
}

python::object getBytesHelper(const FPBReader *self, unsigned int which) {
  boost::uint8_t *bv = self->getBytes(which);
  python::object retval = python::object(python::handle<>(
      PyBytes_FromStringAndSize((const char *)bv, self->nBits() / 8)));
  return retval;
}

double getTaniHelper(const FPBReader *self, unsigned int which,
                     const std::string &bytes) {
  const boost::uint8_t *bv = (const boost::uint8_t *)bytes.c_str();
  return self->getTanimoto(which, bv);
}
}

std::string FPBReaderClassDoc =
    "A class for read-only interactions with FPB files from Andrew Dalke's chemfp.\n\
\n\
\n";

struct FPB_wrapper {
  static void wrap() {
    python::class_<FPBReader, boost::noncopyable>(
        "FPBReader", FPBReaderClassDoc.c_str(), python::init<std::string>())
        .def("Init", &FPBReader::init, "init.\n")
        .def("__len__", &FPBReader::length)
        .def("GetNumBits", &FPBReader::nBits)

        //.def("__getitem__", &FPBReader::operator[]) careful about this
        // leaking
        .def("GetFP", &FPBReader::getFP,
             python::return_value_policy<python::manage_new_object>())
        .def("GetBytes", &getBytesHelper)
        .def("GetId", &FPBReader::getId)
        .def("GetTanimoto", &getTaniHelper)

        .def("GetTanimotoNeighbors", &taniNbrHelper,
             (python::arg("bv"), python::arg("threshold") = 0.7,
              python::arg("topN") = 0))

        ;
  }
};

void wrap_FPB() { FPB_wrapper::wrap(); }
