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
        //.def("__getitem__", &FPBReader::operator[]) careful about this leaking
        .def("GetFP", &FPBReader::getFP,
             python::return_value_policy<python::manage_new_object>())
        .def("GetId", &FPBReader::getId)
        .def("GetTanimoto", &FPBReader::getTanimoto);
  }
};

void wrap_FPB() { FPB_wrapper::wrap(); }
