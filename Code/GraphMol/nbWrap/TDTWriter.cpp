//
//  Copyright (C) 2005-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <string>

// ours
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/python_streambuf_nb.h>

namespace nb = nanobind;
using namespace nb::literals;
using boost_adaptbx::python::streambuf;

namespace RDKit {
namespace {
template <typename T>
T *MolIOEnter(T *self) {
  return self;
}

template <typename T>
bool MolIOExit(T *self, nb::object exc_type, nb::object exc_val,
               nb::object traceback) {
  RDUNUSED_PARAM(exc_type);
  RDUNUSED_PARAM(exc_val);
  RDUNUSED_PARAM(traceback);
  self->close();
  return false;
}

class LocalTDTWriter : public TDTWriter {
 public:
  LocalTDTWriter(nb::object fileObj)
      : TDTWriter(new streambuf::ostream(new streambuf(fileObj, 't')), true) {}

  LocalTDTWriter(std::string fileName) : TDTWriter(fileName) {}
};
}  // namespace

struct tdtwriter_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalTDTWriter>(
        m, "TDTWriter", R"DOC(A class for writing molecules to TDT files.)DOC")
        .def(nb::init<std::string>(), "fileName"_a,
             R"DOC(Constructor.

   If a string argument is provided, it will be treated as the name of the
   output file. If a file-like object is provided, output will be sent there.

)DOC")
        .def(nb::init<nb::object>(), "fileObj"_a)
        .def("__enter__", &MolIOEnter<LocalTDTWriter>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<LocalTDTWriter>)
        .def("SetProps", &TDTWriter::setProps, "props"_a,
             R"DOC(Sets the properties to be written to the output file

  ARGUMENTS:

    - props: a list or tuple of property names

)DOC")
        .def("write", &TDTWriter::write, "mol"_a, "confId"_a = -1,
             R"DOC(Writes a molecule to the output file.

  ARGUMENTS:

    - mol: the Mol to be written
    - confId: (optional) ID of the conformation to write

)DOC")
        .def(
            "flush", &TDTWriter::flush,
            R"DOC(Flushes the output file (forces the disk file to be updated).)DOC")
        .def(
            "close", &TDTWriter::close,
            R"DOC(Flushes the output file and closes it. The Writer cannot be used after this.)DOC")
        .def("NumMols", &TDTWriter::numMols,
             R"DOC(Returns the number of molecules written so far.)DOC")
        .def(
            "SetWrite2D", &TDTWriter::setWrite2D, "state"_a = true,
            R"DOC(Causes 2D conformations to be written (default is 3D conformations).)DOC")
        .def("GetWrite2D", &TDTWriter::getWrite2D)
        .def(
            "SetWriteNames", &TDTWriter::setWriteNames, "state"_a = true,
            R"DOC(Causes names to be written to the output file as NAME records.)DOC")
        .def("GetWriteNames", &TDTWriter::getWriteNames)
        .def(
            "SetNumDigits", &TDTWriter::setNumDigits, "numDigits"_a,
            R"DOC(Sets the number of digits to be written for coordinates.)DOC")
        .def("GetNumDigits", &TDTWriter::getNumDigits);
  };
};
}  // namespace RDKit

void wrap_tdtwriter(nb::module_ &m) { RDKit::tdtwriter_wrap::wrap(m); }
