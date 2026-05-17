//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
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
#include <nanobind/stl/wstring.h>
#include <nanobind/stl/vector.h>
#include <string>

// ours
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/python_streambuf_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
using boost_adaptbx::python::streambuf;
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

class LocalSmilesWriter : public SmilesWriter {
 public:
  LocalSmilesWriter(nb::object fileObj, std::string delimiter,
                    std::string nameHeader, bool includeHeader,
                    bool isomericSmiles, bool kekuleSmiles)
      : SmilesWriter(new streambuf::ostream(new streambuf(fileObj, 't')),
                     delimiter, nameHeader, includeHeader, true, isomericSmiles,
                     kekuleSmiles) {}

  LocalSmilesWriter(std::string fileName, std::string delimiter,
                    std::string nameHeader, bool includeHeader,
                    bool isomericSmiles, bool kekuleSmiles)
      : SmilesWriter(fileName, delimiter, nameHeader, includeHeader,
                     isomericSmiles, kekuleSmiles) {}
};
}  // namespace

std::string swDocStr =
    R"DOC(Constructor.

   ARGUMENTS:

     - fileName: name of the output file. ('-' to write to stdout)
     - delimiter: (optional) delimiter to be used to separate entries on each line.
     - nameHeader: (optional) text to use for the name column in the header line.
                   If this is blank, names will not be included in the output.
     - includeHeader: (optional) toggles inclusion of a header line in the output file.
     - isomericSmiles: (optional) toggles output of isomeric smiles
       (includes stereochem information).
     - kekuleSmiles: (optional) toggles output of kekule smiles (no aromatic
       bonds for molecules that have been kekulized).

)DOC";
struct smiwriter_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalSmilesWriter>(
        m, "SmilesWriter",
        R"DOC(A class for writing molecules to text files.)DOC")
        .def(
            nb::init<std::string, std::string, std::string, bool, bool, bool>(),
            "fileName"_a, "delimiter"_a = " ", "nameHeader"_a = "Name",
            "includeHeader"_a = true, "isomericSmiles"_a = true,
            "kekuleSmiles"_a = false, swDocStr.c_str())
        .def(nb::init<nb::object, std::string, std::string, bool, bool, bool>(),
             "fileObj"_a, "delimiter"_a = " ", "nameHeader"_a = "Name",
             "includeHeader"_a = true, "isomericSmiles"_a = true,
             "kekuleSmiles"_a = false)
        .def("__enter__", &MolIOEnter<LocalSmilesWriter>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<LocalSmilesWriter>,
             "excType"_a = nb::none(), "excValue"_a = nb::none(),
             "traceback"_a = nb::none())
        .def("SetProps", &LocalSmilesWriter::setProps, "props"_a,
             R"DOC(Sets the properties to be written to the output file

  ARGUMENTS:

    - props: a list or tuple of property names

)DOC")
        .def("write", &SmilesWriter::write, "mol"_a, "confId"_a = -1,
             R"DOC(Writes a molecule to the output file.

  ARGUMENTS:

    - mol: the Mol to be written
    - confId: (optional) ignored

)DOC")
        .def(
            "flush", &SmilesWriter::flush,
            R"DOC(Flushes the output file (forces the disk file to be updated).)DOC")
        .def(
            "close", &SmilesWriter::close,
            R"DOC(Flushes the output file and closes it. The Writer cannot be used after this.)DOC")
        .def("NumMols", &SmilesWriter::numMols,
             R"DOC(Returns the number of molecules written so far.)DOC");
  };
};
}  // namespace RDKit

void wrap_smiwriter(nb::module_ &m) { RDKit::smiwriter_wrap::wrap(m); }
