//
//  Copyright (C) 2013-2026 Greg Landrum and other RDKit contributors
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
#include <string>

// ours
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RDKitBase.h>
#include "ContextManagers.h"
#include <RDBoost/python_streambuf_nb.h>

namespace nb = nanobind;
using namespace nb::literals;
using boost_adaptbx::python::streambuf;

namespace RDKit {
namespace {
class LocalPDBWriter : public PDBWriter {
 public:
  LocalPDBWriter(nb::object fileObj, unsigned int flavor = 0)
      : PDBWriter(new streambuf::ostream(new streambuf(fileObj, 't')), true,
                  flavor) {}

  LocalPDBWriter(std::string fileName, unsigned int flavor = 0)
      : PDBWriter(fileName, flavor) {}
};
}  // namespace

struct pdbwriter_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalPDBWriter>(
        m, "PDBWriter", R"DOC(A class for writing molecules to PDB files.)DOC")
        .def(nb::init<std::string, unsigned int>(), "fileName"_a,
             "flavor"_a = 0,
             R"DOC(Constructor.

  ARGUMENTS:

    - fileName: name of the output file. ('-' to write to stdout)
    - flavor: (optional)

)DOC")
        .def(nb::init<nb::object, unsigned int>(), "fileObj"_a, "flavor"_a = 0)
        .def("__enter__", &MolIOEnter<LocalPDBWriter>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<LocalPDBWriter>, "excType"_a = nb::none(),
             "excValue"_a = nb::none(), "traceback"_a = nb::none())
        .def("write", &PDBWriter::write, "mol"_a, "confId"_a = -1,
             R"DOC(Writes a molecule to the output file.

  ARGUMENTS:

    - mol: the Mol to be written
    - confId: (optional) ignored

)DOC")
        .def(
            "flush", &PDBWriter::flush,
            R"DOC(Flushes the output file (forces the disk file to be updated).)DOC")
        .def(
            "close", &PDBWriter::close,
            R"DOC(Flushes the output file and closes it. The Writer cannot be used after this.)DOC")
        .def("NumMols", &PDBWriter::numMols,
             R"DOC(Returns the number of molecules written so far.)DOC");
  };
};
}  // namespace RDKit

void wrap_pdbwriter(nb::module_ &m) { RDKit::pdbwriter_wrap::wrap(m); }
