//
//  Copyright (C) 2013-2021 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#include <string>

// ours
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/RDKitBase.h>
#include "rdchem.h"
#include "ContextManagers.h"
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/python_streambuf.h>

namespace python = boost::python;

namespace RDKit {
using boost_adaptbx::python::streambuf;
namespace {
PDBWriter *getPDBWriter(python::object &fileobj, unsigned int flavor = 0) {
  // FIX: minor leak here
  auto *sb = new streambuf(fileobj, 't');
  auto *ost = new streambuf::ostream(*sb);
  return new PDBWriter(ost, true, flavor);
}
}  // namespace
std::string pdbwDocStr =
    "Constructor.\n\n"
    "   ARGUMENTS:\n\n"
    "     - fileName: name of the output file. ('-' to write to stdout)\n"
    "     - flavor: (optional) \n\n";
struct pdbwriter_wrap {
  static void wrap() {
    python::class_<PDBWriter>("PDBWriter",
                              "A class for writing molecules to PDB files.",
                              python::no_init)
        .def("__init__",
             python::make_constructor(
                 &getPDBWriter, python::default_call_policies(),
                 (python::arg("fileObj"), python::arg("flavor") = 0)))
        .def(python::init<std::string, unsigned int>(
            (python::arg("fileName"), python::arg("flavor") = 0),
            pdbwDocStr.c_str()))
        .def("__enter__", &MolIOEnter<PDBWriter>,
             python::return_internal_reference<>())
        .def("__exit__", &MolIOExit<PDBWriter>)

        .def("write", &PDBWriter::write,
             (python::arg("self"), python::arg("mol"),
              python::arg("confId") = -1),
             "Writes a molecule to the output file.\n\n"
             "  ARGUMENTS:\n\n"
             "    - mol: the Mol to be written\n"
             "    - confId: (optional) ignored \n\n")
        .def("flush", &PDBWriter::flush,
             "Flushes the output file (forces the disk file to be "
             "updated).\n\n")
        .def("close", &PDBWriter::close,
             "Flushes the output file and closes it. The Writer cannot be used "
             "after this.\n\n")
        .def("NumMols", &PDBWriter::numMols,
             "Returns the number of molecules written so far.\n\n");
  };
};
}  // namespace RDKit

void wrap_pdbwriter() { RDKit::pdbwriter_wrap::wrap(); }
