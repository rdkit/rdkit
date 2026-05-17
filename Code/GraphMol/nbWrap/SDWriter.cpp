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

class LocalSDWriter : public SDWriter {
 public:
  LocalSDWriter(nb::object fileObj)
      : SDWriter(new streambuf::ostream(new streambuf(fileObj, 't')), true) {}

  LocalSDWriter(std::string fileName) : SDWriter(fileName) {}
};
}  // namespace

std::string getSDTextHelper(const ROMol &mol, int confId, bool kekulize,
                            bool force_V3000, int molid) {
  return SDWriter::getText(mol, confId, kekulize, force_V3000, molid);
}

struct sdwriter_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalSDWriter>(m, "SDWriter",
                              R"DOC(A class for writing molecules to SD files.

  Usage examples:

    1) writing to a named file:

       >>> writer = SDWriter('out.sdf')
       >>> for mol in list_of_mols:
       ...    writer.write(mol)

    2) writing to a file-like object:

       >>> import gzip
       >>> outf=gzip.open('out.sdf.gz','wt+')
       >>> writer = SDWriter(outf)
       >>> for mol in list_of_mols:
       ...   writer.write(mol)
       >>> writer.close()
       >>> outf.close()

  By default all non-private molecular properties are written to the SD file.
  This can be changed using the SetProps method:

       >>> writer = SDWriter('out.sdf')
       >>> writer.SetProps(['prop1','prop2'])

)DOC")
        .def(nb::init<std::string>(), "fileName"_a,
             R"DOC(Constructor.

If a string argument is provided, it will be treated as the name of the
output file. If a file-like object is provided, output will be sent there.

)DOC")
        .def(nb::init<nb::object>(), "fileObj"_a)
        .def("__enter__", &MolIOEnter<LocalSDWriter>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<LocalSDWriter>, "excType"_a = nb::none(),
             "excValue"_a = nb::none(), "traceback"_a = nb::none())
        .def("SetProps", &SDWriter::setProps, "props"_a,
             R"DOC(Sets the properties to be written to the output file

  ARGUMENTS:

    - props: a list or tuple of property names

)DOC")
        .def("write", &SDWriter::write, "mol"_a, "confId"_a = -1,
             R"DOC(Writes a molecule to the output file.

  ARGUMENTS:

    - mol: the Mol to be written
    - confId: (optional) ID of the conformation to write

)DOC")
        .def(
            "flush", &SDWriter::flush,
            R"DOC(Flushes the output file (forces the disk file to be updated).)DOC")
        .def(
            "close", &SDWriter::close,
            R"DOC(Flushes the output file and closes it. The Writer cannot be used after this.)DOC")
        .def("NumMols", &SDWriter::numMols,
             R"DOC(Returns the number of molecules written so far.)DOC")
        .def(
            "SetForceV3000", &SDWriter::setForceV3000, "val"_a,
            R"DOC(Sets whether or not V3000 mol file writing is being forced.)DOC")
        .def(
            "GetForceV3000", &SDWriter::getForceV3000,
            R"DOC(Returns whether or not V3000 mol file writing is being forced.)DOC")
        .def("SetKekulize", &SDWriter::setKekulize, "val"_a,
             R"DOC(Sets whether or not molecules are kekulized on writing.)DOC")
        .def(
            "GetKekulize", &SDWriter::getKekulize,
            R"DOC(Returns whether or not molecules are kekulized on writing.)DOC")
        .def_static("GetText", &getSDTextHelper, "mol"_a, "confId"_a = -1,
                    "kekulize"_a = true, "force_v3000"_a = false,
                    "molid"_a = -1,
                    R"DOC(returns the SD text for a molecule)DOC");
  };
};
}  // namespace RDKit

void wrap_sdwriter(nb::module_ &m) { RDKit::sdwriter_wrap::wrap(m); }
