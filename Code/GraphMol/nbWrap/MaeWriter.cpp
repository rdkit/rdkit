//
//  Copyright (C) 2023-2026 Schrodinger, LLC and other RDKit contributors
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

#include <memory>
#include <string>
#include <vector>

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

class LocalMaeWriter : public MaeWriter {
 public:
  LocalMaeWriter(nb::object fileobj)
      : dp_streambuf(new streambuf(fileobj, 't')) {
    dp_ostream.reset(new streambuf::ostream(*dp_streambuf));
  }

  LocalMaeWriter(const std::string &fname) : RDKit::MaeWriter(fname) {}

 private:
  std::unique_ptr<streambuf> dp_streambuf = nullptr;
};

}  // namespace

struct wrap_maewriter {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalMaeWriter>(
        m, "MaeWriter",
        R"DOC(An experimental class for writing molecules to Maestro files.

  Usage examples:

    1) writing to a named file:

       >>> writer = MaeWriter('out.mae')
       >>> for mol in list_of_mols:
       ...    writer.write(mol)

    2) writing to a file-like object:

       >>> import gzip
       >>> outf=gzip.open('out.mae.gz','wt+')
       >>> writer = MaeWriter(outf)
       >>> for mol in list_of_mols:
       ...   writer.write(mol)
       >>> writer.close()
       >>> outf.close()

  By default all non-private molecule, atom and bond properties are written
  to the Maestro file. This can be changed using the SetProps method:

       >>> writer = MaeWriter('out.mae')
       >>> writer.SetProps(['prop1','prop2'])

  Properties that are specified, but are not present will be ignored.

  Kekulization is mandatory, as the Maestro format does not have
  the concept of an aromatic bond.

  As this is an experimental writer, many features are not supported yet,
  e.g. chirality and bond stereo labels, stereo groups, substance groups,
  isotopes, or even dummy atoms. Note that these are not supported by
  MaeMolSupplier either.
)DOC")
        .def(nb::init<nb::object>(), "fileobj"_a)
        .def(nb::init<std::string>(), "filename"_a)
        .def("__enter__", &MolIOEnter<LocalMaeWriter>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<LocalMaeWriter>, "excType"_a = nb::none(),
             "excValue"_a = nb::none(), "traceback"_a = nb::none())
        .def(
            "SetProps", &LocalMaeWriter::setProps, "props_list"_a,
            R"DOC(Sets the atom and molecule properties to be written to the output file.

  ARGUMENTS:

    - props_list: a list of atom and molecule property names
)DOC")
        .def("write",
             (void (LocalMaeWriter::*)(const ROMol &,
                                       int))&LocalMaeWriter::write,
             "mol"_a, "confId"_a = defaultConfId,
             R"DOC(Writes a molecule to the output file.

  ARGUMENTS:

    - mol: the molecule to be written
    - confId: (optional) ID of the conformation to write
)DOC")
        .def(
            "flush", &LocalMaeWriter::flush,
            R"DOC(Flushes the output file (forces the disk file to be updated).)DOC")
        .def(
            "close", &LocalMaeWriter::close,
            R"DOC(Flushes the output file and closes it. The writer cannot be used after this.)DOC")
        .def("NumMols", &LocalMaeWriter::numMols,
             R"DOC(Returns the number of molecules written so far.)DOC")
        .def_static(
            "GetText", &LocalMaeWriter::getText, "mol"_a, "confId"_a = -1,
            "props_list"_a = std::vector<std::string>(),
            R"DOC(Returns the Maestro CT block text for a molecule.)DOC");
  };
};
}  // namespace RDKit

void wrap_maewriter(nb::module_ &m) { RDKit::wrap_maewriter::wrap(m); }
