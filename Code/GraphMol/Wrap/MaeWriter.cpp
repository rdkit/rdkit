//
//
//  Copyright (C) 2023 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY

#include <memory>
#include <string>

#include <maeparser/Writer.hpp>

#include <GraphMol/FileParsers/MolWriters.h>
#include <RDBoost/python.h>
#include <RDBoost/python_streambuf.h>
#include <RDBoost/Wrap.h>

#include "ContextManagers.h"

namespace python = boost::python;
using boost_adaptbx::python::streambuf;

namespace RDKit {

class LocalMaeWriter : public MaeWriter {
 public:
  LocalMaeWriter(python::object &fileobj)
      : dp_streambuf(new streambuf(fileobj, 't')) {
    dp_ostream.reset(new streambuf::ostream(*dp_streambuf));
  }

  LocalMaeWriter(streambuf &output) {
    dp_ostream.reset(new streambuf::ostream(output));
  }

  LocalMaeWriter(const std::string &fname) : RDKit::MaeWriter(fname) {}

 private:
  std::unique_ptr<streambuf> dp_streambuf = nullptr;
};

struct wrap_maewriter {
  static void wrap() {
    std::string docStr =
        "An experimental class for writing molecules to Maestro files.\n\
\n\
  Usage examples:\n\
\n\
    1) writing to a named file:\n\n\
       >>> writer = MaeWriter('out.mae')\n\
       >>> for mol in list_of_mols:\n\
       ...    writer.write(mol)\n\
\n\
    2) writing to a file-like object: \n\n\
       >>> import gzip\n\
       >>> outf=gzip.open('out.mae.gz','wt+')\n\
       >>> writer = MaeWriter(outf)\n\
       >>> for mol in list_of_mols:\n\
       ...   writer.write(mol)\n\
       >>> writer.close()\n\
       >>> outf.close()\n\
\n\
  By default all non-private molecule, atom and bond properties are written\n\
  to the Maestro file. This can be changed using the SetProps method:\n\n\
       >>> writer = MaeWriter('out.mae')\n\
       >>> writer.SetProps(['prop1','prop2'])\n\n\
  Properties that are specified, but are not present will be ignored.\n\n\
  Kekulization is mandatory, as the Maestro format does not have\n\
  the concept of an aromatic bond\n\n\
  As this is an experimental writer, many features are not supported yet,\n\
  e.g. chirality and bond stereo labels, stereo groups, substance groups,\n\
  isotopes, or even dummy atoms. Note that these aren't supported by\n\
  MaeMolSupplier either.\n\
\n ";
    python::class_<LocalMaeWriter, boost::noncopyable>(
        "MaeWriter", docStr.c_str(), python::no_init)
        .def(python::init<python::object &>(python::arg("fileobj")))
        .def(python::init<streambuf &>(python::arg("streambuf")))
        .def(python::init<std::string>(python::arg("filename")))
        .def("__enter__", &MolIOEnter<LocalMaeWriter>,
             python::return_internal_reference<>())
        .def("__exit__", &MolIOExit<LocalMaeWriter>)
        .def(
            "SetProps", &LocalMaeWriter::setProps,
            (python::arg("self"), python::args("props_list")),
            "Sets the atom and mol properties to be written to the output file\n\n"
            "  ARGUMENTS:\n\n"
            "    - props: a list or tuple of atom and mol property names\n\n")
        .def(
            "write",
            (void(LocalMaeWriter::*)(const ROMol &, const std::string &, int)) &
                LocalMaeWriter::write,
            (python::arg("self"), python::arg("mol"),
             python::arg("heavyAtomColor") = defaultMaeHeavyAtomColor,
             python::arg("confId") = defaultConfId),
            "Writes a molecule to the output file.\n\n"
            "  ARGUMENTS:\n\n"
            "    - mol: the Mol to be written\n"
            "    - heavyAtomColor: (optional) color which heavy atoms will have in Maestro\n"
            "    - confId: (optional) ID of the conformation to write\n\n")

        .def("flush", &LocalMaeWriter::flush,
             "Flushes the output file (forces the disk file to be "
             "updated).\n\n")
        .def("close", &LocalMaeWriter::close,
             "Flushes the output file and closes it. The Writer cannot be used "
             "after this.\n\n")
        .def("NumMols", &LocalMaeWriter::numMols,
             "Returns the number of molecules written so far.\n\n")
        .def("GetText", &LocalMaeWriter::getText,
             (python::arg("mol"),
              python::arg("heavyAtomColor") = defaultMaeHeavyAtomColor,
              python::arg("confId") = -1,
              python::arg("props_list") = std::vector<std::string>()),
             "returns the Maestro ct block text for a molecule")
        .staticmethod("GetText");

    iterable_converter().from_python<std::vector<std::string>>();
  };
};
}  // namespace RDKit

void wrap_maewriter() { RDKit::wrap_maewriter::wrap(); }
