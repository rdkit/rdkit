//  Copyright (C) 2018  Lorton
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
#include <fstream>
#include <memory>

#include <maeparser/MaeConstants.hpp>
#include <maeparser/Reader.hpp>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/python_streambuf.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

#include "MolSupplier.h"
#include "ContextManagers.h"

namespace python = boost::python;

using namespace schrodinger;
using boost_adaptbx::python::streambuf;
namespace {

bool streamIsGoodOrExhausted(std::istream *stream) {
  PRECONDITION(stream, "bad stream");
  return stream->good() || (stream->eof() && stream->fail() && !stream->bad());
}

class LocalMaeMolSupplier : public RDKit::MaeMolSupplier {
 public:
  LocalMaeMolSupplier() : RDKit::MaeMolSupplier() {}

  LocalMaeMolSupplier(python::object &input, bool sanitize, bool removeHs)
      : dp_streambuf(new streambuf(input)) {
    auto inStream = new streambuf::istream(*dp_streambuf);
    bool owner = true;
    RDKit::v2::FileParsers::MaeMolSupplierParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    dp_supplier.reset(
        new RDKit::v2::FileParsers::MaeMolSupplier(inStream, owner, params));
  }
  LocalMaeMolSupplier(streambuf &input, bool sanitize, bool removeHs) {
    auto inStream = new streambuf::istream(input);

    bool owner = true;
    RDKit::v2::FileParsers::MaeMolSupplierParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    dp_supplier.reset(
        new RDKit::v2::FileParsers::MaeMolSupplier(inStream, owner, params));
  }

  LocalMaeMolSupplier(const std::string &fname, bool sanitize = true,
                      bool removeHs = true)
      : RDKit::MaeMolSupplier(fname, sanitize, removeHs) {}

 private:
  std::unique_ptr<streambuf> dp_streambuf = nullptr;
};

LocalMaeMolSupplier *FwdMolSupplIter(LocalMaeMolSupplier *self) { return self; }

}  // namespace

namespace RDKit {

std::string maeMolSupplierClassDoc =
    "A class which supplies molecules from file-like object containing Maestro data.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n\n\
       >>> suppl = MaeMolSupplier(file('in.mae'))\n\
       >>> for mol in suppl:\n\
       ...    if mol is not None: mol.GetNumAtoms()\n\
\n\
    2) we can also read from compressed files: \n\n\
       >>> import gzip\n\
       >>> suppl = MaeMolSupplier(gzip.open('in.maegz'))\n\
       >>> for mol in suppl:\n\
       ...   if mol is not None: print mol.GetNumAtoms()\n\
\n\
  Properties in the Maestro file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
struct maemolsup_wrap {
  static void wrap() {
    python::class_<LocalMaeMolSupplier, boost::noncopyable>(
        "MaeMolSupplier", maeMolSupplierClassDoc.c_str(),
        python::init<>(python::args("self")))
        .def(python::init<python::object &, bool, bool>(
            (python::arg("self"), python::arg("fileobj"),
             python::arg("sanitize") = true,
             python::arg("removeHs") =
                 true))[python::with_custodian_and_ward_postcall<0, 2>()])
        .def(python::init<streambuf &, bool, bool>(
            (python::arg("self"), python::arg("streambuf"),
             python::arg("sanitize") = true,
             python::arg("removeHs") =
                 true))[python::with_custodian_and_ward_postcall<0, 2>()])
        .def(python::init<std::string, bool, bool>(
            (python::arg("self"), python::arg("filename"),
             python::arg("sanitize") = true, python::arg("removeHs") = true)))
        .def("__enter__", &MolIOEnter<LocalMaeMolSupplier>,
             python::return_internal_reference<>())
        .def("__exit__", &MolIOExit<LocalMaeMolSupplier>)
        .def("__iter__", &FwdMolSupplIter,
             python::return_internal_reference<1>(), python::args("self"))
        .def("__next__", &MolSupplNext<LocalMaeMolSupplier>,
             "Returns the next molecule in the file.  Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>(),
             python::args("self"))
        .def("__getitem__", &MolSupplGetItem<LocalMaeMolSupplier>,
             python::return_value_policy<python::manage_new_object>(),
             python::args("self", "idx"))
        .def("reset", &MaeMolSupplier::reset, python::args("self"),
             "Resets our position in the file to the beginning.\n")
        .def("__len__", &MaeMolSupplier::length, python::args("self"))
        .def("SetData", &MaeMolSupplier::setData, "Sets the text to be parsed",
             ((python::arg("self"), python::arg("data")),
              python::arg("sanitize") = true, python::arg("removeHs") = true))
        .def("atEnd", &MaeMolSupplier::atEnd, python::args("self"),
             "Returns whether or not we have hit EOF.\n");
  };
};
}  // namespace RDKit

void wrap_maesupplier() { RDKit::maemolsup_wrap::wrap(); }
