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

// ours
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/python_streambuf.h>
#include <RDBoost/iterator_next.h>
#include <maeparser/Reader.hpp>

#include "MolSupplier.h"

namespace python = boost::python;

using boost_adaptbx::python::streambuf;
namespace {

class LocalMaeMolSupplier : public RDKit::MaeMolSupplier {
 public:
  LocalMaeMolSupplier(python::object &input, bool sanitize, bool removeHs) {
    // FIX: minor leak here
    auto *sb = new streambuf(input);
    dp_inStream = new streambuf::istream(*sb);
    dp_sInStream.reset(dp_inStream);
    df_owner = true;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    d_reader.reset(new schrodinger::mae::Reader(dp_sInStream));
    d_next_struct = d_reader->next("f_m_ct");
    POSTCONDITION(dp_inStream, "bad instream");
  }
  LocalMaeMolSupplier(streambuf &input, bool sanitize, bool removeHs) {
    dp_inStream = new streambuf::istream(input);
    dp_sInStream.reset(dp_inStream);
    df_owner = true;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    d_reader.reset(new schrodinger::mae::Reader(dp_sInStream));
    d_next_struct = d_reader->next("f_m_ct");
    POSTCONDITION(dp_inStream, "bad instream");
  }

  LocalMaeMolSupplier(const std::string &fname, bool sanitize = true,
                      bool removeHs = true) {
    df_owner = true;
    auto *ifs = new std::ifstream(fname.c_str(), std::ios_base::binary);
    if (!ifs || !(*ifs) || ifs->bad()) {
      std::ostringstream errout;
      errout << "Bad input file " << fname;
      throw RDKit::BadFileException(errout.str());
    }
    dp_inStream = (std::istream *)ifs;
    dp_sInStream.reset(dp_inStream);
    df_sanitize = sanitize;
    df_removeHs = removeHs;

    d_reader.reset(new schrodinger::mae::Reader(dp_sInStream));
    d_next_struct = d_reader->next("f_m_ct");
    POSTCONDITION(dp_inStream, "bad instream");
  };
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
       >>> for mol in suppl:\n \
       ...   if mol is not None: print mol.GetNumAtoms()\n\
\n\
  Properties in the Maestro file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
struct maemolsup_wrap {
  static void wrap() {
    python::class_<LocalMaeMolSupplier, boost::noncopyable>(
        "MaeMolSupplier", maeMolSupplierClassDoc.c_str(), python::no_init)
        .def(python::init<python::object &, bool, bool>(
            (python::arg("fileobj"), python::arg("sanitize") = true,
             python::arg("removeHs") =
                 true))[python::with_custodian_and_ward_postcall<0, 2>()])
        .def(python::init<streambuf &, bool, bool>(
            (python::arg("streambuf"), python::arg("sanitize") = true,
             python::arg("removeHs") =
                 true))[python::with_custodian_and_ward_postcall<0, 2>()])
        .def(python::init<std::string, bool, bool>(
            (python::arg("filename"), python::arg("sanitize") = true,
             python::arg("removeHs") = true)))
        .def(NEXT_METHOD, (ROMol * (*)(LocalMaeMolSupplier *)) & MolSupplNext,
             "Returns the next molecule in the file.  Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("atEnd", &MaeMolSupplier::atEnd,
             "Returns whether or not we have hit EOF.\n")
        .def("__iter__", &FwdMolSupplIter,
             python::return_internal_reference<1>());
  };
};
}  // namespace RDKit

void wrap_maesupplier() { RDKit::maemolsup_wrap::wrap(); }
