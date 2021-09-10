//
//  Copyright (C) 2011-2019  Greg Landrum
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
#include "ContextManagers.h"

#include "MolSupplier.h"

namespace python = boost::python;

using boost_adaptbx::python::streambuf;
namespace {

class LocalForwardSDMolSupplier : public RDKit::ForwardSDMolSupplier {
 public:
  LocalForwardSDMolSupplier(python::object &input, bool sanitize, bool removeHs,
                            bool strictParsing) {
    // FIX: minor leak here
    auto *sb = new streambuf(input, 'b');
    dp_inStream = new streambuf::istream(*sb);
    df_owner = true;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    df_strictParsing = strictParsing;
    POSTCONDITION(dp_inStream, "bad instream");
  }
  LocalForwardSDMolSupplier(streambuf &input, bool sanitize, bool removeHs,
                            bool strictParsing) {
    dp_inStream = new streambuf::istream(input);
    df_owner = true;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    df_strictParsing = strictParsing;
    POSTCONDITION(dp_inStream, "bad instream");
  }
  LocalForwardSDMolSupplier(std::string filename, bool sanitize, bool removeHs,
                            bool strictParsing) {
    std::istream *tmpStream = nullptr;
    tmpStream = static_cast<std::istream *>(
        new std::ifstream(filename.c_str(), std::ios_base::binary));
    if (!(*tmpStream) || tmpStream->bad()) {
      delete tmpStream;
      std::ostringstream errout;
      errout << "Bad input file " << filename;
      throw RDKit::BadFileException(errout.str());
    }
    dp_inStream = tmpStream;
    df_owner = true;
    df_sanitize = sanitize;
    df_removeHs = removeHs;
    df_strictParsing = strictParsing;
    POSTCONDITION(dp_inStream, "bad instream");
  }
};

LocalForwardSDMolSupplier *FwdMolSupplIter(LocalForwardSDMolSupplier *self) {
  return self;
}
}  // namespace

namespace RDKit {

std::string fsdMolSupplierClassDoc =
    "A class which supplies molecules from file-like object containing SD data.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n\n\
       >>> suppl = ForwardSDMolSupplier(file('in.sdf'))\n\
       >>> for mol in suppl:\n\
       ...    if mol is not None: mol.GetNumAtoms()\n\
\n\
    2) we can also read from compressed files: \n\n\
       >>> import gzip\n\
       >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))\n\
       >>> for mol in suppl:\n\
       ...   if mol is not None: print mol.GetNumAtoms()\n\
\n\
  Properties in the SD file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
struct forwardsdmolsup_wrap {
  static void wrap() {
    python::class_<LocalForwardSDMolSupplier, boost::noncopyable>(
        "ForwardSDMolSupplier", fsdMolSupplierClassDoc.c_str(), python::no_init)
        .def(python::init<python::object &, bool, bool, bool>(
            (python::arg("fileobj"), python::arg("sanitize") = true,
             python::arg("removeHs") = true,
             python::arg("strictParsing") =
                 true))[python::with_custodian_and_ward_postcall<0, 2>()])
        .def(python::init<streambuf &, bool, bool, bool>(
            (python::arg("streambuf"), python::arg("sanitize") = true,
             python::arg("removeHs") = true,
             python::arg("strictParsing") =
                 true))[python::with_custodian_and_ward_postcall<0, 2>()])
        .def(python::init<std::string, bool, bool, bool>(
            (python::arg("filename"), python::arg("sanitize") = true,
             python::arg("removeHs") = true,
             python::arg("strictParsing") = true)))
        .def("__enter__",
             (LocalForwardSDMolSupplier * (*)(LocalForwardSDMolSupplier *)) &
                 MolIOEnter,
             python::return_internal_reference<>())
        .def("__exit__", (bool (*)(LocalForwardSDMolSupplier *, python::object,
                                   python::object, python::object)) &
                             MolIOExit)
        .def("__next__",
             (ROMol * (*)(LocalForwardSDMolSupplier *)) & MolForwardSupplNext,
             "Returns the next molecule in the file.  Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("atEnd", &ForwardSDMolSupplier::atEnd,
             "Returns whether or not we have hit EOF.\n")
        .def("GetEOFHitOnRead", &ForwardSDMolSupplier::getEOFHitOnRead,
             "Returns whether or EOF was hit while parsing the previous "
             "entry.\n")
        .def("__iter__", &FwdMolSupplIter,
             python::return_internal_reference<1>())
        .def("GetProcessPropertyLists",
             &ForwardSDMolSupplier::getProcessPropertyLists,
             "returns whether or not any property lists that are present will "
             "be processed when reading molecules")
        .def("SetProcessPropertyLists",
             &ForwardSDMolSupplier::setProcessPropertyLists,
             "sets whether or not any property lists that are present will be "
             "processed when reading molecules");
  };
};
}  // namespace RDKit

void wrap_forwardsdsupplier() { RDKit::forwardsdmolsup_wrap::wrap(); }
