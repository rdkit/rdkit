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

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <fstream>
#include <memory>
#include <string>

// ours
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/FileParseException.h>
#include <RDBoost/python_streambuf_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

using boost_adaptbx::python::streambuf;
namespace {

class LocalForwardSDMolSupplier : public RDKit::ForwardSDMolSupplier {
 private:
  std::unique_ptr<streambuf> dp_streambuf;

 public:
  LocalForwardSDMolSupplier(nb::object input, bool sanitize, bool removeHs,
                            bool strictParsing) {
    dp_streambuf.reset(new streambuf(input, 'b'));
    auto sbis = new streambuf::istream(*dp_streambuf);
    bool owner = true;

    RDKit::v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    dp_supplier.reset(
        new RDKit::v2::FileParsers::ForwardSDMolSupplier(sbis, owner, params));
    POSTCONDITION(sbis, "bad instream");
  }
  LocalForwardSDMolSupplier(streambuf &input, bool sanitize, bool removeHs,
                            bool strictParsing) {
    auto sbis = new streambuf::istream(input);
    bool owner = true;

    RDKit::v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    dp_supplier.reset(
        new RDKit::v2::FileParsers::ForwardSDMolSupplier(sbis, owner, params));
    POSTCONDITION(sbis, "bad instream");
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
    bool owner = true;
    RDKit::v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = strictParsing;
    dp_supplier.reset(new RDKit::v2::FileParsers::ForwardSDMolSupplier(
        tmpStream, owner, params));
    POSTCONDITION(tmpStream, "bad instream");
  }
};

LocalForwardSDMolSupplier *FwdMolSupplIter(LocalForwardSDMolSupplier *self) {
  return self;
}

template <typename T>
RDKit::ROMol *MolForwardSupplNext(T *suppl) {
  RDKit::ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (const RDKit::FileParseException &) {
      throw;
    } catch (...) {
      res = nullptr;
    }
  } else {
    throw nb::stop_iteration();
  }
  if (suppl->atEnd() && suppl->getEOFHitOnRead()) {
    throw nb::stop_iteration();
  }
  return res;
}

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
}  // namespace

namespace RDKit {

std::string fsdMolSupplierClassDoc =
    R"DOC(A class which supplies molecules from a file-like object containing SD data.

  Usage examples:

    1) Lazy evaluation: the molecules are not constructed until we ask for them:

       >>> suppl = ForwardSDMolSupplier(file('in.sdf'))
       >>> for mol in suppl:
       ...    if mol is not None: mol.GetNumAtoms()

    2) we can also read from compressed files:

       >>> import gzip
       >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))
       >>> for mol in suppl:
       ...   if mol is not None: print mol.GetNumAtoms()

  Properties in the SD file are used to set properties on each molecule.
  The properties are accessible using the mol.GetProp(propName) method.

)DOC";
struct forwardsdmolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalForwardSDMolSupplier>(m, "ForwardSDMolSupplier",
                                          fsdMolSupplierClassDoc.c_str())
        .def(nb::init<std::string, bool, bool, bool>(), "filename"_a,
             "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true)
        .def(nb::init<nb::object, bool, bool, bool>(), "fileobj"_a,
             "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true)
        .def("__enter__",
             (LocalForwardSDMolSupplier * (*)(LocalForwardSDMolSupplier *)) &
                 MolIOEnter,
             nb::rv_policy::reference_internal)
        .def("__exit__",
             (bool (*)(LocalForwardSDMolSupplier *, nb::object, nb::object,
                       nb::object))&MolIOExit,
             "excType"_a = nb::none(), "excValue"_a = nb::none(),
             "traceback"_a = nb::none())
        .def(
            "__next__",
            (ROMol * (*)(LocalForwardSDMolSupplier *)) & MolForwardSupplNext,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC")
        .def("atEnd", &ForwardSDMolSupplier::atEnd,
             R"DOC(Returns whether or not we have hit EOF.
)DOC")
        .def("GetEOFHitOnRead", &ForwardSDMolSupplier::getEOFHitOnRead,
             R"DOC(Returns whether EOF was hit while parsing the previous entry.
)DOC")
        .def("__iter__", &FwdMolSupplIter, nb::rv_policy::reference_internal)
        .def(
            "GetProcessPropertyLists",
            &ForwardSDMolSupplier::getProcessPropertyLists,
            R"DOC(Returns whether or not any property lists that are present will be processed when reading molecules.)DOC")
        .def(
            "SetProcessPropertyLists",
            &ForwardSDMolSupplier::setProcessPropertyLists, "val"_a,
            R"DOC(Sets whether or not any property lists that are present will be processed when reading molecules.)DOC");
  };
};
}  // namespace RDKit

void wrap_forwardsdsupplier(nb::module_ &m) {
  RDKit::forwardsdmolsup_wrap::wrap(m);
}
