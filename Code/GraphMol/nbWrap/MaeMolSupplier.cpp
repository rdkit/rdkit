//  Copyright (C) 2018-2026 Lorton and other RDKit contributors
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

#include <memory>
#include <string>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDBoost/python_streambuf_nb.h>
#include <RDGeneral/FileParseException.h>
#include "ContextManagers.h"

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;
using boost_adaptbx::python::streambuf;

namespace {

class LocalMaeMolSupplier : public RDKit::MaeMolSupplier {
 public:
  LocalMaeMolSupplier() : RDKit::MaeMolSupplier() {}

  LocalMaeMolSupplier(nb::object input, bool sanitize, bool removeHs)
      : dp_streambuf(new streambuf(input)) {
    auto inStream = new streambuf::istream(*dp_streambuf);
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

template <typename T>
ROMol *MolSupplNext(T *suppl) {
  ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    try {
      res = suppl->next();
    } catch (const FileParseException &) {
      throw;
    } catch (...) {
      res = nullptr;
    }
  } else {
    throw nb::stop_iteration();
  }
  return res;
}

template <typename T>
ROMol *MolSupplGetItem(T *suppl, int idx) {
  ROMol *res = nullptr;
  if (idx < 0) {
    idx = static_cast<int>(suppl->length()) + idx;
    if (idx < 0) {
      throw nb::index_error("invalid index");
    }
  }
  try {
    res = (*suppl)[static_cast<unsigned int>(idx)];
  } catch (...) {
    if (suppl->atEnd()) {
      throw nb::index_error("invalid index");
    } else {
      res = nullptr;
    }
  }
  return res;
}

}  // namespace

namespace RDKit {

std::string maeMolSupplierClassDoc =
    R"DOC(A class which supplies molecules from a file-like object containing Maestro data.

  Usage examples:

    1) Lazy evaluation: the molecules are not constructed until we ask for them:

       >>> suppl = MaeMolSupplier(file('in.mae'))
       >>> for mol in suppl:
       ...    if mol is not None: mol.GetNumAtoms()

    2) we can also read from compressed files:

       >>> import gzip
       >>> suppl = MaeMolSupplier(gzip.open('in.maegz'))
       >>> for mol in suppl:
       ...   if mol is not None: print mol.GetNumAtoms()

  Properties in the Maestro file are used to set properties on each molecule.
  The properties are accessible using the mol.GetProp(propName) method.

)DOC";

struct maemolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<LocalMaeMolSupplier>(m, "MaeMolSupplier",
                                    maeMolSupplierClassDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<std::string, bool, bool>(), "filename"_a,
             "sanitize"_a = true, "removeHs"_a = true)
        .def(nb::init<nb::object, bool, bool>(), "fileobj"_a,
             "sanitize"_a = true, "removeHs"_a = true)
        .def("__enter__", &MolIOEnter<LocalMaeMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<LocalMaeMolSupplier>)
        .def("__iter__", &FwdMolSupplIter, nb::rv_policy::reference_internal)
        .def(
            "__next__", &MolSupplNext<LocalMaeMolSupplier>,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC")
        .def("__getitem__", &MolSupplGetItem<LocalMaeMolSupplier>,
             nb::rv_policy::take_ownership, "idx"_a)
        .def("reset", &MaeMolSupplier::reset,
             R"DOC(Resets our position in the file to the beginning.
)DOC")
        .def("__len__", &MaeMolSupplier::length)
        .def("SetData", &MaeMolSupplier::setData, "data"_a, "sanitize"_a = true,
             "removeHs"_a = true, R"DOC(Sets the text to be parsed.)DOC")
        .def("atEnd", &MaeMolSupplier::atEnd,
             R"DOC(Returns whether or not we have hit EOF.
)DOC");
  };
};
}  // namespace RDKit

void wrap_maesupplier(nb::module_ &m) { RDKit::maemolsup_wrap::wrap(m); }
