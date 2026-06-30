//
//  Copyright (C) 2005-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

// ours
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include "ContextManagers.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {
template <typename T>
T *MolSupplIter(T *suppl) {
  suppl->reset();
  return suppl;
}

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

std::string tdtMolSupplierClassDoc =
    R"DOC(A class which supplies molecules from a TDT file.

  Usage examples:

    1) Lazy evaluation: the molecules are not constructed until we ask for them:

       >>> suppl = TDTMolSupplier('in.smi')
       >>> for mol in suppl:
       ...    mol.GetNumAtoms()

    2) Lazy evaluation 2:

       >>> suppl = TDTMolSupplier('in.smi')
       >>> mol1 = next(suppl)
       >>> mol2 = next(suppl)
       >>> suppl.reset()
       >>> mol3 = next(suppl)
       # mol3 and mol1 are the same:
       >>> MolToSmiles(mol3)==MolToSmiles(mol1)

    3) Random Access:  all molecules are constructed as soon as we ask for the
       length:

       >>> suppl = TDTMolSupplier('in.smi')
       >>> nMols = len(suppl)
       >>> for i in range(nMols):
       ...   suppl[i].GetNumAtoms()

  Properties in the file are used to set properties on each molecule.
  The properties are accessible using the mol.GetProp(propName) method.
)DOC";

struct tdtmolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<TDTMolSupplier>(m, "TDTMolSupplier",
                               tdtMolSupplierClassDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<std::string, std::string, int, int, bool>(), "fileName"_a,
             "nameRecord"_a = "", "confId2D"_a = -1, "confId3D"_a = -1,
             "sanitize"_a = true)
        .def("__enter__", &MolIOEnter<TDTMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<TDTMolSupplier>, "excType"_a = nb::none(),
             "excValue"_a = nb::none(), "traceback"_a = nb::none())
        .def("__iter__", &MolSupplIter<TDTMolSupplier>,
             nb::rv_policy::reference_internal)
        .def(
            "__next__", &MolSupplNext<TDTMolSupplier>,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC")
        .def("__getitem__", &MolSupplGetItem<TDTMolSupplier>,
             nb::rv_policy::take_ownership, "idx"_a)
        .def("reset", &TDTMolSupplier::reset,
             R"DOC(Resets our position in the file to the beginning.
)DOC")
        .def("__len__", &TDTMolSupplier::length)
        .def("SetData", &TDTMolSupplier::setData, "data"_a, "nameRecord"_a = "",
             "confId2D"_a = -1, "confId3D"_a = -1, "sanitize"_a = true,
             R"DOC(Sets the text to be parsed.)DOC")
        .def("GetItemText", &TDTMolSupplier::getItemText, "index"_a,
             R"DOC(Returns the text for an item.)DOC");
  };
};
}  // namespace RDKit

void wrap_tdtsupplier(nb::module_ &m) { RDKit::tdtmolsup_wrap::wrap(m); }
