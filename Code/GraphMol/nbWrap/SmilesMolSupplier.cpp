//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
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
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>

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

SmilesMolSupplier *SmilesSupplierFromText(
    std::string text, std::string delimiter = " ", int smilesColumn = 0,
    int nameColumn = 1, bool titleLine = true, bool sanitize = true) {
  auto *res = new SmilesMolSupplier();
  res->setData(text, delimiter, smilesColumn, nameColumn, titleLine, sanitize);
  return res;
}

std::string smilesMolSupplierClassDoc =
    R"DOC(A class which supplies molecules from a text file.

  Usage examples:

    1) Lazy evaluation: the molecules are not constructed until we ask for them:

       >>> suppl = SmilesMolSupplier('in.smi')
       >>> for mol in suppl:
       ...    mol.GetNumAtoms()

    2) Lazy evaluation 2:

       >>> suppl = SmilesMolSupplier('in.smi')
       >>> mol1 = next(suppl)
       >>> mol2 = next(suppl)
       >>> suppl.reset()
       >>> mol3 = next(suppl)
       # mol3 and mol1 are the same:
       >>> MolToSmiles(mol3)==MolToSmiles(mol1)

    3) Random Access: all molecules are constructed as soon as we ask for the
       length:

       >>> suppl = SmilesMolSupplier('in.smi')
       >>> nMols = len(suppl)
       >>> for i in range(nMols):
       ...   suppl[i].GetNumAtoms()

  If the input file has a title line and more than two columns (smiles and id), the
  additional columns will be used to set properties on each molecule. The properties
  are accessible using the mol.GetProp(propName) method.
)DOC";

std::string smsDocStr =
    R"DOC(Constructor

  ARGUMENTS:

    - fileName: name of the file to be read

    - delimiter: (optional) text delimiter (a string). Defauts to ' '.

    - smilesColumn: (optional) index of the column containing the SMILES
      data. Defaults to 0.

    - nameColumn: (optional) index of the column containing molecule names.
      Defaults to 1.

    - titleLine: (optional) set this toggle if the file contains a title line.
      Defaults to 1.

    - sanitize: (optional) toggles sanitization of molecules as they are read.
      Defaults to 1.
)DOC";
struct smimolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<SmilesMolSupplier>(m, "SmilesMolSupplier",
                                  smilesMolSupplierClassDoc.c_str())
        .def(nb::init<std::string, std::string, int, int, bool, bool>(),
             "data"_a, "delimiter"_a = " ", "smilesColumn"_a = 0,
             "nameColumn"_a = 1, "titleLine"_a = true, "sanitize"_a = true,
             smsDocStr.c_str())
        .def(nb::init<>())
        .def("__enter__", &MolIOEnter<SmilesMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<SmilesMolSupplier>,
             "excType"_a = nb::none(), "excValue"_a = nb::none(),
             "traceback"_a = nb::none())
        .def("__iter__", &MolSupplIter<SmilesMolSupplier>,
             nb::rv_policy::reference_internal)
        .def(
            "__next__", &MolSupplNext<SmilesMolSupplier>,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC")
        .def("__getitem__", &MolSupplGetItem<SmilesMolSupplier>,
             nb::rv_policy::take_ownership, "idx"_a)
        .def("reset", &SmilesMolSupplier::reset,
             R"DOC(Resets our position in the file to the beginning.
)DOC")
        .def("__len__", &SmilesMolSupplier::length)
        .def("SetData", &SmilesMolSupplier::setData, "data"_a,
             "delimiter"_a = " ", "smilesColumn"_a = 0, "nameColumn"_a = 1,
             "titleLine"_a = true, "sanitize"_a = true,
             R"DOC(Sets the text to be parsed.)DOC")
        .def("GetItemText", &SmilesMolSupplier::getItemText, "index"_a,
             R"DOC(Returns the text for an item.)DOC");

    m.def("SmilesMolSupplierFromText", SmilesSupplierFromText, "text"_a,
          "delimiter"_a = " ", "smilesColumn"_a = 0, "nameColumn"_a = 1,
          "titleLine"_a = true, "sanitize"_a = true,
          nb::rv_policy::take_ownership);
  }
};
}  // namespace RDKit

void wrap_smisupplier(nb::module_ &m) { RDKit::smimolsup_wrap::wrap(m); }
