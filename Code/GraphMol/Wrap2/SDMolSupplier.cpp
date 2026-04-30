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
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

// ours
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>

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

void setDataHelper(SDMolSupplier &self, const std::string &text, bool sanitize,
                   bool removeHs, bool strictParsing) {
  self.setData(text, sanitize, removeHs, strictParsing);
}

}  // namespace

void setStreamIndices(SDMolSupplier &self, const std::vector<int> &arg) {
  std::vector<std::streampos> loc;
  loc.reserve(arg.size());
  for (const auto &idx : arg) {
    loc.push_back(static_cast<std::streampos>(idx));
  }
  self.setStreamIndices(loc);
}

std::string sdMolSupplierClassDoc =
    R"DOC(A class which supplies molecules from an SD file.

     Usage examples:

     1) Lazy evaluation: the molecules are not constructed until we ask for them:

          >>> suppl = SDMolSupplier('in.sdf')
          >>> for mol in suppl:
          ...    mol.GetNumAtoms()

     2) Lazy evaluation 2:

          >>> suppl = SDMolSupplier('in.sdf')
          >>> mol1 = next(suppl)
          >>> mol2 = next(suppl)
          >>> suppl.reset()
          >>> mol3 = next(suppl)
          # mol3 and mol1 are the same:
          >>> MolToSmiles(mol3)==MolToSmiles(mol1)

     3) Random Access:

          >>> suppl = SDMolSupplier('in.sdf')
          >>> mol1 = suppl[0]
          >>> mol2 = suppl[1]
          # NOTE: this will generate an IndexError if the supplier doesn't have that many
          molecules.

     4) Random Access 2: looping over all molecules

          >>> suppl = SDMolSupplier('in.sdf')
          >>> nMols = len(suppl)
          >>> for i in range(nMols):
          ...   suppl[i].GetNumAtoms()

  Properties in the SD file are used to set properties on each molecule.
  The properties are accessible using the mol.GetProp(propName) method.
)DOC";
struct sdmolsup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<SDMolSupplier>(m, "SDMolSupplier", sdMolSupplierClassDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<std::string, bool, bool, bool>(), "fileName"_a,
             "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true)
        .def("__enter__", &MolIOEnter<SDMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<SDMolSupplier>)
        .def("__iter__", &MolSupplIter<SDMolSupplier>,
             nb::rv_policy::reference_internal)
        .def(
            "__next__", &MolSupplNext<SDMolSupplier>,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC")
        .def("__getitem__", &MolSupplGetItem<SDMolSupplier>,
             nb::rv_policy::take_ownership, "idx"_a)
        .def("reset", &SDMolSupplier::reset,
             R"DOC(Resets our position in the file to the beginning.
)DOC")
        .def("__len__", &SDMolSupplier::length)
        .def("SetData", setDataHelper, "data"_a, "sanitize"_a = true,
             "removeHs"_a = true, "strictParsing"_a = true,
             R"DOC(Sets the text to be parsed.)DOC")
        .def(
            "SetData",
            [](SDMolSupplier &suppl, nb::bytes data, bool sanitize,
               bool removeHs, bool strictParsing) {
              // Convert nb::bytes to std::string
              std::string text =
                  std::string(static_cast<const char *>(data.data()),
                              static_cast<size_t>(data.size()));
              suppl.setData(text, sanitize, removeHs, strictParsing);
            },
            "data"_a, "sanitize"_a = true, "removeHs"_a = true,
            "strictParsing"_a = true, R"DOC(Sets the text to be parsed.)DOC")
        .def(
            "_SetStreamIndices", setStreamIndices, "locs"_a,
            R"DOC(Sets the locations of mol beginnings in the input stream. Be *very* careful with this method.)DOC")
        .def("GetItemText", &SDMolSupplier::getItemText, "index"_a,
             R"DOC(Returns the text for an item.)DOC")
        .def("atEnd", &SDMolSupplier::atEnd,
             R"DOC(Returns whether or not we have hit EOF.
)DOC")
        .def(
            "GetProcessPropertyLists", &SDMolSupplier::getProcessPropertyLists,
            R"DOC(Returns whether or not any property lists that are present will be processed when reading molecules.)DOC")
        .def(
            "SetProcessPropertyLists", &SDMolSupplier::setProcessPropertyLists,
            "val"_a,
            R"DOC(Sets whether or not any property lists that are present will be processed when reading molecules.)DOC");
  };
};
}  // namespace RDKit

void wrap_sdsupplier(nb::module_ &m) { RDKit::sdmolsup_wrap::wrap(m); }
