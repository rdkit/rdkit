//
//  Copyright (C) 2020-2026 Shrey Aryan and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_BUILD_THREADSAFE_SSS

#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

// ours
#include <GraphMol/FileParsers/MultithreadedSDMolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>
#include "ContextManagers.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {
template <typename T>
T *MTMolSupplIter(T *suppl) {
  return suppl;
}

template <typename T>
ROMol *MolForwardSupplNext(T *suppl) {
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
  if (suppl->atEnd() && suppl->getEOFHitOnRead()) {
    throw nb::stop_iteration();
  }
  return res;
}

template <typename T>
unsigned int MTMolSupplLastId(T *suppl) {
  return suppl->getLastRecordId();
}

template <typename T>
std::string MTMolSupplLastItem(T *suppl) {
  return suppl->getLastItemText();
}

}  // namespace

std::string multiSDMolSupplierClassDoc =
    R"DOC(A class which concurrently supplies molecules from an SD file.
  Please note that this class is still a bit experimental and the API may
  change in future releases.

  Usage examples:

    1) Lazy evaluation: the molecules might not be constructed until we ask for them:

       >>> suppl = MultithreadedSDMolSupplier('in.sdf')
       >>> for mol in suppl:
       ...    if(mol):
       ...      mol.GetNumAtoms()

    2) Lazy evaluation 2:

       >>> suppl = MultithreadedSDMolSupplier('in.sdf')
       >>> while (!suppl.atEnd()):
       ...    mol = next(mol)
       ...    if(mol):
       ...      mol.GetNumAtoms()

)DOC";

std::string multiSdsDocStr =
    R"DOC(Constructor

  ARGUMENTS:

    - fileName: name of the file to be read

    - sanitize: (optional) toggles sanitization of molecules as they are read.
      Defaults to true.

    - removeHs: (optional) removes Hs. Defaults to true.

    - strictParsing: (optional) allows strict or lax parsing. Defaults to true.

    - numWriterThreads: (optional) number of writer threads. Defaults to 1.

    - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.

    - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.

)DOC";

struct multiSDMolSup_wrap {
  static void wrap(nb::module_ &m) {
    nb::class_<MultithreadedSDMolSupplier>(m, "MultithreadedSDMolSupplier",
                                           multiSDMolSupplierClassDoc.c_str())
        .def(nb::init<>())
        .def(nb::init<std::string, bool, bool, bool, unsigned int, size_t,
                      size_t>(),
             "fileName"_a, "sanitize"_a = true, "removeHs"_a = true,
             "strictParsing"_a = true, "numWriterThreads"_a = 1,
             "sizeInputQueue"_a = 5, "sizeOutputQueue"_a = 5,
             multiSdsDocStr.c_str())
        .def("__iter__", &MTMolSupplIter<MultithreadedSDMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__enter__", &MolIOEnter<MultithreadedSDMolSupplier>,
             nb::rv_policy::reference_internal)
        .def("__exit__", &MolIOExit<MultithreadedSDMolSupplier>,
             "excType"_a = nb::none(), "excValue"_a = nb::none(),
             "traceback"_a = nb::none())
        .def(
            "__next__", &MolForwardSupplNext<MultithreadedSDMolSupplier>,
            nb::rv_policy::take_ownership,
            R"DOC(Returns the next molecule in the file. Raises _StopIteration_ on EOF.
)DOC")
        .def("atEnd", &MultithreadedSDMolSupplier::atEnd,
             R"DOC(Returns true if we have read all records else false.)DOC")
        .def("GetLastRecordId", &MTMolSupplLastId<MultithreadedSDMolSupplier>,
             R"DOC(Returns the record id for the last extracted item.)DOC")
        .def("GetLastItemText", &MTMolSupplLastItem<MultithreadedSDMolSupplier>,
             R"DOC(Returns the text for the last extracted item.)DOC")
        .def(
            "GetProcessPropertyLists",
            &MultithreadedSDMolSupplier::getProcessPropertyLists,
            R"DOC(Returns whether or not any property lists that are present will be processed when reading molecules.)DOC")
        .def(
            "SetProcessPropertyLists",
            &MultithreadedSDMolSupplier::setProcessPropertyLists, "val"_a,
            R"DOC(Sets whether or not any property lists that are present will be processed when reading molecules.)DOC");
  };
};
}  // namespace RDKit

void wrap_multiSDSupplier(nb::module_ &m) {
  RDKit::multiSDMolSup_wrap::wrap(m);
}
#endif
