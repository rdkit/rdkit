//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_BUILD_THREADSAFE_SSS

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>

#include <string>

// ours
#include <GraphMol/FileParsers/MultithreadedSmilesMolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/Wrap/ContextManagers.h>

#include "MolSupplier.h"
#include "MultithreadedMolSupplier.h"

namespace python = boost::python;

namespace RDKit {
std::string multiSmilesMolSupplierClassDoc =
    "A class which concurrently supplies molecules from a text file.\n\
  Please note that this class is still a bit experimental and the API may\n\
  change in future releases.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules might not be constructed until we ask for them:\n\n\
       >>> suppl = MultithreadedSmilesMolSupplier('in.smi')\n\
       >>> for mol in suppl:\n\
       ...    if(mol):\n\
       ...      mol.GetNumAtoms()\n\
\n\
    2) Lazy evaluation 2:\n\n\
       >>> suppl = MultithreadedSmilesMolSupplier('in.smi')\n\
       >>> while (!suppl.atEnd()):\n\
       ...    mol = next(mol)\n\
       ...    if(mol):\n\
       ...      mol.GetNumAtoms()\n\
\n";

std::string multiSmsDocStr =
    "Constructor\n\n\
  ARGUMENTS: \n\
\n\
    - fileName: name of the file to be read\n\
\n\
    - delimiter: (optional) text delimiter (a string).  Defauts to ' \t'.\n\
\n\
    - smilesColumn: (optional) index of the column containing the SMILES\n\
      data.  Defaults to 0.\n\
\n\
    - nameColumn: (optional) index of the column containing molecule names.\n\
      Defaults to 1.\n\
\n\
    - titleLine: (optional) set this toggle if the file contains a title line.\n\
      Defaults to true.\n\
\n\
    - sanitize: (optional) toggles sanitization of molecules as they are read.\n\
      Defaults to true.\n\
\n\
    - numWriterThreads: (optional) number of writer threads. Defaults to 1.\n\
\n\
    - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.\n\
\n\
    - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.\n\
\n";

struct multiSmiMolSup_wrap {
  static void wrap() {
    python::class_<MultithreadedSmilesMolSupplier, boost::noncopyable>(
        "MultithreadedSmilesMolSupplier",
        multiSmilesMolSupplierClassDoc.c_str(), python::init<>())
        .def(python::init<std::string, std::string, int, int, bool, bool,
                          unsigned int, size_t, size_t>(
            (python::arg("fileName"), python::arg("delimiter") = " \t",
             python::arg("smilesColumn") = 0, python::arg("nameColumn") = 1,
             python::arg("titleLine") = true, python::arg("sanitize") = true,
             python::arg("numWriterThreads") = 1,
             python::arg("sizeInputQueue") = 5,
             python::arg("sizeOutputQueue") = 5),
            multiSmsDocStr.c_str()))
        .def("__iter__",
             (MultithreadedSmilesMolSupplier *
              (*)(MultithreadedSmilesMolSupplier *)) &
                 MTMolSupplIter,
             python::return_internal_reference<1>())
        .def("__enter__",
             (MultithreadedSmilesMolSupplier *
              (*)(MultithreadedSmilesMolSupplier *)) &
                 MolIOEnter,
             python::return_internal_reference<>())
        .def("__exit__",
             (bool (*)(MultithreadedSmilesMolSupplier *, python::object,
                       python::object, python::object)) &
                 MolIOExit)
        .def("__next__",
             (ROMol * (*)(MultithreadedSmilesMolSupplier *)) &
                 MolForwardSupplNext,
             "Returns the next molecule in the file. Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("atEnd", &MultithreadedSmilesMolSupplier::atEnd,
             "Returns true if we have read all records else false.\n")
        .def("GetLastRecordId",
             (unsigned int (*)(MultithreadedSmilesMolSupplier *)) &
                 MTMolSupplLastId,
             "Returns the record id for the last extracted item.\n")
        .def("GetLastItemText",
             (std::string(*)(MultithreadedSmilesMolSupplier *)) &
                 MTMolSupplLastItem,
             "Returns the text for the last extracted item.\n");
  };
};
}  // namespace RDKit

void wrap_multiSmiSupplier() { RDKit::multiSmiMolSup_wrap::wrap(); }
#endif
