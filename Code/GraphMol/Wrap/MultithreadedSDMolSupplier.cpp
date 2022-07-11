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

#include <fstream>
#include <string>
// ours
#include <GraphMol/FileParsers/MultithreadedSDMolSupplier.h>
#include <GraphMol/Wrap/ContextManagers.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>

#include "MolSupplier.h"
#include "MultithreadedMolSupplier.h"

namespace python = boost::python;
using boost_adaptbx::python::streambuf;

namespace RDKit {

std::string multiSDMolSupplierClassDoc =
    "A class which concurrently supplies molecules from a text file.\n\
  Please note that this class is still a bit experimental and the API may\n\
  change in future releases.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules might not be constructed until we ask for them:\n\n\
       >>> suppl = MultithreadedSDMolSupplier('in.sdf')\n\
       >>> for mol in suppl:\n\
       ...    if(mol):\n\
       ...      mol.GetNumAtoms()\n\
\n\
    2) Lazy evaluation 2:\n\n\
       >>> suppl = MultithreadedSDMolSupplier('in.sdf')\n\
       >>> while (!suppl.atEnd()):\n\
       ...    mol = next(mol)\n\
       ...    if(mol):\n\
       ...      mol.GetNumAtoms()\n\
\n";

std::string multiSdsDocStr =
    "Constructor\n\n\
  ARGUMENTS: \n\
\n\
    - fileName: name of the file to be read\n\
\n\
    - sanitize: (optional) toggles sanitization of molecules as they are read.\n\
      Defaults to true.\n\
\n\
    - removeHs: (optional) removes Hs. Defaults to true.\n\
\n\
    - strictParsing: (optional) allows strict or lax parsing. Defaults to true.\n\
\n\
    - numWriterThreads: (optional) number of writer threads. Defaults to 1.\n\
\n\
    - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.\n\
\n\
    - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.\n\
\n";
#if 0
// FIX: disabled until we figure out how to make this stable
MultithreadedSDMolSupplier* MTMolSupplStream(
    python::object& input, bool sanitize = true, bool removeHs = true,
    bool strictParsing = true, unsigned int numWriterThreads = 1,
    size_t sizeInputQueue = 5, size_t sizeOutputQueue = 5) {
  auto* sb = new streambuf(input, 'b');
  auto* inStream = new streambuf::istream(*sb);
  MultithreadedSDMolSupplier* sup = new MultithreadedSDMolSupplier(
      inStream, true, sanitize, removeHs, strictParsing, numWriterThreads,
      sizeInputQueue, sizeOutputQueue);
  return sup;
}
#endif
struct multiSDMolSup_wrap {
  static void wrap() {
#if 0
    // FIX: disabled until we make it stable and figure out an API that we're happy with
    python::def(
        "SDMolSupplierFromStream", MTMolSupplStream,
        "Returns MultithreadedSDMolSupplier object constructed from a file object or stream",
        (python::arg("fileobj"), python::arg("sanitize") = true,
         python::arg("removeHs") = true, python::arg("strictParsing") = true,
         python::arg("numWriterThreads") = 1, python::arg("sizeInputQueue") = 5,
         python::arg("sizeOutputQueue") = 5),
        python::with_custodian_and_ward_postcall<
            0, 1, python::return_value_policy<python::manage_new_object>>());
#endif

    python::class_<MultithreadedSDMolSupplier, boost::noncopyable>(
        "MultithreadedSDMolSupplier", multiSDMolSupplierClassDoc.c_str(),
        python::init<>())
        .def(python::init<std::string, bool, bool, bool, unsigned int, size_t,
                          size_t>(
            (python::arg("fileName"), python::arg("sanitize") = true,
             python::arg("removeHs") = true,
             python::arg("strictParsing") = true,
             python::arg("numWriterThreads") = 1,
             python::arg("sizeInputQueue") = 5,
             python::arg("sizeOutputQueue") = 5),
            multiSdsDocStr.c_str()))
        .def("__iter__",
             (MultithreadedSDMolSupplier * (*)(MultithreadedSDMolSupplier *)) &
                 MTMolSupplIter,
             python::return_internal_reference<1>())
        .def("__enter__",
             (MultithreadedSDMolSupplier * (*)(MultithreadedSDMolSupplier *)) &
                 MolIOEnter,
             python::return_internal_reference<>())
        .def("__exit__", (bool (*)(MultithreadedSDMolSupplier *, python::object,
                                   python::object, python::object)) &
                             MolIOExit)
        .def("__next__",
             (ROMol * (*)(MultithreadedSDMolSupplier *)) & MolForwardSupplNext,
             "Returns the next molecule in the file. Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("atEnd", &MultithreadedSDMolSupplier::atEnd,
             "Returns true if we have read all records else false.\n")
        .def(
            "GetLastRecordId",
            (unsigned int (*)(MultithreadedSDMolSupplier *)) & MTMolSupplLastId,
            "Returns the record id for the last extracted item.\n")
        .def(
            "GetLastItemText",
            (std::string(*)(MultithreadedSDMolSupplier *)) & MTMolSupplLastItem,
            "Returns the text for the last extracted item.\n")
        .def("GetProcessPropertyLists",
             &MultithreadedSDMolSupplier::getProcessPropertyLists,
             "returns whether or not any property lists that are present will "
             "be processed when reading molecules")
        .def("SetProcessPropertyLists",
             &MultithreadedSDMolSupplier::setProcessPropertyLists,
             "sets whether or not any property lists that are present will be "
             "processed when reading molecules");
  };
};
}  // namespace RDKit

void wrap_multiSDSupplier() { RDKit::multiSDMolSup_wrap::wrap(); }
#endif