//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/iterator_next.h>
#include <RDBoost/python.h>

#include <string>

// ours
#include <GraphMol/FileParsers/MultithreadedSmilesMolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>

#include "MultithreadedMolSupplier.h"

namespace python = boost::python;

namespace RDKit {

std::string multiSmilesMolSupplierClassDoc =
    "A class which concurrently supplies molecules from a text file.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules might not be constructed until we ask for them:\n\n\
       >>> suppl = MultithreadedSmilesMolSupplier('in.smi')\n\
       >>> while not suppl.atEnd():\n\
			 ...    mol = next(suppl)\n\
       ...    mol.GetNumAtoms()\n\
\n\
  If the input file has a title line and more than two columns (smiles and id), the\n\
  additional columns will be used to set properties on each molecule.  The properties\n\
  are accessible using the mol.GetProp(propName) method.\n\
\n";

std::string multiSmsDocStr =
    "Constructor\n\n\
  ARGUMENTS: \n\
\n\
    - fileName: name of the file to be read\n\
\n\
    - delimiter: (optional) text delimiter (a string).  Defauts to ' '.\n\
\n\
    - smilesColumn: (optional) index of the column containing the SMILES\n\
      data.  Defaults to 0.\n\
\n\
    - nameColumn: (optional) index of the column containing molecule names.\n\
      Defaults to 1.\n\
\n\
    - titleLine: (optional) set this toggle if the file contains a title line.\n\
      Defaults to 1.\n\
\n\
    - sanitize: (optional) toggles sanitization of molecules as they are read.\n\
      Defaults to 1.\n\
\n\
    - numWriterThreads: (optional) number of writer threads\n\
\n\
    - sizeInputQueue: (optional) size of input queue\n\
\n\
    - sizeOutputQueue: (optional) size of output queue\n\
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
             python::arg("numWriterThreads") = 2,
             python::arg("sizeInputQueue") = 5,
             python::arg("sizeOutputQueue") = 5),
            multiSmsDocStr.c_str()))
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
                 MolSupplLastId,
             "Returns the record id for the last extracted item.\n")
        .def("GetItemText",
             (std::string(*)(MultithreadedSmilesMolSupplier *)) &
                 MolSupplLastItem,
             "Returns the text for the last extracted item.\n");
  };
};
}  // namespace RDKit

void wrap_multiSmiSupplier() { RDKit::multiSmiMolSup_wrap::wrap(); }
