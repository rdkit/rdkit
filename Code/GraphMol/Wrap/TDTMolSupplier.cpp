//
//  Copyright (C) 2005-2021 Greg Landrum and Rational Discovery LLC
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

// ours
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include "MolSupplier.h"
#include "ContextManagers.h"

namespace python = boost::python;

namespace RDKit {

std::string tdtMolSupplierClassDoc =
    "A class which supplies molecules from a TDT file.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules are not constructed until we ask for them:\n\n\
       >>> suppl = TDTMolSupplier('in.smi')\n\
       >>> for mol in suppl:\n\
       ...    mol.GetNumAtoms()\n\
\n\
    2) Lazy evaluation 2:\n\n\
       >>> suppl = TDTMolSupplier('in.smi')\n\
       >>> mol1 = next(suppl)\n\
       >>> mol2 = next(suppl)\n\
       >>> suppl.reset()\n\
       >>> mol3 = next(suppl)\n\n\
       # mol3 and mol1 are the same:\
       >>> MolToSmiles(mol3)==MolToSmiles(mol1)\n\
\n\
    3) Random Access:  all molecules are constructed as soon as we ask for the\n\
       length:\n\n\
       >>> suppl = TDTMolSupplier('in.smi')\n\
       >>> nMols = len(suppl)\n\
       >>> for i in range(nMols):\n\
       ...   suppl[i].GetNumAtoms()\n\
\n\
  Properties in the file are used to set properties on each molecule.\n\
  The properties are accessible using the mol.GetProp(propName) method.\n\
\n";
struct tdtmolsup_wrap {
  static void wrap() {
    python::class_<TDTMolSupplier, boost::noncopyable>(
        "TDTMolSupplier", tdtMolSupplierClassDoc.c_str(), python::init<>())
        .def(python::init<std::string, std::string, int, int, bool>(
            (python::arg("fileName"), python::arg("nameRecord") = "",
             python::arg("confId2D") = -1, python::arg("confId3D") = -1,
             python::arg("sanitize") = true)))
        .def("__enter__", &MolIOEnter<TDTMolSupplier>,
             python::return_internal_reference<>())
        .def("__exit__", &MolIOExit<TDTMolSupplier>)
        .def("__iter__", &MolSupplIter<TDTMolSupplier>,
             python::return_internal_reference<1>())
        .def("__next__", &MolSupplNext<TDTMolSupplier>,
             "Returns the next molecule in the file.  Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("__getitem__", &MolSupplGetItem<TDTMolSupplier>,
             python::return_value_policy<python::manage_new_object>())
        .def("reset", &TDTMolSupplier::reset,
             "Resets our position in the file to the beginning.\n")
        .def("__len__", &TDTMolSupplier::length)
        .def("SetData", &TDTMolSupplier::setData, "Sets the text to be parsed",
             (python::arg("self"), python::arg("data"),
              python::arg("nameRecord") = "", python::arg("confId2D") = -1,
              python::arg("confId3D") = -1, python::arg("sanitize") = true))
        .def("GetItemText", &TDTMolSupplier::getItemText,
             "returns the text for an item",
             (python::arg("self"), python::arg("index")));
  };
};
}  // namespace RDKit

void wrap_tdtsupplier() { RDKit::tdtmolsup_wrap::wrap(); }
