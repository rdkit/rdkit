// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
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

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
namespace RDKit {

PeriodicTable *GetTable() { return PeriodicTable::getTable(); }

std::string periodicTableClassDoc =
    "A class which stores information from the Periodic Table.\n\
\n\
  It is not possible to create a PeriodicTable object directly from Python,\n\
  use GetPeriodicTable() to get the global table.\n\
\n\
  The PeriodicTable object can be queried for a variety of properties:\n\
\n\
    - GetAtomicWeight\n\
\n\
    - GetAtomicNumber\n\
\n\
    - GetElementSymbol\n\
\n\
    - GetElementName\n\
\n\
    - GetRvdw (van der Waals radius)\n\
\n\
    - GetRCovalent (covalent radius)\n\
\n\
    - GetDefaultValence\n\
\n\
    - GetValenceList\n\
\n\
    - GetNOuterElecs (number of valence electrons)\n\
\n\
    - GetMostCommonIsotope\n\
\n\
    - GetMostCommonIsotopeMass\n\
\n\
    - GetRb0\n\
\n\
    - GetAbundanceForIsotope\n\
\n\
    - GetMassForIsotope\n\
\n\
  When it makes sense, these can be queried using either an atomic number (integer)\n\
  or an atomic symbol (string)\n\
\n";

struct table_wrapper {
  static void wrap() {
    python::class_<PeriodicTable>(
        "PeriodicTable", periodicTableClassDoc.c_str(), python::no_init)
        .def("GetAtomicWeight",
             (double(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getAtomicWeight,
             python::args("self", "elementSymbol"))
        .def("GetAtomicWeight",
             (double(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getAtomicWeight,
             python::args("self", "elementSymbol"))
        .def("GetAtomicNumber",
             (int(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getAtomicNumber,
             python::args("self", "elementSymbol"))
        .def("GetElementSymbol",
             (std::string(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getElementSymbol,
             python::args("self", "atomicNumber"))
        .def("GetElementName",
             (std::string(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getElementName,
             python::args("self", "atomicNumber"))
        .def("GetRvdw",
             (double(PeriodicTable::*)(UINT) const) & PeriodicTable::getRvdw,
             python::args("self", "elementSymbol"))
        .def("GetRvdw",
             (double(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getRvdw,
             python::args("self", "elementSymbol"))
        .def("GetRcovalent",
             (double(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getRcovalent,
             python::args("self", "elementSymbol"))
        .def("GetRcovalent",
             (double(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getRcovalent,
             python::args("self", "elementSymbol"))
        .def("GetDefaultValence",
             (int(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getDefaultValence,
             python::args("self", "elementSymbol"))
        .def("GetDefaultValence",
             (int(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getDefaultValence,
             python::args("self", "elementSymbol"))
        .def("GetValenceList",
             (const INT_VECT &(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getValenceList,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self", "elementSymbol"))
        .def("GetValenceList",
             (const INT_VECT &(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getValenceList,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self", "elementSymbol"))
        .def(
            "GetNOuterElecs",
            (int(PeriodicTable::*)(UINT) const) & PeriodicTable::getNouterElecs,
            python::args("self", "elementSymbol"))
        .def("GetNOuterElecs",
             (int(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getNouterElecs,
             python::args("self", "elementSymbol"))
        .def("GetMostCommonIsotope",
             (int(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getMostCommonIsotope,
             python::args("self", "elementSymbol"))
        .def("GetMostCommonIsotope",
             (int(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getMostCommonIsotope,
             python::args("self", "elementSymbol"))
        .def("GetMostCommonIsotopeMass",
             (double(PeriodicTable::*)(UINT) const) &
                 PeriodicTable::getMostCommonIsotopeMass,
             python::args("self", "elementSymbol"))
        .def("GetMostCommonIsotopeMass",
             (double(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getMostCommonIsotopeMass,
             python::args("self", "elementSymbol"))
        .def("GetRb0",
             (double(PeriodicTable::*)(UINT) const) & PeriodicTable::getRb0,
             python::args("self", "elementSymbol"))
        .def("GetRb0",
             (double(PeriodicTable::*)(const std::string &) const) &
                 PeriodicTable::getRb0,
             python::args("self", "elementSymbol"))
        .def("GetAbundanceForIsotope",
             (double(PeriodicTable::*)(UINT, UINT) const) &
                 PeriodicTable::getAbundanceForIsotope,
             python::args("self", "elementSymbol", "isotope"))
        .def("GetAbundanceForIsotope",
             (double(PeriodicTable::*)(const std::string &, UINT) const) &
                 PeriodicTable::getAbundanceForIsotope,
             python::args("self", "elementSymbol", "isotope"))
        .def("GetMassForIsotope",
             (double(PeriodicTable::*)(UINT, UINT) const) &
                 PeriodicTable::getMassForIsotope,
             python::args("self", "elementSymbol", "isotope"))
        .def("GetMassForIsotope",
             (double(PeriodicTable::*)(const std::string &, UINT) const) &
                 PeriodicTable::getMassForIsotope,
             python::args("self", "elementSymbol", "isotope"))
        .def("GetMaxAtomicNumber", &PeriodicTable::getMaxAtomicNumber,
             python::args("self"));

    python::def(
        "GetPeriodicTable", GetTable,
        "Returns the application's PeriodicTable instance.\n\n",
        python::return_value_policy<python::reference_existing_object>());
  };
};
}  // namespace RDKit
void wrap_table() { RDKit::table_wrapper::wrap(); }
