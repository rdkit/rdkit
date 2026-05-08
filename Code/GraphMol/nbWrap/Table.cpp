//
//  Copyright (C) 2026 Greg Landrm
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

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
    - GetRow\n\
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
  static void wrap(nb::module_ &m) {
    nb::class_<PeriodicTable>(m, "PeriodicTable")
        .def("GetAtomicWeight",
             nb::overload_cast<unsigned int>(&PeriodicTable::getAtomicWeight,
                                             nb::const_),
             "atomicNumber"_a)
        .def("GetAtomicWeight",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getAtomicWeight, nb::const_),
             "elementSymbol"_a)
        .def("GetAtomicNumber",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getAtomicNumber, nb::const_),
             "elementSymbol"_a)
        .def("GetElementSymbol",
             nb::overload_cast<unsigned int>(&PeriodicTable::getElementSymbol,
                                             nb::const_),
             "atomicNumber"_a)
        .def("GetElementName",
             nb::overload_cast<unsigned int>(&PeriodicTable::getElementName,
                                             nb::const_),
             "atomicNumber"_a)
        .def(
            "GetRow",
            nb::overload_cast<unsigned int>(&PeriodicTable::getRow, nb::const_),
            "atomicNumber"_a)
        .def("GetRow",
             nb::overload_cast<const std::string &>(&PeriodicTable::getRow,
                                                    nb::const_),
             "elementSymbol"_a)
        .def("GetRvdw",
             nb::overload_cast<unsigned int>(&PeriodicTable::getRvdw,
                                             nb::const_),
             "atomicNumber"_a)
        .def("GetRvdw",
             nb::overload_cast<const std::string &>(&PeriodicTable::getRvdw,
                                                    nb::const_),
             "elementSymbol"_a)
        .def("GetRcovalent",
             nb::overload_cast<unsigned int>(&PeriodicTable::getRcovalent,
                                             nb::const_),
             "atomicNumber"_a)
        .def("GetRcovalent",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getRcovalent, nb::const_),
             "elementSymbol"_a)
        .def("GetDefaultValence",
             nb::overload_cast<unsigned int>(&PeriodicTable::getDefaultValence,
                                             nb::const_),
             "atomicNumber"_a)
        .def("GetDefaultValence",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getDefaultValence, nb::const_),
             "elementSymbol"_a)
        .def("GetValenceList",
             nb::overload_cast<unsigned int>(&PeriodicTable::getValenceList,
                                             nb::const_),
             nb::rv_policy::reference_internal, "atomicNumber"_a)
        .def("GetValenceList",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getValenceList, nb::const_),
             nb::rv_policy::reference_internal, "elementSymbol"_a)
        .def("GetNOuterElecs",
             nb::overload_cast<unsigned int>(&PeriodicTable::getNouterElecs,
                                             nb::const_),
             "atomicNumber"_a)
        .def("GetNOuterElecs",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getNouterElecs, nb::const_),
             "elementSymbol"_a)
        .def("GetMostCommonIsotope",
             nb::overload_cast<unsigned int>(
                 &PeriodicTable::getMostCommonIsotope, nb::const_),
             "atomicNumber"_a)
        .def("GetMostCommonIsotope",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getMostCommonIsotope, nb::const_),
             "elementSymbol"_a)
        .def("GetMostCommonIsotopeMass",
             nb::overload_cast<unsigned int>(
                 &PeriodicTable::getMostCommonIsotopeMass, nb::const_),
             "atomicNumber"_a)
        .def("GetMostCommonIsotopeMass",
             nb::overload_cast<const std::string &>(
                 &PeriodicTable::getMostCommonIsotopeMass, nb::const_),
             "elementSymbol"_a)
        .def(
            "GetRb0",
            nb::overload_cast<unsigned int>(&PeriodicTable::getRb0, nb::const_),
            "atomicNumber"_a)
        .def("GetRb0",
             nb::overload_cast<const std::string &>(&PeriodicTable::getRb0,
                                                    nb::const_),
             "elementSymbol"_a)
        .def("GetAbundanceForIsotope",
             nb::overload_cast<unsigned int, unsigned int>(
                 &PeriodicTable::getAbundanceForIsotope, nb::const_),
             "atomicNumber"_a, "isotope"_a)
        .def("GetAbundanceForIsotope",
             nb::overload_cast<const std::string &, unsigned int>(
                 &PeriodicTable::getAbundanceForIsotope, nb::const_),
             "elementSymbol"_a, "isotope"_a)
        .def("GetMassForIsotope",
             nb::overload_cast<unsigned int, unsigned int>(
                 &PeriodicTable::getMassForIsotope, nb::const_),
             "atomicNumber"_a, "isotope"_a)
        .def("GetMassForIsotope",
             nb::overload_cast<const std::string &, unsigned int>(
                 &PeriodicTable::getMassForIsotope, nb::const_),
             "elementSymbol"_a, "isotope"_a)
        .def("GetMaxAtomicNumber", &PeriodicTable::getMaxAtomicNumber)
        .doc() = periodicTableClassDoc.c_str();

    m.def("GetPeriodicTable", PeriodicTable::getTable,
          "Returns the application's PeriodicTable instance.\n\n",
          nb::rv_policy::reference);
  };
};
}  // namespace RDKit
void wrap_table(nb::module_ &m) { RDKit::table_wrapper::wrap(m); }
