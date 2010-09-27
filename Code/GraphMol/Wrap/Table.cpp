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
#include <boost/python.hpp>
#include <string>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
namespace RDKit{

  PeriodicTable *GetTable(){
    return PeriodicTable::getTable();
  }

  std::string periodicTableClassDoc="A class which stores information from the Periodic Table.\n\
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
  When it makes sense, these can be queried using either an atomic number (integer)\n\
  or an atomic symbol (string)\n\
\n";


struct table_wrapper {
  static void wrap(){
    python::class_<PeriodicTable>("PeriodicTable",periodicTableClassDoc.c_str(),
				  python::no_init)
      .def("GetAtomicWeight",(double (PeriodicTable::*)(UINT) const)&PeriodicTable::getAtomicWeight)
      .def("GetAtomicWeight",(double (PeriodicTable::*)(const std::string &) const)&PeriodicTable::getAtomicWeight)
      .def("GetAtomicNumber",(int (PeriodicTable::*)(const std::string &) const)&PeriodicTable::getAtomicNumber)
      .def("GetElementSymbol",(std::string (PeriodicTable::*)(UINT) const)&PeriodicTable::getElementSymbol)
      .def("GetRvdw",(double (PeriodicTable::*)(UINT) const)&PeriodicTable::getRvdw)
      .def("GetRvdw",(double (PeriodicTable::*)(const std::string &) const)&PeriodicTable::getRvdw)
      .def("GetRcovalent",(double (PeriodicTable::*)(UINT) const)&PeriodicTable::getRcovalent)
      .def("GetRcovalent",(double (PeriodicTable::*)(const std::string &) const)&PeriodicTable::getRcovalent)
      .def("GetDefaultValence",(int (PeriodicTable::*)(UINT) const)&PeriodicTable::getDefaultValence)
      .def("GetDefaultValence",(int (PeriodicTable::*)(const std::string &) const)&PeriodicTable::getDefaultValence)
      .def("GetValenceList",(const INT_VECT &(PeriodicTable::*)(UINT) const)&PeriodicTable::getValenceList,
	   python::return_value_policy<python::copy_const_reference>())
      .def("GetValenceList",(const INT_VECT &(PeriodicTable::*)(const std::string &) const)&PeriodicTable::getValenceList,
	   python::return_value_policy<python::copy_const_reference>())
      .def("GetNOuterElecs",(int (PeriodicTable::*)(UINT) const)&PeriodicTable::getNouterElecs)
      .def("GetNOuterElecs",(int (PeriodicTable::*)(const std::string &) const)&PeriodicTable::getNouterElecs)
      ;

    python::def("GetPeriodicTable",GetTable,
		"Returns the application's PeriodicTable instance.\n\n",
		python::return_value_policy<python::reference_existing_object>());


  };
};
} // end of namespace
void wrap_table() {
  RDKit::table_wrapper::wrap();
}
