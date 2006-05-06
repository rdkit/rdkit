// $Id: FragCatParams.cpp 4968 2006-02-18 00:27:15Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>
#include <string>

#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FragCatalog/FragCatalogEntry.h>
#include <GraphMol/FragCatalog/FragCatGenerator.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <Catalogs/Catalog.h>
#include <Catalogs/CatalogParams.h>


namespace python = boost::python;
namespace RDKit{
  struct fragparams_wrapper {
    static void wrap() {
      // this exposed to be read only 
      // i.e. one constructed there can be no changes done to the parameter object
      python::class_<FragCatParams>("FragCatParams",
				    python::init<int, int, std::string>())
	.def(python::init<int, int, std::string,double>())
	.def("GetTypeString", &FragCatParams::getTypeStr)
	.def("GetUpperFragLength", &FragCatParams::getUpperFragLength)
	.def("GetLowerFragLength", &FragCatParams::getLowerFragLength)
	.def("GetTolerance", &FragCatParams::getTolerance)
	.def("GetNumFuncGroups", &FragCatParams::getNumFuncGroups)
	.def("GetFuncGroup", (const ROMol* (FragCatParams::*)(int) const)&FragCatParams::getFuncGroup,
	     python::return_value_policy<python::reference_existing_object>())
	.def("Serialize", &FragCatParams::Serialize)
	;
    };
  };
}

void wrap_fragparams() {
  RDKit::fragparams_wrapper::wrap();
}
