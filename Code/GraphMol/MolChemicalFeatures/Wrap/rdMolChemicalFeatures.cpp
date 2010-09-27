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
#include <RDBoost/Wrap.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <boost/python/register_ptr_to_python.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <GraphMol/MolChemicalFeatures/FeatureParser.h>
namespace python = boost::python;
using namespace RDKit;

void wrap_MolChemicalFeat();
void wrap_factory();
void wrap_ChemicalFeatureUtils();

namespace RDKit {
  MolChemicalFeatureFactory *buildFeatFactory(std::string fileName) {
    std::ifstream inStream(fileName.c_str());
    if(!inStream.is_open()){
      std::string errorstring="File: " + fileName + " could not be opened.";
      PyErr_SetString(PyExc_IOError, errorstring.c_str());
      python::throw_error_already_set();
    }
    std::istream &instrm = static_cast<std::istream &>(inStream);
    return buildFeatureFactory(instrm);
  }

  MolChemicalFeatureFactory *buildFeatFactoryFromString(std::string fdefString) {
    std::istringstream inStream(fdefString);
    std::istream &instrm = static_cast<std::istream &>(inStream);
    return buildFeatureFactory(instrm);
  }
}

void translate_FeatureFileParse_error(RDKit::FeatureFileParseException const&e){
  std::stringstream err;
  err << "Error parsing feature file at line " << e.lineNo() << ":" << std::endl;
  err << e.message() << std::endl;
  PyErr_SetString(PyExc_ValueError, err.str().c_str());
  python::throw_error_already_set();
}


BOOST_PYTHON_MODULE(rdMolChemicalFeatures)
{

  python::scope().attr("__doc__")=
    "Module containing from chemical feature and functions to generate the"
    ;
  python::register_exception_translator<RDKit::FeatureFileParseException>(&translate_FeatureFileParse_error);
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
    
  python::def("BuildFeatureFactory", RDKit::buildFeatFactory,
              "Construct a feature factory given a feature definition in a file",
              python::return_value_policy<python::manage_new_object>());
  python::def("BuildFeatureFactoryFromString", RDKit::buildFeatFactoryFromString,
              "Construct a feature factory given a feature definition block",
              python::return_value_policy<python::manage_new_object>());

  wrap_MolChemicalFeat();
  wrap_factory();
  wrap_ChemicalFeatureUtils();
}
