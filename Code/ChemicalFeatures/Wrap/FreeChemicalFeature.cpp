// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include<RDGeneral/Invariant.h>
#include <RDBoost/PySequenceHolder.h>
#include <ChemicalFeatures/FreeChemicalFeature.h>

namespace ChemicalFeatures {

  // allows BitVects to be pickled
  struct chemfeat_pickle_suite : python::pickle_suite
  {
    static python::tuple
    getinitargs(const FreeChemicalFeature& self)
    {
      return python::make_tuple(self.toString());
    };
  };

  
  std::string featClassDoc="Class to represent a free chemical  feature.\n\
    These chemical features are not derived from a molecule, though they can be matched \n\
    to features from a molecule\n";
  struct freefeat_wrapper {
    static void wrap() {
      python::class_<FreeChemicalFeature>("FreeChemicalFeature", featClassDoc.c_str(),
					  python::init<const std::string &>()) //<>("Default Constructor"))
        .def(python::init<>("Default Constructor"))
        .def(python::init<std::string, std::string, const RDGeom::Point3D &>
             (python::args("family", "type", "loc"),
              "Constructor with family, type and location specified"))
        .def(python::init<std::string, const RDGeom::Point3D &>
             (python::args("family", "loc"), 
              "constructor with family and location specified, empty type"))
        .def("SetFamily", &FreeChemicalFeature::setFamily,
             "Set the family of the feature")
        .def("SetType", &FreeChemicalFeature::setType,
             "Set the sepcific type for the feature")
        .def("GetFamily", &FreeChemicalFeature::getFamily,
             "Get the family of the feature",
             python::return_value_policy<python::copy_const_reference>())
        .def("GetType", &FreeChemicalFeature::getType,
             "Get the sepcific type for the feature",
             python::return_value_policy<python::copy_const_reference>())
        .def("SetPos", &FreeChemicalFeature::setPos,
             "Set the feature position")
        .def("GetPos", &FreeChemicalFeature::getPos,
             "Get the psotion for the feature")
	.def_pickle(chemfeat_pickle_suite())
        ;
      //python::def("FreeChemicalFeature", ChemicalFeatures::makeFreeChemFeat,
      //            (python::arg("family")="", python::arg("type")="",
      //             python::arg("pos")=python::list()),
      //            "This a kind of faking a constructor for FreeChemicalFeature");


    };
  }; 
}

void wrap_freefeat() {
  ChemicalFeatures::freefeat_wrapper::wrap();
}

