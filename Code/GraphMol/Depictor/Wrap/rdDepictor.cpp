// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>

#include <GraphMol/Depictor/RDDepictor.h>

using namespace RDDepict;


namespace python = boost::python;
 
namespace RDDepict {
  
  unsigned int Compute2DCoords(RDKit::ROMol &mol, bool canonOrient,
			       bool clearConfs, python::dict &coordMap){
    RDKit::INT_POINT2D_MAP cMap;
    cMap.clear();
    python::list ks = coordMap.keys();
    for(unsigned int i=0;
	i<python::extract<unsigned int>(ks.attr("__len__")());
	i++){
      unsigned int id = python::extract<unsigned int>(ks[i]);
      if(id>=mol.getNumAtoms()){
	throw_value_error("atom index out of range");
      }
      cMap[id] = python::extract<RDGeom::Point2D>(coordMap[id]);
    }
    return RDDepict::compute2DCoords(mol,&cMap,canonOrient, clearConfs);
  }
}

BOOST_PYTHON_MODULE(rdDepictor)
{
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  python::scope().attr("__doc__") =
    "Module containing the functionality to compute 2D coordinates for a molecule"
    ;
  
  std::string docString;
  
  docString = "Compute 2D coordinates for a molecule. \n\
  The resulting coordinates are stored on each atom of the molecule \n\n\
  ARGUMENTS: \n\n\
     mol - the molecule of interest\n\
     canonOrient - orient the molecule in a canonical way\n\
     clearConf - if true, all existing conformations on the molecule will be cleared\n\
     coordMap - a dictionary mapping atom Ids -> Point2D objects with starting coordinates\n\
                for atoms that should have forced positions.\n\n\
  RETURNS: \n\n\
     ID of the conformation added to the molecule\n";
  python::def("Compute2DCoords", RDDepict::Compute2DCoords,
	      (python::arg("mol"),
	       python::arg("canonOrient")=false,
	       python::arg("clearConfs")=true,
	       python::arg("coordMap")=python::dict()),
	      docString.c_str());
}
