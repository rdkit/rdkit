// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <boost/python.hpp>

#define PY_ARRAY_UNIQUE_SYMBOL Depictor_array_API
#include "Numeric/arrayobject.h"
#include <RDBoost/Wrap.h>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Depictor/EmbeddedFrag.h>

using namespace RDDepict;


namespace python = boost::python;
 
namespace RDDepict {
  
  unsigned int Compute2DCoords(RDKit::ROMol &mol, bool canonOrient,
			       bool clearConfs, python::dict &coordMap,
                               unsigned int nFlipsPerSample=3, unsigned int nSamples=100,
                               int sampleSeed=100, bool permuteDeg4Nodes=false){
    RDGeom::INT_POINT2D_MAP cMap;
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
    return RDDepict::compute2DCoords(mol,&cMap,canonOrient, clearConfs, nFlipsPerSample, 
                                     nSamples, sampleSeed, permuteDeg4Nodes);
  }
  
  unsigned int Compute2DCoordsMimicDistmat(RDKit::ROMol &mol, python::object distMat, 
                                           bool canonOrient, bool clearConfs, double weightDistMat,
                                           unsigned int nFlipsPerSample, unsigned int nSamples,
                                           int sampleSeed, bool permuteDeg4Nodes) {
    PyObject *distMatPtr = distMat.ptr();
    if(!PyArray_Check(distMatPtr)){
      throw_value_error("Argument isn't an array");
    }

    PyArrayObject *dmatrix = reinterpret_cast<PyArrayObject *>(distMatPtr);
    unsigned int nitems = dmatrix->dimensions[0];
    unsigned int na = mol.getNumAtoms();

    if (nitems != na*(na-1)/2) {
      throw_value_error("The array size does not match the number of atoms in the molecule");
    }
    double *inData = reinterpret_cast<double *>(dmatrix->data);
    double *cData = new double[nitems];
    
    memcpy(static_cast<void *>(cData), 
           static_cast<const void *>(inData),
           nitems*sizeof(double));
    
    DOUBLE_SMART_PTR dmat(cData);
    return RDDepict::compute2DCoordsMimicDistMat(mol, &dmat, canonOrient, clearConfs, weightDistMat,
                                                 nFlipsPerSample, nSamples, sampleSeed, permuteDeg4Nodes);
  }
    
}

BOOST_PYTHON_MODULE(rdDepictor)
{
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  python::scope().attr("__doc__") =
    "Module containing the functionality to compute 2D coordinates for a molecule"
    ;

  import_array();

  std::string docString;
  
  docString = "Compute 2D coordinates for a molecule. \n\
  The resulting coordinates are stored on each atom of the molecule \n\n\
  ARGUMENTS: \n\n\
     mol - the molecule of interest\n\
     canonOrient - orient the molecule in a canonical way\n\
     clearConf - if true, all existing conformations on the molecule will be cleared\n\
     coordMap - a dictionary mapping atom Ids -> Point2D objects with starting coordinates\n\
                for atoms that should have forced positions.\n\
     nFlipsPerSample - number of rotatable bonds that are flipped at random at a time.\n\
     nSample - Number of random samplings of rotatable bond.\n\
     sampleSeed - see for the random sampling process.\n\
     permuteDeg4Nodes - allow permutation of bonds at a degree node node during the sampling process \n\n\
  RETURNS: \n\n\
     ID of the conformation added to the molecule\n";
  python::def("Compute2DCoords", RDDepict::Compute2DCoords,
	      (python::arg("mol"),
	       python::arg("canonOrient")=false,
	       python::arg("clearConfs")=true,
	       python::arg("coordMap")=python::dict(),
               python::arg("nFlipsPerSample")=0,
               python::arg("nSample")=0,
               python::arg("sampleSeed")=0,
               python::arg("permuteDeg4Nodes")=false),
	      docString.c_str());

  docString = "Compute 2D coordinates such that inter-atom distance mimic those in a distance matrix. \n\
  The resulting coordinates are stored on each atom of the molecule \n\n\
  ARGUMENTS: \n\n\
     mol - the molecule of interest\n\
     distMat - distance matrix that we want the 2D structure to mimic\n\
     canonOrient - orient the molecule in a canonical way\n\
     clearConf - if true, all existing conformations on the molecule will be cleared\n\
     weightDistMat - weight assigned in the cost function to mimicing the distance matrix.\n\
                     This is must be between (0.0,1.0). (1.0-weightDistMat) is then the\n\
                     weght assined to improving the density of the 2D structure \n\
                     i.e. try to make it spread out\n\
     nFlipsPerSample - number of rotatable bonds that are flipped at random at a time.\n\
     nSample - Number of random samplings of rotatable bond.\n\
     sampleSeed - see for the random sampling process.\n\
     permuteDeg4Nodes - allow permutation of bonds at a degree node node during the sampling process \n\n\
  RETURNS: \n\n\
     ID of the conformation added to the molecule\n";
  python::def("Compute2DCoordsMimicDistmat", RDDepict::Compute2DCoordsMimicDistmat,
	      (python::arg("mol"),
               python::arg("distMat"),
	       python::arg("canonOrient")=false,
	       python::arg("clearConfs")=true,
               python::arg("weightDistMat")=0.5,
               python::arg("nFlipsPerSample")=3,
               python::arg("nSample")=100,
               python::arg("sampleSeed")=100,
               python::arg("permuteDeg4Nodes")=true),
	      docString.c_str());
}
