// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <DistGeom/BoundsMatrix.h>
#include "rdDistGeom.h"

#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace python = boost::python;

namespace RDKit {
  int EmbedMolecule(ROMol &mol, unsigned int maxAttempts=30,
                    int seed=-1, bool clearConfs=true,
                    bool randNegEig=true, unsigned int numZeroFail=1){
   
    int res = DGeomHelpers::EmbedMolecule(mol, maxAttempts, 
                                          seed, clearConfs, randNegEig,
                                          numZeroFail);
    return res;
  }

  INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs=10, unsigned int maxAttempts=30, 
                              int seed=-1, bool clearConfs=true, 
                              bool randNegEig=true, unsigned int numZeroFail=1) {
    INT_VECT res = DGeomHelpers::EmbedMultipleConfs(mol, numConfs, maxAttempts,
                                                    seed, clearConfs,
                                                    randNegEig, numZeroFail);
    return res;
  } 

  PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds=true,
                               bool scaleVDW=false) {
    int dims[2];
    unsigned int nats = mol.getNumAtoms();
    dims[0] = nats;
    dims[1] = nats;

    DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
    DGeomHelpers::initBoundsMat(mat);
    DGeomHelpers::setTopolBounds(mol,mat, set15bounds, scaleVDW);
    PyArrayObject *res = (PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE);
    memcpy(static_cast<void *>(res->data),
	   static_cast<void *>(mat->getData()),
	   nats*nats*sizeof(double));
	   
    return (PyObject *) res;
  }

}

BOOST_PYTHON_MODULE(rdDistGeom) {
  python::scope().attr("__doc__") =
    "Module containing functions to compute starting atomic coordinates in 3D using distance geometry"
    ;

  import_array();

  //RegisterListConverter<RDKit::Atom*>();

  std::string docString = "Use distance geometry to obtain intial coordinates for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - maxAttempts : the maximum number of attempts to try embedding \n\
    - randomSeed : provide a seed for the random number generator so that the same coordinates can be obtained \n\
                   for a molecule on multiple runs. The default (-1) produces a random embedding\n\
    - clearConfs : clear all existing conformations on the molecule\n\
    - randNegEig : If the embedding yields a negative eigen value, pick coordinates that correspond \n\
                   to this component at random \n\
    - numZeroFail : fail embedding is we have this more zero eigen values \n\
                  \n\
 RETURNS:\n\n\
    ID of the new conformation added to the molecule \n\
\n";
  python::def("EmbedMolecule", RDKit::EmbedMolecule,
              (python::arg("mol"), python::arg("maxAttempts")=30,
               python::arg("randomSeed")=-1, python::arg("clearConfs")=true,
               python::arg("randNegEig")=true, python::arg("numZeroFail")=1),
              docString.c_str());

  docString = "Use distance geometry to obtain intial coordinates for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - numConfs : the number of conformers to generate \n\
    - maxAttempts : the maximum number of attempts to try embedding \n\
    - randomSeed : provide a seed for the random number generator so that the same coordinates can be obtained \n\
                   for a molecule on multiple runs. The default (-1) produces a random embedding\n\
    - clearConfs : clear all existing conformations on the molecule\n\
    - randNegEig : If the embedding yields a negative eigen value, pick coordinates that correspond \n\
                   to this component at random \n\
    - numZeroFail : fail embedding is we have this more zero eigen values \n\
                  \n\
 RETURNS:\n\n\
    List of new conformation IDs \n\
\n";
  python::def("EmbedMultipleConfs", RDKit::EmbedMultipleConfs,
              (python::arg("mol"), python::arg("numConfs")=10, 
               python::arg("maxAttempts")=10,
               python::arg("randomSeed")=-1, python::arg("clearConfs")=true,
               python::arg("randNegEig")=true, python::arg("numZeroFail")=1),
              docString.c_str());

  docString = "Returns the distance bounds matrix for a molecule\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - set15bounds : set bounds 15 atom distance bounds based on topology (otherwise stop at 14s)\n\
    - scaleVDW : scale down the sum of VDW radii when setting the lower bounds for \n\
                 atoms less 5 bonds aparts \n\
 RETURNS:\n\n\
    the bounds matrix as a Numeric array with lower bounds in the lower triangle\n\
    and upper bounds in the upper triangle\n\
\n";
  python::def("GetMoleculeBoundsMatrix", RDKit::getMolBoundsMatrix,
              (python::arg("mol"), python::arg("set15bounds")=true,
               python::arg("scaleVDW")=false),
              docString.c_str());


}
