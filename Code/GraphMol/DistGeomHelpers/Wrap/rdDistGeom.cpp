// $Id$
//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
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
		    bool useRandomCoords=false,double boxSizeMult=2.0,
                    bool randNegEig=true, unsigned int numZeroFail=1){
    int res = DGeomHelpers::EmbedMolecule(mol, maxAttempts, 
                                          seed, clearConfs,
					  useRandomCoords,boxSizeMult,
					  randNegEig,
                                          numZeroFail);
    return res;
  }

  INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs=10,
			      unsigned int maxAttempts=30, 
                              int seed=-1, bool clearConfs=true, 
			      bool useRandomCoords=false,double boxSizeMult=2.0,
                              bool randNegEig=true, unsigned int numZeroFail=1,
			      double pruneRmsThresh=-1.0) {
    INT_VECT res = DGeomHelpers::EmbedMultipleConfs(mol, numConfs, maxAttempts,
                                                    seed, clearConfs,
						    useRandomCoords,boxSizeMult,
                                                    randNegEig, numZeroFail, 1e-3,
						    5.0, pruneRmsThresh);
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
    "Module containing functions to compute atomic coordinates in 3D using distance geometry"
    ;

  import_array();

  //RegisterListConverter<RDKit::Atom*>();

  std::string docString = "Use distance geometry to obtain intial \n\
 coordinates for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - maxAttempts : the maximum number of attempts to try embedding \n\
    - randomSeed : provide a seed for the random number generator \n\
                   so that the same coordinates can be obtained \n\
                   for a molecule on multiple runs. The default \n\
                   (-1) uses a random seed \n\
    - clearConfs : clear all existing conformations on the molecule\n\
    - useRandomCoords : Start the embedding from random coordinates instead of\n\
                        using eigenvalues of the distance matrix.\n\
    - boxSizeMult    Determines the size of the box that is used for\n\
                     random coordinates. If this is a positive number, the \n\
                     side length will equal the largest element of the distance\n\
                     matrix times boxSizeMult. If this is a negative number,\n\
                     the side length will equal -boxSizeMult (i.e. independent\n\
                     of the elements of the distance matrix).\n\
    - randNegEig : If the embedding yields a negative eigenvalue, \n\
                   pick coordinates that correspond \n\
                   to this component at random \n\
    - numZeroFail : fail embedding is we have this more zero eigenvalues \n\
\n\
 RETURNS:\n\n\
    ID of the new conformation added to the molecule \n\
\n";
  python::def("EmbedMolecule", RDKit::EmbedMolecule,
              (python::arg("mol"), python::arg("maxAttempts")=30,
               python::arg("randomSeed")=-1, python::arg("clearConfs")=true,
               python::arg("useRandomCoords")=false,
	       python::arg("boxSizeMult")=2.0,
               python::arg("randNegEig")=true, python::arg("numZeroFail")=1),
              docString.c_str());

  docString = "Use distance geometry to obtain multiple sets of \n\
 coordinates for a molecule\n\
 \n\
 ARGUMENTS:\n\n\
  - mol : the molecule of interest\n\
  - numConfs : the number of conformers to generate \n\
  - maxAttempts : the maximum number of attempts to try embedding \n\
  - randomSeed : provide a seed for the random number generator \n\
                 so that the same coordinates can be obtained \n\
                 for a molecule on multiple runs. The default \n\
                 (-1) uses a random seed \n\
  - clearConfs : clear all existing conformations on the molecule\n\
  - useRandomCoords : Start the embedding from random coordinates instead of\n\
                      using eigenvalues of the distance matrix.\n\
  - boxSizeMult    Determines the size of the box that is used for\n\
                   random coordinates. If this is a positive number, the \n\
                   side length will equal the largest element of the distance\n\
                   matrix times boxSizeMult. If this is a negative number,\n\
                   the side length will equal -boxSizeMult (i.e. independent\n\
                   of the elements of the distance matrix).\n\
  - randNegEig : If the embedding yields a negative eigenvalue, \n\
                 pick coordinates that correspond \n\
                 to this component at random \n\
  - numZeroFail : fail embedding is we have this more zero eigenvalues \n\
  - pruneRmsThresh : Retain only the conformations out of 'numConfs' \n\
                    after embedding that are at least \n\
                    this far apart from each other. \n\
          RMSD is computed on the heavy atoms. \n\
          Pruning is greedy; i.e. the first embedded conformation\n\
          is retained and from then on only those that are at\n\
          least pruneRmsThresh away from all retained conformations\n\
          are kept. The pruning is done after embedding and \n\
          bounds violation minimization. No pruning by default.\n\
 RETURNS:\n\n\
    List of new conformation IDs \n\
\n";
  python::def("EmbedMultipleConfs", RDKit::EmbedMultipleConfs,
              (python::arg("mol"), python::arg("numConfs")=10, 
               python::arg("maxAttempts")=10,
               python::arg("randomSeed")=-1, python::arg("clearConfs")=true,
               python::arg("useRandomCoords")=false,
	       python::arg("boxSizeMult")=2.0,
               python::arg("randNegEig")=true, python::arg("numZeroFail")=1,
	       python::arg("pruneRmsThresh")=-1.0),
              docString.c_str());

  docString = "Returns the distance bounds matrix for a molecule\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - set15bounds : set bounds for 1-5 atom distances based on \n\
                    topology (otherwise stop at 1-4s)\n\
    - scaleVDW : scale down the sum of VDW radii when setting the \n\
                 lower bounds for atoms less than 5 bonds apart \n\
 RETURNS:\n\n\
    the bounds matrix as a Numeric array with lower bounds in \n\
    the lower triangle and upper bounds in the upper triangle\n\
\n";
  python::def("GetMoleculeBoundsMatrix", RDKit::getMolBoundsMatrix,
              (python::arg("mol"), python::arg("set15bounds")=true,
               python::arg("scaleVDW")=false),
              docString.c_str());


}
