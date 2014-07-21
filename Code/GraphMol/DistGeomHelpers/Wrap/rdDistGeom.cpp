// $Id$
//
//  Copyright (C) 2004-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>
#define PY_ARRAY_UNIQUE_SYMBOL rdDistGeom_array_API
#include "numpy/arrayobject.h"
#include <DistGeom/BoundsMatrix.h>

#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace python = boost::python;

namespace RDKit {
  int EmbedMolecule(ROMol &mol, unsigned int maxAttempts,
                    int seed, bool clearConfs,
		    bool useRandomCoords,double boxSizeMult,
                    bool randNegEig, unsigned int numZeroFail,
                    python::dict &coordMap,double forceTol,
                    bool ignoreSmoothingFailures){
    std::map<int,RDGeom::Point3D> pMap;
    python::list ks = coordMap.keys();
    unsigned int nKeys=python::extract<unsigned int>(ks.attr("__len__")());
    for(unsigned int i=0;i<nKeys;++i){
      unsigned int id = python::extract<unsigned int>(ks[i]);
      pMap[id]= python::extract<RDGeom::Point3D>(coordMap[id]);
    }
    std::map<int,RDGeom::Point3D> *pMapPtr=0;
    if(nKeys){
      pMapPtr=&pMap;
    }

    int res = DGeomHelpers::EmbedMolecule(mol, maxAttempts, 
                                          seed, clearConfs,
					  useRandomCoords,boxSizeMult,
					  randNegEig,
                                          numZeroFail,
                                          pMapPtr,forceTol,
                                          ignoreSmoothingFailures);
    return res;
  }

  INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs,
			      unsigned int maxAttempts,
                              int seed, bool clearConfs,
			      bool useRandomCoords,double boxSizeMult,
                              bool randNegEig, unsigned int numZeroFail,
			      double pruneRmsThresh,python::dict &coordMap,
                              double forceTol,
                              bool ignoreSmoothingFailures) {

    std::map<int,RDGeom::Point3D> pMap;
    python::list ks = coordMap.keys();
    unsigned int nKeys=python::extract<unsigned int>(ks.attr("__len__")());
    for(unsigned int i=0;i<nKeys;++i){
      unsigned int id = python::extract<unsigned int>(ks[i]);
      pMap[id]= python::extract<RDGeom::Point3D>(coordMap[id]);
    }
    std::map<int,RDGeom::Point3D> *pMapPtr=0;
    if(nKeys){
      pMapPtr=&pMap;
    }

    INT_VECT res = DGeomHelpers::EmbedMultipleConfs(mol, numConfs, maxAttempts,
                                                    seed, clearConfs,
						    useRandomCoords,boxSizeMult, 
                                                    randNegEig, numZeroFail,
                                                    pruneRmsThresh,pMapPtr,forceTol,
                                                    ignoreSmoothingFailures);

    return res;
  } 

  PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds=true,
                               bool scaleVDW=false) {
    unsigned int nats = mol.getNumAtoms();
    npy_intp dims[2];
    dims[0] = nats;
    dims[1] = nats;

    DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
    DGeomHelpers::initBoundsMat(mat);
    DGeomHelpers::setTopolBounds(mol,mat, set15bounds, scaleVDW);
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    memcpy(static_cast<void *>(res->data),
	   static_cast<void *>(mat->getData()),
	   nats*nats*sizeof(double));
	   
    return PyArray_Return(res);
  }
}

BOOST_PYTHON_MODULE(rdDistGeom) {
  python::scope().attr("__doc__") =
    "Module containing functions to compute atomic coordinates in 3D using distance geometry"
    ;

  rdkit_import_array();

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
    - coordMap : a dictionary mapping atom IDs->coordinates. Use this to \n\
                 require some atoms to have fixed coordinates in the resulting \n\
                 conformation.\n\
    - forceTol : tolerance to be used during the force-field minimization with \n\
                 the distance geometry force field.\n\
    - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing\n\
                 of the bounds matrix fails.\n\
\n\
 RETURNS:\n\n\
    ID of the new conformation added to the molecule \n\
\n";
  python::def("EmbedMolecule", RDKit::EmbedMolecule,
              (python::arg("mol"), python::arg("maxAttempts")=0,
               python::arg("randomSeed")=-1, python::arg("clearConfs")=true,
               python::arg("useRandomCoords")=false,
	       python::arg("boxSizeMult")=2.0,
               python::arg("randNegEig")=true, python::arg("numZeroFail")=1,
               python::arg("coordMap")=python::dict(),
               python::arg("forceTol")=1e-3,
               python::arg("ignoreSmoothingFailures")=false),
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
    - coordMap : a dictionary mapping atom IDs->coordinates. Use this to \n\
                 require some atoms to have fixed coordinates in the resulting \n\
                 conformation.\n\
    - forceTol : tolerance to be used during the force-field minimization with \n\
                 the distance geometry force field.\n\
    - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing\n\
                 of the bounds matrix fails.\n\
 RETURNS:\n\n\
    List of new conformation IDs \n\
\n";
  python::def("EmbedMultipleConfs", RDKit::EmbedMultipleConfs,
              (python::arg("mol"), python::arg("numConfs")=10, 
               python::arg("maxAttempts")=0,
               python::arg("randomSeed")=-1, python::arg("clearConfs")=true,
               python::arg("useRandomCoords")=false,
	       python::arg("boxSizeMult")=2.0,
               python::arg("randNegEig")=true, python::arg("numZeroFail")=1,
	       python::arg("pruneRmsThresh")=-1.0,
               python::arg("coordMap")=python::dict(),
               python::arg("forceTol")=1e-3,
               python::arg("ignoreSmoothingFailures")=false),
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
