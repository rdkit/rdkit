//
//  Copyright (C) 2004-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#define PY_ARRAY_UNIQUE_SYMBOL rdDistGeom_array_API
#include <RDBoost/import_array.h>
#include "numpy/arrayobject.h"
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>

#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace python = boost::python;

namespace RDKit {
int EmbedMolecule(ROMol &mol, unsigned int maxAttempts, int seed,
                  bool clearConfs, bool useRandomCoords, double boxSizeMult,
                  bool randNegEig, unsigned int numZeroFail,
                  python::dict &coordMap, double forceTol,
                  bool ignoreSmoothingFailures, bool enforceChirality,
                  bool useExpTorsionAnglePrefs, bool useBasicKnowledge,
                  bool printExpTorsionAngles, bool useSmallRingTorsions,
                  bool useMacrocycleTorsions, unsigned int ETversion) {
  std::map<int, RDGeom::Point3D> pMap;
  python::list ks = coordMap.keys();
  unsigned int nKeys = python::extract<unsigned int>(ks.attr("__len__")());
  for (unsigned int i = 0; i < nKeys; ++i) {
    unsigned int id = python::extract<unsigned int>(ks[i]);
    pMap[id] = python::extract<RDGeom::Point3D>(coordMap[id]);
  }
  std::map<int, RDGeom::Point3D> *pMapPtr = nullptr;
  if (nKeys) {
    pMapPtr = &pMap;
  }

  bool verbose = printExpTorsionAngles;
  int numThreads = 1;
  double pruneRmsThresh = -1.;
  const double basinThresh = DGeomHelpers::EmbedParameters().basinThresh;
  bool onlyHeavyAtomsForRMS = false;
  DGeomHelpers::EmbedParameters params(
      maxAttempts, numThreads, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, pMapPtr, forceTol, ignoreSmoothingFailures,
      enforceChirality, useExpTorsionAnglePrefs, useBasicKnowledge, verbose,
      basinThresh, pruneRmsThresh, onlyHeavyAtomsForRMS, ETversion, nullptr,
      true, useSmallRingTorsions, useMacrocycleTorsions);

  int res;
  {
    NOGIL gil;
    res = DGeomHelpers::EmbedMolecule(mol, params);
  }
  return res;
}

int EmbedMolecule2(ROMol &mol, DGeomHelpers::EmbedParameters &params) {
  int res;
  {
    NOGIL gil;
    res = DGeomHelpers::EmbedMolecule(mol, params);
  }
  return res;
}

INT_VECT EmbedMultipleConfs(
    ROMol &mol, unsigned int numConfs, unsigned int maxAttempts, int seed,
    bool clearConfs, bool useRandomCoords, double boxSizeMult, bool randNegEig,
    unsigned int numZeroFail, double pruneRmsThresh, python::dict &coordMap,
    double forceTol, bool ignoreSmoothingFailures, bool enforceChirality,
    int numThreads, bool useExpTorsionAnglePrefs, bool useBasicKnowledge,
    bool printExpTorsionAngles, bool useSmallRingTorsions,
    bool useMacrocycleTorsions, unsigned int ETversion) {
  std::map<int, RDGeom::Point3D> pMap;
  python::list ks = coordMap.keys();
  unsigned int nKeys = python::extract<unsigned int>(ks.attr("__len__")());
  for (unsigned int i = 0; i < nKeys; ++i) {
    unsigned int id = python::extract<unsigned int>(ks[i]);
    pMap[id] = python::extract<RDGeom::Point3D>(coordMap[id]);
  }
  std::map<int, RDGeom::Point3D> *pMapPtr = nullptr;
  if (nKeys) {
    pMapPtr = &pMap;
  }
  bool verbose = printExpTorsionAngles;
  const double basinThresh = DGeomHelpers::EmbedParameters().basinThresh;
  bool onlyHeavyAtomsForRMS = false;
  DGeomHelpers::EmbedParameters params(
      maxAttempts, numThreads, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, pMapPtr, forceTol, ignoreSmoothingFailures,
      enforceChirality, useExpTorsionAnglePrefs, useBasicKnowledge, verbose,
      basinThresh, pruneRmsThresh, onlyHeavyAtomsForRMS, ETversion, nullptr,
      true, useSmallRingTorsions, useMacrocycleTorsions);

  INT_VECT res;
  {
    NOGIL gil;
    DGeomHelpers::EmbedMultipleConfs(mol, res, numConfs, params);
  }
  return res;
}

INT_VECT EmbedMultipleConfs2(ROMol &mol, unsigned int numConfs,
                             DGeomHelpers::EmbedParameters &params) {
  INT_VECT res;
  {
    NOGIL gil;
    DGeomHelpers::EmbedMultipleConfs(mol, res, numConfs, params);
  }
  return res;
}

PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds = true,
                             bool scaleVDW = false,
                             bool doTriangleSmoothing = true,
                             bool useMacrocycle14config = false) {
  unsigned int nats = mol.getNumAtoms();
  npy_intp dims[2];
  dims[0] = nats;
  dims[1] = nats;

  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
  DGeomHelpers::initBoundsMat(mat);
  DGeomHelpers::setTopolBounds(mol, mat, set15bounds, scaleVDW,
                               useMacrocycle14config);
  if (doTriangleSmoothing) {
    DistGeom::triangleSmoothBounds(mat);
  }
  auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  memcpy(static_cast<void *>(PyArray_DATA(res)),
         static_cast<void *>(mat->getData()), nats * nats * sizeof(double));

  return PyArray_Return(res);
}
DGeomHelpers::EmbedParameters *getETKDG() {  // ET version 1
  return new DGeomHelpers::EmbedParameters(DGeomHelpers::ETKDG);
}
DGeomHelpers::EmbedParameters *getETKDGv2() {  // ET version 2
  return new DGeomHelpers::EmbedParameters(DGeomHelpers::ETKDGv2);
}
DGeomHelpers::EmbedParameters *
getETKDGv3() {  //! Parameters corresponding improved ETKDG by Wang, Witek,
                //! Landrum and Riniker (10.1021/acs.jcim.0c00025) - the
                //! macrocycle part
  return new DGeomHelpers::EmbedParameters(DGeomHelpers::ETKDGv3);
}
DGeomHelpers::EmbedParameters *
getsrETKDGv3() {  //! Parameters corresponding improved ETKDG by Wang, Witek,
                  //! Landrum and Riniker (10.1021/acs.jcim.0c00025) - the
                  //! macrocycle part
  return new DGeomHelpers::EmbedParameters(DGeomHelpers::srETKDGv3);
}
DGeomHelpers::EmbedParameters *getKDG() {
  return new DGeomHelpers::EmbedParameters(DGeomHelpers::KDG);
}
DGeomHelpers::EmbedParameters *getETDG() {
  return new DGeomHelpers::EmbedParameters(DGeomHelpers::ETDG);
}

void setCPCI(DGeomHelpers::EmbedParameters *self, python::dict &CPCIdict) {
  // CPCI has the atom pair tuple as key and charge product as value
  std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>> CPCI(
      new std::map<std::pair<unsigned int, unsigned int>, double>);

  python::list ks = CPCIdict.keys();
  unsigned int nKeys = python::extract<unsigned int>(ks.attr("__len__")());

  for (unsigned int i = 0; i < nKeys; ++i) {
    python::tuple id = python::extract<python::tuple>(ks[i]);
    unsigned int a = python::extract<unsigned int>(id[0]);
    unsigned int b = python::extract<unsigned int>(id[1]);
    (*CPCI)[std::make_pair(a, b)] = python::extract<double>(CPCIdict[id]);
  }

  self->CPCI = CPCI;
}

void setBoundsMatrix(DGeomHelpers::EmbedParameters *self,
                     python::object boundsMatArg) {
  PyObject *boundsMatObj = boundsMatArg.ptr();
  if (!PyArray_Check(boundsMatObj)) {
    throw_value_error("Argument isn't an array");
  }

  auto *boundsMat = reinterpret_cast<PyArrayObject *>(boundsMatObj);
  // get the dimensions of the array
  int nrows = PyArray_DIM(boundsMat, 0);
  int ncols = PyArray_DIM(boundsMat, 1);
  if (nrows != ncols) {
    throw_value_error("The array has to be square");
  }
  if (nrows <= 0) {
    throw_value_error("The array has to have a nonzero size");
  }
  if (PyArray_DESCR(boundsMat)->type_num != NPY_DOUBLE) {
    throw_value_error("Only double arrays are currently supported");
  }

  unsigned int dSize = nrows * nrows;
  auto *cData = new double[dSize];
  auto *inData = reinterpret_cast<double *>(PyArray_DATA(boundsMat));
  memcpy(static_cast<void *>(cData), static_cast<const void *>(inData),
         dSize * sizeof(double));
  DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
  self->boundsMat = boost::shared_ptr<const DistGeom::BoundsMatrix>(
      new DistGeom::BoundsMatrix(nrows, sdata));
}

}  // namespace RDKit

BOOST_PYTHON_MODULE(rdDistGeom) {
  python::scope().attr("__doc__") =
      "Module containing functions to compute atomic coordinates in 3D using "
      "distance geometry";

  rdkit_import_array();

  // RegisterListConverter<RDKit::Atom*>();

  std::string docString =
      "Use distance geometry to obtain initial \n\
 coordinates for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - maxAttempts : the maximum number of attempts to try embedding \n\
    - randomSeed : provide a seed for the random number generator \n\
                   so that the same coordinates can be obtained \n\
                   for a molecule on multiple runs. If -1, the \n\
                   RNG will not be seeded. \n\
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
    - numZeroFail : fail embedding if we have at least this many zero eigenvalues \n\
    - coordMap : a dictionary mapping atom IDs->coordinates. Use this to \n\
                 require some atoms to have fixed coordinates in the resulting \n\
                 conformation.\n\
    - forceTol : tolerance to be used during the force-field minimization with \n\
                 the distance geometry force field.\n\
    - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing\n\
                 of the bounds matrix fails.\n\
    - enforceChirality : enforce the correct chirality if chiral centers are present.\n\
    - useExpTorsionAnglePrefs : impose experimental torsion angle preferences\n\
    - useBasicKnowledge : impose basic knowledge such as flat rings\n\
    - printExpTorsionAngles : print the output from the experimental torsion angles\n\
\n\
 RETURNS:\n\n\
    ID of the new conformation added to the molecule \n\
\n";
  python::def(
      "EmbedMolecule", RDKit::EmbedMolecule,
      (python::arg("mol"), python::arg("maxAttempts") = 0,
       python::arg("randomSeed") = -1, python::arg("clearConfs") = true,
       python::arg("useRandomCoords") = false, python::arg("boxSizeMult") = 2.0,
       python::arg("randNegEig") = true, python::arg("numZeroFail") = 1,
       python::arg("coordMap") = python::dict(), python::arg("forceTol") = 1e-3,
       python::arg("ignoreSmoothingFailures") = false,
       python::arg("enforceChirality") = true,
       python::arg("useExpTorsionAnglePrefs") = true,
       python::arg("useBasicKnowledge") = true,
       python::arg("printExpTorsionAngles") = false,
       python::arg("useSmallRingTorsions") = false,
       python::arg("useMacrocycleTorsions") = false,
       python::arg("ETversion") = 1),
      docString.c_str());

  docString =
      "Use distance geometry to obtain multiple sets of \n\
 coordinates for a molecule\n\
 \n\
 ARGUMENTS:\n\n\
  - mol : the molecule of interest\n\
  - numConfs : the number of conformers to generate \n\
  - maxAttempts : the maximum number of attempts to try embedding \n\
  - randomSeed : provide a seed for the random number generator \n\
                 so that the same coordinates can be obtained \n\
                 for a molecule on multiple runs. If -1, the \n\
                 RNG will not be seeded. \n\
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
  - numZeroFail : fail embedding if we have at least this many zero eigenvalues \n\
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
  - enforceChirality : enforce the correct chirality if chiral centers are present.\n\
  - numThreads : number of threads to use while embedding. This only has an effect if the RDKit\n\
               was built with multi-thread support.\n\
              If set to zero, the max supported by the system will be used.\n\
  - useExpTorsionAnglePrefs : impose experimental torsion angle preferences\n\
  - useBasicKnowledge : impose basic knowledge such as flat rings\n\
  - printExpTorsionAngles : print the output from the experimental torsion angles\n\
 RETURNS:\n\n\
    List of new conformation IDs \n\
\n";
  python::def(
      "EmbedMultipleConfs", RDKit::EmbedMultipleConfs,
      (python::arg("mol"), python::arg("numConfs") = 10,
       python::arg("maxAttempts") = 0, python::arg("randomSeed") = -1,
       python::arg("clearConfs") = true, python::arg("useRandomCoords") = false,
       python::arg("boxSizeMult") = 2.0, python::arg("randNegEig") = true,
       python::arg("numZeroFail") = 1, python::arg("pruneRmsThresh") = -1.0,
       python::arg("coordMap") = python::dict(), python::arg("forceTol") = 1e-3,
       python::arg("ignoreSmoothingFailures") = false,
       python::arg("enforceChirality") = true, python::arg("numThreads") = 1,
       python::arg("useExpTorsionAnglePrefs") = true,
       python::arg("useBasicKnowledge") = true,
       python::arg("printExpTorsionAngles") = false,
       python::arg("useSmallRingTorsions") = false,
       python::arg("useMacrocycleTorsions") = false,
       python::arg("ETversion") = 1),
      docString.c_str());

  python::class_<RDKit::DGeomHelpers::EmbedParameters, boost::noncopyable>(
      "EmbedParameters", "Parameters controlling embedding")
      .def_readwrite("maxIterations",
                     &RDKit::DGeomHelpers::EmbedParameters::maxIterations,
                     "maximum number of embedding attempts to use for a "
                     "single conformation")
      .def_readwrite(
          "numThreads", &RDKit::DGeomHelpers::EmbedParameters::numThreads,
          "number of threads to use when embedding multiple conformations")
      .def_readwrite("randomSeed",
                     &RDKit::DGeomHelpers::EmbedParameters::randomSeed,
                     "seed for the random number generator")
      .def_readwrite("clearConfs",
                     &RDKit::DGeomHelpers::EmbedParameters::clearConfs,
                     "clear all existing conformations on the molecule")
      .def_readwrite("useRandomCoords",
                     &RDKit::DGeomHelpers::EmbedParameters::useRandomCoords,
                     "start the embedding from random coordinates instead of "
                     "using eigenvalues of the distance matrix")
      .def_readwrite(
          "boxSizeMult", &RDKit::DGeomHelpers::EmbedParameters::boxSizeMult,
          "determines the size of the box used for random coordinates")
      .def_readwrite("randNegEig",
                     &RDKit::DGeomHelpers::EmbedParameters::randNegEig,
                     "if the embedding yields a negative eigenvalue, pick "
                     "coordinates that correspond to this component at random")
      .def_readwrite(
          "numZeroFail", &RDKit::DGeomHelpers::EmbedParameters::numZeroFail,
          "fail embedding if we have at least this many zero eigenvalues")
      .def_readwrite("optimizerForceTol",
                     &RDKit::DGeomHelpers::EmbedParameters::optimizerForceTol,
                     "the tolerance to be used during the distance-geometry "
                     "force field minimization")
      .def_readwrite(
          "ignoreSmoothingFailures",
          &RDKit::DGeomHelpers::EmbedParameters::ignoreSmoothingFailures,
          "try and embed the molecule if if triangle smoothing of "
          "the bounds matrix fails")
      .def_readwrite("enforceChirality",
                     &RDKit::DGeomHelpers::EmbedParameters::enforceChirality,
                     "enforce correct chirilaty if chiral centers are present")
      .def_readwrite(
          "useExpTorsionAnglePrefs",
          &RDKit::DGeomHelpers::EmbedParameters::useExpTorsionAnglePrefs,
          "impose experimental torsion angle preferences")
      .def_readwrite("useBasicKnowledge",
                     &RDKit::DGeomHelpers::EmbedParameters::useBasicKnowledge,
                     "impose basic-knowledge constraints such as flat rings")
      .def_readwrite("ETversion",
                     &RDKit::DGeomHelpers::EmbedParameters::ETversion,
                     "version of the experimental torsion-angle preferences")
      .def_readwrite("verbose", &RDKit::DGeomHelpers::EmbedParameters::verbose,
                     "be verbose about configuration")
      .def_readwrite("pruneRmsThresh",
                     &RDKit::DGeomHelpers::EmbedParameters::pruneRmsThresh,
                     "used to filter multiple conformations: keep only "
                     "conformations that are at least this far apart from each "
                     "other")
      .def_readwrite(
          "onlyHeavyAtomsForRMS",
          &RDKit::DGeomHelpers::EmbedParameters::onlyHeavyAtomsForRMS,
          "Only consider heavy atoms when doing RMS filtering")
      .def_readwrite(
          "embedFragmentsSeparately",
          &RDKit::DGeomHelpers::EmbedParameters::embedFragmentsSeparately,
          "split the molecule into fragments and embed them separately")
      .def_readwrite(
          "useSmallRingTorsions",
          &RDKit::DGeomHelpers::EmbedParameters::useSmallRingTorsions,
          "impose small ring torsion angle preferences")
      .def_readwrite(
          "useMacrocycleTorsions",
          &RDKit::DGeomHelpers::EmbedParameters::useMacrocycleTorsions,
          "impose macrocycle torsion angle preferences")
      .def_readwrite(
          "boundsMatForceScaling",
          &RDKit::DGeomHelpers::EmbedParameters::boundsMatForceScaling,
          "scale the weights of the atom pair distance restraints relative to "
          "the other types of restraints")
      .def_readwrite(
          "useSymmetryForPruning",
          &RDKit::DGeomHelpers::EmbedParameters::useSymmetryForPruning,
          "use molecule symmetry when doing the RMSD pruning. Note that this "
          "option automatically also sets onlyHeavyAtomsForRMS to true.")
      .def("SetBoundsMat", &RDKit::setBoundsMatrix,
           "set the distance-bounds matrix to be used (no triangle smoothing "
           "will be done on this) from a Numpy array")
      .def("SetCPCI", &RDKit::setCPCI,
           "set the customised pairwise Columb-like interaction to atom pairs."
           "used during structural minimisation stage")
      .def_readwrite("forceTransAmides",
                     &RDKit::DGeomHelpers::EmbedParameters::forceTransAmides,
                     "constrain amide bonds to be trans");
  docString =
      "Use distance geometry to obtain multiple sets of \n\
 coordinates for a molecule\n\
 \n\
 ARGUMENTS:\n\n\
  - mol : the molecule of interest\n\
  - numConfs : the number of conformers to generate \n\
  - params : an EmbedParameters object \n\
 RETURNS:\n\n\
    List of new conformation IDs \n\
\n";
  python::def(
      "EmbedMultipleConfs", RDKit::EmbedMultipleConfs2,
      (python::arg("mol"), python::arg("numConfs"), python::arg("params")),
      docString.c_str());

  docString =
      "Use distance geometry to obtain intial \n\
 coordinates for a molecule\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - params : an EmbedParameters object \n\
\n\
 RETURNS:\n\n\
    ID of the new conformation added to the molecule \n\
\n";
  python::def("EmbedMolecule", RDKit::EmbedMolecule2,
              (python::arg("mol"), python::arg("params")), docString.c_str());
  python::def(
      "ETKDG", RDKit::getETKDG,
      "Returns an EmbedParameters object for the ETKDG method - version 1.",
      python::return_value_policy<python::manage_new_object>());
  python::def(
      "ETKDGv2", RDKit::getETKDGv2,
      "Returns an EmbedParameters object for the ETKDG method - version 2.",
      python::return_value_policy<python::manage_new_object>());
  python::def("srETKDGv3", RDKit::getsrETKDGv3,
              "Returns an EmbedParameters object for the ETKDG method - "
              "version 3 (small rings).",
              python::return_value_policy<python::manage_new_object>());
  python::def("ETKDGv3", RDKit::getETKDGv3,
              "Returns an EmbedParameters object for the ETKDG method - "
              "version 3 (macrocycles).",
              python::return_value_policy<python::manage_new_object>());
  python::def("ETDG", RDKit::getETDG,
              "Returns an EmbedParameters object for the ETDG method.",
              python::return_value_policy<python::manage_new_object>());
  python::def("KDG", RDKit::getKDG,
              "Returns an EmbedParameters object for the KDG method.",
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the distance bounds matrix for a molecule\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - set15bounds : set bounds for 1-5 atom distances based on \n\
                    topology (otherwise stop at 1-4s)\n\
    - scaleVDW : scale down the sum of VDW radii when setting the \n\
                 lower bounds for atoms less than 5 bonds apart \n\
    - doTriangleSmoothing : run triangle smoothing on the bounds \n\
                 matrix before returning it \n\
 RETURNS:\n\n\
    the bounds matrix as a Numeric array with lower bounds in \n\
    the lower triangle and upper bounds in the upper triangle\n\
\n";
  python::def("GetMoleculeBoundsMatrix", RDKit::getMolBoundsMatrix,
              (python::arg("mol"), python::arg("set15bounds") = true,
               python::arg("scaleVDW") = false,
               python::arg("doTriangleSmoothing") = true,
               python::arg("useMacrocycle14config") = false),
              docString.c_str());
}
