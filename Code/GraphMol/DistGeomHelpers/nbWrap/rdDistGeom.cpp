//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <RDBoost/Wrap_nb.h>

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>

#include <GraphMol/GraphMol.h>
#include <RDGeneral/ControlCHandler.h>

#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
struct PyEmbedParameters : public RDKit::DGeomHelpers::EmbedParameters {
 public:
  PyEmbedParameters() : RDKit::DGeomHelpers::EmbedParameters() {}
  PyEmbedParameters(const RDKit::DGeomHelpers::EmbedParameters &other)
      : RDKit::DGeomHelpers::EmbedParameters(other) {}

  void setCoordMap(const nb::dict &cmap) {
    d_coordMap.reset(new std::map<int, RDGeom::Point3D>());
    for (auto item : cmap) {
      (*d_coordMap)[nb::cast<int>(item.first)] =
          nb::cast<RDGeom::Point3D>(item.second);
    }
    coordMap = d_coordMap.get();
  }

  nb::tuple getFailureCounts() const {
    nb::list lst;
    for (auto failure : failures) {
      lst.append(failure);
    }
    return nb::tuple(lst);
  }

  void setCPCI(const nb::dict &CPCIdict) {
    CPCI = std::make_shared<
        std::map<std::pair<unsigned int, unsigned int>, double>>();
    for (auto item : CPCIdict) {
      auto id = nb::cast<nb::tuple>(item.first);
      unsigned int a = nb::cast<unsigned int>(id[0]);
      unsigned int b = nb::cast<unsigned int>(id[1]);
      (*CPCI)[{a, b}] = nb::cast<double>(item.second);
    }
  }

  void setBoundsMatrix(
      nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig> bm) {
    if (bm.shape(0) != bm.shape(1)) {
      throw nb::value_error("The array has to be square");
    }
    if (bm.shape(0) == 0) {
      throw nb::value_error("The array has to have a nonzero size");
    }
    unsigned int nrows = (unsigned int)bm.shape(0);
    unsigned int dSize = nrows * nrows;
    auto *cData = new double[dSize];
    memcpy(cData, bm.data(), dSize * sizeof(double));
    DistGeom::BoundsMatrix::DATA_SPTR sdata(cData);
    boundsMat = boost::shared_ptr<const DistGeom::BoundsMatrix>(
        new DistGeom::BoundsMatrix(nrows, sdata));
  }

 private:
  std::unique_ptr<std::map<int, RDGeom::Point3D>> d_coordMap;
};
}  // namespace

namespace RDKit {

static int EmbedMolecule(ROMol &mol, unsigned int maxAttempts, int seed,
                         bool clearConfs, bool useRandomCoords,
                         double boxSizeMult, bool randNegEig,
                         unsigned int numZeroFail, nb::dict coordMap,
                         double forceTol, bool ignoreSmoothingFailures,
                         bool enforceChirality, bool useExpTorsionAnglePrefs,
                         bool useBasicKnowledge, bool printExpTorsionAngles,
                         bool useSmallRingTorsions, bool useMacrocycleTorsions,
                         unsigned int ETversion, bool useMacrocycle14config) {
  std::map<int, RDGeom::Point3D> pMap;
  for (auto item : coordMap) {
    pMap[nb::cast<int>(item.first)] = nb::cast<RDGeom::Point3D>(item.second);
  }
  std::map<int, RDGeom::Point3D> *pMapPtr = pMap.empty() ? nullptr : &pMap;

  const double basinThresh = DGeomHelpers::EmbedParameters().basinThresh;
  DGeomHelpers::EmbedParameters params(
      maxAttempts, 1, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, pMapPtr, forceTol, ignoreSmoothingFailures,
      enforceChirality, useExpTorsionAnglePrefs, useBasicKnowledge,
      printExpTorsionAngles, basinThresh, -1., false, ETversion, nullptr, true,
      useSmallRingTorsions, useMacrocycleTorsions, useMacrocycle14config);

  int res;
  {
    nb::gil_scoped_release release;
    res = DGeomHelpers::EmbedMolecule(mol, params);
  }
  if (ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Embedding cancelled");
    throw nb::python_error();
  }
  return res;
}

static int EmbedMolecule2(ROMol &mol, PyEmbedParameters &params) {
  int res;
  {
    nb::gil_scoped_release release;
    res = DGeomHelpers::EmbedMolecule(mol, params);
  }
  if (ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Embedding cancelled");
    throw nb::python_error();
  }
  return res;
}

static INT_VECT EmbedMultipleConfs(
    ROMol &mol, unsigned int numConfs, unsigned int maxAttempts, int seed,
    bool clearConfs, bool useRandomCoords, double boxSizeMult, bool randNegEig,
    unsigned int numZeroFail, double pruneRmsThresh, nb::dict coordMap,
    double forceTol, bool ignoreSmoothingFailures, bool enforceChirality,
    int numThreads, bool useExpTorsionAnglePrefs, bool useBasicKnowledge,
    bool printExpTorsionAngles, bool useSmallRingTorsions,
    bool useMacrocycleTorsions, unsigned int ETversion,
    bool useMacrocycle14config) {
  std::map<int, RDGeom::Point3D> pMap;
  for (auto item : coordMap) {
    pMap[nb::cast<int>(item.first)] = nb::cast<RDGeom::Point3D>(item.second);
  }
  std::map<int, RDGeom::Point3D> *pMapPtr = pMap.empty() ? nullptr : &pMap;

  const double basinThresh = DGeomHelpers::EmbedParameters().basinThresh;
  DGeomHelpers::EmbedParameters params(
      maxAttempts, numThreads, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, pMapPtr, forceTol, ignoreSmoothingFailures,
      enforceChirality, useExpTorsionAnglePrefs, useBasicKnowledge,
      printExpTorsionAngles, basinThresh, pruneRmsThresh, false, ETversion,
      nullptr, true, useSmallRingTorsions, useMacrocycleTorsions,
      useMacrocycle14config);

  INT_VECT res;
  {
    nb::gil_scoped_release release;
    DGeomHelpers::EmbedMultipleConfs(mol, res, numConfs, params);
  }
  if (ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Embedding cancelled");
    throw nb::python_error();
  }
  return res;
}

static INT_VECT EmbedMultipleConfs2(ROMol &mol, unsigned int numConfs,
                                    PyEmbedParameters &params) {
  INT_VECT res;
  {
    nb::gil_scoped_release release;
    DGeomHelpers::EmbedMultipleConfs(mol, res, numConfs, params);
  }
  if (ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Embedding cancelled");
    throw nb::python_error();
  }
  return res;
}

static nb::ndarray<nb::numpy, double, nb::ndim<2>> getMolBoundsMatrix(
    const ROMol &mol, bool set15bounds, bool scaleVDW,
    bool doTriangleSmoothing, bool useMacrocycle14config) {
  unsigned int nats = mol.getNumAtoms();
  DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
  DGeomHelpers::initBoundsMat(mat);
  DGeomHelpers::setTopolBounds(mol, mat, set15bounds, scaleVDW,
                               useMacrocycle14config);
  if (doTriangleSmoothing) {
    DistGeom::triangleSmoothBounds(mat);
  }
  auto *resData = new double[nats * nats];
  memcpy(resData, mat->getData(), nats * nats * sizeof(double));
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(resData, {nats, nats},
                                                     owner);
}

static nb::list getExpTorsHelper(const ROMol &mol, bool useExpTorsions,
                                 bool useSmallRingTorsions,
                                 bool useMacrocycleTorsions,
                                 bool useBasicKnowledge, unsigned int version,
                                 bool verbose) {
  ForceFields::CrystalFF::CrystalFFDetails details;
  std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                         const ForceFields::CrystalFF::ExpTorsionAngle *>>
      torsionBonds;
  ForceFields::CrystalFF::getExperimentalTorsions(
      mol, details, torsionBonds, useExpTorsions, useSmallRingTorsions,
      useMacrocycleTorsions, useBasicKnowledge, version, verbose);
  nb::list result;
  for (const auto &pr : torsionBonds) {
    nb::dict d;
    d["bondIndex"] = std::get<0>(pr);
    d["torsionIndex"] = std::get<2>(pr)->torsionIdx;
    d["smarts"] = std::get<2>(pr)->smarts;
    d["V"] = std::get<2>(pr)->V;
    d["signs"] = std::get<2>(pr)->signs;
    d["atomIndices"] = std::get<1>(pr);
    result.append(d);
  }
  return result;
}

}  // namespace RDKit

NB_MODULE(rdDistGeom, m) {
  m.doc() =
      R"DOC(Module containing functions to compute atomic coordinates in 3D using
distance geometry)DOC";

  m.def(
      "GetExperimentalTorsions",
      [](const RDKit::ROMol &mol, bool useExpTorsionAnglePrefs,
         bool useSmallRingTorsions, bool useMacrocycleTorsions,
         bool useBasicKnowledge, unsigned int ETversion,
         bool printExpTorsionAngles) {
        return RDKit::getExpTorsHelper(mol, useExpTorsionAnglePrefs,
                                      useSmallRingTorsions,
                                      useMacrocycleTorsions, useBasicKnowledge,
                                      ETversion, printExpTorsionAngles);
      },
      "mol"_a, "useExpTorsionAnglePrefs"_a = true,
      "useSmallRingTorsions"_a = false, "useMacrocycleTorsions"_a = true,
      "useBasicKnowledge"_a = true, "ETversion"_a = 2,
      "printExpTorsionAngles"_a = false,
      "returns information about the bonds corresponding to experimental torsions");

  m.def(
      "GetExperimentalTorsions",
      [](const RDKit::ROMol &mol, const PyEmbedParameters &ps) {
        return RDKit::getExpTorsHelper(
            mol, ps.useExpTorsionAnglePrefs, ps.useSmallRingTorsions,
            ps.useMacrocycleTorsions, ps.useBasicKnowledge, ps.ETversion,
            ps.verbose);
      },
      "mol"_a, "embedParams"_a,
      "returns information about the bonds corresponding to experimental torsions");

  m.def(
      "EmbedMolecule", &RDKit::EmbedMolecule,
      "mol"_a, "maxAttempts"_a = 0, "randomSeed"_a = -1,
      "clearConfs"_a = true, "useRandomCoords"_a = false,
      "boxSizeMult"_a = 2.0, "randNegEig"_a = true, "numZeroFail"_a = 1,
      "coordMap"_a = nb::dict(), "forceTol"_a = 1e-3,
      "ignoreSmoothingFailures"_a = false, "enforceChirality"_a = true,
      "useExpTorsionAnglePrefs"_a = true, "useBasicKnowledge"_a = true,
      "printExpTorsionAngles"_a = false, "useSmallRingTorsions"_a = false,
      "useMacrocycleTorsions"_a = true, "ETversion"_a = 2,
      "useMacrocycle14config"_a = true,
      R"DOC(Use distance geometry to obtain initial
coordinates for a molecule

ARGUMENTS:

   - mol : the molecule of interest
   - maxAttempts : maximum number of embedding attempts to use for a single conformation
   - randomSeed : provide a seed for the random number generator
                  so that the same coordinates can be obtained
                  for a molecule on multiple runs. If -1, the
                  RNG will not be seeded.
   - clearConfs : clear all existing conformations on the molecule
   - useRandomCoords : Start the embedding from random coordinates instead of
                       using eigenvalues of the distance matrix.
   - boxSizeMult :  Determines the size of the box that is used for
                    random coordinates. If this is a positive number, the
                    side length will equal the largest element of the distance
                    matrix times boxSizeMult. If this is a negative number,
                    the side length will equal -boxSizeMult (i.e. independent
                    of the elements of the distance matrix).
   - randNegEig : If the embedding yields a negative eigenvalue,
                  pick coordinates that correspond
                  to this component at random
   - numZeroFail : fail embedding if we have at least this many zero eigenvalues
   - coordMap : a dictionary mapping atom IDs->coordinates. Use this to
                require some atoms to have fixed coordinates in the resulting
                conformation.
   - forceTol : tolerance to be used during the force-field minimization with
                the distance geometry force field.
   - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                of the bounds matrix fails.
   - enforceChirality : enforce the correct chirality if chiral centers are present.
   - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
   - useBasicKnowledge : impose basic knowledge such as flat rings
   - printExpTorsionAngles : print the output from the experimental torsion angles
   - useMacrocycleTorsions : use additional torsion profiles for macrocycles
   - ETversion : version of the standard torsion definitions to use. NOTE for both
                 ETKDGv2 and ETKDGv3 this should be 2 since ETKDGv3 uses the ETKDGv2
                 definitions for standard torsions
   - useMacrocycle14config : use the 1-4 distance bounds from ETKDGv3

RETURNS:

   ID of the new conformation added to the molecule or -1 if the embedding fails.
)DOC");

  m.def(
      "EmbedMultipleConfs", &RDKit::EmbedMultipleConfs,
      "mol"_a, "numConfs"_a = 10, "maxAttempts"_a = 0, "randomSeed"_a = -1,
      "clearConfs"_a = true, "useRandomCoords"_a = false,
      "boxSizeMult"_a = 2.0, "randNegEig"_a = true, "numZeroFail"_a = 1,
      "pruneRmsThresh"_a = -1.0, "coordMap"_a = nb::dict(),
      "forceTol"_a = 1e-3, "ignoreSmoothingFailures"_a = false,
      "enforceChirality"_a = true, "numThreads"_a = 1,
      "useExpTorsionAnglePrefs"_a = true, "useBasicKnowledge"_a = true,
      "printExpTorsionAngles"_a = false, "useSmallRingTorsions"_a = false,
      "useMacrocycleTorsions"_a = true, "ETversion"_a = 2,
      "useMacrocycle14config"_a = true,
      R"DOC(Use distance geometry to obtain multiple sets of
coordinates for a molecule

ARGUMENTS:

  - mol : the molecule of interest
  - numConfs : the number of conformers to generate
  - maxAttempts : maximum number of embedding attempts to use for a single conformation
  - randomSeed : provide a seed for the random number generator
                 so that the same coordinates can be obtained
                 for a molecule on multiple runs. If -1, the
                 RNG will not be seeded.
  - clearConfs : clear all existing conformations on the molecule
  - useRandomCoords : Start the embedding from random coordinates instead of
                      using eigenvalues of the distance matrix.
  - boxSizeMult    Determines the size of the box that is used for
                   random coordinates. If this is a positive number, the
                   side length will equal the largest element of the distance
                   matrix times boxSizeMult. If this is a negative number,
                   the side length will equal -boxSizeMult (i.e. independent
                   of the elements of the distance matrix).
  - randNegEig : If the embedding yields a negative eigenvalue,
                 pick coordinates that correspond
                 to this component at random
  - numZeroFail : fail embedding if we have at least this many zero eigenvalues
  - pruneRmsThresh : Retain only the conformations out of 'numConfs'
                    after embedding that are at least
                    this far apart from each other.
                    RMSD is computed on the heavy atoms.
                    Pruning is greedy; i.e. the first embedded conformation
                    is retained and from then on only those that are at
                    least pruneRmsThresh away from all retained conformations
                    are kept. The pruning is done after embedding and
                    bounds violation minimization. No pruning by default.
  - coordMap : a dictionary mapping atom IDs->coordinates. Use this to
               require some atoms to have fixed coordinates in the resulting
               conformation.
  - forceTol : tolerance to be used during the force-field minimization with
               the distance geometry force field.
  - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
               of the bounds matrix fails.
  - enforceChirality : enforce the correct chirality if chiral centers are present.
  - numThreads : number of threads to use while embedding. This only has an effect if the RDKit
               was built with multi-thread support.
              If set to zero, the max supported by the system will be used.
  - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
  - useBasicKnowledge : impose basic knowledge such as flat rings
  - printExpTorsionAngles : print the output from the experimental torsion angles

RETURNS:

   Iterator which yields new conformation IDs
)DOC");

  nb::enum_<RDKit::DGeomHelpers::EmbedFailureCauses>(m, "EmbedFailureCauses")
      .value("INITIAL_COORDS",
             RDKit::DGeomHelpers::EmbedFailureCauses::INITIAL_COORDS)
      .value("FIRST_MINIMIZATION",
             RDKit::DGeomHelpers::EmbedFailureCauses::FIRST_MINIMIZATION)
      .value("CHECK_TETRAHEDRAL_CENTERS",
             RDKit::DGeomHelpers::EmbedFailureCauses::CHECK_TETRAHEDRAL_CENTERS)
      .value("CHECK_CHIRAL_CENTERS",
             RDKit::DGeomHelpers::EmbedFailureCauses::CHECK_CHIRAL_CENTERS)
      .value(
          "MINIMIZE_FOURTH_DIMENSION",
          RDKit::DGeomHelpers::EmbedFailureCauses::MINIMIZE_FOURTH_DIMENSION)
      .value("ETK_MINIMIZATION",
             RDKit::DGeomHelpers::EmbedFailureCauses::ETK_MINIMIZATION)
      .value("FINAL_CHIRAL_BOUNDS",
             RDKit::DGeomHelpers::EmbedFailureCauses::FINAL_CHIRAL_BOUNDS)
      .value("FINAL_CENTER_IN_VOLUME",
             RDKit::DGeomHelpers::EmbedFailureCauses::FINAL_CENTER_IN_VOLUME)
      .value("LINEAR_DOUBLE_BOND",
             RDKit::DGeomHelpers::EmbedFailureCauses::LINEAR_DOUBLE_BOND)
      .value("BAD_DOUBLE_BOND_STEREO",
             RDKit::DGeomHelpers::EmbedFailureCauses::BAD_DOUBLE_BOND_STEREO)
      .value("CHECK_CHIRAL_CENTERS2",
             RDKit::DGeomHelpers::EmbedFailureCauses::CHECK_CHIRAL_CENTERS2)
      .value("EXCEEDED_TIMEOUT",
             RDKit::DGeomHelpers::EmbedFailureCauses::EXCEEDED_TIMEOUT)
      .export_values();

  nb::class_<PyEmbedParameters>(m, "EmbedParameters",
                                "Parameters controlling embedding")
      .def(nb::init<>())
      .def_rw("maxIterations", &PyEmbedParameters::maxIterations,
              R"DOC(maximum number of embedding attempts to use for a
single conformation)DOC")
      .def_rw(
          "numThreads", &PyEmbedParameters::numThreads,
          "number of threads to use when embedding multiple conformations")
      .def_rw("timeout", &PyEmbedParameters::timeout,
              R"DOC(maximum time in seconds to generate a conformer for a
single molecule fragment. If set to 0, no timeout is set)DOC")
      .def_rw("randomSeed", &PyEmbedParameters::randomSeed,
              "seed for the random number generator")
      .def_rw("clearConfs", &PyEmbedParameters::clearConfs,
              "clear all existing conformations on the molecule")
      .def_rw("useRandomCoords", &PyEmbedParameters::useRandomCoords,
              R"DOC(start the embedding from random coordinates instead of
using eigenvalues of the distance matrix)DOC")
      .def_rw("boxSizeMult", &PyEmbedParameters::boxSizeMult,
              "determines the size of the box used for random coordinates")
      .def_rw("randNegEig", &PyEmbedParameters::randNegEig,
              R"DOC(if the embedding yields a negative eigenvalue, pick
coordinates that correspond to this component at random)DOC")
      .def_rw(
          "numZeroFail", &PyEmbedParameters::numZeroFail,
          "fail embedding if we have at least this many zero eigenvalues")
      .def_rw("optimizerForceTol", &PyEmbedParameters::optimizerForceTol,
              R"DOC(the tolerance to be used during the distance-geometry
force field minimization)DOC")
      .def_rw("basinThresh", &PyEmbedParameters::basinThresh,
              "set the basin threshold for the DGeom force field.")
      .def_rw("ignoreSmoothingFailures",
              &PyEmbedParameters::ignoreSmoothingFailures,
              R"DOC(try and embed the molecule if if triangle smoothing of
the bounds matrix fails)DOC")
      .def_rw("enforceChirality", &PyEmbedParameters::enforceChirality,
              "enforce correct chirilaty if chiral centers are present")
      .def_rw("useExpTorsionAnglePrefs",
              &PyEmbedParameters::useExpTorsionAnglePrefs,
              "impose experimental torsion angle preferences")
      .def_rw("useBasicKnowledge", &PyEmbedParameters::useBasicKnowledge,
              "impose basic-knowledge constraints such as flat rings")
      .def_rw("ETversion", &PyEmbedParameters::ETversion,
              "version of the experimental torsion-angle preferences")
      .def_rw("verbose", &PyEmbedParameters::verbose,
              "be verbose about configuration")
      .def_rw("pruneRmsThresh", &PyEmbedParameters::pruneRmsThresh,
              R"DOC(used to filter multiple conformations: keep only
conformations that are at least this far apart from each other)DOC")
      .def_rw("onlyHeavyAtomsForRMS", &PyEmbedParameters::onlyHeavyAtomsForRMS,
              "Only consider heavy atoms when doing RMS filtering")
      .def_rw(
          "embedFragmentsSeparately",
          &PyEmbedParameters::embedFragmentsSeparately,
          "split the molecule into fragments and embed them separately")
      .def_rw("useSmallRingTorsions", &PyEmbedParameters::useSmallRingTorsions,
              "impose small ring torsion angle preferences")
      .def_rw("useMacrocycleTorsions",
              &PyEmbedParameters::useMacrocycleTorsions,
              "impose macrocycle torsion angle preferences")
      .def_rw("useMacrocycle14config",
              &PyEmbedParameters::useMacrocycle14config,
              "use the 1-4 distance bounds from ETKDGv3")
      .def_rw(
          "boundsMatForceScaling", &PyEmbedParameters::boundsMatForceScaling,
          R"DOC(scale the weights of the atom pair distance restraints relative to
the other types of restraints)DOC")
      .def_rw(
          "useSymmetryForPruning", &PyEmbedParameters::useSymmetryForPruning,
          R"DOC(use molecule symmetry when doing the RMSD pruning. Note that this
option automatically also sets onlyHeavyAtomsForRMS to true.)DOC")
      .def("SetBoundsMat", &PyEmbedParameters::setBoundsMatrix,
           "boundsMatArg"_a,
           R"DOC(set the distance-bounds matrix to be used (no triangle smoothing
will be done on this) from a Numpy array)DOC")
      .def("SetCPCI", &PyEmbedParameters::setCPCI, "CPCIdict"_a,
           R"DOC(set the customised pairwise Columb-like interaction to atom pairs.
used during structural minimisation stage)DOC")
      .def_rw("forceTransAmides", &PyEmbedParameters::forceTransAmides,
              "constrain amide bonds to be trans")
      .def_rw(
          "trackFailures", &PyEmbedParameters::trackFailures,
          "keep track of which checks during the embedding process fail")
      .def("GetFailureCounts", &PyEmbedParameters::getFailureCounts,
           "returns the counts of each failure type")
      .def_rw(
          "enableSequentialRandomSeeds",
          &PyEmbedParameters::enableSequentialRandomSeeds,
          "handle random number seeds so that conformer generation can be restarted")
      .def_rw(
          "symmetrizeConjugatedTerminalGroupsForPruning",
          &PyEmbedParameters::symmetrizeConjugatedTerminalGroupsForPruning,
          "symmetrize terminal conjugated groups for RMSD pruning")
      .def("SetCoordMap", &PyEmbedParameters::setCoordMap,
           "sets the coordmap to be used")
      .def("__setattr__", &safeSetattr);

  m.def(
      "EmbedMultipleConfs", &RDKit::EmbedMultipleConfs2,
      "mol"_a, "numConfs"_a, "params"_a,
      R"DOC(Use distance geometry to obtain multiple sets of
coordinates for a molecule

ARGUMENTS:

  - mol : the molecule of interest
  - numConfs : the number of conformers to generate
  - params : an EmbedParameters object

RETURNS:

   Iterator which yields new conformation IDs
)DOC");

  m.def(
      "EmbedMolecule", &RDKit::EmbedMolecule2, "mol"_a, "params"_a,
      R"DOC(Use distance geometry to obtain initial
coordinates for a molecule

ARGUMENTS:

   - mol : the molecule of interest
   - params : an EmbedParameters object

RETURNS:

   ID of the new conformation added to the molecule or -1 if the embedding fails.
)DOC");

  m.def("ETKDG", []() { return PyEmbedParameters(RDKit::DGeomHelpers::ETKDG); },
        "Returns an EmbedParameters object for the ETKDG method - version 1.");
  m.def("ETKDGv2",
        []() { return PyEmbedParameters(RDKit::DGeomHelpers::ETKDGv2); },
        "Returns an EmbedParameters object for the ETKDG method - version 2.");
  m.def("srETKDGv3",
        []() { return PyEmbedParameters(RDKit::DGeomHelpers::srETKDGv3); },
        R"DOC(Returns an EmbedParameters object for the ETKDG method -
version 3 (small rings).)DOC");
  m.def("ETKDGv3",
        []() { return PyEmbedParameters(RDKit::DGeomHelpers::ETKDGv3); },
        R"DOC(Returns an EmbedParameters object for the ETKDG method -
version 3 (macrocycles).)DOC");
  m.def("ETDG", []() { return PyEmbedParameters(RDKit::DGeomHelpers::ETDG); },
        "Returns an EmbedParameters object for the ETDG method.");
  m.def("ETDGv2",
        []() { return PyEmbedParameters(RDKit::DGeomHelpers::ETDGv2); },
        "Returns an EmbedParameters object for the ETDG method - version 2.");
  m.def("KDG", []() { return PyEmbedParameters(RDKit::DGeomHelpers::KDG); },
        "Returns an EmbedParameters object for the KDG method.");

  m.def(
      "GetMoleculeBoundsMatrix", &RDKit::getMolBoundsMatrix,
      "mol"_a, "set15bounds"_a = true, "scaleVDW"_a = false,
      "doTriangleSmoothing"_a = true, "useMacrocycle14config"_a = false,
      R"DOC(Returns the distance bounds matrix for a molecule

ARGUMENTS:

   - mol : the molecule of interest
   - set15bounds : set bounds for 1-5 atom distances based on
                   topology (otherwise stop at 1-4s)
   - scaleVDW : scale down the sum of VDW radii when setting the
                lower bounds for atoms less than 5 bonds apart
   - doTriangleSmoothing : run triangle smoothing on the bounds
                matrix before returning it

RETURNS:

   the bounds matrix as a Numeric array with lower bounds in
   the lower triangle and upper bounds in the upper triangle
)DOC");

  m.def(
      "EmbedParametersToJSON",
      [](const PyEmbedParameters &ps) {
        return embedParametersToJSON(ps);
      },
      "embedParameters"_a,
      R"DOC(Returns json string containing embedParameters attributes

ARGUMENTS:

  - embedParameters : the Params object you want serialized

RETURNS:

  The Params object as json string
)DOC");
}
