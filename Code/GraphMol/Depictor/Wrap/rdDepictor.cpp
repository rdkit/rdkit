// $Id$
//
//  Copyright (C) 2003-2010 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#define PY_ARRAY_UNIQUE_SYMBOL Depictor_array_API
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Depictor/EmbeddedFrag.h>
#include <GraphMol/Depictor/DepictUtils.h>

using namespace RDDepict;

namespace python = boost::python;

void rdDepictExceptionTranslator(DepictException const &e) {
  std::ostringstream oss;
  oss << "Depict error: " << e.what();
  PyErr_SetString(PyExc_ValueError, oss.str().c_str());
}

namespace RDDepict {

unsigned int Compute2DCoords(RDKit::ROMol &mol, bool canonOrient,
                             bool clearConfs, python::dict &coordMap,
                             unsigned int nFlipsPerSample = 3,
                             unsigned int nSamples = 100, int sampleSeed = 100,
                             bool permuteDeg4Nodes = false,
                             double bondLength = -1.0, bool forceRDKit = false,
                             bool useRingTemplates = false) {
  RDGeom::INT_POINT2D_MAP cMap;
  cMap.clear();
  python::list ks = coordMap.keys();
  for (unsigned int i = 0; i < python::len(ks); ++i) {
    unsigned int id = python::extract<unsigned int>(ks[i]);
    if (id >= mol.getNumAtoms()) {
      throw_value_error("atom index out of range");
    }
    cMap[id] = python::extract<RDGeom::Point2D>(coordMap[id]);
  }
  double oBondLen = RDDepict::BOND_LEN;
  if (bondLength > 0) {
    RDDepict::BOND_LEN = bondLength;
  }
  unsigned int res;
  res = RDDepict::compute2DCoords(
      mol, &cMap, canonOrient, clearConfs, nFlipsPerSample, nSamples,
      sampleSeed, permuteDeg4Nodes, forceRDKit, useRingTemplates);
  if (bondLength > 0) {
    RDDepict::BOND_LEN = oBondLen;
  }
  return res;
}

unsigned int Compute2DCoordsMimicDistmat(
    RDKit::ROMol &mol, python::object distMat, bool canonOrient,
    bool clearConfs, double weightDistMat, unsigned int nFlipsPerSample,
    unsigned int nSamples, int sampleSeed, bool permuteDeg4Nodes,
    double bondLength = -1.0, bool forceRDKit = false) {
  PyObject *distMatPtr = distMat.ptr();
  if (!PyArray_Check(distMatPtr)) {
    throw_value_error("Argument isn't an array");
  }

  auto *dmatrix = reinterpret_cast<PyArrayObject *>(distMatPtr);
  unsigned int nitems = PyArray_DIM(dmatrix, 0);
  unsigned int na = mol.getNumAtoms();

  if (nitems != na * (na - 1) / 2) {
    throw_value_error(
        "The array size does not match the number of atoms in the molecule");
  }
  auto *inData = reinterpret_cast<double *>(PyArray_DATA(dmatrix));
  auto *cData = new double[nitems];

  memcpy(static_cast<void *>(cData), static_cast<const void *>(inData),
         nitems * sizeof(double));

  DOUBLE_SMART_PTR dmat(cData);
  double oBondLen = RDDepict::BOND_LEN;
  if (bondLength > 0) {
    RDDepict::BOND_LEN = bondLength;
  }
  unsigned int res;
  res = RDDepict::compute2DCoordsMimicDistMat(
      mol, &dmat, canonOrient, clearConfs, weightDistMat, nFlipsPerSample,
      nSamples, sampleSeed, permuteDeg4Nodes, forceRDKit);
  if (bondLength > 0) {
    RDDepict::BOND_LEN = oBondLen;
  }
  return res;
}

python::tuple GenerateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
    const python::object &refPatt, const ConstrainedDepictionParams &params) {
  RDKit::ROMol *referencePattern = nullptr;
  if (!refPatt.is_none()) {
    referencePattern = python::extract<RDKit::ROMol *>(refPatt);
  }
  auto matchVect = RDDepict::generateDepictionMatching2DStructure(
      mol, reference, confId, referencePattern, params);
  python::list atomMap;
  for (const auto &pair : matchVect) {
    atomMap.append(python::make_tuple(pair.first, pair.second));
  }
  return python::tuple(atomMap);
}

python::tuple GenerateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
    const python::object &refPatt, bool acceptFailure, bool forceRDKit,
    bool allowRGroups) {
  ConstrainedDepictionParams params;
  params.acceptFailure = acceptFailure;
  params.forceRDKit = forceRDKit;
  params.allowRGroups = allowRGroups;
  return GenerateDepictionMatching2DStructure(mol, reference, confId, refPatt,
                                              params);
}

python::tuple GenerateDepictionMatching2DStructure(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
    const python::object &refPatt, const python::object &pyParams) {
  ConstrainedDepictionParams params;
  if (pyParams) {
    params = python::extract<ConstrainedDepictionParams>(pyParams);
  }
  return GenerateDepictionMatching2DStructure(mol, reference, confId, refPatt,
                                              params);
}

void GenerateDepictionMatching2DStructureAtomMap(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const python::object &atomMap, int confId,
    const ConstrainedDepictionParams &params) {
  std::unique_ptr<RDKit::MatchVectType> matchVect(translateAtomMap(atomMap));
  RDDepict::generateDepictionMatching2DStructure(mol, reference, *matchVect,
                                                 confId, params);
}

void GenerateDepictionMatching2DStructureAtomMap(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const python::object &atomMap, int confId, const python::object &pyParams) {
  ConstrainedDepictionParams params;
  if (pyParams) {
    params = python::extract<ConstrainedDepictionParams>(pyParams);
  }
  GenerateDepictionMatching2DStructureAtomMap(mol, reference, atomMap, confId,
                                              params);
}

void GenerateDepictionMatching2DStructureAtomMap(RDKit::ROMol &mol,
                                                 const RDKit::ROMol &reference,
                                                 const python::object &atomMap,
                                                 int confId, bool forceRDKit) {
  ConstrainedDepictionParams params;
  params.forceRDKit = forceRDKit;

  GenerateDepictionMatching2DStructureAtomMap(mol, reference, atomMap, confId,
                                              params);
}

void GenerateDepictionMatching3DStructure(RDKit::ROMol &mol,
                                          RDKit::ROMol &reference, int confId,
                                          python::object refPatt,
                                          bool acceptFailure,
                                          bool forceRDKit = false) {
  RDKit::ROMol *referencePattern = nullptr;
  if (refPatt) {
    referencePattern = python::extract<RDKit::ROMol *>(refPatt);
  }

  RDDepict::generateDepictionMatching3DStructure(
      mol, reference, confId, referencePattern, acceptFailure, forceRDKit);
}

bool isCoordGenSupportAvailable() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return true;
#else
  return false;
#endif
}

void setPreferCoordGen(bool value) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = value;
#endif
}
bool getPreferCoordGen() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return RDDepict::preferCoordGen;
#else
  return false;
#endif
}

class UsingCoordGen : public boost::noncopyable {
 public:
  UsingCoordGen() = delete;
  UsingCoordGen(bool temp_state)
      : m_initial_state{getPreferCoordGen()}, m_temp_state(temp_state) {}
  ~UsingCoordGen() = default;

  void enter() { setPreferCoordGen(m_temp_state); }

  void exit(python::object exc_type, python::object exc_val,
            python::object traceback) {
    RDUNUSED_PARAM(exc_type);
    RDUNUSED_PARAM(exc_val);
    RDUNUSED_PARAM(traceback);
    setPreferCoordGen(m_initial_state);
  }

 private:
  bool m_initial_state;
  bool m_temp_state;
};

}  // namespace RDDepict

BOOST_PYTHON_MODULE(rdDepictor) {
  python::scope().attr("__doc__") =
      "Module containing the functionality to compute 2D coordinates for a "
      "molecule";
  python::register_exception_translator<RDDepict::DepictException>(
      &rdDepictExceptionTranslator);

  rdkit_import_array();

  python::def("IsCoordGenSupportAvailable", isCoordGenSupportAvailable,
              "Returns whether RDKit was built with CoordGen support.");

  python::def("SetPreferCoordGen", setPreferCoordGen, python::arg("val"),
#ifdef RDK_BUILD_COORDGEN_SUPPORT
              "Sets whether or not the CoordGen library should be preferred to "
              "the RDKit depiction library."
#else
              "Has no effect (CoordGen support not enabled)"
#endif
  );
  python::def(
      "SetRingSystemTemplates", RDDepict::setRingSystemTemplates,
      (python::arg("templatePath")),
      "Loads the ring system templates from the specified file to be "
      "used in 2D coordinate generation. Each template must be a single "
      "line in the file represented using CXSMILES, and the structure should "
      "be a single ring system. Throws a DepictException if any templates "
      "are invalid.");
  python::def(
      "AddRingSystemTemplates", RDDepict::addRingSystemTemplates,
      (python::arg("templatePath")),
      "Adds the ring system templates from the specified file to be "
      "used in 2D coordinate generation. If there are duplicates, the most "
      "recently added template will be used. Each template must be a single "
      "line in the file represented using CXSMILES, and the structure should "
      "be a single ring system. Throws a DepictException if any templates "
      "are invalid.");
  python::def(
      "LoadDefaultRingSystemTemplates",
      RDDepict::loadDefaultRingSystemTemplates,
      "Loads the default ring system templates and removes existing ones, if present.");
  python::def(
      "GetPreferCoordGen", getPreferCoordGen,
#ifdef RDK_BUILD_COORDGEN_SUPPORT
      "Return whether or not the CoordGen library is used for coordinate "
      "generation in the RDKit depiction library."
#else
      "Always returns False (CoordGen support not enabled)"
#endif
  );

  python::class_<UsingCoordGen, boost::noncopyable>(
      "UsingCoordGen",
      "Context manager to temporarily set CoordGen library preference in RDKit depiction.",
      python::init<bool>(python::args("self", "temp_state"), "Constructor"))
      .def("__enter__", &UsingCoordGen::enter)
      .def("__exit__", &UsingCoordGen::exit);

  python::class_<ConstrainedDepictionParams, boost::noncopyable>(
      "ConstrainedDepictionParams",
      "Parameters controlling constrained depiction")
      .def_readwrite(
          "acceptFailure", &ConstrainedDepictionParams::acceptFailure,
          "if False (default), a DepictException is thrown if the molecule "
          "does not have a substructure match to the reference; "
          "if True, an unconstrained depiction will be generated")
      .def_readwrite(
          "forceRDKit", &ConstrainedDepictionParams::forceRDKit,
          "if True, use RDKit to generate coordinates even if preferCoordGen "
          "is set to True; defaults to False")
      .def_readwrite(
          "allowRGroups", &ConstrainedDepictionParams::allowRGroups,
          "if True, terminal dummy atoms in the reference are ignored "
          "if they match an implicit hydrogen in the molecule or if they are "
          "attached top a query atom; defaults to False")
      .def_readwrite(
          "alignOnly", &ConstrainedDepictionParams::alignOnly,
          "if False (default), a part of the molecule is hard-constrained "
          "to have the same coordinates as the reference, and the rest of "
          "the molecule is built around it; if True, coordinates "
          "from conformation existingConfId are preserved (if they exist) "
          "or generated without constraints (if they do not exist), then "
          "the conformation is rigid-body aligned to the reference")
      .def_readwrite("adjustMolBlockWedging",
                     &ConstrainedDepictionParams::adjustMolBlockWedging,
                     "if True (default), existing wedging information "
                     "will be updated or cleared as required; if False, "
                     "existing molblock wedging information will "
                     "always be preserved")
      .def_readwrite(
          "existingConfId", &ConstrainedDepictionParams::existingConfId,
          "conformation id whose 2D coordinates should be "
          "rigid-body aligned to the reference (if alignOnly is True), "
          "or used to determine whether existing molblock wedging information "
          "can be preserved following the constrained depiction (if "
          "adjustMolBlockWedging is True")
      .def_readwrite(
          "useRingTemplates", &ConstrainedDepictionParams::useRingTemplates,
          "use templates to generate coordinates of complex ring systems");

  std::string docString;
  docString =
      "Compute 2D coordinates for a molecule. \n\
  The resulting coordinates are stored on each atom of the molecule \n\n\
  ARGUMENTS: \n\n\
     mol - the molecule of interest\n\
     canonOrient - orient the molecule in a canonical way\n\
     clearConfs - if true, all existing conformations on the molecule\n\
             will be cleared\n\
     coordMap - a dictionary mapping atom Ids -> Point2D objects \n\
                with starting coordinates for atoms that should\n\
                have their positions locked.\n\
     nFlipsPerSample - number of rotatable bonds that are\n\
                flipped at random at a time.\n\
     nSample - Number of random samplings of rotatable bonds.\n\
     sampleSeed - seed for the random sampling process.\n\
     permuteDeg4Nodes - allow permutation of bonds at a degree 4\n\
                 node during the sampling process \n\
     bondLength - change the default bond length for depiction \n\
     forceRDKit - use RDKit to generate coordinates even if \n\
                  preferCoordGen is set to true\n\
     useRingTemplates - use templates to generate coordinates of complex\n\
                  ring systems\n\n\
  RETURNS: \n\n\
     ID of the conformation added to the molecule\n";
  python::def(
      "Compute2DCoords", RDDepict::Compute2DCoords,
      (python::arg("mol"), python::arg("canonOrient") = true,
       python::arg("clearConfs") = true,
       python::arg("coordMap") = python::dict(),
       python::arg("nFlipsPerSample") = 0, python::arg("nSample") = 0,
       python::arg("sampleSeed") = 0, python::arg("permuteDeg4Nodes") = false,
       python::arg("bondLength") = -1.0, python::arg("forceRDKit") = false,
       python::arg("useRingTemplates") = false),
      docString.c_str());

  docString =
      "Compute 2D coordinates for a molecule such \n\
  that the inter-atom distances mimic those in a user-provided\n\
  distance matrix. \n\
  The resulting coordinates are stored on each atom of the molecule \n\n\
  ARGUMENTS: \n\n\
     mol - the molecule of interest\n\
     distMat - distance matrix that we want the 2D structure to mimic\n\
     canonOrient - orient the molecule in a canonical way\n\
     clearConfs - if true, all existing conformations on the molecule\n\
             will be cleared\n\
     weightDistMat - weight assigned in the cost function to mimicking\n\
                     the distance matrix.\n\
                     This must be between (0.0,1.0). (1.0-weightDistMat)\n\
                     is then the weight assigned to improving \n\
                     the density of the 2D structure i.e. try to\n\
                     make it spread out\n\
     nFlipsPerSample - number of rotatable bonds that are\n\
                flipped at random at a time.\n\
     nSample - Number of random samplings of rotatable bonds.\n\
     sampleSeed - seed for the random sampling process.\n\
     permuteDeg4Nodes - allow permutation of bonds at a degree 4\n\
                 node during the sampling process \n\
     bondLength - change the default bond length for depiction \n\
     forceRDKit - use RDKit to generate coordinates even if \n\
                  preferCoordGen is set to true\n\n\
  RETURNS: \n\n\
     ID of the conformation added to the molecule\n";
  python::def(
      "Compute2DCoordsMimicDistmat", RDDepict::Compute2DCoordsMimicDistmat,
      (python::arg("mol"), python::arg("distMat"),
       python::arg("canonOrient") = false, python::arg("clearConfs") = true,
       python::arg("weightDistMat") = 0.5, python::arg("nFlipsPerSample") = 3,
       python::arg("nSample") = 100, python::arg("sampleSeed") = 100,
       python::arg("permuteDeg4Nodes") = true, python::arg("bondLength") = -1.0,
       python::arg("forceRDKit") = false),
      docString.c_str());

  docString =
      "Generate a depiction for a molecule where a piece of the \n\
  molecule is constrained to have the same coordinates as a reference. \n\n\
  The constraint can be hard (default) or soft. \n\
\n\
  Hard (default, ConstrainedDepictionParams.alignOnly=False): \n\
  Existing molecule coordinates, if present, are discarded; \n\
  new coordinates are generated constraining a piece of the molecule \n\
  to have the same coordinates as the reference, while the rest of \n\
  the molecule is built around it. \n\
  If ConstrainedDepictionParams.adjustMolBlockWedging is False \n\
  (default), existing molblock wedging information is always preserved. \n\
  If ConstrainedDepictionParams.adjustMolBlockWedging is True, \n\
  existing molblock wedging information is preserved in case it \n\
  only involves the invariant core and the core conformation has not \n\
  changed, while it is cleared in case the wedging is also outside \n\
  the invariant core, or core coordinates were changed. \n\
  If ConstrainedDepictionParams.acceptFailure is set to True and no \n\
  substructure match is found, coordinates will be recomputed from \n\
  scratch, hence molblock wedging information will be cleared. \n\
 \n\
  Soft (ConstrainedDepictionParams.alignOnly=True): \n\
  Existing coordinates in the conformation identified by \n\
  ConstrainedDepictionParams.existingConfId are preserved if present, \n\
  otherwise unconstrained new coordinates are generated. \n\
  Subsequently, coodinates undergo a rigid-body alignment to the reference. \n\
  If ConstrainedDepictionParams.adjustMolBlockWedging is False \n\
  (default), existing molblock wedging information is always preserved. \n\
  If ConstrainedDepictionParams.adjustMolBlockWedging is True, \n\
  existing molblock wedging information is inverted in case the rigid-body \n\
  alignment involved a flip around the Z axis. \n\
 \n\
  This is useful, for example, for generating depictions of SAR data \n\
  sets such that the cores of the molecules are all oriented the same way. \n\
  ARGUMENTS: \n\n\
  mol -    the molecule to be aligned, this will come back \n\
           with a single conformer. \n\
  reference -    a molecule with the reference atoms to align to; \n\
                 this should have a depiction. \n\
  confId -       (optional) the id of the reference conformation to use \n\
  refPatt -      (optional) a query molecule to be used to generate \n\
                 the atom mapping between the molecule and the reference \n\
  params - (optional) a ConstrainedDepictionParams instance\n\n\
  RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom \n\
           indices in mol constrained to have the same coordinates as atom \n\
           indices in reference.\n";
  python::def(
      "GenerateDepictionMatching2DStructure",
      static_cast<python::tuple (*)(RDKit::ROMol &, const RDKit::ROMol &,
                                    int confId, const python::object &,
                                    const python::object &)>(
          GenerateDepictionMatching2DStructure),
      (python::arg("mol"), python::arg("reference"), python::arg("confId") = -1,
       python::arg("refPatt") = python::object(),
       python::arg("params") = python::object()),
      docString.c_str());

  docString =
      "Generate a depiction for a molecule where a piece of the \n\
  molecule is constrained to have the same coordinates as a reference. \n\n\
  This is useful, for example, for generating depictions of SAR data \n\
  sets such that the cores of the molecules are all oriented the same way. \n\
  ARGUMENTS: \n\n\
  mol -    the molecule to be aligned, this will come back \n\
           with a single conformer. \n\
  reference -    a molecule with the reference atoms to align to; \n\
                 this should have a depiction. \n\
  confId -       the id of the reference conformation to use \n\
  refPatt -      a query molecule to be used to generate \n\
                 the atom mapping between the molecule and the reference \n\
  acceptFailure - if True, standard depictions will be generated \n\
                  for molecules that don't have a substructure match to the \n\
                  reference; if False, throws a DepictException.\n\
  forceRDKit -    (optional) use RDKit to generate coordinates even if \n\
                  preferCoordGen is set to true\n\
  allowRGroups -  (optional) if True, terminal dummy atoms in the \n\
                  reference are ignored if they match an implicit \n\
                  hydrogen in the molecule, and a constrained \n\
                  depiction is still attempted\n\n\
  RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom \n\
           indices in mol constrained to have the same coordinates as atom \n\
           indices in reference.\n";
  python::def(
      "GenerateDepictionMatching2DStructure",
      static_cast<python::tuple (*)(RDKit::ROMol &, const RDKit::ROMol &,
                                    int confId, const python::object &, bool,
                                    bool, bool)>(
          GenerateDepictionMatching2DStructure),
      (python::arg("mol"), python::arg("reference"), python::arg("confId") = -1,
       python::arg("refPatt") = python::object(),
       python::arg("acceptFailure") = false, python::arg("forceRDKit") = false,
       python::arg("allowRGroups") = false),
      docString.c_str());

  docString =
      "Generate a depiction for a molecule where a piece of the \n\
  molecule is constrained to have the same coordinates as a reference. \n\n\
  This is useful for, for example, generating depictions of SAR data \n\
  sets so that the cores of the molecules are all oriented the same way. \n\
  ARGUMENTS: \n\n\
  mol -    the molecule to be aligned, this will come back \n\
           with a single conformer. \n\
  reference -    a molecule with the reference atoms to align to; \n\
                 this should have a depiction. \n\
  atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will \n\
                 be used to generate the atom mapping between the molecule \n\
                 and the reference. Note that this sequence can be shorter \n\
                 than the number of atoms in the reference.\n\
  confId -       (optional) the id of the reference conformation to use \n\
  params -       (optional) an instance of ConstrainedDepictionParams\n";
  python::def(
      "GenerateDepictionMatching2DStructure",
      static_cast<void (*)(RDKit::ROMol &, const RDKit::ROMol &,
                           const python::object &, int,
                           const python::object &)>(
          RDDepict::GenerateDepictionMatching2DStructureAtomMap),
      (python::arg("mol"), python::arg("reference"), python::arg("atomMap"),
       python::arg("confId") = -1, python::arg("params") = python::object()),
      docString.c_str());

  docString =
      "Generate a depiction for a molecule where a piece of the \n\
  molecule is constrained to have the same coordinates as a reference. \n\n\
  This is useful for, for example, generating depictions of SAR data \n\
  sets so that the cores of the molecules are all oriented the same way. \n\
  ARGUMENTS: \n\n\
  mol -    the molecule to be aligned, this will come back \n\
           with a single conformer. \n\
  reference -    a molecule with the reference atoms to align to; \n\
                 this should have a depiction. \n\
  atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will \n\
                 be used to generate the atom mapping between the molecule \n\
                 and the reference. Note that this sequence can be shorter \n\
                 than the number of atoms in the reference.\n\
  confId -       the id of the reference conformation to use \n\
  forceRDKit -   use RDKit to generate coordinates even if \n\
                 preferCoordGen is set to true\n";
  python::def(
      "GenerateDepictionMatching2DStructure",
      static_cast<void (*)(RDKit::ROMol &, const RDKit::ROMol &,
                           const python::object &, int, bool)>(
          RDDepict::GenerateDepictionMatching2DStructureAtomMap),
      (python::arg("mol"), python::arg("reference"), python::arg("atomMap"),
       python::arg("confId"), python::arg("forceRDKit")),
      docString.c_str());

  docString =
      "Generate a depiction for a molecule where a piece of the molecule \n\
  is constrained to have coordinates similar to those of a 3D reference \n\
  structure.\n\
  ARGUMENTS: \n\n\
  mol -    the molecule to be aligned, this will come back \n\
           with a single conformer containing the 2D coordinates. \n\
  reference -    a molecule with the reference atoms to align to. \n\
                 By default this should be the same as mol, but with \n\
                 3D coordinates \n\
  confId -       (optional) the id of the reference conformation to use \n\
  referencePattern -  (optional) a query molecule to map a subset of \n\
                      the reference onto the mol, so that only some of the \n\
                      atoms are aligned. \n\
  acceptFailure - (optional) if True, standard depictions will be generated \n\
                  for molecules that don't match the reference or the\n\
                  referencePattern; if False, throws a DepictException.\n\
  forceRDKit -    (optional) use RDKit to generate coordinates even if \n\
                  preferCoordGen is set to true";

  python::def(
      "GenerateDepictionMatching3DStructure",
      RDDepict::GenerateDepictionMatching3DStructure,
      (python::arg("mol"), python::arg("reference"), python::arg("confId") = -1,
       python::arg("refPatt") = python::object(),
       python::arg("acceptFailure") = false, python::arg("forceRDKit") = false),
      docString.c_str());

  docString =
      "Rotate the 2D depiction such that the majority of bonds have a\n\
  30-degree angle with the X axis.\n\
  ARGUMENTS:\n\n\
  mol              - the molecule to be rotated.\n\
  confId           - (optional) the id of the reference conformation to use.\n\
  minimizeRotation - (optional) if False (the default), the molecule\n\
                     is rotated such that the majority of bonds have an angle\n\
                     with the X axis of 30 or 90 degrees. If True, the minimum\n\
                     rotation is applied such that the majority of bonds have\n\
                     an angle with the X axis of 0, 30, 60, or 90 degrees,\n\
                     with the goal of altering the initial orientation as\n\
                     little as possible .";

  python::def("StraightenDepiction", RDDepict::straightenDepiction,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("minimizeRotation") = false),
              docString.c_str());

  docString =
      "Normalizes the 2D depiction.\n\
If canonicalize is != 0, the depiction is subjected to a canonical\n\
transformation such that its main axis is aligned along the X axis\n\
(canonicalize >0, the default) or the Y axis (canonicalize <0).\n\
If canonicalize is 0, no canonicalization takes place.\n\
If scaleFactor is <0.0 (the default) the depiction is scaled such\n\
that bond lengths conform to RDKit standards. The applied scaling\n\
factor is returned.\n\n\
ARGUMENTS:\n\n\
mol          - the molecule to be normalized\n\
confId       - (optional) the id of the reference conformation to use\n\
canonicalize - (optional) if != 0, a canonical transformation is\n\
               applied: if >0 (the default), the main molecule axis is\n\
               aligned to the X axis, if <0 to the Y axis.\n\
               If 0, no canonical transformation is applied.\n\
scaleFactor  - (optional) if >0.0, the scaling factor to apply. The default\n\
               (-1.0) means that the depiction is automatically scaled\n\
               such that bond lengths are the standard RDKit ones.\n\n\
RETURNS: the applied scaling factor.";

  python::def(
      "NormalizeDepiction", RDDepict::normalizeDepiction,
      (python::arg("mol"), python::arg("confId") = -1,
       python::arg("canonicalize") = 1, python::arg("scaleFactor") = -1.),
      docString.c_str());
}
