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
                             double bondLength = -1.0,
                             bool forceRDKit = false) {
  RDGeom::INT_POINT2D_MAP cMap;
  cMap.clear();
  python::list ks = coordMap.keys();
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(ks.attr("__len__")()); i++) {
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
  res = RDDepict::compute2DCoords(mol, &cMap, canonOrient, clearConfs,
                                  nFlipsPerSample, nSamples, sampleSeed,
                                  permuteDeg4Nodes, forceRDKit);
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
    RDKit::ROMol &mol, RDKit::ROMol &reference, int confId,
    python::object refPatt, bool acceptFailure, bool forceRDKit,
    bool allowRGroups) {
  RDKit::ROMol *referencePattern = nullptr;
  if (refPatt != python::object()) {
    referencePattern = python::extract<RDKit::ROMol *>(refPatt);
  }
  auto matchVect = RDDepict::generateDepictionMatching2DStructure(
      mol, reference, confId, referencePattern, acceptFailure, forceRDKit,
      allowRGroups);
  python::list atomMap;
  for (const auto &pair : matchVect) {
    atomMap.append(python::make_tuple(pair.first, pair.second));
  }
  return python::tuple(atomMap);
}

void GenerateDepictionMatching2DStructureAtomMap(RDKit::ROMol &mol,
                                                 RDKit::ROMol &reference,
                                                 python::object atomMap,
                                                 int confId, bool forceRDKit) {
  std::unique_ptr<RDKit::MatchVectType> matchVect(translateAtomMap(atomMap));

  RDDepict::generateDepictionMatching2DStructure(mol, reference, *matchVect,
                                                 confId, forceRDKit);
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
void setPreferCoordGen(bool value) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = value;
#endif
}
}  // namespace RDDepict

BOOST_PYTHON_MODULE(rdDepictor) {
  python::scope().attr("__doc__") =
      "Module containing the functionality to compute 2D coordinates for a "
      "molecule";
  python::register_exception_translator<RDDepict::DepictException>(
      &rdDepictExceptionTranslator);

  rdkit_import_array();

  python::def("SetPreferCoordGen", setPreferCoordGen, python::arg("val"),
#ifdef RDK_BUILD_COORDGEN_SUPPORT
              "Sets whether or not the CoordGen library should be preferred to "
              "the RDKit depiction library."
#else
              "Has no effect (CoordGen support not enabled)"
#endif
  );
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
                  preferCoordGen is set to true\n\n\
  RETURNS: \n\n\
     ID of the conformation added to the molecule\n";
  python::def(
      "Compute2DCoords", RDDepict::Compute2DCoords,
      (python::arg("mol"), python::arg("canonOrient") = true,
       python::arg("clearConfs") = true,
       python::arg("coordMap") = python::dict(),
       python::arg("nFlipsPerSample") = 0, python::arg("nSample") = 0,
       python::arg("sampleSeed") = 0, python::arg("permuteDeg4Nodes") = false,
       python::arg("bondLength") = -1.0, python::arg("forceRDKit") = false),
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
  This is useful for, for example, generating depictions of SAR data \n\
  sets so that the cores of the molecules are all oriented the same way. \n\
  ARGUMENTS: \n\n\
  mol -    the molecule to be aligned, this will come back \n\
           with a single conformer. \n\
  reference -    a molecule with the reference atoms to align to; \n\
                 this should have a depiction. \n\
  confId -       (optional) the id of the reference conformation to use \n\
  refPatt -      (optional) a query molecule to be used to generate \n\
                 the atom mapping between the molecule and the reference \n\
  acceptFailure - (optional) if True, standard depictions will be generated \n\
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
      RDDepict::GenerateDepictionMatching2DStructure,
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
                 and the reference. \n\
  confId -       (optional) the id of the reference conformation to use \n\
  forceRDKit -   (optional) use RDKit to generate coordinates even if \n\
                 preferCoordGen is set to true\n";
  python::def(
      "GenerateDepictionMatching2DStructure",
      RDDepict::GenerateDepictionMatching2DStructureAtomMap,
      (python::arg("mol"), python::arg("reference"), python::arg("atomMap"),
       python::arg("confId") = -1, python::arg("forceRDKit") = false),
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
