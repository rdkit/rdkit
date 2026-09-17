//
//  Copyright (C) 2003-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Depictor/EmbeddedFrag.h>
#include <GraphMol/Depictor/DepictUtils.h>
#include <Geometry/point.h>
#include <RDBoost/Wrap_nb.h>

using namespace RDDepict;

namespace nb = nanobind;
using namespace nb::literals;

static std::unique_ptr<RDKit::MatchVectType> translateAtomMap(
    const nb::object &atomMap) {
  if (atomMap.is_none()) {
    return nullptr;
  }
  auto res = std::make_unique<RDKit::MatchVectType>();
  for (auto item : atomMap) {
    auto pair = nb::cast<nb::sequence>(item);
    int v1 = nb::cast<int>(pair[0]);
    int v2 = nb::cast<int>(pair[1]);
    res->push_back(std::make_pair(v1, v2));
  }
  return res;
}

static unsigned int compute2DCoordsHelper(RDKit::ROMol &mol, bool canonOrient,
                                          bool clearConfs, nb::dict coordMap,
                                          unsigned int nFlipsPerSample,
                                          unsigned int nSamples, int sampleSeed,
                                          int permuteDeg4Nodes, double bondLength,
                                          bool forceRDKit, bool useRingTemplates) {
  RDGeom::INT_POINT2D_MAP cMap;
  cMap.clear();
  for (auto item : coordMap) {
    unsigned int id = nb::cast<unsigned int>(item.first);
    if (id >= mol.getNumAtoms()) {
      throw ValueErrorException("atom index out of range");
    }
    cMap[id] = nb::cast<RDGeom::Point2D>(item.second);
  }
  double oBondLen = RDDepict::BOND_LEN;
  if (bondLength > 0) {
    RDDepict::BOND_LEN = bondLength;
  }
  unsigned int res =
      RDDepict::compute2DCoords(mol, &cMap, canonOrient, clearConfs,
                                nFlipsPerSample, nSamples, sampleSeed,
                                (bool)permuteDeg4Nodes, forceRDKit, useRingTemplates);
  if (bondLength > 0) {
    RDDepict::BOND_LEN = oBondLen;
  }
  return res;
}

static unsigned int compute2DCoordsMimicDistmatHelper(
    RDKit::ROMol &mol,
    nb::ndarray<nb::numpy, double, nb::ndim<1>> distMat, bool canonOrient,
    bool clearConfs, double weightDistMat, unsigned int nFlipsPerSample,
    unsigned int nSamples, int sampleSeed, bool permuteDeg4Nodes,
    double bondLength, bool forceRDKit) {
  unsigned int na = mol.getNumAtoms();
  unsigned int nitems = na * (na - 1) / 2;

  if (distMat.size() != nitems) {
    throw ValueErrorException(
        "The array size does not match the number of atoms in the molecule");
  }
  auto *inData = distMat.data();
  auto *cData = new double[nitems];
  memcpy(static_cast<void *>(cData), static_cast<const void *>(inData),
         nitems * sizeof(double));

  DOUBLE_SMART_PTR dmat(cData);
  double oBondLen = RDDepict::BOND_LEN;
  if (bondLength > 0) {
    RDDepict::BOND_LEN = bondLength;
  }
  unsigned int res = RDDepict::compute2DCoordsMimicDistMat(
      mol, &dmat, canonOrient, clearConfs, weightDistMat, nFlipsPerSample,
      nSamples, sampleSeed, permuteDeg4Nodes, forceRDKit);
  if (bondLength > 0) {
    RDDepict::BOND_LEN = oBondLen;
  }
  return res;
}

static nb::tuple generate2DStructureHelper(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
    const nb::object &refPatt, const ConstrainedDepictionParams &params) {
  RDKit::ROMol *referencePattern = nullptr;
  if (!refPatt.is_none()) {
    referencePattern = nb::cast<RDKit::ROMol *>(refPatt);
  }
  auto matchVect = RDDepict::generateDepictionMatching2DStructure(
      mol, reference, confId, referencePattern, params);
  nb::list atomMap;
  for (const auto &pair : matchVect) {
    atomMap.append(nb::make_tuple(pair.first, pair.second));
  }
  return nb::tuple(atomMap);
}

static nb::tuple generate2DStructureWithParamsHelper(
    RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
    const nb::object &refPatt, const nb::object &pyParams) {
  ConstrainedDepictionParams params;
  if (!pyParams.is_none()) {
    params = nb::cast<ConstrainedDepictionParams>(pyParams);
  }
  return generate2DStructureHelper(mol, reference, confId, refPatt, params);
}

static void generate2DStructureAtomMapHelper(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const nb::object &atomMap, int confId,
    const ConstrainedDepictionParams &params) {
  std::unique_ptr<RDKit::MatchVectType> matchVect(translateAtomMap(atomMap));
  RDDepict::generateDepictionMatching2DStructure(mol, reference, *matchVect,
                                                 confId, params);
}

static void generate2DStructureAtomMapWithParamsHelper(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const nb::object &atomMap, int confId, const nb::object &pyParams) {
  ConstrainedDepictionParams params;
  if (!pyParams.is_none()) {
    params = nb::cast<ConstrainedDepictionParams>(pyParams);
  }
  generate2DStructureAtomMapHelper(mol, reference, atomMap, confId, params);
}

static void generate2DStructureAtomMapForceRDKitHelper(
    RDKit::ROMol &mol, const RDKit::ROMol &reference,
    const nb::object &atomMap, int confId, bool forceRDKit) {
  ConstrainedDepictionParams params;
  params.forceRDKit = forceRDKit;
  generate2DStructureAtomMapHelper(mol, reference, atomMap, confId, params);
}

static void generateDepictionMatching3DStructureHelper(RDKit::ROMol &mol,
                                                 RDKit::ROMol &reference,
                                                 int confId,
                                                 nb::object refPatt,
                                                 bool acceptFailure,
                                                 bool forceRDKit) {
  RDKit::ROMol *referencePattern = nullptr;
  if (!refPatt.is_none()) {
    referencePattern = nb::cast<RDKit::ROMol *>(refPatt);
  }
  RDDepict::generateDepictionMatching3DStructure(
      mol, reference, confId, referencePattern, acceptFailure, forceRDKit);
}

static bool isCoordGenSupportAvailable() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return true;
#else
  return false;
#endif
}

static void setPreferCoordGen(bool value) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = value;
#endif
}

static bool getPreferCoordGen() {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  return RDDepict::preferCoordGen;
#else
  return false;
#endif
}

class UsingCoordGen {
 public:
  UsingCoordGen() = delete;
  UsingCoordGen(bool temp_state)
      : m_initial_state{getPreferCoordGen()}, m_temp_state(temp_state) {}
  ~UsingCoordGen() = default;

  void enter() { setPreferCoordGen(m_temp_state); }

  void exit(nb::object /*exc_type*/, nb::object /*exc_val*/,
            nb::object /*traceback*/) {
    setPreferCoordGen(m_initial_state);
  }

 private:
  bool m_initial_state;
  bool m_temp_state;
};

NB_MODULE(rdDepictor, m) {
  m.doc() =
      R"DOC(Module containing the functionality to compute 2D coordinates for a molecule)DOC";

  nb::exception<RDDepict::DepictException>(m, "DepictException", PyExc_ValueError);

  m.def("IsCoordGenSupportAvailable", isCoordGenSupportAvailable,
        "Returns whether RDKit was built with CoordGen support.");

  m.def("SetPreferCoordGen", setPreferCoordGen, "val"_a,
#ifdef RDK_BUILD_COORDGEN_SUPPORT
        "Sets whether or not the CoordGen library should be preferred to "
        "the RDKit depiction library."
#else
        "Has no effect (CoordGen support not enabled)"
#endif
  );

  m.def(
      "SetRingSystemTemplates", RDDepict::setRingSystemTemplates,
      "templatePath"_a,
      R"DOC(Loads the ring system templates from the specified file to be
used in 2D coordinate generation. Each template must be a single
line in the file represented using CXSMILES, and the structure should
be a single ring system. Throws a DepictException if any templates
are invalid.)DOC");

  m.def(
      "AddRingSystemTemplates", RDDepict::addRingSystemTemplates,
      "templatePath"_a,
      R"DOC(Adds the ring system templates from the specified file to be
used in 2D coordinate generation. If there are duplicates, the most
recently added template will be used. Each template must be a single
line in the file represented using CXSMILES, and the structure should
be a single ring system. Throws a DepictException if any templates
are invalid.)DOC");

  m.def(
      "LoadDefaultRingSystemTemplates",
      RDDepict::loadDefaultRingSystemTemplates,
      "Loads the default ring system templates and removes existing ones, if present.");

  m.def("GetPreferCoordGen", getPreferCoordGen,
#ifdef RDK_BUILD_COORDGEN_SUPPORT
        "Return whether or not the CoordGen library is used for coordinate "
        "generation in the RDKit depiction library."
#else
        "Always returns False (CoordGen support not enabled)"
#endif
  );

  nb::class_<UsingCoordGen>(
      m, "UsingCoordGen",
      "Context manager to temporarily set CoordGen library preference in RDKit depiction.")
      .def(nb::init<bool>(), "temp_state"_a, "Constructor")
      .def("__enter__", &UsingCoordGen::enter)
      .def("__exit__", &UsingCoordGen::exit,
           "exc_type"_a = nb::none(), "exc_val"_a = nb::none(),
           "traceback"_a = nb::none());

  nb::class_<ConstrainedDepictionParams>(m, "ConstrainedDepictionParams",
                                         "Parameters controlling constrained depiction")
      .def(nb::init<>())
      .def_rw(
          "acceptFailure", &ConstrainedDepictionParams::acceptFailure,
          R"DOC(if False (default), a DepictException is thrown if the molecule
does not have a substructure match to the reference;
if True, an unconstrained depiction will be generated)DOC")
      .def_rw(
          "forceRDKit", &ConstrainedDepictionParams::forceRDKit,
          R"DOC(if True, use RDKit to generate coordinates even if preferCoordGen
is set to True; defaults to False)DOC")
      .def_rw(
          "allowRGroups", &ConstrainedDepictionParams::allowRGroups,
          R"DOC(if True, terminal dummy atoms in the reference are ignored
if they match an implicit hydrogen in the molecule or if they are
attached top a query atom; defaults to False)DOC")
      .def_rw(
          "alignOnly", &ConstrainedDepictionParams::alignOnly,
          R"DOC(if False (default), a part of the molecule is hard-constrained
to have the same coordinates as the reference, and the rest of
the molecule is built around it; if True, coordinates
from conformation existingConfId are preserved (if they exist)
or generated without constraints (if they do not exist), then
the conformation is rigid-body aligned to the reference)DOC")
      .def_rw("adjustMolBlockWedging",
               &ConstrainedDepictionParams::adjustMolBlockWedging,
               R"DOC(if True (default), existing wedging information
will be updated or cleared as required; if False,
existing molblock wedging information will
always be preserved)DOC")
      .def_rw(
          "existingConfId", &ConstrainedDepictionParams::existingConfId,
          R"DOC(conformation id whose 2D coordinates should be
rigid-body aligned to the reference (if alignOnly is True),
or used to determine whether existing molblock wedging information
can be preserved following the constrained depiction (if
adjustMolBlockWedging is True)DOC")
      .def_rw(
          "useRingTemplates", &ConstrainedDepictionParams::useRingTemplates,
          "use templates to generate coordinates of complex ring systems")
      .def("__setattr__", &safeSetattr);

  m.def(
      "Compute2DCoords", compute2DCoordsHelper,
      "mol"_a, "canonOrient"_a = true, "clearConfs"_a = true,
      "coordMap"_a = nb::dict(), "nFlipsPerSample"_a = 0, "nSample"_a = 0,
      "sampleSeed"_a = 0, "permuteDeg4Nodes"_a = 0, "bondLength"_a = -1.0,
      "forceRDKit"_a = false, "useRingTemplates"_a = false,
      R"DOC(Compute 2D coordinates for a molecule.
  The resulting coordinates are stored on each atom of the molecule

  ARGUMENTS:

     mol - the molecule of interest
     canonOrient - orient the molecule in a canonical way
     clearConfs - if true, all existing conformations on the molecule
             will be cleared
     coordMap - a dictionary mapping atom Ids -> Point2D objects
                with starting coordinates for atoms that should
                have their positions locked.
     nFlipsPerSample - number of rotatable bonds that are
                flipped at random at a time.
     nSample - Number of random samplings of rotatable bonds.
     sampleSeed - seed for the random sampling process.
     permuteDeg4Nodes - allow permutation of bonds at a degree 4
                 node during the sampling process
     bondLength - change the default bond length for depiction
     forceRDKit - use RDKit to generate coordinates even if
                  preferCoordGen is set to true
     useRingTemplates - use templates to generate coordinates of complex
                  ring systems

  RETURNS:

     ID of the conformation added to the molecule)DOC");

  m.def(
      "Compute2DCoordsMimicDistmat", compute2DCoordsMimicDistmatHelper,
      "mol"_a, "distMat"_a, "canonOrient"_a = false, "clearConfs"_a = true,
      "weightDistMat"_a = 0.5, "nFlipsPerSample"_a = 3, "nSample"_a = 100,
      "sampleSeed"_a = 100, "permuteDeg4Nodes"_a = true, "bondLength"_a = -1.0,
      "forceRDKit"_a = false,
      R"DOC(Compute 2D coordinates for a molecule such
  that the inter-atom distances mimic those in a user-provided
  distance matrix.
  The resulting coordinates are stored on each atom of the molecule

  ARGUMENTS:

     mol - the molecule of interest
     distMat - distance matrix that we want the 2D structure to mimic
     canonOrient - orient the molecule in a canonical way
     clearConfs - if true, all existing conformations on the molecule
             will be cleared
     weightDistMat - weight assigned in the cost function to mimicking
                     the distance matrix.
                     This must be between (0.0,1.0). (1.0-weightDistMat)
                     is then the weight assigned to improving
                     the density of the 2D structure i.e. try to
                     make it spread out
     nFlipsPerSample - number of rotatable bonds that are
                flipped at random at a time.
     nSample - Number of random samplings of rotatable bonds.
     sampleSeed - seed for the random sampling process.
     permuteDeg4Nodes - allow permutation of bonds at a degree 4
                 node during the sampling process
     bondLength - change the default bond length for depiction
     forceRDKit - use RDKit to generate coordinates even if
                  preferCoordGen is set to true

  RETURNS:

     ID of the conformation added to the molecule)DOC");

  m.def(
      "GenerateDepictionMatching2DStructure",
      generate2DStructureWithParamsHelper,
      "mol"_a, "reference"_a, "confId"_a = -1, "refPatt"_a = nb::none(),
      "params"_a = nb::none(),
      R"DOC(Generate a depiction for a molecule where a piece of the
  molecule is constrained to have the same coordinates as a reference.

  The constraint can be hard (default) or soft.

  Hard (default, ConstrainedDepictionParams.alignOnly=False):
  Existing molecule coordinates, if present, are discarded;
  new coordinates are generated constraining a piece of the molecule
  to have the same coordinates as the reference, while the rest of
  the molecule is built around it.
  If ConstrainedDepictionParams.adjustMolBlockWedging is False
  (default), existing molblock wedging information is always preserved.
  If ConstrainedDepictionParams.adjustMolBlockWedging is True,
  existing molblock wedging information is preserved in case it
  only involves the invariant core and the core conformation has not
  changed, while it is cleared in case the wedging is also outside
  the invariant core, or core coordinates were changed.
  If ConstrainedDepictionParams.acceptFailure is set to True and no
  substructure match is found, coordinates will be recomputed from
  scratch, hence molblock wedging information will be cleared.

  Soft (ConstrainedDepictionParams.alignOnly=True):
  Existing coordinates in the conformation identified by
  ConstrainedDepictionParams.existingConfId are preserved if present,
  otherwise unconstrained new coordinates are generated.
  Subsequently, coodinates undergo a rigid-body alignment to the reference.
  If ConstrainedDepictionParams.adjustMolBlockWedging is False
  (default), existing molblock wedging information is always preserved.
  If ConstrainedDepictionParams.adjustMolBlockWedging is True,
  existing molblock wedging information is inverted in case the rigid-body
  alignment involved a flip around the Z axis.

  This is useful, for example, for generating depictions of SAR data
  sets such that the cores of the molecules are all oriented the same way.
  ARGUMENTS:

  mol -    the molecule to be aligned, this will come back
           with a single conformer.
  reference -    a molecule with the reference atoms to align to;
                 this should have a depiction.
  confId -       (optional) the id of the reference conformation to use
  refPatt -      (optional) a query molecule to be used to generate
                 the atom mapping between the molecule and the reference
  params - (optional) a ConstrainedDepictionParams instance

  RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom
           indices in mol constrained to have the same coordinates as atom
           indices in reference.)DOC");

  m.def(
      "GenerateDepictionMatching2DStructure",
      [](RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
         const nb::object &refPatt, bool acceptFailure, bool forceRDKit,
         bool allowRGroups) {
        ConstrainedDepictionParams params;
        params.acceptFailure = acceptFailure;
        params.forceRDKit = forceRDKit;
        params.allowRGroups = allowRGroups;
        return generate2DStructureHelper(mol, reference, confId, refPatt,
                                        params);
      },
      "mol"_a, "reference"_a, "confId"_a = -1, "refPatt"_a = nb::none(),
      "acceptFailure"_a = false, "forceRDKit"_a = false,
      "allowRGroups"_a = false,
      R"DOC(Generate a depiction for a molecule where a piece of the
  molecule is constrained to have the same coordinates as a reference.

  This is useful, for example, for generating depictions of SAR data
  sets such that the cores of the molecules are all oriented the same way.
  ARGUMENTS:

  mol -    the molecule to be aligned, this will come back
           with a single conformer.
  reference -    a molecule with the reference atoms to align to;
                 this should have a depiction.
  confId -       the id of the reference conformation to use
  refPatt -      a query molecule to be used to generate
                 the atom mapping between the molecule and the reference
  acceptFailure - if True, standard depictions will be generated
                  for molecules that don't have a substructure match to the
                  reference; if False, throws a DepictException.
  forceRDKit -    (optional) use RDKit to generate coordinates even if
                  preferCoordGen is set to true
  allowRGroups -  (optional) if True, terminal dummy atoms in the
                  reference are ignored if they match an implicit
                  hydrogen in the molecule, and a constrained
                  depiction is still attempted

  RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom
           indices in mol constrained to have the same coordinates as atom
           indices in reference.)DOC");

  m.def(
      "GenerateDepictionMatching2DStructure",
      generate2DStructureAtomMapWithParamsHelper,
      "mol"_a, "reference"_a, "atomMap"_a, "confId"_a = -1,
      "params"_a = nb::none(),
      R"DOC(Generate a depiction for a molecule where a piece of the
  molecule is constrained to have the same coordinates as a reference.

  This is useful for, for example, generating depictions of SAR data
  sets so that the cores of the molecules are all oriented the same way.
  ARGUMENTS:

  mol -    the molecule to be aligned, this will come back
           with a single conformer.
  reference -    a molecule with the reference atoms to align to;
                 this should have a depiction.
  atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will
                 be used to generate the atom mapping between the molecule
                 and the reference. Note that this sequence can be shorter
                 than the number of atoms in the reference.
  confId -       (optional) the id of the reference conformation to use
  params -       (optional) an instance of ConstrainedDepictionParams)DOC");

  m.def(
      "GenerateDepictionMatching2DStructure",
      generate2DStructureAtomMapForceRDKitHelper,
      "mol"_a, "reference"_a, "atomMap"_a, "confId"_a, "forceRDKit"_a,
      R"DOC(Generate a depiction for a molecule where a piece of the
  molecule is constrained to have the same coordinates as a reference.

  This is useful for, for example, generating depictions of SAR data
  sets so that the cores of the molecules are all oriented the same way.
  ARGUMENTS:

  mol -    the molecule to be aligned, this will come back
           with a single conformer.
  reference -    a molecule with the reference atoms to align to;
                 this should have a depiction.
  atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will
                 be used to generate the atom mapping between the molecule
                 and the reference. Note that this sequence can be shorter
                 than the number of atoms in the reference.
  confId -       the id of the reference conformation to use
  forceRDKit -   use RDKit to generate coordinates even if
                 preferCoordGen is set to true)DOC");

  m.def(
      "GenerateDepictionMatching3DStructure",
      generateDepictionMatching3DStructureHelper,
      "mol"_a, "reference"_a, "confId"_a = -1, "refPatt"_a = nb::none(),
      "acceptFailure"_a = false, "forceRDKit"_a = false,
      R"DOC(Generate a depiction for a molecule where a piece of the molecule
  is constrained to have coordinates similar to those of a 3D reference
  structure.
  ARGUMENTS:

  mol -    the molecule to be aligned, this will come back
           with a single conformer containing the 2D coordinates.
  reference -    a molecule with the reference atoms to align to.
                 By default this should be the same as mol, but with
                 3D coordinates
  confId -       (optional) the id of the reference conformation to use
  referencePattern -  (optional) a query molecule to map a subset of
                      the reference onto the mol, so that only some of the
                      atoms are aligned.
  acceptFailure - (optional) if True, standard depictions will be generated
                  for molecules that don't match the reference or the
                  referencePattern; if False, throws a DepictException.
  forceRDKit -    (optional) use RDKit to generate coordinates even if
                  preferCoordGen is set to true)DOC");

  m.def(
      "StraightenDepiction", RDDepict::straightenDepiction,
      "mol"_a, "confId"_a = -1, "minimizeRotation"_a = false,
      R"DOC(Rotate the 2D depiction such that the majority of bonds have a
  30-degree angle with the X axis.
  ARGUMENTS:

  mol              - the molecule to be rotated.
  confId           - (optional) the id of the reference conformation to use.
  minimizeRotation - (optional) if False (the default), the molecule
                     is rotated such that the majority of bonds have an angle
                     with the X axis of 30 or 90 degrees. If True, the minimum
                     rotation is applied such that the majority of bonds have
                     an angle with the X axis of 0, 30, 60, or 90 degrees,
                     with the goal of altering the initial orientation as
                     little as possible .)DOC");

  m.def(
      "NormalizeDepiction", RDDepict::normalizeDepiction,
      "mol"_a, "confId"_a = -1, "canonicalize"_a = 1, "scaleFactor"_a = -1.,
      R"DOC(Normalizes the 2D depiction.
If canonicalize is != 0, the depiction is subjected to a canonical
transformation such that its main axis is aligned along the X axis
(canonicalize >0, the default) or the Y axis (canonicalize <0).
If canonicalize is 0, no canonicalization takes place.
If scaleFactor is <0.0 (the default) the depiction is scaled such
that bond lengths conform to RDKit standards. The applied scaling
factor is returned.

ARGUMENTS:

mol          - the molecule to be normalized
confId       - (optional) the id of the reference conformation to use
canonicalize - (optional) if != 0, a canonical transformation is
               applied: if >0 (the default), the main molecule axis is
               aligned to the X axis, if <0 to the Y axis.
               If 0, no canonical transformation is applied.
scaleFactor  - (optional) if >0.0, the scaling factor to apply. The default
               (-1.0) means that the depiction is automatically scaled
               such that bond lengths are the standard RDKit ones.

RETURNS: the applied scaling factor.)DOC");
}
