//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <string>

#include <boost/python.hpp>

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>
#include <RDBoost/Wrap.h>

namespace python = boost::python;

namespace RDKit {

namespace helpers {
void set_customFeatures(GaussianShape::ShapeInputOptions &shp,
                        const python::object &s) {
  shp.customFeatures.clear();
  auto numVecs = python::len(s);
  shp.customFeatures.reserve(numVecs);
  for (auto i = 0u; i < numVecs; ++i) {
    const auto outVec = s[i];
    auto numFeats = python::len(outVec);
    std::vector<GaussianShape::CustomFeature> feats;
    feats.reserve(numFeats);
    for (auto j = 0u; j < numFeats; ++j) {
      const auto feat = outVec[j];
      unsigned int featType = python::extract<unsigned int>(feat[0]);
      RDGeom::Point3D pos = python::extract<RDGeom::Point3D>(feat[1]);
      double radius = python::extract<double>(feat[2]);
      std::vector<unsigned int> atoms;
      if (len(feat) == 4) {
        for (unsigned int k = 0; k < len(feat[3]); ++k) {
          atoms.push_back(python::extract<unsigned int>(feat[3][k]));
        }
      }
      feats.emplace_back(featType, pos, radius, atoms);
    }
    shp.customFeatures.emplace_back(std::move(feats));
  }
}

python::tuple get_customFeatures(const GaussianShape::ShapeInputOptions &shp) {
  python::list allFeatLists;
  for (const auto &feats : shp.customFeatures) {
    python::list featList;
    for (const auto &feat : feats) {
      python::list elem;
      elem.append(static_cast<int>(feat.type));
      elem.append(feat.pos);
      elem.append(feat.rad);
      elem.append(feat.atoms);
      featList.append(elem);
    }
    allFeatLists.append(featList);
  }
  return python::tuple(allFeatLists);
}

python::tuple alignMol1(const ROMol &ref, ROMol &fit,
                        const python::object &py_refOpts,
                        const python::object &py_fitOpts,
                        const python::object &py_overlayOpts, int refConfId,
                        int fitConfId) {
  GaussianShape::ShapeInputOptions refOpts, fitOpts;
  if (!py_refOpts.is_none()) {
    refOpts = python::extract<GaussianShape::ShapeInputOptions>(py_refOpts);
  }
  if (!py_fitOpts.is_none()) {
    fitOpts = python::extract<GaussianShape::ShapeInputOptions>(py_fitOpts);
  }
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  auto results = GaussianShape::AlignMolecule(
      ref, fit, refOpts, fitOpts, nullptr, overlayOpts, refConfId, fitConfId);
  return python::make_tuple(results[0], results[1], results[2]);
}

python::tuple alignMol2(const GaussianShape::ShapeInput &refShape, ROMol &fit,
                        const python::object &py_fitOpts,
                        const python::object &py_overlayOpts, int fitConfId) {
  GaussianShape::ShapeInputOptions fitOpts;
  if (!py_fitOpts.is_none()) {
    fitOpts = python::extract<GaussianShape::ShapeInputOptions>(py_fitOpts);
  }
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  auto results = GaussianShape::AlignMolecule(refShape, fit, fitOpts, nullptr,
                                              overlayOpts, fitConfId);
  return python::make_tuple(results[0], results[1], results[2]);
}

python::tuple alignShapes(const GaussianShape::ShapeInput &refShape,
                          GaussianShape::ShapeInput &fitShape,
                          const python::object &py_overlayOpts) {
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  RDGeom::Transform3D xform;
  auto results =
      GaussianShape::AlignShape(refShape, fitShape, &xform, overlayOpts);
  python::list pyMatrix;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      pyMatrix.append(xform.getValUnchecked(i, j));
    }
  }
  return python::make_tuple(results[0], results[1], results[2], pyMatrix);
}

python::tuple scoreMol1(const ROMol &ref, const ROMol &fit,
                        const python::object &py_refOpts,
                        const python::object &py_fitOpts,
                        const python::object &py_overlayOpts, int refConfId,
                        int fitConfId) {
  GaussianShape::ShapeInputOptions refOpts, fitOpts;
  if (!py_refOpts.is_none()) {
    refOpts = python::extract<GaussianShape::ShapeInputOptions>(py_refOpts);
  }
  if (!py_fitOpts.is_none()) {
    fitOpts = python::extract<GaussianShape::ShapeInputOptions>(py_fitOpts);
  }
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  std::array<double, 2> ovVols;
  auto results = GaussianShape::ScoreMolecule(
      ref, fit, refOpts, fitOpts, overlayOpts, refConfId, fitConfId, &ovVols);
  return python::make_tuple(results[0], results[1], results[2], ovVols[0],
                            ovVols[1]);
}

python::tuple scoreMol2(const GaussianShape::ShapeInput &refShape,
                        const ROMol &fit, const python::object &py_fitOpts,
                        const python::object &py_overlayOpts, int fitConfId) {
  GaussianShape::ShapeInputOptions fitOpts;
  if (!py_fitOpts.is_none()) {
    fitOpts = python::extract<GaussianShape::ShapeInputOptions>(py_fitOpts);
  }
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  std::array<double, 2> ovVols;
  auto results = GaussianShape::ScoreMolecule(refShape, fit, fitOpts,
                                              overlayOpts, fitConfId, &ovVols);
  return python::make_tuple(results[0], results[1], results[2], ovVols[0],
                            ovVols[1]);
}

python::tuple scoreShape(const GaussianShape::ShapeInput &refShape,
                         const GaussianShape::ShapeInput &fitShape,
                         const python::object &py_overlayOpts) {
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  std::array<double, 2> ovVols;
  auto results =
      GaussianShape::ScoreShape(refShape, fitShape, overlayOpts, &ovVols);
  return python::make_tuple(results[0], results[1], results[2], ovVols[0],
                            ovVols[1]);
}

void set_atomSubset(GaussianShape::ShapeInputOptions &opts,
                    const python::object &as) {
  pythonObjectToVect<unsigned int>(as, opts.atomSubset);
}

python::tuple get_atomSubset(const GaussianShape::ShapeInputOptions &opts) {
  python::list py_list;
  for (const auto &val : opts.atomSubset) {
    py_list.append(val);
  }
  return python::tuple(py_list);
}

void set_atomRadii(GaussianShape::ShapeInputOptions &opts,
                   const python::object &ar) {
  int len = python::len(ar);
  opts.atomRadii.resize(len);
  for (int i = 0; i < len; i++) {
    unsigned int atomIdx = python::extract<unsigned int>(ar[i][0]);
    double radius = python::extract<double>(ar[i][1]);
    opts.atomRadii[i] = std::make_pair(atomIdx, radius);
  }
}

python::tuple get_atomRadii(const GaussianShape::ShapeInputOptions &opts) {
  python::list py_list;
  for (const auto &val : opts.atomRadii) {
    py_list.append(python::make_tuple(static_cast<int>(val.first), val.second));
  }
  return python::tuple(py_list);
}

double getShapeVolume_helper(const GaussianShape::ShapeInput &shape) {
  return shape.getShapeVolume();
}

double getColorVolume_helper(const GaussianShape::ShapeInput &shape) {
  return shape.getColorVolume();
}

python::tuple bestSimilarity_helper(GaussianShape::ShapeInput &refShape,
                                    const GaussianShape::ShapeInput &fitShape,
                                    double threshold,
                                    const python::object &py_overlayOpts) {
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  unsigned int bestFitShape, bestThisShape;
  RDGeom::Transform3D bestXform;
  auto bestSim = refShape.bestSimilarity(fitShape, bestThisShape, bestFitShape,
                                         bestXform, threshold, overlayOpts);
  python::list results;
  results.append(python::make_tuple(bestSim[0], bestSim[1], bestSim[2]));
  results.append(bestThisShape);
  results.append(bestFitShape);
  python::list pyMatrix;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      pyMatrix.append(bestXform.getValUnchecked(i, j));
    }
  }
  results.append(pyMatrix);
  return python::tuple(results);
}

double maxPossibleSimilarity_helper(GaussianShape::ShapeInput &refShape,
                                    GaussianShape::ShapeInput &fitShape,
                                    const python::object &py_overlayOpts) {
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  return refShape.maxPossibleSimilarity(fitShape, overlayOpts);
}

ROMol *shapeToMol_helper(GaussianShape::ShapeInput &shape, bool includeColors,
                         bool withBonds) {
  auto mol = shape.shapeToMol(includeColors, withBonds);
  return static_cast<ROMol *>(mol.release());
}
}  // namespace helpers

void wrap_rdGaussianShape() {
  python::scope().attr("__doc__") =
      "Module containing implementation of Gaussian-based shape overlay and"
      " scoring."
      "NOTE: This functionality is experimental and the API"
      " and/or results may change in future releases.";

  python::enum_<RDKit::GaussianShape::StartMode>("StartMode")
      .value("ROTATE_0", RDKit::GaussianShape::StartMode::ROTATE_0)
      .value("ROTATE_180", RDKit::GaussianShape::StartMode::ROTATE_180)
      .value("ROTATE_180_WIGGLE",
             RDKit::GaussianShape::StartMode::ROTATE_180_WIGGLE)
      .value("ROTATE_45", RDKit::GaussianShape::StartMode::ROTATE_45)
      .value("ROTATE_0_FRAGMENT",
             RDKit::GaussianShape::StartMode::ROTATE_0_FRAGMENT)
      .value("ROTATE_180_FRAGMENT",
             RDKit::GaussianShape::StartMode::ROTATE_180_FRAGMENT)
      .value("ROTATE_45_FRAGMENT",
             RDKit::GaussianShape::StartMode::ROTATE_45_FRAGMENT)
      .value("A_LA_PUBCHEM", GaussianShape::StartMode::A_LA_PUBCHEM)
      .export_values();

  python::enum_<RDKit::GaussianShape::OptimMode>("OptimMode")
      .value("SHAPE_ONLY", RDKit::GaussianShape::OptimMode::SHAPE_ONLY)
      .value("SHAPE_PLUS_COLOR_SCORE",
             RDKit::GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE)
      .value("SHAPE_PLUS_COLOR",
             RDKit::GaussianShape::OptimMode::SHAPE_PLUS_COLOR)
      .export_values();

  python::class_<GaussianShape::ShapeInputOptions, boost::noncopyable>(
      "ShapeInputOptions",
      "ShapeInputOptions - options for setting up ShapeInput objects.")
      .def_readwrite("useColors", &GaussianShape::ShapeInputOptions::useColors,
                     "Whether to use color features in overlay.  Default=True.")
      .def_readwrite(
          "allCarbonRadii", &GaussianShape::ShapeInputOptions::allCarbonRadii,
          "Whether to use the same radius, appropriate for Carbon, for all atoms.  There is a"
          " slight accuracy penalty but significant speed gain if used.  Default=True.")
      .add_property(
          "atomSubset", &helpers::get_atomSubset, &helpers::set_atomSubset,
          "If not empty, use just these atoms in the molecule to form the ShapeInput object.")
      .add_property(
          "customFeatures", &helpers::get_customFeatures,
          &helpers::set_customFeatures,
          "Custom features for the shape.  Requires a list of lists of tuples of"
          " int (the feature type), Point3D (the coordinates), float (the radius)"
          " and optionally a list of indices of the atoms that the feature was derived from.")
      .add_property(
          "atomRadii", &helpers::get_atomRadii, &helpers::set_atomRadii,
          "Non-standard radii to use for the atoms specified by their indices"
          " in the molecule.  Not all atoms need have a radius specified."
          "  A list of tuples of [int, float].")
      .def("__setattr__", &safeSetattr);

  python::class_<GaussianShape::ShapeOverlayOptions, boost::noncopyable>(
      "ShapeOverlayOptions",
      "ShapeOverlayOptions - options for controlling the shape overlay process.")
      .def_readwrite(
          "startMode", &RDKit::GaussianShape::ShapeOverlayOptions::startMode,
          "Start modes for optimisation.  Default is A_LA_PUBCHEM - as used by the"
          " PubChem code - either ROTATE_180_WIGGLE or ROTATE_45 depending on the shape"
          " of the two molecules.  ROTATE_180_WIGGLE means 180 rotations about"
          " the x, y and z axes, then a small"
          " rotation about each axis from that point, using the best scoring one of"
          " those. ROTATE_180 uses 180 degree rotations for 4 start points,"
          " ROTATE_45 uses 45 degree rotations for 9 start points and ROTATE_0"
          " leaves the relative orientations of the 2 molecules as passed in before"
          " optimisation.  There are also ROTATE_0_FRAGMENT, ROTATE_45_FRAGMENT"
          " and ROTATE_180_FRAGMENT that as well as the above move the fit"
          " molecule to the ends of each of the principal axes and then does"
          " the appropriate rotations.  This is useful when the fit molecule is"
          " a lot smaller than the reference molecule, but requires a large number"
          " of optimisations so is relatively slow.")
      .def_readwrite(
          "optimMode", &GaussianShape::ShapeOverlayOptions::optimMode,
          "Optimisation mode, controlling what parameters are used"
          " to drive the overlay.  Default=SHAPE_PLUS_COLOR_SCORE which"
          " optimises using just the overlap of shape, but uses the"
          " color to decide which is the best overlay.  Other options"
          " are SHAPE_ONLY and SHAPE_AND_COLOR with the latter using"
          " the overlap of color features as well. ")
      .def_readwrite(
          "simAlpha", &GaussianShape::ShapeOverlayOptions::simAlpha,
          "When doing a Tversky similarity, the alpha value.  If alpha and"
          " beta are both the default 1.0, it's a Tanimoto similarity.  A"
          " high alpha and low beta emphasize the fit volume in the"
          " similarity and vice versa. Tversky is O / (A * (R - O) + B * (F"
          " - O) + O) where O is the overlap volume, R is the reference's"
          " volume and F is the fit's volume.  This is different from that"
          " used by OpenEye (O / (A * R + B * F)).")
      .def_readwrite("simBeta", &GaussianShape::ShapeOverlayOptions::simBeta,
                     "When doing a Tversky similarity, the beta value.")
      .def_readwrite(
          "optParam", &GaussianShape::ShapeOverlayOptions::optParam,
          "If using colors, the relative weights of the shape and color scores,"
          " as a fraction of 1.  Default=0.5.")
      .def_readwrite(
          "nSteps", &GaussianShape::ShapeOverlayOptions::nSteps,
          "Maximum number of steps for the shape overlay process. Default=100.")
      .def_readwrite(
          "normalize", &GaussianShape::ShapeOverlayOptions::normalize,
          "Whether to normalize the shapes before overlay by putting them into their"
          " canonical orientation (centred on the origin, aligned along its"
          " principal axes.  Default=True.")
      .def_readwrite(
          "useDistCutoff", &GaussianShape::ShapeOverlayOptions::useDistCutoff,
          "Whether to use distance cutoff when calculating the shape volumes.  If used,"
          " there will be a small penalty in accuracy but a significant increase in speed."
          "  Default=True.")
      .def_readwrite(
          "distCutoff", &GaussianShape::ShapeOverlayOptions::distCutoff,
          "If using a distance cutoff, this is the value used.  Default=4.5 of whatever"
          " units the coordinates are in.")
      .def_readwrite(
          "shapeConvergenceCriterion",
          &GaussianShape::ShapeOverlayOptions::shapeConvergenceCriterion,
          "Optimisation stops when the shape Tversky score changes by less"
          " than this amount after an optimisation step.  A larger number is"
          " faster but gives less precise overlays.  Default=0.001.")
      .def("__setattr__", &safeSetattr);

  std::string docString("ShapeInput object");
  python::class_<GaussianShape::ShapeInput, boost::noncopyable>(
      "ShapeInput", docString.c_str(),
      python::init<const ROMol &, int, const GaussianShape::ShapeInputOptions &,
                   const GaussianShape::ShapeOverlayOptions &>(
          python::args("self", "confId", "shapeOpt", "overlayOpts")))
      .add_property(
          "GetSmiles", &GaussianShape::ShapeInput::getSmiles,
          "Get the SMILES string for the molecule that the shape relates to.")
      .add_property(
          "setActiveShape", &GaussianShape::ShapeInput::setActiveShape,
          "Set the active shape, the one that will be used for overlays etc.")
      .add_property("getActiveShape",
                    &GaussianShape::ShapeInput::getActiveShape,
                    "Return the number of the active shape.")
      .add_property("NumAtoms", &GaussianShape::ShapeInput::getNumAtoms,
                    "Get the number of atoms defining the shape.")
      .add_property("NumFeatures", &GaussianShape::ShapeInput::getNumFeatures,
                    "Get the number of features in the shape.")
      .add_property(
          "NumShapes", &GaussianShape::ShapeInput::getNumShapes,
          "Get the number of shapes.  There will be a shape for each conformation "
          "of the input molecule, unless shape pruning was carried out in which case"
          " there may be fewer.")
      .add_property("ShapeVolume", &helpers::getShapeVolume_helper,
                    "Get the volume due to the atoms for the active shape.")
      .add_property(
          "ColorVolume", &helpers::getColorVolume_helper,
          "Get the volume of the shape's color features for the active shap.")
      .def(
          "NormalizeCoords", &GaussianShape::ShapeInput::normalizeCoords,
          "Align the principal axes to the cartesian axes and centre on the origin."
          " Doesn't require that the shape was created from a molecule.  Creates"
          " the necessary transformation if not already done.")
      .def(
          "ShapeToMol", &helpers::shapeToMol_helper,
          (python::arg("self"), python::arg("includeColors") = false,
           python::arg("withBonds") = true),
          "Return a molecule with coordinates of the current active shape."
          "  If includeColors is True, (default is False) the color features"
          " will be added as xenon atoms.  If withBonds is True (the default)"
          " a molecule with bonds will be created, if not then just atoms at the"
          " appropriate positions will be produced.",
          python::return_value_policy<python::manage_new_object>())
      .def(
          "BestSimilarity", &helpers::bestSimilarity_helper,
          (python::arg("self"), python::arg("fitShape"),
           python::arg("threshold") = -1.0,
           python::arg("overlayOpts") = python::object()),
          "Find the best similarity score between all shapes in this shape and the"
          " other one. Stops as soon as it gets something above the threshold."
          " The score runs between 0.0 and 1.0, so the default threshold of -1.0"
          " means no threshold. Fills in the shape numbers of the two that were"
          " responsible if there is something above the threshold, and the"
          " transformation that did it. Returns a tuple of the similarity scores"
          " ((-1.0, -1.0, -1.0) if there was nothing above the threshold), the number of the"
          " shape for this object and the shape number of the fitShape that gave"
          " the best similarity and the transformation matrix (as a list of 16 floats)"
          " that will reproduce the best overlay.  The shapes won't necessarily"
          " be left in the state that gave the best similarity.  Note that the"
          " shape numbers are not necessarily the same as the original molecule"
          " conformation numbers.")
      .def("MaxPossibleSimilarity", &helpers::maxPossibleSimilarity_helper,
           (python::arg("self"), python::arg("fitShape"),
            python::arg("overlayOpts") = python::object()),
           "Get the maximum possible similarity score between all shapes in"
           " this shape and all shapes in the fitShape.  The maximum similarity"
           " is when one shape is entirely inside the other.  This returns"
           " the similarity in that case, which is the upper bound on what"
           " is achievable between these 2 shapes.")
      .def("__setattr__", &safeSetattr);

  python::def(
      "AlignMol", &helpers::alignMol1,
      (python::arg("ref"), python::arg("fit"),
       python::arg("refOpts") = python::object(),
       python::arg("fitOpts") = python::object(),
       python::arg("overlayOpts") = python::object(),
       python::arg("refConfId") = -1, python::arg("fitConfId") = -1),
      R"DOC(Aligns a fit molecule onto a reference molecule.  The fit is modified.

Parameters
----------
ref: RDKit.ROMol
    Reference molecule
fit: RDKit.ROMol
    Fit molecule that will be overlaid
refOpts: ShapeInputOptions, optional
    Options for building the ref shape
fitOpts: ShapeInputOptions, optional
    Options for building the fit shape
overlayOpts: ShapeOverlayOptions, optional
    Options for controlling the overlay
refConfId : int, optional
    Reference conformer ID (default is -1)
fitConfId : int, optional
    fit conformer ID (default is -1)

Returns
-------
3-tuple of floats
    The results are (combo_score, shape_score, color_score).  The color_score is
    0.0 if color features not used, in which case combo_score and shape_score will
    be the same.
)DOC");

  python::def(
      "AlignMol", &helpers::alignMol2,
      (python::arg("refShape"), python::arg("fit"),
       python::arg("fitOpts") = python::object(),
       python::arg("overlayOpts") = python::object(),
       python::arg("fitConfId") = -1),
      R"DOC(Aligns a fit molecule onto a reference shape.  The fit is modified.

Parameters
----------
refShape: ShapeInput
    Reference shape
fit: RDKit.ROMol
    Fit molecule that will be overlaid
fitOpts: ShapeInputOptions, optional
    Options for building the fit shape
overlayOpts: ShapeOverlayOptions, optional
    Options for controlling the overlay
fitConfId : int, optional
    Fit conformer ID (default is -1)

Returns
-------
3-tuple of floats
    The results are (combo_score, shape_score, color_score).  The color_score is
    0.0 if color features not used, in which case combo_score and shape_score will
    be the same.)DOC");

  python::def(
      "AlignShapes", &helpers::alignShapes,
      (python::arg("refShape"), python::arg("fitShape"),
       python::arg("overlayOpts") = python::object()),
      R"DOC(Aligns a fit shape to a reference shape. The fit is modified.

Parameters
----------
refShape : ShapeInput
    Reference shape
fitShape : ShapeInput
    fit shape
overlayOpts: ShapeOverlayOptions, optional
    Options for controlling the overlay


Returns
-------
 4-tuple of float, float, list of floats
    The results are (combo_score, shape_score, color_score, matrix)
    The matrix is a 16-float list giving the transformation matrix that
    overlays the fit onto the reference.)DOC");

  python::def("ScoreMol", &helpers::scoreMol1,
              (python::arg("ref"), python::arg("fit"),
               python::arg("refOpts") = python::object(),
               python::arg("fitOpts") = python::object(),
               python::arg("overlayOpts") = python::object(),
               python::arg("refConfId") = -1, python::arg("fitConfId") = -1),
              R"DOC(Calculates the scores between a reference molecule and a fit
molecule without overlay.

Parameters
----------
ref: RDKit.ROMol
    Reference molecule
fit: RDKit.ROMol
    Fit molecule that will be scored
refOpts: ShapeInputOptions, optional
    Options for building the ref shape
fitOpts: ShapeInputOptions, optional
    Options for building the fit shape
overlayOpts: ShapeOverlayOptions, optional
    Options for controlling the volume calculation
refConfId : int, optional
    Reference conformer ID (default is -1)
fitConfId : int, optional
    fit conformer ID (default is -1)

Returns
-------
3-tuple of floats
    The results are (combo_score, shape_score, color_score).  The color_score is
    0.0 if color features not used, in which case combo_score and shape_score will
    be the same.
)DOC");

  python::def(
      "ScoreMol", &helpers::scoreMol2,
      (python::arg("refShape"), python::arg("fit"),
       python::arg("fitOpts") = python::object(),
       python::arg("overlayOpts") = python::object(),
       python::arg("fitConfId") = -1),
      R"DOC(Calculates the scores between a reference shape and a fit molecule
without overlay.

Parameters
----------
refShape: ShapeInput
    Reference shape
fit: RDKit.ROMol
    Fit molecule that will be scored
fitOpts: ShapeInputOptions, optional
    Options for building the fit shape
overlayOpts: ShapeOverlayOptions, optional
    Options for controlling the volume calculation
fitConfId : int, optional
    fit conformer ID (default is -1)

Returns
-------
3-tuple of floats
    The results are (combo_score, shape_score, color_score).  The color_score is
    0.0 if color features not used, in which case combo_score and shape_score will
    be the same.
)DOC");

  python::def(
      "ScoreShape", &helpers::scoreShape,
      (python::arg("refShape"), python::arg("fitShape"),
       python::arg("overlayOpts") = python::object()),
      R"DOC(Calculates the scores between a reference shape and a fit shape without
overlay.

Parameters
----------
refShape: ShapeInput
    Reference shape
fitShape: ShapeInput
    Fit shape
fitOpts: ShapeInputOptions, optional
    Options for building the fit shape
overlayOpts: ShapeOverlayOptions, optional
    Options for controlling the volume calculation

Returns
-------
3-tuple of floats
    The results are (combo_score, shape_score, color_score).  The color_score is
    0.0 if color features not used, in which case combo_score and shape_score will
    be the same.
)DOC");
}

BOOST_PYTHON_MODULE(rdGaussianShape) { wrap_rdGaussianShape(); }

}  // namespace RDKit