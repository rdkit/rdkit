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

#include "pubchem_shape/PubChemShape.hpp"

#include <string>

#include <boost/python.hpp>

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
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
  auto len = python::len(s);
  shp.customFeatures.reserve(len);
  for (auto i = 0u; i < len; ++i) {
    const auto elem = s[i];
    unsigned int featType = python::extract<unsigned int>(elem[0]);
    RDGeom::Point3D pos = python::extract<RDGeom::Point3D>(elem[1]);
    double radius = python::extract<double>(elem[2]);
    shp.customFeatures.emplace_back(featType, pos, radius);
  }
}
python::tuple get_customFeatures(const GaussianShape::ShapeInputOptions &shp) {
  python::list py_list;
  for (const auto &val : shp.customFeatures) {
    python::list elem;
    elem.append(static_cast<int>(std::get<0>(val)));
    elem.append(std::get<1>(val));
    elem.append(std::get<2>(val));
    py_list.append(python::tuple(elem));
  }
  return python::tuple(py_list);
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
  auto results = GaussianShape::ScoreMolecule(
      ref, fit, refOpts, fitOpts, overlayOpts, refConfId, fitConfId);
  return python::make_tuple(results[0], results[1], results[2]);
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
  auto results = GaussianShape::ScoreMolecule(refShape, fit, fitOpts,
                                              overlayOpts, fitConfId);
  return python::make_tuple(results[0], results[1], results[2]);
}

python::tuple scoreShape(const GaussianShape::ShapeInput &refShape,
                         const GaussianShape::ShapeInput &fitShape,
                         const python::object &py_overlayOpts) {
  GaussianShape::ShapeOverlayOptions overlayOpts;
  if (!py_overlayOpts.is_none()) {
    overlayOpts =
        python::extract<GaussianShape::ShapeOverlayOptions>(py_overlayOpts);
  }
  auto results = GaussianShape::ScoreShape(refShape, fitShape, overlayOpts);
  return python::make_tuple(results[0], results[1], results[2]);
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
      .value("ROTATE_90", RDKit::GaussianShape::StartMode::ROTATE_90)
      .value("ROTATE_0_FRAGMENT",
             RDKit::GaussianShape::StartMode::ROTATE_0_FRAGMENT)
      .value("ROTATE_180_FRAGMENT",
             RDKit::GaussianShape::StartMode::ROTATE_180_FRAGMENT)
      .value("ROTATE_90_FRAGMENT",
             RDKit::GaussianShape::StartMode::ROTATE_90_FRAGMENT)
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
          "customFeatures", &helpers::get_customFeatures,
          &helpers::set_customFeatures,
          "Custom features for the shape.  Requires a list of tuples of"
          " int (the feature type), Point3D (the coordinates) and float (the radius).")
      .def("__setattr__", &safeSetattr);

  python::class_<GaussianShape::ShapeOverlayOptions, boost::noncopyable>(
      "ShapeOverlayOptions",
      "ShapeOverlayOptions - options for controlling the shape overlay process.")
      .def_readwrite(
          "startMode", &RDKit::GaussianShape::ShapeOverlayOptions::startMode,
          "Start modes for optimisation.  Default is ROTATE_180_WIGGLE - is as used by the"
          " PubChem code - 180 rotations about the x, y and z axes, then a small"
          " rotation about each axis from that point, using the best scoring one of"
          " those. ROTATE_180 uses 180 degree rotations for 4 start points,"
          " ROTATE_90 uses 90 degree rotations for 9 start points and ROTATE_0"
          " leaves the relative orientations of the 2 molecules as passed in before"
          " optimisation.  There are also ROTATE_0_FRAGMENT, ROTATE_90_FRAGMENT"
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
          "there will be a small penalty in accuracy but a significant increase in speed."
          "  Default=True.")
      .def_readwrite(
          "distCutoff", &GaussianShape::ShapeOverlayOptions::distCutoff,
          "If using a distance cutoff, this is the value used.  Default=4.5 of whatever"
          " units the coordinates are in.")
      .def("__setattr__", &safeSetattr);

  std::string docString("ShapeInput object");
  python::class_<GaussianShape::ShapeInput, boost::noncopyable>(
      "ShapeInput", docString.c_str(),
      python::init<const ROMol &, int, const GaussianShape::ShapeInputOptions &,
                   const GaussianShape::ShapeOverlayOptions &>(
          python::args("self", "confId", "shapeOpt", "overlayOpts")))
      .add_property("NumAtoms", &GaussianShape::ShapeInput::getNumAtoms,
                    "Get the number of atoms defining the shape.")
      .add_property("NumFeatures", &GaussianShape::ShapeInput::getNumFeatures,
                    "Get the number of features in the shape.")
      .add_property("ShapeVolume", &GaussianShape::ShapeInput::getShapeVolume,
                    "Get the shape's volume due to the atoms.")
      .add_property("ColorVolume", &GaussianShape::ShapeInput::getColorVolume,
                    "Get the volume of the shape's color features.")
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
    0.0 if color features not used, in which case combo_score and shape_score should
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
    0.0 if color features not used, in which case combo_score and shape_score should
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
    0.0 if color features not used, in which case combo_score and shape_score should
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
    0.0 if color features not used, in which case combo_score and shape_score should
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
    0.0 if color features not used, in which case combo_score and shape_score should
    be the same.
)DOC");
}

BOOST_PYTHON_MODULE(rdGaussianShape) { wrap_rdGaussianShape(); }

}  // namespace RDKit