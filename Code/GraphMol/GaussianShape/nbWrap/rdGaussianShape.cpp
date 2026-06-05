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

#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>
#include <RDBoost/Wrap_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace helpers {

void set_customFeatures(GaussianShape::ShapeInputOptions &shp, nb::object s) {
  shp.customFeatures.clear();
  auto len = nb::len(s);
  shp.customFeatures.reserve(len);
  for (auto i = 0u; i < len; ++i) {
    const auto elem = s[i];
    unsigned int featType = nb::cast<unsigned int>(elem[0]);
    RDGeom::Point3D pos = nb::cast<RDGeom::Point3D>(elem[1]);
    double radius = nb::cast<double>(elem[2]);
    shp.customFeatures.emplace_back(featType, pos, radius);
  }
}

nb::tuple get_customFeatures(const GaussianShape::ShapeInputOptions &shp) {
  nb::list py_list;
  for (const auto &val : shp.customFeatures) {
    nb::list elem;
    elem.append(static_cast<int>(std::get<0>(val)));
    elem.append(std::get<1>(val));
    elem.append(std::get<2>(val));
    py_list.append(nb::tuple(elem));
  }
  return nb::tuple(py_list);
}

}  // namespace helpers

NB_MODULE(rdGaussianShape, m) {
  m.doc() =
      R"DOC(Module containing implementation of Gaussian-based shape overlay and scoring.
NOTE: This functionality is experimental and the API and/or results may change in future releases.)DOC";

  nb::enum_<RDKit::GaussianShape::StartMode>(m, "StartMode")
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

  nb::enum_<RDKit::GaussianShape::OptimMode>(m, "OptimMode")
      .value("SHAPE_ONLY", RDKit::GaussianShape::OptimMode::SHAPE_ONLY)
      .value("SHAPE_PLUS_COLOR_SCORE",
             RDKit::GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE)
      .value("SHAPE_PLUS_COLOR",
             RDKit::GaussianShape::OptimMode::SHAPE_PLUS_COLOR)
      .export_values();

  nb::class_<GaussianShape::ShapeInputOptions>(
      m, "ShapeInputOptions",
      "ShapeInputOptions - options for setting up ShapeInput objects.")
      .def(nb::init<>())
      .def_rw("useColors", &GaussianShape::ShapeInputOptions::useColors,
              "Whether to use color features in overlay.  Default=True.")
      .def_rw(
          "allCarbonRadii", &GaussianShape::ShapeInputOptions::allCarbonRadii,
          R"DOC(Whether to use the same radius, appropriate for Carbon, for all atoms.  There is a
slight accuracy penalty but significant speed gain if used.  Default=True.)DOC")
      .def_prop_rw(
          "atomSubset",
          [](const GaussianShape::ShapeInputOptions &opts) {
            nb::list py_list;
            for (const auto &val : opts.atomSubset) {
              py_list.append(val);
            }
            return nb::tuple(py_list);
          },
          [](GaussianShape::ShapeInputOptions &opts, nb::object as) {
            pythonObjectToVect<unsigned int>(as, opts.atomSubset);
          },
          "If not empty, use just these atoms in the molecule to form the ShapeInput object.")
      .def_prop_rw(
          "customFeatures", &helpers::get_customFeatures,
          &helpers::set_customFeatures,
          R"DOC(Custom features for the shape.  Requires a list of tuples of
int (the feature type), Point3D (the coordinates) and float (the radius).)DOC")
      .def_prop_rw(
          "atomRadii",
          [](const GaussianShape::ShapeInputOptions &opts) {
            nb::list py_list;
            for (const auto &val : opts.atomRadii) {
              py_list.append(
                  nb::make_tuple(static_cast<int>(val.first), val.second));
            }
            return nb::tuple(py_list);
          },
          [](GaussianShape::ShapeInputOptions &opts, nb::object ar) {
            int len = nb::len(ar);
            opts.atomRadii.resize(len);
            for (int i = 0; i < len; i++) {
              unsigned int atomIdx = nb::cast<unsigned int>(ar[i][0]);
              double radius = nb::cast<double>(ar[i][1]);
              opts.atomRadii[i] = std::make_pair(atomIdx, radius);
            }
          },
          R"DOC(Non-standard radii to use for the atoms specified by their indices
in the molecule.  Not all atoms need have a radius specified.
A list of tuples of [int, float].)DOC");

  nb::class_<GaussianShape::ShapeOverlayOptions>(
      m, "ShapeOverlayOptions",
      "ShapeOverlayOptions - options for controlling the shape overlay process.")
      .def(nb::init<>())
      .def_rw(
          "startMode", &RDKit::GaussianShape::ShapeOverlayOptions::startMode,
          R"DOC(Start modes for optimisation.  Default is A_LA_PUBCHEM - as used by the
PubChem code - either ROTATE_180_WIGGLE or ROTATE_45 depending on the shape
of the two molecules.  ROTATE_180_WIGGLE means 180 rotations about
the x, y and z axes, then a small
rotation about each axis from that point, using the best scoring one of
those. ROTATE_180 uses 180 degree rotations for 4 start points,
ROTATE_45 uses 45 degree rotations for 9 start points and ROTATE_0
leaves the relative orientations of the 2 molecules as passed in before
optimisation.  There are also ROTATE_0_FRAGMENT, ROTATE_45_FRAGMENT
and ROTATE_180_FRAGMENT that as well as the above move the fit
molecule to the ends of each of the principal axes and then does
the appropriate rotations.  This is useful when the fit molecule is
a lot smaller than the reference molecule, but requires a large number
of optimisations so is relatively slow.)DOC")
      .def_rw("optimMode", &GaussianShape::ShapeOverlayOptions::optimMode,
              R"DOC(Optimisation mode, controlling what parameters are used
to drive the overlay.  Default=SHAPE_PLUS_COLOR_SCORE which
optimises using just the overlap of shape, but uses the
color to decide which is the best overlay.  Other options
are SHAPE_ONLY and SHAPE_AND_COLOR with the latter using
the overlap of color features as well. )DOC")
      .def_rw(
          "simAlpha", &GaussianShape::ShapeOverlayOptions::simAlpha,
          R"DOC(When doing a Tversky similarity, the alpha value.  If alpha and
beta are both the default 1.0, it's a Tanimoto similarity.  A
high alpha and low beta emphasize the fit volume in the
similarity and vice versa. Tversky is O / (A * (R - O) + B * (F
- O) + O) where O is the overlap volume, R is the reference's
volume and F is the fit's volume.  This is different from that
used by OpenEye (O / (A * R + B * F)).)DOC")
      .def_rw("simBeta", &GaussianShape::ShapeOverlayOptions::simBeta,
              "When doing a Tversky similarity, the beta value.")
      .def_rw(
          "optParam", &GaussianShape::ShapeOverlayOptions::optParam,
          R"DOC(If using colors, the relative weights of the shape and color scores,
as a fraction of 1.  Default=0.5.)DOC")
      .def_rw(
          "nSteps", &GaussianShape::ShapeOverlayOptions::nSteps,
          "Maximum number of steps for the shape overlay process. Default=100.")
      .def_rw(
          "normalize", &GaussianShape::ShapeOverlayOptions::normalize,
          R"DOC(Whether to normalize the shapes before overlay by putting them into their
canonical orientation (centred on the origin, aligned along its
principal axes.  Default=True.)DOC")
      .def_rw(
          "useDistCutoff", &GaussianShape::ShapeOverlayOptions::useDistCutoff,
          R"DOC(Whether to use distance cutoff when calculating the shape volumes.  If used,
there will be a small penalty in accuracy but a significant increase in speed.
Default=True.)DOC")
      .def_rw(
          "distCutoff", &GaussianShape::ShapeOverlayOptions::distCutoff,
          R"DOC(If using a distance cutoff, this is the value used.  Default=4.5 of whatever
units the coordinates are in.)DOC")
      .def_rw(
          "shapeConvergenceCriterion",
          &GaussianShape::ShapeOverlayOptions::shapeConvergenceCriterion,
          R"DOC(Optimisation stops when the shape Tversky score changes by less
than this amount after an optimisation step.  A larger number is
faster but gives less precise overlays.  Default=0.001.)DOC");

  nb::class_<GaussianShape::ShapeInput>(m, "ShapeInput", "ShapeInput object")
      .def(
          nb::init<const ROMol &, int, const GaussianShape::ShapeInputOptions &,
                   const GaussianShape::ShapeOverlayOptions &>(),
          "mol"_a, "confId"_a = -1,
          "shapeOpt"_a = GaussianShape::ShapeInputOptions(),
          "overlayOpts"_a = GaussianShape::ShapeOverlayOptions())
      .def_prop_ro("NumAtoms", &GaussianShape::ShapeInput::getNumAtoms,
                   "Get the number of atoms defining the shape.")
      .def_prop_ro("NumFeatures", &GaussianShape::ShapeInput::getNumFeatures,
                   "Get the number of features in the shape.")
      .def_prop_ro("ShapeVolume", &GaussianShape::ShapeInput::getShapeVolume,
                   "Get the shape's volume due to the atoms.")
      .def_prop_ro("ColorVolume", &GaussianShape::ShapeInput::getColorVolume,
                   "Get the volume of the shape's color features.")
      .def("__setattr__", &safeSetattr);

  m.def(
      "AlignMol",
      [](const ROMol &ref, ROMol &fit,
         std::optional<GaussianShape::ShapeInputOptions> py_refOpts,
         std::optional<GaussianShape::ShapeInputOptions> py_fitOpts,
         std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts,
         int refConfId, int fitConfId) {
        auto results = GaussianShape::AlignMolecule(
            ref, fit, py_refOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_fitOpts.value_or(GaussianShape::ShapeInputOptions()), nullptr,
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            refConfId, fitConfId);
        return nb::make_tuple(results[0], results[1], results[2]);
      },
      "ref"_a, "fit"_a, "refOpts"_a = nb::none(), "fitOpts"_a = nb::none(),
      "overlayOpts"_a = nb::none(), "refConfId"_a = -1, "fitConfId"_a = -1,
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

  m.def(
      "AlignMol",
      [](const GaussianShape::ShapeInput &refShape, ROMol &fit,
         std::optional<GaussianShape::ShapeInputOptions> py_fitOpts,
         std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts,
         int fitConfId) {
        auto results = GaussianShape::AlignMolecule(
            refShape, fit,
            py_fitOpts.value_or(GaussianShape::ShapeInputOptions()), nullptr,
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            fitConfId);
        return nb::make_tuple(results[0], results[1], results[2]);
      },
      "refShape"_a, "fit"_a, "fitOpts"_a = nb::none(),
      "overlayOpts"_a = nb::none(), "fitConfId"_a = -1,
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

  m.def(
      "AlignShapes",
      [](const GaussianShape::ShapeInput &refShape,
         GaussianShape::ShapeInput &fitShape,
         std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts) {
        RDGeom::Transform3D xform;
        auto results = GaussianShape::AlignShape(
            refShape, fitShape, &xform,
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()));
        nb::list pyMatrix;
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            pyMatrix.append(xform.getValUnchecked(i, j));
          }
        }
        return nb::make_tuple(results[0], results[1], results[2], pyMatrix);
      },
      "refShape"_a, "fitShape"_a, "overlayOpts"_a = nb::none(),
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

  m.def(
      "ScoreMol",
      [](const ROMol &ref, const ROMol &fit,
         std::optional<GaussianShape::ShapeInputOptions> py_refOpts,
         std::optional<GaussianShape::ShapeInputOptions> py_fitOpts,
         std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts,
         int refConfId, int fitConfId) {
        auto results = GaussianShape::ScoreMolecule(
            ref, fit, py_refOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_fitOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            refConfId, fitConfId);
        return nb::make_tuple(results[0], results[1], results[2]);
      },
      "ref"_a, "fit"_a, "refOpts"_a = nb::none(), "fitOpts"_a = nb::none(),
      "overlayOpts"_a = nb::none(), "refConfId"_a = -1, "fitConfId"_a = -1,
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

  m.def(
      "ScoreMol",
      [](const GaussianShape::ShapeInput &refShape, const ROMol &fit,
         std::optional<GaussianShape::ShapeInputOptions> py_fitOpts,
         std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts,
         int fitConfId) {
        auto results = GaussianShape::ScoreMolecule(
            refShape, fit,
            py_fitOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            fitConfId);
        return nb::make_tuple(results[0], results[1], results[2]);
      },
      "refShape"_a, "fit"_a, "fitOpts"_a = nb::none(),
      "overlayOpts"_a = nb::none(), "fitConfId"_a = -1,
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

  m.def(
      "ScoreShape",
      [](const GaussianShape::ShapeInput &refShape,
         const GaussianShape::ShapeInput &fitShape,
         std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts) {
        auto results = GaussianShape::ScoreShape(
            refShape, fitShape,
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()));
        return nb::make_tuple(results[0], results[1], results[2]);
      },
      "refShape"_a, "fitShape"_a, "overlayOpts"_a = nb::none(),
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

}  // namespace RDKit
