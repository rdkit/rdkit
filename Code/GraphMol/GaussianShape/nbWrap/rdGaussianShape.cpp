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
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
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
  auto numVecs = nb::len(s);
  shp.customFeatures.reserve(numVecs);
  for (auto i = 0u; i < numVecs; ++i) {
    nb::object outVec = s[i];
    auto numFeats = nb::len(outVec);
    std::vector<GaussianShape::CustomFeature> feats;
    feats.reserve(numFeats);
    for (auto j = 0u; j < numFeats; ++j) {
      nb::object feat = outVec[j];
      unsigned int featType = nb::cast<unsigned int>(feat[0]);
      RDGeom::Point3D pos = nb::cast<RDGeom::Point3D>(feat[1]);
      double radius = nb::cast<double>(feat[2]);
      std::vector<unsigned int> atoms;
      if (nb::len(feat) == 4) {
        nb::object featAtoms = feat[3];
        for (unsigned int k = 0; k < nb::len(featAtoms); ++k) {
          atoms.push_back(nb::cast<unsigned int>(featAtoms[k]));
        }
      }
      feats.emplace_back(featType, pos, radius, atoms);
    }
    shp.customFeatures.emplace_back(std::move(feats));
  }
}

nb::tuple get_customFeatures(const GaussianShape::ShapeInputOptions &shp) {
  nb::list allFeatLists;
  for (const auto &feats : shp.customFeatures) {
    nb::list featList;
    for (const auto &feat : feats) {
      nb::list elem;
      elem.append(static_cast<int>(feat.type));
      elem.append(feat.pos);
      elem.append(feat.rad);
      elem.append(feat.atoms);
      featList.append(elem);
    }
    allFeatLists.append(featList);
  }
  return nb::tuple(allFeatLists);
}

double getShapeVolume_helper(const GaussianShape::ShapeInput &shape) {
  return shape.getShapeVolume();
}

double getColorVolume_helper(const GaussianShape::ShapeInput &shape) {
  return shape.getColorVolume();
}

nb::tuple bestSimilarity_helper(
    GaussianShape::ShapeInput &refShape,
    const GaussianShape::ShapeInput &fitShape, double threshold,
    std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts) {
  unsigned int bestFitShape, bestThisShape;
  RDGeom::Transform3D bestXform;
  auto bestSim = refShape.bestSimilarity(
      fitShape, bestThisShape, bestFitShape, bestXform, threshold,
      py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()));
  nb::list pyMatrix;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      pyMatrix.append(bestXform.getValUnchecked(i, j));
    }
  }
  return nb::make_tuple(nb::make_tuple(bestSim[0], bestSim[1], bestSim[2]),
                        bestThisShape, bestFitShape, pyMatrix);
}

double maxPossibleSimilarity_helper(
    GaussianShape::ShapeInput &refShape, GaussianShape::ShapeInput &fitShape,
    std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts) {
  return refShape.maxPossibleSimilarity(
      fitShape, py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()));
}

ROMol *shapeToMol_helper(GaussianShape::ShapeInput &shape, bool includeColors,
                         bool withBonds) {
  auto mol = shape.shapeToMol(includeColors, withBonds);
  return static_cast<ROMol *>(mol.release());
}

nb::tuple scoreMolAllConfs_helper(
    const ROMol &ref, const ROMol &fit,
    std::optional<GaussianShape::ShapeInputOptions> py_refOpts,
    std::optional<GaussianShape::ShapeInputOptions> py_fitOpts,
    std::optional<GaussianShape::ShapeOverlayOptions> py_overlayOpts) {
  std::vector<std::vector<double>> combScores;
  int bestRefConf, bestFitConf;
  RDGeom::Transform3D bestXform;
  GaussianShape::ScoreMoleculeAllConformers(
      ref, fit, bestRefConf, bestFitConf, combScores,
      py_refOpts.value_or(GaussianShape::ShapeInputOptions()),
      py_fitOpts.value_or(GaussianShape::ShapeInputOptions()),
      py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
      &bestXform);

  nb::list pyScores;
  for (const auto &scores : combScores) {
    nb::list s;
    for (const auto &score : scores) {
      s.append(score);
    }
    pyScores.append(s);
  }

  nb::list pyMatrix;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      pyMatrix.append(bestXform.getValUnchecked(i, j));
    }
  }
  return nb::make_tuple(pyScores, bestRefConf, bestFitConf, pyMatrix);
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
          R"DOC(Custom features for the shape.  Requires a list of lists of tuples of
int (the feature type), Point3D (the coordinates), float (the radius)
and optionally a list of indices of the atoms that the feature was derived from.)DOC")
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
    A list of tuples of [int, float].)DOC")
      .def_rw(
          "shapePruneThreshold",
          &GaussianShape::ShapeInputOptions::shapePruneThreshold,
          "If there is more than 1 conformer for the input molecule, prune the"
          " shapes so that none of them are more similar to each other than the"
          " threshold.  Default -1.0 means no pruning.")
      .def_rw(
          "sortShapes", &GaussianShape::ShapeInputOptions::sortShapes,
          "If True (the default), the shapes are sorted into descending order"
          " of total volume.")
      .def_rw(
          "includeDummies", &GaussianShape::ShapeInputOptions::includeDummies,
          "Whether to include dummy atoms in the shape or not.  Default=True.")
      .def("__setattr__", &safeSetattr);

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
faster but gives less precise overlays.  Default=0.001.)DOC")
      .def("__setattr__", &safeSetattr);

  nb::class_<GaussianShape::ShapeInput>(m, "ShapeInput", "ShapeInput object")
      .def(
          nb::init<const ROMol &, int, const GaussianShape::ShapeInputOptions &,
                   const GaussianShape::ShapeOverlayOptions &>(),
          "mol"_a, "confId"_a = -1,
          "shapeOpt"_a = GaussianShape::ShapeInputOptions(),
          "overlayOpts"_a = GaussianShape::ShapeOverlayOptions())
      .def_prop_ro(
          "GetSmiles", &GaussianShape::ShapeInput::getSmiles,
          "Get the SMILES string for the molecule that the shape relates to.")
      .def("setActiveShape", &GaussianShape::ShapeInput::setActiveShape,
           "Set the active shape, the one that will be used for overlays etc.")
      .def("getActiveShape", &GaussianShape::ShapeInput::getActiveShape,
           "Return the number of the active shape.")
      .def_prop_ro("NumAtoms", &GaussianShape::ShapeInput::getNumAtoms,
                   "Get the number of atoms defining the shape.")
      .def_prop_ro("NumFeatures", &GaussianShape::ShapeInput::getNumFeatures,
                   "Get the number of features in the shape.")
      .def_prop_ro(
          "NumShapes", &GaussianShape::ShapeInput::getNumShapes,
          "Get the number of shapes.  There will be a shape for each conformation "
          "of the input molecule, unless shape pruning was carried out in which case"
          " there may be fewer.")
      .def_prop_ro("ShapeVolume", &helpers::getShapeVolume_helper,
                   "Get the volume due to the atoms for the active shape.")
      .def_prop_ro(
          "ColorVolume", &helpers::getColorVolume_helper,
          "Get the volume of the shape's color features for the active shape.")
      .def(
          "NormalizeCoords", &GaussianShape::ShapeInput::normalizeCoords,
          "Align the principal axes to the cartesian axes and centre on the origin."
          " Doesn't require that the shape was created from a molecule.  Creates"
          " the necessary transformation if not already done.")
      .def(
          "ShapeToMol", &helpers::shapeToMol_helper, "includeColors"_a = false,
          "withBonds"_a = true,
          "Return a molecule with coordinates of the current active shape."
          "  If includeColors is True, (default is False) the color features"
          " will be added as xenon atoms.  If withBonds is True (the default)"
          " a molecule with bonds will be created, if not then just atoms at the"
          " appropriate positions will be produced.",
          nb::rv_policy::take_ownership)
      .def(
          "BestSimilarity", &helpers::bestSimilarity_helper, "fitShape"_a,
          "threshold"_a = -1.0, "overlayOpts"_a = nb::none(),
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
           "fitShape"_a, "overlayOpts"_a = nb::none(),
           "Get the maximum possible similarity score between all shapes in"
           " this shape and all shapes in the fitShape.  The maximum similarity"
           " is when one shape is entirely inside the other.  This returns"
           " the similarity in that case, which is the upper bound on what"
           " is achievable between these 2 shapes.")
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
        std::pair<double, double> ovVols;
        auto results = GaussianShape::ScoreMolecule(
            ref, fit, py_refOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_fitOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            refConfId, fitConfId, &ovVols);
        return nb::make_tuple(results[0], results[1], results[2], ovVols.first,
                              ovVols.second);
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
        std::pair<double, double> ovVols;
        auto results = GaussianShape::ScoreMolecule(
            refShape, fit,
            py_fitOpts.value_or(GaussianShape::ShapeInputOptions()),
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            fitConfId, &ovVols);
        return nb::make_tuple(results[0], results[1], results[2], ovVols.first,
                              ovVols.second);
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
        std::pair<double, double> ovVols;
        auto results = GaussianShape::ScoreShape(
            refShape, fitShape,
            py_overlayOpts.value_or(GaussianShape::ShapeOverlayOptions()),
            &ovVols);
        return nb::make_tuple(results[0], results[1], results[2], ovVols.first,
                              ovVols.second);
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

  m.def("ScoreMoleculeAllConformers", &helpers::scoreMolAllConfs_helper,
        "ref"_a, "fit"_a, "refOpts"_a = nb::none(), "fitOpts"_a = nb::none(),
        "overlayOpts"_a = nb::none(),
        R"DOC(Calculate the scores for the alignment of all conformers
 of the fit molecule onto the reference.  The molecules themselves are not
 altered.

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

Returns
-------
A complex tuple containing:
    A tuple of tuples containing the scores from aligning the fit conformations
    onto the reference conformations.  scores[0][1] is the score of aligning
    fit conformation 1 onto ref conformation 0.
    The ID of the ref conformer from the best-scoring alignment
    The ID of the fit conformer from the best-scoring alignment
    The transformation that gives the best-scoring alignment for those
    conformers as a 16-float tuple.
)DOC");
}

}  // namespace RDKit
