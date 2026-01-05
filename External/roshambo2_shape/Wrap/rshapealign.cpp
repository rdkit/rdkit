//
//  Copyright (C) 2021-2025 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <boost/python.hpp>

#include <GraphMol/ROMol.h>
#include "numpy/arrayobject.h"
#include <RDBoost/Wrap.h>
#include "../Roshambo2Shape.hpp"
#include "../ShapeInput.h"
#include "../ShapeOverlayOptions.h"

namespace python = boost::python;
using DTYPE = float;

namespace helpers {

python::tuple alignMol1(const RDKit::ROMol &ref, RDKit::ROMol &probe,
                        const python::object &py_opts, int refConfId,
                        int fitConfId) {
  RDKit::ShapeAlign::ShapeOverlayOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RDKit::ShapeAlign::ShapeOverlayOptions>(py_opts);
  }

  auto [st, ct] = RDKit::ShapeAlign::AlignMolecule(ref, probe, nullptr, opts,
                                                   refConfId, fitConfId);
  return python::make_tuple(st, ct);
}

python::tuple alignMol2(const RDKit::ShapeAlign::ShapeInput &refShape,
                        RDKit::ROMol &probe, const python::object &py_opts,
                        int fitConfId) {
  RDKit::ShapeAlign::ShapeOverlayOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RDKit::ShapeAlign::ShapeOverlayOptions>(py_opts);
  }

  auto [st, ct] = RDKit::ShapeAlign::AlignMolecule(refShape, probe, nullptr,
                                                   opts, fitConfId);
  return python::make_tuple(st, ct);
}

python::tuple alignShape(const RDKit::ShapeAlign::ShapeInput &refShape,
                         RDKit::ShapeAlign::ShapeInput &probeShape,
                         const python::object &py_opts) {
  RDKit::ShapeAlign::ShapeOverlayOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RDKit::ShapeAlign::ShapeOverlayOptions>(py_opts);
  }

  RDGeom::Transform3D xform;
  auto [st, ct] =
      RDKit::ShapeAlign::AlignShape(refShape, probeShape, &xform, opts);
  python::list pyMatrix;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      pyMatrix.append(xform.getValUnchecked(i, j));
    }
  }
  return python::make_tuple(st, ct, pyMatrix);
}

python::tuple scoreShape(const RDKit::ShapeAlign::ShapeInput &refShape,
                         const RDKit::ShapeAlign::ShapeInput &probeShape,
                         const python::object &py_opts) {
  RDKit::ShapeAlign::ShapeOverlayOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RDKit::ShapeAlign::ShapeOverlayOptions>(py_opts);
  }

  auto [st, ct] = RDKit::ShapeAlign::ScoreShape(refShape, probeShape, opts);
  return python::make_tuple(st, ct);
}

python::tuple scoreMol1(const RDKit::ROMol &ref, const RDKit::ROMol &probe,
                        const python::object &py_opts, int refConfId,
                        int probeConfId) {
  RDKit::ShapeAlign::ShapeOverlayOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RDKit::ShapeAlign::ShapeOverlayOptions>(py_opts);
  }

  auto [st, ct] = RDKit::ShapeAlign::ScoreMolecule(ref, probe, opts, refConfId,
                                                   probeConfId);
  return python::make_tuple(st, ct);
}

python::tuple scoreMol2(const RDKit::ShapeAlign::ShapeInput &refShape,
                        const RDKit::ROMol &probe,
                        const python::object &py_opts, int probeConfId) {
  RDKit::ShapeAlign::ShapeOverlayOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RDKit::ShapeAlign::ShapeOverlayOptions>(py_opts);
  }

  auto [st, ct] =
      RDKit::ShapeAlign::ScoreMolecule(refShape, probe, opts, probeConfId);
  return python::make_tuple(st, ct);
}

python::list getShapeType_helper(const RDKit::ShapeAlign::ShapeInput &shape) {
  python::list types;
  for (auto &t : shape.getTypes()) {
    types.append(t);
  }
  return types;
}
}  // namespace helpers

void wrap_roshambo2shape() {
  python::class_<RDKit::ShapeAlign::ShapeOverlayOptions, boost::noncopyable>(
      "ShapeOverlayOptions", "Shape Overlay Options")
      .def_readwrite(
          "useColors", &RDKit::ShapeAlign::ShapeOverlayOptions::d_useColors,
          "Whether to use colors (pharmacophore features) in the score.  Default=True.")
      .def_readwrite("normalize",
                     &RDKit::ShapeAlign::ShapeOverlayOptions::d_normalize,
                     "Whether to normalise the shape by putting it into"
                     "its inertial frame.  Default=True.")
      .def_readwrite(
          "optParam", &RDKit::ShapeAlign::ShapeOverlayOptions::d_optParam,
          "If using colors, the relative weights fo shape and color scores. Default=0.5.")
      .def_readwrite("nSteps",
                     &RDKit::ShapeAlign::ShapeOverlayOptions::d_nSteps,
                     "Number of steps of the optimiser to take.  Default=100.")
      .def("__setattr__", &safeSetattr);

  std::string docString =
      "Shape object for use with Roshambo2-based shape alignment.";
  python::class_<RDKit::ShapeAlign::ShapeInput, boost::noncopyable>(
      "ShapeInput", docString.c_str(),
      python::init<RDKit::ROMol &, int,
                   RDKit::ShapeAlign::ShapeOverlayOptions &>(
          (python::arg("self"), python::arg("mol"), python::arg("confId") = -1,
           python::arg("opts") = python::object())))
      .def("GetNumAtoms", &RDKit::ShapeAlign::ShapeInput::getNumAtoms,
           "Get number of atoms in the shape.")
      .def("GetNumFeatures", &RDKit::ShapeAlign::ShapeInput::getNumFeatures,
           "Get number of features in the shape.")
      .def(
          "GetTypes", &helpers::getShapeType_helper,
          "Get the types of the atoms and features.  Currently, atom types are all 0.")
      .def("GetVolume", &RDKit::ShapeAlign::ShapeInput::getSelfOverlapVol,
           "Get the shape volume.")
      .def("GetColorVolume",
           &RDKit::ShapeAlign::ShapeInput::getSelfOverlapColor,
           "Get the color volume.");

  docString =
      R"DOC(Aligns a probe molecule to a reference molecule. The probe is modified.

Parameters
----------
ref : RDKit.ROMol
    Reference molecule
probe : RDKit.ROMol
    Probe molecule
opts : ShapeOverlayOptions, optional
    Options for the overlay
refConfId : int, optional
    Reference conformer ID (default is -1)
probeConfId : int, optional
    Probe conformer ID (default is -1)


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if color optimisation not selected is 1.0.)DOC";

  python::def("AlignMol", &helpers::alignMol1,
              (python::arg("ref"), python::arg("probe"),
               python::arg("opts") = python::object(),
               python::arg("refConfId") = -1, python::arg("probeConfId") = -1),
              docString.c_str());

  docString =
      R"DOC(Aligns a probe molecule to a reference shape. The probe is modified.

Parameters
----------
ref : RDKit.ROMol
    Reference molecule
probe : RDKit.ROMol
    Probe molecule
opts : ShapeOverlayOptions, optional
    Options for the overlay
probeConfId : int, optional
    Probe conformer ID (default is -1)


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if color optimisation not selected is 1.0.)DOC";

  python::def(
      "AlignMol", &helpers::alignMol2,
      (python::arg("refShape"), python::arg("probe"),
       python::arg("opts") = python::object(), python::arg("probeConfId") = -1),
      docString.c_str());

  docString =
      R"DOC(Aligns a probe shape to a reference shape. The probe is modified.

Parameters
----------
refShape : ShapeInput
    Reference shape
probeShape : ShapeInput
    Probe shape
opts : ShapeOverlayOptions, optional
    Options for the overlay


Returns
-------
 3-tuple of double, double, list of doubles
    The results are (shape_score, color_score, matrix)
    The matrix is a 16-float list giving the transformation matrix that
    overlays the probe onto the reference.  To turn it into a numpy array
    that can be used to transform conformations, do something like:
    ttrans = [tpl[2][0:4], tpl[2][4:8], tpl[2][8:12], tpl[2][12:16]]
    nptrans = np.array(ttrans)
)DOC";

  python::def("AlignShape", &helpers::alignShape,
              (python::arg("refShape"), python::arg("probeShape"),
               python::arg("opts") = python::object()),
              docString.c_str());

  docString =
      R"DOC(Score the overlap of a shape to a reference shape without moving
either.  Note that if you take the output from one of the Align...
functions and feed it into a Score... function you won't get
exactly the same answer.  This is because the formula for the
tanimoto uses the fit volume and that is calculated once at the
start before the rotations and translations that form the
optimisation.  Floating point cruft moves the atoms by small
amounts relative to each other which means that the final
calculated volume differs slightly from the one calculated at
the start.  The fixed scoring obviously doesn't have this effect.

Parameters
----------
refShape : ShapeInput
    Reference shape
probeShape : ShapeInput
    Probe shape
opts : ShapeOverlayOptions, optional
    Options for the scoring


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if color optimisation not selected is 1.0.)DOC";
  python::def("ScoreShape", &helpers::scoreShape,
              (python::arg("refShape"), python::arg("probeShape"),
               python::arg("opts") = python::object()),
              docString.c_str());

  docString =
      R"DOC(Score the overlap of a molecule to a reference shape without moving
either.
Parameters
----------
refShape : ShapeInput
    Reference shape
probe : RDKit.ROMol
    Probe molecule
opts : ShapeOverlayOptions, optional
    Options for the scoring
probeConfId : int, optional
    Probe conformer ID (default is -1)


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if color optimisation not selected is 1.0.)DOC";

  python::def(
      "ScoreMol", &helpers::scoreMol2,
      (python::arg("refShape"), python::arg("probe"),
       python::arg("opts") = python::object(), python::arg("probeConfId") = -1),
      docString.c_str());

  docString =
      R"DOC(Score the overlap of a molecule to a reference molecule without moving
either.
Parameters
----------
ref : RDKit.ROMol
    Reference molecule
probe : RDKit.ROMol
    Probe molecule
opts : ShapeOverlayOptions, optional
    Options for the scoring
refConfId : int, optional
    Reference conformer ID (default is -1)
probeConfId : int, optional
    Probe conformer ID (default is -1)


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if color optimisation not selected is 1.0.)DOC";

  python::def("ScoreMol", &helpers::scoreMol1,
              (python::arg("ref"), python::arg("probe"),
               python::arg("opts") = python::object(),
               python::arg("refConfId") = -1, python::arg("probeConfId") = -1),
              docString.c_str());
}

BOOST_PYTHON_MODULE(rdShapeAlign2) { wrap_roshambo2shape(); }
