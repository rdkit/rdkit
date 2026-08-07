/*******************************************************************************

Copyright 2024-2026 by Greg Landrum and the pubchem_shape contributors and
other RDKit contributors

This file is part of pubchem_shape

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********************************************************************/

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include <vector>

#include "../PubChemShape.hpp"
#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/Wrap_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(rdShapeAlign, m) {
  m.doc() = "Module containing PubChem shape alignment functionality.";

  nb::class_<ShapeInputOptions>(m, "ShapeInputOptions", "Shape Input Options")
      .def(nb::init<>())
      .def_rw(
          "useColors", &ShapeInputOptions::useColors,
          "Whether to use colors (pharmacophore features) in the score.  Default=True.")
      .def_rw(
          "includeDummies", &ShapeInputOptions::includeDummies,
          "Whether to use dummy atoms in the alignment. Default=False.")
      .def_rw(
          "dummyRadius", &ShapeInputOptions::dummyRadius,
          R"DOC(If using dummy atoms in the alignment, what radius to use for them.
  Default=2.16 (the radius of Xe).)DOC")
      .def_prop_rw(
          "atomSubset",
          [](const ShapeInputOptions &opts) {
            nb::list py_list;
            for (const auto &val : opts.atomSubset) {
              py_list.append(val);
            }
            return nb::tuple(py_list);
          },
          [](ShapeInputOptions &opts, nb::object as) {
            pythonObjectToVect<unsigned int>(as, opts.atomSubset);
          },
          "If not empty, use just these atoms in the molecule to form the ShapeInput object.")
      .def_prop_rw(
          "notColorAtoms",
          [](const ShapeInputOptions &opts) {
            nb::list py_list;
            for (const auto &val : opts.notColorAtoms) {
              py_list.append(val);
            }
            return nb::tuple(py_list);
          },
          [](ShapeInputOptions &opts, nb::object nca) {
            pythonObjectToVect<unsigned int>(nca, opts.notColorAtoms);
          },
          "Any atoms mentioned here by index should not be used in a color feature.")
      .def_prop_rw(
          "atomRadii",
          [](const ShapeInputOptions &opts) {
            nb::list py_list;
            for (const auto &val : opts.atomRadii) {
              py_list.append(nb::make_tuple(static_cast<int>(val.first), val.second));
            }
            return nb::tuple(py_list);
          },
          [](ShapeInputOptions &opts, nb::object ar) {
            int len = nb::len(ar);
            opts.atomRadii.resize(len);
            for (int i = 0; i < len; i++) {
              unsigned int atomIdx = nb::cast<unsigned int>(ar[i][0]);
              double radius = nb::cast<double>(ar[i][1]);
              opts.atomRadii[i] = std::make_pair(atomIdx, radius);
            }
          },
          "Non-standard radii to use for the atoms specified by their indices"
          " in the molecule.  A list of tuples of [int, float].")
      .def_prop_rw(
          "customFeatures",
          [](const ShapeInputOptions &shp) {
            nb::list py_list;
            for (const auto &val : shp.customFeatures) {
              nb::list elem;
              elem.append(static_cast<int>(std::get<0>(val)));
              elem.append(std::get<1>(val));
              elem.append(std::get<2>(val));
              py_list.append(nb::tuple(elem));
            }
            return nb::tuple(py_list);
          },
          [](ShapeInputOptions &shp, nb::object s) {
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
          },
          "Custom features for the shape.")
      .def_rw("normalize", &ShapeInputOptions::normalize,
              R"DOC(Whether to normalise the shape by putting into
its inertial frame.  Default=True.)DOC")
      .def("__setattr__", &safeSetattr);

  nb::class_<ShapeInput>(m, "ShapeInput")
      .def_rw("coord", &ShapeInput::coord)
      .def_rw("alpha_vector", &ShapeInput::alpha_vector)
      .def_rw("atom_type_vector", &ShapeInput::atom_type_vector)
      .def_rw("volumeAtomIndexVector", &ShapeInput::volumeAtomIndexVector)
      .def_prop_rw(
          "shift",
          [](const ShapeInput &shp) {
            nb::list py_list;
            for (const auto &val : shp.shift) {
              py_list.append(val);
            }
            return py_list;
          },
          [](ShapeInput &shp, nb::object s) {
            pythonObjectToVect<double>(s, shp.shift);
          },
          "Translation of centre of shape coordinates to origin.")
      .def_prop_ro(
          "inertialRot",
          [](const ShapeInput &shp) {
            nb::list py_list;
            for (const auto &val : shp.inertialRot) {
              py_list.append(val);
            }
            return py_list;
          },
          "Rotation applied to put the shape into its principal axes frame of reference.")
      .def_rw("sov", &ShapeInput::sov)
      .def_rw("sof", &ShapeInput::sof);

  m.def(
      "AlignMol",
      [](const RDKit::ROMol &ref, RDKit::ROMol &probe, int refConfId,
         int probeConfId, bool useColors, double opt_param,
         unsigned int max_preiters, unsigned int max_postiters) {
        std::vector<float> matrix(12, 0.0);
        auto [nbr_st, nbr_ct] =
            AlignMolecule(ref, probe, matrix, refConfId, probeConfId, useColors,
                          opt_param, max_preiters, max_postiters);
        return nb::make_tuple(nbr_st, nbr_ct);
      },
      "ref"_a, "probe"_a, "refConfId"_a = -1, "probeConfId"_a = -1,
      "useColors"_a = true, "opt_param"_a = 1.0, "max_preiters"_a = 10,
      "max_postiters"_a = 30,
      R"DOC(Aligns a probe molecule to a reference molecule. The probe is modified.

Parameters
----------
ref : RDKit.ROMol
    Reference molecule
probe : RDKit.ROMol
    Probe molecule
refConfId : int, optional
    Reference conformer ID (default is -1)
probeConfId : int, optional
    Probe conformer ID (default is -1)
useColors : bool, optional
    Whether or not to use colors in the scoring (default is True)
opt_param : float, optional
    Balance of shape and color for optimization.
    0 is only color, 0.5 is equal weight, and 1.0 is only shape.
    Default is 1.0.
max_preiters : int, optional
    In the two phase optimization, the maximum iterations done on all poses.
max_postiters : int, optional
    In the two phase optimization, the maximum iterations during the second phase on
    only the best poses from the first phase


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if opt_param is 1.0.)DOC");

  m.def(
      "AlignMol",
      [](const RDKit::ROMol &ref, RDKit::ROMol &probe,
         const ShapeInputOptions &refShapeOpts,
         const ShapeInputOptions &probeShapeOpts, int refConfId, int probeConfId,
         double opt_param, unsigned int max_preiters,
         unsigned int max_postiters) {
        std::vector<float> matrix(12, 0.0);
        auto [nbr_st, nbr_ct] =
            AlignMolecule(ref, probe, matrix, refShapeOpts, probeShapeOpts,
                          refConfId, probeConfId, opt_param, max_preiters,
                          max_postiters);
        return nb::make_tuple(nbr_st, nbr_ct);
      },
      "ref"_a, "probe"_a, "refShapeOpts"_a, "probeShapeOpts"_a,
      "refConfId"_a = -1, "probeConfId"_a = -1, "opt_param"_a = 1.0,
      "max_preiters"_a = 10, "max_postiters"_a = 30,
      R"DOC(Aligns a probe molecule to a reference molecule. The probe is modified.

Parameters
----------
ref : RDKit.ROMol
    Reference molecule
probe : RDKit.ROMol
    Probe molecule
refShapeOpts : ShapeInputOptions
    Options for constructing the shape for the reference molecule
probeShapeOpts : ShapeInputOptions
    Options for constructing the shape for the probe molecule
refConfId : int, optional
    Reference conformer ID (default is -1)
probeConfId : int, optional
    Probe conformer ID (default is -1)
opt_param : float, optional
    Balance of shape and color for optimization.
    0 is only color, 0.5 is equal weight, and 1.0 is only shape.
    Default is 1.0.
max_preiters : int, optional
    In the two phase optimization, the maximum iterations done on all poses.
max_postiters : int, optional
    In the two phase optimization, the maximum iterations during the second phase on
    only the best poses from the first phase


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if opt_param is 1.0.)DOC");

  m.def(
      "AlignMol",
      [](const ShapeInput &refShape, RDKit::ROMol &probe, int probeConfId,
         bool useColors, double opt_param, unsigned int max_preiters,
         unsigned int max_postiters, bool applyRefShift) {
        std::vector<float> matrix(12, 0.0);
        auto [nbr_st, nbr_ct] =
            AlignMolecule(refShape, probe, matrix, probeConfId, useColors,
                          opt_param, max_preiters, max_postiters, applyRefShift);
        return nb::make_tuple(nbr_st, nbr_ct);
      },
      "refShape"_a, "probe"_a, "probeConfId"_a = -1, "useColors"_a = true,
      "opt_param"_a = 1.0, "max_preiters"_a = 10, "max_postiters"_a = 30,
      "applyRefShift"_a = false,
      R"DOC(Aligns a probe molecule to a reference shape. The probe is modified.
Assumes the shapes are both centred on the origin.

Parameters
----------
refShape : ShapeInput
    Reference shape
probe : RDKit.ROMol
    Probe molecule
probeConfId : int, optional
    Probe conformer ID (default is -1)
useColors : bool, optional
    Whether or not to use colors in the scoring (default is True)
opt_param : float, optional
    Balance of shape and color for optimization.
    0 is only color, 0.5 is equal weight, and 1.0 is only shape.
    Default is 1.0.
max_preiters : int, optional
    In the two phase optimization, the maximum iterations done on all poses.
max_postiters : int, optional
    In the two phase optimization, the maximum iterations during the second phase on
    only the best poses from the first phase
applyRefShift : bool, optional
    If True, apply the reference shape's shift translation to the final
    coordinates.


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if opt_param is 1.0.)DOC");

  m.def(
      "AlignShapes",
      [](const ShapeInput &refShape, ShapeInput &fitShape, double opt_param,
         unsigned int max_preiters, unsigned int max_postiters) {
        std::vector<float> matrix(12, 0.0);
        auto [nbr_st, nbr_ct] = AlignShape(refShape, fitShape, matrix,
                                           opt_param, max_preiters, max_postiters);
        nb::list pyMatrix;
        for (auto v : matrix) {
          pyMatrix.append(v);
        }
        return nb::make_tuple(nbr_st, nbr_ct, pyMatrix);
      },
      "refShape"_a, "probeShape"_a, "opt_param"_a = 1.0,
      "max_preiters"_a = 10, "max_postiters"_a = 30,
      R"DOC(Aligns a probe shape to a reference shape. The probe is modified.

Parameters
----------
refShape : ShapeInput
    Reference shape
probeShape : ShapeInput
    Probe shape
opt_param : float, optional
    Balance of shape and color for optimization.
    0 is only color, 0.5 is equal weight, and 1.0 is only shape
max_preiters : int, optional
    In the two phase optimization, the maximum iterations done on all poses.
max_postiters : int, optional
    In the two phase optimization, the maximum iterations during the second phase on
    only the best poses from the first phase


Returns
-------
 3-tuple of double, double, list of doubles
    The results are (shape_score, color_score, matrix)
    The matrix is a 12-float list giving the transformation matrix that
    overlays the probe onto the reference.)DOC");

  m.def(
      "TransformConformer",
      [](nb::object pyFinalTrans, nb::object pyFinalRot, nb::object pyMatrix,
         ShapeInput probeShape, RDKit::Conformer &probeConf) {
        std::vector<float> matrix;
        pythonObjectToVect<float>(pyMatrix, matrix);
        if (matrix.size() != 12) {
          throw nb::value_error(
              ("The transformation matrix must have 12 values.  It had " +
               std::to_string(matrix.size()) + ".")
                  .c_str());
        }
        std::vector<double> finalTrans;
        pythonObjectToVect<double>(pyFinalTrans, finalTrans);
        if (finalTrans.size() != 3) {
          throw nb::value_error(
              ("The final translation vector must have 3 values.  It had " +
               std::to_string(finalTrans.size()) + ".")
                  .c_str());
        }
        std::vector<double> finalRot;
        pythonObjectToVect<double>(pyFinalRot, finalRot);
        if (finalRot.size() != 9) {
          throw nb::value_error(
              ("The final rotation vector must have 9 values.  It had " +
               std::to_string(finalRot.size()) + ".")
                  .c_str());
        }
        TransformConformer(finalTrans, finalRot, matrix, probeShape, probeConf);
      },
      "finalTrans"_a, "finalRot"_a, "matrix"_a, "probeShape"_a, "probeConformer"_a,
      R"DOC(Assuming that probeShape has been overlaid onto refShape to give
the supplied transformation matrix, applies that transformation to the
 given conformer.

Parameters
----------
finalTrans : list[float * 3]
    The final translation to apply to conformer.
matrix: list[float * 12]
    The transformation matrix
probeShape : ShapeInput
    Probe shape
probeConformer : Conformer
    Probe conformer
)DOC");

  m.def(
      "PrepareConformer",
      [](const RDKit::ROMol &mol, int confId, nb::object py_opts) {
        ShapeInputOptions opts;
        if (!py_opts.is_none()) {
          opts = nb::cast<ShapeInputOptions>(py_opts);
        }
        return new ShapeInput(PrepareConformer(mol, confId, opts));
      },
      "mol"_a, "confId"_a = -1, "opts"_a = nb::none(),
      R"DOC(Generates a ShapeInput object for a molecule

Parameters
----------
mol : RDKit.ROMol
    Reference molecule
confId : int, optional
    Conformer ID to use (default is -1)
opts : ShapeInputOptions, optional
    Options for Shapeinput

Returns
-------
 a ShapeInput for the molecule)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "ScoreMol",
      [](const RDKit::ROMol &mol1, RDKit::ROMol &mol2,
         const ShapeInputOptions &mol1ShapeOpts,
         const ShapeInputOptions &mol2ShapeOpts, int mol1ConfId,
         int mol2ConfId) {
        auto [nbr_st, nbr_ct] = ScoreMolecule(mol1, mol2, mol1ShapeOpts,
                                              mol2ShapeOpts, mol1ConfId, mol2ConfId);
        return nb::make_tuple(nbr_st, nbr_ct);
      },
      "mol1"_a, "mol2"_a, "mol1ShapeOpts"_a, "mol2ShapeOpts"_a,
      "mol1ConfId"_a = -1, "mol2ConfId"_a = -1,
      R"DOC(Calculate the scores between a shape and a molecule without moving them.

Parameters
----------
mol1 : RDKit.ROMol
    First molecule
mol2 : RDKit.ROMol
    Second molecule
mol1ShapeOptions:
    Options for constructing the shape for molecule 1
mol2ShapeOptions:
    Options for constructing the shape for molecule 2
mol1ConfId : int, optional
    First molecule conformer ID (default is -1)
mol2ConfId : int, optional
    Second conformer ID (default is -1)


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if useColors is False for either of the
shape options)DOC");

  m.def(
      "ScoreMol",
      [](const ShapeInput &shape, RDKit::ROMol &mol,
         const ShapeInputOptions &molShapeOpts, int molConfId) {
        auto [nbr_st, nbr_ct] =
            ScoreMolecule(shape, mol, molShapeOpts, molConfId);
        return nb::make_tuple(nbr_st, nbr_ct);
      },
      "shape"_a, "mol"_a, "molShapeOpts"_a, "molConfId"_a = -1,
      R"DOC(Calculate the scores between 2 molecules without moving them.

Parameters
----------
shape : ShapeInput
    Shape
mol : RDKit.ROMol
    Molecule
molShapeOptions:
    Options for constructing the shape for molecule
molConfId : int, optional
    Molecule conformer ID (default is -1)

Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if shape.useColors is False)DOC");

  m.def(
      "ScoreShape",
      [](const ShapeInput &shape1, ShapeInput &shape2, bool useColors) {
        auto [nbr_st, nbr_ct] = ScoreShape(shape1, shape2, useColors);
        return nb::make_tuple(nbr_st, nbr_ct);
      },
      "shape1"_a, "shape2"_a, "useColors"_a = false,
      R"DOC(Calculate the scores between 2 shapes without moving them.

Parameters
----------
shape1 : ShapeInput
    Shape
shape2 : ShapeInput
    Shape
useColors : bool
    Whether to use colors for the score or not.
Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if useColors is False)DOC");
}
