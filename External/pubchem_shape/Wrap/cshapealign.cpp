/*******************************************************************************

Copyright 2024 by Greg Landrum and the pubchem_shape contributors

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

#include <boost/python.hpp>

#include <vector>

#include "../PubChemShape.hpp"
#include <GraphMol/RDKitBase.h>
#include <RDBoost/Wrap.h>
namespace python = boost::python;

namespace helpers {
python::tuple alignMol(const RDKit::ROMol &ref, RDKit::ROMol &probe,
                       int refConfId, int probeConfId, bool useColors,
                       double opt_param, unsigned int max_preiters,
                       unsigned int max_postiters) {
  std::vector<float> matrix(12, 0.0);
  auto [nbr_st, nbr_ct] =
      AlignMolecule(ref, probe, matrix, refConfId, probeConfId, useColors,
                    opt_param, max_preiters, max_postiters);
  return python::make_tuple(nbr_st, nbr_ct);
}
python::tuple alignMol2(const ShapeInput &ref, RDKit::ROMol &probe,
                        int probeConfId, bool useColors, double opt_param,
                        unsigned int max_preiters, unsigned int max_postiters,
                        bool applyRefShift) {
  std::vector<float> matrix(12, 0.0);
  auto [nbr_st, nbr_ct] =
      AlignMolecule(ref, probe, matrix, probeConfId, useColors, opt_param,
                    max_preiters, max_postiters, applyRefShift);
  return python::make_tuple(nbr_st, nbr_ct);
}
python::tuple alignShapes(const ShapeInput &refShape, ShapeInput &fitShape,
                          double opt_param, unsigned int max_preiters,
                          unsigned int max_postiters) {
  std::vector<float> matrix(12, 0.0);
  auto [nbr_st, nbr_ct] = AlignShape(refShape, fitShape, matrix, opt_param,
                                     max_preiters, max_postiters);
  python::list pyMatrix;
  for (auto m : matrix) {
    pyMatrix.append(m);
  }
  return python::make_tuple(nbr_st, nbr_ct, pyMatrix);
}
void transformConformer(const python::list &pyFinalTrans,
                        const python::list &pyMatrix, ShapeInput probeShape,
                        RDKit::Conformer &probeConf) {
  std::vector<float> matrix;
  pythonObjectToVect<float>(pyMatrix, matrix);
  if (matrix.size() != 12) {
    throw_value_error(
        "The transformation matrix must have 12 values.  It had " +
        std::to_string(matrix.size()) + ".");
  }
  std::vector<double> finalTrans;
  pythonObjectToVect<double>(pyFinalTrans, finalTrans);
  if (finalTrans.size() != 3) {
    throw_value_error(
        "The final translation vector must have 3 values.  It had " +
        std::to_string(finalTrans.size()) + ".");
  }
  TransformConformer(finalTrans, matrix, probeShape, probeConf);
}
ShapeInput *prepConf(const RDKit::ROMol &mol, int confId,
                     const python::object &py_opts) {
  ShapeInputOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<ShapeInputOptions>(py_opts);
  }
  return new ShapeInput(PrepareConformer(mol, confId, opts));
}
void set_atomSubset(ShapeInputOptions &opts, const python::list &as) {
  pythonObjectToVect<unsigned int>(as, opts.atomSubset);
}

python::list get_atomSubset(const ShapeInputOptions &opts) {
  python::list py_list;
  for (const auto &val : opts.atomSubset) {
    py_list.append(val);
  }
  return py_list;
}

void set_notColorAtoms(ShapeInputOptions &opts, const python::list &nca) {
  pythonObjectToVect<unsigned int>(nca, opts.notColorAtoms);
}

python::list get_notColorAtoms(const ShapeInputOptions &opts) {
  python::list py_list;
  for (const auto &val : opts.notColorAtoms) {
    py_list.append(val);
  }
  return py_list;
}

void set_atomRadii(ShapeInputOptions &opts, const python::list &ar) {
  int len = python::len(ar);
  opts.atomRadii.resize(len);
  for (int i = 0; i < len; i++) {
    unsigned int atomIdx = python::extract<unsigned int>(ar[i][0]);
    double radius = python::extract<double>(ar[i][1]);
    opts.atomRadii[i] = std::make_pair(atomIdx, radius);
  }
}

python::list get_atomRadii(const ShapeInputOptions &opts) {
  python::list py_list;
  for (const auto &val : opts.atomRadii) {
    py_list.append(python::make_tuple(static_cast<int>(val.first), val.second));
  }
  return py_list;
}

void set_shapeShift(ShapeInput &shp, const python::list &s) {
  pythonObjectToVect<double>(s, shp.shift);
}
python::list get_shapeShift(const ShapeInput &shp) {
  python::list py_list;
  for (const auto &val : shp.shift) {
    py_list.append(val);
  }
  return py_list;
}
}  // namespace helpers

void wrap_pubchemshape() {
  RegisterVectorConverter<float>("FloatVector");
  RegisterVectorConverter<double>("DoubleVector");
  RegisterVectorConverter<unsigned int>("UnsignedIntVector");

  python::class_<ShapeInputOptions, boost::noncopyable>("ShapeInputOptions",
                                                        "Shape Input Options")
      .def_readwrite(
          "useColors", &ShapeInputOptions::useColors,
          "Whether to use colors (pharmacophore features) in the score.  Default=True.")
      .def_readwrite(
          "includeDummies", &ShapeInputOptions::includeDummies,
          "Whether to use dummy atoms in the alignment. Default=False.")
      .def_readwrite(
          "dummyRadius", &ShapeInputOptions::dummyRadius,
          "If using dummy atoms in the alignment, what radius to use for them."
          "  Default=2.16 (the radius of Xe).")
      .add_property(
          "atomSubset", &helpers::get_atomSubset, &helpers::set_atomSubset,
          "If not empty, use just these atoms in the molecule to form the ShapeInput object.")
      .add_property(
          "notColorAtoms", &helpers::get_notColorAtoms,
          &helpers::set_notColorAtoms,
          "Any atoms mentioned here by index should not be used in a color feature.")
      .add_property(
          "atomRadii", &helpers::get_atomRadii, &helpers::set_atomRadii,
          "Non-standard radii to use for the atoms specified by their indices"
          " in the molecule.  A list of tuples of [int, float].")
      .def("__setattr__", &safeSetattr);

  python::def(
      "AlignMol", &helpers::alignMol,
      (python::arg("ref"), python::arg("probe"), python::arg("refConfId") = -1,
       python::arg("probeConfId") = -1, python::arg("useColors") = true,
       python::arg("opt_param") = 1.0, python::arg("max_preiters") = 10,
       python::arg("max_postiters") = 30),
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
max_preiters : int, optional
max_postiters : int, optional


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if useColors is False)DOC");

  python::def(
      "AlignMol", &helpers::alignMol2,
      (python::arg("refShape"), python::arg("probe"),
       python::arg("probeConfId") = -1, python::arg("useColors") = true,
       python::arg("opt_param") = 1.0, python::arg("max_preiters") = 10,
       python::arg("max_postiters") = 30, python::arg("applyRefShift") = false),
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
optParam : float, optional
max_preiters : int, optional
max_postiters : int, optional
applyRefShift : bool, optional
    If True, apply the reference shape's shift translation to the final
    coordinates.


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if useColors is False)DOC");

  python::def(
      "AlignShapes", &helpers::alignShapes,
      (python::arg("refShape"), python::arg("probeShape"),
       python::arg("opt_param") = 1.0, python::arg("max_preiters") = 10,
       python::arg("max_postiters") = 30),
      R"DOC(Aligns a probe shape to a reference shape. The probe is modified.

Parameters
----------
refShape : ShapeInput
    Reference shape
probeShape : ShapeInput
    Probe shape
opt_param : float, optional
max_preiters : int, optional
max_postiters : int, optional


Returns
-------
 3-tuple of double, double, list of doubles
    The results are (shape_score, color_score, matrix)
    The matrix is a 12-float list giving the transformation matrix that
    overlays the probe onto the reference.)DOC");

  python::def(
      "TransformConformer", &helpers::transformConformer,
      (python::arg("finalTrans"), python::arg("matrix"),
       python::arg("probeShape"), python::arg("probeConformer")),
      R"DOC(Assuming that probeShape has been overlaid onto refShape to give
the supplied transformation matrix, applies that transformation to the
 given conformer.

Parameters
----------
refShape : ShapeInput
    Reference shape
matrix: list[float * 12]
    The transformation matrix
probeShape : ShapeInput
    Probe shape
probeConformer : Conformer
    Probe conformer
)DOC");

  python::def("PrepareConformer", &helpers::prepConf,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("opts") = python::object()),
              R"DOC(Generates a ShapeInput object for a molecule

Parameters
----------
mol : RDKit.ROMol
    Reference molecule
confId : int, optional
    Conformer ID to use (default is -1)
useColors : bool, optional
    Whether or not to assign chemical features (colors) (default is True)

Returns
-------
 a ShapeInput for the molecule)DOC",
              python::return_value_policy<python::manage_new_object>());
  python::class_<ShapeInput, boost::noncopyable>("ShapeInput", python::no_init)
      .def_readwrite("coord", &ShapeInput::coord)
      .def_readwrite("alpha_vector", &ShapeInput::alpha_vector)
      .def_readwrite("atom_type_vector", &ShapeInput::atom_type_vector)
      .def_readwrite("volumeAtomIndexVector",
                     &ShapeInput::volumeAtomIndexVector)
      .add_property("shift", &helpers::get_shapeShift, &helpers::set_shapeShift,
                    "Translation of centre of shape coordinates to origin.")
      .def_readwrite("sov", &ShapeInput::sov)
      .def_readwrite("sof", &ShapeInput::sof);
}

BOOST_PYTHON_MODULE(rdShapeAlign) { wrap_pubchemshape(); }