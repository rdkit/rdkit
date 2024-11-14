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
                        unsigned int max_preiters, unsigned int max_postiters) {
  std::vector<float> matrix(12, 0.0);
  auto [nbr_st, nbr_ct] =
      AlignMolecule(ref, probe, matrix, probeConfId, useColors, opt_param,
                    max_preiters, max_postiters);
  return python::make_tuple(nbr_st, nbr_ct);
}
ShapeInput *prepConf(const RDKit::ROMol &mol, int confId, bool useColors) {
  return new ShapeInput(PrepareConformer(mol, confId, useColors));
}
}  // namespace helpers

void wrap_pubchemshape() {
  RegisterVectorConverter<float>("FloatVector");
  RegisterVectorConverter<double>("DoubleVector");
  RegisterVectorConverter<unsigned int>("UnsignedIntVector");

  python::def(
      "AlignMol", &helpers::alignMol,
      (python::arg("ref"), python::arg("probe"), python::arg("refConfId") = -1,
       python::arg("probeConfId") = -1, python::arg("useColors") = true,
       python::arg("opt_param") = 0.5, python::arg("max_preiters") = 3,
       python::arg("max_postiters") = 16),
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
optParam : float, optional
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
       python::arg("opt_param") = 0.5, python::arg("max_preiters") = 3,
       python::arg("max_postiters") = 16),
      R"DOC(Aligns a probe molecule to a reference shape. The probe is modified.

Parameters
----------
refShape : ShapeInput
    Reference molecule
probe : RDKit.ROMol
    Probe molecule
probeConfId : int, optional
    Probe conformer ID (default is -1)
useColors : bool, optional
    Whether or not to use colors in the scoring (default is True)
optParam : float, optional
max_preiters : int, optional
max_postiters : int, optional


Returns
-------
 2-tuple of doubles
    The results are (shape_score, color_score)
    The color_score is zero if useColors is False)DOC");

  python::def("PrepareConformer", &helpers::prepConf,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("useColors") = true),
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
      .def_readwrite("shift", &ShapeInput::shift)
      .def_readwrite("sov", &ShapeInput::sov)
      .def_readwrite("sof", &ShapeInput::sof);
}

BOOST_PYTHON_MODULE(rdShapeAlign) { wrap_pubchemshape(); }