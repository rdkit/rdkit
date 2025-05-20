//
//  Copyright (C) 2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>

#include <RDBoost/python.h>
#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/DetermineBonds/DetermineBonds.h>

namespace python = boost::python;
using namespace RDKit;

namespace {

void determineConnectivityHelper(ROMol &mol, bool useHueckel, int charge,
                                 double covFactor, bool useVdw) {
  auto &wmol = static_cast<RWMol &>(mol);
  determineConnectivity(wmol, useHueckel, charge, covFactor, useVdw);
}
void determineBondOrdersHelper(ROMol &mol, int charge,
                               bool allowChargedFragments, bool embedChiral,
                               bool useAtomMap, size_t maxIterations) {
  auto &wmol = static_cast<RWMol &>(mol);
  determineBondOrders(wmol, charge, allowChargedFragments, embedChiral,
                      useAtomMap, maxIterations);
}
void determineBondsHelper(ROMol &mol, bool useHueckel, int charge,
                          double covFactor, bool allowChargedFragments,
                          bool embedChiral, bool useAtomMap, bool useVdw,
                          size_t maxIterations) {
  auto &wmol = static_cast<RWMol &>(mol);
  determineBonds(wmol, useHueckel, charge, covFactor, allowChargedFragments,
                 embedChiral, useAtomMap, useVdw, maxIterations);
}
bool hueckelSupportEnabled() {
#ifdef RDK_BUILD_YAEHMOP_SUPPORT
  return true;
#else
  return false;
#endif
}
}  // namespace

void rdMaxFindBondOrdersItersExceededTranslator(
    const RDKit::MaxFindBondOrdersItersExceeded &x) {
  PyErr_SetString(PyExc_RuntimeError, x.what());
}

BOOST_PYTHON_MODULE(rdDetermineBonds) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of the xyz2mol algorithm. This is based on xyz2mol: https://github.com/jensengroup/xyz2mol";

  python::register_exception_translator<RDKit::MaxFindBondOrdersItersExceeded>(
      &rdMaxFindBondOrdersItersExceededTranslator);

  std::string docs;
  docs =
      R"DOC(Assigns atomic connectivity to a molecule using atomic coordinates,
disregarding pre-existing bonds

Args:
   mol : the molecule of interest; it must have a 3D conformer
   useHueckel : (optional) if this is  \c true, extended Hueckel theory
       will be used to determine connectivity rather than the van der Waals
       or connect-the-dots methods
   charge : (optional) the charge of the molecule; it must be provided if
       the Hueckel method is used and charge is non-zero
   covFactor : (optional) the factor with which to multiply each covalent
       radius if the van der Waals method is used
   useVdw: (optional) if this is false, the connect-the-dots method
       will be used instead of the van der Waals method
)DOC";
  python::def("DetermineConnectivity", &determineConnectivityHelper,
              (python::arg("mol"), python::arg("useHueckel") = false,
               python::arg("charge") = 0, python::arg("covFactor") = 1.3,
               python::arg("useVdw") = false),
              docs.c_str());

  docs =
      R"DOC(Assigns atomic connectivity to a molecule using atomic coordinates,
disregarding pre-existing bonds

Args:
   mol : the molecule of interest; it must have a 3D conformer
   charge : (optional) the charge of the molecule; it must be provided if
       the Hueckel method is used and charge is non-zero
   allowChargedFragments : (optional) if this is true, formal charges
       will be placed on atoms according to their valency; otherwise, radical
       electrons will be placed on the atoms
   embedChiral : (optional) if this is true,
      chirality information will be embedded into the molecule; the function calls
      sanitizeMol() when this is true
   useAtomMap : (optional) if this is true, an atom map will be created for the
      molecule
)DOC";
  python::def(
      "DetermineBondOrders", &determineBondOrdersHelper,
      (python::arg("mol"), python::arg("charge") = 0,
       python::arg("allowChargedFragments") = true,
       python::arg("embedChiral") = true, python::arg("useAtomMap") = false,
       python::arg("maxIterations") = 0),
      docs.c_str());

  docs =
      R"DOC(Assigns atomic connectivity to a molecule using atomic coordinates,
disregarding pre-existing bonds

Args:
   mol : the molecule of interest; it must have a 3D conformer
   useHueckel : (optional) if this is true, extended Hueckel theory
       will be used to determine connectivity rather than the van der Waals
       or connect-the-dots methods
   charge : (optional) the charge of the molecule; it must be provided if
       the Hueckel method is used and charge is non-zero
   covFactor : (optional) the factor with which to multiply each covalent
       radius if the van der Waals method is used
   allowChargedFragments : (optional) if this is true, formal charges
       will be placed on atoms according to their valency; otherwise, radical
       electrons will be placed on the atoms
   embedChiral : (optional) if this is true,
      chirality information will be embedded into the molecule; the function calls
      sanitizeMol() when this is true
   useAtomMap : (optional) if this is true, an atom map will be created for the
      molecule
   useVdw: (optional) if this is false, the connect-the-dots method
       will be used instead of the van der Waals method
)DOC";
  python::def(
      "DetermineBonds", &determineBondsHelper,
      (python::arg("mol"), python::arg("useHueckel") = false,
       python::arg("charge") = 0, python::arg("covFactor") = 1.3,
       python ::arg("allowChargedFragments") = true,
       python::arg("embedChiral") = true, python::arg("useAtomMap") = false,
       python::arg("useVdw") = false, python::arg("maxIterations") = 0),
      docs.c_str());
  python::def("hueckelEnabled", &hueckelSupportEnabled,
              "whether or not the RDKit was compiled with YAeHMOP support");
}
