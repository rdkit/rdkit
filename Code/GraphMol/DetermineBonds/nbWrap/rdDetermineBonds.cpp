//
//  Copyright (C) 2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/DetermineBonds/DetermineBonds.h>
#include <RDGeneral/ControlCHandler.h>

namespace nb = nanobind;
using namespace nb::literals;
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
  if (ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Determine Bond Orders cancelled");
    throw nb::python_error();
  }
}

void determineBondsHelper(ROMol &mol, bool useHueckel, int charge,
                          double covFactor, bool allowChargedFragments,
                          bool embedChiral, bool useAtomMap, bool useVdw,
                          size_t maxIterations) {
  auto &wmol = static_cast<RWMol &>(mol);
  determineBonds(wmol, useHueckel, charge, covFactor, allowChargedFragments,
                 embedChiral, useAtomMap, useVdw, maxIterations);
  if (ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Determine Bond Orders cancelled");
    throw nb::python_error();
  }
}

bool hueckelSupportEnabled() {
#ifdef RDK_BUILD_YAEHMOP_SUPPORT
  return true;
#else
  return false;
#endif
}
}  // namespace

NB_MODULE(rdDetermineBonds, m) {
  m.doc() =
      "Module containing a C++ implementation of the xyz2mol algorithm. This "
      "is based on xyz2mol: https://github.com/jensengroup/xyz2mol";

  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *) {
        try {
          std::rethrow_exception(p);
        } catch (const RDKit::MaxFindBondOrdersItersExceeded &e) {
          PyErr_SetString(PyExc_RuntimeError, e.what());
        }
      });

  m.def("DetermineConnectivity", &determineConnectivityHelper, "mol"_a,
        "useHueckel"_a = false, "charge"_a = 0, "covFactor"_a = 1.3,
        "useVdw"_a = false,
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
)DOC");

  m.def("DetermineBondOrders", &determineBondOrdersHelper, "mol"_a,
        "charge"_a = 0, "allowChargedFragments"_a = true,
        "embedChiral"_a = true, "useAtomMap"_a = false,
        "maxIterations"_a = size_t(0),
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
   maxIterations: (optional) maximum number of iterations to run in the bond order
   determination algorithm, after which a MaxFindBondOrdersItersExceeded
   exception will be thrown. Defaults to 0 (no limit)
)DOC");

  m.def("DetermineBonds", &determineBondsHelper, "mol"_a,
        "useHueckel"_a = false, "charge"_a = 0, "covFactor"_a = 1.3,
        "allowChargedFragments"_a = true, "embedChiral"_a = true,
        "useAtomMap"_a = false, "useVdw"_a = false,
        "maxIterations"_a = size_t(0),
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
   maxIterations: (optional) maximum number of iterations to run in the bond order
   determination algorithm, after which a MaxFindBondOrdersItersExceeded
   exception will be thrown. Defaults to 0 (no limit)
)DOC");

  m.def("hueckelEnabled", &hueckelSupportEnabled,
        "whether or not the RDKit was compiled with YAeHMOP support");
}
