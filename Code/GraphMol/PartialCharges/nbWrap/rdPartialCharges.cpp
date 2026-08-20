//
//  Copyright (C) 2003-2026 Rational Discovery LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(rdPartialCharges, m) {
  m.doc() =
      "Module containing functions to set partial charges - currently "
      "Gasteiger Charges";

  m.def(
      "ComputeGasteigerCharges",
      [](const RDKit::ROMol &mol, int nIter, bool throwOnParamFailure) {
        RDKit::computeGasteigerCharges(&mol, nIter, throwOnParamFailure);
      },
      "mol"_a, "nIter"_a = 12, "throwOnParamFailure"_a = false,
      R"DOC(Compute Gasteiger partial charges for molecule

The charges are computed using an iterative procedure presented in

Ref : J.Gasteiger, M. Marseli, Iterative Equalization of Oribital Electronegatiity
A Rapid Access to Atomic Charges, Tetrahedron Vol 36 p3219 1980

The computed charges are stored on each atom are stored a computed property (under the name
_GasteigerCharge). In addition, each atom also stored the total charge for the implicit hydrogens
on the atom (under the property name _GasteigerHCharge)

ARGUMENTS:

   - mol : the molecule of interrest
   - nIter : number of iteration (defaults to 12)
   - throwOnParamFailure : toggles whether or not an exception should be raised if parameters
     for an atom cannot be found.  If this is false (the default), all parameters for unknown
     atoms will be set to zero.  This has the effect of removing that atom from the iteration.

)DOC");
}
