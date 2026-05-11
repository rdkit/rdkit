//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include "RDGeneral/ControlCHandler.h"

namespace nb = nanobind;
using namespace nb::literals;
using RDKit::CIPLabeler::assignCIPLabels;

namespace {
boost::dynamic_bitset<> pyObjToBitset(nb::object obj, size_t n) {
  boost::dynamic_bitset<> result(n);
  if (!obj.is_none()) {
    for (nb::handle h : nb::iter(obj)) {
      result.set(nb::cast<size_t>(h));
    }
  }
  return result;
}

void assignCIPLabelsHelper(RDKit::ROMol &mol, nb::object atomsToLabel,
                           nb::object bondsToLabel,
                           unsigned int maxRecursiveIterations) {
  auto atoms = pyObjToBitset(atomsToLabel, mol.getNumAtoms());
  auto bonds = pyObjToBitset(bondsToLabel, mol.getNumBonds());

  // If both atoms and bonds are None, assign all the mol.
  if (atomsToLabel.is_none() && bondsToLabel.is_none()) {
    atoms.set();
    bonds.set();
  }

  assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
  if (RDKit::ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Assign CIP labels cancelled");
    throw nb::python_error();
  }
}
}  // namespace

NB_MODULE(rdCIPLabeler, m) {
  m.doc() =
      "Module containing a function to assign stereochemical labels based "
      "on an accurate CIP rules implementation. This algoritm is a port "
      "of https://github.com/SiMolecule/centres, which was originally "
      "written by John Mayfield. The original algorithm is described in:\n\n"
      "Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, "
      "A., Redkin, D.\nAlgorithmic Analysis of Cahn--Ingold--Prelog Rules of "
      "Stereochemistry:\nProposals for Revised Rules and a Guide for Machine "
      "Implementation.\nJ. Chem. Inf. Model. 2018, 58, 1755-1765.\n";

  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *) {
        try {
          std::rethrow_exception(p);
        } catch (const RDKit::CIPLabeler::MaxIterationsExceeded &e) {
          std::string msg = e.what();
          PyErr_SetString(PyExc_RuntimeError, msg.c_str());
        }
      });

  m.def(
      "AssignCIPLabels", assignCIPLabelsHelper, "mol"_a,
      "atomsToLabel"_a = nb::none(), "bondsToLabel"_a = nb::none(),
      "maxRecursiveIterations"_a = 0u,
      R"DOC(New implementation of Stereo assignment using a true CIP ranking.
On return:  The molecule to contains CIP flags
Errors:  when maxRecursiveIterations is exceeded, throws a MaxIterationsExceeded error
ARGUMENTS:

 - mol: the molecule
 - atomsToLabel: (optional) list of atoms to label
 - bondsToLabel: (optional) list of bonds to label
 - maxRecursiveIterations: (optional) protects against pseudo-infinite
recursion for highly symmetrical structures.
 A value of 1,250,000 take about 1 second.  Most structures requires less than 10,000iterations.
 A peptide with MW~3000 took about 100 iterations, and a 20,000 mw protein took about 600 iterations
(0 = default - no limit)
)DOC");
}
