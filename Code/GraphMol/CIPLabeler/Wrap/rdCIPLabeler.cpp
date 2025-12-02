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
#include <string>

#include <RDBoost/Wrap.h>
#include <RDBoost/python.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include "RDGeneral/ControlCHandler.h"

namespace python = boost::python;
using RDKit::CIPLabeler::assignCIPLabels;

void rdMaxIterationsExceededTranslator(
    RDKit::CIPLabeler::MaxIterationsExceeded const &x) {
  std::ostringstream ss;
  ss << x.what();
  PyErr_SetString(PyExc_RuntimeError, ss.str().c_str());
}

void assignCIPLabelsWrapHelper(RDKit::ROMol &mol,
                               const python::object &atomsToLabel,
                               const python::object &bondsToLabel,
                               unsigned int maxRecursiveIterations) {
  auto atoms = pythonObjectToDynBitset(atomsToLabel, mol.getNumAtoms());
  auto bonds = pythonObjectToDynBitset(bondsToLabel, mol.getNumBonds());

  // If both atoms and bonds are None, assign all the mol.
  if (!atomsToLabel && !bondsToLabel) {
    atoms.set();
    bonds.set();
  }

  assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
  if (RDKit::ControlCHandler::getGotSignal()) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "Assign CIP labels cancelled");
    boost::python::throw_error_already_set();
  }
}

BOOST_PYTHON_MODULE(rdCIPLabeler) {
  python::scope().attr("__doc__") =
      "Module containing a function to assign stereochemical labels based "
      "on an accurate CIP rules implementation. This algoritm is a port "
      "of https://github.com/SiMolecule/centres, which was originally "
      "written by John Mayfield. The original algorithm is described in:\n\n"
      "Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, "
      "A., Redkin, D.\nAlgorithmic Analysis of Cahn--Ingold--Prelog Rules of "
      "Stereochemistry:\nProposals for Revised Rules and a Guide for Machine "
      "Implementation.\nJ. Chem. Inf. Model. 2018, 58, 1755-1765.\n";

  python::register_exception_translator<
      RDKit::CIPLabeler::MaxIterationsExceeded>(
      &rdMaxIterationsExceededTranslator);

  std::string docString =
      "New implementation of Stereo assignment using a true CIP ranking.\n"
      "On return:  The molecule to contains CIP flags\n"
      "Errors:  when maxRecursiveIterations is exceeded, throws a "
      "MaxIterationsExceeded error\nARGUMENTS:\n\n"
      " - mol: the molecule\n"
      " - atomsToLabel: (optional) list of atoms to label\n"
      " - bondsToLabel: (optional) list of bonds to label\n"
      " - maxRecursiveIterations: (optional) protects against pseudo-infinite\n"
      "recursion for highly symmetrical structures.\n A value of 1,250,000 take"
      " about 1 second.  Most structures requires less than 10,000"
      "iterations.\n A peptide with MW~3000 took about 100 iterations, and a "
      "20,000 mw protein took about 600 iterations\n(0 = default - no limit)\n";

  python::def(
      "AssignCIPLabels", assignCIPLabelsWrapHelper,
      (python::arg("mol"), python::arg("atomsToLabel") = python::object(),
       python::arg("bondsToLabel") = python::object(),
       python::arg("maxRecursiveIterations") = 0),
      docString.c_str());
}
