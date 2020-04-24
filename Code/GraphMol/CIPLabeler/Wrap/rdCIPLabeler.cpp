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

namespace python = boost::python;
using RDKit::CIPLabeler::assignCIPLabels;

BOOST_PYTHON_MODULE(rdCIPLabeler) {
  python::scope().attr("__doc__") =
      "Module containing a function to assign stereochemical labels based "
      "on an accurate CIP rules implementation:\n\nHanson, R. M., Musacchio, "
      "S., Mayfield, J. W., Vainio, M. J., Yerin, A., Redkin, D.\nAlgorithmic "
      "Analysis of Cahn−Ingold−Prelog Rules of Stereochemistry:\nProposals "
      "for Revised Rules and a Guide for Machine Implementation.\nJ. Chem. "
      "Inf. Model. 2018, 58, 1755-1765.\n";

  std::string docString =
      "New implementation of Stereo assignment using a true CIP ranking";

  python::def("AssignCIPLabels", assignCIPLabels, python::arg("mol"),
              docString.c_str());
}
