// $Id$
//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/FMCS/FMCS.h>

namespace python = boost::python;

namespace RDKit {
}

BOOST_PYTHON_MODULE(rdFMCS) {
  
  python::scope().attr("__doc__") =
    "Module containing a C++ implementation of the FMCS algorithm";
   
  /*
  std::string docString = "Compute the centroid of the conformation - hydrogens are ignored and no attention\n\
                           if paid to the difference in sizes of the heavy atoms\n";
  python::def("ComputeCentroid", MolTransforms::computeCentroid,
              (python::arg("conf"), python::arg("ignoreHs")=true),
              docString.c_str());
  */
}
