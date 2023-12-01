// $Id$
//
//  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <boost/python.hpp>

#include <GraphMol/GraphMol.h>
#include <RDBoost/Wrap.h>

#include "../ConformerParser.h"

namespace python = boost::python;

namespace RDKit {

INT_VECT AddConformersFromAmberTrajectory(ROMol &mol, std::string fName,
                                          int numConfs, bool clearConfs) {
  if (clearConfs) {
    mol.clearConformers();
  }
  std::vector<std::vector<double>> coords;
  ConformerParser::readAmberTrajectory(fName, coords, mol.getNumAtoms());
  INT_VECT res = ConformerParser::addConformersFromList(mol, coords, numConfs);
  if (numConfs < 0) {
    numConfs = coords.size();
  }
  return res;
}

}  // namespace RDKit

BOOST_PYTHON_MODULE(rdConformerParser) {
  python::scope().attr("__doc__") =
      "Module containing functions to read conformations of a molecule from MD trajectories";

  // import_array();

  std::string docString =
      "Read conformations of a molecule from \n\
 an Amber trajectory\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - traj : the filename of the trajectory \n\
    - numConfs : number of conformations to read \n\
                The default (-1) reads all. \n\
    - clearConfs : clear all existing conformations on the molecule\n\
                   The default is true. \n\
\n\
 RETURNS:\n\n\
    IDs of the new conformations added to the molecule \n\
\n";
  python::def("AddConformersFromAmberTrajectory",
              RDKit::AddConformersFromAmberTrajectory,
              (python::arg("mol"), python::arg("traj"),
               python::arg("numConfs") = -1, python::arg("clearConfs") = true),
              docString.c_str());
}
