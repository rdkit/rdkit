//
//  Copyright (c) 2013-2026 Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
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
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <GraphMol/GraphMol.h>

#include "../ConformerParser.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(rdConformerParser, m) {
  m.doc() =
      "Module containing functions to read conformations of a molecule from MD trajectories";

  m.def(
      "AddConformersFromAmberTrajectory",
      [](RDKit::ROMol &mol, std::string fName, int numConfs, bool clearConfs) {
        if (clearConfs) {
          mol.clearConformers();
        }
        std::vector<std::vector<double>> coords;
        RDKit::ConformerParser::readAmberTrajectory(fName, coords,
                                                    mol.getNumAtoms());
        INT_VECT res =
            RDKit::ConformerParser::addConformersFromList(mol, coords, numConfs);
        if (numConfs < 0) {
          numConfs = coords.size();
        }
        return res;
      },
      "mol"_a, "traj"_a, "numConfs"_a = -1, "clearConfs"_a = true,
      R"DOC(Read conformations of a molecule from
an Amber trajectory

ARGUMENTS:

   - mol : the molecule of interest
   - traj : the filename of the trajectory
   - numConfs : number of conformations to read
               The default (-1) reads all.
   - clearConfs : clear all existing conformations on the molecule
                  The default is true.

RETURNS:

   IDs of the new conformations added to the molecule
)DOC");
}
