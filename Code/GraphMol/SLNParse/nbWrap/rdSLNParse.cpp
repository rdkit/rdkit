// $Id$
//
//  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior
//       written permission.
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
// Created by Greg Landrum, September 2006
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/SLNParse/SLNParse.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/FileParseException.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
static ROMol *MolFromSLN(std::string sln, bool sanitize = true,
                         bool debugParser = false) {
  RWMol *newM = SLNToMol(sln, sanitize, debugParser);
  return static_cast<ROMol *>(newM);
}
static ROMol *MolFromQuerySLN(std::string sln, bool mergeHs = true,
                               bool debugParser = false) {
  RWMol *newM = SLNQueryToMol(sln, mergeHs, debugParser);
  return static_cast<ROMol *>(newM);
}
}  // namespace RDKit

NB_MODULE(rdSLNParse, m) {
  m.doc() =
      "Module containing classes and functions for working with Sybyl line "
      "notation (SLN).";

  nb::register_exception_translator(
      [](const std::exception_ptr &p, void *) {
        try {
          std::rethrow_exception(p);
        } catch (const RDKit::SLNParseException &e) {
          std::string msg = std::string("SLNParseException: ") + e.what();
          PyErr_SetString(PyExc_ValueError, msg.c_str());
        }
      });

  m.def("MolFromSLN", &RDKit::MolFromSLN, "SLN"_a, "sanitize"_a = true,
        "debugParser"_a = false, nb::rv_policy::take_ownership,
        R"DOC(Construct a molecule from an SLN string.

    ARGUMENTS:

    - SLN: the SLN string

    - sanitize: (optional) toggles sanitization of the molecule.
      Defaults to True.

  RETURNS:

    a Mol object, None on failure.

  NOTE: the SLN should not contain query information or properties. To build a
    query from SLN, use MolFromQuerySLN.
)DOC");

  m.def("MolFromQuerySLN", &RDKit::MolFromQuerySLN, "SLN"_a,
        "mergeHs"_a = true, "debugParser"_a = false,
        nb::rv_policy::take_ownership,
        R"DOC(Construct a query molecule from an SLN string.

  ARGUMENTS:

    - SLN: the SLN string

    - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
      heavy atoms. Defaults to False.

  RETURNS:

    a Mol object suitable for using in substructure queries, None on failure.
)DOC");
}
