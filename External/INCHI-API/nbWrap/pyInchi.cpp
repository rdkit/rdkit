//
//  Copyright (c) 2011-2026 Novartis Institutes for BioMedical Research Inc.
//    and other RDKit contributors
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
#include <nanobind/stl/tuple.h>

#include <GraphMol/GraphMol.h>
#include <RDBoost/boost_shared_ptr.h>
#include "../inchi.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(rdinchi, m) {
  m.def(
      "InchiToMol",
      [](const std::string &inchi, bool sanitize, bool removeHs) {
        RDKit::ExtraInchiReturnValues rv;
        RDKit::ROMol *mol = RDKit::InchiToMol(inchi, rv, sanitize, removeHs);
        if (mol == nullptr) {
          return std::make_tuple(RDKit::ROMOL_SPTR(), rv.returnCode,
                                 rv.messagePtr, rv.logPtr);
        } else {
          return std::make_tuple(RDKit::ROMOL_SPTR(mol), rv.returnCode,
                                 rv.messagePtr, rv.logPtr);
        }
      },
      "inchi"_a, "sanitize"_a = true, "removeHs"_a = true,
      R"DOC(return a ROMol for a InChI string
  Returns:
    a tuple with:
      - the molecule
      - the return code from the InChI conversion
      - a string with any messages from the InChI conversion
      - a string with any log messages from the InChI conversion)DOC");

  m.def(
      "MolToInchi",
      [](const RDKit::ROMol &mol, std::string options) {
        RDKit::ExtraInchiReturnValues rv;
        const char *_options = nullptr;
        if (options.size()) {
          _options = options.c_str();
        }
        std::string inchi = RDKit::MolToInchi(mol, rv, _options);
        return std::make_tuple(inchi, rv.returnCode, rv.messagePtr, rv.logPtr,
                               rv.auxInfoPtr);
      },
      "mol"_a, "options"_a = std::string(),
      R"DOC(return the InChI for a ROMol molecule.

  Arguments:
    - mol: the molecule to use.
    - options: the InChI generation options.
      Options should be prefixed with either a - or a /
      Available options are explained in the InChI technical FAQ:
      https://www.inchi-trust.org/technical-faq/#15.14
      and the User Guide:
      https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
  Returns:
    a tuple with:
      - the InChI
      - the return code from the InChI conversion
      - a string with any messages from the InChI conversion
      - a string with any log messages from the InChI conversion
      - a string with the InChI AuxInfo)DOC");

  m.def(
      "MolBlockToInchi",
      [](const std::string &molblock, std::string options) {
        RDKit::ExtraInchiReturnValues rv;
        const char *_options = nullptr;
        if (options.size()) {
          _options = options.c_str();
        }
        std::string inchi = RDKit::MolBlockToInchi(molblock, rv, _options);
        return std::make_tuple(inchi, rv.returnCode, rv.messagePtr, rv.logPtr,
                               rv.auxInfoPtr);
      },
      "molblock"_a, "options"_a = std::string(),
      R"DOC(return the InChI for a ROMol molecule.

  Arguments:
    - molblock: the mol block to use.
    - options: the InChI generation options.
      Options should be prefixed with either a - or a /
      Available options are explained in the InChI technical FAQ:
      https://www.inchi-trust.org/technical-faq/#15.14
      and the User Guide:
      https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
  Returns:
    a tuple with:
      - the InChI
      - the return code from the InChI conversion
      - a string with any messages from the InChI conversion
      - a string with any log messages from the InChI conversion
      - a string with the InChI AuxInfo)DOC");

  m.def("InchiToInchiKey", RDKit::InchiToInchiKey, "inchi"_a,
        "return the InChI key for an InChI string");

  m.def(
      "MolToInchiKey", RDKit::MolToInchiKey, "mol"_a,
      "options"_a = std::string(),
      R"DOC(return the InChI key for a ROMol molecule.

  Arguments:
    - mol: the molecule to use.
    - options: the InChI generation options.
      Options should be prefixed with either a - or a /
      Available options are explained in the InChI technical FAQ:
      https://www.inchi-trust.org/technical-faq/#15.14
      and the User Guide available from:
      https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
  Returns: the InChI key)DOC");

  m.def("GetInchiVersion", RDKit::getInchiVersion,
        "returns the version of the InChI software being used");
}
