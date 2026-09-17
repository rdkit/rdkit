//
//  Copyright (c) 2025-2026 Glysade Inc. and other RDKit contributors
//
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
//     * Neither the name of Novartis Institutues for BioMedical Research Inc.
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
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <ChemDraw/chemdraw.h>
#include <ChemDraw/chemdrawreaction.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

std::string pyObjectToString(nb::object input) {
  if (nb::isinstance<nb::str>(input)) {
    return nb::cast<std::string>(input);
  }
  std::wstring ws = nb::cast<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

nb::tuple MolsFromChemDrawBlockHelper(
    const std::string &block, bool sanitize, bool removeHs,
    RDKit::v2::NeedsCleanPolicy needsCleanPolicy =
        RDKit::v2::NeedsCleanPolicy::TrustSource) {
  std::vector<std::unique_ptr<RWMol>> mols;
  try {
    mols = RDKit::v2::MolsFromChemDrawBlock(
        block,
        {sanitize, removeHs, RDKit::v2::CDXFormat::CDXML, needsCleanPolicy});
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    std::shared_ptr<ROMol> sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::tuple MolsFromChemDrawFileHelper(
    nb::object cdxml, bool sanitize, bool removeHs,
    RDKit::v2::NeedsCleanPolicy needsCleanPolicy =
        RDKit::v2::NeedsCleanPolicy::TrustSource) {
  auto mols = RDKit::v2::MolsFromChemDrawFile(
      pyObjectToString(cdxml),
      {sanitize, removeHs, RDKit::v2::CDXFormat::CDXML, needsCleanPolicy});
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    std::shared_ptr<ROMol> sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::tuple ReactionsFromChemDrawFileHelper(const std::string &filename,
                                          bool sanitize, bool removeHs) {
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  try {
    rxns = RDKit::v2::ChemDrawFileToChemicalReactions(filename, sanitize,
                                                      removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &rxn : rxns) {
    // take ownership of the data from the unique_ptr
    res.append(std::shared_ptr<ChemicalReaction>(rxn.release()));
  }
  return nb::tuple(res);
}

nb::tuple ReactionsFromChemDrawBlockHelper(nb::object imolBlock, bool sanitize,
                                           bool removeHs) {
  std::istringstream inStream(pyObjectToString(imolBlock));
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  try {
    rxns = RDKit::v2::ChemDrawDataStreamToChemicalReactions(inStream, sanitize,
                                                            removeHs);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &rxn : rxns) {
    // take ownership of the data from the unique_ptr
    res.append(std::shared_ptr<ChemicalReaction>(rxn.release()));
  }
  return nb::tuple(res);
}

}  // namespace

NB_MODULE(rdChemDraw, m) {
  m.doc() =
      "Module containing classes and functions for working with ChemDraw files.";

  nb::enum_<v2::CDXFormat>(m, "CDXFormat")
      .value("CDX", v2::CDXFormat::CDX)
      .value("CDXML", v2::CDXFormat::CDXML);

  nb::enum_<v2::NeedsCleanPolicy>(m, "NeedsCleanPolicy")
      .value("TrustSource", v2::NeedsCleanPolicy::TrustSource)
      .value("TrustExplicitHydrogens",
             v2::NeedsCleanPolicy::TrustExplicitHydrogens);

  m.def("MolsFromChemDrawFile", MolsFromChemDrawFileHelper, "filename"_a,
        "sanitize"_a = true, "removeHs"_a = true,
        "needsCleanPolicy"_a = v2::NeedsCleanPolicy::TrustSource,
        R"DOC(Extract all molecules from a ChemDraw file.

Note that the ChemDraw format is large and complex, the RDKit doesn't support
full functionality, just the base ones required for molecule and
reaction parsing.

ARGUMENTS:

  - filename: the chemdraw filename (.cdx/.cdxml)

  - sanitize: if True, sanitize the molecules [default True]

  - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

  - needsCleanPolicy: how to handle `NeedsClean` hydrogen metadata.
  `TrustSource` honors `NeedsClean` by allowing sanitization to
  recompute hydrogens. `TrustExplicitHydrogens` preserves the literal
  source metadata when sanitize is True. [default TrustSource]

RETURNS:
  a tuple of parsed Mol objects.)DOC");

  m.def("MolsFromChemDrawBlock", MolsFromChemDrawBlockHelper, "block"_a,
        "sanitize"_a = true, "removeHs"_a = true,
        "needsCleanPolicy"_a = v2::NeedsCleanPolicy::TrustSource,
        R"DOC(Extract all molecules from a ChemDraw block.

Note that the ChemDraw format is large and complex, the RDKit doesn't support
full functionality, just the base ones required for molecule and
reaction parsing.

ARGUMENTS:

  - block: the CDX/CDXML block

  - sanitize: if True, sanitize the molecules [default True]

  - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

  - needsCleanPolicy: how to handle `NeedsClean` hydrogen metadata.
  `TrustSource` honors `NeedsClean` by allowing sanitization to
  recompute hydrogens. `TrustExplicitHydrogens` preserves the literal
  source metadata when sanitize is True. [default TrustSource]

RETURNS:
  a tuple of parsed Mol objects.)DOC");

  m.def("ReactionsFromChemDrawFile", ReactionsFromChemDrawFileHelper,
        "filename"_a, "sanitize"_a = false, "removeHs"_a = false,
        R"DOC(Extract all reactions from a ChemDraw file.

Note that the ChemDraw format is large and complex, the RDKit doesn't support
full functionality, just the base ones required for molecule and
reaction parsing.

ARGUMENTS:

  - filename: the chemdraw filename (.cdx/.cdxml)

  - sanitize: if True, sanitize the molecules [default True]

  - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

RETURNS:
  a tuple of parsed ChemicalReaction objects.)DOC");

  m.def("ReactionsFromChemDrawBlock", ReactionsFromChemDrawBlockHelper,
        "rxnblock"_a, "sanitize"_a = false, "removeHs"_a = false,
        R"DOC(Extract all reactions from a ChemDraw text block.

Note that the ChemDraw format is large and complex, the RDKit doesn't support
full functionality, just the base ones required for molecule and
reaction parsing.

ARGUMENTS:

  - rxnblock: the ChemDraw text block

  - sanitize: if True, sanitize the molecules [default True]

  - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

RETURNS:
  a tuple of parsed ChemicalReaction objects.)DOC");

  m.def(
      "MolToChemDrawBlock", v2::MolToChemDrawBlock, "mol"_a,
      "format"_a = v2::CDXFormat::CDXML,
      R"DOC(Convert a molecule into a chemdraw string using the specified format

ARGUMENTS:

  - mol: the molecule to convert

  - format: The ChemDraw format to use, CDXML/CDX [default CDXML]

RETURNS:
  the ChemDraw string.)DOC");
}
