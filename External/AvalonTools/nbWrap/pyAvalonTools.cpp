//
//  Created by Greg Landrum, July 2008
//  Copyright (C) 2008-2026 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include <RDBoost/boost_shared_ptr.h>

#include <GraphMol/GraphMol.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseIntVect.h>
#include <AvalonTools.h>

extern "C" {
#include "struchk.h"
}

namespace nb = nanobind;
using namespace nb::literals;

namespace {

enum StruChkFlag {
  bad_molecule = BAD_MOLECULE,
  alias_conversion_failed = ALIAS_CONVERSION_FAILED,
  transformed = TRANSFORMED,
  fragments_found = FRAGMENTS_FOUND,
  either_warning = EITHER_WARNING,
  stereo_error = STEREO_ERROR,
  dubious_stereo_removed = DUBIOUS_STEREO_REMOVED,
  atom_clash = ATOM_CLASH,
  atom_check_failed = ATOM_CHECK_FAILED,
  size_check_failed = SIZE_CHECK_FAILED,
  recharged = RECHARGED,
  stereo_forced_bad = STEREO_FORCED_BAD,
  stereo_transformed = STEREO_TRANSFORMED,
  template_transformed = TEMPLATE_TRANSFORMED,
};

enum StruChkResult {
  success = 0,
  bad_set = BAD_SET,
  transformed_set = TRANSFORMED_SET,
};

}  // namespace

NB_MODULE(pyAvalonTools, m) {
  m.doc() = R"DOC(Module containing functionality from the Avalon toolkit.

The functions currently exposed are:
  - GetCanonSmiles()   : return the canonical smiles for a molecule
  - GetAvalonFP()      : return the Avalon fingerprint for a molecule as
                         an RDKit ExplicitBitVector
  - GetAvalonCountFP()      : return the Avalon fingerprint for a molecule as
                              an RDKit SparseIntVector
  - Generate2DCoords() : use the Avalon coordinate generator to create
                         a set of 2D coordinates for a molecule
Each function can be called with either an RDKit molecule or some
molecule data as text (e.g. a SMILES or an MDL mol block).

See the individual docstrings for more information.
)DOC";

  m.def(
      "GetCanonSmiles",
      (std::string(*)(RDKit::ROMol &, int))AvalonTools::getCanonSmiles,
      "mol"_a, "flags"_a = -1,
      R"DOC(returns canonical smiles for an RDKit molecule)DOC");

  m.def(
      "GetCanonSmiles",
      (std::string(*)(const std::string &, bool, int))AvalonTools::getCanonSmiles,
      "molData"_a, "isSmiles"_a, "flags"_a = -1,
      R"DOC(Returns canonical smiles for some molecule data.
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
MDL mol data is assumed.)DOC");

  m.def(
      "GetAvalonFP",
      [](const RDKit::ROMol &mol, unsigned int nBits, bool isQuery,
         bool resetVect, unsigned int bitFlags) {
        auto *res = new ExplicitBitVect(nBits);
        AvalonTools::getAvalonFP(mol, *res, nBits, isQuery, resetVect,
                                 bitFlags);
        return res;
      },
      "mol"_a, "nBits"_a = 512, "isQuery"_a = false, "resetVect"_a = false,
      "bitFlags"_a = AvalonTools::avalonSimilarityBits,
      R"DOC(returns the Avalon fingerprint for an RDKit molecule)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetAvalonFP",
      [](const std::string &data, bool isSmiles, unsigned int nBits,
         bool isQuery, bool resetVect, unsigned int bitFlags) {
        auto *res = new ExplicitBitVect(nBits);
        AvalonTools::getAvalonFP(data, isSmiles, *res, nBits, isQuery,
                                 resetVect, bitFlags);
        return res;
      },
      "molData"_a, "isSmiles"_a, "nBits"_a = 512, "isQuery"_a = false,
      "resetVect"_a = false,
      "bitFlags"_a = AvalonTools::avalonSimilarityBits,
      R"DOC(returns the Avalon fingerprint for some molecule data.
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
MDL mol data is assumed.)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "Generate2DCoords",
      (unsigned int (*)(RDKit::ROMol &, bool))AvalonTools::set2DCoords,
      "mol"_a, "clearConfs"_a = true,
      R"DOC(Generates 2d coordinates for an RDKit molecule)DOC");

  m.def(
      "Generate2DCoords",
      (std::string(*)(const std::string &, bool))AvalonTools::set2DCoords,
      "molData"_a, "isSmiles"_a,
      R"DOC(returns an MDL mol block with 2D coordinates for some molecule data.
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
MDL mol data is assumed.)DOC");

  m.def(
      "GetAvalonFPAsWords",
      [](const RDKit::ROMol &mol, unsigned int nBits, bool isQuery,
         bool resetVect, unsigned int bitFlags) {
        std::vector<uint32_t> words;
        AvalonTools::getAvalonFP(mol, words, nBits, isQuery, resetVect,
                                 bitFlags);
        nb::list res;
        for (auto w : words) {
          res.append(static_cast<unsigned long>(w));
        }
        return res;
      },
      "mol"_a, "nBits"_a = 512, "isQuery"_a = false, "resetVect"_a = false,
      "bitFlags"_a = AvalonTools::avalonSimilarityBits,
      R"DOC(returns the Avalon fingerprint for an RDKit molecule as a list of ints)DOC");

  m.def(
      "GetAvalonCountFP",
      [](const RDKit::ROMol &mol, unsigned int nBits, bool isQuery,
         unsigned int bitFlags) {
        auto *res = new RDKit::SparseIntVect<uint32_t>(nBits);
        AvalonTools::getAvalonCountFP(mol, *res, nBits, isQuery, bitFlags);
        return res;
      },
      "mol"_a, "nBits"_a = 512, "isQuery"_a = false,
      "bitFlags"_a = AvalonTools::avalonSimilarityBits,
      R"DOC(returns the Avalon count fingerprint for an RDKit molecule)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetAvalonCountFP",
      [](const std::string &data, bool isSmiles, unsigned int nBits,
         bool isQuery, unsigned int bitFlags) {
        auto *res = new RDKit::SparseIntVect<uint32_t>(nBits);
        AvalonTools::getAvalonCountFP(data, isSmiles, *res, nBits, isQuery,
                                     bitFlags);
        return res;
      },
      "molData"_a, "isSmiles"_a, "nBits"_a = 512, "isQuery"_a = false,
      "bitFlags"_a = AvalonTools::avalonSimilarityBits,
      R"DOC(returns the Avalon count fingerprint for some molecule data.
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
MDL mol data is assumed.)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetAvalonFPAsWords",
      [](const std::string &data, bool isSmiles, unsigned int nBits,
         bool isQuery, bool resetVect, unsigned int bitFlags) {
        std::vector<uint32_t> words;
        AvalonTools::getAvalonFP(data, isSmiles, words, nBits, isQuery,
                                 resetVect, bitFlags);
        nb::list res;
        for (auto w : words) {
          res.append(static_cast<unsigned long>(w));
        }
        return res;
      },
      "molData"_a, "isSmiles"_a, "nBits"_a = 512, "isQuery"_a = false,
      "resetVect"_a = false,
      "bitFlags"_a = AvalonTools::avalonSimilarityBits,
      R"DOC(returns the Avalon fingerprint for some molecule data as a list of ints.
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
MDL mol data is assumed.)DOC");

  m.def(
      "InitializeCheckMol",
      (int (*)(const std::string &))AvalonTools::initCheckMol,
      "options"_a = "",
      R"DOC(initializes the structure checker.
The argument should contain option lines separated by embedded newlines.An empty string will be used if the argument is omitted.An non-zero error code is returned in case of failure.)DOC");

  m.def("CloseCheckMolFiles", AvalonTools::closeCheckMolFiles,
        R"DOC(close open files used by molecule-checking functions.)DOC");

  m.def(
      "CheckMolecule",
      [](const std::string &data, bool isSmiles) {
        int errs = 0;
        RDKit::ROMOL_SPTR rMol = AvalonTools::checkMol(errs, data, isSmiles);
        return nb::make_tuple(errs, rMol);
      },
      "molstring"_a, "isSmiles"_a,
      R"DOC(check a molecule passed in as a string.
If the isSmiles argument is true, the string should represent the SMILES encoding
of the molecule, otherwise it should be encoded as an MDL molfile.
The first member of the return tuple contains the bit-encoded corrections made to the molecule.
If possible, the molecule (corrected when appropriate) is returned as the second member of
the return tuple. Otherwise, None is returned.)DOC");

  m.def(
      "CheckMolecule",
      [](RDKit::ROMol &mol) {
        int errs = 0;
        RDKit::ROMOL_SPTR rMol = AvalonTools::checkMol(errs, mol);
        return nb::make_tuple(errs, rMol);
      },
      "mol"_a,
      R"DOC(check a molecule passed in as an RDKit molecule.
The first member of the return tuple contains the bit-encoded corrections made to the molecule.
If possible, the molecule (corrected when appropriate) is returned as the second member of
the return tuple. Otherwise, None is returned.)DOC");

  m.def(
      "CheckMoleculeString",
      [](const std::string &data, bool isSmiles) {
        std::pair<std::string, int> res =
            AvalonTools::checkMolString(data, isSmiles);
        return nb::make_tuple(res.second, res.first);
      },
      "molstring"_a, "isSmiles"_a,
      R"DOC(check a molecule passed in as a string and returns the result as a string.
If the isSmiles argument is true, the string should represent the SMILES encoding
of the molecule, otherwise it should be encoded as an MDL molfile.
The first member of the return tuple contains the bit-encoded corrections made to the molecule.
If possible, a corrected CTAB for the molecule is returned as the second member of
the return tuple.)DOC");

  m.def("GetCheckMolLog", AvalonTools::getCheckMolLog,
        R"DOC(Returns the Struchk log for the last molecules processed.)DOC");

  m.attr("avalonSSSBits") = AvalonTools::avalonSSSBits;
  m.attr("avalonSimilarityBits") = AvalonTools::avalonSimilarityBits;

  nb::enum_<StruChkFlag>(m, "StruChkFlag")
      .value("bad_molecule", bad_molecule)
      .value("alias_conversion_failed", alias_conversion_failed)
      .value("transformed", transformed)
      .value("fragments_found", fragments_found)
      .value("either_warning", either_warning)
      .value("stereo_error", stereo_error)
      .value("dubious_stereo_removed", dubious_stereo_removed)
      .value("atom_clash", atom_clash)
      .value("atom_check_failed", atom_check_failed)
      .value("size_check_failed", size_check_failed)
      .value("recharged", recharged)
      .value("stereo_forced_bad", stereo_forced_bad)
      .value("stereo_transformed", stereo_transformed)
      .value("template_transformed", template_transformed);

  nb::enum_<StruChkResult>(m, "StruChkResult")
      .value("success", success)
      .value("bad_set", bad_set)
      .value("transformed_set", transformed_set);
}
