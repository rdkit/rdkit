//  Copyright (c) 2017-2026, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
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
#include <nanobind/stl/shared_ptr.h>
#include <chrono>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/RGroupDecomposition/RGroupUtils.h>
#include <RDBoost/Wrap_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

class RGroupDecompositionHelper {
  std::unique_ptr<RGroupDecomposition> decomp;

 public:
  RGroupDecompositionHelper(nb::object cores,
                            const RGroupDecompositionParameters &params =
                                RGroupDecompositionParameters()) {
    try {
      const ROMol &mol = nb::cast<const ROMol &>(cores);
      decomp.reset(new RGroupDecomposition(mol, params));
    } catch (const nb::cast_error &) {
      MOL_SPTR_VECT coreMols;
      for (nb::handle h : nb::iter(cores)) {
        auto *mol_ptr = nb::cast<ROMol *>(h);
        if (!mol_ptr) {
          throw nb::value_error("reaction called with None reactants");
        }
        ROMOL_SPTR sptr(mol_ptr, [](ROMol *) {});
        coreMols.push_back(sptr);
      }
      decomp.reset(new RGroupDecomposition(coreMols, params));
    }
  }

  int Add(const ROMol &mol) {
    NOGIL gil;
    return decomp->add(mol);
  }

  int GetMatchingCoreIdx(const ROMol &mol, nb::object matches) {
    std::vector<MatchVectType> matchVect;
    int coreIdx;
    {
      NOGIL gil;
      coreIdx = decomp->getMatchingCoreIdx(mol, &matchVect);
    }
    if (!matches.is_none() && nb::isinstance<nb::list>(matches)) {
      for (const auto &match : matchVect) {
        nb::list atomMap;
        for (const auto &pair : match) {
          atomMap.append(nb::make_tuple(pair.first, pair.second));
        }
        matches.attr("append")(nb::tuple(atomMap));
      }
    }
    return coreIdx;
  }

  bool Process() {
    NOGIL gil;
    return decomp->process();
  }

  nb::tuple ProcessAndScore() {
    NOGIL gil;
    auto result = decomp->processAndScore();
    return nb::make_tuple(result.success, result.score);
  }

  nb::list GetRGroupLabels() {
    nb::list result;
    std::vector<std::string> labels = decomp->getRGroupLabels();
    for (const auto &label : labels) {
      result.append(label);
    }
    return result;
  }

  nb::list GetRGroupsAsRows(bool asSmiles = false) {
    const RGroupRows &groups = decomp->getRGroupsAsRows();
    nb::list result;

    for (const auto &side_chains : groups) {
      nb::dict dict;
      for (const auto &[lbl, mol] : side_chains) {
        if (asSmiles) {
          dict[nb::cast(lbl)] = nb::cast(MolToSmiles(*mol, true));
        } else {
          dict[nb::cast(lbl)] = nb::cast(toStd(mol));
        }
      }
      result.append(dict);
    }
    return result;
  }

  nb::dict GetRGroupsAsColumns(bool asSmiles = false) {
    nb::dict result;

    RGroupColumns groups = decomp->getRGroupsAsColumns();

    for (RGroupColumns::const_iterator it = groups.begin(); it != groups.end();
         ++it) {
      nb::list col;

      for (const auto &cit : it->second) {
        if (asSmiles) {
          col.append(nb::cast(MolToSmiles(*cit, true)));
        } else {
          col.append(nb::cast(toStd(cit)));
        }
      }
      result[nb::cast(it->first)] = col;
    }
    return result;
  }
};

nb::object RGroupDecomp(nb::object cores, nb::object mols,
                        bool asSmiles = false, bool asRows = true,
                        const RGroupDecompositionParameters &options =
                            RGroupDecompositionParameters()) {
  auto t0 = std::chrono::steady_clock::now();
  RGroupDecompositionHelper decomp(cores, options);
  nb::list unmatched;

  unsigned int idx = 0;
  for (nb::handle h : nb::iter(mols)) {
    auto *mol_ptr = nb::cast<ROMol *>(h);
    if (!mol_ptr) {
      throw nb::value_error("reaction called with None reactants");
    }
    ROMOL_SPTR sptr(mol_ptr, [](ROMol *) {});
    if (decomp.Add(*sptr) == -1) {
      unmatched.append(idx);
    }
    ++idx;
    checkForTimeout(t0, options.timeout);
  }

  decomp.Process();
  if (asRows) {
    return nb::cast(
        nb::make_tuple(decomp.GetRGroupsAsRows(asSmiles), unmatched));
  } else {
    return nb::cast(
        nb::make_tuple(decomp.GetRGroupsAsColumns(asSmiles), unmatched));
  }
}

void relabelMappedDummiesHelper(ROMol &mol, unsigned int inputLabels,
                                unsigned int outputLabels) {
  relabelMappedDummies(mol, static_cast<RGroupLabelling>(inputLabels),
                       static_cast<RGroupLabelling>(outputLabels));
}

}  // namespace RDKit

NB_MODULE(rdRGroupDecomposition, m) {
  m.doc() =
      R"DOC(TEST!!! Module containing RGroupDecomposition classes and functions.)DOC";

  nb::enum_<RDKit::RGroupLabels>(m, "RGroupLabels", nb::is_arithmetic())
      .value("IsotopeLabels", RDKit::IsotopeLabels)
      .value("AtomMapLabels", RDKit::AtomMapLabels)
      .value("AtomIndexLabels", RDKit::AtomIndexLabels)
      .value("RelabelDuplicateLabels", RDKit::RelabelDuplicateLabels)
      .value("MDLRGroupLabels", RDKit::MDLRGroupLabels)
      .value("DummyAtomLabels", RDKit::DummyAtomLabels)
      .value("AutoDetect", RDKit::AutoDetect);

  nb::enum_<RDKit::RGroupMatching>(m, "RGroupMatching", nb::is_arithmetic())
      .value("Greedy", RDKit::Greedy)
      .value("GreedyChunks", RDKit::GreedyChunks)
      .value("Exhaustive", RDKit::Exhaustive)
      .value("NoSymmetrization", RDKit::NoSymmetrization)
      .value("GA", RDKit::GA);

  nb::enum_<RDKit::RGroupLabelling>(m, "RGroupLabelling", nb::is_arithmetic())
      .value("AtomMap", RDKit::AtomMap)
      .value("Isotope", RDKit::Isotope)
      .value("MDLRGroup", RDKit::MDLRGroup);

  nb::enum_<RDKit::RGroupCoreAlignment>(m, "RGroupCoreAlignment",
                                        nb::is_arithmetic())
      // DEPRECATED, remove the following line in release 2021.03
      .value("None", RDKit::NoAlignment)
      .value("NoAlignment", RDKit::NoAlignment)
      .value("MCS", RDKit::MCS);

  nb::enum_<RDKit::RGroupScore>(m, "RGroupScore", nb::is_arithmetic())
      .value("Match", RDKit::Match)
      .value("FingerprintVariance", RDKit::FingerprintVariance);

  nb::class_<RDKit::RGroupDecompositionParameters>(
      m, "RGroupDecompositionParameters",
      R"DOC(RGroupDecompositionParameters controls how the RGroupDecomposition
sets labelling and matches structures
OPTIONS:
  - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or
RGroupCoreAlignment.MCS
                         If set to MCS, cores labels are mapped to
each other using their
                         Maximum common substructure overlap.
  - RGroupLabels: optionally set where the rgroup labels to use are
encoded.
                   RGroupLabels.IsotopeLabels - labels are stored
on isotopes
                   RGroupLabels.AtomMapLabels - labels are stored
on atommaps
                   RGroupLabels.MDLRGroupLabels - labels are stored
on MDL R-groups
                   RGroupLabels.DummyAtomLabels - labels are stored
on dummy atoms
                   RGroupLabels.AtomIndexLabels - use the atom index
as the label
                   RGroupLabels.RelabelDuplicateLabels - fix any
duplicate labels
                   RGroupLabels.AutoDetect - auto detect the label
[default]
     Note: in all cases, any rgroups found on unlabelled atoms will
be automatically
            labelled.
  - RGroupLabelling: choose where the rlabels are stored on the
decomposition
                      RGroupLabelling.AtomMap - store rgroups as atom
maps (for smiles)
                      RGroupLabelling.Isotope - store rgroups on the
isotope
                      RGroupLabelling.MDLRGroup - store rgroups as mdl
rgroups (for molblocks)
                     default: AtomMap | MDLRGroup
  - onlyMatchAtRGroups: only allow rgroup decomposition at the
specified rgroups
  - removeAllHydrogenRGroups: remove all user-defined rgroups that
only have hydrogens
  - removeAllHydrogenRGroupsAndLabels: remove all user-defined
rgroups that only have hydrogens, and also remove the corresponding
labels from the core
  - removeHydrogensPostMatch: remove all hydrogens from the output
molecules
  - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or
more
  - doTautomers: match all tautomers of a core against each
input structure
  - doEnumeration: expand input cores into enumerated mol bundles
  - allowMultipleRGroupsOnUnlabelled: permit more than one rgroup to
be attached to an unlabelled core atom
  - allowMultipleCoresInSameMol: permit a core to match more than
once in the same molecule if the sets of matched atoms are not equal
(default=False))DOC")
      .def(nb::init<>(), "Constructor, takes no arguments")
      .def_rw("labels", &RDKit::RGroupDecompositionParameters::labels)
      .def_rw("matchingStrategy",
              &RDKit::RGroupDecompositionParameters::matchingStrategy)
      .def_rw("scoreMethod", &RDKit::RGroupDecompositionParameters::scoreMethod)
      .def_rw("rgroupLabelling",
              &RDKit::RGroupDecompositionParameters::rgroupLabelling)
      .def_rw("alignment", &RDKit::RGroupDecompositionParameters::alignment)
      .def_rw("chunkSize", &RDKit::RGroupDecompositionParameters::chunkSize)
      .def_rw("onlyMatchAtRGroups",
              &RDKit::RGroupDecompositionParameters::onlyMatchAtRGroups)
      .def_rw("removeAllHydrogenRGroups",
              &RDKit::RGroupDecompositionParameters::removeAllHydrogenRGroups)
      .def_rw("removeHydrogensPostMatch",
              &RDKit::RGroupDecompositionParameters::removeHydrogensPostMatch)
      .def_rw("timeout", &RDKit::RGroupDecompositionParameters::timeout)
      .def_rw("gaPopulationSize",
              &RDKit::RGroupDecompositionParameters::gaPopulationSize)
      .def_rw("gaMaximumOperations",
              &RDKit::RGroupDecompositionParameters::gaMaximumOperations)
      .def_rw("gaNumberOperationsWithoutImprovement",
              &RDKit::RGroupDecompositionParameters::
                  gaNumberOperationsWithoutImprovement)
      .def_rw("gaRandomSeed",
              &RDKit::RGroupDecompositionParameters::gaRandomSeed)
      .def_rw("gaNumberRuns",
              &RDKit::RGroupDecompositionParameters::gaNumberRuns)
      .def_rw("gaParallelRuns",
              &RDKit::RGroupDecompositionParameters::gaParallelRuns)
      .def_rw("allowNonTerminalRGroups",
              &RDKit::RGroupDecompositionParameters::allowNonTerminalRGroups)
      .def_rw("removeAllHydrogenRGroupsAndLabels",
              &RDKit::RGroupDecompositionParameters::
                  removeAllHydrogenRGroupsAndLabels)
      .def_rw("allowMultipleRGroupsOnUnlabelled",
              &RDKit::RGroupDecompositionParameters::
                  allowMultipleRGroupsOnUnlabelled)
      .def_rw(
          "allowMultipleCoresInSameMol",
          &RDKit::RGroupDecompositionParameters::allowMultipleCoresInSameMol)
      .def_rw("doTautomers", &RDKit::RGroupDecompositionParameters::doTautomers)
      .def_rw("doEnumeration",
              &RDKit::RGroupDecompositionParameters::doEnumeration)
      .def_ro("substructMatchParams",
              &RDKit::RGroupDecompositionParameters::substructmatchParams)
      .def_rw("includeTargetMolInResults",
              &RDKit::RGroupDecompositionParameters::includeTargetMolInResults)
      .def("__setattr__", &safeSetattr);

  nb::class_<RDKit::RGroupDecompositionHelper>(m, "RGroupDecomposition")
      .def(nb::init<nb::object>(), "cores"_a,
           "Construct from a molecule or sequence of molecules")
      .def(nb::init<nb::object, const RDKit::RGroupDecompositionParameters &>(),
           "cores"_a, "params"_a,
           "Construct from a molecule or sequence of molecules and a "
           "parameters object")
      .def("Add", &RDKit::RGroupDecompositionHelper::Add, "mol"_a)
      .def("GetMatchingCoreIdx",
           &RDKit::RGroupDecompositionHelper::GetMatchingCoreIdx, "mol"_a,
           "matches"_a = nb::none())
      .def("Process", &RDKit::RGroupDecompositionHelper::Process,
           "Process the rgroups (must be done prior to "
           "GetRGroupsAsRows/Columns and GetRGroupLabels)")
      .def("ProcessAndScore",
           &RDKit::RGroupDecompositionHelper::ProcessAndScore,
           "Process the rgroups and returns the score (must be done prior to "
           "GetRGroupsAsRows/Columns and GetRGroupLabels)")
      .def("GetRGroupLabels",
           &RDKit::RGroupDecompositionHelper::GetRGroupLabels,
           R"DOC(Return the current list of found rgroups.
Note, Process() should be called first)DOC")
      .def("GetRGroupsAsRows",
           &RDKit::RGroupDecompositionHelper::GetRGroupsAsRows,
           "asSmiles"_a = false,
           R"DOC(Return the rgroups as rows (note: can be fed directly into a
pandas datatable)
ARGUMENTS:
 - asSmiles: if True return smiles strings, otherwise return
molecules [default: False]
  Row structure:
     rows[idx] = {rgroup_label: molecule_or_smiles}
)DOC")
      .def("GetRGroupsAsColumns",
           &RDKit::RGroupDecompositionHelper::GetRGroupsAsColumns,
           "asSmiles"_a = false,
           R"DOC(Return the rgroups as columns (note: can be fed directly into a
pandas datatable)
ARGUMENTS:
 - asSmiles: if True return smiles strings, otherwise return
molecules [default: False]
  Column structure:
     columns[rgroup_label] = [ mols_or_smiles ]
)DOC");

  m.def("RGroupDecompose", &RDKit::RGroupDecomp, "cores"_a, "mols"_a,
        "asSmiles"_a = false, "asRows"_a = true,
        "options"_a = RDKit::RGroupDecompositionParameters(),
        R"DOC(Decompose a collection of molecules into their Rgroups
ARGUMENTS:
  - cores: a set of cores from most to least specific.
           See RGroupDecompositionParameters for more details
           on how the cores can be labelled
  - mols: the molecules to be decomposed
  - asSmiles: if True return smiles strings, otherwise return
molecules [default: False]
  - asRows: return the results as rows (default) otherwise return
columns
  - options: RGroupDecompositionParameters object that defines
the parameters for the decomposition.
           See RGroupDecompositionParameters for defaults

RETURNS: row_or_column_results, unmatched

  Row structure:
     rows[idx] = {rgroup_label: molecule_or_smiles}
  Column structure:
     columns[rgroup_label] = [ mols_or_smiles ]

  unmatched is a vector of indices in the input mols that were not
matched.
)DOC");

  m.def("RelabelMappedDummies", &RDKit::relabelMappedDummiesHelper, "mol"_a,
        "inputLabels"_a = static_cast<unsigned int>(
            RDKit::AtomMap | RDKit::Isotope | RDKit::MDLRGroup),
        "outputLabels"_a = static_cast<unsigned int>(RDKit::MDLRGroup),
        R"DOC(Relabel dummy atoms bearing an R-group mapping (as
atom map number, isotope or MDLRGroup label) such that
they will be displayed by the rendering code as R# rather
than #*, *:#, #*:#, etc. By default, only the MDLRGroup label
is retained on output; this may be configured through the
outputLabels parameter.
In case there are multiple potential R-group mappings,
the priority on input is Atom map number > Isotope > MDLRGroup.
The inputLabels parameter allows to configure which mappings
are taken into consideration.
)DOC");
}
