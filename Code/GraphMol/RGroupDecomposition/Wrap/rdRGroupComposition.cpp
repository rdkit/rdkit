//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
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
#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <boost/python/list.hpp>
#include <string>
#include <cmath>
#include <chrono>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/python_streambuf.h>

namespace python = boost::python;
using boost_adaptbx::python::streambuf;

namespace RDKit {

class RGroupDecompositionHelper {
  RGroupDecomposition *decomp;

 public:
  ~RGroupDecompositionHelper() { delete decomp; }

  RGroupDecompositionHelper(python::object cores,
                            const RGroupDecompositionParameters &params =
                                RGroupDecompositionParameters()) {
    python::extract<ROMol> isROMol(cores);
    if (isROMol.check()) {
      decomp = new RGroupDecomposition(isROMol(), params);
    } else {
      MOL_SPTR_VECT coreMols;
      python::stl_input_iterator<ROMOL_SPTR> iter(cores), end;
      while (iter != end) {
        if (!*iter) {
          throw_value_error("reaction called with None reactants");
        }
        coreMols.push_back(*iter);
        ++iter;
      }
      decomp = new RGroupDecomposition(coreMols, params);
    }
  }

  int Add(const ROMol &mol) {
    NOGIL gil;
    return decomp->add(mol);
  }
  bool Process() {
    NOGIL gil;
    return decomp->process();
  }
  python::tuple ProcessAndScore() {
    NOGIL gil;
    auto result = decomp->processAndScore();
    return python::make_tuple(result.success, result.score);
  }

  python::list GetRGroupLabels() {
    python::list result;
    std::vector<std::string> labels = decomp->getRGroupLabels();
    for (auto label : labels) {
      result.append(label);
    }
    return result;
  }
  python::list GetRGroupsAsRows(bool asSmiles = false) {
    const RGroupRows &groups = decomp->getRGroupsAsRows();
    python::list result;

    for (const auto &side_chains : groups) {
      python::dict dict;
      for (const auto &side_chain : side_chains) {
        if (asSmiles) {
          dict[side_chain.first] = MolToSmiles(*side_chain.second, true);
        } else {
          dict[side_chain.first] = side_chain.second;
        }
      }
      result.append(dict);
    }
    return result;
  }

  python::dict GetRGroupsAsColumn(bool asSmiles = false) {
    python::dict result;

    RGroupColumns groups = decomp->getRGroupsAsColumns();

    for (RGroupColumns::const_iterator it = groups.begin(); it != groups.end();
         ++it) {
      python::list col;

      for (const auto &cit : it->second) {
        if (asSmiles) {
          col.append(MolToSmiles(*cit, true));
        } else {
          col.append(cit);
        }
      }
      result[it->first] = col;
    }
    return result;
  }
};

python::object RGroupDecomp(python::object cores, python::object mols,
                            bool asSmiles = false, bool asRows = true,
                            const RGroupDecompositionParameters &options =
                                RGroupDecompositionParameters()) {
  auto t0 = std::chrono::steady_clock::now();
  RGroupDecompositionHelper decomp(cores, options);
  python::list unmatched;

  python::stl_input_iterator<ROMOL_SPTR> iter(mols), end;
  unsigned int idx = 0;
  while (iter != end) {
    if (!*iter) {
      throw_value_error("reaction called with None reactants");
    }
    if (decomp.Add(*(*iter)) == -1) {
      unmatched.append(idx);
    }
    ++iter;
    ++idx;
    checkForTimeout(t0, options.timeout);
  }

  decomp.Process();
  if (asRows) {
    return make_tuple(decomp.GetRGroupsAsRows(asSmiles), unmatched);
  } else {
    return make_tuple(decomp.GetRGroupsAsColumn(asSmiles), unmatched);
  }
}  // namespace RDKit

struct rgroupdecomp_wrapper {
  static void wrap() {
    bool noproxy = true;
    RegisterVectorConverter<RDKit::ROMOL_SPTR>("MOL_SPTR_VECT", noproxy);

    std::string docString = "";
    python::enum_<RDKit::RGroupLabels>("RGroupLabels")
        .value("IsotopeLabels", RDKit::IsotopeLabels)
        .value("AtomMapLabels", RDKit::AtomMapLabels)
        .value("AtomIndexLabels", RDKit::AtomIndexLabels)
        .value("RelabelDuplicateLabels", RDKit::RelabelDuplicateLabels)
        .value("MDLRGroupLabels", RDKit::MDLRGroupLabels)
        .value("DummyAtomLabels", RDKit::DummyAtomLabels)
        .value("AutoDetect", RDKit::AutoDetect)
        .export_values();

    python::enum_<RDKit::RGroupMatching>("RGroupMatching")
        .value("Greedy", RDKit::Greedy)
        .value("GreedyChunks", RDKit::GreedyChunks)
        .value("Exhaustive", RDKit::Exhaustive)
        .value("NoSymmetrization", RDKit::NoSymmetrization)
        .value("GA", RDKit::GA)
        .export_values();

    python::enum_<RDKit::RGroupLabelling>("RGroupLabelling")
        .value("AtomMap", RDKit::AtomMap)
        .value("Isotope", RDKit::Isotope)
        .value("MDLRGroup", RDKit::MDLRGroup)
        .export_values();

    python::enum_<RDKit::RGroupCoreAlignment>("RGroupCoreAlignment")
        // DEPRECATED, remove the folowing line in release 2021.03
        .value("None", RDKit::NoAlignment)
        .value("NoAlignment", RDKit::NoAlignment)
        .value("MCS", RDKit::MCS)
        .export_values();

    python::enum_<RDKit::RGroupScore>("RGroupScore")
        .value("Match", RDKit::Match)
        .value("FingerprintVariance", RDKit::FingerprintVariance)
        .export_values();

    docString =
        "RGroupDecompositionParameters controls how the RGroupDecomposition "
        "sets labelling and matches structures\n"
        "  OPTIONS:\n"
        "    - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or "
        "RGroupCoreAlignment.MCS\n"
        "                           If set to MCS, cores labels are mapped to "
        "each other using their\n"
        "                           Maximum common substructure overlap.\n"
        "    - RGroupLabels: optionally set where the rgroup labels to use are "
        "encoded.\n"
        "                     RGroupLabels.IsotopeLabels - labels are stored "
        "on isotopes\n"
        "                     RGroupLabels.AtomMapLabels - labels are stored "
        "on atommaps\n"
        "                     RGroupLabels.MDLRGroupLabels - labels are stored "
        "on MDL R-groups\n"
        "                     RGroupLabels.DummyAtomLabels - labels are stored "
        "on dummy atoms\n"
        "                     RGroupLabels.AtomIndexLabels - use the atom "
        "index "
        "as the label\n"
        "                     RGroupLabels.RelabelDuplicateLabels - fix any "
        "duplicate labels\n"
        "                     RGroupLabels.AutoDetect - auto detect the label "
        "[default]\n"
        "       Note: in all cases, any rgroups found on unlabelled atoms will "
        "be automatically\n"
        "              labelled.\n"
        "    - RGroupLabelling: choose where the rlabels are stored on the "
        "decomposition\n"
        "                        RGroupLabels.AtomMap - store rgroups as atom "
        "maps (for smiles)\n"
        "                        RGroupLabels.Isotope - store rgroups on the "
        "isotope\n"
        "                        RGroupLabels.MDLRGroup - store rgroups as mdl "
        "rgroups (for molblocks)\n"
        "                       default: AtomMap | MDLRGroup\n"
        "    - onlyMatchAtRGroups: only allow rgroup decomposition at the "
        "specified rgroups\n"
        "    - removeAllHydrogenRGroups: remove all user-defined rgroups that "
        "only have hydrogens\n"
        "    - removeAllHydrogenRGroupsAndLabels: remove all user-defined "
        "rgroups that only have hydrogens, and also remove the corresponding "
        "labels from the core\n"
        "    - removeHydrogensPostMatch: remove all hydrogens from the output "
        "molecules\n"
        "    - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or "
        "more\n";
    python::class_<RDKit::RGroupDecompositionParameters>(
        "RGroupDecompositionParameters", docString.c_str(),
        python::init<>("Constructor, takes no arguments"))

        .def_readwrite("labels", &RDKit::RGroupDecompositionParameters::labels)
        .def_readwrite("matchingStrategy",
                       &RDKit::RGroupDecompositionParameters::matchingStrategy)
        .def_readwrite("scoreMethod",
                       &RDKit::RGroupDecompositionParameters::scoreMethod)
        .def_readwrite("rgroupLabelling",
                       &RDKit::RGroupDecompositionParameters::rgroupLabelling)
        .def_readwrite("alignment",
                       &RDKit::RGroupDecompositionParameters::alignment)
        .def_readwrite("chunkSize",
                       &RDKit::RGroupDecompositionParameters::chunkSize)
        .def_readwrite(
            "onlyMatchAtRGroups",
            &RDKit::RGroupDecompositionParameters::onlyMatchAtRGroups)
        .def_readwrite(
            "removeAllHydrogenRGroups",
            &RDKit::RGroupDecompositionParameters::removeAllHydrogenRGroups)
        .def_readwrite(
            "removeHydrogensPostMatch",
            &RDKit::RGroupDecompositionParameters::removeHydrogensPostMatch)
        .def_readwrite("timeout",
                       &RDKit::RGroupDecompositionParameters::timeout)
        .def_readwrite("gaPopulationSize",
                       &RDKit::RGroupDecompositionParameters::gaPopulationSize)
        .def_readwrite(
            "gaMaximumOperations",
            &RDKit::RGroupDecompositionParameters::gaMaximumOperations)
        .def_readwrite("gaNumberOperationsWithoutImprovement",
                       &RDKit::RGroupDecompositionParameters::
                           gaNumberOperationsWithoutImprovement)
        .def_readwrite("gaRandomSeed",
                       &RDKit::RGroupDecompositionParameters::gaRandomSeed)
        .def_readwrite("gaNumberRuns",
                       &RDKit::RGroupDecompositionParameters::gaNumberRuns)
        .def_readwrite("gaParallelRuns",
                       &RDKit::RGroupDecompositionParameters::gaParallelRuns)
        .def_readwrite(
            "allowNonTerminalRGroups",
            &RDKit::RGroupDecompositionParameters::allowNonTerminalRGroups)
        .def_readwrite("removeAllHydrogenRGroupsAndLabels",
                       &RDKit::RGroupDecompositionParameters::
                           removeAllHydrogenRGroupsAndLabels)
        .def_readwrite("allowMultipleRGroupsOnUnlabelled",
                       &RDKit::RGroupDecompositionParameters::
                           allowMultipleRGroupsOnUnlabelled)
        .def_readonly(
            "substructMatchParams",
            &RDKit::RGroupDecompositionParameters::substructmatchParams);

    python::class_<RDKit::RGroupDecompositionHelper, boost::noncopyable>(
        "RGroupDecomposition", docString.c_str(),
        python::init<python::object>(
            "Construct from a molecule or sequence of molecules"))
        .def(
            python::init<python::object, const RGroupDecompositionParameters &>(
                "Construct from a molecule or sequence of molecules and a "
                "parameters object"))
        .def("Add", &RGroupDecompositionHelper::Add)
        .def("Process", &RGroupDecompositionHelper::Process,
             "Process the rgroups (must be done prior to "
             "GetRGroupsAsRows/Columns and GetRGroupLabels)")
        .def("ProcessAndScore", &RGroupDecompositionHelper::ProcessAndScore,
             "Process the rgroups and returns the score (must be done prior to "
             "GetRGroupsAsRows/Columns and GetRGroupLabels)")
        .def("GetRGroupLabels", &RGroupDecompositionHelper::GetRGroupLabels,
             "Return the current list of found rgroups.\n"
             "Note, Process() should be called first")
        .def("GetRGroupsAsRows", &RGroupDecompositionHelper::GetRGroupsAsRows,
             python::arg("asSmiles") = false,
             "Return the rgroups as rows (note: can be fed directrly into a "
             "pandas datatable)\n"
             "  ARGUMENTS:\n"
             "   - asSmiles: if True return smiles strings, otherwise return "
             "molecules [default: False]\n"
             "    Row structure:\n"
             "       rows[idx] = {rgroup_label: molecule_or_smiles}\n")
        .def("GetRGroupsAsColumns",
             &RGroupDecompositionHelper::GetRGroupsAsColumn,
             python::arg("asSmiles") = false,
             "Return the rgroups as columns (note: can be fed directrly into a "
             "pandas datatable)\n"
             "  ARGUMENTS:\n"
             "   - asSmiles: if True return smiles strings, otherwise return "
             "molecules [default: False]\n"
             "    Column structure:\n"
             "       columns[rgroup_label] = [ mols_or_smiles ]\n");

    docString =
        "Decompose a collecion of molecules into their Rgroups\n"
        "  ARGUMENTS:\n"
        "    - cores: a set of cores from most to least specific.\n"
        "             See RGroupDecompositionParameters for more details\n"
        "             on how the cores can be labelled\n"
        "    - mols: the molecules to be decomposed\n"
        "    - asSmiles: if True return smiles strings, otherwise return "
        "molecules [default: False]\n"
        "    - asRows: return the results as rows (default) otherwise return "
        "columns\n"
        "\n"
        "  RETURNS: row_or_column_results, unmatched\n"
        "\n"
        "    Row structure:\n"
        "       rows[idx] = {rgroup_label: molecule_or_smiles}\n"
        "    Column structure:\n"
        "       columns[rgroup_label] = [ mols_or_smiles ]\n"
        "\n"
        "    unmatched is a vector of indices in the input mols that were not "
        "matched.\n";
    python::def("RGroupDecompose", RDKit::RGroupDecomp,
                (python::arg("cores"), python::arg("mols"),
                 python::arg("asSmiles") = false, python::arg("asRows") = true,
                 python::arg("options") = RGroupDecompositionParameters()),
                docString.c_str());
  };
};
}  // namespace RDKit

BOOST_PYTHON_MODULE(rdRGroupDecomposition) {
  python::scope().attr("__doc__") =
      "Module containing RGroupDecomposition classes and functions.";
  RDKit::rgroupdecomp_wrapper::wrap();
}
