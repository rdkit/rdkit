//  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
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
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>
#include <math.h>

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
        if (!*iter) throw_value_error("reaction called with None reactants");
        coreMols.push_back(*iter);
        ++iter;
      }
      decomp = new RGroupDecomposition(coreMols, params);
    }
  }

  int Add(const ROMol &mol) { return decomp->add(mol); }
  bool Process() { return decomp->process(); }

  python::list GetRGroupsAsRows(bool asSmiles=false) {
    const RGroupRows &groups = decomp->getRGroupsAsRows();
    python::list result;

    for (RGroupRows::const_iterator it = groups.begin(); it != groups.end();
         ++it) {
      python::dict dict;
      const RGroupRow &side_chains = *(it);
      for (RGroupRow::const_iterator sit = side_chains.begin();
           sit != side_chains.end(); ++sit) {
        if (asSmiles) {
          dict[sit->first] = MolToSmiles(*sit->second, true);
        } else {
          dict[sit->first] = sit->second;
        }
      }
      result.append(dict);
    }
    return result;
  }

  python::dict GetRGroupsAsColumn(bool asSmiles=false) {
    python::dict result;

    RGroupColumns groups = decomp->getRGroupsAsColumns();

    for (RGroupColumns::const_iterator it = groups.begin(); it != groups.end();
         ++it) {
      python::list col;

      for (RGroupColumn::const_iterator cit = it->second.begin();
           cit != it->second.end(); ++cit) {
        if (asSmiles) {
          col.append(MolToSmiles(**cit,true));
        } else {
          col.append(*cit);
        }
      }
      result[it->first] = col;
    }
    return result;
  }

  
};

python::object RGroupDecomp(python::object cores,
                            python::object mols,
                            bool asSmiles = false,
                            bool asRows = true,
                            const RGroupDecompositionParameters &options = RGroupDecompositionParameters()
                            ) {
  RGroupDecompositionHelper decomp(cores, options);
  python::list unmatched;

  python::stl_input_iterator<ROMOL_SPTR> iter(mols), end;
  unsigned int idx=0;
  while (iter != end) {
    if (!*iter) throw_value_error("reaction called with None reactants");
    if(decomp.Add(*(*iter)) == -1) {
      unmatched.append(idx);
    }
    ++iter;
    ++idx;
  }

  decomp.Process();
  if ( asRows ) {
    return make_tuple(decomp.GetRGroupsAsRows(asSmiles), unmatched);
  } else {
    return make_tuple(decomp.GetRGroupsAsColumn(asSmiles), unmatched);
  }
}


struct rgroupdecomp_wrapper {
  static void wrap() {
    // logic from https://stackoverflow.com/a/13017303
    boost::python::type_info info =
        boost::python::type_id<RDKit::MOL_SPTR_VECT>();
    const boost::python::converter::registration *reg =
        boost::python::converter::registry::query(info);
    if (reg == NULL || (*reg).m_to_python == NULL) {
      python::class_<RDKit::MOL_SPTR_VECT>("MOL_SPTR_VECT")
          .def(python::vector_indexing_suite<RDKit::MOL_SPTR_VECT, true>());
    }
    std::string docString = "";
    python::enum_<RDKit::RGroupLabels>("RGroupLabels")
        .value("IsotopeLabels", RDKit::IsotopeLabels)
        .value("AtomMapLabels", RDKit::AtomMapLabels)
        .value("AtomIndexLabels", RDKit::AtomIndexLabels)
        .value("RelabelDuplicateLabels", RDKit::RelabelDuplicateLabels)
        .value("AutoDetect", RDKit::AutoDetect)
        .export_values();

    python::enum_<RDKit::RGroupMatching>("RGroupMatching")
        .value("Greedy", RDKit::Greedy)
        .value("GreedyChunks", RDKit::GreedyChunks)
        .value("Exhaustive", RDKit::Exhaustive)
        .export_values();

    python::enum_<RDKit::RGroupLabelling>("RGroupLabelling")
        .value("AtomMap", RDKit::AtomMap)
        .value("Isotope", RDKit::Isotope)
        .value("MDLRGroup", RDKit::MDLRGroup)
        .export_values();

    python::enum_<RDKit::RGroupCoreAlignment>("RGroupCoreAlignment")
        .value("None", RDKit::None)
        .value("MCS", RDKit::MCS)
        .export_values();

    docString = "RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures\n"
        "  OPTIONS:\n"
        "    - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None or RGroupCoreAlignment.MCS\n"
        "                           If set to MCS, cores labels are mapped to each other using their\n"
        "                           Maximum common substructure overlap.\n"
        "    - RGroupLabels: optionally set where the rgroup labels to use are encoded.\n"
        "                     RroupLabels.IsotopeLabels - labels are stored on isotopes\n"
        "                     RroupLabels.AtomMapLabels - labels are stored on atommaps\n"
        "                     RroupLabels.AtomIndexLabels - use the atom index as the label\n"
        "                     RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels\n"
        "                     RGroupLabels.AutoDetect - auto detect the label [default]\n"
        "       Note: in all cases, any rgroups found on unlabelled atoms will be automatically\n"
        "              labelled.\n"
        "    - RGroupLabelling: choose where the rlabels are stored on the decomposition\n"
        "                        RGroupLabels.AtomMap - store rgroups as atom maps (for smiles)\n"
        "                        RGroupLabels.Isotope - stroe rgroups on the isotope\n"
        "                        RGroupLabels.MDLRGroup - store rgroups as mdl rgroups (for molblocks)\n"
        "                       default: AtomMap | MDLRGroup\n"
        "    - matchOnlyAtRGroups: only allow rgroup decomposition at the specified rgroups\n"
        "    - setRemoveRGroupsThatAreAllHydrogen: remove all rgroups that only have hydrogens\n"
        "    - removeHydrogensPostMatch: remove all hydrogens from the output molecules\n";
    python::class_<RDKit::RGroupDecompositionParameters>(
        "RGroupDecompositionParameters", docString.c_str(),
        python::init<>("Constructor, takes no arguments"))

        .def(python::init<RGroupLabels, RGroupMatching, RGroupLabelling,
                          RGroupCoreAlignment, unsigned int, bool, bool>())
        .def_readwrite("labels", &RDKit::RGroupDecompositionParameters::labels)
        .def_readwrite("matchingStrategy",
                       &RDKit::RGroupDecompositionParameters::matchingStrategy)
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
            &RDKit::RGroupDecompositionParameters::removeHydrogensPostMatch);

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
             "GetRGroupsAsRows/Columns)")
        .def("GetRGroupsAsRows", &RGroupDecompositionHelper::GetRGroupsAsRows,
             python::arg("asSmiles")=false,
             "Return the rgroups as rows (note: can be fed directrly into a pandas datatable)\n"
             "  ARGUMENTS:\n"
             "   - asSmiles: if True return smiles strings, otherwise return molecules [default: False]\n"
             "    Row structure:\n"
             "       rows[idx] = {rgroup_label: molecule_or_smiles}\n"
             )
        .def("GetRGroupsAsColumns",
             &RGroupDecompositionHelper::GetRGroupsAsColumn,
             python::arg("asSmiles")=false,
             "Return the rgroups as columns (note: can be fed directrly into a pandas datatable)\n"
             "  ARGUMENTS:\n"
             "   - asSmiles: if True return smiles strings, otherwise return molecules [default: False]\n"
             "    Column structure:\n"
             "       columns[rgroup_label] = [ mols_or_smiles ]\n"             
             );
    
    docString = "Decompose a collecion of molecules into their Rgroups\n"
        "  ARGUMENTS:\n"
        "    - cores: a set of cores from most to least specific.\n"
        "             See RGroupDecompositionParameters for more details\n"
        "             on how the cores can be labelled\n"
        "    - mols: the molecules to be decomposed\n"
        "    - asSmiles: if True return smiles strings, otherwise return molecules [default: False]\n"
        "    - asRows: return the results as rows (default) otherwise return columns\n"
        "\n"
        "  RETURNS: row_or_column_results, unmatched\n"
        "\n"
        "    Row structure:\n"
        "       rows[idx] = {rgroup_label: molecule_or_smiles}\n"
        "    Column structure:\n"
        "       columns[rgroup_label] = [ mols_or_smiles ]\n"
        "\n"
        "    unmatched is a vector of indices in the input mols that were not matched.\n";
    python::def(
        "RGroupDecompose", RDKit::RGroupDecomp,
        (python::arg("cores"), python::arg("mols"),
         python::arg("asSmiles") = false,
         python::arg("asRows") = true,
         python::arg("options") = RGroupDecompositionParameters()),
        docString.c_str());
  };  
};
}

BOOST_PYTHON_MODULE(rdRGroupDecomposition) {
  python::scope().attr("__doc__") =
      "Module containing RGroupDecomposition classes and functions.";
  RDKit::rgroupdecomp_wrapper::wrap();
}
