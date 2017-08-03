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
#include <DataStructs/ExplicitBitVect.h>
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

  python::list GetRGroupsAsRows() {
    const RGroupRows &groups = decomp->getRGroupsAsRows();
    python::list result;

    for (RGroupRows::const_iterator it = groups.begin(); it != groups.end();
         ++it) {
      python::dict dict;
      const RGroupRow &side_chains = *(it);
      for (RGroupRow::const_iterator sit = side_chains.begin();
           sit != side_chains.end(); ++sit) {
        dict[sit->first] = sit->second;
      }
      result.append(dict);
    }
    return result;
  }

  python::dict GetRGroupsAsColumn() {
    python::dict result;

    RGroupColumns groups = decomp->getRGroupsAsColumns();

    for (RGroupColumns::const_iterator it = groups.begin(); it != groups.end();
         ++it) {
      python::list col;

      for (RGroupColumn::const_iterator cit = it->second.begin();
           cit != it->second.end(); ++cit) {
        col.append(*cit);
      }
      result[it->first] = col;
    }
    return result;
  }
};

struct rgroupdecomp_wrapper {
  static void wrap() {
    python::class_<RDKit::MOL_SPTR_VECT>("MOL_SPTR_VECT")
        .def(python::vector_indexing_suite<RDKit::MOL_SPTR_VECT, true>());

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

    python::class_<RDKit::RGroupDecompositionParameters>(
        "RGroupDecompositionParameters", docString.c_str(),
        python::init<>("Constructor, takes no arguments"))

        .def(python::init<RGroupLabels, RGroupMatching, RGroupLabelling,
                          RGroupCoreAlignment, unsigned int, bool, bool>())
        .def("SetRGroupLabels",
             &RDKit::RGroupDecompositionParameters::SetRGroupLabels)
        .def("GetRGroupLabels",
             &RDKit::RGroupDecompositionParameters::GetRGroupLabels)
        .def("SetRGroupLabelling",
             &RDKit::RGroupDecompositionParameters::SetRGroupLabelling)
        .def("GetRGroupLabelling",
             &RDKit::RGroupDecompositionParameters::GetRGroupLabelling)
        .def("SetRGroupMatching",
             &RDKit::RGroupDecompositionParameters::SetRGroupMatching)
        .def("GetRGroupMatching",
             &RDKit::RGroupDecompositionParameters::GetRGroupMatching)
        .def("SetRGroupCoreAlignment",
             &RDKit::RGroupDecompositionParameters::SetRGroupCoreAlignment)
        .def("GetRGroupCoreAlignment",
             &RDKit::RGroupDecompositionParameters::GetRGroupCoreAlignment)

        ;

    python::class_<RDKit::RGroupDecompositionHelper, boost::noncopyable>(
        "RGroupDecomposition", docString.c_str(),
        python::init<python::object>(
            "Construct from a molecule or sequence of molecules"))
        .def("Add", &RGroupDecompositionHelper::Add)
        .def("Process", &RGroupDecompositionHelper::Process,
             "Process the rgroups (must be done prior to "
             "GetRGroupsAsRows/Columns)")
        .def("GetRGroupsAsRows", &RGroupDecompositionHelper::GetRGroupsAsRows)
        .def("GetRGroupsAsColumns",
             &RGroupDecompositionHelper::GetRGroupsAsColumn);
  };
};
}

BOOST_PYTHON_MODULE(rdRGroupDecomposition) {
  python::scope().attr("__doc__") =
      "Module containing RGroupDecomposition classes and functions.";
  RDKit::rgroupdecomp_wrapper::wrap();
}
