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
#include <cmath>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <RDFreeSASA.h>
#include <RDBoost/Wrap.h>

namespace python = boost::python;

namespace RDKit {
namespace {
python::object classifyAtomsHelper(RDKit::ROMol &mol,
                                   const FreeSASA::SASAOpts &opts) {
  std::vector<double> radii;
  python::list l;
  if (FreeSASA::classifyAtoms(mol, radii, opts)) {
    for (double &i : radii) {
      l.append(i);
    }
    return std::move(l);
  }
  return std::move(l);
}

double calcSASAHelper(const RDKit::ROMol &mol, python::object radii,
                      int confIdx, const RDKit::Atom *query,
                      const FreeSASA::SASAOpts &opts) {
  const RDKit::QueryAtom *atom = nullptr;
  if (query) {
    atom = dynamic_cast<const RDKit::QueryAtom *>(query);
    if (!atom) {
      throw ValueErrorException("Query is not a query atom!");
    }
  }

  std::vector<double> vradii;

  unsigned int sz = python::extract<unsigned int>(radii.attr("__len__")());
  for (unsigned int i = 0; i < sz; ++i) {
    vradii.push_back(python::extract<double>(radii[i])());
  }

  return FreeSASA::calcSASA(mol, vradii, confIdx, atom, opts);
}
}  // namespace

struct freesasa_wrapper {
  static void wrap() {
    std::string docString = "";
    python::enum_<FreeSASA::SASAOpts::Algorithm>("SASAAlgorithm")
        .value("LeeRichards", FreeSASA::SASAOpts::LeeRichards)
        .value("ShrakeRupley", FreeSASA::SASAOpts::ShrakeRupley)
        .export_values();

    python::enum_<FreeSASA::SASAOpts::Classifier>("SASAClassifier")
        .value("Protor", FreeSASA::SASAOpts::Protor)
        .value("NACCESS", FreeSASA::SASAOpts::NACCESS)
        .value("OONS", FreeSASA::SASAOpts::OONS)
        .export_values();

    python::enum_<FreeSASA::SASAOpts::Classes>("SASAClass")
        .value("Unclassified", FreeSASA::SASAOpts::Unclassified)
        .value("APolar", FreeSASA::SASAOpts::APolar)
        .value("Polar", FreeSASA::SASAOpts::Polar)
        .export_values();

    python::class_<FreeSASA::SASAOpts>(
        "SASAOpts", docString.c_str(),
        python::init<>("Constructor takes no arguments"))
        .def(python::init<FreeSASA::SASAOpts::Algorithm,
                          FreeSASA::SASAOpts::Classifier>())
        .def(python::init<FreeSASA::SASAOpts::Algorithm,
                          FreeSASA::SASAOpts::Classifier, double>())
        .def_readwrite("algorithm", &FreeSASA::SASAOpts::algorithm)
        .def_readwrite("classifier", &FreeSASA::SASAOpts::classifier)
        .def_readwrite("probeRadius", &FreeSASA::SASAOpts::probeRadius);

    docString =
        "Classify the atoms in the molecule returning their radii if "
        "possible.\n"
        "ARGUMENTS:\n"
        "   - mol: molecule to classify\n"
        "   - options: FreeSASA options class specifying the classification "
        "method.\n"
        "               Current classifiers are Protor, NACCESS and OONS\n"
        "               classification is stored as atom property 'SASAClass' "
        "for the integer value\n"
        "                and 'SASAClassName' for the string name of the class, "
        "Polar, APolar...\n"
        "\n"
        "RETURNS:\n"
        "  list of radii where radii[atom.GetIdx()] is the radii of the atom.\n"
        "  If classification fails, NONE is returned\n";

    python::def(
        "classifyAtoms", classifyAtomsHelper,
        (python::arg("mol"), python::arg("options") = FreeSASA::SASAOpts()),
        docString.c_str());

    docString =
        "Compute the Solvent Accessible Surface Area using the FreeSASA "
        "library\n"
        "ARGUMENTS:\n"
        "  - mol: The molecule to compute.\n"
        "  - radii:  A list of atom raddii where radii[atom.GetIdx()] is the "
        "radius of the atom\n"
        "            These can be passed in or calculated with classifyAtoms "
        "for some proteins\n"
        "  - confIdx: Specify the conformer to use for the 3D geometry  "
        "[default -1]\n"
        "  - query: Pass along a query atom to compute the SASA for a subset "
        "of atoms.\n"
        "           precanned query atoms can be made with "
        "MakeFreeSasaPolarAtomQuery and\n"
        "           MakeFreeSasaAPolarAtomQuery for classified polar and "
        "apolar atoms respectively.\n"
        "  - opts: a SASAOpts class specifying the algorithm to use\n"
        "\n"
        "RETURNS:\n"
        "The computed solvent accessible surface area.\n";

    python::def(
        "CalcSASA", calcSASAHelper,
        (python::arg("mol"), python::arg("radii"), python::arg("confIdx") = -1,
         python::arg("query") = python::object(),
         python::arg("opts") = FreeSASA::SASAOpts()),
        docString.c_str());

    python::def(
        "MakeFreeSasaAPolarAtomQuery", FreeSASA::makeFreeSasaAPolarAtomQuery,
        python::return_value_policy<python::manage_new_object>(),
        "Returns an APolar atom query for use with CalcSASA.  An apolar atom "
        "has the SASAClass\n"
        "and SASAClassName set to the APOLAR class.  (see classifyAtoms)");

    python::def(
        "MakeFreeSasaPolarAtomQuery", FreeSASA::makeFreeSasaPolarAtomQuery,
        python::return_value_policy<python::manage_new_object>(),
        "Returns a polar atom query for use with CalcSASA.  An polar atom has "
        "the SASAClass\n"
        "and SASAClassName set to the POLAR class.  (see classifyAtoms)");
  }
};
}  // namespace RDKit

BOOST_PYTHON_MODULE(rdFreeSASA) {
  python::scope().attr("__doc__") =
      "Module containing rdFreeSASA classes and functions.";
  RDKit::freesasa_wrapper::wrap();
}
