//
//  Copyright (C) 2017-2026 Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>

#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <RDFreeSASA.h>
#include <RDBoost/Wrap_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(rdFreeSASA, m) {
  m.doc() = "Module containing rdFreeSASA classes and functions.";

  nb::enum_<FreeSASA::SASAOpts::Algorithm>(m, "SASAAlgorithm")
      .value("LeeRichards", FreeSASA::SASAOpts::LeeRichards)
      .value("ShrakeRupley", FreeSASA::SASAOpts::ShrakeRupley)
      .export_values();

  nb::enum_<FreeSASA::SASAOpts::Classifier>(m, "SASAClassifier")
      .value("Protor", FreeSASA::SASAOpts::Protor)
      .value("NACCESS", FreeSASA::SASAOpts::NACCESS)
      .value("OONS", FreeSASA::SASAOpts::OONS)
      .export_values();

  nb::enum_<FreeSASA::SASAOpts::Classes>(m, "SASAClass")
      .value("Unclassified", FreeSASA::SASAOpts::Unclassified)
      .value("APolar", FreeSASA::SASAOpts::APolar)
      .value("Polar", FreeSASA::SASAOpts::Polar)
      .export_values();

  nb::class_<FreeSASA::SASAOpts>(m, "SASAOpts")
      .def(nb::init<>(), "Constructor takes no arguments")
      .def(nb::init<FreeSASA::SASAOpts::Algorithm,
                    FreeSASA::SASAOpts::Classifier>(),
           "alg"_a, "cls"_a)
      .def(nb::init<FreeSASA::SASAOpts::Algorithm,
                    FreeSASA::SASAOpts::Classifier, double>(),
           "alg"_a, "cls"_a, "pr"_a)
      .def_rw("algorithm", &FreeSASA::SASAOpts::algorithm)
      .def_rw("classifier", &FreeSASA::SASAOpts::classifier)
      .def_rw("probeRadius", &FreeSASA::SASAOpts::probeRadius)
      .def("__setattr__", &safeSetattr);

  m.def(
      "classifyAtoms",
      [](RDKit::ROMol &mol, const FreeSASA::SASAOpts &opts) {
        std::vector<double> radii;
        nb::list l;
        if (FreeSASA::classifyAtoms(mol, radii, opts)) {
          for (double r : radii) {
            l.append(r);
          }
        }
        return l;
      },
      "mol"_a, "options"_a = FreeSASA::SASAOpts(),
      R"DOC(Classify the atoms in the molecule returning their radii if possible.
ARGUMENTS:
   - mol: molecule to classify
   - options: FreeSASA options class specifying the classification method.
               Current classifiers are Protor, NACCESS and OONS
               classification is stored as atom property 'SASAClass' for the integer value
                and 'SASAClassName' for the string name of the class, Polar, APolar...

RETURNS:
  list of radii where radii[atom.GetIdx()] is the radii of the atom.
  If classification fails, NONE is returned)DOC");

  m.def(
      "CalcSASA",
      [](const RDKit::ROMol &mol, nb::object radii, int confIdx,
         const RDKit::Atom *query, const FreeSASA::SASAOpts &opts) {
        const RDKit::QueryAtom *atom = nullptr;
        if (query) {
          atom = dynamic_cast<const RDKit::QueryAtom *>(query);
          if (!atom) {
            throw ValueErrorException("Query is not a query atom!");
          }
        }
        std::vector<double> vradii;
        unsigned int sz = nb::len(radii);
        for (unsigned int i = 0; i < sz; ++i) {
          vradii.push_back(nb::cast<double>(radii[i]));
        }
        return FreeSASA::calcSASA(mol, vradii, confIdx, atom, opts);
      },
      "mol"_a, "radii"_a, "confIdx"_a = -1,
      "query"_a = nb::none(), "opts"_a = FreeSASA::SASAOpts(),
      R"DOC(Compute the Solvent Accessible Surface Area using the FreeSASA library
ARGUMENTS:
  - mol: The molecule to compute.
  - radii:  A list of atom raddii where radii[atom.GetIdx()] is the radius of the atom
            These can be passed in or calculated with classifyAtoms for some proteins
  - confIdx: Specify the conformer to use for the 3D geometry  [default -1]
  - query: Pass along a query atom to compute the SASA for a subset of atoms.
           precanned query atoms can be made with MakeFreeSasaPolarAtomQuery and
           MakeFreeSasaAPolarAtomQuery for classified polar and apolar atoms respectively.
  - opts: a SASAOpts class specifying the algorithm to use

RETURNS:
The computed solvent accessible surface area.)DOC");

  m.def(
      "MakeFreeSasaAPolarAtomQuery", FreeSASA::makeFreeSasaAPolarAtomQuery,
      nb::rv_policy::take_ownership,
      R"DOC(Returns an APolar atom query for use with CalcSASA.  An apolar atom has the SASAClass
and SASAClassName set to the APOLAR class.  (see classifyAtoms))DOC");

  m.def(
      "MakeFreeSasaPolarAtomQuery", FreeSASA::makeFreeSasaPolarAtomQuery,
      nb::rv_policy::take_ownership,
      R"DOC(Returns a polar atom query for use with CalcSASA.  An polar atom has the SASAClass
and SASAClassName set to the POLAR class.  (see classifyAtoms))DOC");
}
