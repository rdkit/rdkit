//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmmpa_array_API
#include <boost/python.hpp>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/MMPA/MMPA.h>
#include <boost/foreach.hpp>

namespace python = boost::python;

namespace {
python::tuple fragmentMolHelper(const RDKit::ROMol& mol, unsigned int maxCuts,
                                unsigned int maxCutBonds,
                                const std::string& pattern,
                                bool resultsAsMols) {
  std::vector<std::pair<RDKit::ROMOL_SPTR, RDKit::ROMOL_SPTR> > tres;
  bool ok = RDKit::MMPA::fragmentMol(mol, tres, maxCuts, maxCutBonds, pattern);
  python::list pyres;
  if (ok) {
    for (std::vector<std::pair<RDKit::ROMOL_SPTR,
                               RDKit::ROMOL_SPTR> >::const_iterator pr =
             tres.begin();
         pr != tres.end(); ++pr) {
      python::list lres;
      if (resultsAsMols) {
        lres.append(pr->first);
        lres.append(pr->second);
      } else {
        if (pr->first) {
          lres.append(RDKit::MolToSmiles(*(pr->first), true));
        } else {
          lres.append("");
        }
        lres.append(RDKit::MolToSmiles(*(pr->second), true));
      }
      pyres.append(python::tuple(lres));
    }
  }
  return python::tuple(pyres);
}
}

BOOST_PYTHON_MODULE(rdMMPA) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of code for doing MMPA";

  std::string docString =
      "Does the fragmentation necessary for an MMPA analysis";
  python::def("FragmentMol", fragmentMolHelper,
              (python::arg("mol"), python::arg("maxCuts") = 3,
               python::arg("maxCutBonds") = 20,
               python::arg("pattern") = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]",
               python::arg("resultsAsMols") = true),
              docString.c_str());
}
