//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolHash/MolHash.h>
#include <RDBoost/Wrap.h>

namespace python = boost::python;
using namespace RDKit;

namespace {
std::string GenMolHashString(const ROMol &mol, python::object atomsToUse,
                             python::object bondsToUse) {
  std::unique_ptr<std::vector<unsigned> > avect;
  if (atomsToUse) {
    avect = pythonObjectToVect(atomsToUse,
                               static_cast<unsigned>(mol.getNumAtoms()));
  }
  std::unique_ptr<std::vector<unsigned> > bvect;
  if (bondsToUse) {
    bvect = pythonObjectToVect(bondsToUse,
                               static_cast<unsigned>(mol.getNumBonds()));
  }
  std::string res =
      MolHash::generateMoleculeHashSet(mol, avect.get(), bvect.get());
  return res;
}
}

BOOST_PYTHON_MODULE(rdMolHash) {
  python::scope().attr("__doc__") =
      "Module containing functions to generate a hash/key for molecules";

  std::string docString = "Generates a hash string for a molecule";
  python::def("GenerateMoleculeHashString", GenMolHashString,
              (python::arg("mol"), python::arg("atomsToUse") = python::list(),
               python::arg("bondsToUse") = python::list()),
              docString.c_str());
}
