//
//  Copyright (C) 2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <sstream>
#include <vector>
#include <string>
namespace pt = boost::property_tree;

namespace RDKit {
namespace {
void ptreeToMol(RWMol *mol, const pt::ptree &molE) {
  PRECONDITION(mol, "no molecule");
  for (const auto &atE : molE) {
    if (atE.first == "rdkit:atom") {
      std::string asmi = atE.second.get<std::string>("<xmlattr>.atom-smiles");
      Atom *atom = SmilesToAtom(asmi);
      bool updateLabel = false, takeOwnership = true;
      mol->addAtom(atom, updateLabel, takeOwnership);
    }
  }
  for (const auto &atE : molE) {
    if (atE.first == "rdkit:bond") {
      std::string asmi = atE.second.get<std::string>("<xmlattr>.bond-smiles");
      Bond *bond = SmilesToBond(asmi);
      bond->setBeginAtomIdx(atE.second.get<int>("<xmlattr>.begin-atom-idx") -
                            1);
      bond->setEndAtomIdx(atE.second.get<int>("<xmlattr>.end-atom-idx") - 1);
      bool takeOwnership = true;
      mol->addBond(bond, takeOwnership);
    }
  }
}
}  // namespace

RWMol *RDKitSVGToMol(const std::string &svg, bool sanitize, bool removeHs) {
  std::stringstream iss(svg);
  pt::ptree tree;
  pt::read_xml(iss, tree);
  RWMol *res = nullptr;
  pt::ptree empty_ptree;
  const pt::ptree &molsE = tree.get_child("svg", empty_ptree);
  for (const auto &molE : molsE) {
    if (molE.first == "rdkit:mol") {
      res = new RWMol();
      ptreeToMol(res, molE.second);
      if (res->getNumAtoms()) {
        if (removeHs) {
          bool implicitOnly = false, updateExplicitCount = false;
          MolOps::removeHs(*res, implicitOnly, updateExplicitCount, sanitize);
        } else if (sanitize) {
          MolOps::sanitizeMol(*res);
        }
      }
    }
  }
  return res;
}
}  // namespace RDKit
