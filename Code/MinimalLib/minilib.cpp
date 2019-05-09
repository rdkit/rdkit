//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <iostream>
#include "minilib.h"

#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <INCHI-API/inchi.h>

using namespace RDKit;

namespace {
ROMol *mol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    res = MolBlockToMol(input, sanitize);
  } else {
    SmilesParserParams ps;
    ps.sanitize = false;
    res = SmilesToMol(input, ps);
  }
  if (res) {
    try {
      MolOps::sanitizeMol(*res);
      MolOps::assignStereochemistry(*res, true, true, true);
    } catch (...) {
      delete res;
      res = nullptr;
    }
  }
  return res;
}

ROMol *qmol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    res = MolBlockToMol(input, sanitize);
  } else {
    res = SmartsToMol(input);
  }
  return res;
}

std::string smiles_(const ROMol &m) { return MolToSmiles(m); }
std::string molblock_(const ROMol &m) { return MolToMolBlock(m); }
std::string inchi_(const ROMol &m) {
  ExtraInchiReturnValues rv;
  return MolToInchi(m, rv);
}
std::string svg_(const ROMol &m) {
  MolDraw2DSVG drawer(250, 200);
  MolDraw2DUtils::prepareAndDrawMolecule(drawer, m);
  drawer.finishDrawing();

  return drawer.getDrawingText();
}
}  // namespace

std::string JSMol::get_smiles() const {
  if (!d_mol) return "";
  return smiles_(*d_mol);
}
std::string JSMol::get_svg() const {
  if (!d_mol) return "";
  return svg_(*d_mol);
}
std::string JSMol::get_inchi() const {
  if (!d_mol) return "";
  return inchi_(*d_mol);
}
std::string JSMol::get_molblock() const {
  if (!d_mol) return "";
  return molblock_(*d_mol);
}
std::vector<unsigned int> JSMol::get_substruct_match(const JSMol &q) const {
  std::vector<unsigned int> res;
  MatchVectType match;
  if (SubstructMatch(*d_mol, *(q.d_mol), match)) {
    for (const auto &pr : match) {
      res.push_back(pr.second);
    }
  }
  return res;
}

std::string get_smiles(const std::string &input) {
  std::unique_ptr<ROMol> mol(mol_from_input(input));
  if (!mol) return "";
  return smiles_(*mol);
}
std::string get_svg(const std::string &input) {
  std::unique_ptr<ROMol> mol(mol_from_input(input));
  if (!mol) return "";
  return svg_(*mol);
}

std::string get_inchi(const std::string &input) {
  std::unique_ptr<ROMol> mol(mol_from_input(input));
  if (!mol) return "";
  return inchi_(*mol);
}

std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

JSMol *get_mol(const std::string &input) {
  ROMol *mol = mol_from_input(input);
  return new JSMol(mol);
}

JSMol *get_qmol(const std::string &input) {
  ROMol *mol = qmol_from_input(input);
  return new JSMol(mol);
}

std::string version() { return std::string(rdkitVersion); }

int ping() {
  std::cerr << "alive" << std::endl;
  return 1;
}
