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
#include <INCHI-API/inchi.h>

using namespace RDKit;

namespace {
std::unique_ptr<ROMol> mol_from_input(const std::string &input) {
  RWMol *res = nullptr;
#if 0
  try {  // pickle?
    res = new RWMol(ROMol(input));
  } catch (MolPicklerException &) {
    res = nullptr;
  }
#endif
  if (!res) {
    SmilesParserParams ps;
    ps.sanitize = false;
    res = SmilesToMol(input, ps);
  }
  if (!res) {
    bool sanitize = false;
    res = MolBlockToMol(input, sanitize);
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
  return std::unique_ptr<ROMol>(res);
}
}  // namespace

std::string get_smiles(const std::string &input) {
  std::unique_ptr<ROMol> mol = mol_from_input(input);
  if (!mol) return "";
  return MolToSmiles(*mol);
}

std::string get_svg(const std::string &input) {
  std::unique_ptr<ROMol> mol = mol_from_input(input);
  if (!mol) return "";

  MolDraw2DSVG drawer(250, 200);
  MolDraw2DUtils::prepareAndDrawMolecule(drawer, *mol);
  drawer.finishDrawing();

  return drawer.getDrawingText();
}

std::string get_inchi(const std::string &input) {
  std::unique_ptr<ROMol> mol = mol_from_input(input);
  if (!mol) return "";
  ExtraInchiReturnValues rv;
  return MolToInchi(*mol, rv);
}

std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

std::string get_pkl(const std::string &input) {
  std::unique_ptr<ROMol> mol = mol_from_input(input);
  if (!mol) return "";
  std::string res;
  MolPickler::pickleMol(*mol, res);
  return res;
}

std::string version() { return std::string(rdkitVersion); }

int ping() {
  std::cerr << "alive" << std::endl;
  return 1;
}
