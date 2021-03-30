//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <cstring>
#include <iostream>

#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <DataStructs/BitOps.h>

#include "common.h"

#include <INCHI-API/inchi.h>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

using namespace RDKit;

namespace {
char *str_to_c(const std::string &str, size_t *len = nullptr) {
  if (len) {
    *len = str.size();
  }
  char *res;
  res = (char *)malloc(str.size() + 1);
  memcpy((void *)res, (const void *)str.c_str(), str.size());
  res[str.size()] = 0;
  return res;
}
char *str_to_c(const char *str) {
  char *res;
  res = (char *)malloc(strlen(str) + 1);
  strcpy(res, str);
  return res;
}
}  // namespace

#if 0
std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}
#endif

#define MOL_FROM_PKL(pkl, len) \
  if (!pkl || !len) {          \
    return str_to_c("");       \
  }                            \
  std::string lpkl(pkl, len);  \
  ROMol mol(lpkl);

extern "C" char *get_smiles(const char *pkl, size_t len) {
  MOL_FROM_PKL(pkl, len);
  auto data = MolToSmiles(mol);
  return str_to_c(data);
}
extern "C" char *get_cxsmiles(const char *pkl, size_t len) {
  MOL_FROM_PKL(pkl, len);
  auto data = MolToCXSmiles(mol);
  return str_to_c(data);
}
extern "C" char *get_molblock(const char *pkl, size_t len) {
  MOL_FROM_PKL(pkl, len);
  auto data = MolToMolBlock(mol);
  return str_to_c(data);
}
extern "C" char *get_v3kmolblock(const char *pkl, size_t len) {
  MOL_FROM_PKL(pkl, len);
  auto data = MolToV3KMolBlock(mol);
  return str_to_c(data);
}
extern "C" char *get_json(const char *pkl, size_t len) {
  MOL_FROM_PKL(pkl, len);
  auto data = MolInterchange::MolToJSONData(mol);
  return str_to_c(data);
}

#if 0
std::string get_cxsmiles(const char *pkl, size_t len) const {
  if (!d_mol) return "";
  return MolToCXSmiles(*d_mol);
}

std::string JSMol::get_svg(unsigned int w, unsigned int h) const {
  if (!d_mol) return "";
  return svg_(*d_mol, w, h);
}
std::string JSMol::get_svg_with_highlights(const std::string &details) const {
  if (!d_mol) return "";

  unsigned int w = d_defaultWidth;
  unsigned int h = d_defaultHeight;
  return svg_(*d_mol, w, h, details);
}

std::string JSMol::get_inchi() const {
  if (!d_mol) return "";
  ExtraInchiReturnValues rv;
  return MolToInchi(*d_mol, rv);
}
std::string JSMol::get_molblock() const {
  if (!d_mol) return "";
  return MolToMolBlock(*d_mol);
}
std::string JSMol::get_v3Kmolblock() const {
  if (!d_mol) return "";
  return MolToV3KMolBlock(*d_mol);
}
#endif

extern "C" char *get_mol(const char *input, size_t *len) {
  RWMol *mol = MinimalLib::mol_from_input(input);
  if (!mol) {
    *len = 0;
    return str_to_c("Error!");
  }
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl);
  return str_to_c(pkl, len);
}
#if 0
JSMol *get_qmol(const std::string &input) {
  RWMol *mol = qmol_from_input(input);
  return new JSMol(mol);
}
#endif
extern "C" char *version() { return str_to_c(rdkitVersion); }

extern "C" void prefer_coordgen(bool useCoordGen) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = useCoordGen;
#endif
}
