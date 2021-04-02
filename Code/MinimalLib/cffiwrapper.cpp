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
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
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

#define MOL_FROM_PKL(mol, pkl, pkl_sz)       \
  if (!(pkl) || !(pkl_sz)) {                 \
    return NULL;                             \
  }                                          \
  std::string mol##inp_str((pkl), (pkl_sz)); \
  ROMol mol(mol##inp_str);                   \
  mol.setProp(common_properties::_StereochemDone, 1, true);

extern "C" char *get_smiles(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  auto data = MolToSmiles(mol);
  return str_to_c(data);
}
extern "C" char *get_smarts(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  auto data = MolToSmarts(mol);
  return str_to_c(data);
}
extern "C" char *get_cxsmiles(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  auto data = MolToCXSmiles(mol);
  return str_to_c(data);
}
extern "C" char *get_molblock(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  auto data = MolToMolBlock(mol);
  return str_to_c(data);
}
extern "C" char *get_v3kmolblock(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  auto data = MolToV3KMolBlock(mol);
  return str_to_c(data);
}
extern "C" char *get_json(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  auto data = MolInterchange::MolToJSONData(mol);
  return str_to_c(data);
}
extern "C" void free_ptr(char *ptr) {
  if (ptr) {
    free(ptr);
  }
}

extern "C" char *get_svg(const char *pkl, size_t pkl_sz, unsigned int width,
                         unsigned int height) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  return str_to_c(MinimalLib::mol_to_svg(mol, width, height));
}
extern "C" char *get_svg_with_highlights(const char *pkl, size_t pkl_sz,
                                         const char *details_json) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  unsigned int width = MinimalLib::d_defaultWidth;
  unsigned int height = MinimalLib::d_defaultHeight;
  return str_to_c(MinimalLib::mol_to_svg(mol, width, height, details_json));
}

extern "C" char *get_inchi(const char *pkl, size_t pkl_sz) {
  MOL_FROM_PKL(mol, pkl, pkl_sz);
  ExtraInchiReturnValues rv;
  return str_to_c(MolToInchi(mol, rv));
}

extern "C" char *get_inchi_for_molblock(const char *ctab) {
  if (!ctab) {
    return str_to_c("");
  }
  ExtraInchiReturnValues rv;
  return str_to_c(MolBlockToInchi(ctab, rv));
}

extern "C" char *get_inchikey_for_inchi(const char *inchi) {
  if (!inchi) {
    return str_to_c("");
  }
  return str_to_c(InchiToInchiKey(inchi));
}

extern "C" char *get_mol(const char *input, size_t *pkl_sz) {
  RWMol *mol = MinimalLib::mol_from_input(input);
  if (!mol) {
    *pkl_sz = 0;
    return NULL;
  }
  unsigned int propFlags = PicklerOps::PropertyPickleOptions::AllProps ^
                           PicklerOps::PropertyPickleOptions::ComputedProps;
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl, propFlags);
  return str_to_c(pkl, pkl_sz);
}
extern "C" char *get_qmol(const char *input, size_t *pkl_sz) {
  RWMol *mol = MinimalLib::qmol_from_input(input);
  if (!mol) {
    *pkl_sz = 0;
    return str_to_c("Error!");
  }
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl);
  return str_to_c(pkl, pkl_sz);
}
extern "C" char *version() { return str_to_c(rdkitVersion); }

extern "C" char *get_substruct_match(const char *mol_pkl, size_t mol_pkl_sz,
                                     const char *query_pkl, size_t query_pkl_sz,
                                     const char *options_json) {
  MOL_FROM_PKL(mol, mol_pkl, mol_pkl_sz)
  MOL_FROM_PKL(query, query_pkl, query_pkl_sz)

  SubstructMatchParameters params;
  if (options_json) {
    std::string json(options_json);
    updateSubstructMatchParamsFromJSON(params, json);
  }
  params.maxMatches = 1;

  std::string res = "{}";
  auto matches = SubstructMatch(mol, query, params);
  if (!matches.empty()) {
    auto match = matches[0];
    rj::Document doc;
    doc.SetObject();
    MinimalLib::get_sss_json(mol, query, match, doc, doc);
    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return str_to_c(res);
}
extern "C" char *get_substruct_matches(const char *mol_pkl, size_t mol_pkl_sz,
                                       const char *query_pkl,
                                       size_t query_pkl_sz,
                                       const char *options_json) {
  MOL_FROM_PKL(mol, mol_pkl, mol_pkl_sz)
  MOL_FROM_PKL(query, query_pkl, query_pkl_sz)
  SubstructMatchParameters params;
  if (options_json) {
    std::string json(options_json);
    updateSubstructMatchParamsFromJSON(params, json);
  }

  std::string res = "{}";
  auto matches = SubstructMatch(mol, query, params);
  if (!matches.empty()) {
    rj::Document doc;
    doc.SetArray();

    for (const auto &match : matches) {
      rj::Value rjMatch(rj::kObjectType);
      MinimalLib::get_sss_json(mol, query, match, rjMatch, doc);
      doc.PushBack(rjMatch, doc.GetAllocator());
    }

    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return str_to_c(res);
}

extern "C" char *get_descriptors(const char *mol_pkl, size_t mol_pkl_sz) {
  MOL_FROM_PKL(mol, mol_pkl, mol_pkl_sz)
  return str_to_c(MinimalLib::get_descriptors(mol));
}

extern "C" char *get_morgan_fp(const char *mol_pkl, size_t mol_pkl_sz,
                               unsigned int radius, size_t fplen) {
  MOL_FROM_PKL(mol, mol_pkl, mol_pkl_sz)
  auto fp = MorganFingerprints::getFingerprintAsBitVect(mol, radius, fplen);
  auto res = BitVectToText(*fp);
  delete fp;

  return str_to_c(res);
}

extern "C" char *get_rdkit_fp(const char *mol_pkl, size_t mol_pkl_sz,
                              size_t fplen) {
  MOL_FROM_PKL(mol, mol_pkl, mol_pkl_sz)
  auto fp = RDKFingerprintMol(mol, 1, 7, fplen);
  auto res = BitVectToText(*fp);
  delete fp;

  return str_to_c(res);
}

extern "C" void prefer_coordgen(short val) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = val;
#endif
};

extern "C" short set_2d_coords(char **mol_pkl, size_t *mol_pkl_sz) {
  MOL_FROM_PKL(mol, *mol_pkl, *mol_pkl_sz)
  RDDepict::compute2DCoords(mol);

  unsigned int propFlags = PicklerOps::PropertyPickleOptions::AllProps ^
                           PicklerOps::PropertyPickleOptions::ComputedProps;
  std::string pkl;
  MolPickler::pickleMol(mol, pkl, propFlags);

  free(*mol_pkl);
  *mol_pkl = str_to_c(pkl, mol_pkl_sz);
  return 1;
}
