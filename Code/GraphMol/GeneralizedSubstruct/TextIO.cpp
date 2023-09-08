//
//  Copyright (c) 2023, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include <DataStructs/base64.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <RDGeneral/StreamOps.h>
#include "XQMol.h"

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <cstdint>
#include <sstream>
#include <variant>

namespace bpt = boost::property_tree;

namespace RDKit {
namespace GeneralizedSubstruct {
namespace detail {
constexpr std::uint16_t recognition = 0xbe73;
constexpr std::uint16_t version = 1000;

}  // namespace detail
std::string ExtendedQueryMol::toBinary() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out);
  streamWrite(ss, detail::recognition);
  streamWrite(ss, detail::version);
  std::string pkl;
  if (std::holds_alternative<RWMol_T>(xqmol)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_MOL);
    MolPickler::pickleMol(*std::get<RWMol_T>(xqmol), pkl);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<MolBundle_T>(xqmol)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_MOLBUNDLE);
    pkl = std::get<MolBundle_T>(xqmol)->serialize();
  } else if (std::holds_alternative<TautomerQuery_T>(xqmol)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_TAUTOMERQUERY);
    pkl = std::get<TautomerQuery_T>(xqmol)->serialize();
  } else if (std::holds_alternative<TautomerBundle_T>(xqmol)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE);
    const auto &itm = std::get<TautomerBundle_T>(xqmol);
    std::uint16_t nTauts = itm->size();
    streamWrite(ss, nTauts);
    for (const auto &taut : *itm) {
      streamWrite(ss, taut->serialize());
    }
#endif
  } else {
    UNDER_CONSTRUCTION("unrecognized type in ExtendedQueryMol");
  }
  if (!pkl.empty()) {
    streamWrite(ss, pkl);
  }
  return ss.str();
}

namespace {

ExtendedQueryMol::TautomerBundle_T readTautomerQueries(std::stringstream &ss) {
  ExtendedQueryMol::TautomerBundle_T res{
      new std::vector<ExtendedQueryMol::TautomerQuery_T>()};
  std::uint16_t nTauts = 0;
  streamRead(ss, nTauts);
  res->reserve(nTauts);
  for (auto i = 0u; i < nTauts; ++i) {
    std::string pkl;
    streamRead(ss, pkl, 0);
    res->emplace_back(std::make_unique<TautomerQuery>(pkl));
  }
  return res;
}
}  // namespace

void ExtendedQueryMol::initFromBinary(const std::string &pickle) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::in |
                       std::ios_base::out);
  ss.write(pickle.c_str(), pickle.length());
  std::uint16_t recog;
  streamRead(ss, recog);

  if (recog != detail::recognition) {
    throw ValueErrorException("bad pickle format: bad endian ID");
  }
  std::uint16_t vers;
  streamRead(ss, vers);
  if (vers > detail::version) {
    throw ValueErrorException(
        "attempt to depickle from more recent pickle version");
  }
  unsigned char readType;
  streamRead(ss, readType);
  std::string pkl;
  switch (readType) {
    case ExtendedQueryMolTypes::XQM_MOL:
      streamRead(ss, pkl, 0);
      xqmol = std::make_unique<RWMol>(pkl);
      break;
#ifdef RDK_USE_BOOST_SERIALIZATION
    case ExtendedQueryMolTypes::XQM_MOLBUNDLE:
      streamRead(ss, pkl, 0);
      xqmol = std::make_unique<MolBundle>(pkl);
      break;
    case ExtendedQueryMolTypes::XQM_TAUTOMERQUERY:
      streamRead(ss, pkl, 0);
      xqmol = std::make_unique<TautomerQuery>(pkl);
      break;
    case ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE:
      xqmol = readTautomerQueries(ss);
      break;
#endif
    default:
      UNDER_CONSTRUCTION("unrecognized type in ExtendedQueryMol");
  }
}

namespace {
bool has_query_feature(const ROMol &mol) {
  for (const auto atom : mol.atoms()) {
    if (atom->hasQuery()) {
      return true;
    }
  }
  for (const auto bond : mol.bonds()) {
    if (bond->hasQuery()) {
      return true;
    }
  }
  return false;
}

void add_mol_to_elem(bpt::ptree &elem, const ROMol &mol) {
  std::string pkl;
  MolPickler::pickleMol(mol, pkl);
  std::unique_ptr<char[]> b64(Base64Encode(pkl.c_str(), pkl.length()));
  elem.put("pkl", b64.get());
  if (has_query_feature(mol)) {
    elem.put("smarts", MolToCXSmarts(mol));
  } else {
    elem.put("smiles", MolToCXSmiles(mol));
  }
}

void to_pt(bpt::ptree &pt, const ROMol &mol) {
  bpt::ptree elem;
  add_mol_to_elem(elem, mol);
  pt.add_child("mol", elem);
}

void to_pt(bpt::ptree &pt, const MolBundle &bndl) {
  {
    bpt::ptree children;
    for (const auto &mol : bndl.getMols()) {
      bpt::ptree elem;
      RWMol mcp(*mol);
      mcp.clearComputedProps();
      add_mol_to_elem(elem, mcp);
      children.push_back({"", elem});
    }
    pt.add_child("mols", children);
  }
}
void to_pt(bpt::ptree &pt, const TautomerQuery &tq) {
  {
    bpt::ptree children;
    for (const auto &taut : tq.getTautomers()) {
      bpt::ptree elem;
      add_mol_to_elem(elem, *taut);
      children.push_back({"", elem});
    }
    pt.add_child("tautomers", children);
  }
  {
    bpt::ptree elem;
    add_mol_to_elem(elem, tq.getTemplateMolecule());
    pt.add_child("template", elem);
  }
  {
    bpt::ptree children;
    for (const auto idx : tq.getModifiedAtoms()) {
      bpt::ptree elem;
      elem.put_value(idx);
      children.push_back({"", elem});
    }
    pt.add_child("modifiedAtoms", children);
  }
  {
    bpt::ptree children;
    for (const auto idx : tq.getModifiedBonds()) {
      bpt::ptree elem;
      elem.put_value(idx);
      children.push_back({"", elem});
    }
    pt.add_child("modifiedBonds", children);
  }
}

void to_pt(bpt::ptree &pt,
           const ExtendedQueryMol::TautomerBundle_T &tautQueries) {
  {
    bpt::ptree children;
    for (const auto &tq : *tautQueries) {
      bpt::ptree elem;
      to_pt(elem, *tq);
      children.push_back({"", elem});
    }
    pt.add_child("tautomerQueries", children);
  }
}

RWMol *pt_to_mol(bpt::ptree &pt) {
  auto b64pkl = pt.get<std::string>("pkl", "");
  if (!b64pkl.empty()) {
    unsigned int len;
    std::unique_ptr<char[]> cpkl(Base64Decode(b64pkl.c_str(), &len));
    std::string pkl(cpkl.get(), len);
    return new RWMol(pkl);
  }
  auto smi = pt.get<std::string>("smiles", "");
  if (!smi.empty()) {
    return SmilesToMol(smi);
  } else {
    auto sma = pt.get<std::string>("smarts");
    return SmartsToMol(sma);
  }
}

template <typename T>
T from_pt(bpt::ptree &) {}

template <>
ExtendedQueryMol::RWMol_T from_pt(bpt::ptree &pt) {
  return ExtendedQueryMol::RWMol_T(pt_to_mol(pt.get_child("mol")));
}

template <>
ExtendedQueryMol::TautomerQuery_T from_pt(bpt::ptree &pt) {
  std::vector<ROMOL_SPTR> tautomers;
  for (auto &child : pt.get_child("tautomers")) {
    tautomers.push_back(ROMOL_SPTR(pt_to_mol(child.second)));
  }

  ROMol *templ = pt_to_mol(pt.get_child("template"));
  templ->updatePropertyCache(false);

  std::vector<size_t> modifiedAtoms;
  for (auto &child : pt.get_child("modifiedAtoms")) {
    modifiedAtoms.push_back(child.second.get<size_t>(""));
  }

  std::vector<size_t> modifiedBonds;
  for (auto &child : pt.get_child("modifiedBonds")) {
    modifiedBonds.push_back(child.second.get<size_t>(""));
  }

  auto res = ExtendedQueryMol::TautomerQuery_T(
      new TautomerQuery(tautomers, templ, modifiedAtoms, modifiedBonds));
  return res;
}

template <>
ExtendedQueryMol::MolBundle_T from_pt(bpt::ptree &pt) {
  auto res = ExtendedQueryMol::MolBundle_T(new MolBundle);
  for (auto &child : pt.get_child("mols")) {
    res->addMol(ROMOL_SPTR(pt_to_mol(child.second)));
  }
  return res;
}

template <>
ExtendedQueryMol::TautomerBundle_T from_pt(bpt::ptree &pt) {
  ExtendedQueryMol::TautomerBundle_T res{
      new std::vector<std::unique_ptr<TautomerQuery>>()};
  for (auto &child : pt.get_child("tautomerQueries")) {
    res->emplace_back(from_pt<ExtendedQueryMol::TautomerQuery_T>(child.second));
  }
  return res;
}
}  // namespace

void ExtendedQueryMol::initFromJSON(const std::string &json) {
  std::istringstream ss;
  ss.str(json);
  try {
    bpt::ptree pt;
    bpt::read_json(ss, pt);
    auto xqmType = pt.get<unsigned char>("xqm_type");
    switch (xqmType) {
      case ExtendedQueryMol::ExtendedQueryMolTypes::XQM_MOL:
        xqmol = from_pt<ExtendedQueryMol::RWMol_T>(pt);
        break;
      case ExtendedQueryMol::ExtendedQueryMolTypes::XQM_MOLBUNDLE:
        xqmol = from_pt<ExtendedQueryMol::MolBundle_T>(pt);
        break;
      case ExtendedQueryMol::ExtendedQueryMolTypes::XQM_TAUTOMERQUERY:
        xqmol = from_pt<ExtendedQueryMol::TautomerQuery_T>(pt);
        break;
      case ExtendedQueryMol::ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE:
        xqmol = from_pt<ExtendedQueryMol::TautomerBundle_T>(pt);
        break;
      default:
        UNDER_CONSTRUCTION("unrecognized type in JSON");
    }
  } catch (const bpt::ptree_error &) {
    throw ValueErrorException("problems parsing JSON");
  }
}

std::string ExtendedQueryMol::toJSON() const {
  bpt::ptree pt;

  if (std::holds_alternative<RWMol_T>(xqmol)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_MOL);
    to_pt(pt, *std::get<RWMol_T>(xqmol));
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<MolBundle_T>(xqmol)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_MOLBUNDLE);
    to_pt(pt, *std::get<MolBundle_T>(xqmol));
  } else if (std::holds_alternative<TautomerQuery_T>(xqmol)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_TAUTOMERQUERY);
    to_pt(pt, *std::get<TautomerQuery_T>(xqmol));
  } else if (std::holds_alternative<TautomerBundle_T>(xqmol)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE);
    const auto &itm = std::get<TautomerBundle_T>(xqmol);
    pt.put("num_entries", itm->size());
    to_pt(pt, itm);
#endif
  } else {
    UNDER_CONSTRUCTION("unrecognized type in ExtendedQueryMol");
  }

  std::stringstream ss;
  bpt::json_parser::write_json(ss, pt);
  return ss.str();
}
}  // namespace GeneralizedSubstruct
}  // namespace RDKit
