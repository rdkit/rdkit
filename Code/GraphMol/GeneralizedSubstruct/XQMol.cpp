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

#include <cstdint>
#include <variant>
#include <sstream>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

#include <GraphMol/TautomerQuery/TautomerQuery.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "XQMol.h"

namespace RDKit {

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

void ExtendedQueryMol::initFromJSON(const std::string &) {
  UNDER_CONSTRUCTION("not yet implemented");
}

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

ExtendedQueryMol::ExtendedQueryMol(const std::string &text, bool isJSON) {
  if (!isJSON) {
    initFromBinary(text);
  } else {
    initFromJSON(text);
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

void to_pt(boost::property_tree::ptree &pt, const MolBundle &bndl) {
  {
    boost::property_tree::ptree children;
    for (const auto &mol : bndl.getMols()) {
      boost::property_tree::ptree elem;
      if (has_query_feature(*mol)) {
        elem.put_value(MolToCXSmarts(*mol));
      } else {
        elem.put_value(MolToCXSmiles(*mol));
      }
      children.push_back({"", elem});
    }
    pt.add_child("mols", children);
  }
}
void to_pt(boost::property_tree::ptree &pt, const TautomerQuery &tq) {
  {
    boost::property_tree::ptree children;
    for (const auto &taut : tq.getTautomers()) {
      boost::property_tree::ptree elem;
      if (has_query_feature(*taut)) {
        elem.put_value(MolToCXSmarts(*taut));
      } else {
        elem.put_value(MolToCXSmiles(*taut));
      }

      children.push_back({"", elem});
    }
    pt.add_child("tautomers", children);
  }
  {
    boost::property_tree::ptree children;
    for (const auto idx : tq.getModifiedAtoms()) {
      boost::property_tree::ptree elem;
      elem.put_value(idx);
      children.push_back({"", elem});
    }
    pt.add_child("modifiedAtoms", children);
  }
  {
    boost::property_tree::ptree children;
    for (const auto idx : tq.getModifiedBonds()) {
      boost::property_tree::ptree elem;
      elem.put_value(idx);
      children.push_back({"", elem});
    }
    pt.add_child("modifiedBonds", children);
  }
  pt.put("template", MolToCXSmarts(tq.getTemplateMolecule()));
}

void to_pt(boost::property_tree::ptree &pt,
           const ExtendedQueryMol::TautomerBundle_T &tautQueries) {
  {
    boost::property_tree::ptree children;
    for (const auto &tq : *tautQueries) {
      boost::property_tree::ptree elem;
      to_pt(elem, *tq);
      children.push_back({"", elem});
    }
    pt.add_child("tautomerQueries", children);
  }
}
}  // namespace

std::string ExtendedQueryMol::toJSON() const {
  boost::property_tree::ptree pt;

  if (std::holds_alternative<RWMol_T>(xqmol)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_MOL);
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
  boost::property_tree::json_parser::write_json(ss, pt);
  return ss.str();
}

std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  if (std::holds_alternative<ExtendedQueryMol::RWMol_T>(query.xqmol)) {
    res = RDKit::SubstructMatch(
        mol, *std::get<ExtendedQueryMol::RWMol_T>(query.xqmol), params);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<ExtendedQueryMol::MolBundle_T>(
                 query.xqmol)) {
    res = RDKit::SubstructMatch(
        mol, *std::get<ExtendedQueryMol::MolBundle_T>(query.xqmol), params);
  } else if (std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(
                 query.xqmol)) {
    res = std::get<ExtendedQueryMol::TautomerQuery_T>(query.xqmol)
              ->substructOf(mol, params);
  } else if (std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(
                 query.xqmol)) {
    const auto &vect =
        std::get<ExtendedQueryMol::TautomerBundle_T>(query.xqmol);
    for (const auto &tq : *vect) {
      res = tq->substructOf(mol, params);
      if (!res.empty()) {
        break;
      }
    }
#endif
  } else {
    UNDER_CONSTRUCTION("unrecognized type in ExtendedQueryMol");
  }
  return res;
}
}  // namespace RDKit
