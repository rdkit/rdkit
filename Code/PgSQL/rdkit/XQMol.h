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
#include <RDGeneral/export.h>
#ifndef XQMOL_H_MAY2023
#define XQMOL_H_MAY2023

#include <cstdint>
#include <variant>
#include <sstream>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
enum ExtendedQueryMolTypes : unsigned char {
  XQM_MOL = 1,
  XQM_MOLBUNDLE = 2,
  XQM_TAUTOMERQUERY = 3,
  XQM_TAUTOMERBUNDLE = 4
};
using ExtendedQueryMol =
    std::variant<std::unique_ptr<RWMol>, std::unique_ptr<MolBundle>,
                 std::unique_ptr<TautomerQuery>,
                 std::vector<std::unique_ptr<TautomerQuery>>>;

namespace detail {
constexpr std::uint16_t recognition = 0xbe73;
constexpr std::uint16_t version = 1000;

}  // namespace detail
std::string pickle(const ExtendedQueryMol &xqm) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out);
  streamWrite(ss, detail::recognition);
  streamWrite(ss, detail::version);
  std::string pkl;
  if (std::holds_alternative<std::unique_ptr<RWMol>>(xqm)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_MOL);
    MolPickler::pickleMol(*std::get<std::unique_ptr<RWMol>>(xqm), pkl);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<std::unique_ptr<MolBundle>>(xqm)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_MOLBUNDLE);
    pkl = std::get<std::unique_ptr<MolBundle>>(xqm)->serialize();
  } else if (std::holds_alternative<std::unique_ptr<TautomerQuery>>(xqm)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_TAUTOMERQUERY);
    pkl = std::get<std::unique_ptr<TautomerQuery>>(xqm)->serialize();
  } else if (std::holds_alternative<
                 std::vector<std::unique_ptr<TautomerQuery>>>(xqm)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE);
    const auto &itm =
        std::get<std::vector<std::unique_ptr<TautomerQuery>>>(xqm);
    std::uint16_t nTauts = itm.size();
    streamWrite(ss, nTauts);
    for (const auto &taut : itm) {
      streamWrite(ss, taut->serialize());
    }
#endif
  } else {
    throw ValueErrorException("unrecognized type in ExtendedQueryMol");
  }
  if (!pkl.empty()) {
    streamWrite(ss, pkl);
  }
  return ss.str();
}

namespace {

std::vector<std::unique_ptr<TautomerQuery>> readTautomerQueries(
    std::stringstream &ss) {
  std::vector<std::unique_ptr<TautomerQuery>> res;
  std::uint16_t nTauts = 0;
  streamRead(ss, nTauts);

  for (auto i = 0u; i < nTauts; ++i) {
    std::string pkl;
    streamRead(ss, pkl, 0);
    res.emplace_back(std::make_unique<TautomerQuery>(pkl));
  }
  return res;
}
}  // namespace

ExtendedQueryMol *depickle(const std::string &pickle) {
  ExtendedQueryMol *res = nullptr;
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
      res = new ExtendedQueryMol(std::make_unique<RWMol>(pkl));
      break;
#ifdef RDK_USE_BOOST_SERIALIZATION
    case ExtendedQueryMolTypes::XQM_MOLBUNDLE:
      streamRead(ss, pkl, 0);
      res = new ExtendedQueryMol(std::make_unique<MolBundle>(pkl));
      break;
    case ExtendedQueryMolTypes::XQM_TAUTOMERQUERY:
      streamRead(ss, pkl, 0);
      res = new ExtendedQueryMol(std::make_unique<TautomerQuery>(pkl));
      break;
    case ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE: {
      res = new ExtendedQueryMol(std::move(readTautomerQueries(ss)));
    } break;
#endif
    default:
      throw ValueErrorException("unknown type in pickle");
  }
  return res;
}

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
           const std::vector<std::unique_ptr<TautomerQuery>> &tautQueries) {
  {
    boost::property_tree::ptree children;
    for (const auto &tq : tautQueries) {
      boost::property_tree::ptree elem;
      to_pt(elem, *tq);
      children.push_back({"", elem});
    }
    pt.add_child("tautomerQueries", children);
  }
}
std::string to_text(const ExtendedQueryMol &xqm) {
  boost::property_tree::ptree pt;

  if (std::holds_alternative<std::unique_ptr<RWMol>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_MOL);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<std::unique_ptr<MolBundle>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_MOLBUNDLE);
    to_pt(pt, *std::get<std::unique_ptr<MolBundle>>(xqm));
  } else if (std::holds_alternative<std::unique_ptr<TautomerQuery>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_TAUTOMERQUERY);
    to_pt(pt, *std::get<std::unique_ptr<TautomerQuery>>(xqm));
  } else if (std::holds_alternative<
                 std::vector<std::unique_ptr<TautomerQuery>>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_TAUTOMERBUNDLE);
    const auto &itm =
        std::get<std::vector<std::unique_ptr<TautomerQuery>>>(xqm);
    pt.put("num_entries", itm.size());
    to_pt(pt, itm);
#endif
  } else {
    throw ValueErrorException("unrecognized type in ExtendedQueryMol");
  }

  std::stringstream ss;
  boost::property_tree::json_parser::write_json(ss, pt);
  return ss.str();
}

}  // namespace RDKit
#endif
