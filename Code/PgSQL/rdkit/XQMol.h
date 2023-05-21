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
  XQM_BUNDLE = 2,
  XQM_TAUTOMERQUERY = 3
};
using ExtendedQueryMol =
    std::variant<std::unique_ptr<RWMol>, std::unique_ptr<MolBundle>,
                 std::unique_ptr<TautomerQuery>>;

namespace detail {
constexpr std::uint16_t recognition = 0xfeed;
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
  } else if (std::holds_alternative<std::unique_ptr<MolBundle>>(xqm)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_BUNDLE);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<std::unique_ptr<TautomerQuery>>(xqm)) {
    streamWrite(ss, ExtendedQueryMolTypes::XQM_TAUTOMERQUERY);
    pkl = std::get<std::unique_ptr<TautomerQuery>>(xqm)->serialize();
#endif
  } else {
    throw ValueErrorException("unrecognized type in ExtendedQueryMol");
  }
  streamWrite(ss, pkl);
  return ss.str();
}

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
      res = new ExtendedQueryMol(std::make_unique<RWMol>(pickle));
      break;
    case ExtendedQueryMolTypes::XQM_BUNDLE:
      break;
#ifdef RDK_USE_BOOST_SERIALIZATION
    case ExtendedQueryMolTypes::XQM_TAUTOMERQUERY:
      streamRead(ss, pkl, 0);
      res = new ExtendedQueryMol(std::make_unique<TautomerQuery>(pkl));
      break;
#endif
    default:
      throw ValueErrorException("unknown type in pickle");
  }
  return res;
}

std::string to_text(const ExtendedQueryMol &xqm) {
  boost::property_tree::ptree pt;

  if (std::holds_alternative<std::unique_ptr<RWMol>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_MOL);
  } else if (std::holds_alternative<std::unique_ptr<MolBundle>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_BUNDLE);
#ifdef RDK_USE_BOOST_SERIALIZATION
  } else if (std::holds_alternative<std::unique_ptr<TautomerQuery>>(xqm)) {
    pt.put("xqm_type", (int)ExtendedQueryMolTypes::XQM_TAUTOMERQUERY);
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
