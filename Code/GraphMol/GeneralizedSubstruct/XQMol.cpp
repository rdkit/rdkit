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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>

#include "XQMol.h"

namespace RDKit {

ExtendedQueryMol::ExtendedQueryMol(const std::string &text, bool isJSON) {
  if (!isJSON) {
    initFromBinary(text);
  } else {
    initFromJSON(text);
  }
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

ExtendedQueryMol createExtendedQueryMol(const RWMol &mol) {
  auto bndl = MolEnumerator::enumerate(mol);
  if (bndl.empty()) {
    // nothing enumerated
    auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(mol));
    if (tq->getTautomers().size() == 1) {
      // no tautomers, return just the mol
      return {std::make_unique<RWMol>(mol)};
    } else {
      // return the tautomers
      return {std::move(tq)};
    }
    bndl.addMol(boost::shared_ptr<ROMol>(new ROMol(mol)));
  }

  if (bndl.size() == 1) {
    auto tq =
        std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*bndl.getMol(0)));
    if (tq->getTautomers().size() == 1) {
      // no tautomers, just one molecule, return the molecule:
      return {std::make_unique<RWMol>(*bndl.getMol(0))};
    } else {
      // return the tautomers
      return {std::move(tq)};
    }
  } else {
    bool hadTautomers = false;
    auto tautomerBundle =
        std::make_unique<std::vector<std::unique_ptr<TautomerQuery>>>(0);
    for (const auto &bmol : bndl.getMols()) {
      auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*bmol));
      if (tq->getTautomers().size() > 1) {
        hadTautomers = true;
      }
      tautomerBundle->emplace_back(std::move(tq));
    }
    if (!hadTautomers) {
      // no tautomers, just return the bundle
      return {std::make_unique<MolBundle>(bndl)};
    } else {
      // return the tautomer bundle
      return {std::move(tautomerBundle)};
    }
  }
}

}  // namespace RDKit
