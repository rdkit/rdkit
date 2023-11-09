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
namespace GeneralizedSubstruct {

ExtendedQueryMol::ExtendedQueryMol(const std::string &text, bool isJSON) {
  if (!isJSON) {
    initFromBinary(text);
  } else {
    initFromJSON(text);
  }
}

void ExtendedQueryMol::initFromOther(const ExtendedQueryMol &other) {
  if (std::holds_alternative<ExtendedQueryMol::RWMol_T>(other.xqmol)) {
    xqmol = std::make_unique<RWMol>(
        *std::get<ExtendedQueryMol::RWMol_T>(other.xqmol));
  } else if (std::holds_alternative<ExtendedQueryMol::MolBundle_T>(
                 other.xqmol)) {
    xqmol = std::make_unique<MolBundle>(
        *std::get<ExtendedQueryMol::MolBundle_T>(other.xqmol));
  } else if (std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(
                 other.xqmol)) {
    xqmol = std::make_unique<TautomerQuery>(
        *std::get<ExtendedQueryMol::TautomerQuery_T>(other.xqmol));
  } else if (std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(
                 other.xqmol)) {
    auto tb = std::make_unique<std::vector<std::unique_ptr<TautomerQuery>>>();
    for (const auto &tqp :
         *std::get<ExtendedQueryMol::TautomerBundle_T>(other.xqmol)) {
      tb->emplace_back(std::make_unique<TautomerQuery>(*tqp));
    }
    xqmol = std::move(tb);
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

ExtendedQueryMol createExtendedQueryMol(const RWMol &mol, bool doEnumeration,
                                        bool doTautomers,
                                        bool adjustQueryProperties,
                                        MolOps::AdjustQueryParameters params) {
  MolBundle bndl;
  if (doEnumeration) {
    bndl = MolEnumerator::enumerate(mol);
  }
  if (bndl.empty()) {
    // nothing enumerated, just add the input molecule
    bndl.addMol(boost::shared_ptr<ROMol>(new ROMol(mol)));
  }

  if (bndl.size() == 1) {
    const ROMol *lmol = nullptr;
    std::unique_ptr<ROMol> holder;
    if (adjustQueryProperties) {
      holder.reset(MolOps::adjustQueryProperties(*bndl.getMol(0), &params));
      lmol = holder.get();
    } else {
      lmol = bndl.getMol(0).get();
    }
    if (doTautomers) {
      auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*lmol));
      if (tq->getTautomers().size() == 1) {
        // no tautomers, just one molecule, return the molecule:
        return {std::make_unique<RWMol>(*lmol)};
      } else {
        // return the tautomers
        return {std::move(tq)};
      }
    } else {
      return {std::make_unique<RWMol>(*lmol)};
    }
  } else {
    MolBundle lbndl;
    for (auto &bmol : bndl.getMols()) {
      if (adjustQueryProperties) {
        boost::shared_ptr<ROMol> lmol(
            MolOps::adjustQueryProperties(*bmol, &params));
        lbndl.addMol(lmol);
      } else {
        lbndl.addMol(bmol);
      }
    }
    bool hadTautomers = false;
    auto tautomerBundle =
        std::make_unique<std::vector<std::unique_ptr<TautomerQuery>>>(0);
    if (doTautomers) {
      for (const auto &bmol : lbndl.getMols()) {
        auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*bmol));
        if (tq->getTautomers().size() > 1) {
          hadTautomers = true;
        }
        tautomerBundle->emplace_back(std::move(tq));
      }
    }
    if (!hadTautomers) {
      // no tautomers, just return the bundle
      return {std::make_unique<MolBundle>(lbndl)};
    } else {
      // return the tautomer bundle
      return {std::move(tautomerBundle)};
    }
  }
}
}  // namespace GeneralizedSubstruct
}  // namespace RDKit
