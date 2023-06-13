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

ExtendedQueryMol createExtendedQueryMol(const RWMol &mol,
                                        bool adjustQueryProperties,
                                        MolOps::AdjustQueryParameters params) {
  auto bndl = MolEnumerator::enumerate(mol);
  if (bndl.empty()) {
    const ROMol *lmol = nullptr;
    std::unique_ptr<ROMol> holder;
    if (adjustQueryProperties) {
      holder.reset(MolOps::adjustQueryProperties(mol, &params));
      lmol = holder.get();
    } else {
      lmol = &mol;
    }
    // nothing enumerated
    auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*lmol));
    if (tq->getTautomers().size() == 1) {
      // no tautomers, return just the mol
      return {std::make_unique<RWMol>(*lmol)};
    } else {
      // return the tautomers
      return {std::move(tq)};
    }
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
    auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*lmol));
    if (tq->getTautomers().size() == 1) {
      // no tautomers, just one molecule, return the molecule:
      return {std::make_unique<RWMol>(*lmol)};
    } else {
      // return the tautomers
      return {std::move(tq)};
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
    for (const auto &bmol : lbndl.getMols()) {
      auto tq = std::unique_ptr<TautomerQuery>(TautomerQuery::fromMol(*bmol));
      if (tq->getTautomers().size() > 1) {
        hadTautomers = true;
      }
      tautomerBundle->emplace_back(std::move(tq));
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
