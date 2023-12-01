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

#include <variant>
#include <memory>
#include <string>
#include <vector>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/core/noncopyable.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/TautomerQuery/TautomerQuery.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit {
namespace GeneralizedSubstruct {
struct RDKIT_GENERALIZEDSUBSTRUCT_EXPORT ExtendedQueryMol {
  enum ExtendedQueryMolTypes : unsigned char {
    XQM_MOL = 1,
    XQM_MOLBUNDLE = 2,
    XQM_TAUTOMERQUERY = 3,
    XQM_TAUTOMERBUNDLE = 4
  };
  using RWMol_T = std::unique_ptr<RWMol>;
  using MolBundle_T = std::unique_ptr<MolBundle>;
  using TautomerQuery_T = std::unique_ptr<TautomerQuery>;
  using TautomerBundle_T =
      std::unique_ptr<std::vector<std::unique_ptr<TautomerQuery>>>;
  using ContainedType =
      std::variant<RWMol_T, MolBundle_T, TautomerQuery_T, TautomerBundle_T>;
  ExtendedQueryMol(std::unique_ptr<RWMol> mol) : xqmol(std::move(mol)) {}
  ExtendedQueryMol(std::unique_ptr<MolBundle> bundle)
      : xqmol(std::move(bundle)) {}
  ExtendedQueryMol(std::unique_ptr<TautomerQuery> tq) : xqmol(std::move(tq)) {}
  ExtendedQueryMol(
      std::unique_ptr<std::vector<std::unique_ptr<TautomerQuery>>> tqs)
      : xqmol(std::move(tqs)) {}
  ExtendedQueryMol(const ExtendedQueryMol &other) { initFromOther(other); }
  ExtendedQueryMol &operator=(const ExtendedQueryMol &other) {
    if (this == &other) {
      return *this;
    }
    initFromOther(other);
    return *this;
  }

  ExtendedQueryMol(ExtendedQueryMol &&o) noexcept : xqmol(std::move(o.xqmol)) {}
  ExtendedQueryMol(const std::string &text, bool isJSON = false);

  void initFromBinary(const std::string &pkl);
  void initFromJSON(const std::string &text);
  void initFromOther(const ExtendedQueryMol &other);

  ContainedType xqmol;
  std::string toBinary() const;
  std::string toJSON() const;
};

//! Creates an ExtendedQueryMol from the input molecule
/*!
  This takes a query molecule and, conceptually, performs the following steps to
  produce an ExtendedQueryMol:

    1. Enumerates features like Link Nodes and SRUs
    2. Converts everything into TautomerQueries
    3. Runs adjustQueryProperties()

  Each step is optional

    \param mol the molecule to start with
    \param doEnumeration  enumerate features like Link Nodes and SRUs
    \param doTautomers generate TautomerQueries
    \param adjustQueryProperties call adjustQueryProperties on each of the
       results
    \param params  AdjustQueryParameters object controlling the operation of
       adjustQueryProperties

    \return The new ExtendedQueryMol

*/
RDKIT_GENERALIZEDSUBSTRUCT_EXPORT ExtendedQueryMol createExtendedQueryMol(
    const RWMol &mol, bool doEnumeration = true, bool doTautomers = true,
    bool adjustQueryProperties = false,
    MolOps::AdjustQueryParameters params = {});

//! does a substructure search with an ExtendedQueryMol
RDKIT_GENERALIZEDSUBSTRUCT_EXPORT std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParameters &params = SubstructMatchParameters());

//! checks if a molecule has a match to an ExtendedQueryMol
inline bool hasSubstructMatch(
    const ROMol &mol, const ExtendedQueryMol &query,
    const SubstructMatchParameters &params = SubstructMatchParameters()) {
  SubstructMatchParameters lparams = params;
  lparams.maxMatches = 1;
  return !SubstructMatch(mol, query, lparams).empty();
}
}  // namespace GeneralizedSubstruct
}  // namespace RDKit
#endif
