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

namespace RDKit {
class RWMol;
class MolBundle;
class TautomerQuery;

struct RDKIT_GENERALIZEDSUBSTRUCT_EXPORT ExtendedQueryMol
    : private boost::noncopyable {
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

  ExtendedQueryMol(ExtendedQueryMol &&o) noexcept : xqmol(std::move(o.xqmol)) {}
  ExtendedQueryMol(const std::string &text, bool isJSON = false);

  void initFromBinary(const std::string &pkl);
  void initFromJSON(const std::string &text);

  ContainedType xqmol;
  std::string toBinary() const;
  std::string toJSON() const;
};

}  // namespace RDKit
#endif
