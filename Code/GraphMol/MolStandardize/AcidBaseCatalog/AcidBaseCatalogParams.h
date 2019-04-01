//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_ACIDBASE_CATALOG_PARAMS_H__
#define __RD_ACIDBASE_CATALOG_PARAMS_H__

#include <Catalogs/CatalogParams.h>
#include "AcidBaseCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <string>
#include <vector>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class RDKIT_MOLSTANDARDIZE_EXPORT AcidBaseCatalogParams
    : public RDCatalog::CatalogParams {
 public:
  AcidBaseCatalogParams() {
    d_typeStr = "AcidBase Catalog Parameters";
    d_pairs.clear();
  }

  AcidBaseCatalogParams(const std::string &acidBaseFile);
  AcidBaseCatalogParams(std::istream &acidBaseFile);
  // copy constructor
  AcidBaseCatalogParams(const AcidBaseCatalogParams &other);

  ~AcidBaseCatalogParams() override;

  unsigned int getNumPairs() const {
    return static_cast<unsigned int>(d_pairs.size());
  }

  const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &getPairs() const;

  const std::pair<ROMOL_SPTR, ROMOL_SPTR> getPair(unsigned int fid) const;

  void toStream(std::ostream &) const override;
  std::string Serialize() const override;
  void initFromStream(std::istream &ss) override;
  void initFromString(const std::string &text) override;

 private:
  std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> d_pairs;

};  // class AcidBaseCatalogParams

}  // namespace MolStandardize
}  // namespace RDKit

#endif
