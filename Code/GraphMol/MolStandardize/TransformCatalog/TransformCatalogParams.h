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
#ifndef __RD_TRANSFORM_CATALOG_PARAMS_H__
#define __RD_TRANSFORM_CATALOG_PARAMS_H__

#include <Catalogs/CatalogParams.h>
#include "TransformCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <string>
#include <vector>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class RDKIT_MOLSTANDARDIZE_EXPORT TransformCatalogParams
    : public RDCatalog::CatalogParams {
 public:
  TransformCatalogParams() {
    d_typeStr = "Transform Catalog Parameters";
    d_transformations.clear();
  }

  TransformCatalogParams(const std::string &transformFile);
  TransformCatalogParams(std::istream &transformStream);
  // copy constructor
  TransformCatalogParams(const TransformCatalogParams &other);

  ~TransformCatalogParams() override;

  unsigned int getNumTransformations() const {
    return static_cast<unsigned int>(d_transformations.size());
  }

  const std::vector<std::shared_ptr<ChemicalReaction>> &getTransformations()
      const;

  const ChemicalReaction *getTransformation(unsigned int fid) const;

  void toStream(std::ostream &) const override;
  std::string Serialize() const override;
  void initFromStream(std::istream &ss) override;
  void initFromString(const std::string &text) override;

 private:
  std::vector<std::shared_ptr<ChemicalReaction>> d_transformations;

};  // class TransformCatalogParams

}  // namespace MolStandardize
}  // namespace RDKit

#endif
