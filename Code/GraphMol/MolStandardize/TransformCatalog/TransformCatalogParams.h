//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_TRANSFORM_CATALOG_PARAMS_H
#define RD_TRANSFORM_CATALOG_PARAMS_H

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

  // if the string here is empty the default transforms will be used
  TransformCatalogParams(const std::string &transformFile);
  TransformCatalogParams(std::istream &transformStream);
  TransformCatalogParams(
      const std::vector<std::pair<std::string, std::string>> &data);
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
