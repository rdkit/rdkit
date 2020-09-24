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
#ifndef __RD_TAUTOMER_CATALOG_PARAMS_H__
#define __RD_TAUTOMER_CATALOG_PARAMS_H__

#include <Catalogs/CatalogParams.h>
#include "TautomerCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <string>
#include <vector>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class TautomerTransform;

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerCatalogParams
    : public RDCatalog::CatalogParams {
 public:
  TautomerCatalogParams() {
    d_typeStr = "Tautomer Catalog Parameters";
    d_transforms.clear();
  }

  TautomerCatalogParams(const std::string &tautomerFile);
  // copy constructor
  TautomerCatalogParams(const TautomerCatalogParams &other);

  ~TautomerCatalogParams() override;

  // DEPRECATED: remove in release 2021.01
  // the function name is misleading and getTransforms().size()
  // yields the same information
  unsigned int getNumTautomers() const {
    return static_cast<unsigned int>(d_transforms.size());
  }

  const std::vector<TautomerTransform> &getTransforms() const;

  const TautomerTransform getTransform(unsigned int fid) const;

  void toStream(std::ostream &) const override;
  std::string Serialize() const override;
  void initFromStream(std::istream &ss) override;
  void initFromString(const std::string &text) override;

 private:
  //		std::vector<std::pair<ROMol*, ROMol*>> d_pairs;
  std::vector<TautomerTransform> d_transforms;

};  // class TautomerCatalogParams

}  // namespace MolStandardize
}  // namespace RDKit

#endif
