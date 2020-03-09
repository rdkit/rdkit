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
#ifndef __RD_FRAGMENT_CATALOG_PARAMS_H__
#define __RD_FRAGMENT_CATALOG_PARAMS_H__

#include <Catalogs/CatalogParams.h>
#include "FragmentCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <string>
#include <vector>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class RDKIT_MOLSTANDARDIZE_EXPORT FragmentCatalogParams
    : public RDCatalog::CatalogParams {
 public:
  FragmentCatalogParams() {
    d_typeStr = "Fragment Catalog Parameters";
    d_funcGroups.clear();
  }

  FragmentCatalogParams(const std::string &fgroupFile);
  FragmentCatalogParams(std::istream &fgroupStream);
  // copy constructor
  FragmentCatalogParams(const FragmentCatalogParams &other);

  ~FragmentCatalogParams() override;

  unsigned int getNumFuncGroups() const {
    return static_cast<unsigned int>(d_funcGroups.size());
  }

  const std::vector<std::shared_ptr<ROMol>> &getFuncGroups() const;

  const ROMol *getFuncGroup(unsigned int fid) const;

  void toStream(std::ostream &) const override;
  std::string Serialize() const override;
  void initFromStream(std::istream &ss) override;
  void initFromString(const std::string &text) override;

 private:
  std::vector<std::shared_ptr<ROMol>> d_funcGroups;

};  // class FragmentCatalogParams

}  // namespace MolStandardize
}  // namespace RDKit

#endif
