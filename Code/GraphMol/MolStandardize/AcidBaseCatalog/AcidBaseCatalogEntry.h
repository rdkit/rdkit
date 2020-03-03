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
#ifndef __RD_ACIDBASE_CATALOG_ENTRY_H__
#define __RD_ACIDBASE_CATALOG_ENTRY_H__

#include <Catalogs/CatalogEntry.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include "AcidBaseCatalogParams.h"
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace MolStandardize {

class RDKIT_MOLSTANDARDIZE_EXPORT AcidBaseCatalogEntry
    : public RDCatalog::CatalogEntry {
 public:
  AcidBaseCatalogEntry() {
    d_descrip = "";
    dp_props = new Dict();
    setBitId(-1);
  }

  ~AcidBaseCatalogEntry() override {
    delete dp_props;
    dp_props = nullptr;
  }

  // TODO Catalog.h requires a getOrder function
  unsigned int getOrder() const { return 0; }  // dp_mol->getNumBonds(); }

  void toStream(std::ostream &ss) const override;
  std::string Serialize() const override;
  void initFromStream(std::istream &ss) override;
  void initFromString(const std::string &text) override;

 private:
  std::shared_ptr<std::pair<ROMOL_SPTR, ROMOL_SPTR>> dp_pair;
  Dict *dp_props;
  std::string d_descrip;

};  // class AcidBaseCatalogEntry

}  // namespace MolStandardize
}  // namespace RDKit

#endif
