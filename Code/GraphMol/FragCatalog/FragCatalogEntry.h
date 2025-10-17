//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_FRAGCATALOGENTRY_H_
#define _RD_FRAGCATALOGENTRY_H_

#include "FragCatParams.h"
#include <RDGeneral/utils.h>
#include <Catalogs/CatalogEntry.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <map>
#include <sstream>
#include <string_view>

namespace RDKit {

class RDKIT_FRAGCATALOG_EXPORT FragCatalogEntry
    : public RDCatalog::CatalogEntry {
 public:
  FragCatalogEntry() : d_descrip("") {
    dp_props = new Dict();
    setBitId(-1);
  }

  FragCatalogEntry(const ROMol *omol, const PATH_TYPE &path,
                   const MatchVectType &aidToFid);
  FragCatalogEntry(const std::string &pickle);

  ~FragCatalogEntry() override {
    delete dp_mol;
    dp_mol = nullptr;
    if (dp_props) {
      delete dp_props;
      dp_props = nullptr;
    }
  }

  std::string getDescription() const override { return d_descrip; }

  void setDescription(const std::string &val) { d_descrip = val; }

  void setDescription(const FragCatParams *params);

  // check if this fragment matches the one specified
  //

  bool match(const FragCatalogEntry *other, double tol) const;

  Subgraphs::DiscrimTuple getDiscrims() const;

  unsigned int getOrder() const { return dp_mol->getNumBonds(); }

  const INT_INT_VECT_MAP &getFuncGroupMap() const { return d_aToFmap; }

  // REVIEW: this should be removed?
  std::string getSmarts() { return ""; }

  // FUnctions on the property dictionary
  template <typename T>
  void setProp(const std::string_view &key, T &val) const {
    dp_props->setVal(key, val);
  }

  void setProp(const std::string_view &key, int val) const {
    dp_props->setVal(key, val);
  }

  void setProp(const std::string_view &key, float val) const {
    dp_props->setVal(key, val);
  }

  void setProp(const std::string_view &key, std::string &val) const {
    dp_props->setVal(key, val);
  }

  template <typename T>
  void getProp(const std::string_view key, T &res) const {
    getProp(std::string(key).c_str(), res);
  }

  bool hasProp(const std::string_view key) const {
    return hasProp(std::string(key).c_str());
  }

  void clearProp(const std::string_view key) const {
    clearProp(std::string(key).c_str());
  }

  void toStream(std::ostream &ss) const override;
  std::string Serialize() const override;
  void initFromStream(std::istream &ss) override;
  void initFromString(const std::string &text) override;

 private:
  ROMol *dp_mol{nullptr};
  Dict *dp_props;

  std::string d_descrip;

  unsigned int d_order{0};

  // a map between the atom ids in mol that connect to
  // a functional group and the corresponding functional
  // group ID
  INT_INT_VECT_MAP d_aToFmap;
};
}  // namespace RDKit

#endif
