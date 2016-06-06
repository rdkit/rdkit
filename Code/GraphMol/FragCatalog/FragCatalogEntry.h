//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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

namespace RDKit {

class FragCatalogEntry : public RDCatalog::CatalogEntry, public RDKit::RDProps {
 public:
  FragCatalogEntry() : dp_mol(0), d_descrip(""), d_order(0) {
    setBitId(-1);
  }

  FragCatalogEntry(const ROMol *omol, const PATH_TYPE &path,
                   const MatchVectType &aidToFid);
  FragCatalogEntry(const std::string &pickle);

  ~FragCatalogEntry() {
    delete dp_mol;
    dp_mol = 0;
  }

  std::string getDescription() const { return d_descrip; }

  void setDescription(const std::string &val) { d_descrip = val; }

  void setDescription(const FragCatParams *params);

  // check if this fragment macthes the one specified
  //

  bool match(const FragCatalogEntry *other, double tol) const;

  Subgraphs::DiscrimTuple getDiscrims() const;

  unsigned int getOrder() const { return dp_mol->getNumBonds(); }

  const INT_INT_VECT_MAP &getFuncGroupMap() const { return d_aToFmap; }

  // REVIEW: this should be removed?
  std::string getSmarts() { return ""; }

  void toStream(std::ostream &ss) const;
  std::string Serialize() const;
  void initFromStream(std::istream &ss);
  void initFromString(const std::string &text);

 private:
  ROMol *dp_mol;
  std::string d_descrip;

  unsigned int d_order;

  // a map between the atom ids in mol that connect to
  // a functional group and the corresponding functional
  // group ID
  INT_INT_VECT_MAP d_aToFmap;
};
}

#endif
