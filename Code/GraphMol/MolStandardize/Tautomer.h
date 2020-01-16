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
#ifndef RD_TAUTOMER_H
#define RD_TAUTOMER_H

#include <boost/function.hpp>
#include <string>
#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogEntry.h>
#include <GraphMol/MolStandardize/TautomerCatalog/TautomerCatalogParams.h>

namespace RDKit {
class ROMol;
class RWMol;

namespace MolStandardize {

typedef RDCatalog::HierarchCatalog<TautomerCatalogEntry, TautomerCatalogParams,
                                   int>
    TautomerCatalog;

namespace TautomerScoringFunctions {
RDKIT_MOLSTANDARDIZE_EXPORT int scoreRings(const ROMol &mol);
RDKIT_MOLSTANDARDIZE_EXPORT int scoreSubstructs(const ROMol &mol);
RDKIT_MOLSTANDARDIZE_EXPORT int scoreHeteroHs(const ROMol &mol);

inline int scoreTautomer(const ROMol &mol) {
  return scoreRings(mol) + scoreSubstructs(mol) + scoreHeteroHs(mol);
}
}  // namespace TautomerScoringFunctions

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumerator {
 public:
  TautomerEnumerator() = delete;
  TautomerEnumerator(TautomerCatalog *tautCat) : dp_catalog(tautCat){};
  TautomerEnumerator(const TautomerEnumerator &other)
      : dp_catalog(other.dp_catalog){};
  TautomerEnumerator &operator=(const TautomerEnumerator &other) {
    if (this == &other) return *this;
    dp_catalog = other.dp_catalog;
    return *this;
  }

  std::vector<ROMOL_SPTR> enumerate(const ROMol &mol) const;
  ROMol *pickCanonical(const std::vector<ROMOL_SPTR> &tautomers,
                       boost::function<int(const ROMol &mol)> scoreFunc =
                           TautomerScoringFunctions::scoreTautomer) const;
  ROMol *canonicalize(const ROMol &mol,
                      boost::function<int(const ROMol &mol)> scoreFunc =
                          TautomerScoringFunctions::scoreTautomer) const {
    auto tautomers = enumerate(mol);
    if (!tautomers.size()) {
      BOOST_LOG(rdWarningLog)
          << "no tautomers found, returning input molecule" << std::endl;
      return new ROMol(mol);
    }
    return pickCanonical(tautomers, scoreFunc);
  };

 private:
  std::shared_ptr<TautomerCatalog> dp_catalog;
};  // TautomerEnumerator class

}  // namespace MolStandardize
}  // namespace RDKit

#endif
