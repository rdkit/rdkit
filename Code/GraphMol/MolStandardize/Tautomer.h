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
const std::string tautomerScoringVersion = "1.0.0";

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

  //! returns all tautomers for the input molecule
  std::vector<ROMOL_SPTR> enumerate(const ROMol &mol) const;

  //! returns the canonical tautomer from a set of possible tautomers
  /*!
    Note that the canonical tautomer is very likely not the most stable tautomer
    for any given conditions. The default scoring rules are designed to produce
    "reasonable" tautomers, but the primary concern is that the results are
    canonical: you always get the same canonical tautomer for a molecule
    regardless of what the input tautomer or atom ordering were.
  */
  ROMol *pickCanonical(const std::vector<ROMOL_SPTR> &tautomers,
                       boost::function<int(const ROMol &mol)> scoreFunc =
                           TautomerScoringFunctions::scoreTautomer) const;

  //! returns the canonical tautomer for a molecule
  /*!
    Note that the canonical tautomer is very likely not the most stable tautomer
    for any given conditions. The default scoring rules are designed to produce
    "reasonable" tautomers, but the primary concern is that the results are
    canonical: you always get the same canonical tautomer for a molecule
    regardless of what the input tautomer or atom ordering were.
  */
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
