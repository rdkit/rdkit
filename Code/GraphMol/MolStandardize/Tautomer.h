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
#include <boost/dynamic_bitset.hpp>

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
  /*!
    The enumeration rules are inspired by the publication:
    M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
    https://doi.org/10.1007/s10822-010-9346-4

    \param mol: the molecule to be enumerated
    \param modifiedAtoms: if provided this is used to return which atoms are
    modified during the tautomerization
    \param modifiedBonds: if provided this is used to return which bonds are
    modified during the tautomerization

    Note: the definitions used here are that the atoms modified during
    tautomerization are the atoms at the beginning and end of each tautomer
    transform (the H "donor" and H "acceptor" in the transform) and the bonds
    modified during transformation are any bonds whose order is changed during
    the tautomer transform (these are the bonds between the "donor" and the
    "acceptor")

  */
  std::vector<ROMOL_SPTR> enumerate(
      const ROMol &mol, boost::dynamic_bitset<> *modifiedAtoms = nullptr,
      boost::dynamic_bitset<> *modifiedBonds = nullptr) const;

  //! returns the canonical tautomer from a set of possible tautomers
  /*!
    Note that the canonical tautomer is very likely not the most stable tautomer
    for any given conditions. The default scoring rules are designed to produce
    "reasonable" tautomers, but the primary concern is that the results are
    canonical: you always get the same canonical tautomer for a molecule
    regardless of what the input tautomer or atom ordering were.

    The default scoring scheme is inspired by the publication:
    M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
    https://doi.org/10.1007/s10822-010-9346-4

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

    The default scoring scheme is inspired by the publication:
    M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
    https://doi.org/10.1007/s10822-010-9346-4

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
