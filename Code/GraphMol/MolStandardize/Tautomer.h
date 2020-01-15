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

namespace TautomerScoringFunctions{
  RDKIT_MOLSTANDARDIZE_EXPORT int scoreRings(const ROMol &mol);
  RDKIT_MOLSTANDARDIZE_EXPORT int scoreSubstructs(const ROMol &mol);
  RDKIT_MOLSTANDARDIZE_EXPORT int scoreHeteroHs(const ROMol &mol);

  inline int scoreTautomer(const ROMol &mol){
    return scoreRings(mol)+scoreSubstructs(mol)+scoreHeteroHs(mol);
  }
}


class RDKIT_MOLSTANDARDIZE_EXPORT TautomerCanonicalizer {
 public:
  ROMol *canonicalize(const ROMol &mol, TautomerCatalog *tautcat);
};  // TautomerCanonicalizer class

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumerator {
 public:
  std::vector<ROMOL_SPTR> enumerate(const ROMol &mol, TautomerCatalog *tautcat);
};  // TautomerEnumerator class

}  // namespace MolStandardize
}  // namespace RDKit

#endif
