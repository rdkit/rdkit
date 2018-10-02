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
#ifndef __RD_TAUTOMER_H__
#define __RD_TAUTOMER_H__

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

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerCanonicalizer {
 public:
  //	TautomerCanonicalizer(unsigned int max_tautomers)
  //			: MAX_TAUTOMERS(max_tautomers) {};
  //	TautomerCanonicalizer(const TautomerCanonicalizer &other) {
  //		MAX_TAUTOMERS = other.MAX_TAUTOMERS;
  //	};
  //	~TautomerCanonicalizer() {};

  ROMol *canonicalize(const ROMol &mol, TautomerCatalog *tautcat);

  //	private:
  //		unsigned int MAX_TAUTOMERS;
};  // TautomerCanonicalizer class

class RDKIT_MOLSTANDARDIZE_EXPORT TautomerEnumerator {
 public:
  std::vector<ROMOL_SPTR> enumerate(const ROMol &mol, TautomerCatalog *tautcat);

  //		struct Tautomer {
  //			std::string Smiles;
  //			boost::shared_ptr<ROMol> Mol;
  //			Tautomer(std::string smiles, boost::shared_ptr<ROMol>
  // mol) 				: Smiles(smiles), Mol(mol) {}
  //
  //			// sorting products alphabetically by SMILES
  //			bool operator < (const Tautomer &tautomer) const {
  //				return (Smiles < tautomer.Smiles);
  //			}
  //
  //		};

  //		TautomerEnumerator(unsigned int max_tautomers)
  //			: MAX_TAUTOMERS(max_tautomers) {};
  //		TautomerEnumerator(const TautomerEnumerator &other) {
  //			MAX_TAUTOMERS = other.MAX_TAUTOMERS;
  //		};
  //		~TautomerEnumerator() {};
  //	private:
  //		unsigned int MAX_TAUTOMERS;
};  // TautomerEnumerator class

std::vector<std::pair<unsigned int, unsigned int>> pairwise(
    const std::vector<int> vect);
}  // namespace MolStandardize
}  // namespace RDKit

#endif
