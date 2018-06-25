#ifndef __RD_NORMALIZE_H__
#define __RD_NORMALIZE_H__

#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogEntry.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

typedef RDCatalog::HierarchCatalog<TransformCatalogEntry, TransformCatalogParams, int>
    TransformCatalog;

class Normalizer {
	public:
		ROMol* normalize(const ROMol &mol, TransformCatalog *tcat);
		struct Product {
			std::string Smiles;
			boost::shared_ptr<ROMol> Mol;
			Product(std::string smiles, boost::shared_ptr<ROMol> &mol): Smiles(smiles), Mol(mol) {}

			// sorting products alphabetically by SMILES
			bool operator < (const Product &pdt) const {
							return (Smiles < pdt.Smiles);
			}
		};
	private:
		ROMol* normalizeFragment(const ROMol &mol, const std::vector<std::shared_ptr<ChemicalReaction>> &transforms);
		boost::shared_ptr<ROMol> applyTransform(const ROMol &mol, ChemicalReaction &rule);

}; // Normalizer class
}
}

#endif
