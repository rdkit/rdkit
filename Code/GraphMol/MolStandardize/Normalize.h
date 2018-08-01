#ifndef __RD_NORMALIZE_H__
#define __RD_NORMALIZE_H__

#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogEntry.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>
#include <GraphMol/MolStandardize/MolStandardize.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {
extern const CleanupParameters defaultCleanupParameters;

typedef RDCatalog::HierarchCatalog<TransformCatalogEntry,
                                   TransformCatalogParams, int>
    TransformCatalog;

class Normalizer {
 public:
	Normalizer();
	Normalizer(const std::string normalizeFile, const unsigned int maxRestarts);
	Normalizer(const Normalizer &other);
	~Normalizer();

  ROMol *normalize(const ROMol &mol);
  struct Product {
    std::string Smiles;
    boost::shared_ptr<ROMol> Mol;
    Product(std::string smiles, boost::shared_ptr<ROMol> &mol)
        : Smiles(smiles), Mol(mol) {}

    // sorting products alphabetically by SMILES
    bool operator<(const Product &pdt) const { return (Smiles < pdt.Smiles); }
  };

 private:
	TransformCatalog *d_tcat;
	unsigned int MAX_RESTARTS;

  ROMol *normalizeFragment(
      const ROMol &mol,
      const std::vector<std::shared_ptr<ChemicalReaction>> &transforms);
  boost::shared_ptr<ROMol> applyTransform(const ROMol &mol,
                                          ChemicalReaction &rule);

};  // Normalizer class
}  // namespace MolStandardize
}  // namespace RDKit

#endif
