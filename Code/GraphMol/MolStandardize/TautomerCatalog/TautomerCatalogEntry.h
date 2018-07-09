#ifndef __RD_TAUTOMER_CATALOG_ENTRY_H__
#define __RD_TAUTOMER_CATALOG_ENTRY_H__

#include <Catalogs/CatalogEntry.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include "TautomerCatalogParams.h"
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace MolStandardize {

class TautomerCatalogEntry : public RDCatalog::CatalogEntry {
	public:
		TautomerCatalogEntry(){
				dp_transform = nullptr;
			 	d_descrip = "";
				dp_props = new Dict();
				setBitId(-1);
		}

		~TautomerCatalogEntry() override {
			delete dp_transform;
			dp_transform = nullptr;
			delete dp_props;
			dp_props = nullptr;
		}

		// TODO Catalog.h requires a getOrder function
		unsigned int getOrder() const { return 0; }//dp_mol->getNumBonds(); }

		void toStream(std::ostream &ss) const override;
		std::string Serialize() const override;
		void initFromStream(std::istream &ss) override;
		void initFromString(const std::string &text) override;
	private:
		TautomerTransform* dp_transform;
		Dict *dp_props;
		std::string d_descrip;

}; // class TautomerCatalogEntry

} // namespace MolStandardize 
} // namespace RDKit

#endif
