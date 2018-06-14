#ifndef __RD_FRAGMENT_CATALOG_ENTRY_H__
#define __RD_FRAGMENT_CATALOG_ENTRY_H__

#include <Catalogs/CatalogEntry.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "FragmentCalaogParams.h"
#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace MolStandardize {

class FragmentCatalogEntry : public RDCatalog::CatalogEntry {
	public:
		FragmentCatalogEntry() :
			dp_mol(0), d_descrip("") {
				dp_props = new Dict();
				setBitId(-1);
			}

		FragmentCatalogEntry(const ROMol *omol, const PATH_TYPE &path);
		FragmentCatalogEntry(const std::string &pickle);

		~FragmentCatalogEntry() override {
			delete dp_mol;
			dp_mol = 0;
			if (dp_props) {
				delete dp_props;
				dp_props = 0;
			}
		}

		std::string getDescription() const override { return d_descrip; }

		unsigned int getOrder() const { return dp_mol->getNumBonds(); }

		void toStream(std::ostream &ss) const override;
		std::string Serialize() const override;
		void initFromStream(std::istream &ss) override;
		void initFromString(const std::string &text) override;
	private:
		ROMol *dp_mol;
		Dict *dp_props;
		std::string d_descrip;

}; // class FragmentCatalogEntry

} // namespace MolStandardize 
} // namespace RDKit

#endif
