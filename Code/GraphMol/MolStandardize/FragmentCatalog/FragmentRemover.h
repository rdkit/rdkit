#ifndef __RD_FRAGMENT_REMOVER_H__
#define __RD_FRAGMENT_REMOVER_H__

#include <Catalogs/Catalog.h>
#include "FragmentCalaogEntry.h"
#include "FragmentCalaogParams.h"

namespace RDKit {
class ROMol;

namespace MolStandardize {

typedef RDCatalog::HierarchCatalog<FragmentCatalogEntry, FragmentCatalogParams, int>
    FragmentCatalog;

class FragmentRemover {
	public:
		FragmentRemover() {};
		FragmentRemover(bool leave_last) :
			LEAVE_LAST(leave_last) {};
		FragmentRemover(const FragmentRemover &other);
		~FragmentRemover() {};

		ROMol* remove(const ROMol &mol, FragmentCatalog *fcat);
	private:
		// Setting leave_last to True will ensure at least one fragment
		//  is left in the molecule, even if it is matched by a 
		//  FragmentPattern
		bool LEAVE_LAST = true;

}; // class FragmentRemover
} // namespace MolStandardize
} // namespace RDKit

#endif
