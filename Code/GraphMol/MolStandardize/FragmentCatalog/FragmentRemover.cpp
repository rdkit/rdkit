#include "FragmentRemover.h"
#include "FragmentCalaogUtils.h"
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace RDKit {
namespace MolStandardize {

FragmentRemover::FragmentRemover(const FragmentRemover &other) {
	LEAVE_LAST = other.LEAVE_LAST;
};

ROMol* FragmentRemover::remove(const ROMol &mol, FragmentCatalog *fcat) {

	PRECONDITION(fcat, "");
	const FragmentCatalogParams *fparams = fcat->getCatalogParams();

	PRECONDITION(fparams, "");

	const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
	std::vector<std::shared_ptr<ROMol>>::const_iterator fgci;
	auto *removed =  new ROMol(mol) ;

	for (fgci = fgrps.begin(); fgci != fgrps.end(); fgci++) {

		std::vector<boost::shared_ptr<ROMol>> frags =
		       	MolOps::getMolFrags(*removed);
		// If nothing is left or leave_last and only one fragment, end here
		if (removed->getNumAtoms() == 0 || ( this->LEAVE_LAST 
		  && frags.size() <= 1 ) ) { break; }

		const ROMol *fgrp = fgci->get();
		std::string fname;
		(*fgci)->getProp(common_properties::_Name, fname);
		ROMol *tmp = RDKit::deleteSubstructs(*removed, *fgrp, true);

		if (tmp->getNumAtoms() != removed->getNumAtoms()) {
			std::cout << "Removed fragment: " << fname << std::endl;
		}

		if ( this->LEAVE_LAST && tmp->getNumAtoms() == 0 ) {
			// All the remaining fragments match this pattern - leave them all
			break;
		}
		removed = tmp;
	}
	return removed;

//	ROMol *nm = removeFrags(mol, fparams);
//	return nm;
}

} // namespace MolStandardize
} // namespace RDKit
