#include "FragmentRemover.h"
#include "FragmentCalaogUtils.h"
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

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
}

bool isOrganic(const ROMol &frag) {
	// Returns true if fragment contains at least one carbon atom.
	ROMol::ConstAtomIterator ai;
	for(ai = frag.beginAtoms(); ai != frag.endAtoms() ; ++ai) {
    if ( (*ai)->getAtomicNum() == 6 ) {
					return true;
				}
	}
	return false;
}

LargestFragmentChooser::LargestFragmentChooser(const LargestFragmentChooser &other) {
	PREFER_ORGANIC = other.PREFER_ORGANIC;
}

boost::shared_ptr<ROMol> LargestFragmentChooser::choose(const ROMol &mol) {
//	auto m = new ROMol(mol);

	std::vector<boost::shared_ptr<ROMol>> frags = 
		MolOps::getMolFrags(mol);
	LargestFragmentChooser::Largest l;

	for (const auto &frag : frags) {
		std::string smiles = MolToSmiles(*frag);
		std::cout << "Fragment: " << smiles << std::endl;
		bool organic = isOrganic(*frag);
  	if ( this->PREFER_ORGANIC ) {
  		// Skip this fragment if not organic and we already have an organic fragment as the largest so far
  		if (l.Fragment != NULL && l.Organic && ! organic ) continue;
			// Reset largest if it wasn't organic and this fragment is organic
			// if largest and organic and not largest['organic']:
			if (l.Fragment != NULL && organic && ! l.Organic) {
				l.Fragment = nullptr;	
			}
		}
		ROMol::AtomIterator ai;
		unsigned int atoms = 0;
		for (ai = (*frag).beginAtoms(); ai != (*frag).endAtoms() ; ++ai) {
						atoms += 1 + (*ai)->getTotalNumHs();
		}
		// Skip this fragment if fewer atoms than the largest
		if ( l.Fragment != NULL && (atoms < l.Atoms) ) continue;

		// Skip this fragment if equal number of atoms but weight is lower
		double weight = Descriptors::calcExactMW(*frag);
		if ( l.Fragment != NULL && (atoms == l.Atoms) && (weight < l.Weight) ) continue;
			 
 		// Skip this fragment if equal atoms and equal weight but smiles 
		// comes last alphabetically
		if ( l.Fragment != NULL && (atoms == l.Atoms) && (weight == l.Weight) &&
						(smiles > l.Smiles) ) continue;

		//Otherwise this is the largest so far
		l.Smiles = smiles;
		l.Fragment = frag;
		l.Atoms = atoms;
		l.Weight = weight;
		l.Organic = organic;
	}
	
	return l.Fragment;
//	return m;

}

LargestFragmentChooser::Largest::Largest() 
	: Smiles(""), Fragment(nullptr), Atoms(0), Weight(0), Organic(false) {}

LargestFragmentChooser::Largest::Largest(std::string &smiles,
							 	const boost::shared_ptr<ROMol> &fragment, 
								unsigned int &atoms, double &weight, bool &organic)
					: Smiles(smiles), Fragment(fragment), Atoms(atoms),
Weight(weight), Organic(organic) {}

} // namespace MolStandardize
} // namespace RDKit
