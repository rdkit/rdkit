#include "Normalize.h"
#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

using namespace std;
using namespace RDKit;

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

int MAX_RESTARTS = 200;

ROMol* Normalizer::normalize(const ROMol &mol, TransformCatalog *tcat) {

	PRECONDITION(tcat, "");
	const TransformCatalogParams *tparams = tcat->getCatalogParams();

	PRECONDITION(tparams, "");
	const std::vector<std::shared_ptr<ChemicalReaction>> &transforms = tparams->getTransformations();

	std::vector<boost::shared_ptr<ROMol>> frags =	MolOps::getMolFrags(mol);
	std::vector<ROMol*> nfrags;//( frags.size() );
	for (const auto &frag : frags) {	
		ROMol* nfrag = this->normalizeFragment(*frag, transforms);
		nfrags.push_back(nfrag);
	}
	ROMol* outmol = nfrags.back();
	nfrags.pop_back();
	for (const auto &nfrag : nfrags) {
		outmol = combineMols(*outmol, *nfrag);
	}
	return outmol;
}

ROMol* Normalizer::normalizeFragment(const ROMol &mol, 
								const std::vector<std::shared_ptr<ChemicalReaction>> &transforms) {

	ROMol *nfrag = new ROMol(mol);
	bool loop_brake = false;
	for (int i = 0; i < MAX_RESTARTS; ++i) {
		// Iterate through Normalization transforms and apply each in order
		for (auto &transform : transforms) {
			std::string tname;
			transform->getProp(common_properties::_Name, tname);
			boost::shared_ptr<ROMol> product = this->applyTransform(mol, *transform);
			if (product != nullptr) {
//				std::cout << "Rule applied: " << tname << std::endl;
				nfrag = new ROMol(*product);
				loop_brake = true;
				break;
			}			
		}
		// For loop finishes normally, all applicable transforms have been applied
		if (loop_brake == false) { return nfrag; }
	}
	std::cout << "Gave up normalization after " << MAX_RESTARTS << " restarts." << 
				 std::endl;	
	return nfrag;
}

boost::shared_ptr<ROMol> Normalizer::applyTransform(const ROMol &mol, ChemicalReaction &transform) {
	// Repeatedly apply normalization transform to molecule until no changes occur.
	// 
	// It is possible for multiple products to be produced when a rule is applied. 
	// The rule is applied repeatedly to each of the products, until no further 
	// changes occur or after 20 attempts.
	//
	// If there are multiple unique products after the final application, the 
	// first product (sorted alphabetically by SMILES) is chosen.
	
	boost::shared_ptr<ROMol> tmp( new ROMol(mol) );
	MOL_SPTR_VECT mols;
	mols.push_back(tmp);

	transform.initReactantMatchers();
	for (int i = 0; i < 20; ++i) {
		std::vector<Normalizer::Product> pdts;
		for (auto &m : mols) {
			std::vector<MOL_SPTR_VECT> products = transform.runReactants( {m} );	
			for (const auto &pdt : products) {
				shared_ptr<RWMol> p0( new RWMol(*pdt[0]) );
//				std::cout << MolToSmiles(*p0) << std::endl;
				unsigned int failed;
				MolOps::sanitizeMol(*p0, failed);
				if (failed == 0) {
					boost::shared_ptr<ROMol> p0_ROMol( new ROMol(*p0) );
					Normalizer::Product np(MolToSmiles(*p0), p0_ROMol);
					pdts.push_back(np);
				}	else { std::cout << "FAILED sanitizeMol." << std::endl; }	
			}
		}
		if (pdts.size() != 0) {
			std::sort(pdts.begin(), pdts.end());
			mols.clear();
			for (const auto pdt : pdts) {
				mols.push_back(pdt.Mol);
			}
		} else {
			if (i > 0) {
				return mols[0];
			} else {
				return nullptr;
			}
		}
	}
}

} // namespace MolStandardize
} // namespace RDKit

