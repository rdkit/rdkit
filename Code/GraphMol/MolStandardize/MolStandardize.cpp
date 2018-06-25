#include "MolStandardize.h"
#include "Metal.h"
#include "Normalize.h"
#include <GraphMol/RDKitBase.h>
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentRemover.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
using namespace std;

namespace RDKit{
namespace MolStandardize{

bool cleanup(RWMol &mol, const CleanupParameters &params){	
	bool passedOp = false;

	auto *newM = new RWMol(mol);
	RDKit::MolOps::sanitizeMol(*newM);	
	RDKit::MolOps::removeHs(*newM);

	// TODO
	MolStandardize::MetalDisconnector md;
	md.disconnect(mol);
	MolStandardize::normalize(mol, params);
	// reionize(*newM)
	RDKit::MolOps::assignStereochemistry(*newM);

	return passedOp;
}

void tautomerParent(RWMol &mol, const CleanupParameters &params){
}

// Return the fragment parent of a given molecule.
// The fragment parent is the largest organic covalent unit in the molecule.
//
void fragmentParent(RWMol &mol, const CleanupParameters &params, bool skip_standardize) {
	
	if (!skip_standardize) {
		cleanup(mol, params);
	}

	// TODO
	// largest fragment 
	LargestFragmentChooser lfragchooser(params.preferOrganic);
	std::shared_ptr<ROMol> nm( new ROMol(mol) );
	boost::shared_ptr<ROMol> lfrag = lfragchooser.choose(*nm);
	mol = RWMol(*lfrag);
}

void stereoParent(RWMol &mol, const CleanupParameters &params){
}

void isotopeParent(RWMol &mol, const CleanupParameters &params){
}

void chargeParent(RWMol &mol, const CleanupParameters &params){
}

void superParent(RWMol &mol, const CleanupParameters &params){
}

void normalize(RWMol &mol, const CleanupParameters &params){
	
	auto *tparams = new TransformCatalogParams(params.normalizations);
	TransformCatalog tcat(tparams);
	Normalizer normalizer;

	std::shared_ptr<ROMol> m( new ROMol(mol) );
	ROMol* normalized = normalizer.normalize(*m, &tcat);

	mol = RWMol(*normalized);
}

} // end of namespace MolStandardize
} // end of namespace RDKit

