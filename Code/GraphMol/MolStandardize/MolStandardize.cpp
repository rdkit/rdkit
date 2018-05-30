#include "MolStandardize.h"
#include "Metal.h"
#include <GraphMol/RDKitBase.h>
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>

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
        // normalize(*newM)
	// reionize(*newM)
	RDKit::MolOps::assignStereochemistry(*newM);

	return passedOp;
}

void tautomerParent(RWMol &mol, const CleanupParameters &params){
}

void stereoParent(RWMol &mol, const CleanupParameters &params){
}

void isotopeParent(RWMol &mol, const CleanupParameters &params){
}

void chargeParent(RWMol &mol, const CleanupParameters &params){
}

void superParent(RWMol &mol, const CleanupParameters &params){
}


} // end of namespace MolStandardize
} // end of namespace RDKit

