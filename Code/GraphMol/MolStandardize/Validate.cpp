#include "Validate.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace RDKit;

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

std::vector<ValidationErrorInfo> RDKitValidation::validate(const ROMol &mol, bool reportAllFailures) {
	
	ROMol molCopy = mol;
	std::vector<ValidationErrorInfo> errors;

	unsigned int na = mol.getNumAtoms();
	
	if (!na){
		errors.push_back(ValidationErrorInfo("Molecule has no atoms"));
	}
	
	// loop over atoms
	for (size_t i = 0; i < na; ++i) {
		if (!reportAllFailures) {
			if (errors.size() >= 1) {break;}
		}
		Atom* atom = molCopy.getAtomWithIdx(i);
		try{
			int explicitValence = atom->calcExplicitValence();
		} catch (const MolSanitizeException &e) {
			//std::cout << e.message() << std::endl;
			errors.push_back(ValidationErrorInfo(e.message()));
		}
	}
	return errors;
}

} // namespace MolStandardize
} // namespace RDKit
