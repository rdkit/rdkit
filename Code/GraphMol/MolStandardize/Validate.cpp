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

Validator::Validator() {};

Validator::Validator(const Validator &other) {};

Validator::~Validator() {};

std::vector<ValidationErrorInfo> Validator::validate(const ROMol &mol, const ValidationType &vtype, bool reportAllFailures) {
	std::vector<ValidationErrorInfo> errors;

	switch (vtype) {
		case MolStandardize::RDKitDefault:
			errors  = Validator::validateRDKitDefault(mol);
			break;
//		case MolStandardize::AllowedAtoms:
//			errout << "AllowedAtoms mode";
//                        ValidationErrorInfo res = Validator::validateAllowedAtoms(const ROMol &mol);
//			break;
//		case MolStandardize::DisallowedAtoms:
//			errout << "DisallowedAtoms mode";
//                        ValidationErrorInfo res = Validator::validateDisallowedAtoms(const ROMol &mol);
//			break;
	}
	return errors;
}

std::vector<ValidationErrorInfo> Validator::validateRDKitDefault(const ROMol &mol) {
	
	ROMol molCopy = mol;
	std::vector<ValidationErrorInfo> errors;

	unsigned int na = mol.getNumAtoms();
	
	if (!na){
		errors.push_back(ValidationErrorInfo("Molecule has no atoms"));
	}
	
	// loop over atoms
	for (size_t i = 0; i < na; ++i) {
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

//ValidationErrorInfo Validator::validateAllowedAtoms(const ROMol &mol) {
//        return ValidationErrorInfo("RDKitDefault mode");
//}
//
//ValidationErrorInfo Validator::validateDisallowedAtoms(const ROMol &mol) {
//        return ValidationErrorInfo("RDKitDefault mode");
//}

} // namespace MolStandardize
} // namespace RDKit
