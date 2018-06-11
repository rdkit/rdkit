#include "Validate.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace RDKit;

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

std::vector<ValidationErrorInfo> RDKitValidation::validate(const ROMol &mol, bool reportAllFailures) const {
	
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


std::vector<ValidationErrorInfo> MolVSValidation::validate(const ROMol &mol, bool reportAllFailures) const {

	std::vector<ValidationErrorInfo> errors;

	this->noAtomValidation(mol, reportAllFailures, errors);
	this->fragmentValidation(mol, reportAllFailures, errors); //TODO filterscatalog stuff
	this->neutralValidation(mol, reportAllFailures, errors);
	this->isotopeValidation(mol, reportAllFailures, errors);

	return errors;
}

void MolVSValidation::noAtomValidation(const ROMol &mol, bool reportAllFailures, std::vector<ValidationErrorInfo> &errors) const {
	unsigned int na = mol.getNumAtoms();
	
	if (!na){
		errors.push_back(ValidationErrorInfo("Molecule has no atoms"));
	}	
}

void MolVSValidation::fragmentValidation(const ROMol &mol, bool reportAllFailures, std::vector<ValidationErrorInfo> &errors) const {

}

void MolVSValidation::neutralValidation(const ROMol &mol, bool reportAllFailures, std::vector<ValidationErrorInfo> &errors) const {
	int charge = RDKit::MolOps::getFormalCharge(mol);
	if (charge != 0) {
		std::string msg = "Not an overall neutral system (" + std::to_string(charge) + ')';
//		std::cout << msg << std::endl;
		errors.push_back(ValidationErrorInfo(msg));
	}
}

void MolVSValidation::isotopeValidation(const ROMol &mol, bool reportAllFailures, std::vector<ValidationErrorInfo> &errors) const {

	unsigned int na = mol.getNumAtoms();
	std::set<string> isotopes;
	
	// loop over atoms
	for (size_t i = 0; i < na; ++i) {
		if (!reportAllFailures) {
			if (errors.size() >= 1) {break;}
		}
		const Atom* atom = mol.getAtomWithIdx(i) ;
		unsigned int isotope = atom->getIsotope();
		if (isotope != 0) {
			std::string symbol = atom->getSymbol();
			isotopes.insert( std::to_string(isotope) + symbol );
		}
	}

	for (auto &isotope : isotopes) {
		errors.push_back(ValidationErrorInfo("Molecule contains isotope " + isotope));
	}
}

std::vector<ValidationErrorInfo> AllowedAtomsValidation::validate(const ROMol &mol, bool reportAllFailures) const {

	std::vector<ValidationErrorInfo> errors;
	unsigned int na = mol.getNumAtoms();

	for (size_t i = 0; i < na; ++i) {
		if (!reportAllFailures) {
			if (errors.size() >= 1) {break;}
		}
		const Atom* qatom = mol.getAtomWithIdx(i);
		bool match = false;
		// checks to see qatom matches one of list of allowedAtoms
		for (const auto &allowedAtom : this->d_allowedList) {
			if ( allowedAtom->Match(qatom) ) {
			       match = true;
			       break;
			}		
		}
		// if no match, append to list of errors.
		if (! match) {
			std::string symbol = qatom->getSymbol();
			errors.push_back(ValidationErrorInfo("Atom " + symbol + " is not in allowedAtoms list"));
		}
	}
	return errors;
}

std::vector<ValidationErrorInfo> DisallowedAtomsValidation::validate(const ROMol &mol, bool reportAllFailures) const {

	std::vector<ValidationErrorInfo> errors;
	unsigned int na = mol.getNumAtoms();

	for (size_t i = 0; i < na; ++i) {
		if (!reportAllFailures) {
			if (errors.size() >= 1) {break;}
		}
		const Atom* qatom = mol.getAtomWithIdx(i);
		bool match = false;
		// checks to see qatom matches one of list of allowedAtoms
		for (const auto &disallowedAtom : this->d_disallowedList) {
			if ( disallowedAtom->Match(qatom) ) {
			       match = true;
			}		
		}
		// if no match, append to list of errors.
		if (match) {
			std::string symbol = qatom->getSymbol();
			errors.push_back(ValidationErrorInfo("Atom " + symbol + " is in disallowedAtoms list"));
		}
	}
	return errors;
}

} // namespace MolStandardize
} // namespace RDKit
