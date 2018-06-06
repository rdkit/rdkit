#ifndef __RD_VALIDATE_H__
#define __RD_VALIDATE_H__

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <iostream>
#include <exception>
#include <string>
#include <vector>

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

enum ValidationType {
	Unspecified = 0,
	RDKitDefault, 
	AllowedAtoms,
	DisallowedAtoms
};

//struct ValidationParams{
//	ValidationParams()
//	       : valType(Unspecified) {}	
//	ValidationParams(ValidationType valType)
//		: valType(valType) {}
//	ValidationType valType;
//};

class ValidationErrorInfo: public std::exception {
	public:
		ValidationErrorInfo(const std::string &msg): _msg(msg){};
		const char* message() const { return _msg.c_str(); };
		~ValidationErrorInfo() throw() {};
	private:
		std::string _msg;
}; // class ValidationErrorInfo

class Validator{
	public:
		Validator();
		Validator(const Validator &other);
		~Validator();
		
		std::vector<ValidationErrorInfo> validate(const ROMol &mol, const ValidationType &vtype, bool reportAllFailures = true);

		std::vector<ValidationErrorInfo> validateRDKitDefault(const ROMol &mol);
//		ValidationErrorInfo validateAllowedAtoms(const ROMol &mol); 
//		ValidationErrorInfo validateDisallowedAtoms(const ROMol &mol); 
}; // class Validator
} // namespace MolStandardize
} // namespace RDKit

#endif
