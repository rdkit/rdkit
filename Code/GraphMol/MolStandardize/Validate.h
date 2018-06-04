#ifndef __RD_VALIDATE_H__
#define __RD_VALIDATE_H__

#include <iostream>

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

class Validator{
	public:
		Validator();
		Validator(const Validator &other);
		~Validator();
		
		std::ostringstream validate(const ROMol &mol, const ValidationType &vtype);

}; // class Validator
} // namespace MolStandardize
} // namespace RDKit

#endif
