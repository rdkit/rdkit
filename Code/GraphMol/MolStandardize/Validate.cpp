#include "Validate.h"
#include <GraphMol/RDKitBase.h>
#include <iostream>

using namespace std;
using namespace RDKit;

namespace RDKit{
class RWMol;
class ROMol;

namespace MolStandardize{

Validator::Validator() {};

Validator::Validator(const Validator &other) {};

Validator::~Validator() {};

ostringstream Validator::validate(const ROMol &mol, const ValidationType &vtype) {
	ostringstream errout;

	switch (vtype) {
		case MolStandardize::RDKitDefault:
			errout << "RDKitDefault mode";
			break;
		case MolStandardize::AllowedAtoms:
			errout << "AllowedAtoms mode";
			break;
		case MolStandardize::DisallowedAtoms:
			errout << "DisallowedAtoms mode";
			break;
	}
	return errout;
}

} // namespace MolStandardize
} // namespace RDKit
