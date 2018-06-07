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

class ValidationErrorInfo: public std::exception {
	public:
		ValidationErrorInfo(const std::string &msg): _msg(msg){};
		const char* message() const { return _msg.c_str(); };
		~ValidationErrorInfo() throw() {};
	private:
		std::string _msg;
}; // class ValidationErrorInfo

class ValidationMethod{
	public:
		virtual std::vector<ValidationErrorInfo> validate(const ROMol &mol, bool reportAllFailures) = 0;
};

class RDKitValidation : public ValidationMethod {
	public:
		std::vector<ValidationErrorInfo> validate(const ROMol &mol, bool reportAllFailures) override;
};

} // namespace MolStandardize
} // namespace RDKit

#endif
