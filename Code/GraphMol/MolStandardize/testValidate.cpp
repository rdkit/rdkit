#include "Validate.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <iostream>

using namespace RDKit;
using namespace std;

void testValidate(){
	string smi1;
	MolStandardize::Validator v;

	smi1 = "ClCCCl.c1ccccc1O";
	unique_ptr<ROMol> m1( SmilesToMol(smi1) );
	std::ostringstream errout = v.validate(*m1, MolStandardize::RDKitDefault);
	std::string msg = errout.str();
	TEST_ASSERT(msg == "RDKitDefault mode");
	std::cout << msg << std::endl;

	std::ostringstream errout2 = v.validate(*m1, MolStandardize::AllowedAtoms);    
        std::string msg2 = errout2.str();                                               
        TEST_ASSERT(msg2 == "AllowedAtoms mode");
	std::cout << msg2 << std::endl;

	std::ostringstream errout3 = v.validate(*m1, MolStandardize::DisallowedAtoms);
	std::string msg3 = errout3.str();
	TEST_ASSERT(msg3 == "DisallowedAtoms mode");
	std::cout << msg3 << std::endl;
}

int main() {
	testValidate();
	return 0;
}
