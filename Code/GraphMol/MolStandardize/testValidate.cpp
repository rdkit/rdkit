#include "Validate.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <iostream>

using namespace RDKit;
using namespace std;
using namespace MolStandardize;

void testValidate(){
	string smi1, smi2, smi3;
	Validator v;

	// testing RDKitDefault
	smi1 = "CO(C)C";
	unique_ptr<ROMol> m1( SmilesToMol(smi1, 0, false) );
	vector<ValidationErrorInfo> errout1 = v.validate(*m1, MolStandardize::RDKitDefault);
	for (auto &query : errout1) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "Explicit valence for atom # 1 O, 3, is greater than permitted");
	}

	smi2 = "";
	unique_ptr<ROMol> m2( SmilesToMol(smi2, 0, false) );
	vector<ValidationErrorInfo> errout2 = v.validate(*m2, MolStandardize::RDKitDefault);
	for (auto &query : errout2) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "Molecule has no atoms");
	}

	smi3 = "CO(C)CCN(=O)=O";
	unique_ptr<ROMol> m3( SmilesToMol(smi3, 0, false) );
	vector<ValidationErrorInfo> errout3 = v.validate(*m3, MolStandardize::RDKitDefault);
	std::vector<string> msgs;
	std::vector<string> ans = {"Explicit valence for atom # 1 O, 3, is greater than permitted", 
			"Explicit valence for atom # 5 N, 5, is greater than permitted"};
	for (auto &query : errout3) {
		msgs.push_back(query.message());
	}
	TEST_ASSERT(msgs == ans);
	
//	MolStandardize::ValidationErrorInfo errout2 = v.validate(*m1, MolStandardize::AllowedAtoms);    
//        std::string msg2 = errout2.message();                                               
//        TEST_ASSERT(msg2 == "AllowedAtoms mode");
//	std::cout << msg2 << std::endl;
//
//	MolStandardize::ValidationErrorInfo errout3 = v.validate(*m1, MolStandardize::DisallowedAtoms);
//	std::string msg3 = errout3.message();
//	TEST_ASSERT(msg3 == "DisallowedAtoms mode");
//	std::cout << msg3 << std::endl;
}

int main() {
	testValidate();
	return 0;
}
