#include "Validate.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <iostream>

using namespace RDKit;
using namespace std;
using namespace MolStandardize;

void testRDKitValidation(){
	string smi1, smi2, smi3, smi4;
	RDKitValidation vm;

	// testing RDKitDefault
	smi1 = "CO(C)C";
	unique_ptr<ROMol> m1( SmilesToMol(smi1, 0, false) );
	vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
	for (auto &query : errout1) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "Explicit valence for atom # 1 O, 3, is greater than permitted");
	}

	// testing for molecule with no atoms
	smi2 = "";
	unique_ptr<ROMol> m2( SmilesToMol(smi2, 0, false) );
	vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
	for (auto &query : errout2) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "Molecule has no atoms");
	}

	// testing molecule with multiple valency errors
	smi3 = "CO(C)CCN(=O)=O";
	unique_ptr<ROMol> m3( SmilesToMol(smi3, 0, false) );
	vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
	std::vector<string> msgs1;
	std::vector<string> ans1 = {"Explicit valence for atom # 1 O, 3, is greater than permitted", 
			"Explicit valence for atom # 5 N, 5, is greater than permitted"};
	for (auto &query : errout3) {
		msgs1.push_back(query.message());
	}
	TEST_ASSERT(msgs1 == ans1);
	
	// testing molecule with multiple valency errors and only outputting
	// first error
	smi4 = "CO(C)CCN(=O)=O";
	unique_ptr<ROMol> m4( SmilesToMol(smi4, 0, false) );
	vector<ValidationErrorInfo> errout4 = vm.validate(*m4, false);
	std::vector<string> msgs2;
	std::vector<string> ans2 = {"Explicit valence for atom # 1 O, 3, is greater than permitted"};
	for (auto &query : errout4) {
		msgs2.push_back(query.message());
	}
	TEST_ASSERT(msgs2 == ans2);
}

void testMolVSValidation() {
	string file_root = getenv("RDBASE");
	file_root += "/Docs/Book/data/";
	MolVSValidation vm;

	string invalid_mol = file_root + "invalid.mol";
	unique_ptr<RWMol> m1( RDKit::MolFileToMol(invalid_mol, false) );
	unique_ptr<ROMol> res( new ROMol(*m1) );
	vector<ValidationErrorInfo> errout1 = vm.validate(*res, true);

	for (auto &query : errout1) {
	std::string msg = query.message();
	std::cout << msg << std::endl;
	TEST_ASSERT(msg == "Molecule is None.");
	}

}

int main() {
	//testRDKitValidation();
	testMolVSValidation();
	return 0;
}
