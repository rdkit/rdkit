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
	TEST_ASSERT(msg == "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is greater than permitted");
	}

	// testing for molecule with no atoms
	smi2 = "";
	unique_ptr<ROMol> m2( SmilesToMol(smi2, 0, false) );
	vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
	for (auto &query : errout2) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
	}

	// testing molecule with multiple valency errors
	smi3 = "CO(C)CCN(=O)=O";
	unique_ptr<ROMol> m3( SmilesToMol(smi3, 0, false) );
	vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
	std::vector<string> msgs1;
	std::vector<string> ans1 = {"INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is greater than permitted", 
			"INFO: [ValenceValidation] Explicit valence for atom # 5 N, 5, is greater than permitted"};
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
	std::vector<string> ans2 = {"INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is greater than permitted"};
	for (auto &query : errout4) {
		msgs2.push_back(query.message());
	}
	TEST_ASSERT(msgs2 == ans2);
}

void testMolVSValidation() {
	string smi1, smi2, smi3, smi4, smi5, smi6;
	MolVSValidation vm;

	// testing MolVSDefault
	// testing for molecule with no atoms
	smi1 = "";
	unique_ptr<ROMol> m1( SmilesToMol(smi1, 0, false) );
	vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
	for (auto &query : errout1) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
	}

	smi2 = "O=C([O-])c1ccccc1";
	unique_ptr<ROMol> m2( SmilesToMol(smi2, 0, false) );
	vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
	for (auto &query : errout2) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [NeutralValidation] Not an overall neutral system (-1)");
	}

	smi3 = "CN=[NH+]CN=N";
	unique_ptr<ROMol> m3( SmilesToMol(smi3, 0, false) );
	vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
	for (auto &query : errout3) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [NeutralValidation] Not an overall neutral system (+1)"); // fix to show + sign
	}

	smi4 = "[13CH4]";
	unique_ptr<ROMol> m4( SmilesToMol(smi4, 0, false) );
	vector<ValidationErrorInfo> errout4 = vm.validate(*m4, true);
	for (auto &query : errout4) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [IsotopeValidation] Molecule contains isotope 13C");
	}

	smi5 = "[2H]C(Cl)(Cl)Cl";
	unique_ptr<ROMol> m5( SmilesToMol(smi5, 0, false) );
	vector<ValidationErrorInfo> errout5 = vm.validate(*m5, true);
	for (auto &query : errout5) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [IsotopeValidation] Molecule contains isotope 2H");
	}

	smi6 = "[2H]OC([2H])([2H])[2H]";
	unique_ptr<ROMol> m6( SmilesToMol(smi6, 0, false) );
	vector<ValidationErrorInfo> errout6 = vm.validate(*m6, true);
	for (auto &query : errout6) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [IsotopeValidation] Molecule contains isotope 2H");
	}

	std::string smi7 = "COc1cccc(C=N[N-]C(N)=O)c1[O-].O.O.O.O=[U+2]=O";
	unique_ptr<ROMol> m7( SmilesToMol(smi7, 0, false) );
	vector<ValidationErrorInfo> errout7 = vm.validate(*m7, true);
	TEST_ASSERT(errout7.size() != 0);
	for (auto &query : errout7) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [FragmentValidation] water/hydroxide is present");
	}

	std::string smi8 = "CC(=O)O.NCC(=O)NCCCCCCCCCCNC(=O)CN";
	unique_ptr<ROMol> m8( SmilesToMol(smi8, 0, false) );
	vector<ValidationErrorInfo> errout8 = vm.validate(*m8, true);
	TEST_ASSERT(errout8.size() != 0);
	for (auto &query : errout8) {
	std::string msg = query.message();
	TEST_ASSERT(msg == "INFO: [FragmentValidation] acetate/acetic acid is present");
	}

	std::string smi9 = "N#CC(Br)(Br)C#N.[Br-].[K+]";
	unique_ptr<ROMol> m9( SmilesToMol(smi9, 0, false) );
	vector<ValidationErrorInfo> errout9 = vm.validate(*m9, true);
	std::vector<std::string> ans = {"INFO: [FragmentValidation] bromine is present", 
					"INFO: [FragmentValidation] potassium is present"};
	TEST_ASSERT(errout9.size() == ans.size());
	for (size_t i=0; i < errout9.size(); ++i) {
		TEST_ASSERT(errout9[i].message() == ans[i]);
	}

	std::string smi10 = "C1COCCO1.O=C(NO)NO";
	unique_ptr<ROMol> m10( SmilesToMol(smi10, 0, false) );
	vector<ValidationErrorInfo> errout10 = vm.validate(*m10, true);
	std::vector<std::string> ans10 = {"INFO: [FragmentValidation] 1,2-dimethoxyethane is present", 
					"INFO: [FragmentValidation] 1,4-dioxane is present"};
	TEST_ASSERT(errout10.size() == ans10.size());
	for (size_t i=0; i < errout10.size(); ++i) {
		TEST_ASSERT(errout10[i].message() == ans10[i]);
	}
}

void testAllowedAtomsValidation() {
//	std::vector<string> atoms = {"C", "N", "O"};
	std::vector<unsigned int> atoms = {6, 7, 8};
	std::vector<shared_ptr<Atom>> atomList;	

	for (auto &atom : atoms) {
		shared_ptr<Atom> a( new Atom(atom) );
		atomList.push_back(a);
	}

	AllowedAtomsValidation vm(atomList);
	std::string smi1;

	smi1 = "CC(=O)CF";
	unique_ptr<ROMol> m1( SmilesToMol(smi1) );
	vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
	for (auto &query : errout1) {
	std::string msg = query.message();
	std::cout << msg << std::endl;
	TEST_ASSERT(msg == "INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list");
	}

}

void testDisallowedAtomsValidation() {
//	std::vector<string> atoms = {"F", "Cl", "Br"};
	std::vector<unsigned int> atoms = {9, 17, 35};
	std::vector<shared_ptr<Atom>> atomList;	

	for (auto &atom : atoms) {
		shared_ptr<Atom> a( new Atom(atom) );
		atomList.push_back(a);
	}

	DisallowedAtomsValidation vm(atomList);
	std::string smi1;

	smi1 = "CC(=O)CF";
	unique_ptr<ROMol> m1( SmilesToMol(smi1) );
	vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
	for (auto &query : errout1) {
	std::string msg = query.message();
	std::cout << msg << std::endl;
	TEST_ASSERT(msg == "INFO: [DisallowedAtomsValidation] Atom F is in disallowedAtoms list");
	}
}

void testFragment() {
	string smi1, smi2, smi3, smi4, smi5, smi6;
	MolVSValidation vm;

	// testing MolVSValidation fragmentValidation
	// FragmentValidation should identify 1,2-dichloroethane.
	smi1 = "ClCCCl.c1ccccc1O";
	unique_ptr<ROMol> m1( SmilesToMol(smi1, 0, false) );
	vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
	for (auto &query : errout1) {
	std::string msg = query.message();
	std::cout << msg << std::endl;
	TEST_ASSERT(msg == "INFO: [FragmentValidation] 1,2-dichloroethane is present");
	}

	smi2 = "COCCOC.CCCBr";
	unique_ptr<ROMol> m2( SmilesToMol(smi2, 0, false) );
	vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
	for (auto &query : errout2) {
	std::string msg = query.message();
	std::cout << msg << std::endl;
	TEST_ASSERT(msg == "INFO: [FragmentValidation] 1,2-dimethoxyethane is present");
	}

}

int main() {
	testRDKitValidation();
	testMolVSValidation();
	testAllowedAtomsValidation();
	testDisallowedAtomsValidation();
	testFragment();
	return 0;
}
