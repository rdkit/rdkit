#include "MolStandardize.h"
#include "Metal.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>

#include <iostream>

using namespace RDKit;
using namespace std;

void testCleaup(){
	string smi1, smi2, smi3, smi4;
	CleanupParameters *params;

	// Test covalent metal is disconnected during standardize.
	smi1 = "CCC(=O)O[Na]";
	unique_ptr<RWMol> m1( SmilesToMol(smi1) );
	MolStandardize::cleanup(*m1, *params);
	TEST_ASSERT(MolToSmiles(*m1) == "CCC(=O)[O-].[Na+]");

	// Test metal ion is untouched during standardize.
	smi2 = "CCC(=O)[O-].[Na+]";
	unique_ptr<RWMol> m2( SmilesToMol(smi2) );
	MolStandardize::cleanup(*m2, *params);
	TEST_ASSERT(MolToSmiles(*m2) == "CCC(=O)[O-].[Na+]");

	// Test Hg is disconnected from O during standardize.
	smi3 = "CCC(=O)O[Hg]";
	unique_ptr<RWMol> m3( SmilesToMol(smi3) );
	MolStandardize::cleanup(*m3, *params);
	TEST_ASSERT(MolToSmiles(*m3) == "CCC(=O)[O-].[Hg+]")

	// Test dimethylmercury is not disconnected during standardize.
	smi4 = "C[Hg]C";
	unique_ptr<RWMol> m4( SmilesToMol(smi4) );
	MolStandardize::cleanup(*m4, *params);
	TEST_ASSERT(MolToSmiles(*m4) == "C[Hg]C")

}

void testStandardize(){
	string smi;
	CleanupParameters *params;

	smi = "CCCC(=O)O";
	unique_ptr<RWMol> m( SmilesToMol(smi) );
	unique_ptr<RWMol> m2( SmilesToMol(smi) );
	TEST_ASSERT(m);

	// cleanup function
	MolStandardize::cleanup(*m, *params);
	// testing nothing has changed
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty tautomer parent function
	MolStandardize::tautomerParent(*m, *params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty isotope parent function
	MolStandardize::isotopeParent(*m, *params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty charge parent function
	MolStandardize::chargeParent(*m, *params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty super parent function
	MolStandardize::superParent(*m, *params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));
}

void testMetalDisconnector(){

	MolStandardize::MetalDisconnector md;

	string smi1 = "CCC(=O)O[Na]";
	unique_ptr<RWMol> m1( SmilesToMol(smi1) );
	TEST_ASSERT(m1);
	md.disconnect(*m1);
	TEST_ASSERT(MolToSmiles(*m1) == "CCC(=O)[O-].[Na+]");

	string smi2 = "[Na]OC(=O)CCC(=O)O[Na]";
	unique_ptr<RWMol> m2( SmilesToMol(smi2) );
	TEST_ASSERT(m2);
	md.disconnect(*m2);
	TEST_ASSERT(MolToSmiles(*m2) == "O=C([O-])CCC(=O)[O-].[Na+].[Na+]");

	string smi3 = "c1ccccc1[Mg]Br";
	unique_ptr<RWMol> m3( SmilesToMol(smi3) );
	TEST_ASSERT(m3);
	md.disconnect(*m3);
	TEST_ASSERT(MolToSmiles(*m3) == "Br[Mg]c1ccccc1");

	string smi4 = "Br[Mg]c1ccccc1CCC(=O)O[Na]";
        unique_ptr<RWMol> m4( SmilesToMol(smi4) );
        TEST_ASSERT(m4);
        md.disconnect(*m4);
        TEST_ASSERT(MolToSmiles(*m4) == "O=C([O-])CCc1ccccc1[Mg]Br.[Na+]");

}
int main() {
	testCleaup();
	testStandardize();
	testMetalDisconnector();
	return 0;
}
