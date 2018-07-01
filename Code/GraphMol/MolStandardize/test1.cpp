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

void testCleanup(){
	string smi1, smi2, smi3, smi4;
	CleanupParameters params;

	// Test covalent metal is disconnected during standardize.
	smi1 = "CCC(=O)O[Na]";
	unique_ptr<RWMol> m1( SmilesToMol(smi1) );
	MolStandardize::cleanup(*m1, params);
	std::cout << MolToSmiles(*m1) << std::endl;
	TEST_ASSERT(MolToSmiles(*m1) == "CCC(=O)[O-].[Na+]");

	// Test metal ion is untouched during standardize.
	smi2 = "CCC(=O)[O-].[Na+]";
	unique_ptr<RWMol> m2( SmilesToMol(smi2) );
	MolStandardize::cleanup(*m2, params);
	TEST_ASSERT(MolToSmiles(*m2) == "CCC(=O)[O-].[Na+]");

	// Test Hg is disconnected from O during standardize.
	smi3 = "CCC(=O)O[Hg]";
	unique_ptr<RWMol> m3( SmilesToMol(smi3) );
	MolStandardize::cleanup(*m3, params);
	TEST_ASSERT(MolToSmiles(*m3) == "CCC(=O)[O-].[Hg+]")

	// Test dimethylmercury is not disconnected during standardize.
	smi4 = "C[Hg]C";
	unique_ptr<RWMol> m4( SmilesToMol(smi4) );
	MolStandardize::cleanup(*m4, params);
	TEST_ASSERT(MolToSmiles(*m4) == "C[Hg]C")

}

void testStandardize(){
	string smi;
	CleanupParameters params;

	smi = "CCCC(=O)O";
	unique_ptr<RWMol> m( SmilesToMol(smi) );
	unique_ptr<RWMol> m2( SmilesToMol(smi) );
	TEST_ASSERT(m);

	// cleanup function
	MolStandardize::cleanup(*m, params);
	// testing nothing has changed
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty tautomer parent function
	MolStandardize::tautomerParent(*m, params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty isotope parent function
	MolStandardize::isotopeParent(*m, params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty charge parent function
	MolStandardize::chargeParent(*m, params);
	TEST_ASSERT(MolToSmiles(*m) == MolToSmiles(*m2));

	// empty super parent function
	MolStandardize::superParent(*m, params);
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

void testNormalize() { 
	string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10, smi11, 
				smi12, smi13, smi14, smi15, smi16, smi17; 
	CleanupParameters params;

	std::string rdbase = getenv("RDBASE");
	std::string transformFile = 
		rdbase + "/Code/GraphMol/MolStandardize/TransformCatalog/test_data/normalizations.txt";
	params.normalizations = transformFile;
	std::cout << params.normalizations << std::endl;

	// Normalize nitro group.
	smi1 = "C1(=CC=CC=C1)[N+](=O)[O-]";
	unique_ptr<RWMol> m1( SmilesToMol(smi1) );
	MolStandardize::cleanup(*m1, params);
	std::cout << MolToSmiles(*m1) << std::endl;
	TEST_ASSERT(MolToSmiles(*m1) == "O=[N+]([O-])c1ccccc1");

	// Normalize nitro group.
	smi2 = "O=[N](=O)c1ccccc1";
	unique_ptr<RWMol> m2( SmilesToMol(smi2) );
	MolStandardize::cleanup(*m2, params);
	std::cout << MolToSmiles(*m2) << std::endl;
	TEST_ASSERT(MolToSmiles(*m2) == "O=[N+]([O-])c1ccccc1");

	// Normalize nitro group.
	smi3 = "[O-][N+](=O)c1ccccc1";
	unique_ptr<RWMol> m3( SmilesToMol(smi3) );
	MolStandardize::cleanup(*m3, params);
	std::cout << MolToSmiles(*m3) << std::endl;
	TEST_ASSERT(MolToSmiles(*m3) == "O=[N+]([O-])c1ccccc1");

	// Normalize nitro group.
	smi4 = "[N](=O)(=O)O";
	unique_ptr<RWMol> m4( SmilesToMol(smi4) );
	MolStandardize::cleanup(*m4, params);
	std::cout << MolToSmiles(*m4) << std::endl;
	TEST_ASSERT(MolToSmiles(*m4) == "O=[N+]([O-])O");

	// Normalize nitro group.
	smi5 = "O[N+](=O)[O-]";
	unique_ptr<RWMol> m5( SmilesToMol(smi5) );
	MolStandardize::cleanup(*m5, params);
	std::cout << MolToSmiles(*m5) << std::endl;
	TEST_ASSERT(MolToSmiles(*m5) == "O=[N+]([O-])O");

	// Normalize pyridine oxide.
	smi6 = "C1=[N](C=CC=C1)=O";
	unique_ptr<RWMol> m6( SmilesToMol(smi6) );
	MolStandardize::cleanup(*m6, params);
	std::cout << MolToSmiles(*m6) << std::endl;
	TEST_ASSERT(MolToSmiles(*m6) == "[O-][n+]1ccccc1");

	// Normalize pyridine oxide.
	smi7 = "O=n1ccccc1";
	unique_ptr<RWMol> m7( SmilesToMol(smi7) );
	MolStandardize::cleanup(*m7, params);
	std::cout << MolToSmiles(*m7) << std::endl;
	TEST_ASSERT(MolToSmiles(*m7) == "[O-][n+]1ccccc1");

	// normalize sulfone.
	smi8 = "C[S+2]([O-])([O-])C";
	unique_ptr<RWMol> m8( SmilesToMol(smi8) );
	MolStandardize::cleanup(*m8, params);
	std::cout << MolToSmiles(*m8) << std::endl;
	TEST_ASSERT(MolToSmiles(*m8) == "CS(C)(=O)=O");

	// normalize sulfone.
	smi9 = "C[S+2]([O-])([O-])O";
	unique_ptr<RWMol> m9( SmilesToMol(smi9) );
	MolStandardize::cleanup(*m9, params);
	std::cout << MolToSmiles(*m9) << std::endl;
	TEST_ASSERT(MolToSmiles(*m9) == "CS(=O)(=O)O");

	// normalize sulfoxide..
	smi10 = "CS(=O)C";
	unique_ptr<RWMol> m10( SmilesToMol(smi10) );
	MolStandardize::cleanup(*m10, params);
	std::cout << MolToSmiles(*m10) << std::endl;
	TEST_ASSERT(MolToSmiles(*m10) == "C[S+](C)[O-]");

	// normalize sulfoxide.
	smi11 = "COC1=CC2=C(C=C1)[N]C(=N2)[S](=O)CC3=C(C(=C(C=N3)C)OC)C";
	unique_ptr<RWMol> m11( SmilesToMol(smi11) );
	MolStandardize::cleanup(*m11, params);
	std::cout << MolToSmiles(*m11) << std::endl;
	TEST_ASSERT(MolToSmiles(*m11) == 
									"COc1ccc2c(c1)N=C([S+]([O-])Cc1ncc(C)c(OC)c1C)[N]2");

	// Normalize sulfoxide.
	smi12 = "COc1ccc2c(c1)nc([nH]2)S(=O)Cc1ncc(c(c1C)OC)C";
	unique_ptr<RWMol> m12( SmilesToMol(smi12) );
	MolStandardize::cleanup(*m12, params);
	std::cout << MolToSmiles(*m12) << std::endl;
	TEST_ASSERT(MolToSmiles(*m12) == 
									"COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1");

	// Normalize azide.
	smi13 = "C1(=CC=C(C=C1)N)N=[N]#N";
	unique_ptr<RWMol> m13( SmilesToMol(smi13) );
	MolStandardize::cleanup(*m13, params);
	std::cout << MolToSmiles(*m13) << std::endl;
	TEST_ASSERT(MolToSmiles(*m13) == 
									"[N-]=[N+]=Nc1ccc(N)cc1");

	// Normalize diazo.
	smi14 = "[N](#N)=C1C(NC(N=C1)=O)=O";
	unique_ptr<RWMol> m14( SmilesToMol(smi14) );
	MolStandardize::cleanup(*m14, params);
	std::cout << MolToSmiles(*m14) << std::endl;
	TEST_ASSERT(MolToSmiles(*m14) == "[N-]=[N+]=C1C=NC(=O)NC1=O");

	// Normalize phosphate.
	smi15 = "C1=NC=C([N]1)CO[P+]([O-])([O-])[O-]";
	unique_ptr<RWMol> m15( SmilesToMol(smi15) );
	MolStandardize::cleanup(*m15, params);
	std::cout << MolToSmiles(*m15) << std::endl;
	TEST_ASSERT(MolToSmiles(*m15) == "O=P([O-])([O-])OCC1=CN=C[N]1");

	// Normalize hydrazine-diazonium.
	smi16 = "CNNC[N+]#N";
	unique_ptr<RWMol> m16( SmilesToMol(smi16) );
	MolStandardize::cleanup(*m16, params);
	std::cout << MolToSmiles(*m16) << std::endl;
	TEST_ASSERT(MolToSmiles(*m16) == "CN=[NH+]CN=N");

	// Normalize amidinium.
	smi17 = "[C+](C)(N)N";
	unique_ptr<RWMol> m17( SmilesToMol(smi17) );
	MolStandardize::cleanup(*m17, params);
	std::cout << MolToSmiles(*m17) << std::endl;
	TEST_ASSERT(MolToSmiles(*m17) == "CC(N)=[NH2+]");
}

void testNormalizeMultiFrags() {
	string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10, smi11, 
				smi12, smi13, smi14, smi15, smi16, smi17; 
	CleanupParameters params;

	std::string rdbase = getenv("RDBASE");
	std::string transformFile = 
		rdbase + "/Code/GraphMol/MolStandardize/TransformCatalog/test_data/normalizations.txt";
	params.normalizations = transformFile;

	// All fragments should stay if one gets transformed by normalization.
	smi1 = "[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1";
	unique_ptr<RWMol> m1( SmilesToMol(smi1) );
	MolStandardize::cleanup(*m1, params);
	std::cout << MolToSmiles(*m1) << std::endl;
	TEST_ASSERT(MolToSmiles(*m1) == "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]");

	// All fragments should stay if one gets transformed by normalization.
	smi2 = "[Na+].[O-]C(=O)c1ccc(C[S+2]([O-])([O-]))cc1";
	unique_ptr<RWMol> m2( SmilesToMol(smi2) );
	MolStandardize::cleanup(*m2, params);
	std::cout << MolToSmiles(*m2) << std::endl;
	TEST_ASSERT(MolToSmiles(*m2) == "O=C([O-])c1ccc(C[S](=O)=O)cc1.[Na+]");

	// Recombine non-aromatic 1,3-separated charges.
	smi3 = "C[N-]C(C)=[N+](C)C";
	unique_ptr<RWMol> m3( SmilesToMol(smi3) );
	MolStandardize::cleanup(*m3, params);
	std::cout << MolToSmiles(*m3) << std::endl;
	TEST_ASSERT(MolToSmiles(*m3) == "CN=C(C)N(C)C");

	// Recombine aromatic 1,3-separated charges.
	smi4 = "[n-]1c(=[N+](C)C)cccc1";
	unique_ptr<RWMol> m4( SmilesToMol(smi4) );
	MolStandardize::cleanup(*m4, params);
	std::cout << MolToSmiles(*m4) << std::endl;
	TEST_ASSERT(MolToSmiles(*m4) == "CN(C)c1ccccn1");

	// Recombine aromatic 1,3-separated charges.
	smi5 = "C[n+]1c([N-](C))cccc1";
	unique_ptr<RWMol> m5( SmilesToMol(smi5) );
	MolStandardize::cleanup(*m5, params);
	std::cout << MolToSmiles(*m5) << std::endl;
	TEST_ASSERT(MolToSmiles(*m5) == "CN=c1ccccn1C");

	// Recombine aromatic 1,3-separated charges to form pyrimidone.
	smi6 = "[O-]c1[n+](C)cccc1";
	unique_ptr<RWMol> m6( SmilesToMol(smi6) );
	MolStandardize::cleanup(*m6, params);
	std::cout << MolToSmiles(*m6) << std::endl;
	TEST_ASSERT(MolToSmiles(*m6) == "Cn1ccccc1=O");

	// Recombine aromatic 1,3-separated charges to form pyrimidone.
	smi7 = "COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c([O-])[n+](C)c3c2cc1OC";
	unique_ptr<RWMol> m7( SmilesToMol(smi7) );
	MolStandardize::cleanup(*m7, params);
	std::cout << MolToSmiles(*m7) << std::endl;
	TEST_ASSERT(MolToSmiles(*m7) == 
									"COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c(=O)n(C)c3c2cc1OC");

	// Recombine non-aromatic 1,5-separated charges.
	smi8 = "C[N-]C=CC=[N+](C)C";
	unique_ptr<RWMol> m8( SmilesToMol(smi8) );
	MolStandardize::cleanup(*m8, params);
	std::cout << MolToSmiles(*m8) << std::endl;
	TEST_ASSERT(MolToSmiles(*m8) == "CN=CC=CN(C)C");

	// Recombine aromatic 1,5-separated charges.
	smi9 = "[n-]1ccc(=[N+](C)C)cc1";
	unique_ptr<RWMol> m9( SmilesToMol(smi9) );
	MolStandardize::cleanup(*m9, params);
	std::cout << MolToSmiles(*m9) << std::endl;
	TEST_ASSERT(MolToSmiles(*m9) == "CN(C)c1ccncc1");

	// Recombine aromatic 1,5-separated charges.
	smi10 = "C[n+]1ccc([N-]C)cc1";
	unique_ptr<RWMol> m10( SmilesToMol(smi10) );
	MolStandardize::cleanup(*m10, params);
	std::cout << MolToSmiles(*m10) << std::endl;
	TEST_ASSERT(MolToSmiles(*m10) == "CN=c1ccn(C)cc1");

	// Shift positive charge from nonprotonated to protonated atom.
	smi11 = "CNC=[N+](C)C";
	unique_ptr<RWMol> m11( SmilesToMol(smi11) );
	MolStandardize::cleanup(*m11, params);
	std::cout << MolToSmiles(*m11) << std::endl;
	TEST_ASSERT(MolToSmiles(*m11) == "C[NH+]=CN(C)C");

	// Shift positive charge from nonprotonated to protonated atom."
	smi12 = "CNC=CC=[N+](C)C";
	unique_ptr<RWMol> m12( SmilesToMol(smi12) );
	MolStandardize::cleanup(*m12, params);
	std::cout << MolToSmiles(*m12) << std::endl;
	TEST_ASSERT(MolToSmiles(*m12) == "C[NH+]=CC=CN(C)C");

	// Shift positive charge from nonprotonated to protonated atom."
	smi13 = "[nH]1ccc(=[N+](C)C)cc1";
	unique_ptr<RWMol> m13( SmilesToMol(smi13) );
	MolStandardize::cleanup(*m13, params);
	std::cout << MolToSmiles(*m13) << std::endl;
	TEST_ASSERT(MolToSmiles(*m13) == "CN(C)c1cc[nH+]cc1");

	// Ensure no transforms inadvertently breaks open rings.
	smi14 = "[O-]C1=CC=CC2=CC=CC=[N+]12";
	unique_ptr<RWMol> m14( SmilesToMol(smi14) );
	MolStandardize::cleanup(*m14, params);
	std::cout << MolToSmiles(*m14) << std::endl;
	TEST_ASSERT(MolToSmiles(*m14) == "O=c1cccc2ccccn12");

	// Shift positive charge from nonprotonated to protonated atom.
	smi15 = "[nH]1c(=[N+](C)C)cccc1";
	unique_ptr<RWMol> m15( SmilesToMol(smi15) );
	MolStandardize::cleanup(*m15, params);
	std::cout << MolToSmiles(*m15) << std::endl;
	TEST_ASSERT(MolToSmiles(*m15) == "CN(C)c1cccc[nH+]1");

}

void testCharge() {
	std::string smi1, smi2, smi3, smi4;

	// Reionization should not infinitely loop forever on these molecules.
	smi1 = "CCCCCCCCCCCCCCCCCC(=O)CC(=C)C(=O)O[Ti](=O)(OC(C)C)C(C)C";
	std::string ss =	MolStandardize::standardizeSmiles(smi1);
	std::cout << ss << std::endl;
	TEST_ASSERT(ss == "C=C(CC(=O)[CH-]CCCCCCCCCCCCCCCC)C(=O)[O-].CC(C)[O-].CCC.[O-2].[Ti+5]");

	// Reionization should not infinitely loop forever on these molecules.
	smi2 = "OP(=O)(O)[O-].OP(=O)([O-])[O-].[O-]S(=O)(=O)[O-].[Na+].[Na+].[Na+].[Mg+2].[Cl-].[Cl-].[K+].[K+]";
	std::string ss2 = MolStandardize::standardizeSmiles(smi2);
	TEST_ASSERT(ss2 == "O=P([O-])(O)O.O=P([O-])([O-])O.O=S(=O)([O-])[O-].[Cl-].[Cl-].[K+].[K+].[Mg+2].[Na+].[Na+].[Na+]");

  // Charge parent!!
  // Test reionizer moves proton to weaker acid.
  smi3 = "[Na].[Na]";
	std::string ss3 = MolStandardize::standardizeSmiles(smi3);
  TEST_ASSERT(ss3 == "[Na+].[Na+]");

	// TODO: Arguably should become selenite ion... O=[Se]([O-])[O-]. 
	// Need an AcidBasePair?
  smi4 = "[Na].[Na].O[Se](O)=O";
	std::string ss4 = MolStandardize::standardizeSmiles(smi4);
  TEST_ASSERT(ss4 == "O=[Se](O)O.[Na+].[Na+]");

}

int main() {
	testCleanup();
	testStandardize();
	testMetalDisconnector();
	testNormalize();
	testNormalizeMultiFrags();
	testCharge();
	return 0;
}
