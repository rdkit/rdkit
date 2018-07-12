#include "MolStandardize.h"
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogParams.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogUtils.h>
#include "Charge.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;
using namespace MolStandardize;

void testReionizer() {

	std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7;

	// testing parsing of acid base catalog
	std::string rdbase = getenv("RDBASE");
	std::string acidbaseFile = 
		rdbase + "/Code/GraphMol/MolStandardize/AcidBaseCatalog/test_data/acid_base_pairs.txt";
	std::shared_ptr<AcidBaseCatalogParams> abparams( new AcidBaseCatalogParams(acidbaseFile) );
	unsigned int npairs = abparams->getNumPairs();
	TEST_ASSERT(npairs == 33);

	AcidBaseCatalog abcat(abparams.get());
	Reionizer reionizer;

	// Test table salt.
	smi1 = "[Na].[Cl]";
	std::shared_ptr<ROMol> m1( SmilesToMol(smi1) );
	std::cout << "Before reionizing: " << MolToSmiles(*m1) << std::endl;
	ROMOL_SPTR reionized( reionizer.reionize(*m1, &abcat) );
	std::cout << MolToSmiles(*reionized) << std::endl;
	TEST_ASSERT(MolToSmiles(*reionized) == "[Cl-].[Na+]");

	// Test forced charge correction maintaining overall neutral charge.
	smi2 = "[Na].O=C(O)c1ccccc1";
	std::shared_ptr<ROMol> m2( SmilesToMol(smi2) );
	std::cout << "Before reionizing: " << MolToSmiles(*m2) << std::endl;
	ROMOL_SPTR reionized2( reionizer.reionize(*m2, &abcat) );
	std::cout << MolToSmiles(*reionized2) << std::endl;
	TEST_ASSERT(MolToSmiles(*reionized2) == "O=C([O-])c1ccccc1.[Na+]");
	
	// Test reionizer moves proton to weaker acid.
	smi3 = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O";
	std::shared_ptr<ROMol> m3( SmilesToMol(smi3) );
	std::cout << "Before reionizing: " << MolToSmiles(*m3) << std::endl;
	ROMOL_SPTR reionized3( reionizer.reionize(*m3, &abcat) );
	std::cout << MolToSmiles(*reionized3) << std::endl;
	TEST_ASSERT(MolToSmiles(*reionized3) == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");

	// Test reionizer moves proton to weaker acid.
	smi5 = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O";
	std::shared_ptr<ROMol> m5( SmilesToMol(smi5) );
	std::cout << "Before reionizing: " << MolToSmiles(*m5) << std::endl;
	ROMOL_SPTR reionized5( reionizer.reionize(*m5, &abcat) );
	std::cout << MolToSmiles(*reionized5) << std::endl;
	TEST_ASSERT(MolToSmiles(*reionized3) == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1");

	// Test charged carbon doesn't get recognised as alpha-carbon-hydrogen-keto.
	smi6 = "CCOC(=O)C(=O)[CH-]C#N";
	std::shared_ptr<ROMol> m6( SmilesToMol(smi6) );
	std::cout << "Before reionizing: " << MolToSmiles(*m6) << std::endl;
	ROMOL_SPTR reionized6( reionizer.reionize(*m6, &abcat) );
	std::cout << MolToSmiles(*reionized6) << std::endl;
	TEST_ASSERT(MolToSmiles(*reionized6) == "CCOC(=O)C(=O)[CH-]C#N");

	// TODO... can't make this work. Python SanitizeMol looks to correct... 
	// what is different with MolOps::sanitizeMol? 
	smi7 = "C[N+]1=C[CH-]N(C(=N)N)/C1=C/[N+](=O)[O-]";
	std::shared_ptr<ROMol> m7( SmilesToMol(smi7) );
	std::cout << "Before reionizing: " << MolToSmiles(*m7) << std::endl;
	ROMOL_SPTR reionized7( reionizer.reionize(*m7, &abcat) );
	std::cout << MolToSmiles(*reionized7) << std::endl;
//	TEST_ASSERT(MolToSmiles(*reionized7) == "C[N+]1=CCN(C(=N)N)C1=[C-][N+](=O)[O-]");

}

void testChargeParent() {
	std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10, 
					smi11, smi12;
	CleanupParameters params;
	// initialize CleanupParameters with preferOrganic=true
  CleanupParameters params_preferorg;
	params_preferorg.preferOrganic = true;

	// Test neutralization of ionized acids and bases.
	smi1 = "C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])";
	std::unique_ptr<RWMol> m1( SmilesToMol(smi1) );
	MolStandardize::chargeParent(*m1, params);
	std::cout << MolToSmiles(*m1) << std::endl;
	TEST_ASSERT(MolToSmiles(*m1) == "NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O");
/*
	// Test preservation of zwitterion.
	smi2 = "n(C)1cc[n+]2cccc([O-])c12";
	std::unique_ptr<RWMol> m2( SmilesToMol(smi2) );
	MolStandardize::chargeParent(*m2, params);
	std::cout << MolToSmiles(*m2) << std::endl;
	TEST_ASSERT(MolToSmiles(*m2) == "Cn1cc[n+]2cccc([O-])c12");

	// Choline should be left with a positive charge.
	smi3 = "C[N+](C)(C)CCO";
	std::unique_ptr<RWMol> m3( SmilesToMol(smi3) );
	MolStandardize::chargeParent(*m3, params);
	std::cout << MolToSmiles(*m3) << std::endl;
	TEST_ASSERT(MolToSmiles(*m3) == "C[N+](C)(C)CCO");

	// Hydrogen should be removed to give deanol as a charge parent.
	smi4 = "C[NH+](C)CCO";
	std::unique_ptr<RWMol> m4( SmilesToMol(smi4) );
	MolStandardize::chargeParent(*m4, params);
	std::cout << MolToSmiles(*m4) << std::endl;
	TEST_ASSERT(MolToSmiles(*m4) == "CN(C)CCO");

	// Sodium benzoate to benzoic acid.
	smi5 = "[Na+].O=C([O-])c1ccccc1";
	std::unique_ptr<RWMol> m5( SmilesToMol(smi5) );
	MolStandardize::chargeParent(*m5, params);
	std::cout << MolToSmiles(*m5) << std::endl;
	TEST_ASSERT(MolToSmiles(*m5) == "O=C(O)c1ccccc1");

	// Benzoate ion to benzoic acid.
	smi6 = "O=C([O-])c1ccccc1";
	std::unique_ptr<RWMol> m6( SmilesToMol(smi6) );
	MolStandardize::chargeParent(*m6, params);
	std::cout << MolToSmiles(*m6) << std::endl;
	TEST_ASSERT(MolToSmiles(*m6) == "O=C(O)c1ccccc1");

	// Charges in histidine should be neutralized.
	smi7 = "[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]";
	std::unique_ptr<RWMol> m7( SmilesToMol(smi7) );
	MolStandardize::chargeParent(*m7, params);
	std::cout << MolToSmiles(*m7) << std::endl;
	TEST_ASSERT(MolToSmiles(*m7) == "NC(Cc1cnc[nH]1)C(=O)O");

	// 
	smi8 = "C[NH+](C)(C).[Cl-]";
	std::unique_ptr<RWMol> m8( SmilesToMol(smi8) );
	MolStandardize::chargeParent(*m8, params);
	std::cout << MolToSmiles(*m8) << std::endl;
	TEST_ASSERT(MolToSmiles(*m8) == "CN(C)C");

	// No organic fragments.
	smi9 = "[N+](=O)([O-])[O-]";
	std::unique_ptr<RWMol> m9( SmilesToMol(smi9) );
	MolStandardize::chargeParent(*m9, params);
	std::cout << MolToSmiles(*m9) << std::endl;
	TEST_ASSERT(MolToSmiles(*m9) == "O=[N+]([O-])[O-]");

	// TODO switch prefer_organic=true
	// No organic fragments.
	smi10 = "[N+](=O)([O-])[O-]";
	std::unique_ptr<RWMol> m10( SmilesToMol(smi10) );
	MolStandardize::chargeParent(*m10, params_preferorg);
	std::cout << MolToSmiles(*m10) << std::endl;
	TEST_ASSERT(MolToSmiles(*m10) == "O=[N+]([O-])[O-]");

	// Larger inorganic fragment should be chosen.
	smi11 = "[N+](=O)([O-])[O-].[CH2]";
	std::unique_ptr<RWMol> m11( SmilesToMol(smi11) );
	MolStandardize::chargeParent(*m11, params);
	std::cout << MolToSmiles(*m11) << std::endl;
	TEST_ASSERT(MolToSmiles(*m11) == "O=[N+]([O-])[O-]");

	// TODO prefer_organic=true
	// Smaller organic fragment should be chosen over larger inorganic fragment.
	smi12 = "[N+](=O)([O-])[O-].[CH2]";
	std::unique_ptr<RWMol> m12( SmilesToMol(smi12) );
	MolStandardize::chargeParent(*m12, params_preferorg);
	std::cout << MolToSmiles(*m12) << std::endl;
	TEST_ASSERT(MolToSmiles(*m12) == "[CH2]");
*/
}

int main() {
	testReionizer();
	testChargeParent();
	return 0;
}
