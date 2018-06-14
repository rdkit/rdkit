#include "FragmentCalaogParams.h"
#include "FragmentCalaogUtils.h"
#include "FragmentRemover.h"
//#include "FragmentCalaogEntry.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <iostream>
#include <fstream>

using namespace RDKit;
using namespace MolStandardize;

void test2() {

	std::string smi1, smi2, smi3, smi4, smi5, smi6, smi8, smi9, smi10, smi11,
		smi12;

	// testing parsing of fragment catalog
	std::string rdbase = getenv("RDBASE");
	std::string fgrpFile = 
		rdbase + "/Code/GraphMol/MolStandardize/FragmentCatalog/test_data/fragmentPatterns.txt";
	auto *fparams = new FragmentCatalogParams(fgrpFile);
	unsigned int numfg = fparams->getNumFuncGroups();
	TEST_ASSERT(fparams->getNumFuncGroups() == 61);

	FragmentCatalog fcat(fparams);
	FragmentRemover fragremover;

	// single salt removal
	smi1 = "CN(C)C.Cl";
	std::shared_ptr<ROMol> m1( SmilesToMol(smi1) );
	ROMol* remove = fragremover.remove(*m1, &fcat);
	std::cout << MolToSmiles(*remove) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove) == "CN(C)C");

	// double salt removal
	smi2 = "CN(C)C.Cl.Cl.Br";
	std::shared_ptr<ROMol> m2( SmilesToMol(smi2) );
	ROMol* remove2 = fragremover.remove(*m2, &fcat);
	std::cout << MolToSmiles(*remove2) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove2) == "CN(C)C");
	
	// FragmentPatterns should match entire fragments only, 
	// matches within larger fragments should be left
	smi3 = "CN(Br)Cl";
	std::shared_ptr<ROMol> m3( SmilesToMol(smi3) );
	ROMol* remove3 = fragremover.remove(*m3, &fcat);
	std::cout << MolToSmiles(*remove3) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove3) == "CN(Cl)Br");

	// FragmentPatterns should match entire fragments only, 
	// matches within larger fragments should be left
	smi4 = "CN(Br)Cl.Cl";
	std::shared_ptr<ROMol> m4( SmilesToMol(smi4) );
	ROMol* remove4 = fragremover.remove(*m4, &fcat);
	std::cout << MolToSmiles(*remove4) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove4) == "CN(Cl)Br");

	// charged salts
	smi5 = "C[NH+](C)(C).[Cl-]";
	std::shared_ptr<ROMol> m5( SmilesToMol(smi5) );
	ROMol* remove5 = fragremover.remove(*m5, &fcat);
	std::cout << MolToSmiles(*remove5) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove5) == "C[NH+](C)C");

	// Last match should be left.
	smi6 = "CC(=O)O.[Na]";
	std::shared_ptr<ROMol> m6( SmilesToMol(smi6) );
	ROMol* remove6 = fragremover.remove(*m6, &fcat);
	std::cout << MolToSmiles(*remove6) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove6) == "CC(=O)O");

	// Last match should be removed.
	FragmentRemover fr_noleavelast(false);
	ROMol* remove7 = fr_noleavelast.remove(*m6, &fcat);
	std::cout << MolToSmiles(*remove7) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove7) == "");

	// Multiple identical last fragments should all be left.
	smi8 = "Cl.Cl";
	std::shared_ptr<ROMol> m8( SmilesToMol(smi8) );
	ROMol* remove8 = fragremover.remove(*m8, &fcat);
	std::cout << MolToSmiles(*remove8) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove8) == "Cl.Cl");

	// Last match should be left.
	smi9 = "[Na+].OC(=O)Cc1ccc(CN)cc1.OS(=O)(=O)C(F)(F)F";
	std::shared_ptr<ROMol> m9( SmilesToMol(smi9) );
	ROMol* remove9 = fragremover.remove(*m9, &fcat);
	std::cout << MolToSmiles(*remove9) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove9) == "NCc1ccc(cc1)CC(=O)O");

	// 1,4-Dioxane should be removed..
	smi10 = "c1ccccc1O.O1CCOCC1";
	std::shared_ptr<ROMol> m10( SmilesToMol(smi10) );
	ROMol* remove10 = fragremover.remove(*m10, &fcat);
	std::cout << MolToSmiles(*remove10) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove10) == "Oc1ccccc1");

	// Benzene should be removed.
	smi11 = "c1ccccc1.CCCBr";
	std::shared_ptr<ROMol> m11( SmilesToMol(smi11) );
	ROMol* remove11 = fragremover.remove(*m11, &fcat);
	std::cout << MolToSmiles(*remove11) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove11) == "CCCBr");

	// Various fragments should be removed.
	smi12 = "CC(NC1=CC=C(O)C=C1)=O.CCCCC.O.CCO.CCCO.C1CCCCC1.C1CCCCCC1";
	std::shared_ptr<ROMol> m12( SmilesToMol(smi12) );
	ROMol* remove12 = fragremover.remove(*m12, &fcat);
	std::cout << MolToSmiles(*remove12) << std::endl;
	TEST_ASSERT(MolToSmiles(*remove12) == "CC(=O)Nc1ccc(O)cc1");
}

int main() {
//	test1();
	test2();
	return 0;
}
