#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogParams.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogUtils.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentRemover.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolStandardize/MolStandardize.h>

#include <iostream>
#include <fstream>

using namespace RDKit;
using namespace MolStandardize;

void test2() {
  std::string smi1, smi2, smi3, smi4, smi5, smi6, smi8, smi9, smi10, smi11,
      smi12;

  // testing parsing of fragment catalog
  std::string rdbase = getenv("RDBASE");
  std::string fgrpFile = rdbase +
                         "/Code/GraphMol/MolStandardize/FragmentCatalog/"
                         "data/fragmentPatterns.txt";
  auto* fparams = new FragmentCatalogParams(fgrpFile);
  unsigned int numfg = fparams->getNumFuncGroups();
  TEST_ASSERT(fparams->getNumFuncGroups() == 61);

  FragmentCatalog fcat(fparams);
  FragmentRemover fragremover;

  // single salt removal
  smi1 = "CN(C)C.Cl";
  std::shared_ptr<ROMol> m1(SmilesToMol(smi1));
  ROMol* remove = fragremover.remove(*m1, &fcat);
  std::cout << MolToSmiles(*remove) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove) == "CN(C)C");

  // double salt removal
  smi2 = "CN(C)C.Cl.Cl.Br";
  std::shared_ptr<ROMol> m2(SmilesToMol(smi2));
  ROMol* remove2 = fragremover.remove(*m2, &fcat);
  std::cout << MolToSmiles(*remove2) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove2) == "CN(C)C");

  // FragmentPatterns should match entire fragments only,
  // matches within larger fragments should be left
  smi3 = "CN(Br)Cl";
  std::shared_ptr<ROMol> m3(SmilesToMol(smi3));
  ROMol* remove3 = fragremover.remove(*m3, &fcat);
  std::cout << MolToSmiles(*remove3) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove3) == "CN(Cl)Br");

  // FragmentPatterns should match entire fragments only,
  // matches within larger fragments should be left
  smi4 = "CN(Br)Cl.Cl";
  std::shared_ptr<ROMol> m4(SmilesToMol(smi4));
  ROMol* remove4 = fragremover.remove(*m4, &fcat);
  std::cout << MolToSmiles(*remove4) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove4) == "CN(Cl)Br");

  // charged salts
  smi5 = "C[NH+](C)(C).[Cl-]";
  std::shared_ptr<ROMol> m5(SmilesToMol(smi5));
  ROMol* remove5 = fragremover.remove(*m5, &fcat);
  std::cout << MolToSmiles(*remove5) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove5) == "C[NH+](C)C");

  // Last match should be left.
  smi6 = "CC(=O)O.[Na]";
  std::shared_ptr<ROMol> m6(SmilesToMol(smi6));
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
  std::shared_ptr<ROMol> m8(SmilesToMol(smi8));
  ROMol* remove8 = fragremover.remove(*m8, &fcat);
  std::cout << MolToSmiles(*remove8) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove8) == "Cl.Cl");

  // Last match should be left.
  smi9 = "[Na+].OC(=O)Cc1ccc(CN)cc1.OS(=O)(=O)C(F)(F)F";
  std::shared_ptr<ROMol> m9(SmilesToMol(smi9));
  ROMol* remove9 = fragremover.remove(*m9, &fcat);
  std::cout << MolToSmiles(*remove9) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove9) == "NCc1ccc(CC(=O)O)cc1");

  // 1,4-Dioxane should be removed..
  smi10 = "c1ccccc1O.O1CCOCC1";
  std::shared_ptr<ROMol> m10(SmilesToMol(smi10));
  ROMol* remove10 = fragremover.remove(*m10, &fcat);
  std::cout << MolToSmiles(*remove10) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove10) == "Oc1ccccc1");

  // Benzene should be removed.
  smi11 = "c1ccccc1.CCCBr";
  std::shared_ptr<ROMol> m11(SmilesToMol(smi11));
  ROMol* remove11 = fragremover.remove(*m11, &fcat);
  std::cout << MolToSmiles(*remove11) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove11) == "CCCBr");

  // Various fragments should be removed.
  smi12 = "CC(NC1=CC=C(O)C=C1)=O.CCCCC.O.CCO.CCCO.C1CCCCC1.C1CCCCCC1";
  std::shared_ptr<ROMol> m12(SmilesToMol(smi12));
  ROMol* remove12 = fragremover.remove(*m12, &fcat);
  std::cout << MolToSmiles(*remove12) << std::endl;
  TEST_ASSERT(MolToSmiles(*remove12) == "CC(=O)Nc1ccc(O)cc1");
}

void test_largest_fragment() {
  std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7, smi8, smi9, smi10;
  LargestFragmentChooser lfragchooser;
  LargestFragmentChooser lfrag_preferOrg(true);
  CleanupParameters* params = new CleanupParameters;

  // Multiple organic fragments of different sizes.
  smi2 = "O=C(O)c1ccccc1.O=C(O)c1ccccc1.O=C(O)c1ccccc1";
  //	std::shared_ptr<ROMol> m2( SmilesToMol(smi2) );
  //	boost::shared_ptr<ROMol> lfrag2 = lfragchooser.choose(*m2);
  //	std::cout << MolToSmiles(*lfrag2) << std::endl;
  //	TEST_ASSERT(MolToSmiles(*lfrag2) == "O=C(O)c1ccccc1");
  std::shared_ptr<RWMol> m2(SmilesToMol(smi2));
  MolStandardize::fragmentParent(*m2, *params);
  TEST_ASSERT(MolToSmiles(*m2) == "O=C(O)c1ccccc1");

  // No organic fragments
  smi3 = "[N+](=O)([O-])[O-]";
  std::shared_ptr<ROMol> m3(SmilesToMol(smi3));
  boost::shared_ptr<ROMol> lfrag3 = lfragchooser.choose(*m3);
  std::cout << MolToSmiles(*lfrag3) << std::endl;
  TEST_ASSERT(MolToSmiles(*lfrag3) == "O=[N+]([O-])[O-]");

  // Larger inorganic should be chosen
  smi4 = "[N+](=O)([O-])[O-].[CH3+]";
  std::shared_ptr<ROMol> m4(SmilesToMol(smi4));
  boost::shared_ptr<ROMol> lfrag4 = lfragchooser.choose(*m4);
  std::cout << MolToSmiles(*lfrag4) << std::endl;
  TEST_ASSERT(MolToSmiles(*lfrag4) == "O=[N+]([O-])[O-]");

  // Smaller organic fragment should be chosen over larger inorganic fragment.
  smi5 = "[N+](=O)([O-])[O-].[CH3+]";
  std::shared_ptr<ROMol> m5(SmilesToMol(smi5));
  boost::shared_ptr<ROMol> lfrag5 = lfrag_preferOrg.choose(*m5);
  std::cout << MolToSmiles(*lfrag5) << std::endl;
  TEST_ASSERT(MolToSmiles(*lfrag5) == "[CH3+]");

  // Salt without charges.
  smi1 = "[Na].O=C(O)c1ccccc1";
  std::shared_ptr<RWMol> m1(SmilesToMol(smi1));
  MolStandardize::fragmentParent(*m1, *params);
  std::cout << MolToSmiles(*m1) << std::endl;
  //	TEST_ASSERT(MolToSmiles(*lfrag) == "CN(C)C");

  smi6 = "[Na]OC(=O)c1ccccc1";
  std::shared_ptr<RWMol> m6(SmilesToMol(smi6));
  //	MolStandardize::cleanup(*m6, *params);
  MolStandardize::fragmentParent(*m6, *params);
  std::cout << MolToSmiles(*m6) << std::endl;
  TEST_ASSERT(MolToSmiles(*m6) == "O=C([O-])c1ccccc1");

  smi7 = "c1ccccc1C(=O)O[Ca]OC(=O)c1ccccc1";
  std::shared_ptr<RWMol> m7(SmilesToMol(smi7));
  MolStandardize::fragmentParent(*m7, *params);
  std::cout << MolToSmiles(*m7) << std::endl;
  TEST_ASSERT(MolToSmiles(*m7) == "O=C([O-])c1ccccc1");

  smi8 = "[Pt](Cl)(Cl)(O)(O)(NC(C)C)NC(C)C";
  std::shared_ptr<RWMol> m8(SmilesToMol(smi8));
  MolStandardize::fragmentParent(*m8, *params);
  std::cout << MolToSmiles(*m8) << std::endl;
  TEST_ASSERT(MolToSmiles(*m8) == "CC(C)[NH-]");

  // Mercury containing compound.
  smi9 = "CC[Hg]SC1=C(C=CC=C1)C(=O)[O][Na]";
  std::shared_ptr<RWMol> m9(SmilesToMol(smi9));
  MolStandardize::fragmentParent(*m9, *params);
  std::cout << MolToSmiles(*m9) << std::endl;
  TEST_ASSERT(MolToSmiles(*m9) == "CC[Hg]Sc1ccccc1C(=O)[O-]");

  // Covalent bond with metal.
  smi10 = "[Ag]OC(=O)O[Ag]";
  std::shared_ptr<RWMol> m10(SmilesToMol(smi10));
  MolStandardize::fragmentParent(*m10, *params);
  std::cout << MolToSmiles(*m10) << std::endl;
  TEST_ASSERT(MolToSmiles(*m10) == "O=C([O-])[O-]");

  // params->preferOrganic = true;
}

int main() {
  test2();
  test_largest_fragment();
  return 0;
}
