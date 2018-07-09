#include "Tautomer.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;
using namespace MolStandardize;

void test1() {
	std::string rdbase = getenv("RDBASE");
	std::string tautomerFile = 
		rdbase + "/Code/GraphMol/MolStandardize/TautomerCatalog/test_data/tautomerTransforms.in";
	auto *tautparams = new TautomerCatalogParams(tautomerFile);
	unsigned int ntautomers = tautparams->getNumTautomers();
	TEST_ASSERT(ntautomers == 34);

	TautomerCatalog tautcat(tautparams);
	TautomerEnumerator te;

	// Enumerate 1,3 keto/enol tautomer.
	std::string smi1 = "C1(=CCCCC1)O";
	std::shared_ptr<ROMol> m1( SmilesToMol(smi1) );
	std::vector<std::string> res = te.enumerate(*m1, &tautcat);
	std::vector<std::string> ans = {"O=C1CCCCC1", "OC1=CCCCC1"};
	TEST_ASSERT(res == ans);

	// Enumerate 1,3 keto/enol tautomer.
	std::string smi2 = "C1(CCCCC1)=O";
	std::shared_ptr<ROMol> m2( SmilesToMol(smi2) );
	std::vector<std::string> res2 = te.enumerate(*m2, &tautcat);
	std::vector<std::string> ans2 = {"O=C1CCCCC1", "OC1=CCCCC1"};
	TEST_ASSERT(res2 == ans2);

	// Enumerate acetophenone keto/enol tautomer.
	std::string smi3 = "C(=C)(O)C1=CC=CC=C1";
	std::shared_ptr<ROMol> m3( SmilesToMol(smi3) );
	std::vector<std::string> res3 = te.enumerate(*m3, &tautcat);
	std::vector<std::string> ans3 = {"C=C(O)c1ccccc1", "CC(=O)c1ccccc1"};
	TEST_ASSERT(res3 == ans3);

	// Enumerate acetone keto/enol tautomer.
	std::string smi4 = "CC(C)=O";
	std::shared_ptr<ROMol> m4( SmilesToMol(smi4) );
	std::vector<std::string> res4 = te.enumerate(*m4, &tautcat);
	std::vector<std::string> ans4 = {"C=C(C)O", "CC(C)=O"};
	TEST_ASSERT(res4 == ans4);

	// keto/enol tautomer
	std::string smi5 = "OC(C)=C(C)C";
	std::shared_ptr<ROMol> m5( SmilesToMol(smi5) );
	std::vector<std::string> res5 = te.enumerate(*m5, &tautcat);
	std::vector<std::string> ans5 = {"C=C(O)C(C)C", "CC(=O)C(C)C", "CC(C)=C(C)O"};
	TEST_ASSERT(res5 == ans5);

	// 1-phenyl-2-propanone enol/keto
	std::string smi6 = "c1(ccccc1)CC(=O)C";
	std::shared_ptr<ROMol> m6( SmilesToMol(smi6) );
	std::vector<std::string> res6 = te.enumerate(*m6, &tautcat);
	std::vector<std::string> ans6 = {"C=C(O)Cc1ccccc1", "CC(=O)Cc1ccccc1", "CC(O)=Cc1ccccc1"};
	TEST_ASSERT(res6 == ans6);
	
	// 1,5 keto/enol tautomer
	std::string smi7 = "Oc1nccc2cc[nH]c(=N)c12";
	std::shared_ptr<ROMol> m7( SmilesToMol(smi7) );
	std::vector<std::string> res7 = te.enumerate(*m7, &tautcat);
	std::vector<std::string> ans7 = {
					"N=c1[nH]ccc2cc[nH]c(=O)c12", 
					"N=c1[nH]ccc2ccnc(O)c12", 
					"N=c1nccc2cc[nH]c(O)c1-2",
					"Nc1[nH]ccc2ccnc(=O)c1-2", 
					"Nc1nccc2cc[nH]c(=O)c12",
					"Nc1nccc2ccnc(O)c12"};
	TEST_ASSERT(res7 == ans7);
	
	// 1,5 keto/enol tautomer
	std::string smi8 = "C1(C=CCCC1)=O";
	std::shared_ptr<ROMol> m8( SmilesToMol(smi8) );
	std::vector<std::string> res8 = te.enumerate(*m8, &tautcat);
	std::vector<std::string> ans8 = {
					"O=C1C=CCCC1", 
					"O=C1CC=CCC1", 
					"OC1=CC=CCC1", 
					"OC1=CCC=CC1", 
					"OC1=CCCC=C1"};
	TEST_ASSERT(res8 == ans8);
	
	// 1,5 keto/enol tautomer
	std::string smi9 = "C1(=CC=CCC1)O";
	std::shared_ptr<ROMol> m9( SmilesToMol(smi9) );
	std::vector<std::string> res9 = te.enumerate(*m9, &tautcat);
	std::vector<std::string> ans9 = {
					"O=C1C=CCCC1", 
					"O=C1CC=CCC1", 
					"OC1=CC=CCC1", 
					"OC1=CCC=CC1", 
					"OC1=CCCC=C1"};
	TEST_ASSERT(res9 == ans9);
	
	// aliphatic imine tautomer
	std::string smi10 = "C1(CCCCC1)=N";
	std::shared_ptr<ROMol> m10( SmilesToMol(smi10) );
	std::vector<std::string> res10 = te.enumerate(*m10, &tautcat);
	std::vector<std::string> ans10 = {"N=C1CCCCC1", "NC1=CCCCC1"};
	TEST_ASSERT(res10 == ans10);
	
	// aliphatic imine tautomer
	std::string smi11 = "C1(=CCCCC1)N";
	std::shared_ptr<ROMol> m11( SmilesToMol(smi11) );
	std::vector<std::string> res11 = te.enumerate(*m11, &tautcat);
	std::vector<std::string> ans11 = {"N=C1CCCCC1", "NC1=CCCCC1"};
	TEST_ASSERT(res11 == ans11);

	// special imine tautomer
	std::string smi12 = "C1(C=CC=CN1)=CC";
	std::shared_ptr<ROMol> m12( SmilesToMol(smi12) );
	std::vector<std::string> res12 = te.enumerate(*m12, &tautcat);
	std::vector<std::string> ans12 = {
					"CC=C1C=CC=CN1", 
					"CC=C1C=CCC=N1",
					"CCc1ccccn1"};
	TEST_ASSERT(res12 == ans12);
	
	// special imine tautomer
	std::string smi13 = "C1(=NC=CC=C1)CC";
	std::shared_ptr<ROMol> m13( SmilesToMol(smi13) );
	std::vector<std::string> res13 = te.enumerate(*m13, &tautcat);
	std::vector<std::string> ans13 = {
					"CC=C1C=CC=CN1", 
					"CC=C1C=CCC=N1",
					"CCc1ccccn1"};
	TEST_ASSERT(res13 == ans13);

	// 1,3 aromatic heteroatom H shift
	std::string smi14 = "O=c1cccc[nH]1";
	std::shared_ptr<ROMol> m14( SmilesToMol(smi14) );
	std::vector<std::string> res14 = te.enumerate(*m14, &tautcat);
	std::vector<std::string> ans14 = {"O=c1cccc[nH]1", "Oc1ccccn1"};
	TEST_ASSERT(res14 == ans14);

	// 1,3 aromatic heteroatom H shift
	std::string smi15 = "Oc1ccccn1";
	std::shared_ptr<ROMol> m15( SmilesToMol(smi15) );
	std::vector<std::string> res15 = te.enumerate(*m15, &tautcat);
	std::vector<std::string> ans15 = {"O=c1cccc[nH]1", "Oc1ccccn1"};
	TEST_ASSERT(res15 == ans15);

	// 1,3 aromatic heteroatom H shift
	std::string smi16 = "Oc1ncc[nH]1";
	std::shared_ptr<ROMol> m16( SmilesToMol(smi16) );
	std::vector<std::string> res16 = te.enumerate(*m16, &tautcat);
	std::vector<std::string> ans16 = {"O=c1[nH]cc[nH]1", "Oc1ncc[nH]1"};
	TEST_ASSERT(res16 == ans16);

	// 1,3 heteroatom H shift
	std::string smi17 = "OC(C)=NC";
	std::shared_ptr<ROMol> m17( SmilesToMol(smi17) );
	std::vector<std::string> res17 = te.enumerate(*m17, &tautcat);
	std::vector<std::string> ans17 = {"C=C(O)NC", "CN=C(C)O", "CNC(C)=O"};
	TEST_ASSERT(res17 == ans17);

	// 1,3 heteroatom H shift
	std::string smi18 = "CNC(C)=O";
	std::shared_ptr<ROMol> m18( SmilesToMol(smi18) );
	std::vector<std::string> res18 = te.enumerate(*m18, &tautcat);
	std::vector<std::string> ans18 = {"C=C(O)NC", "CN=C(C)O", "CNC(C)=O"};
	TEST_ASSERT(res18 == ans18);

	// 1,3 heteroatom H shift
	std::string smi19 = "S=C(N)N";
	std::shared_ptr<ROMol> m19( SmilesToMol(smi19) );
	std::vector<std::string> res19 = te.enumerate(*m19, &tautcat);
	std::vector<std::string> ans19 = {"N=C(N)S", "NC(N)=S"};
	TEST_ASSERT(res19 == ans19);

	// 1,3 heteroatom H shift
	std::string smi20 = "SC(N)=N";
	std::shared_ptr<ROMol> m20( SmilesToMol(smi20) );
	std::vector<std::string> res20 = te.enumerate(*m20, &tautcat);
	std::vector<std::string> ans20 = {"N=C(N)S", "NC(N)=S"};
	TEST_ASSERT(res20 == ans20);

	// 1,3 heteroatom H shift
	std::string smi21 = "N=c1[nH]ccn(C)1";
	std::shared_ptr<ROMol> m21( SmilesToMol(smi21) );
	std::vector<std::string> res21 = te.enumerate(*m21, &tautcat);
	std::vector<std::string> ans21 = {"Cn1cc[nH]c1=N", "Cn1ccnc1N"};
	TEST_ASSERT(res21 == ans21);

	// 1,3 heteroatom H shift
	std::string smi22 = "CN=c1[nH]cncc1";
	std::shared_ptr<ROMol> m22( SmilesToMol(smi22) );
	std::vector<std::string> res22 = te.enumerate(*m22, &tautcat);
	std::vector<std::string> ans22 = {"CN=c1cc[nH]cn1", "CN=c1ccnc[nH]1", "CNc1ccncn1"};
	TEST_ASSERT(res22 == ans22);

	// 1,5 aromatic heteroatom H shift
	std::string smi23 = "Oc1cccc2ccncc12";
	std::shared_ptr<ROMol> m23( SmilesToMol(smi23) );
	std::vector<std::string> res23 = te.enumerate(*m23, &tautcat);
	std::vector<std::string> ans23 = {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"};
	TEST_ASSERT(res23 == ans23);

	// 1,5 aromatic heteroatom H shift
	std::string smi24 = "O=c1cccc2cc[nH]cc1-2";
	std::shared_ptr<ROMol> m24( SmilesToMol(smi24) );
	std::vector<std::string> res24 = te.enumerate(*m24, &tautcat);
	std::vector<std::string> ans24 = {"O=c1cccc2cc[nH]cc1-2", "Oc1cccc2ccncc12"};
	TEST_ASSERT(res24 == ans24);

	// 1,5 aromatic heteroatom H shift
	std::string smi25 = "Cc1n[nH]c2ncnn12";
	std::shared_ptr<ROMol> m25( SmilesToMol(smi25) );
	std::vector<std::string> res25 = te.enumerate(*m25, &tautcat);
	std::vector<std::string> ans25 = {
					"C=C1NN=C2N=CNN12", 
					"C=C1NN=C2NC=NN12",
					"C=C1NNc2ncnn21", 
					"Cc1n[nH]c2ncnn12", 
					"Cc1nnc2[nH]cnn12", 
					"Cc1nnc2nc[nH]n12"}; 
	TEST_ASSERT(res25 == ans25);

	// 1,5 aromatic heteroatom H shift
	std::string smi26 = "Cc1nnc2nc[nH]n12";
	std::shared_ptr<ROMol> m26( SmilesToMol(smi26) );
	std::vector<std::string> res26 = te.enumerate(*m26, &tautcat);
	std::vector<std::string> ans26 = {
					"C=C1NN=C2N=CNN12", 
					"C=C1NN=C2NC=NN12",
					"C=C1NNc2ncnn21", 
					"Cc1n[nH]c2ncnn12", 
					"Cc1nnc2[nH]cnn12", 
					"Cc1nnc2nc[nH]n12"}; 
	TEST_ASSERT(res26 == ans26);

	// 1,5 aromatic heteroatom H shift
	std::string smi29 = "Oc1ccncc1";
	std::shared_ptr<ROMol> m29( SmilesToMol(smi29) );
	std::vector<std::string> res29 = te.enumerate(*m29, &tautcat);
	std::vector<std::string> ans29 = {"O=c1cc[nH]cc1", "Oc1ccncc1"};
	TEST_ASSERT(res29 == ans29);

	// 1,5 aromatic heteroatom H shift
	std::string smi27 = "Oc1c(cccc3)c3nc2ccncc12";
	std::shared_ptr<ROMol> m27( SmilesToMol(smi27) );
	std::vector<std::string> res27 = te.enumerate(*m27, &tautcat);
	std::vector<std::string> ans27 = {
					"O=c1c2c[nH]ccc-2nc2ccccc12",
					"O=c1c2ccccc2[nH]c2ccncc12", 
					"Oc1c2ccccc2nc2ccncc12"};
	TEST_ASSERT(res27 == ans27);

	// 1,3 and 1,5 aromatic heteroatom H shift
	std::string smi28 = "Oc1ncncc1";
	std::shared_ptr<ROMol> m28( SmilesToMol(smi28) );
	std::vector<std::string> res28 = te.enumerate(*m28, &tautcat);
	std::vector<std::string> ans28 = {
					"O=c1cc[nH]cn1",
					"O=c1ccnc[nH]1", 
					"Oc1ccncn1"};
	TEST_ASSERT(res28 == ans28);

	// 1,5 aromatic heteroatom H shift
	std::string smi30 = "C2(=C1C(=NC=N1)[NH]C(=N2)N)O";
	std::shared_ptr<ROMol> m30( SmilesToMol(smi30) );
	std::vector<std::string> res30 = te.enumerate(*m30, &tautcat);
	std::vector<std::string> ans30 = {
					"N=c1[nH]c(=O)c2[nH]cnc2[nH]1", 
					"N=c1[nH]c(=O)c2nc[nH]c2[nH]1", 
					"N=c1[nH]c2ncnc-2c(O)[nH]1", 
					"N=c1nc(O)c2[nH]cnc2[nH]1", 
					"N=c1nc(O)c2nc[nH]c2[nH]1", 
					"N=c1nc2[nH]cnc2c(O)[nH]1", 
					"N=c1nc2nc[nH]c2c(O)[nH]1", 
					"Nc1nc(=O)c2[nH]cnc2[nH]1", 
					"Nc1nc(=O)c2nc[nH]c2[nH]1", 
					"Nc1nc(O)c2[nH]cnc2n1",
					"Nc1nc(O)c2nc[nH]c2n1", 
					"Nc1nc(O)c2ncnc-2[nH]1", 
					"Nc1nc2[nH]cnc2c(=O)[nH]1", 
					"Nc1nc2nc[nH]c2c(=O)[nH]1", 
					"Nc1nc2ncnc-2c(O)[nH]1"}; 
	TEST_ASSERT(res30 == ans30);

	// 1,5 aromatic heteroatom H shift
	std::string smi31 = "C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O";
	std::shared_ptr<ROMol> m31( SmilesToMol(smi31) );
	std::vector<std::string> res31 = te.enumerate(*m31, &tautcat);
	std::vector<std::string> ans31 = {
					"N=c1[nH]c(=O)c2[nH]cnc2[nH]1", 
					"N=c1[nH]c(=O)c2nc[nH]c2[nH]1", 
					"N=c1[nH]c2ncnc-2c(O)[nH]1", 
					"N=c1nc(O)c2[nH]cnc2[nH]1", 
					"N=c1nc(O)c2nc[nH]c2[nH]1", 
					"N=c1nc2[nH]cnc2c(O)[nH]1", 
					"N=c1nc2nc[nH]c2c(O)[nH]1", 
					"Nc1nc(=O)c2[nH]cnc2[nH]1", 
					"Nc1nc(=O)c2nc[nH]c2[nH]1", 
					"Nc1nc(O)c2[nH]cnc2n1",
					"Nc1nc(O)c2nc[nH]c2n1", 
					"Nc1nc(O)c2ncnc-2[nH]1", 
					"Nc1nc2[nH]cnc2c(=O)[nH]1", 
					"Nc1nc2nc[nH]c2c(=O)[nH]1", 
					"Nc1nc2ncnc-2c(O)[nH]1"}; 
	TEST_ASSERT(res31 == ans31);

	// 1,5 aromatic heteroatom H shift
	std::string smi32 = "Oc1n(C)ncc1";
	std::shared_ptr<ROMol> m32( SmilesToMol(smi32) );
	std::vector<std::string> res32 = te.enumerate(*m32, &tautcat);
	std::vector<std::string> ans32 = {
					"CN1N=CCC1=O",
					"Cn1[nH]ccc1=O", 
					"Cn1nccc1O"};
	TEST_ASSERT(res32 == ans32);

	// 1,5 aromatic heteroatom H shift
	std::string smi33 = "O=c1nc2[nH]ccn2cc1";
	std::shared_ptr<ROMol> m33( SmilesToMol(smi33) );
	std::vector<std::string> res33 = te.enumerate(*m33, &tautcat);
	std::vector<std::string> ans33 = {
					"O=c1ccn2cc[nH]c2n1",
					"O=c1ccn2ccnc2[nH]1", 
					"Oc1ccn2ccnc2n1"};
	TEST_ASSERT(res33 == ans33);

	// 1,5 aromatic heteroatom H shift
	std::string smi34 = "N=c1nc[nH]cc1";
	std::shared_ptr<ROMol> m34( SmilesToMol(smi34) );
	std::vector<std::string> res34 = te.enumerate(*m34, &tautcat);
	std::vector<std::string> ans34 = {
					"N=c1cc[nH]cn1",
					"N=c1ccnc[nH]1", 
					"Nc1ccncn1"};
	TEST_ASSERT(res34 == ans34);

	// 1,5 aromatic heteroatom H shift
	std::string smi35 = "N=c(c1)ccn2cc[nH]c12";
	std::shared_ptr<ROMol> m35( SmilesToMol(smi35) );
	std::vector<std::string> res35 = te.enumerate(*m35, &tautcat);
	std::vector<std::string> ans35 = {
					"N=c1ccn2cc[nH]c2c1",
					"Nc1ccn2ccnc2c1"};
	TEST_ASSERT(res35 == ans35);

	// 1,5 aromatic heteroatom H shift
	std::string smi36 = "CN=c1nc[nH]cc1";
	std::shared_ptr<ROMol> m36( SmilesToMol(smi36) );
	std::vector<std::string> res36 = te.enumerate(*m36, &tautcat);
	std::vector<std::string> ans36 = {
					"CN=c1cc[nH]cn1",
					"CN=c1ccnc[nH]1", 
					"CNc1ccncn1"};
	TEST_ASSERT(res36 == ans36);

	// 1,7 aromatic heteroatom H shift
	std::string smi37 = "c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1";
	std::shared_ptr<ROMol> m37( SmilesToMol(smi37) );
	std::vector<std::string> res37 = te.enumerate(*m37, &tautcat);
	std::vector<std::string> ans37 = {
					"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1",
				 	"c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2", 
					"c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2"};
	TEST_ASSERT(res37 == ans37);

	// 1,7 aromatic heteroatom H shift
	std::string smi38 = "c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2";
	std::shared_ptr<ROMol> m38( SmilesToMol(smi38) );
	std::vector<std::string> res38 = te.enumerate(*m38, &tautcat);
	std::vector<std::string> ans38 = {
					"c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1", 
					"c1ccc2c(c1)=NC(c1nc3ccccc3[nH]1)N=2", 
					"c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2"};
	TEST_ASSERT(res38 == ans38);

	// 1,9 aromatic heteroatom H shift
	std::string smi39 = "CNc1ccnc2ncnn21";
	std::shared_ptr<ROMol> m39( SmilesToMol(smi39) );
	std::vector<std::string> res39 = te.enumerate(*m39, &tautcat);
	std::vector<std::string> ans39 = {
					"CN=c1cc[nH]c2ncnn12", 
				 	"CN=c1ccnc2[nH]cnn12", 
					"CN=c1ccnc2nc[nH]n12",
					"CNc1ccnc2ncnn12"};
	TEST_ASSERT(res39 == ans39);

	// 1,9 aromatic heteroatom H shift
	std::string smi40 = "CN=c1ccnc2nc[nH]n21";
	std::shared_ptr<ROMol> m40( SmilesToMol(smi40) );
	std::vector<std::string> res40 = te.enumerate(*m40, &tautcat);
	std::vector<std::string> ans40 = {
					"CN=c1cc[nH]c2ncnn12", 
					"CN=c1ccnc2[nH]cnn12", 
					"CN=c1ccnc2nc[nH]n12", 
					"CNc1ccnc2ncnn12"};
	TEST_ASSERT(res40 == ans40);

	// 1,11 aromatic heteroatom H shift
	std::string smi41 = "Nc1ccc(C=C2C=CC(=O)C=C2)cc1";
	std::shared_ptr<ROMol> m41( SmilesToMol(smi41) );
	std::vector<std::string> res41 = te.enumerate(*m41, &tautcat);
	std::vector<std::string> ans41 = {
					"N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1", 
					"N=C1C=CC(=Cc2ccc(O)cc2)C=C1", 
					"N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1",
					"Nc1ccc(C=C2C=CC(=O)C=C2)cc1"};
      	TEST_ASSERT(res41 == ans41);

	// 1,11 aromatic heteroatom H shift
	std::string smi42 = "N=C1C=CC(=Cc2ccc(O)cc2)C=C1";
	std::shared_ptr<ROMol> m42( SmilesToMol(smi42) );
	std::vector<std::string> res42 = te.enumerate(*m42, &tautcat);
	std::vector<std::string> ans42 = {
					"N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1", 
					"N=C1C=CC(=Cc2ccc(O)cc2)C=C1", 
					"N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1",
					"Nc1ccc(C=C2C=CC(=O)C=C2)cc1"}; 
	TEST_ASSERT(res42 == ans42);

	// heterocyclic tautomer
	std::string smi43 = "n1ccc2ccc[nH]c12";
	std::shared_ptr<ROMol> m43( SmilesToMol(smi43) );
	std::vector<std::string> res43 = te.enumerate(*m43, &tautcat);
	std::vector<std::string> ans43 = {"c1c[nH]c2nccc-2c1", "c1cnc2[nH]ccc2c1"}; 
	TEST_ASSERT(res43 == ans43);

	// heterocyclic tautomer
	std::string smi44 = "c1cc(=O)[nH]c2nccn12";
	std::shared_ptr<ROMol> m44( SmilesToMol(smi44) );
	std::vector<std::string> res44 = te.enumerate(*m44, &tautcat);
	std::vector<std::string> ans44 = {
					"O=c1ccn2cc[nH]c2n1", 
					"O=c1ccn2ccnc2[nH]1", 
					"Oc1ccn2ccnc2n1"};
	TEST_ASSERT(res44 == ans44);

	// heterocyclic tautomer
	std::string smi45 = "c1cnc2c[nH]ccc12";
	std::shared_ptr<ROMol> m45( SmilesToMol(smi45) );
	std::vector<std::string> res45 = te.enumerate(*m45, &tautcat);
	std::vector<std::string> ans45 = {"c1cc2cc[nH]c2cn1", "c1cc2cc[nH]cc-2n1"}; 
	TEST_ASSERT(res45 == ans45);

	// heterocyclic tautomer
	std::string smi46 = "n1ccc2c[nH]ccc12";
	std::shared_ptr<ROMol> m46( SmilesToMol(smi46) );
	std::vector<std::string> res46 = te.enumerate(*m46, &tautcat);
	std::vector<std::string> ans46 = {
					"c1cc2[nH]ccc2cn1", 
					"c1cc2c[nH]ccc-2n1"};
	TEST_ASSERT(res46 == ans46);

	// heterocyclic tautomer
	std::string smi47 = "c1cnc2ccc[nH]c12";
	std::shared_ptr<ROMol> m47( SmilesToMol(smi47) );
	std::vector<std::string> res47 = te.enumerate(*m47, &tautcat);
	std::vector<std::string> ans47 = {
					"c1c[nH]c2ccnc-2c1", 
					"c1cnc2cc[nH]c2c1"};
	TEST_ASSERT(res47 == ans47);

	// furanone tautomer
	std::string smi48 = "C1=CC=C(O1)O";
	std::shared_ptr<ROMol> m48( SmilesToMol(smi48) );
	std::vector<std::string> res48 = te.enumerate(*m48, &tautcat);
	std::vector<std::string> ans48 = {
					"O=C1CC=CO1", 
					"Oc1ccco1"};
	TEST_ASSERT(res48 == ans48);

	// furanone tautomer
	std::string smi49 = "O=C1CC=CO1";
	std::shared_ptr<ROMol> m49( SmilesToMol(smi49) );
	std::vector<std::string> res49 = te.enumerate(*m49, &tautcat);
	std::vector<std::string> ans49 = {
					"O=C1CC=CO1", 
					"Oc1ccco1"};
	TEST_ASSERT(res49 == ans49);

	// keten/ynol tautomer
	std::string smi50 = "CC=C=O";
	std::shared_ptr<ROMol> m50( SmilesToMol(smi50) );
	std::vector<std::string> res50 = te.enumerate(*m50, &tautcat);
	std::vector<std::string> ans50 = {"CC#CO", "CC=C=O"};
//	TEST_ASSERT(res50 == ans50);

//	// keten/ynol tautomer
//	std::string smi51 = "CC#CO";
//	std::shared_ptr<ROMol> m51( SmilesToMol(smi51) );
//	std::vector<std::string> res51 = te.enumerate(*m51, &tautcat);
//	std::vector<std::string> ans51 = {"CC#CO", "CC=C=O"};
//	TEST_ASSERT(res51 == ans51);




	//	smi2 = "c1(ccccc1)/C=C(/O)\\C";
//	std::shared_ptr<ROMol> m2( SmilesToMol(smi2) );
//	ROMol* res2 = te.enumerate(*m2, &tautcat);

}

void ttt() {
	//ROMol* ref = SmartsToMol("[O,S,Se,Te;X2!H0]-[CH0]=[C]-[C]=[C,N]");
	ROMol* ref = SmilesToMol("C=CC=CO");
	ROMol* mol = SmilesToMol("CC(O)=Cc1ccccc1");
//	ROMol* mol = SmilesToMol("CC(O)=CC1=CC=CC=C1");
	std::vector< MatchVectType > matches;
	unsigned int matched = SubstructMatch( *mol, *ref, matches );
	if (!matched) {
		std::cout << "No match" << std::endl;
	} else {
		std::cout << "Match!" << std::endl;
	}


}

int main() {
	test1();
//	ttt();
	return 0;
}
