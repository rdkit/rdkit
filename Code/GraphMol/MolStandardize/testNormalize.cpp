//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogUtils.h>
#include "Normalize.h"
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
#include <memory>

using namespace RDKit;
using namespace MolStandardize;

void test1() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test1" << std::endl;
  std::string smi1, smi2, smi3, smi4, smi5, smi6, smi7;

  Normalizer normalizer;
  auto normalize = [&normalizer](const std::string & smiles) -> std::string {
    std::unique_ptr<ROMol> molecule(SmilesToMol(smiles));
    std::unique_ptr<ROMol> normalized(normalizer.normalize(*molecule));
    return MolToSmiles(*normalized);
  };

  // Test sulfoxide normalization.
  TEST_ASSERT(normalize("CS(C)=O") == "C[S+](C)[O-]");

  // Test sulfone
  TEST_ASSERT(normalize("C[S+2]([O-])([O-])O") == "CS(=O)(=O)O");

  // Test 1,3-separated charges are recombined.
  TEST_ASSERT(normalize("CC([O-])=[N+](C)C") == "CC(=O)N(C)C");

  // Test 1,3-separated charges are recombined.
  TEST_ASSERT(normalize("C[n+]1ccccc1[O-]") == "Cn1ccccc1=O");

  // Test a case where 1,3-separated charges should not be recombined.
  TEST_ASSERT(
    normalize("CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]")
    == "CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]");

  // Test a case where 1,3-separated charges should not be recombined.
  TEST_ASSERT(
    normalize("[O-][n+]1cccc[n+]1[O-]") == "[O-][n+]1cccc[n+]1[O-]");

  // Test 1,5-separated charges are recombined.
  TEST_ASSERT(normalize(R"(C[N+](C)=C\C=C\[O-])") == "CN(C)C=CC=O");

  // Test a case where 1,5-separated charges should not be recombined.
  TEST_ASSERT(
    normalize("C[N+]1=C2C=[N+]([O-])C=CN2CCC1")
    == "C[N+]1=C2C=[N+]([O-])C=CN2CCC1");

  // Failed on 1k normalize test sanitizeMol step
  TEST_ASSERT(
    normalize("O=c1cc([O-])[n+](C2OC(CO)C(O)C2O)c2sccn12")
    == "O=c1cc([O-])[n+](C2OC(CO)C(O)C2O)c2sccn12");

  // Test normalization of 1,3-conjugated cations
  TEST_ASSERT(normalize("C[N+](C)=CN") == "CN(C)C=[NH2+]");
  // verify it doesn't apply to oximes
  TEST_ASSERT(normalize("C[N+](C)=NO") == "C[N+](C)=NO");
  // verify it doesn't affect diazo groups
  TEST_ASSERT(normalize("[N-]=[N+]=CN") == "[N-]=[N+]=CN");
  // but still applies if the adjacent heteroatom is neutral
  TEST_ASSERT(normalize("C[N+](N)=CN") == "CN(N)C=[NH2+]");

  // Test normalization of 1,5-conjugated cations
  TEST_ASSERT(normalize("C[N+](C)=CC=CN") == "CN(C)C=CC=[NH2+]");
  // verify it doesn't apply to oximes
  TEST_ASSERT(normalize("C[N+](C)=CC=NO") == "C[N+](C)=CC=NO");
  // verify it doesn't affect diazo groups
  TEST_ASSERT(normalize("[N-]=[N+]=CC=CN") == "[N-]=[N+]=CC=CN");
  // but still applies if the adjacent heteroatom is neutral
  TEST_ASSERT(normalize("C[N+](N)=CC=CN") == "CN(N)C=CC=[NH2+]");

  // Test normalization of 1,3-conjugated cations
  TEST_ASSERT(normalize("C[N+](C)=c1cccc[nH]1") == "CN(C)c1cccc[nH+]1");
  // verify it doesn't affect diazo groups
  TEST_ASSERT(normalize("[N-]=[N+]=c1cccc[nH]1") == "[N-]=[N+]=c1cccc[nH]1");
  // but still applies if the adjacent heteroatom is neutral
  TEST_ASSERT(normalize("C[N+](N)=c1cccc[nH]1") == "CN(N)c1cccc[nH+]1");

  // Test normalization of 1,5-conjugated cations
  TEST_ASSERT(normalize("C[N+](C)=c1cc[nH]cc1") == "CN(C)c1cc[nH+]cc1");
  // verify it doesn't affect diazo groups
  TEST_ASSERT(normalize("[N-]=[N+]=c1cc[nH]cc1") == "[N-]=[N+]=c1cc[nH]cc1");
  // but still applies if the adjacent heteroatom is neutral
  TEST_ASSERT(normalize("C[N+](N)=c1cc[nH]cc1") == "CN(N)c1cc[nH+]cc1");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n test2" << std::endl;
  {
    // initialization from string:
    std::string tfdata = R"DATA(//	Name	SMIRKS
Nitro to N+(O-)=O	[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]
Sulfone to S(=O)(=O)	[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])
Pyridine oxide to n+O-	[n:1]=[O:2]>>[n+:1][O-:2]
)DATA";
    std::stringstream sstr(tfdata);
    Normalizer nn(sstr, 10);
    bool debugParse = false;
    bool sanitize = false;
    std::unique_ptr<ROMol> imol(
        SmilesToMol("O=N(=O)CCN=N#N", debugParse, sanitize));
    std::unique_ptr<ROMol> m2(nn.normalize(*imol));
    TEST_ASSERT(MolToSmiles(*m2) == "N#N=NCC[N+](=O)[O-]");
  }
  {
    // initialization from vector:
    std::vector<std::pair<std::string, std::string>> tfdata = {
        {"Nitro to N+(O-)=O",
         "[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:"
         "3]"},
        {"Sulfone to S(=O)(=O)",
         "[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])"},
        {"Pyridine oxide to n+O-", "[n:1]=[O:2]>>[n+:1][O-:2]"}};
    Normalizer nn(tfdata, 10);
    bool debugParse = false;
    bool sanitize = false;
    std::unique_ptr<ROMol> imol(
        SmilesToMol("O=N(=O)CCN=N#N", debugParse, sanitize));
    std::unique_ptr<ROMol> m2(nn.normalize(*imol));
    TEST_ASSERT(MolToSmiles(*m2) == "N#N=NCC[N+](=O)[O-]");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub2414() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing github #2414: "
                          "combinatorial explosion in Normalizer"
                       << std::endl;
  // pubchem sid 7638352
  std::string molb = R"MOL(7638352
  -OEChem-01301907472D

143151  0     0  0  0  0  0  0999 V2000
    1.6643   -1.2092    0.0000 Rh  0  0  0  0  0  0  0  0  0  0  0  0
    1.3553   -0.2582    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    2.4646   -6.2630    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7350    2.0001    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3172    0.8272    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9127   -3.0768    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.2298   -4.2734    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    6.5607    0.3263    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    3.1090    3.8378    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    3.3001    1.7525    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    6.7114   -2.6538    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.7500    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.6154   -0.9002    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    1.9734   -2.1603    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    2.8297   -7.9562    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
   -2.6054    3.4976    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4723    2.1178    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
   -5.3110   -2.0546    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
   -0.6407   -5.7709    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    8.2539    0.6914    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    3.2849    5.5609    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    3.9999    3.3368    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    8.4344   -2.8298    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    3.1359   -7.0042    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6025    2.4976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5222    1.8059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.9069   -2.9693    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6377   -4.7709    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3019    0.9976    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.6995    4.7501    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1927    2.7467    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.6237   -2.2444    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.7234   -6.9343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2058   -5.5917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2376    2.8676    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2324    1.1326    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2959    0.6221    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3384    1.0322    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0202   -4.0710    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8053   -2.0826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7272   -5.1409    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2677   -3.4059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2321   -0.4149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8894    1.0674    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0213    4.2472    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1966    3.4284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.2943    1.8599    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3059    1.6450    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.1208   -3.5662    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3019   -1.7415    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9511   -1.4410    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5878   -2.5590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6154   -1.9002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5665   -1.2092    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6440    0.5156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7824   -1.5725    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3891   -3.1889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4326   -3.8246    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1561   -0.4622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5820   -2.4516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7479   -2.3977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3077   -0.5379    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2346    1.4280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6750   -0.5782    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3015   -2.7795    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1039   -4.5658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6923   -2.1123    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1784   -3.4714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4829   -2.3977    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7715   -2.1880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6383    0.4082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6948   -1.9819    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2817   -4.1831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4538   -4.0297    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7933   -5.5219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1121   -0.1516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1729   -3.2648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7479   -3.4029    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2637   -0.8485    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8255    2.2412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4882    0.0126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1148   -3.3703    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6483   -1.8017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7692   -4.2846    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4829   -3.4029    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7275   -2.4986    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2291    1.2214    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5080   -1.3910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0950   -4.7740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1432   -4.9857    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8113   -5.7366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    2.0104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8631   -0.8197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7695   -4.1855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6154   -3.9106    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.4785   -1.8305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8257    2.1420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4088   -0.3908    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0156   -4.3706    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4408   -2.3568    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1746   -0.2934    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2234   -2.1603    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7328   -0.0038    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.1461    0.7069    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9874   -1.5375    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8826   -1.8964    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.9414    1.1733    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4940    1.6084    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0752    0.1312    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.4819   -1.0388    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8151   -4.1995    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7328   -0.0038    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4829   -3.0901    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8163   -3.5745    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.3482   -1.8964    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4693   -3.3571    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.0437   -0.5060    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7979   -2.9766    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6820   -4.8926    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7153   -2.7275    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3856    2.3732    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3894   -2.4730    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0527   -5.8812    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.3504   -3.9003    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.9326   -3.4774    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.2233    1.1139    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4203   -1.8005    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.9069   -6.5139    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5696   -5.3445    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2706   -7.4010    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    3.0104    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8141   -0.5107    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3573   -4.9945    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6154   -4.9106    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.4295   -2.1395    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.4135    2.9510    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.2178    0.1970    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.8246   -4.9584    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1 12  1  0  0  0  0
  1 13  1  0  0  0  0
  1 14  1  0  0  0  0
  3 24  1  0  0  0  0
  3 33  2  0  0  0  0
  3 34  2  0  0  0  0
  3 78  1  0  0  0  0
  4 25  1  0  0  0  0
  4 35  2  0  0  0  0
  4 36  2  0  0  0  0
  4 79  1  0  0  0  0
  5 26  1  0  0  0  0
  5 37  2  0  0  0  0
  5 38  2  0  0  0  0
  5 80  1  0  0  0  0
  6 27  1  0  0  0  0
  6 39  2  0  0  0  0
  6 40  2  0  0  0  0
  6 81  1  0  0  0  0
  7 28  1  0  0  0  0
  7 41  2  0  0  0  0
  7 42  2  0  0  0  0
  7 82  1  0  0  0  0
  8 29  1  0  0  0  0
  8 43  2  0  0  0  0
  8 44  2  0  0  0  0
  8 83  1  0  0  0  0
  9 30  1  0  0  0  0
  9 45  2  0  0  0  0
  9 46  2  0  0  0  0
  9 84  1  0  0  0  0
 10 31  1  0  0  0  0
 10 47  2  0  0  0  0
 10 48  2  0  0  0  0
 10 85  1  0  0  0  0
 11 32  1  0  0  0  0
 11 49  2  0  0  0  0
 11 50  2  0  0  0  0
 11 86  1  0  0  0  0
 12 51  1  0  0  0  0
 12 52  1  0  0  0  0
 12 53  1  0  0  0  0
 12105  1  0  0  0  0
 13 54  1  0  0  0  0
 13 55  1  0  0  0  0
 13 56  1  0  0  0  0
 13106  1  0  0  0  0
 14 57  1  0  0  0  0
 14 58  1  0  0  0  0
 14 59  1  0  0  0  0
 14107  1  0  0  0  0
 15 24  1  0  0  0  0
 16 25  1  0  0  0  0
 17 26  1  0  0  0  0
 18 27  1  0  0  0  0
 19 28  1  0  0  0  0
 20 29  1  0  0  0  0
 21 30  1  0  0  0  0
 22 31  1  0  0  0  0
 23 32  1  0  0  0  0
 51 60  2  0  0  0  0
 51 69  1  0  0  0  0
 52 61  2  0  0  0  0
 52 70  1  0  0  0  0
 53 62  2  0  0  0  0
 53 71  1  0  0  0  0
 54 63  2  0  0  0  0
 54 72  1  0  0  0  0
 55 64  2  0  0  0  0
 55 73  1  0  0  0  0
 56 65  2  0  0  0  0
 56 74  1  0  0  0  0
 57 66  2  0  0  0  0
 57 75  1  0  0  0  0
 58 67  2  0  0  0  0
 58 76  1  0  0  0  0
 59 68  1  0  0  0  0
 59 77  2  0  0  0  0
 60 79  1  0  0  0  0
 60108  1  0  0  0  0
 61 80  1  0  0  0  0
 61109  1  0  0  0  0
 62 81  1  0  0  0  0
 62110  1  0  0  0  0
 63 82  1  0  0  0  0
 63111  1  0  0  0  0
 64 83  1  0  0  0  0
 64112  1  0  0  0  0
 65 84  1  0  0  0  0
 65113  1  0  0  0  0
 66 85  1  0  0  0  0
 66114  1  0  0  0  0
 67 86  1  0  0  0  0
 67115  1  0  0  0  0
 68 78  2  0  0  0  0
 68116  1  0  0  0  0
 69 87  2  0  0  0  0
 69117  1  0  0  0  0
 70 88  2  0  0  0  0
 70118  1  0  0  0  0
 71 89  2  0  0  0  0
 71119  1  0  0  0  0
 72 90  2  0  0  0  0
 72120  1  0  0  0  0
 73 91  2  0  0  0  0
 73121  1  0  0  0  0
 74 92  2  0  0  0  0
 74122  1  0  0  0  0
 75 93  2  0  0  0  0
 75123  1  0  0  0  0
 76 94  2  0  0  0  0
 76124  1  0  0  0  0
 77 95  1  0  0  0  0
 77125  1  0  0  0  0
 78 96  1  0  0  0  0
 79 97  2  0  0  0  0
 80 98  2  0  0  0  0
 81 99  2  0  0  0  0
 82100  2  0  0  0  0
 83101  2  0  0  0  0
 84102  2  0  0  0  0
 85103  2  0  0  0  0
 86104  2  0  0  0  0
 87 97  1  0  0  0  0
 87126  1  0  0  0  0
 88 98  1  0  0  0  0
 88127  1  0  0  0  0
 89 99  1  0  0  0  0
 89128  1  0  0  0  0
 90100  1  0  0  0  0
 90129  1  0  0  0  0
 91101  1  0  0  0  0
 91130  1  0  0  0  0
 92102  1  0  0  0  0
 92131  1  0  0  0  0
 93103  1  0  0  0  0
 93132  1  0  0  0  0
 94104  1  0  0  0  0
 94133  1  0  0  0  0
 95 96  2  0  0  0  0
 95134  1  0  0  0  0
 96135  1  0  0  0  0
 97136  1  0  0  0  0
 98137  1  0  0  0  0
 99138  1  0  0  0  0
100139  1  0  0  0  0
101140  1  0  0  0  0
102141  1  0  0  0  0
103142  1  0  0  0  0
104143  1  0  0  0  0
M  END)MOL";

  std::string tfdata = R"DATA(//	Name	SMIRKS
Alkaline oxide to ions	[Li,Na,K;+0:1]-[O+0:2]>>([*+1:1].[O-:2])
)DATA";
  std::stringstream sstr(tfdata);
  Normalizer nn(sstr, 10);
  bool sanitize = false;
  bool removeHs = false;
  std::unique_ptr<ROMol> imol(MolBlockToMol(molb, sanitize, removeHs));
  std::unique_ptr<ROMol> m2(nn.normalize(*imol));
  TEST_ASSERT(m2);
  auto p = "[Na]-O"_smarts;
  TEST_ASSERT(p);
  TEST_ASSERT(SubstructMatch(*imol, *p).size() == 9);
  TEST_ASSERT(SubstructMatch(*m2, *p).size() == 0);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testNormalizeMultipleAltSmarts() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing that multiple SMARTS "
         "matching the same group work"
      << std::endl;
  std::string azideNormalizations = R"DATA(// Name	SMARTS
Azide to N=N+=N-	[N:2]=[N:3]#[N:4]>>[N:2]=[N+:3]=[N-:4]
Broken azide to N=N+=N-	[N:2]=[N:3]=[N:4]>>[NH0:2]=[NH0+:3]=[NH0-:4])DATA";
  std::stringstream azideNormalizationsStream(azideNormalizations);
  std::stringstream captureLog;
  auto ostate = rdInfoLog->df_enabled;
  rdInfoLog->df_enabled = true;
  rdInfoLog->SetTee(captureLog);
  Normalizer nn(azideNormalizationsStream, 200);
  const std::string brokenAzideSmi = "CN=[N+]=[NH-]";
  const int debugParse = 0;
  const bool sanitize = false;
  ROMOL_SPTR brokenAzide(SmilesToMol(brokenAzideSmi, debugParse, sanitize));
  ROMOL_SPTR normalizedAzide(nn.normalize(*brokenAzide));
  rdInfoLog->ClearTee();
  std::string line;
  unsigned int count = 0;
  while (std::getline(captureLog, line)) {
    if (line.find("Rule applied: BrokenazidetoN=N+=N-") != std::string::npos) {
      ++count;
    }
  }
  rdInfoLog->df_enabled = ostate;
  TEST_ASSERT(count == 1);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testGithub3460() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Github #3460: "
                          "Normalization rule incorrectly matches sulfones"
                       << std::endl;
  std::stringstream captureLog;
  auto ostate = rdInfoLog->df_enabled;
  rdInfoLog->df_enabled = true;
  rdInfoLog->SetTee(captureLog);
  Normalizer nn;
  auto mol = "[O-][S+]1Nc2c(Cl)cc(Cl)c3c(Cl)cc(Cl)c(c23)N1"_smiles;
  TEST_ASSERT(mol);
  ROMOL_SPTR normalized(nn.normalize(*mol));
  rdInfoLog->ClearTee();
  rdInfoLog->df_enabled = ostate;
  auto logged = captureLog.str();
  TEST_ASSERT(logged.find("Running Normalizer") != std::string::npos);
  TEST_ASSERT(logged.find("Rule applied: C/S+NtoC/S=N+") == std::string::npos);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testEmptyMol() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Test that Normalizer "
                          "does not crash on an empty mol"
                       << std::endl;
  Normalizer nn;
  std::unique_ptr<ROMol> emptyMol(new ROMol());
  std::unique_ptr<ROMol> normalized(nn.normalize(*emptyMol));
  TEST_ASSERT(!normalized->getNumAtoms());
}

void testGithub4281() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing Github #3460: "
         "Normalization rule incorrectly matches [N+]=O rather than n=O"
      << std::endl;
  auto mol = "Cn1c(=O)c2nc[nH][n+](=O)c2n(C)c1=O"_smiles;
  std::stringstream captureLog;
  rdInfoLog->SetTee(captureLog);
  auto ostate = rdInfoLog->df_enabled;
  rdInfoLog->df_enabled = true;
  Normalizer nn;
  std::unique_ptr<ROMol> normalized(nn.normalize(*mol));
  rdInfoLog->ClearTee();
  rdInfoLog->df_enabled = ostate;
  auto logged = captureLog.str();
  TEST_ASSERT(logged.find("FAILED sanitizeMol") == std::string::npos);
}

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
#endif
  testGithub2414();
  testNormalizeMultipleAltSmarts();
  testGithub3460();
  testEmptyMol();
  testGithub4281();
  return 0;
}
