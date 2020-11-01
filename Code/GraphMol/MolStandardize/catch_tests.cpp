//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/MolStandardize/Charge.h>

#include <iostream>
#include <fstream>

using namespace RDKit;

TEST_CASE("SKIP_IF_ALL_MATCH") {
  auto m = "[Na+].[Cl-]"_smiles;
  REQUIRE(m);

  SECTION("default") {
    MolStandardize::FragmentRemover fragRemover;
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[Na+]");
  }
  SECTION("don't remove all") {
    MolStandardize::FragmentRemover fragRemover("", true, true);
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[Cl-].[Na+]");
  }
  SECTION("feel free to remove everything") {
    MolStandardize::FragmentRemover fragRemover("", false, false);
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(outm->getNumAtoms() == 0);
  }
  SECTION("don't remove all 2") {
    MolStandardize::FragmentRemover fragRemover("", true, true);
    auto m = "[Na+].[Cl-].[Na+].[Cl-]"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[Cl-].[Cl-].[Na+].[Na+]");
  }
}

TEST_CASE("symmetry in the uncharger", "[uncharger]") {
  SECTION("case 1") {
    auto m = "C[N+](C)(C)CC(C(=O)[O-])CC(=O)[O-]"_smiles;
    REQUIRE(m);
    {
      bool canonicalOrdering = false;
      MolStandardize::Uncharger uncharger(canonicalOrdering);
      std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)[O-])C(=O)O");
    }
    {
      bool canonicalOrdering = true;
      MolStandardize::Uncharger uncharger(canonicalOrdering);
      std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)O)C(=O)[O-]");
    }
    {
      MolStandardize::CleanupParameters params;
      std::unique_ptr<ROMol> outm(MolStandardize::chargeParent(*m, params));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)O)C(=O)[O-]");
    }
    {
      MolStandardize::CleanupParameters params;
      params.doCanonical = false;
      std::unique_ptr<ROMol> outm(MolStandardize::chargeParent(*m, params));
      REQUIRE(outm);
      CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CC(=O)[O-])C(=O)O");
    }
  }
}

TEST_CASE("uncharger bug with duplicates", "[uncharger]") {
  SECTION("case 1") {
    auto m = "[NH3+]CC([O-])C[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "NCC(O)CO");
  }
  SECTION("case 2") {
    auto m = "CC([O-])C[O-].[Na+]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CC([O-])CO.[Na+]");
  }
  SECTION("acids + others 1, github #2392") {
    auto m = "C[N+](C)(C)CC(C[O-])CC(=O)[O-]"_smiles;
    REQUIRE(m);
    bool doCanonical = false;
    MolStandardize::Uncharger uncharger(doCanonical);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CO)CC(=O)[O-]");
  }
  SECTION("acids + others 2, github #2392") {
    auto m = "C[N+](C)(C)CC(CC(=O)[O-])C[O-]"_smiles;
    REQUIRE(m);
    bool doCanonical = false;
    MolStandardize::Uncharger uncharger(doCanonical);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)CC(CO)CC(=O)[O-]");
  }
}

TEST_CASE(
    "github #2411: MolStandardize: FragmentRemover should not sanitize "
    "[fragments]") {
  SECTION("demo") {
    std::string smi = "CN(C)(C)C.Cl";
    bool debugParse = false;
    bool sanitize = false;
    std::unique_ptr<ROMol> m(SmilesToMol(smi, debugParse, sanitize));
    REQUIRE(m);

    MolStandardize::FragmentRemover fragRemover;
    std::unique_ptr<ROMol> outm(fragRemover.remove(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CN(C)(C)C");
  }
}

TEST_CASE(
    "github #2452: incorrectly removing charge from boron anions"
    "[fragments][uncharger]") {
  SECTION("demo") {
    auto m = "C[B-](C)(C)C"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;

    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(MolToSmiles(*outm) == "C[B-](C)(C)C");
  }
  SECTION("should be removed") {
    auto m = "C[BH-](C)(C)"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;

    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == 0);
    CHECK(MolToSmiles(*outm) == "CB(C)C");
  }
}

TEST_CASE("github #2602: Uncharger ignores dications", "[uncharger]") {
  SECTION("demo") {
    auto m = "[O-]CCC[O-].[Ca+2]"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 2);
    CHECK(outm->getAtomWithIdx(0)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(4)->getFormalCharge() == -1);
    CHECK(MolToSmiles(*outm) == "[Ca+2].[O-]CCC[O-]");
  }
}

TEST_CASE(
    "github #2605: Uncharger incorrectly neutralizes cations when "
    "non-neutralizable anions are present.",
    "[uncharger]") {
  SECTION("demo") {
    auto m = "F[B-](F)(F)F.[NH3+]CCC"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "CCC[NH3+].F[B-](F)(F)F");
  }
  SECTION("multiple positively charged sites") {
    auto m = "F[B-](F)(F)F.[NH3+]CC=C[NH3+]"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(9)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "F[B-](F)(F)F.NCC=C[NH3+]");
  }
  SECTION("make sure we don't go too far") {
    auto m = "F[B-](F)(F)F.[OH3+2]CCC"_smiles;  // totally bogus structure
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "CCC[OH2+].F[B-](F)(F)F");
  }
}

TEST_CASE("github #2610: Uncharger incorrectly modifying a zwitterion.",
          "[uncharger]") {
  SECTION("demo") {
    auto m = "C1=CC=CC[NH+]1-[O-]"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(outm->getAtomWithIdx(6)->getFormalCharge() == -1);
    CHECK(MolToSmiles(*outm) == "[O-][NH+]1C=CC=CC1");
  }
}

TEST_CASE("problems with ringInfo initialization", "[normalizer]") {
  std::string tfs =
      R"TXT(Bad amide tautomer1	[C:1]([OH1;D1:2])=;!@[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]
Bad amide tautomer2	[C:1]([OH1;D1:2])=;!@[NH0:3]>>[C:1](=[OH0:2])-[NH1:3])TXT";
  std::stringstream iss(tfs);
  MolStandardize::Normalizer nrml(iss, 20);
  SECTION("example1") {
    auto m = "Cl.Cl.OC(=N)NCCCCCCCCCCCCNC(O)=N"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(nrml.normalize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "Cl.Cl.NC(=O)NCCCCCCCCCCCCNC(N)=O");
  }
}

TEST_CASE("segfault in normalizer", "[normalizer]") {
  std::string tfs =
      R"TXT(Bad amide tautomer1	[C:1]([OH1;D1:2])=;!@[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]
Bad amide tautomer2	[C:1]([OH1;D1:2])=;!@[NH0:3]>>[C:1](=[OH0:2])-[NH1:3])TXT";
  std::stringstream iss(tfs);
  MolStandardize::Normalizer nrml(iss, 20);
  SECTION("example1") {
    std::string molblock = R"CTAB(molblock = """
  SciTegic12221702182D

 47 51  0  0  0  0            999 V2000
    0.2962    6.2611    0.0000 C   0  0
   -3.9004    4.4820    0.0000 C   0  0
    1.4195    5.2670    0.0000 C   0  0
   -3.8201   -7.4431    0.0000 C   0  0
   -4.9433   -6.4490    0.0000 C   0  0
   -2.3975   -6.9674    0.0000 C   0  0
    3.5921   -3.5947    0.0000 C   0  0
   -3.1475    2.3700    0.0000 C   0  0
    2.1695   -4.0705    0.0000 C   0  0
   -2.0242    1.3759    0.0000 C   0  0
   -4.6440   -4.9792    0.0000 C   0  0
    2.7681   -1.1308    0.0000 C   0  0
   -5.8626    1.1332    0.0000 C   0  0
    3.0674    0.3391    0.0000 C   0  0
    3.6660    3.2787    0.0000 C   0  0
    8.1591   -0.6978    0.0000 C   0  0
    7.3351    1.7662    0.0000 C   0  0
   -6.3876    3.5028    0.0000 C   0  0
   -0.6756   -5.0219    0.0000 C   0  0
    7.0358    0.2964    0.0000 C   0  0
    3.8914   -2.1249    0.0000 C   0  0
   -2.0982   -5.4976    0.0000 C   0  0
   -4.5701    1.8943    0.0000 C   0  0  1  0  0  0
   -6.9859    2.1273    0.0000 C   0  0  1  0  0  0
    4.4900    0.8148    0.0000 C   0  0
    1.3455   -1.6065    0.0000 C   0  0
    4.7893    2.2846    0.0000 C   0  0
    1.9442    1.3332    0.0000 C   0  0
    1.0462   -3.0763    0.0000 C   0  0
    2.2435    2.8030    0.0000 C   0  0
   -0.6017    1.8516    0.0000 C   0  0
    5.6132   -0.1794    0.0000 C   0  0
    0.2223   -0.6124    0.0000 Cl  0  0
    9.2823   -1.6919    0.0000 N   0  0
   -3.2215   -4.5035    0.0000 N   0  0
    6.2119    2.7603    0.0000 N   0  0
    5.3139   -1.6492    0.0000 N   0  0
    0.5216    0.8575    0.0000 N   0  0
   -4.8945    3.3588    0.0000 N   0  0
   -8.2913    2.8662    0.0000 O   0  0
   -0.3024    3.3214    0.0000 O   0  0
    1.1202    3.7971    0.0000 O   0  0
   -0.3763   -3.5520    0.0000 O   0  0
   -2.8482    3.8398    0.0000 H   0  0
   -2.3235   -0.0940    0.0000 H   0  0
   -3.9483    0.5292    0.0000 H   0  0
   -7.8572    0.9063    0.0000 H   0  0
  1  3  1  0
  2 39  1  0
  3 42  1  0
  4  5  2  0
  4  6  1  0
  5 11  1  0
  6 22  2  0
  7  9  2  0
  7 21  1  0
  8 44  1  0
  8 10  2  0
  8 23  1  0
  9 29  1  0
 10 45  1  0
 10 31  1  0
 11 35  2  0
 12 21  2  0
 12 26  1  0
 13 23  1  0
 13 24  1  0
 14 25  2  0
 14 28  1  0
 15 27  2  0
 15 30  1  0
 16 20  1  0
 16 34  3  0
 17 20  2  0
 17 36  1  0
 18 24  1  0
 18 39  1  0
 19 22  1  0
 19 43  1  0
 20 32  1  0
 21 37  1  0
 22 35  1  0
 23 46  1  6
 23 39  1  0
 24 47  1  1
 24 40  1  0
 25 27  1  0
 25 32  1  0
 26 29  2  0
 26 33  1  0
 27 36  1  0
 28 30  2  0
 28 38  1  0
 29 43  1  0
 30 42  1  0
 31 38  2  0
 31 41  1  0
 32 37  2  3
M  END
"""

)CTAB";
    std::unique_ptr<RWMol> m(MolBlockToMol(molblock, false, false));
    REQUIRE(m);
    m->updatePropertyCache();
    MolOps::fastFindRings(*m);
    MolOps::setBondStereoFromDirections(*m);
    MolOps::removeHs(*m, false, false, false);
    std::unique_ptr<RWMol> res((RWMol *)nrml.normalize(*m));
    REQUIRE(res);
    MolOps::sanitizeMol(*res);
    MolOps::assignStereochemistry(*res);
    CHECK(MolToSmiles(*res) ==
          "CCOc1cc2[nH]cc(C#N)c(=Nc3ccc(OCc4ccccn4)c(Cl)c3)c2cc1NC(=O)/C=C/"
          "[C@H]1C[C@H](O)CN1C");
  }
}
TEST_CASE("problems with uncharging HS- from mol file", "[normalizer]") {
  SECTION("example1") {
    std::string mb = R"CTAB(
  SciTegic12231509382D

  1  0  0  0  0  0            999 V2000
   13.0092   -4.9004    0.0000 S   0  5
M  CHG  1   1  -1
M  END)CTAB";
    std::unique_ptr<ROMol> m(MolBlockToMol(mb));
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    CHECK(MolToSmiles(*outm) == "S");
  }
}

TEST_CASE("explicit Hs and Ns when neutralizing", "[normalizer]") {
  SECTION("example1") {
    std::string molblock = R"CTAB(
  Mrv1810 10301909502D          

  2  1  0  0  0  0            999 V2000
   -3.0000    0.6316    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1750    0.6316    0.0000 N   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  1   2  -1
M  END
)CTAB";
    std::unique_ptr<RWMol> m(MolBlockToMol(molblock, false, false));
    REQUIRE(m);
    m->updatePropertyCache();
    MolStandardize::Uncharger uc;
    std::unique_ptr<ROMol> res((ROMol *)uc.uncharge(*m));
    REQUIRE(res);
    CHECK(res->getAtomWithIdx(1)->getFormalCharge() == 0);
    CHECK(res->getAtomWithIdx(1)->getTotalNumHs() == 2);
    auto mb = MolToMolBlock(*res);
    // should be no valence markers in the output mol block:
    CHECK(mb.find("0.0000 N   0  0  0  0  0  0") != std::string::npos);
  }
}

TEST_CASE("fragment remover not considering bond counts", "[fragments][bug]") {
  std::string salts = R"DATA(Benethamine	C(Cc1ccccc1)NCc2ccccc2
Chloride	Cl
)DATA";
  std::istringstream iss(salts);
  bool leave_last = false;
  MolStandardize::FragmentRemover rmv(iss, leave_last);

  SECTION("example that should not be removed") {
    std::string molblock = R"CTAB(
  SciTegic11261411092D

 17 18  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 Cl  0  0
    2.2393    0.5156    0.0000 N   0  0
    3.6682    0.5156    0.0000 C   0  0
    2.9538    0.1031    0.0000 C   0  0
    3.6682    1.3406    0.0000 C   0  0
    2.9538   -0.7219    0.0000 C   0  0
    4.3827    0.1031    0.0000 C   0  0
    2.9538    1.7531    0.0000 C   0  0
    4.3827    1.7531    0.0000 C   0  0
    2.2393    1.3406    0.0000 C   0  0
    3.6682   -1.1344    0.0000 C   0  0
    2.2393   -1.1344    0.0000 C   0  0
    5.0972    0.5156    0.0000 C   0  0
    5.0972    1.3406    0.0000 C   0  0
    3.6682   -1.9594    0.0000 C   0  0
    2.2393   -1.9594    0.0000 C   0  0
    2.9538   -2.3719    0.0000 C   0  0
  2  4  1  0
  2 10  1  0
  3  4  1  0
  3  5  1  0
  3  7  2  0
  4  6  1  0
  5  8  1  0
  5  9  2  0
  6 11  2  0
  6 12  1  0
  7 13  1  0
  8 10  1  0
  9 14  1  0
 11 15  1  0
 12 16  2  0
 13 14  2  0
 15 17  2  0
 16 17  1  0
M  END)CTAB";
    std::unique_ptr<RWMol> m(MolBlockToMol(molblock));
    REQUIRE(m);
    m->updatePropertyCache();

    std::unique_ptr<ROMol> sm(rmv.remove(*m));
    REQUIRE(sm);
    CHECK(sm->getNumAtoms() == 16);
  }

  SECTION("example that should be removed") {
    std::string molblock = R"CTAB(
  Mrv1810 11071914502D          

 17 17  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    2.2393    0.5156    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.6682    0.5156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9538    0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6682    1.3406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9538   -0.7219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3827    0.1031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9538    1.7531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.3827    1.7531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2393    1.3406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6682   -1.1344    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2393   -1.1344    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0972    0.5156    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0972    1.3406    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6682   -1.9594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2393   -1.9594    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9538   -2.3719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  4  1  0  0  0  0
  2 10  1  0  0  0  0
  3  5  1  0  0  0  0
  3  7  2  0  0  0  0
  4  6  1  0  0  0  0
  5  8  1  0  0  0  0
  5  9  2  0  0  0  0
  6 11  2  0  0  0  0
  6 12  1  0  0  0  0
  7 13  1  0  0  0  0
  8 10  1  0  0  0  0
  9 14  1  0  0  0  0
 11 15  1  0  0  0  0
 12 16  2  0  0  0  0
 13 14  2  0  0  0  0
 15 17  2  0  0  0  0
 16 17  1  0  0  0  0
M  END
)CTAB";
    std::unique_ptr<RWMol> m(MolBlockToMol(molblock));
    REQUIRE(m);
    m->updatePropertyCache();

    std::unique_ptr<ROMol> sm(rmv.remove(*m));
    REQUIRE(sm);
    CHECK(sm->getNumAtoms() == 0);
  }
}

TEST_CASE("github #2792: carbon in the uncharger", "[uncharger][bug]") {
  SECTION("carbocation 1") {
    auto m = "C[CH2+]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(1)->getTotalNumHs() == 3);
  }
  SECTION("boron cation") {
    auto m = "C[BH+]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(1)->getTotalNumHs() == 2);
  }
  SECTION("carbanion 1") {
    auto m = "C[CH2-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(1)->getTotalNumHs() == 3);
  }
  SECTION("carbocation 2") {
    auto m = "CN1C=CN[CH+]1"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(5)->getTotalNumHs() == 2);
  }
  SECTION("carbocation 2 without sanitization") {
    SmilesParserParams params;
    params.sanitize = false;
    std::unique_ptr<ROMol> m(SmilesToMol("CN1C=CN[CH+]1", params));
    REQUIRE(m);
    m->updatePropertyCache();
    MolStandardize::Uncharger uncharger;
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 0);
    CHECK(outm->getAtomWithIdx(5)->getTotalNumHs() == 2);
  }
}

TEST_CASE("github #2965: molecules properties not retained after cleanup",
          "[cleanup][bug]") {
  SECTION("example 1") {
    MolStandardize::CleanupParameters params;
    std::unique_ptr<RWMol> m(SmilesToMol("Cl.c1cnc(OCCCC2CCNCC2)cn1"));
    REQUIRE(m);
    m->setProp("testing_prop", "1234");
    std::unique_ptr<RWMol> res(MolStandardize::cleanup(*m, params));
    REQUIRE(res);
    auto x = res->getDict();
    CHECK(x.getVal<std::string>("testing_prop") == "1234");
  }
}

TEST_CASE(
    "github #2970: chargeParent() segmentation fault when standardization is "
    "skipped i.e. skip_standardize is set to true") {
  auto m = "COC=1C=CC(NC=2N=CN=C3NC=NC23)=CC1"_smiles;
  REQUIRE(m);
  MolStandardize::CleanupParameters params;
  std::unique_ptr<RWMol> res(MolStandardize::cleanup(*m, params));

  std::unique_ptr<ROMol> outm(MolStandardize::chargeParent(*res, params, true));

  REQUIRE(outm);
  CHECK(MolToSmiles(*outm) == "COc1ccc(Nc2ncnc3[nH]cnc23)cc1");
}
