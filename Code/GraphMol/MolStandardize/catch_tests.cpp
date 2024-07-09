//
//  Copyright (C) 2019-2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/MolStandardize/Charge.h>
#include <GraphMol/MolStandardize/Tautomer.h>

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

TEST_CASE("uncharger 'force' option") {
  SECTION("force=false (default)") {
    MolStandardize::Uncharger uncharger;
    auto m1 = "C[N+](C)(C)CC([O-])C[O-]"_smiles;
    REQUIRE(m1);
    std::unique_ptr<ROMol> outm1(uncharger.uncharge(*m1));
    REQUIRE(outm1);
    CHECK(MolToSmiles(*outm1) == "C[N+](C)(C)CC([O-])CO");
    auto m2 = "C[B-](C)(C)CC([NH3+])C[NH3+]"_smiles;
    REQUIRE(m2);
    std::unique_ptr<ROMol> outm2(uncharger.uncharge(*m2));
    REQUIRE(outm2);
    CHECK(MolToSmiles(*outm2) == "C[B-](C)(C)CC(N)C[NH3+]");
  }
  SECTION("force=true") {
    MolStandardize::Uncharger uncharger(false, true);
    auto m1 = "C[N+](C)(C)CC([O-])C[O-]"_smiles;
    REQUIRE(m1);
    std::unique_ptr<ROMol> outm1(uncharger.uncharge(*m1));
    REQUIRE(outm1);
    CHECK(MolToSmiles(*outm1) == "C[N+](C)(C)CC(O)CO");
    auto m2 = "C[B-](C)(C)CC([NH3+])C[NH3+]"_smiles;
    REQUIRE(m2);
    std::unique_ptr<ROMol> outm2(uncharger.uncharge(*m2));
    REQUIRE(outm2);
    CHECK(MolToSmiles(*outm2) == "C[B-](C)(C)CC(N)CN");
  }
  SECTION("force=true doesn't alter nitro groups") {
    auto m = "CCC[N+](=O)[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(false, true);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CCC[N+](=O)[O-]");
  }
  SECTION("force=true doesn't alter n-oxides") {
    auto m = "[O-][n+]1ccccc1"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(false, true);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "[O-][n+]1ccccc1");
  }
  SECTION("tetramethylammonium acetate (force=false)") {
    auto m = "C[N+](C)(C)C.CC(=O)[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(true, false);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CC(=O)[O-].C[N+](C)(C)C");
  }
  SECTION("tetramethylammonium acetate (force=true)") {
    auto m = "C[N+](C)(C)C.CC(=O)[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(true, true);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "CC(=O)O.C[N+](C)(C)C");
  }
  SECTION("tetramethylammonium nitrate (force=false)") {
    auto m = "C[N+](C)(C)C.O=[N+]([O-])[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(true, false);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)C.O=[N+]([O-])[O-]");
  }
  SECTION("tetramethylammonium nitrate (force=true)") {
    auto m = "C[N+](C)(C)C.O=[N+]([O-])[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(true, true);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)C.O=[N+]([O-])O");
  }
  SECTION("bookkeeping (force=false)") {
    auto m = "O=[N+]([O-])[O-].O=[N+]([O-])[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(true, false);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "O=[N+]([O-])O.O=[N+]([O-])O");
  }
  SECTION("bookkeeping (force=true)") {
    auto m = "O=[N+]([O-])[O-].O=[N+]([O-])[O-]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uncharger(true, true);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "O=[N+]([O-])O.O=[N+]([O-])O");
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
    v2::SmilesParse::SmilesParserParams ps;
    ps.sanitize = false;
    auto m = v2::SmilesParse::MolFromSmiles("F[B-](F)(F)F.[NH4+2]CCC",
                                            ps);  // totally bogus structure
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == -1);
    CHECK(outm->getAtomWithIdx(5)->getFormalCharge() == 1);
    CHECK(MolToSmiles(*outm) == "CCC[NH3+].F[B-](F)(F)F");
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
  SECTION("zwitterion also including an N-oxide") {
    auto m = "C[N+](C)(C)C(C(=O)[O-])c1cc[n+]([O-])cc1"_smiles;
    REQUIRE(m);
    bool canonicalOrdering = true;
    MolStandardize::Uncharger uncharger(canonicalOrdering);
    std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
    REQUIRE(outm);
    CHECK(MolToSmiles(*outm) == "C[N+](C)(C)C(C(=O)[O-])c1cc[n+]([O-])cc1");
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

TEST_CASE("update parameters from JSON") {
  std::string rdbase = std::getenv("RDBASE");

  // a few tests to make sure the basics work
  MolStandardize::CleanupParameters params;
  CHECK(params.maxRestarts == 200);
  CHECK(params.tautomerReassignStereo == true);

  MolStandardize::updateCleanupParamsFromJSON(params,
                                              R"JSON({"maxRestarts":12,
  "tautomerReassignStereo":false,
  "fragmentFile":"foo.txt"})JSON");
  CHECK(params.maxRestarts == 12);
  CHECK(params.tautomerReassignStereo == false);
  CHECK(params.fragmentFile == "foo.txt");
}

TEST_CASE("provide normalizer parameters as data") {
  std::vector<std::pair<std::string, std::string>> tfs{
      {"Bad amide tautomer1",
       "[C:1]([OH1;D1:2])=;!@[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]"},
      {"Bad amide tautomer2",
       "[C:1]([OH1;D1:2])=;!@[NH0:3]>>[C:1](=[OH0:2])-[NH1:3]"}};
  SECTION("example1") {
    MolStandardize::Normalizer nrml(tfs, 20);
    auto m = "Cl.Cl.OC(=N)NCCCCCCCCCCCCNC(O)=N"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(nrml.normalize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "Cl.Cl.NC(=O)NCCCCCCCCCCCCNC(N)=O");
  }
  SECTION("example2") {
    MolStandardize::Normalizer nrml(tfs, 20);
    auto m = "OC(=N)NCCCCCCCCCCCCNC(O)=N"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(nrml.normalize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "NC(=O)NCCCCCCCCCCCCNC(N)=O");
  }
}

TEST_CASE("provide normalizer parameters as JSON") {
  SECTION("example1") {
    std::string json = R"JSON({"normalizationData":[
      {"name":"silly 1","smarts":"[Cl:1]>>[F:1]"},
      {"name":"silly 2","smarts":"[Br:1]>>[F:1]"}
    ]})JSON";
    MolStandardize::CleanupParameters params;
    MolStandardize::updateCleanupParamsFromJSON(params, json);
    CHECK(params.normalizationData.size() == 2);

    MolStandardize::Normalizer nrml(params.normalizationData, 20);
    auto m = "ClCCCBr"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(nrml.normalize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "FCCCF");
  }
}

TEST_CASE("provide charge parameters as data") {
  std::vector<std::tuple<std::string, std::string, std::string>> params{
      {"-CO2H", "C(=O)[OH]", "C(=O)[O-]"}, {"phenol", "c[OH]", "c[O-]"}};
  SECTION("example1") {
    MolStandardize::Reionizer reion(params);
    auto m = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(reion.reionize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "O=C([O-])c1cccc(O)c1");
  }
  SECTION("example2") {
    MolStandardize::Reionizer reion(params);
    auto m = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(reion.reionize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "O=S([O-])c1ccc(S(=O)(=O)O)cc1");
  }
}

TEST_CASE("provide charge parameters as JSON") {
  SECTION("example1") {
    std::string json = R"JSON({"acidbaseData":[
      {"name":"-CO2H","acid":"C(=O)[OH]","base":"C(=O)[O-]"},
      {"name":"phenol","acid":"c[OH]","base":"c[O-]"}
    ]})JSON";
    MolStandardize::CleanupParameters params;
    MolStandardize::updateCleanupParamsFromJSON(params, json);
    CHECK(params.acidbaseData.size() == 2);

    MolStandardize::Reionizer reion(params.acidbaseData);
    auto m = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(reion.reionize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "O=C([O-])c1cccc(O)c1");
    m = "C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O"_smiles;
    REQUIRE(m);
    res.reset(reion.reionize(*m));
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "O=S([O-])c1ccc(S(=O)(=O)O)cc1");
  }
}

TEST_CASE("provide tautomer parameters as JSON") {
  SECTION("example1") {
    std::string json = R"JSON({"tautomerTransformData":[
      {"name":"1,3 (thio)keto/enol f","smarts":"[CX4!H0]-[C]=[O,S,Se,Te;X1]","bonds":"","charges":""},
      {"name":"1,3 (thio)keto/enol r","smarts":"[O,S,Se,Te;X2!H0]-[C]=[C]"}
    ]})JSON";
    MolStandardize::CleanupParameters params;
    MolStandardize::updateCleanupParamsFromJSON(params, json);
    CHECK(params.tautomerTransformData.size() == 2);
    MolStandardize::TautomerEnumerator te(params);
    auto m = "CCC=O"_smiles;
    REQUIRE(m);
    auto tauts = te.enumerate(*m);
    CHECK(tauts.size() == 2);
    CHECK(MolToSmiles(*tauts[0]) == "CC=CO");
    CHECK(MolToSmiles(*tauts[1]) == "CCC=O");
  }
  SECTION("example 2") {
    std::string json = R"JSON({"tautomerTransformData":[
        {"name":"isocyanide f", "smarts":"[C-0!H0]#[N+0]", "bonds":"#", "charges":"-+"},
        {"name":"isocyanide r", "smarts":"[N+!H0]#[C-]", "bonds":"#", "charges":"-+"}
    ]})JSON";
    MolStandardize::CleanupParameters params;
    MolStandardize::updateCleanupParamsFromJSON(params, json);
    CHECK(params.tautomerTransformData.size() == 2);
    MolStandardize::TautomerEnumerator te(params);
    auto m = "C#N"_smiles;
    REQUIRE(m);
    auto tauts = te.enumerate(*m);
    CHECK(tauts.size() == 2);
    CHECK(MolToSmiles(*tauts[0]) == "C#N");
    CHECK(MolToSmiles(*tauts[1]) == "[C-]#[NH+]");
  }
  SECTION("example3") {
    std::string json = R"JSON({"tautomerTransformData":[
      {"name":"1,3 (thio)keto/enol f","smarts":"[CX4!H0]-[C]=[O,S,Se,Te;X1]","bonds":"","charges":""},
      {"name":"1,3 (thio)keto/enol r","smarts":"[O,S,Se,Te;X2!H0]-[C]=[C]"}
    ]})JSON";
    MolStandardize::CleanupParameters params;
    MolStandardize::updateCleanupParamsFromJSON(params, json);
    CHECK(params.tautomerTransformData.size() == 2);
    auto m = "CCC=O"_smiles;
    REQUIRE(m);
    std::unique_ptr<RWMol> nm{MolStandardize::canonicalTautomer(m.get())};
    CHECK(MolToSmiles(*nm) == "CCC=O");
  }
}

TEST_CASE("provide fragment parameters as JSON") {
  SECTION("example1") {
    std::string json = R"JSON({"fragmentData":[
        {"name":"hydrogen", "smarts":"[H]"}, 
        {"name":"fluorine", "smarts":"[F]"}, 
        {"name":"chlorine", "smarts":"[Cl]"}
    ]})JSON";
    MolStandardize::CleanupParameters params;
    MolStandardize::updateCleanupParamsFromJSON(params, json);
    CHECK(params.fragmentData.size() == 3);
    std::unique_ptr<MolStandardize::FragmentRemover> fm{
        MolStandardize::fragmentRemoverFromParams(params, true)};
    auto m = "[F-].[Cl-].[Br-].CC"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> nm{fm->remove(*m)};
    CHECK(MolToSmiles(*nm) == "CC.[Br-]");
  }
}

TEST_CASE("tautomer parent") {
  SECTION("example1") {
    auto m = "[O-]c1ccc(C(=O)O)cc1CC=CO"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> nm{MolStandardize::tautomerParent(*m)};
    CHECK(MolToSmiles(*nm) == "O=CCCc1cc(C(=O)[O-])ccc1O");
    MolStandardize::tautomerParentInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=CCCc1cc(C(=O)[O-])ccc1O");
  }
}

TEST_CASE("stereo parent") {
  SECTION("example1") {
    auto m = "C[C@](F)(Cl)C/C=C/[C@H](F)Cl"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> nm{MolStandardize::stereoParent(*m)};
    CHECK(MolToSmiles(*nm) == "CC(F)(Cl)CC=CC(F)Cl");
  }
}

TEST_CASE("isotope parent") {
  SECTION("example1") {
    auto m = "[12CH3][13CH3]"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> nm{MolStandardize::isotopeParent(*m)};
    CHECK(MolToSmiles(*nm) == "CC");
  }
  SECTION("attached D") {
    // this behavior - leaving H atoms with no isotope info - is intentional
    // It may be that we're working with molecules which include Hs and we don't
    // want to just automatically remove them.
    auto m = "O[2H]"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> nm{MolStandardize::isotopeParent(*m)};
    CHECK(MolToSmiles(*nm) == "[H]O");
  }
}

TEST_CASE("super parent") {
  SECTION("example1") {
    auto m = "[O-]c1c([12C@H](F)Cl)c(O[2H])c(C(=O)O)cc1CC=CO.[Na+]"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> nm{MolStandardize::superParent(*m)};
    CHECK(MolToSmiles(*nm) == "O=CCCc1cc(C(=O)O)c(O)c(C(F)Cl)c1O");
    MolStandardize::superParentInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=CCCc1cc(C(=O)O)c(O)c(C(F)Cl)c1O");
  }
}

TEST_CASE(
    "Github #4260: Exception thrown by reionizer when dealing with Mg+2") {
  SECTION("reported") {
    auto m = "[Mg].OC(=O)c1ccccc1C"_smiles;
    REQUIRE(m);
    std::unique_ptr<RWMol> m2(MolStandardize::reionize(m.get()));
    REQUIRE(m2);
    CHECK(m2->getAtomWithIdx(0)->getFormalCharge() == 2);
    CHECK(m2->getAtomWithIdx(1)->getFormalCharge() == -1);
  }
}

TEST_CASE("Github #5008: bad tautomers for phosphorous compounds") {
  SECTION("as reported") {
    auto m = "NP(=O)(O)N(CCCl)CCCl"_smiles;
    REQUIRE(m);
    MolStandardize::TautomerEnumerator tenum;
    auto tauts = tenum.enumerate(*m);
    CHECK(tauts.size() == 1);
  }
  SECTION("P which should tautomerize") {
    auto m = "CP(O)C"_smiles;
    REQUIRE(m);
    MolStandardize::TautomerEnumerator tenum;
    auto tauts = tenum.enumerate(*m);
    CHECK(tauts.size() == 2);
  }
  SECTION("Canonical version") {
    auto m = "CP(O)C"_smiles;
    REQUIRE(m);
    std::unique_ptr<RWMol> ct(MolStandardize::canonicalTautomer(m.get()));
    REQUIRE(ct);
    CHECK(MolToSmiles(*ct) == "C[PH](C)=O");
  }
}

TEST_CASE("Github #5169: Standardization via RDKit breaks molecules",
          "[uncharger]") {
  SECTION("basics") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::vector<std::string> smis = {"C[O+](C)C", "[H]/[O+]=C/Cl"};
    for (const auto &smi : smis) {
      std::unique_ptr<RWMol> m{SmilesToMol(smi, ps)};
      REQUIRE(m);
      m->updatePropertyCache(false);
      MolStandardize::Uncharger uncharger;
      std::unique_ptr<ROMol> outm(uncharger.uncharge(*m));
      REQUIRE(outm);
      INFO("failing for smiles " << smi);
      CHECK(outm->getAtomWithIdx(1)->getFormalCharge() == 1);
    }
  }
}

TEST_CASE("asymmetric imine tautomer generation", "[tautomers]") {
  SECTION("basics") {
    MolStandardize::TautomerEnumerator tenum;
    // clang-format off
    std::vector<std::pair<std::string, unsigned>> data = {
        {"C=C1NNC(=O)N1*", 2},
        {"CC1=NN=C(O)N1*", 2},
        {"C-C=NC", 1},
        {"C-C=N", 2},
        {"C-C=Nc1ccccc1", 2},
    };
    // clang-format on
    for (const auto &pr : data) {
      INFO(pr.first);
      std::unique_ptr<RWMol> m(SmilesToMol(pr.first));
      auto res = tenum.enumerate(*m);
      CHECK(res.size() == pr.second);
    }
  }
}

TEST_CASE("Github 5317: standardization failing with zwitterionic sulfone") {
  SECTION("basics") {
    auto m = "C[S+2]([O-])([O-])C([O-])C(=O)O"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uc;
    std::unique_ptr<ROMol> res{uc.uncharge(*m)};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "C[S+2]([O-])([O-])C(O)C(=O)O");
  }
  SECTION("don't overdo it") {
    auto m = "C[S+2]([O-])([O-])C([O-])C(=O)O.[Na+]"_smiles;
    REQUIRE(m);
    MolStandardize::Uncharger uc;
    std::unique_ptr<ROMol> res{uc.uncharge(*m)};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "C[S+2]([O-])([O-])C([O-])C(=O)O.[Na+]");
  }
}

TEST_CASE("Github 5318: standardizing unsanitized molecules should work") {
  SmilesParserParams ps;
  ps.sanitize = false;
  ps.removeHs = false;
  std::unique_ptr<RWMol> m{SmilesToMol("C[S+2]([O-])([O-])C([O-])C(=O)O", ps)};
  REQUIRE(m);
  std::unique_ptr<RWMol> m2{SmilesToMol("Cc1[nH]ncc1.[Cl]", ps)};
  REQUIRE(m2);
  SECTION("reionizer") {
    MolStandardize::Reionizer reion;
    std::unique_ptr<ROMol> res{reion.reionize(*m)};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "C[S+2]([O-])([O-])C(O)C(=O)[O-]");
  }
  SECTION("uncharger") {
    MolStandardize::Uncharger uc;
    std::unique_ptr<ROMol> res{uc.uncharge(*m)};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "C[S+2]([O-])([O-])C(O)C(=O)O");
  }
  SECTION("normalizer") {
    std::unique_ptr<ROMol> res{MolStandardize::normalize(m.get())};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "CS(=O)(=O)C([O-])C(=O)O");
  }
  SECTION("tautomer") {
    std::unique_ptr<ROMol> res{MolStandardize::canonicalTautomer(m2.get())};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "Cc1cc[nH]n1.Cl");
    RWMol cp(*m2);
    MolStandardize::canonicalTautomerInPlace(cp);
    CHECK(MolToSmiles(cp) == "Cc1cc[nH]n1.Cl");
  }
  SECTION("fragments") {
    std::unique_ptr<ROMol> res{MolStandardize::removeFragments(m2.get())};
    REQUIRE(res);
    CHECK(MolToSmiles(*res) == "Cc1ccn[nH]1");
  }
}

TEST_CASE("Github #5320: cleanup() and stereochemistry") {
  SECTION("basics") {
    auto m = "Cl[C@](O)([O-])C(=O)O"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    std::unique_ptr<RWMol> m2{MolStandardize::cleanup(m.get())};
    REQUIRE(m2);
    CHECK(m2->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(MolToSmiles(*m2) == "O=C([O-])C(O)(O)Cl");
  }
}

TEST_CASE("Github #5402: order dependence of tautomer transforms") {
  SECTION("as-reported") {
    MolStandardize::TautomerEnumerator te;

    auto m1 = "c1ccc([C@@H](CC2=NCCN2)c2ccccn2)cc1"_smiles;
    REQUIRE(m1);
    auto m2 = "C([C@H](C1=CC=CC=C1)C2=NC=CC=C2)C3=NCCN3"_smiles;
    REQUIRE(m2);
    std::cerr << " * - * - * - * m1" << std::endl;
    std::unique_ptr<ROMol> res1{te.canonicalize(*m1)};
    REQUIRE(res1);
    std::cerr << " * - * - * - * m2" << std::endl;
    std::unique_ptr<ROMol> res2{te.canonicalize(*m2)};
    REQUIRE(res2);
    CHECK(MolToSmiles(*res1) == MolToSmiles(*res2));
  }
  SECTION("zoom") {
    MolStandardize::CleanupParameters params;
    const std::vector<
        std::tuple<std::string, std::string, std::string, std::string>>
        tTransforms{
            std::make_tuple(std::string("special imine r1"),
                            std::string("[Cz0R0X4!H0]-[c]=[nz0]"),
                            std::string(""), std::string("")),
            std::make_tuple(std::string("special imine r2"),
                            std::string("[Cz0R0X4!H0]-[c](=c)-[nz0]"),
                            std::string("==-"), std::string("")),
        };

    params.tautomerTransformData = tTransforms;
    MolStandardize::TautomerEnumerator te(params);

    auto m1 = "c1ccc([C@@H](CC2=NCCN2)c2ccccn2)cc1"_smiles;
    REQUIRE(m1);
    auto m2 = "C([C@H](C1=CC=CC=C1)C2=NC=CC=C2)C3=NCCN3"_smiles;
    REQUIRE(m2);

    std::cerr << " * - * - * - * m1" << std::endl;
    std::unique_ptr<ROMol> res1{te.canonicalize(*m1)};
    REQUIRE(res1);
    std::cerr << " * - * - * - * m2" << std::endl;
    std::unique_ptr<ROMol> res2{te.canonicalize(*m2)};
    REQUIRE(res2);
    CHECK(MolToSmiles(*res1) == MolToSmiles(*res2));
  }
}

TEST_CASE("Github 5784: kekulization error when enumerating tautomers") {
  std::vector<std::string> smis{"NC1=NC=NC(C)=C1", "CC1N=CN(C)C(=O)C=1",
                                "CC1=CC=CC(=O)N1C"};
  for (const auto &smi : smis) {
    INFO(smi);
    std::unique_ptr<ROMol> m{SmilesToMol(smi)};
    REQUIRE(m);
    MolStandardize::TautomerEnumerator te;
    std::unique_ptr<ROMol> res(te.canonicalize(*m));
    REQUIRE(res);
  }
}

TEST_CASE("in place operations") {
  SECTION("reionizer") {
    MolStandardize::Reionizer reion;
    auto m = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    REQUIRE(m);
    reion.reionizeInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=C([O-])c1cccc(O)c1");
  }
  SECTION("reionize") {
    auto m = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    REQUIRE(m);
    MolStandardize::reionizeInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=C([O-])c1cccc(O)c1");
  }
  SECTION("uncharge") {
    MolStandardize::Uncharger unchg;
    auto m = "c1cc([O-])cc(C(=O)O)c1"_smiles;
    REQUIRE(m);
    unchg.unchargeInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=C(O)c1cccc(O)c1");
  }
  SECTION("normalizer") {
    MolStandardize::Normalizer nrml;
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m{SmilesToMol("O=N(=O)-CC-N(=O)=O", ps)};
    REQUIRE(m);
    nrml.normalizeInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=[N+]([O-])CC[N+](=O)[O-]");
    m.reset(SmilesToMol("OCCN", ps));
    REQUIRE(m);
    nrml.normalizeInPlace(*m);
    CHECK(MolToSmiles(*m) == "NCCO");
  }
  SECTION("normalize") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m{SmilesToMol("O=N(=O)-CC-N(=O)=O", ps)};
    REQUIRE(m);
    MolStandardize::normalizeInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=[N+]([O-])CC[N+](=O)[O-]");
    m.reset(SmilesToMol("OCCN", ps));
    REQUIRE(m);
    MolStandardize::normalizeInPlace(*m);
    CHECK(MolToSmiles(*m) == "NCCO");
  }
  SECTION("FragmentRemover") {
    auto m = "CCCC.Cl.[Na]"_smiles;
    REQUIRE(m);
    MolStandardize::FragmentRemover fragremover;
    RWMol cp1(*m);
    fragremover.removeInPlace(cp1);
    CHECK(MolToSmiles(cp1) == "CCCC");
    RWMol cp2(*m);
    MolStandardize::removeFragmentsInPlace(cp2);
    CHECK(MolToSmiles(cp2) == "CCCC");
  }
  SECTION("FragmentParent") {
    auto m = "CCCC.Cl.[Na]"_smiles;
    REQUIRE(m);
    RWMol cp1(*m);
    MolStandardize::fragmentParentInPlace(cp1);
    // note: this isn't a nice answer, and it should be
    // fixed, but it is what the code currently generates
    CHECK(MolToSmiles(cp1) == "[CH2-]CCC");
  }
  SECTION("ChargeParent") {
    auto m = "[O-]C(=O)CCC.[Na+]"_smiles;
    REQUIRE(m);
    RWMol cp1(*m);
    MolStandardize::chargeParentInPlace(cp1);
    CHECK(MolToSmiles(cp1) == "CCCC(=O)O");
  }
  SECTION("IsotopeParent") {
    auto m = "[13CH3]C"_smiles;
    REQUIRE(m);
    RWMol cp1(*m);
    MolStandardize::isotopeParentInPlace(cp1);
    CHECK(MolToSmiles(cp1) == "CC");
  }
  SECTION("StereoParent") {
    auto m = "F[C@H](O)Cl"_smiles;
    REQUIRE(m);
    RWMol cp1(*m);
    MolStandardize::stereoParentInPlace(cp1);
    CHECK(MolToSmiles(cp1) == "OC(F)Cl");
  }
  SECTION("cleanup") {
    SmilesParserParams ps;
    ps.sanitize = false;
    // silly ugly example which ensures disconnection, normalization, and
    // reionization
    std::unique_ptr<RWMol> m{
        SmilesToMol("O=N(=O)-C(O[Fe])C(C(=O)O)C-N(=O)=O", ps)};
    REQUIRE(m);
    MolStandardize::cleanupInPlace(*m);
    CHECK(MolToSmiles(*m) == "O=C([O-])C(C[N+](=O)[O-])C(O)[N+](=O)[O-].[Fe+]");
  }
  SECTION("disconnect organometallics") {
    auto m("[CH2-](->[K+])c1ccccc1"_smiles);
    TEST_ASSERT(m);
    MolStandardize::disconnectOrganometallicsInPlace(*m);
    TEST_ASSERT(MolToSmiles(*m) == "[CH2-]c1ccccc1.[K+]");
  }
}

TEST_CASE("cleanup with multiple mols") {
  SmilesParserParams ps;
  ps.sanitize = false;
  // silly ugly examples which ensures disconnection, normalization, and
  // reionization
  std::vector<std::pair<std::string, std::string>> data = {
      {"O=N(=O)-C(O[Fe])C(C(=O)O)C-N(=O)=O",
       "O=C([O-])C(C[N+](=O)[O-])C(O)[N+](=O)[O-].[Fe+]"},
      {"O=N(=O)-CC(O[Fe])C(C(=O)O)C-N(=O)=O",
       "O=C([O-])C(C[N+](=O)[O-])C(O)C[N+](=O)[O-].[Fe+]"},
      {"O=N(=O)-CCC(O[Fe])C(C(=O)O)C-N(=O)=O",
       "O=C([O-])C(C[N+](=O)[O-])C(O)CC[N+](=O)[O-].[Fe+]"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &pr : data) {
    mols.emplace_back(SmilesToMol(pr.first, ps));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    MolStandardize::cleanupInPlace(molPtrs);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::cleanupInPlace(molPtrs, numThreads);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("normalize with multiple mols") {
  SmilesParserParams ps;
  ps.sanitize = false;
  std::vector<std::pair<std::string, std::string>> data = {
      {"O=N(=O)-CC-N(=O)=O", "O=[N+]([O-])CC[N+](=O)[O-]"},
      {"O=N(=O)-CCC-N(=O)=O", "O=[N+]([O-])CCC[N+](=O)[O-]"},
      {"O=N(=O)-CCCC-N(=O)=O", "O=[N+]([O-])CCCC[N+](=O)[O-]"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &pr : data) {
    mols.emplace_back(SmilesToMol(pr.first, ps));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    MolStandardize::normalizeInPlace(molPtrs);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::normalizeInPlace(molPtrs, numThreads);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("Reionize with multiple mols") {
  SmilesParserParams ps;
  ps.sanitize = false;
  std::vector<std::pair<std::string, std::string>> data = {
      {"c1cc([O-])cc(C(=O)O)c1", "O=C([O-])c1cccc(O)c1"},
      {"c1cc(C[O-])cc(C(=O)O)c1", "O=C([O-])c1cccc(CO)c1"},
      {"c1cc(CC[O-])cc(C(=O)O)c1", "O=C([O-])c1cccc(CCO)c1"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &pr : data) {
    mols.emplace_back(SmilesToMol(pr.first, ps));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    MolStandardize::reionizeInPlace(molPtrs);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::reionizeInPlace(molPtrs, numThreads);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("RemoveFragments with multiple mols") {
  SmilesParserParams ps;
  ps.sanitize = false;
  std::vector<std::pair<std::string, std::string>> data = {
      {"CCCC.Cl.[Na]", "CCCC"},
      {"CCCCO.Cl.[Na]", "CCCCO"},
      {"CCOC.Cl.[Na]", "CCOC"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &pr : data) {
    mols.emplace_back(SmilesToMol(pr.first, ps));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    MolStandardize::removeFragmentsInPlace(molPtrs);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::removeFragmentsInPlace(molPtrs, numThreads);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("charge with multiple mols") {
  auto params = MolStandardize::defaultCleanupParameters;
  params.preferOrganic = true;

  std::vector<std::pair<std::string, std::string>> data = {
      {"O=C([O-])c1ccccc1", "O=C(O)c1ccccc1"},
      {"C[NH+](C)(C).[Cl-]", "CN(C)C"},
      {"[N+](=O)([O-])[O-].[CH2]", "[CH2]"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &[insmi, outsmi] : data) {
    mols.emplace_back(SmilesToMol(insmi));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    int numThreads = 1;
    MolStandardize::chargeParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::chargeParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("isotope with multiple mols") {
  auto params = MolStandardize::defaultCleanupParameters;

  std::vector<std::pair<std::string, std::string>> data = {
      {"[13CH3]C", "CC"},
      {"[13CH3]C.C", "C.CC"},
      {"[13CH3][12CH3]", "CC"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &[insmi, outsmi] : data) {
    mols.emplace_back(SmilesToMol(insmi));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    int numThreads = 1;
    MolStandardize::isotopeParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::isotopeParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("fragments with multiple mols") {
  auto params = MolStandardize::defaultCleanupParameters;
  params.preferOrganic = true;

  std::vector<std::pair<std::string, std::string>> data = {
      {"O=C([O-])c1ccccc1", "O=C([O-])c1ccccc1"},
      {"C[NH+](C)(C).[Cl-]", "C[NH+](C)C"},
      {"[N+](=O)([O-])[O-].CC", "CC"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &[insmi, outsmi] : data) {
    mols.emplace_back(SmilesToMol(insmi));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    int numThreads = 1;
    MolStandardize::fragmentParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::fragmentParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("stereo with multiple mols") {
  auto params = MolStandardize::defaultCleanupParameters;

  std::vector<std::pair<std::string, std::string>> data = {
      {"F[C@H](O)Cl", "OC(F)Cl"},
      {"F[C@H](CCO)Cl", "OCCC(F)Cl"},
      {"F[C@H](CCO)Cl.F[C@H](O)Cl", "OC(F)Cl.OCCC(F)Cl"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 8; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &[insmi, outsmi] : data) {
    mols.emplace_back(SmilesToMol(insmi));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    int numThreads = 1;
    MolStandardize::stereoParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::stereoParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("tautomerParent with multiple mols") {
  auto params = MolStandardize::defaultCleanupParameters;

  std::vector<std::pair<std::string, std::string>> data = {
      {"[O-]c1ccc(C(=O)O)cc1CC=CO", "O=CCCc1cc(C(=O)[O-])ccc1O"},
      {"[O-]c1ccc(C(=O)O)cc1CC=CO.[Na+]", "O=CCCc1cc(C(=O)[O-])ccc1O.[Na+]"},
      {"[O-]c1ccc(C(=O)O)cc1C[13CH]=CO", "O=C[13CH2]Cc1cc(C(=O)[O-])ccc1O"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 5; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &[insmi, outsmi] : data) {
    mols.emplace_back(SmilesToMol(insmi));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    int numThreads = 1;
    MolStandardize::tautomerParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::tautomerParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}

TEST_CASE("superParent with multiple mols") {
  auto params = MolStandardize::defaultCleanupParameters;

  std::vector<std::pair<std::string, std::string>> data = {
      {"[O-]c1ccc(C(=O)O)cc1CC=CO", "O=CCCc1cc(C(=O)O)ccc1O"},
      {"[O-]c1ccc(C(=O)O)cc1CC=CO.[Na+]", "O=CCCc1cc(C(=O)O)ccc1O"},
      {"[O-]c1ccc(C(=O)O)cc1C[13CH]=CO", "O=CCCc1cc(C(=O)O)ccc1O"},
  };
  // bulk that up a bit
  for (auto iter = 0u; iter < 5; ++iter) {
    auto sz = data.size();
    for (auto i = 0u; i < sz; ++i) {
      data.push_back(data[i]);
    }
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  std::vector<RWMol *> molPtrs;
  for (const auto &[insmi, outsmi] : data) {
    mols.emplace_back(SmilesToMol(insmi));
    REQUIRE(mols.back());
    molPtrs.push_back(mols.back().get());
  }
  SECTION("basics") {
    int numThreads = 1;
    MolStandardize::superParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  SECTION("multithreaded") {
    int numThreads = 4;
    MolStandardize::superParentInPlace(molPtrs, numThreads, params);
    for (auto i = 0u; i < mols.size(); ++i) {
      REQUIRE(mols[i]);
      CHECK(MolToSmiles(*mols[i]) == data[i].second);
    }
  }
#endif
}