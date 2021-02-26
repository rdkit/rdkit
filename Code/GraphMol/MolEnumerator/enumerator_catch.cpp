//
//  Copyright (C) 2020 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include "MolEnumerator.h"

using namespace RDKit;

TEST_CASE("PositionVariation", "[MolEnumerator]") {
  auto mol1 = R"CTAB(
  Mrv2007 06232015292D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7083 2.415 0 0
M  V30 2 C -3.042 1.645 0 0
M  V30 3 C -3.042 0.105 0 0
M  V30 4 N -1.7083 -0.665 0 0
M  V30 5 C -0.3747 0.105 0 0
M  V30 6 C -0.3747 1.645 0 0
M  V30 7 * -0.8192 1.3883 0 0
M  V30 8 O -0.8192 3.6983 0 0
M  V30 9 C 0.5145 4.4683 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(3 1 5 6) ATTACH=ANY
M  V30 8 1 8 9
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(mol1);

  auto mol2 = R"CTAB(
  Mrv2007 06242006032D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7083 2.415 0 0
M  V30 2 C -3.042 1.645 0 0
M  V30 3 C -3.042 0.105 0 0
M  V30 4 N -1.7083 -0.665 0 0
M  V30 5 C -0.3747 0.105 0 0
M  V30 6 C -0.3747 1.645 0 0
M  V30 7 * -3.042 0.875 0 0
M  V30 8 F -5.0434 0.875 0 0
M  V30 9 * -1.0415 2.03 0 0
M  V30 10 Cl -1.0415 4.34 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(2 2 3) ATTACH=ANY
M  V30 8 1 9 10 ENDPTS=(2 1 6) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
  REQUIRE(mol2);

  SECTION("PositionVariationOp unit tests 1") {
    MolEnumerator::PositionVariationOp op(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == 3);

    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "COc1ccccn1");
    }
    {
      std::vector<size_t> elems{2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "COc1cccnc1");
    }
  }
  SECTION("PositionVariationOp unit tests 2") {
    MolEnumerator::PositionVariationOp op(*mol2);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == 2);
    CHECK(vcnts[1] == 2);

    {
      std::vector<size_t> elems{1, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "Fc1cc(Cl)ccn1");
    }
    {
      std::vector<size_t> elems{0, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "Fc1cncc(Cl)c1");
    }
  }

  SECTION("PositionVariationOp unit tests 3") {
    MolEnumerator::PositionVariationOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == 3);
  }

  SECTION("enumeration basics 1") {
    MolEnumerator::MolEnumeratorParams ps;
    ps.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
        new MolEnumerator::PositionVariationOp());
    auto bundle = MolEnumerator::enumerate(*mol1, ps);
    CHECK(bundle.size() == 3);

    CHECK(bundle.getMols()[0]->getAtomWithIdx(0)->getDegree() == 3);
    CHECK(bundle.getMols()[0]->getAtomWithIdx(0)->getImplicitValence() == 0);

    std::vector<std::string> tsmis = {"COc1ccncc1", "COc1ccccn1", "COc1cccnc1"};
    std::vector<std::string> smis;
    for (const auto &molp : bundle.getMols()) {
      smis.push_back(MolToSmiles(*molp));
    }
    CHECK(smis == tsmis);
  }
  SECTION("enumeration basics 2") {
    MolEnumerator::MolEnumeratorParams ps;
    ps.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
        new MolEnumerator::PositionVariationOp());
    auto bundle = MolEnumerator::enumerate(*mol2, ps);
    CHECK(bundle.size() == 4);
    std::vector<std::string> tsmis = {"Fc1cnccc1Cl", "Fc1cncc(Cl)c1",
                                      "Fc1cc(Cl)ccn1", "Fc1ccc(Cl)cn1"};
    std::vector<std::string> smis;
    for (const auto &molp : bundle.getMols()) {
      smis.push_back(MolToSmiles(*molp));
    }
    CHECK(smis == tsmis);
  }
}

TEST_CASE("LINKNODE", "[MolEnumerator]") {
  auto mol1 = R"CTAB(one linknode
  Mrv2007 06222005102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 6.9164 14.4947 0 0
M  V30 4 C 9.5836 14.4947 0 0
M  V30 5 C 9.5836 12.9547 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 END BOND
M  V30 LINKNODE 1 4 2 1 2 1 5
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(mol1);
  SECTION("LinkNode unit tests 1") {
    MolEnumerator::LinkNodeOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == 4);

    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "OC1CCCC1");
    }
    {
      std::vector<size_t> elems{2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "OC1CCCCC(O)C1O");
    }
  }
  SECTION("LinkNode unit tests 2") {
    auto mol2 = R"CTAB(two linknodes
  Mrv2014 07072016412D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 7.2366 14.4611 0 0
M  V30 4 C 8.7681 14.622 0 0
M  V30 5 C 9.3945 13.2151 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 7 F 9.5382 15.9557 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 7 1 4 7
M  V30 END BOND
M  V30 LINKNODE 1 3 2 1 2 1 5
M  V30 LINKNODE 1 4 2 4 3 4 5
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol2);
    MolEnumerator::LinkNodeOp op;
    op.initFromMol(*mol2);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == 3);
    CHECK(vcnts[1] == 4);

    {
      std::vector<size_t> elems{0, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "OC1CCC(F)C1");
    }
    {
      std::vector<size_t> elems{2, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "OC1CCC(F)CC(O)C1O");
    }
    {
      std::vector<size_t> elems{1, 2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "OC1CCC(F)C(F)C(F)CC1O");
    }
  }
  SECTION("LinkNode unit tests 3") {
    auto mol = R"CTAB(neighboring linknodes
  Mrv2014 07082008052D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 7.2366 14.4611 0 0
M  V30 4 C 8.7681 14.622 0 0
M  V30 5 C 9.3945 13.2151 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 7 F 9.5382 15.9557 0 0
M  V30 8 C 9.5837 9.8747 0 0
M  V30 9 C 10.9008 12.8949 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 7 1 4 7
M  V30 8 1 6 8
M  V30 9 1 5 9
M  V30 END BOND
M  V30 LINKNODE 1 3 2 1 2 1 5
M  V30 LINKNODE 1 2 2 5 1 5 4
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    MolEnumerator::LinkNodeOp op;
    op.initFromMol(*mol);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == 3);
    CHECK(vcnts[1] == 2);

    {
      std::vector<size_t> elems{0, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "COC1CCC(F)C1C");
    }
    {
      std::vector<size_t> elems{1, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "COC1CCC(F)C(C)C1OC");
    }
    {
      std::vector<size_t> elems{1, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "COC1CCC(F)C(C)C(C)C1OC");
    }
  }
  SECTION("LinkNode unit tests: isotopes") {
    auto mol = R"CTAB(isotopes and linknodes
  Mrv2014 07082008362D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.5 19.6392 0 0 MASS=14
M  V30 2 C -0.27 18.3054 0 0 MASS=15
M  V30 3 C 1.27 18.3054 0 0 MASS=13
M  V30 4 O 2.6037 17.5354 0 0 MASS=12
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 END BOND
M  V30 LINKNODE 1 3 2 3 1 3 2
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    MolEnumerator::LinkNodeOp op;
    op.initFromMol(*mol);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == 3);

    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "[12OH][13CH]1[14CH2][15CH2][13CH]1[12OH]");
    }
  }
  SECTION("LinkNode unit tests: no linknodes") {
    auto mol = R"CTAB(no linknodes
  Mrv2007 06222005102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 6.9164 14.4947 0 0
M  V30 4 C 9.5836 14.4947 0 0
M  V30 5 C 9.5836 12.9547 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol);
    MolEnumerator::LinkNodeOp op;
    op.initFromMol(*mol);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 0);
  }

  SECTION("enumeration basics 1") {
    MolEnumerator::MolEnumeratorParams ps;
    ps.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
        new MolEnumerator::LinkNodeOp());
    auto bundle = MolEnumerator::enumerate(*mol1, ps);
    std::vector<std::string> tsmis = {"OC1CCCC1", "OC1CCCCC1O",
                                      "OC1CCCCC(O)C1O", "OC1CCCCC(O)C(O)C1O"};
    CHECK(bundle.size() == tsmis.size());
    std::vector<std::string> smis;
    for (const auto &molp : bundle.getMols()) {
      smis.push_back(MolToSmiles(*molp));
    }
    CHECK(smis == tsmis);
  }
  SECTION("enumeration basics 2") {
    auto mol = R"CTAB(no linknodes
  Mrv2007 06222005102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 6.9164 14.4947 0 0
M  V30 4 C 9.5836 14.4947 0 0
M  V30 5 C 9.5836 12.9547 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol);
    MolEnumerator::MolEnumeratorParams ps;
    ps.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
        new MolEnumerator::LinkNodeOp());
    auto bundle = MolEnumerator::enumerate(*mol, ps);
    CHECK(bundle.size() == 0);
  }
}

TEST_CASE("multiple enumeration points", "[MolEnumerator]") {
  auto mol1 = R"CTAB(
  Mrv2014 12212013392D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 12 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6 4.7484 0 0
M  V30 2 C -7.3337 3.9784 0 0
M  V30 3 N -7.3337 2.4383 0 0
M  V30 4 C -6 1.6683 0 0
M  V30 5 C -4.6663 2.4383 0 0
M  V30 6 C -4.6663 3.9784 0 0
M  V30 7 C -5.8773 0.0617 0 0
M  V30 8 C -3.2136 0.2013 0 0
M  V30 9 C -4.5052 -0.6374 0 0
M  V30 10 O -4.4246 -2.1753 0 0
M  V30 11 * -6 4.235 0 0
M  V30 12 C -4.845 6.2355 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 8 9
M  V30 8 1 7 9
M  V30 9 1 5 8
M  V30 10 1 7 4
M  V30 11 1 9 10
M  V30 12 1 11 12 ENDPTS=(3 1 2 6) ATTACH=ANY
M  V30 END BOND
M  V30 LINKNODE 1 3 2 9 7 9 8
M  V30 END CTAB
M  END
)CTAB"_ctab;
  REQUIRE(mol1);
  std::vector<std::string> tsmis = {
      "Cc1cnc2c(c1)CC(O)C2",         "Cc1cnc2c(c1)CC(O)C(O)C2",
      "Cc1cnc2c(c1)CC(O)C(O)C(O)C2", "Cc1ccc2c(n1)CC(O)C2",
      "Cc1ccc2c(n1)CC(O)C(O)C2",     "Cc1ccc2c(n1)CC(O)C(O)C(O)C2",
      "Cc1ccnc2c1CC(O)C2",           "Cc1ccnc2c1CC(O)C(O)C2",
      "Cc1ccnc2c1CC(O)C(O)C(O)C2"};
  SECTION("test1") {
    std::vector<MolEnumerator::MolEnumeratorParams> paramsList;
    MolEnumerator::MolEnumeratorParams posVariationParams;
    posVariationParams.dp_operation =
        std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
            new MolEnumerator::PositionVariationOp());
    paramsList.push_back(posVariationParams);
    MolEnumerator::MolEnumeratorParams linkParams;
    linkParams.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
        new MolEnumerator::LinkNodeOp());
    paramsList.push_back(linkParams);
    auto bundle = MolEnumerator::enumerate(*mol1, paramsList);
    CHECK(bundle.size() == tsmis.size());
    for (const auto &molp : bundle.getMols()) {
      auto smi = MolToSmiles(*molp);
      // std::cerr << smi << std::endl;
      CHECK(std::find(tsmis.begin(), tsmis.end(), smi) != tsmis.end());
    }
  }
  SECTION("test2") {
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == tsmis.size());
    for (const auto &molp : bundle.getMols()) {
      auto smi = MolToSmiles(*molp);
      CHECK(std::find(tsmis.begin(), tsmis.end(), smi) != tsmis.end());
    }
  }
  SECTION("edges1") {
    std::vector<MolEnumerator::MolEnumeratorParams> paramsList;
    auto bundle = MolEnumerator::enumerate(*mol1, paramsList);
    CHECK(bundle.size() == 0);
  }
  SECTION("edges2") {
    std::vector<MolEnumerator::MolEnumeratorParams> paramsList;
    MolEnumerator::MolEnumeratorParams posVariationParams;
    posVariationParams.dp_operation =
        std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
            new MolEnumerator::PositionVariationOp());
    paramsList.push_back(posVariationParams);
    MolEnumerator::MolEnumeratorParams linkParams;
    linkParams.dp_operation = std::shared_ptr<MolEnumerator::MolEnumeratorOp>(
        new MolEnumerator::LinkNodeOp());
    paramsList.push_back(linkParams);
    auto mol = "c1ccccc1"_smiles;
    auto bundle = MolEnumerator::enumerate(*mol, paramsList);
    CHECK(bundle.size() == 0);
  }
}

TEST_CASE("multiple enumeration points 2", "[MolEnumerator][bug]") {
  auto mol1 = R"CTAB(
  Mrv2014 02182116082D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 15 14 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.1474 10.0453 0 0
M  V30 2 C -3.481 9.2753 0 0
M  V30 3 C -3.481 7.7352 0 0
M  V30 4 C -2.1474 6.9652 0 0
M  V30 5 C -0.8137 7.7352 0 0
M  V30 6 C -0.8137 9.2753 0 0
M  V30 7 N 0.6334 9.802 0 0
M  V30 8 C 0.7204 7.3118 0 0
M  V30 9 C 1.5815 8.5885 0 0
M  V30 10 * -3.0365 9.0186 0 0
M  V30 11 O -3.0365 11.3286 0 0
M  V30 12 C 1.0579 11.2824 0 0
M  V30 13 O -0.0119 12.3901 0 0
M  V30 14 * 1.151 7.9501 0 0
M  V30 15 C 1.151 10.2601 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 9
M  V30 8 1 7 6
M  V30 9 1 5 8
M  V30 10 2 8 9
M  V30 11 1 10 11 ENDPTS=(3 2 3 1) ATTACH=ANY
M  V30 12 1 7 12
M  V30 13 1 12 13
M  V30 14 1 14 15 ENDPTS=(2 8 9) ATTACH=ANY
M  V30 END BOND
M  V30 LINKNODE 1 2 2 12 7 12 13
M  V30 END CTAB
M  END
)CTAB"_ctab;
  std::vector<std::string> tsmis = {
      "Cc1cn(CO)c2cc(O)ccc12", "Cc1cn(CCO)c2cc(O)ccc12",
      "Cc1cc2ccc(O)cc2n1CO",   "Cc1cc2ccc(O)cc2n1CCO",
      "Cc1cn(CO)c2ccc(O)cc12", "Cc1cn(CCO)c2ccc(O)cc12",
      "Cc1cc2cc(O)ccc2n1CO",   "Cc1cc2cc(O)ccc2n1CCO",
      "Cc1cn(CO)c2c(O)cccc12", "Cc1cn(CCO)c2c(O)cccc12",
      "Cc1cc2cccc(O)c2n1CO",   "Cc1cc2cccc(O)c2n1CCO"};
  SECTION("test2") {
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == tsmis.size());
    for (const auto &molp : bundle.getMols()) {
      auto smi = MolToSmiles(*molp);
      // std::cerr << "\"" << smi << "\"," << std::endl;
      CHECK(std::find(tsmis.begin(), tsmis.end(), smi) != tsmis.end());
    }
  }
}