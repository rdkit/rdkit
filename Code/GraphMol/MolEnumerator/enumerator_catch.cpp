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

TEST_CASE(
    "github #4382: Enumerate fails on variable attachment points with queries",
    "[MolEnumerator][bug]") {
  SECTION("test1") {
    auto mol1 = R"CTAB(
  Mrv2108 08032115452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.083 5.5632 0 0
M  V30 2 C -1.083 7.1033 0 0
M  V30 3 N 0.3973 7.5278 0 0
M  V30 4 N 0.3104 5.0376 0 0
M  V30 5 C 1.2585 6.251 0 0
M  V30 6 * 0.3539 6.2827 0 0
M  V30 7 A 1.5089 8.2832 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 2
M  V30 3 1 1 4
M  V30 4 1 4 5
M  V30 5 2 3 5
M  V30 6 1 6 7 ENDPTS=(2 4 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == 2);
    std::vector<std::string> tsmas = {"[#6]1-[#6]-[#7]=[#6]-[#7]-1-[!#1]",
                                      "[#6]1-[#6]-[#7](=[#6]-[#7]-1)-[!#1]"};
    for (const auto &molp : bundle.getMols()) {
      auto smarts = MolToSmarts(*molp);
      CHECK(std::find(tsmas.begin(), tsmas.end(), smarts) != tsmas.end());
    }
  }
  SECTION("test1 reversed") {
    // never actually observed one of these, but it should still work
    auto mol1 = R"CTAB(
  Mrv2108 08032115452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.083 5.5632 0 0
M  V30 2 C -1.083 7.1033 0 0
M  V30 3 N 0.3973 7.5278 0 0
M  V30 4 N 0.3104 5.0376 0 0
M  V30 5 C 1.2585 6.251 0 0
M  V30 6 A 1.5089 8.2832 0 0
M  V30 7 * 0.3539 6.2827 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 2
M  V30 3 1 1 4
M  V30 4 1 4 5
M  V30 5 2 3 5
M  V30 6 1 6 7 ENDPTS=(2 4 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == 2);
    std::vector<std::string> tsmas = {"[#6]1-[#6]-[#7]=[#6]-[#7]-1-[!#1]",
                                      "[#6]1-[#6]-[#7](=[#6]-[#7]-1)-[!#1]"};
    for (const auto &molp : bundle.getMols()) {
      auto smarts = MolToSmarts(*molp);
      CHECK(std::find(tsmas.begin(), tsmas.end(), smarts) != tsmas.end());
    }
  }

  SECTION("test3: both are dummy atoms") {
    auto mol1 = R"CTAB(
  Mrv2108 08032115452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.083 5.5632 0 0
M  V30 2 C -1.083 7.1033 0 0
M  V30 3 N 0.3973 7.5278 0 0
M  V30 4 N 0.3104 5.0376 0 0
M  V30 5 C 1.2585 6.251 0 0
M  V30 6 * 0.3539 6.2827 0 0
M  V30 7 * 1.5089 8.2832 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 2
M  V30 3 1 1 4
M  V30 4 1 4 5
M  V30 5 2 3 5
M  V30 6 1 6 7 ENDPTS=(2 4 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == 2);
    std::vector<std::string> tsmas = {"[#6]1-[#6]-[#7]=[#6]-[#7]-1-*",
                                      "[#6]1-[#6]-[#7](=[#6]-[#7]-1)-*"};
    for (const auto &molp : bundle.getMols()) {
      auto smarts = MolToSmarts(*molp);
      CHECK(std::find(tsmas.begin(), tsmas.end(), smarts) != tsmas.end());
    }
  }
}

TEST_CASE("github #4381: need implicit H cleanup after Enumerate",
          "[MolEnumerator][bug]") {
  SECTION("test1") {
    auto mol1 = R"CTAB(
  Mrv2108 08032115452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.083 5.5632 0 0
M  V30 2 C -1.083 7.1033 0 0
M  V30 3 N 0.3973 7.5278 0 0
M  V30 4 N 0.3104 5.0376 0 0
M  V30 5 C 1.2585 6.251 0 0
M  V30 6 * 0.3539 6.2827 0 0
M  V30 7 C 1.5089 8.2832 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 3 2
M  V30 3 1 1 4
M  V30 4 1 4 5
M  V30 5 2 3 5
M  V30 6 1 6 7 ENDPTS=(2 4 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == 2);
    std::vector<std::string> tsmas = {"[#6]1:[#6]:[#7]:[#6]:[#7]:1-[#6]",
                                      "[#6]1:[#6]:[#7](:[#6]:[#7]:1)-[#6]"};
    for (const auto &molp : bundle.getMols()) {
      auto smarts = MolToSmarts(*molp);
      CHECK(std::find(tsmas.begin(), tsmas.end(), smarts) != tsmas.end());
    }
  }
}

TEST_CASE("RepeatUnit", "[MolEnumerator]") {
  SECTION("basic HT enumeration") {
    auto mol1 = R"CTAB(
  ACCLDraw05132106232D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 9.7578 -7.0211 0 0 
M  V30 2 O 10.7757 -7.6201 0 0 
M  V30 3 * 11.8037 -7.0378 0 0 
M  V30 4 * 8.7298 -7.6034 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 1 4 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(2 2 1) XBONDS=(2 3 2) BRKXYZ=(9 9.24 -7.9 0 9.24 -6.72 -
M  V30 0 0 0 0) BRKXYZ=(9 11.29 -6.74 0 11.29 -7.92 0 0 0 0) CONNECT=HT -
M  V30 LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CO*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*COCOCO*");
    }
  }
  SECTION("basic HH enumeration") {
    auto mol1 = R"CTAB(
  ACCLDraw05132106232D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 9.7578 -7.0211 0 0 
M  V30 2 O 10.7757 -7.6201 0 0 
M  V30 3 * 11.8037 -7.0378 0 0 
M  V30 4 * 8.7298 -7.6034 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 1 4 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(2 2 1) XBONDS=(2 3 2) BRKXYZ=(9 9.24 -7.9 0 9.24 -6.72 -
M  V30 0 0 0 0) BRKXYZ=(9 11.29 -6.74 0 11.29 -7.92 0 0 0 0) CONNECT=HH -
M  V30 LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*COOCCO*");
    }
  }
  SECTION("EU enumeration (which we don't do)") {
    auto mol1 = R"CTAB(
  ACCLDraw05132106232D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 9.7578 -7.0211 0 0 
M  V30 2 O 10.7757 -7.6201 0 0 
M  V30 3 * 11.8037 -7.0378 0 0 
M  V30 4 * 8.7298 -7.6034 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 1 4 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(2 2 1) XBONDS=(2 3 2) BRKXYZ=(9 9.24 -7.9 0 9.24 -6.72 -
M  V30 0 0 0 0) BRKXYZ=(9 11.29 -6.74 0 11.29 -7.92 0 0 0 0) CONNECT=EU -
M  V30 LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.empty());
  }

  SECTION("non-overlapping multi-HT enumeration") {
    auto mol1 = R"CTAB(
  Mrv2108 09232105582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 12.719 -9.1518 0 0
M  V30 2 O 14.0458 -9.9326 0 0
M  V30 3 C 15.3857 -9.1735 0 0 MASS=13
M  V30 4 * 11.379 -9.9108 0 0
M  V30 5 N 16.713 -9.9545 0 0
M  V30 6 C 18.053 -9.1955 0 0
M  V30 7 * 19.3803 -9.9764 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 4 1 3 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(2 2 1) XBONDS=(2 2 3) BRKXYZ=(9 12.044 -10.2974 0 -
M  V30 12.044 -8.7593 0 0 0 0) BRKXYZ=(9 14.7161 -8.7854 0 14.7161 -10.3235 0 -
M  V30 0 0 0) CONNECT=HT LABEL=n
M  V30 2 SRU 0 ATOMS=(2 5 6) XBONDS=(2 4 6) BRKXYZ=(9 19.0681 -8.7207 0 -
M  V30 18.131 -10.3134 0 0 0 0) BRKXYZ=(9 15.6979 -10.4293 0 16.6351 -8.8365 -
M  V30 0 0 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    CHECK(vcnts[1] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*[13CH2]*");
    }
    {
      std::vector<size_t> elems{1, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CO[13CH2]*");
    }
    {
      std::vector<size_t> elems{0, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]*");
    }
    {
      std::vector<size_t> elems{1, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]OC*");
    }
    {
      std::vector<size_t> elems{3, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]OCOCOC*");
    }
    {
      std::vector<size_t> elems{1, 3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CNCNCN[13CH2]OC*");
    }
    {
      std::vector<size_t> elems{2, 2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CNCN[13CH2]OCOC*");
    }
  }

  SECTION("non-overlapping multi-HT enumeration reversed example") {
    // same as the previous example, but the head/tail definition of the second
    // SRU is swapped. This shouldn't change the results
    auto mol1 = R"CTAB(
  Mrv2108 09232105582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 12.719 -9.1518 0 0
M  V30 2 O 14.0458 -9.9326 0 0
M  V30 3 C 15.3857 -9.1735 0 0 MASS=13
M  V30 4 * 11.379 -9.9108 0 0
M  V30 5 N 16.713 -9.9545 0 0
M  V30 6 C 18.053 -9.1955 0 0
M  V30 7 * 19.3803 -9.9764 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 4 1 3 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(2 2 1) XBONDS=(2 2 3) BRKXYZ=(9 12.044 -10.2974 0 -
M  V30 12.044 -8.7593 0 0 0 0) BRKXYZ=(9 14.7161 -8.7854 0 14.7161 -10.3235 0 -
M  V30 0 0 0) CONNECT=HT LABEL=n
M  V30 2 SRU 0 ATOMS=(2 6 5) XBONDS=(2 6 4) BRKXYZ=(9 19.0681 -8.7207 0 -
M  V30 18.131 -10.3134 0 0 0 0) BRKXYZ=(9 15.6979 -10.4293 0 16.6351 -8.8365 -
M  V30 0 0 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    CHECK(vcnts[1] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*[13CH2]*");
    }
    {
      std::vector<size_t> elems{1, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CO[13CH2]*");
    }
    {
      std::vector<size_t> elems{0, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]*");
    }
    {
      std::vector<size_t> elems{1, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]OC*");
    }
    {
      std::vector<size_t> elems{3, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]OCOCOC*");
    }
    {
      std::vector<size_t> elems{1, 3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CNCNCN[13CH2]OC*");
    }
    {
      std::vector<size_t> elems{2, 2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CNCN[13CH2]OCOC*");
    }
  }

  SECTION("nested HT enumeration is not supported") {
    auto mol1 = R"CTAB(
  ACCLDraw05132106342D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 2 0 0
M  V30 BEGIN ATOM
M  V30 1 * 9.2115 -7.4169 0 0 
M  V30 2 C 10.2483 -6.8511 0 0 CFG=3 
M  V30 3 C 10.2767 -5.67 0 0 
M  V30 4 O 11.257 -7.4662 0 0 
M  V30 5 * 12.2941 -6.9003 0 0 
M  V30 6 N 9.268 -5.0548 0 0 CFG=3 
M  V30 7 * 9.2964 -3.8737 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 2 4 
M  V30 4 1 4 5 
M  V30 5 1 3 6 
M  V30 6 1 6 7 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(5 3 2 4 6 7) XBONDS=(2 1 4) BRKXYZ=(9 9.72 -7.72 0 9.74 -
M  V30 -6.54 0 0 0 0) BRKXYZ=(9 11.79 -6.59 0 11.76 -7.77 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n 
M  V30 2 SRU 2 ATOMS=(2 3 6) XBONDS=(2 6 2) BRKXYZ=(9 8.78 -4.78 0 9.78 -4.15 -
M  V30 0 0 0 0) BRKXYZ=(9 10.76 -5.95 0 9.76 -6.57 0 0 0 0) CONNECT=HT -
M  V30 LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    CHECK_THROWS_AS(op.initFromMol(*mol1), ValueErrorException);
  }

  SECTION("single atom SRU enumeration") {
    auto mol1 = R"CTAB(
  Mrv2108 09232115452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 1 0 0
M  V30 BEGIN ATOM
M  V30 1 * 9.3333 -3.9167 0 0
M  V30 2 C 10.667 -3.1467 0 0
M  V30 3 * 12.0007 -3.9167 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(1 2) XBONDS=(2 1 2) BRKXYZ=(9 11.6782 -2.6635 0 10.7542 -
M  V30 -4.2639 0 0 0 0) BRKXYZ=(9 10.5799 -4.2639 0 9.6559 -2.6635 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*C*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CCC*");
    }
  }
}

TEST_CASE("SRUs directly bound to each other") {
  SECTION("basics") {
    auto mol1 = R"CTAB(
  Mrv2108 09242104182D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 2 0 0
M  V30 BEGIN ATOM
M  V30 1 * -0.75 -3.2917 0 0
M  V30 2 C 0.5837 -2.5217 0 0
M  V30 3 O 1.9174 -3.2917 0 0
M  V30 4 C 3.251 -2.5217 0 0
M  V30 5 N 4.5847 -3.2917 0 0
M  V30 6 * 5.9184 -2.5217 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(2 2 3) XBONDS=(2 1 3) BRKXYZ=(9 2.0045 -2.1744 0 2.9285 -
M  V30 -3.7748 0 0 0 0) BRKXYZ=(9 0.4965 -3.6389 0 -0.4275 -2.0385 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n
M  V30 2 SRU 0 ATOMS=(2 4 5) XBONDS=(2 3 5) BRKXYZ=(9 4.6719 -2.1744 0 5.5959 -
M  V30 -3.7748 0 0 0 0) BRKXYZ=(9 3.1639 -3.6389 0 2.2399 -2.0385 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    CHECK(vcnts[1] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      // FIX: this result is ugly and is fixable by adding
      // some extra bookkeeping. Somehow doesn't seem
      // super important at the moment though.
      CHECK(MolToSmiles(*newmol) == "*.*");
    }
    {
      std::vector<size_t> elems{1, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CO*");
    }
    {
      std::vector<size_t> elems{0, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN*");
    }
    {
      std::vector<size_t> elems{1, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*COCN*");
    }
    {
      std::vector<size_t> elems{2, 2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*COCOCNCN*");
    }
  }
}

TEST_CASE("RepeatUnit with enumerate()", "[MolEnumerator]") {
  SECTION("basic HT enumeration") {
    auto mol1 = R"CTAB(
  ACCLDraw05132106232D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 9.7578 -7.0211 0 0 
M  V30 2 O 10.7757 -7.6201 0 0 
M  V30 3 * 11.8037 -7.0378 0 0 
M  V30 4 * 8.7298 -7.6034 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 1 4 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(2 2 1) XBONDS=(2 3 2) BRKXYZ=(9 9.24 -7.9 0 9.24 -6.72 -
M  V30 0 0 0 0) BRKXYZ=(9 11.29 -6.74 0 11.29 -7.92 0 0 0 0) CONNECT=HT -
M  V30 LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;

    std::vector<std::string> expected = {"**", "*CO*", "*COCO*", "*COCOCO*"};
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == op.d_defaultRepeatCount);
    for (const auto &molp : bundle.getMols()) {
      auto smiles = MolToSmiles(*molp);
      CHECK(std::find(expected.begin(), expected.end(), smiles) !=
            expected.end());
    }
  }
}

TEST_CASE("RepeatUnit+LINKNODE", "[MolEnumerator]") {
  SECTION("basics") {
    auto mol1 = R"CTAB(
  Mrv2108 09252115192D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 12.719 -9.1518 0 0
M  V30 2 O 14.0458 -9.9326 0 0
M  V30 3 * 15.3857 -9.1735 0 0
M  V30 4 * 11.379 -9.9108 0 0
M  V30 5 C 12.7317 -7.6118 0 0
M  V30 6 C 12.2558 -6.1472 0 0
M  V30 7 C 13.7622 -6.4674 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 5 1 7 6
M  V30 6 1 6 5
M  V30 7 1 5 7
M  V30 END BOND
M  V30 LINKNODE 1 2 2 7 5 7 6
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(5 2 1 5 7 6) XBONDS=(2 2 3) BRKXYZ=(9 12.044 -10.2974 0 -
M  V30 12.044 -8.7593 0 0 0 0) BRKXYZ=(9 14.7161 -8.7854 0 14.7161 -10.3235 0 -
M  V30 0 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    CHECK_THROWS_AS(MolEnumerator::enumerate(*mol1), ValueErrorException);
  }
}

TEST_CASE("ladder") {
  SECTION("basics, no XBCORR") {
    auto mol1 = R"CTAB(
  Mrv2108 09262105082D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 11.375 -0.2516 0 0
M  V30 2 * 10.0413 -1.0216 0 0
M  V30 3 * 10.0413 -2.5617 0 0
M  V30 4 C 11.375 -3.3317 0 0
M  V30 5 C 12.7087 -2.5617 0 0
M  V30 6 N 12.7087 -1.0216 0 0
M  V30 7 * 14.0423 -0.2516 0 0
M  V30 8 * 14.0423 -3.3317 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 4
M  V30 3 1 4 5
M  V30 4 1 5 6
M  V30 5 1 1 6
M  V30 6 1 5 8
M  V30 7 1 7 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(4 1 4 5 6) XBONDS=(4 1 2 6 7) BRKXYZ=(9 13.2958 0.0956 -
M  V30 0 13.2958 -3.679 0 0 0 0) -
M  V30 BRKXYZ=(9 10.8638 -3.8148 0 10.8638 0.2315 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);

    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "**.**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC(*)N(*)C*");
    }
    {
      std::vector<size_t> elems{2};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CN(*)C(*)CN1C*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CN2CC(*)N(*)CC2CN1C*");
    }
  }

  SECTION("basics, with XBCORR") {
    auto mol1 = R"CTAB(
  Mrv2108 09262105082D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 11.375 -0.2516 0 0
M  V30 2 * 10.0413 -1.0216 0 0
M  V30 3 * 10.0413 -2.5617 0 0
M  V30 4 C 11.375 -3.3317 0 0
M  V30 5 C 12.7087 -2.5617 0 0
M  V30 6 N 12.7087 -1.0216 0 0
M  V30 7 * 14.0423 -0.2516 0 0
M  V30 8 * 14.0423 -3.3317 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 4
M  V30 3 1 4 5
M  V30 4 1 5 6
M  V30 5 1 1 6
M  V30 6 1 5 8
M  V30 7 1 7 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(4 1 4 5 6) XBONDS=(4 1 2 6 7) XBCORR=(4 1 7 2 6) -
M  V30 BRKXYZ=(9 13.2958 0.0956 0 13.2958 -3.679 0 0 0 0) -
M  V30 BRKXYZ=(9 10.8638 -3.8148 0 10.8638 0.2315 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);

    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "**.**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC(*)N(*)C*");
    }
    {
      std::vector<size_t> elems{2};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CC(*)N(*)CN1C*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CC2CC(*)N(*)CN2CN1C*");
    }
  }
  SECTION("using enumerate with XBCORR") {
    auto mol1 = R"CTAB(
  Mrv2108 09262105082D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 11.375 -0.2516 0 0
M  V30 2 * 10.0413 -1.0216 0 0
M  V30 3 * 10.0413 -2.5617 0 0
M  V30 4 C 11.375 -3.3317 0 0
M  V30 5 C 12.7087 -2.5617 0 0
M  V30 6 N 12.7087 -1.0216 0 0
M  V30 7 * 14.0423 -0.2516 0 0
M  V30 8 * 14.0423 -3.3317 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 4
M  V30 3 1 4 5
M  V30 4 1 5 6
M  V30 5 1 1 6
M  V30 6 1 5 8
M  V30 7 1 7 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(4 1 4 5 6) XBONDS=(4 1 2 6 7) XBCORR=(4 1 7 2 6) -
M  V30 0 BRKXYZ=(9 13.2958 0.0956 13.2958 -3.679 0 0 0 0) -
M  V30 BRKXYZ=(9 10.8638 -3.8148 0 10.8638 0.2315 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;

    std::vector<std::string> expected = {"**.**", "*CC(*)N(*)C*",
                                         "*CC1CC(*)N(*)CN1C*",
                                         "*CC1CC2CC(*)N(*)CN2CN1C*"};
    auto bundle = MolEnumerator::enumerate(*mol1);
    CHECK(bundle.size() == op.d_defaultRepeatCount);
    for (const auto &molp : bundle.getMols()) {
      auto smiles = MolToSmiles(*molp);
      CHECK(std::find(expected.begin(), expected.end(), smiles) !=
            expected.end());
    }
  }
}

TEST_CASE("starting from CXSMILES") {
  SECTION("basic HT enumeration") {
    auto mol1 = "*-CO-* |$star_e;;;star_e$,Sg:n:2,1::ht|"_smiles;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CO*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*COCOCO*");
    }
  }
  SECTION("basic HH enumeration") {
    auto mol1 = "*-CO-* |$star_e;;;star_e$,Sg:n:2,1::hh|"_smiles;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*COOCCO*");
    }
  }

  SECTION("non-overlapping multi-HT enumeration reversed example") {
    auto mol1 =
        "*-CN[13CH2]OC-* |$star_e;;;;;;star_e$,Sg:n:4,5::ht,Sg:n:1,2::ht|"_smiles;
    REQUIRE(mol1);
    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 2);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    CHECK(vcnts[1] == op.d_defaultRepeatCount);

    {
      std::vector<size_t> elems{0, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*[13CH2]*");
    }
    {
      std::vector<size_t> elems{1, 0};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CO[13CH2]*");
    }
    {
      std::vector<size_t> elems{0, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]*");
    }
    {
      std::vector<size_t> elems{1, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]OC*");
    }
    {
      std::vector<size_t> elems{3, 1};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CN[13CH2]OCOCOC*");
    }
    {
      std::vector<size_t> elems{1, 3};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CNCNCN[13CH2]OC*");
    }
    {
      std::vector<size_t> elems{2, 2};
      std::unique_ptr<ROMol> newmol(op(elems));
      CHECK(MolToSmiles(*newmol) == "*CNCN[13CH2]OCOC*");
    }
  }
  SECTION("ladder basics, no XBCORR") {
    auto mol1 =
        "*-CC(-*)N(-*)C-* |$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::ht:4,2:0,6|"_smiles;
    REQUIRE(mol1);

    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "**.**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC(*)N(*)C*");
    }
    {
      std::vector<size_t> elems{2};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CN(*)C(*)CN1C*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CN2CC(*)N(*)CC2CN1C*");
    }
  }

  SECTION("ladder basics, with XBCORR") {
    auto mol1 =
        "*-CC(-*)N(-*)C-* |$star_e;;;star_e;;star_e;;star_e$,Sg:n:6,1,2,4::ht:6,0:4,2|"_smiles;
    REQUIRE(mol1);

    MolEnumerator::RepeatUnitOp op;
    op.initFromMol(*mol1);
    auto vcnts = op.getVariationCounts();
    REQUIRE(vcnts.size() == 1);
    CHECK(vcnts[0] == op.d_defaultRepeatCount);
    {
      std::vector<size_t> elems{0};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "**.**");
    }
    {
      std::vector<size_t> elems{1};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC(*)N(*)C*");
    }
    {
      std::vector<size_t> elems{2};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CC(*)N(*)CN1C*");
    }
    {
      std::vector<size_t> elems{3};
      std::unique_ptr<ROMol> newmol(op(elems));
      REQUIRE(newmol);
      CHECK(MolToSmiles(*newmol) == "*CC1CC2CC(*)N(*)CN2CN1C*");
    }
  }
}

TEST_CASE(
    "Github #6014: MolEnumerator is not propagating query information to molecules ") {
  SECTION("basics") {
    auto m = "O-[CX4]-[CH3] |LN:1:1.5|"_smarts;
    REQUIRE(m);
    auto bundle = MolEnumerator::enumerate(*m);
    CHECK(bundle.size() == 5);
    auto q0 = bundle[0];
    REQUIRE(q0);
    CHECK(q0->getAtomWithIdx(0)->hasQuery());
    CHECK(MolToSmarts(*q0) == "O-[C&X4]-[C&H3]");
    CHECK(MolToSmarts(*bundle[1]) == "O-[C&X4]-[C&X4]-[C&H3]");
  }
}

TEST_CASE(
    "Github #6432: MolEnumerate should clear the reaction properties on its results") {
  SECTION("basics") {
    auto m = "COC |LN:1:1.3|"_smarts;
    REQUIRE(m);
    auto bundle = MolEnumerator::enumerate(*m);
    CHECK(bundle.size() == 3);
    auto q0 = bundle[0];
    REQUIRE(q0);
    for (const auto atom : q0->atoms()) {
      CHECK(!atom->hasProp(common_properties::reactantAtomIdx));
      CHECK(!atom->hasProp(common_properties::reactionMapNum));
      CHECK(!atom->hasProp("was_dummy"));
    }
  }
}
