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

    CHECK(bundle.getMols()[0]->getAtomWithIdx(0)->getDegree()==3);
    CHECK(bundle.getMols()[0]->getAtomWithIdx(0)->getImplicitValence()==0);    

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
