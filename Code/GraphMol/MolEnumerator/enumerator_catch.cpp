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
  auto mol = R"CTAB(
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
  REQUIRE(mol);

  SECTION("PositionVariationOp unit tests") {
    MolEnumerator::PositionVariationOp op(*mol);
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
}

TEST_CASE("LINKNODE", "[MolEnumerator]") {
  SECTION("basics") {
    auto mol = R"CTAB(
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
    REQUIRE(mol);
  }
}
