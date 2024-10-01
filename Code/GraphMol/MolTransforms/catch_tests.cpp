//
//  Copyright (c) 2024 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
///
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <algorithm>
#include <regex>

using namespace RDKit;
using std::unique_ptr;

TEST_CASE("chirality flip in canonicalizeMolConformer") {
  SECTION("first discovery") {
    auto m1 = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.533948 2.465972 -1.220733 0
M  V30 2 C -1.361222 1.165424 -1.939720 0 CFG=2
M  V30 3 C -2.694145 0.606596 -2.153754 0
M  V30 4 N -3.766061 0.178135 -2.319155 0
M  V30 5 C -0.566772 0.289201 -1.072783 0
M  V30 6 C 0.084234 -0.420637 -0.345400 0
M  V30 7 C 0.867082 -1.237691 0.508367 0
M  V30 8 C 1.549171 -1.919699 1.232613 0
M  V30 9 C 2.339643 -2.731169 2.106347 0
M  V30 10 C 2.993349 -3.411935 2.854815 0
M  V30 11 O 3.747921 -4.178375 3.716269 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 2 3
M  V30 3 3 3 4
M  V30 4 1 2 5
M  V30 5 3 5 6
M  V30 6 1 6 7
M  V30 7 3 7 8
M  V30 8 1 8 9
M  V30 9 3 9 10
M  V30 10 1 10 11
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m1);
    auto m2 = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 10 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.278200 0.712725 -0.066592 0
M  V30 2 C -2.682025 -0.648154 -0.271665 0 CFG=1
M  V30 3 C -2.915612 -1.059871 -1.656844 0
M  V30 4 N -3.080403 -1.390524 -2.739848 0
M  V30 5 C -1.232937 -0.539032 -0.023304 0
M  V30 6 C -0.029030 -0.463778 0.167869 0
M  V30 7 C 1.366183 -0.398546 0.377840 0
M  V30 8 C 2.562372 -0.361417 0.548269 0
M  V30 9 C 3.973502 -0.322846 0.772691 0
M  V30 10 C 5.164230 -0.269523 0.960648 0
M  V30 11 O 6.516239 -0.218493 1.152837 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 3 3 4
M  V30 4 1 2 5
M  V30 5 3 5 6
M  V30 6 1 6 7
M  V30 7 3 7 8
M  V30 8 1 8 9
M  V30 9 3 9 10
M  V30 10 1 10 11
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m2);

    auto m1_cp = ROMol(*m1);
    auto m2_cp = ROMol(*m2);
    MolTransforms::canonicalizeConformer(m1_cp.getConformer());
    MolTransforms::canonicalizeConformer(m2_cp.getConformer());

    MolOps::assignChiralTypesFrom3D(m1_cp);
    MolOps::assignChiralTypesFrom3D(m2_cp);
    CHECK(m1_cp.getAtomWithIdx(1)->getChiralTag() ==
          m1->getAtomWithIdx(1)->getChiralTag());
    CHECK(m2_cp.getAtomWithIdx(1)->getChiralTag() ==
          m2->getAtomWithIdx(1)->getChiralTag());
  }
}
