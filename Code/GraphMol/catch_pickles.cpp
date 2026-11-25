//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <limits>
#include <fstream>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;
TEST_CASE("Github #6312: space overhead of serializing properties") {
  SECTION("Bonds v2000") {
    auto mol = R"CTAB(
  Mrv1810 02111915042D          

  4  3  0  0  0  0            999 V2000
   -1.5625    1.6071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8480    2.0196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2770    2.0196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5625    0.7821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  2  6  0  0  0  0
  1  4  1  1  0  0  0
M  END)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(1)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 6);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 1);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondStereo) == 1);

    std::string basepkl;
    MolPickler::pickleMol(*mol, basepkl);

    std::string pkl;
    MolPickler::pickleMol(*mol, pkl,
                          PicklerOps::PropertyPickleOptions::BondProps |
                              PicklerOps::PropertyPickleOptions::PrivateProps);
    // std::cerr << "!!!! " << pkl.size() << " " << basepkl.size() << std::endl;

    RWMol mol2(pkl);
    CHECK(mol2.getBondWithIdx(1)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 6);
    CHECK(mol2.getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 1);
    CHECK(mol2.getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondStereo) == 1);
  }
  SECTION("bonds-v3k") {
    auto mol = R"CTAB(
  Mrv1810 02111915102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.9167 3 0 0
M  V30 2 C -1.583 3.77 0 0
M  V30 3 C -4.2503 3.77 0 0
M  V30 4 C -2.9167 1.46 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 6 1 2
M  V30 3 1 1 4 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(1)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 6);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 1);
    CHECK(mol->getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondCfg) == 1);

    std::string basepkl;
    MolPickler::pickleMol(*mol, basepkl);

    std::string pkl;
    MolPickler::pickleMol(*mol, pkl,
                          PicklerOps::PropertyPickleOptions::BondProps |
                              PicklerOps::PropertyPickleOptions::PrivateProps);

    CHECK(pkl.size() > basepkl.size());
    // make sure the property names aren't in the pickle
    CHECK(pkl.find(common_properties::_MolFileBondType) == std::string::npos);
    CHECK(pkl.find(common_properties::_MolFileBondCfg) == std::string::npos);
    // std::cerr << "!!!! " << pkl.size() << " " << basepkl.size() << std::endl;

    RWMol mol2(pkl);
    CHECK(mol2.getBondWithIdx(1)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 6);
    CHECK(mol2.getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondType) == 1);
    CHECK(mol2.getBondWithIdx(2)->getProp<unsigned int>(
              common_properties::_MolFileBondCfg) == 1);
  }
  SECTION("atoms") {
    auto mol = R"CTAB(
  Mrv2211 04272306392D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.375 3.125 0 0 CFG=2
M  V30 2 C -6.0413 3.895 0 0
M  V30 3 O -8.7087 3.895 0 0
M  V30 4 F -7.375 1.585 0 0
M  V30 5 Cl -6.1532 2.1875 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 1 1 4
M  V30 3 1 1 2 CFG=1
M  V30 4 1 1 5 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    std::string basepkl;
    MolPickler::pickleMol(*mol, basepkl);

    CHECK(mol->getAtomWithIdx(0)->getProp<int>(common_properties::molParity) ==
          2);
    CHECK(mol->getAtomWithIdx(0)->getProp<int>(
              common_properties::_ChiralityPossible) == 1);
    mol->getAtomWithIdx(0)->setProp(common_properties::_CIPCode,
                                    std::string("S"));
    mol->getAtomWithIdx(1)->setProp(common_properties::molAtomMapNumber, 1);
    mol->getAtomWithIdx(2)->setAtomMapNum(2);
    mol->getAtomWithIdx(3)->setProp(common_properties::dummyLabel,
                                    std::string("foo"));

    std::string pkl;
    MolPickler::pickleMol(*mol, pkl,
                          PicklerOps::PropertyPickleOptions::AtomProps |
                              PicklerOps::PropertyPickleOptions::PrivateProps);

    CHECK(pkl.size() > basepkl.size());

    RWMol mol2(pkl);
    CHECK(mol2.getAtomWithIdx(0)->getProp<int>(common_properties::molParity) ==
          2);
    CHECK(mol2.getAtomWithIdx(0)->getProp<int>(
              common_properties::_ChiralityPossible) == 1);
    CHECK(mol2.getAtomWithIdx(0)->getProp<std::string>(
              common_properties::_CIPCode) == "S");

    CHECK(mol2.getAtomWithIdx(1)->getProp<int>(
              common_properties::molAtomMapNumber) == 1);
    CHECK(mol2.getAtomWithIdx(2)->getProp<int>(
              common_properties::molAtomMapNumber) == 2);
    CHECK(mol2.getAtomWithIdx(3)->getProp<std::string>(
              common_properties::dummyLabel) == "foo");
  }

  SECTION("attachment points") {
    auto mol = R"CTAB(
  Mrv2211 09062306242D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.7917 4.0833 0 0
M  V30 2 C -6.458 4.8533 0 0
M  V30 3 C -5.1243 4.0833 0 0
M  V30 4 * -6.458 4.34 0 0
M  V30 5 C -5.303 6.3405 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 4 5 ENDPTS=(3 1 2 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(2)->getProp<std::string>(
              common_properties::_MolFileBondAttach) == "ANY");
    CHECK(mol->getBondWithIdx(2)->getProp<std::string>(
              common_properties::_MolFileBondEndPts) == "(3 1 2 3)");

    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);

    // make sure the property names aren't in the pickle
    CHECK(pkl.find(common_properties::_MolFileBondAttach) == std::string::npos);
    CHECK(pkl.find(common_properties::_MolFileBondEndPts) == std::string::npos);
    // std::cerr << "!!!! " << pkl.size() << " " << basepkl.size() << std::endl;

    RWMol mol2(pkl);
    CHECK(mol2.getBondWithIdx(2)->getProp<std::string>(
              common_properties::_MolFileBondAttach) == "ANY");
    CHECK(mol2.getBondWithIdx(2)->getProp<std::string>(
              common_properties::_MolFileBondEndPts) == "(3 1 2 3)");
  }
}

TEST_CASE("parsing old pickles with many features") {
  std::string pklName = getenv("RDBASE");
  pklName += "/Code/GraphMol/test_data/mol_with_sgroups_and_stereo.pkl";

  auto m =
      "C/C=C/C[C@H](O)[C@@H](C)F |a:6,o2:4,r,SgD:5:data_pt:4.5::::|"_smiles;
  REQUIRE(m);
  std::ifstream inStream(pklName.c_str(), std::ios_base::binary);
  RWMol m2;
  // if the mol can be read, the primary problem was addressed
  MolPickler::molFromPickle(inStream, m2);
  CHECK(m2.getNumAtoms() == m->getNumAtoms());
  CHECK(MolToCXSmiles(*m) == MolToCXSmiles(m2));
}

TEST_CASE("github #7675 : pickling HasProp queries") {
  SECTION("basics") {
    auto mol = "CC"_smarts;
    REQUIRE(mol);
    mol->getAtomWithIdx(0)->expandQuery(makeHasPropQuery<Atom>("foo"));
    mol->getBondWithIdx(0)->expandQuery(makeHasPropQuery<Bond>("foo"));
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    RWMol mol2(pkl);
    REQUIRE(mol2.getAtomWithIdx(0)->hasQuery());
  }
}

TEST_CASE("pickling HasPropWithValue queries") {
  SECTION("basics") {
    if (0) {
      auto mol = "CC"_smarts;
      REQUIRE(mol);
      mol->getAtomWithIdx(0)->expandQuery(
          makePropQuery<Atom, int>("foo", 1, 2));
      mol->getBondWithIdx(0)->expandQuery(
          makePropQuery<Bond, int>("bar", 1, 0));
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      RWMol mol2(pkl);
      REQUIRE(mol2.getAtomWithIdx(0)->hasQuery());
      REQUIRE(mol2.getAtomWithIdx(1)->hasQuery());
    }
    {
      auto mol = "CC"_smarts;
      REQUIRE(mol);
      mol->getAtomWithIdx(0)->expandQuery(
          makePropQuery<Atom, std::string>("foo", "asdfs"));
      mol->getBondWithIdx(0)->expandQuery(
          makePropQuery<Bond, std::string>("bar", "dsafasdf"));
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      RWMol mol2(pkl);
      REQUIRE(mol2.getAtomWithIdx(0)->hasQuery());
      REQUIRE(mol2.getAtomWithIdx(1)->hasQuery());
    }
    {
      auto mol = "CC"_smarts;
      REQUIRE(mol);
      ExplicitBitVect bv(10);
      mol->getAtomWithIdx(0)->expandQuery(
          makePropQuery<Atom, ExplicitBitVect>("foo", bv, 0.1));
      mol->getBondWithIdx(0)->expandQuery(
          makePropQuery<Bond, ExplicitBitVect>("bar", bv, 0.1));
      std::string pkl;
      MolPickler::pickleMol(*mol, pkl);
      RWMol mol2(pkl);
      REQUIRE(mol2.getAtomWithIdx(0)->hasQuery());
      REQUIRE(mol2.getAtomWithIdx(1)->hasQuery());
    }
  }
}

TEST_CASE(
    "Errors in pickling and depickling mols with < 256 atoms and more than 255 bonds") {
  auto mol = R"CTAB(Generic RDKit S-group binary-round-trip reproducer
  RDKit 3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 129 257 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.500000 0.000000 0.000000 0
M  V30 3 C 3.000000 0.000000 0.000000 0
M  V30 4 C 4.500000 0.000000 0.000000 0
M  V30 5 C 6.000000 0.000000 0.000000 0
M  V30 6 C 7.500000 0.000000 0.000000 0
M  V30 7 C 9.000000 0.000000 0.000000 0
M  V30 8 C 10.500000 0.000000 0.000000 0
M  V30 9 C 12.000000 0.000000 0.000000 0
M  V30 10 C 13.500000 0.000000 0.000000 0
M  V30 11 C 15.000000 0.000000 0.000000 0
M  V30 12 C 16.500000 0.000000 0.000000 0
M  V30 13 C 18.000000 0.000000 0.000000 0
M  V30 14 C 0.000000 1.500000 0.000000 0
M  V30 15 C 1.500000 1.500000 0.000000 0
M  V30 16 C 3.000000 1.500000 0.000000 0
M  V30 17 C 4.500000 1.500000 0.000000 0
M  V30 18 C 6.000000 1.500000 0.000000 0
M  V30 19 C 7.500000 1.500000 0.000000 0
M  V30 20 C 9.000000 1.500000 0.000000 0
M  V30 21 C 10.500000 1.500000 0.000000 0
M  V30 22 C 12.000000 1.500000 0.000000 0
M  V30 23 C 13.500000 1.500000 0.000000 0
M  V30 24 C 15.000000 1.500000 0.000000 0
M  V30 25 C 16.500000 1.500000 0.000000 0
M  V30 26 C 18.000000 1.500000 0.000000 0
M  V30 27 C 0.000000 3.000000 0.000000 0
M  V30 28 C 1.500000 3.000000 0.000000 0
M  V30 29 C 3.000000 3.000000 0.000000 0
M  V30 30 C 4.500000 3.000000 0.000000 0
M  V30 31 C 6.000000 3.000000 0.000000 0
M  V30 32 C 7.500000 3.000000 0.000000 0
M  V30 33 C 9.000000 3.000000 0.000000 0
M  V30 34 C 10.500000 3.000000 0.000000 0
M  V30 35 C 12.000000 3.000000 0.000000 0
M  V30 36 C 13.500000 3.000000 0.000000 0
M  V30 37 C 15.000000 3.000000 0.000000 0
M  V30 38 C 16.500000 3.000000 0.000000 0
M  V30 39 C 18.000000 3.000000 0.000000 0
M  V30 40 C 0.000000 4.500000 0.000000 0
M  V30 41 C 1.500000 4.500000 0.000000 0
M  V30 42 C 3.000000 4.500000 0.000000 0
M  V30 43 C 4.500000 4.500000 0.000000 0
M  V30 44 C 6.000000 4.500000 0.000000 0
M  V30 45 C 7.500000 4.500000 0.000000 0
M  V30 46 C 9.000000 4.500000 0.000000 0
M  V30 47 C 10.500000 4.500000 0.000000 0
M  V30 48 C 12.000000 4.500000 0.000000 0
M  V30 49 C 13.500000 4.500000 0.000000 0
M  V30 50 C 15.000000 4.500000 0.000000 0
M  V30 51 C 16.500000 4.500000 0.000000 0
M  V30 52 C 18.000000 4.500000 0.000000 0
M  V30 53 C 0.000000 6.000000 0.000000 0
M  V30 54 C 1.500000 6.000000 0.000000 0
M  V30 55 C 3.000000 6.000000 0.000000 0
M  V30 56 C 4.500000 6.000000 0.000000 0
M  V30 57 C 6.000000 6.000000 0.000000 0
M  V30 58 C 7.500000 6.000000 0.000000 0
M  V30 59 C 9.000000 6.000000 0.000000 0
M  V30 60 C 10.500000 6.000000 0.000000 0
M  V30 61 C 12.000000 6.000000 0.000000 0
M  V30 62 C 13.500000 6.000000 0.000000 0
M  V30 63 C 15.000000 6.000000 0.000000 0
M  V30 64 C 16.500000 6.000000 0.000000 0
M  V30 65 C 18.000000 6.000000 0.000000 0
M  V30 66 C 0.000000 7.500000 0.000000 0
M  V30 67 C 1.500000 7.500000 0.000000 0
M  V30 68 C 3.000000 7.500000 0.000000 0
M  V30 69 C 4.500000 7.500000 0.000000 0
M  V30 70 C 6.000000 7.500000 0.000000 0
M  V30 71 C 7.500000 7.500000 0.000000 0
M  V30 72 C 9.000000 7.500000 0.000000 0
M  V30 73 C 10.500000 7.500000 0.000000 0
M  V30 74 C 12.000000 7.500000 0.000000 0
M  V30 75 C 13.500000 7.500000 0.000000 0
M  V30 76 C 15.000000 7.500000 0.000000 0
M  V30 77 C 16.500000 7.500000 0.000000 0
M  V30 78 C 18.000000 7.500000 0.000000 0
M  V30 79 C 0.000000 9.000000 0.000000 0
M  V30 80 C 1.500000 9.000000 0.000000 0
M  V30 81 C 3.000000 9.000000 0.000000 0
M  V30 82 C 4.500000 9.000000 0.000000 0
M  V30 83 C 6.000000 9.000000 0.000000 0
M  V30 84 C 7.500000 9.000000 0.000000 0
M  V30 85 C 9.000000 9.000000 0.000000 0
M  V30 86 C 10.500000 9.000000 0.000000 0
M  V30 87 C 12.000000 9.000000 0.000000 0
M  V30 88 C 13.500000 9.000000 0.000000 0
M  V30 89 C 15.000000 9.000000 0.000000 0
M  V30 90 C 16.500000 9.000000 0.000000 0
M  V30 91 C 18.000000 9.000000 0.000000 0
M  V30 92 C 0.000000 10.500000 0.000000 0
M  V30 93 C 1.500000 10.500000 0.000000 0
M  V30 94 C 3.000000 10.500000 0.000000 0
M  V30 95 C 4.500000 10.500000 0.000000 0
M  V30 96 C 6.000000 10.500000 0.000000 0
M  V30 97 C 7.500000 10.500000 0.000000 0
M  V30 98 C 9.000000 10.500000 0.000000 0
M  V30 99 C 10.500000 10.500000 0.000000 0
M  V30 100 C 12.000000 10.500000 0.000000 0
M  V30 101 C 13.500000 10.500000 0.000000 0
M  V30 102 C 15.000000 10.500000 0.000000 0
M  V30 103 C 16.500000 10.500000 0.000000 0
M  V30 104 C 18.000000 10.500000 0.000000 0
M  V30 105 C 0.000000 12.000000 0.000000 0
M  V30 106 C 1.500000 12.000000 0.000000 0
M  V30 107 C 3.000000 12.000000 0.000000 0
M  V30 108 C 4.500000 12.000000 0.000000 0
M  V30 109 C 6.000000 12.000000 0.000000 0
M  V30 110 C 7.500000 12.000000 0.000000 0
M  V30 111 C 9.000000 12.000000 0.000000 0
M  V30 112 C 10.500000 12.000000 0.000000 0
M  V30 113 C 12.000000 12.000000 0.000000 0
M  V30 114 C 13.500000 12.000000 0.000000 0
M  V30 115 C 15.000000 12.000000 0.000000 0
M  V30 116 C 16.500000 12.000000 0.000000 0
M  V30 117 C 18.000000 12.000000 0.000000 0
M  V30 118 C 0.000000 13.500000 0.000000 0
M  V30 119 C 1.500000 13.500000 0.000000 0
M  V30 120 C 3.000000 13.500000 0.000000 0
M  V30 121 C 4.500000 13.500000 0.000000 0
M  V30 122 C 6.000000 13.500000 0.000000 0
M  V30 123 C 7.500000 13.500000 0.000000 0
M  V30 124 C 9.000000 13.500000 0.000000 0
M  V30 125 C 10.500000 13.500000 0.000000 0
M  V30 126 C 12.000000 13.500000 0.000000 0
M  V30 127 C 13.500000 13.500000 0.000000 0
M  V30 128 C 15.000000 13.500000 0.000000 0
M  V30 129 C 16.500000 13.500000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 1 15 16
M  V30 16 1 16 17
M  V30 17 1 17 18
M  V30 18 1 18 19
M  V30 19 1 19 20
M  V30 20 1 20 21
M  V30 21 1 21 22
M  V30 22 1 22 23
M  V30 23 1 23 24
M  V30 24 1 24 25
M  V30 25 1 25 26
M  V30 26 1 26 27
M  V30 27 1 27 28
M  V30 28 1 28 29
M  V30 29 1 29 30
M  V30 30 1 30 31
M  V30 31 1 31 32
M  V30 32 1 32 33
M  V30 33 1 33 34
M  V30 34 1 34 35
M  V30 35 1 35 36
M  V30 36 1 36 37
M  V30 37 1 37 38
M  V30 38 1 38 39
M  V30 39 1 39 40
M  V30 40 1 40 41
M  V30 41 1 41 42
M  V30 42 1 42 43
M  V30 43 1 43 44
M  V30 44 1 44 45
M  V30 45 1 45 46
M  V30 46 1 46 47
M  V30 47 1 47 48
M  V30 48 1 48 49
M  V30 49 1 49 50
M  V30 50 1 50 51
M  V30 51 1 51 52
M  V30 52 1 52 53
M  V30 53 1 53 54
M  V30 54 1 54 55
M  V30 55 1 55 56
M  V30 56 1 56 57
M  V30 57 1 57 58
M  V30 58 1 58 59
M  V30 59 1 59 60
M  V30 60 1 60 61
M  V30 61 1 61 62
M  V30 62 1 62 63
M  V30 63 1 63 64
M  V30 64 1 64 65
M  V30 65 1 65 66
M  V30 66 1 66 67
M  V30 67 1 67 68
M  V30 68 1 68 69
M  V30 69 1 69 70
M  V30 70 1 70 71
M  V30 71 1 71 72
M  V30 72 1 72 73
M  V30 73 1 73 74
M  V30 74 1 74 75
M  V30 75 1 75 76
M  V30 76 1 76 77
M  V30 77 1 77 78
M  V30 78 1 78 79
M  V30 79 1 79 80
M  V30 80 1 80 81
M  V30 81 1 81 82
M  V30 82 1 82 83
M  V30 83 1 83 84
M  V30 84 1 84 85
M  V30 85 1 85 86
M  V30 86 1 86 87
M  V30 87 1 87 88
M  V30 88 1 88 89
M  V30 89 1 89 90
M  V30 90 1 90 91
M  V30 91 1 91 92
M  V30 92 1 92 93
M  V30 93 1 93 94
M  V30 94 1 94 95
M  V30 95 1 95 96
M  V30 96 1 96 97
M  V30 97 1 97 98
M  V30 98 1 98 99
M  V30 99 1 99 100
M  V30 100 1 100 101
M  V30 101 1 101 102
M  V30 102 1 102 103
M  V30 103 1 103 104
M  V30 104 1 104 105
M  V30 105 1 105 106
M  V30 106 1 106 107
M  V30 107 1 107 108
M  V30 108 1 108 109
M  V30 109 1 109 110
M  V30 110 1 110 111
M  V30 111 1 111 112
M  V30 112 1 112 113
M  V30 113 1 113 114
M  V30 114 1 114 115
M  V30 115 1 115 116
M  V30 116 1 116 117
M  V30 117 1 117 118
M  V30 118 1 118 119
M  V30 119 1 119 120
M  V30 120 1 120 121
M  V30 121 1 121 122
M  V30 122 1 122 123
M  V30 123 1 123 124
M  V30 124 1 124 125
M  V30 125 1 125 126
M  V30 126 1 126 127
M  V30 127 1 127 128
M  V30 128 1 128 129
M  V30 129 1 129 1
M  V30 130 1 1 3
M  V30 131 1 2 4
M  V30 132 1 3 5
M  V30 133 1 4 6
M  V30 134 1 5 7
M  V30 135 1 6 8
M  V30 136 1 7 9
M  V30 137 1 8 10
M  V30 138 1 9 11
M  V30 139 1 10 12
M  V30 140 1 11 13
M  V30 141 1 12 14
M  V30 142 1 13 15
M  V30 143 1 14 16
M  V30 144 1 15 17
M  V30 145 1 16 18
M  V30 146 1 17 19
M  V30 147 1 18 20
M  V30 148 1 19 21
M  V30 149 1 20 22
M  V30 150 1 21 23
M  V30 151 1 22 24
M  V30 152 1 23 25
M  V30 153 1 24 26
M  V30 154 1 25 27
M  V30 155 1 26 28
M  V30 156 1 27 29
M  V30 157 1 28 30
M  V30 158 1 29 31
M  V30 159 1 30 32
M  V30 160 1 31 33
M  V30 161 1 32 34
M  V30 162 1 33 35
M  V30 163 1 34 36
M  V30 164 1 35 37
M  V30 165 1 36 38
M  V30 166 1 37 39
M  V30 167 1 38 40
M  V30 168 1 39 41
M  V30 169 1 40 42
M  V30 170 1 41 43
M  V30 171 1 42 44
M  V30 172 1 43 45
M  V30 173 1 44 46
M  V30 174 1 45 47
M  V30 175 1 46 48
M  V30 176 1 47 49
M  V30 177 1 48 50
M  V30 178 1 49 51
M  V30 179 1 50 52
M  V30 180 1 51 53
M  V30 181 1 52 54
M  V30 182 1 53 55
M  V30 183 1 54 56
M  V30 184 1 55 57
M  V30 185 1 56 58
M  V30 186 1 57 59
M  V30 187 1 58 60
M  V30 188 1 59 61
M  V30 189 1 60 62
M  V30 190 1 61 63
M  V30 191 1 62 64
M  V30 192 1 63 65
M  V30 193 1 64 66
M  V30 194 1 65 67
M  V30 195 1 66 68
M  V30 196 1 67 69
M  V30 197 1 68 70
M  V30 198 1 69 71
M  V30 199 1 70 72
M  V30 200 1 71 73
M  V30 201 1 72 74
M  V30 202 1 73 75
M  V30 203 1 74 76
M  V30 204 1 75 77
M  V30 205 1 76 78
M  V30 206 1 77 79
M  V30 207 1 78 80
M  V30 208 1 79 81
M  V30 209 1 80 82
M  V30 210 1 81 83
M  V30 211 1 82 84
M  V30 212 1 83 85
M  V30 213 1 84 86
M  V30 214 1 85 87
M  V30 215 1 86 88
M  V30 216 1 87 89
M  V30 217 1 88 90
M  V30 218 1 89 91
M  V30 219 1 90 92
M  V30 220 1 91 93
M  V30 221 1 92 94
M  V30 222 1 93 95
M  V30 223 1 94 96
M  V30 224 1 95 97
M  V30 225 1 96 98
M  V30 226 1 97 99
M  V30 227 1 98 100
M  V30 228 1 99 101
M  V30 229 1 100 102
M  V30 230 1 101 103
M  V30 231 1 102 104
M  V30 232 1 103 105
M  V30 233 1 104 106
M  V30 234 1 105 107
M  V30 235 1 106 108
M  V30 236 1 107 109
M  V30 237 1 108 110
M  V30 238 1 109 111
M  V30 239 1 110 112
M  V30 240 1 111 113
M  V30 241 1 112 114
M  V30 242 1 113 115
M  V30 243 1 114 116
M  V30 244 1 115 117
M  V30 245 1 116 118
M  V30 246 1 117 119
M  V30 247 1 118 120
M  V30 248 1 119 121
M  V30 249 1 120 122
M  V30 250 1 121 123
M  V30 251 1 122 124
M  V30 252 1 123 125
M  V30 253 1 124 126
M  V30 254 1 125 127
M  V30 255 1 126 128
M  V30 256 1 127 129
M  V30 257 1 128 1
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(1 128) XBONDS=(1 257) LABEL=Generic
M  V30 END SGROUP
M  V30 END CTAB
M  END
$$$$
)CTAB"_ctab;
  REQUIRE(mol);
  std::string pkl;
  MolPickler::pickleMol(*mol, pkl);
  RWMol mol2(pkl);

  // This should not throw
  auto v3kMolBlock = MolToV3KMolBlock(mol2);

  CHECK(v3kMolBlock.find(
            "M  V30 1 SUP 1 ATOMS=(1 128) XBONDS=(1 257) LABEL=Generic") !=
        std::string::npos);
}

TEST_CASE("OOB read in depickling") {
  auto m = "C1CCC1"_smiles;
  REQUIRE(m);
  std::string pkl;
  MolPickler::pickleMol(*m, pkl);
  REQUIRE(pkl.size() == 95);
  // these are bytes in the ring definitions of the pickle:
  REQUIRE(static_cast<unsigned int>(pkl[80]) == 4);
  REQUIRE(static_cast<unsigned int>(pkl[83]) == 2);
  REQUIRE(static_cast<unsigned int>(pkl[84]) == 3);
  SECTION("ring size zero") {
    pkl[80] = 0;
    {
      RWMol m2;
      std::istringstream inStream(pkl);
      REQUIRE_THROWS_AS(MolPickler::molFromPickle(inStream, m2),
                        MolPicklerException);
    }
    {
      RWMol m2;
      std::istringstream inStream(pkl);
      REQUIRE_THROWS_WITH(MolPickler::molFromPickle(inStream, m2),
                          "invalid ring size");
    }
  }
  SECTION("invalid ring member order") {
    pkl[83] = 3;
    pkl[84] = 2;
    {
      RWMol m2;
      std::istringstream inStream(pkl);
      REQUIRE_THROWS_AS(MolPickler::molFromPickle(inStream, m2),
                        MolPicklerException);
    }
    {
      RWMol m2;
      std::istringstream inStream(pkl);
      REQUIRE_THROWS_WITH(MolPickler::molFromPickle(inStream, m2),
                          "ring-bond not found");
    }
  }
}

TEST_CASE("unpickle ring families fast path and RDL fallback") {
  auto compareFamilies = [](const RingInfo &a, const RingInfo &b) {
    CHECK(a.areRingFamiliesInitialized() == b.areRingFamiliesInitialized());
    CHECK(a.numRingFamilies() == b.numRingFamilies());
    CHECK(a.numRelevantCycles() == b.numRelevantCycles());
    CHECK(a.atomRelevantCycles() == b.atomRelevantCycles());

    auto sortedFamilies = [](const RingInfo::VECT_INT_VECT &fams) {
      RingInfo::VECT_INT_VECT res;
      res.reserve(fams.size());
      for (const auto &f : fams) {
        auto copy = f;
        std::ranges::sort(copy);
        res.push_back(std::move(copy));
      }
      std::ranges::sort(res);
      return res;
    };
    CHECK(sortedFamilies(a.atomRingFamilies()) ==
          sortedFamilies(b.atomRingFamilies()));
    CHECK(sortedFamilies(a.bondRingFamilies()) ==
          sortedFamilies(b.bondRingFamilies()));
  };

  // Molecules where stored rings uniquely determine families/relevant cycles
  // (fast path skips RDL, so dp_urfData stays null).
  const std::vector<std::string> verifiedFastPathSmiles = {
      "",
      "CCCCO",
      "c1ccccc1",
      "C1CC1.C1CCC1",
      "C1CCC2(CC1)CCCCC2",
      "c1ccc2ccccc2c1",
      "c1ccc2cc3ccccc3cc2c1",
      "c1ccc2c(c1)ccc1ccccc12",
      "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
      "C1CC2CCC1C2",     // norbornane (bicyclo[2.2.1]): 5-rings share 2 bonds
      "C1CC2CCC(C1)C2",  // bicyclo[3.2.1]octane: 5- and 6-ring share 2 bonds
  };

  for (const auto &smi : verifiedFastPathSmiles) {
    INFO("fast path smiles: " << smi);
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(mol);
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol unpickled(pkl);

    REQUIRE(unpickled.getRingInfo()->areRingFamiliesInitialized());
    CHECK(unpickled.getRingInfo()->dp_urfData == nullptr);
    compareFamilies(*mol->getRingInfo(), *unpickled.getRingInfo());

    // Copy preserves initialized state without RDL data.
    ROMol copied(unpickled);
    CHECK(copied.getRingInfo()->areRingFamiliesInitialized());
    CHECK(copied.getRingInfo()->dp_urfData == nullptr);
    compareFamilies(*mol->getRingInfo(), *copied.getRingInfo());

    // Subsequent findSSSR and symmetrizeSSSR still work cleanly.
    {
      RWMol mForSssr(unpickled);
      VECT_INT_VECT sssr;
      MolOps::findSSSR(mForSssr, sssr);
      CHECK(mForSssr.getRingInfo()->areRingFamiliesInitialized());
      compareFamilies(*mol->getRingInfo(), *mForSssr.getRingInfo());
    }
    {
      RWMol mForSymm(unpickled);
      VECT_INT_VECT symm;
      MolOps::symmetrizeSSSR(mForSymm, symm);
      CHECK(mForSymm.getRingInfo()->areRingFamiliesInitialized());
      compareFamilies(*mol->getRingInfo(), *mForSymm.getRingInfo());
    }
    {
      RWMol mForSymmNoRecalc(unpickled);
      RWMol mFresh(*mol);
      VECT_INT_VECT symmUnpk;
      VECT_INT_VECT symmFresh;
      MolOps::symmetrizeSSSR(mForSymmNoRecalc, symmUnpk,
                             MolOps::SymmetrizeSSSRAlgorithm::RDL, false);
      MolOps::symmetrizeSSSR(mFresh, symmFresh,
                             MolOps::SymmetrizeSSSRAlgorithm::RDL, false);
      auto sortedRings = [](VECT_INT_VECT rings) {
        for (auto &r : rings) {
          std::ranges::sort(r);
        }
        std::ranges::sort(rings);
        return rings;
      };
      CHECK(sortedRings(symmUnpk) == sortedRings(symmFresh));
      CHECK(mForSymmNoRecalc.getRingInfo()->areRingFamiliesInitialized());
      compareFamilies(*mFresh.getRingInfo(), *mForSymmNoRecalc.getRingInfo());
    }
  }

  // Molecules that must fall back to RDL (bridged/cage/cyclic ring-adjacency
  // systems where relevant cycles or URFs are not uniquely determined by a
  // tree-fused strict-inequality basis).
  const std::vector<std::string> rdlFallbackSmiles = {
      "C1CC2CCC1CC2",                   // bicyclo[2.2.2]octane
      "C1C2CC3CC1CC(C2)C3",             // adamantane
      "C12C3C4C1C5C2C3C45",             // cubane
      "c1cc2ccc3cccc4ccc(c1)c2c34",     // pyrene
      "C1CC2CCC1CC1CCC(CC1)CC1CCC(CC1)CC1CCC(CC1)C2",
  };

  for (const auto &smi : rdlFallbackSmiles) {
    INFO("RDL fallback smiles: " << smi);
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(mol);
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol unpickled(pkl);

    REQUIRE(unpickled.getRingInfo()->areRingFamiliesInitialized());
    CHECK(unpickled.getRingInfo()->dp_urfData != nullptr);
    compareFamilies(*mol->getRingInfo(), *unpickled.getRingInfo());
  }

  // Dative bond forming a ring (not included by default findRingFamilies):
  // unpickling must match default findRingFamilies(*mol).
  {
    auto mol = v2::SmilesParse::MolFromSmiles("C1CC->N1");
    REQUIRE(mol);
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol unpickled(pkl);
    REQUIRE(unpickled.getRingInfo()->areRingFamiliesInitialized());
    compareFamilies(*mol->getRingInfo(), *unpickled.getRingInfo());
  }

  // Pickle created from a molecule with FastFindRings must NOT initialize
  // ring families on unpickle (matches existing behavior).
  {
    auto mol = v2::SmilesParse::MolFromSmiles("c1ccccc1");
    REQUIRE(mol);
    MolOps::fastFindRings(*mol);
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol unpickled(pkl);
    CHECK(!unpickled.getRingInfo()->areRingFamiliesInitialized());
  }

  // Depickling into a non-empty molecule must reset existing ring families
  // instead of duplicating them.
  {
    auto m1 = v2::SmilesParse::MolFromSmiles("C1CC1");
    auto m2 = v2::SmilesParse::MolFromSmiles("C1COC1");
    REQUIRE(m1);
    REQUIRE(m2);
    ROMol combined;
    std::string pkl;
    MolPickler::pickleMol(*m1, pkl);
    MolPickler::molFromPickle(pkl, combined);
    MolPickler::pickleMol(*m2, pkl);
    MolPickler::molFromPickle(pkl, combined);

    const auto *ri = combined.getRingInfo();
    REQUIRE(ri->areRingFamiliesInitialized());
    CHECK(ri->numRingFamilies() == ri->numRings());
    CHECK(ri->numRelevantCycles() == ri->atomRelevantCycles().size());

    ROMol fresh(combined);
    fresh.getRingInfo()->resetRingFamilies();
    MolOps::findRingFamilies(fresh);
    compareFamilies(*fresh.getRingInfo(), *ri);
  }

  // Without RDL data, atomRelevantCycles() must derive from ring families
  // rather than drifting when stored rings are mutated via addRing().
  {
    auto mol = v2::SmilesParse::MolFromSmiles("C1CC1.CCC");
    REQUIRE(mol);
    std::string pkl;
    MolPickler::pickleMol(*mol, pkl);
    ROMol unpickled(pkl);
    auto *ri = unpickled.getRingInfo();
    REQUIRE(ri->areRingFamiliesInitialized());
    REQUIRE(ri->dp_urfData == nullptr);
    const auto expectedRC = ri->atomRelevantCycles();
    const RingInfo::INT_VECT extraRing = {3, 4, 5};
    ri->addRing(extraRing, extraRing);
    CHECK(ri->numRelevantCycles() == 1);
    CHECK(ri->atomRelevantCycles() == expectedRC);
  }
}