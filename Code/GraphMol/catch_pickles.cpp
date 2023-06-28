//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

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
}

TEST_CASE("Stereo Group ID preservation: pickle and unpickle",
          "[pickle][StereoGroup]") {
  std::string pkl;

  auto run_checks = [](ROMol &m) {
    const auto &groups = m.getStereoGroups();
    REQUIRE(groups.size() == 3);

    auto g = groups[0];
    CHECK(g.getGroupType() == RDKit::StereoGroupType::STEREO_AND);
    CHECK(g.getId() == 8);
    CHECK(g.getAtoms() ==
          std::vector<Atom *>{m.getAtomWithIdx(3), m.getAtomWithIdx(5)});

    g = groups[1];
    CHECK(g.getGroupType() == RDKit::StereoGroupType::STEREO_OR);
    CHECK(g.getId() == 1);
    CHECK(g.getAtoms() == std::vector<Atom *>{m.getAtomWithIdx(7)});

    g = groups[2];
    CHECK(g.getGroupType() == RDKit::StereoGroupType::STEREO_AND);
    CHECK(g.getId() == 7);
    CHECK(g.getAtoms() == std::vector<Atom *>{m.getAtomWithIdx(1)});
  };

  {
    auto mol = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 O 1.299038 2.250000 0.000000 0
M  V30 4 C 2.598076 -0.000000 0.000000 0
M  V30 5 C 3.897114 0.750000 0.000000 0
M  V30 6 C 2.598076 -1.500000 0.000000 0
M  V30 7 C 3.897114 -2.250000 0.000000 0
M  V30 8 C 1.299038 -2.250000 0.000000 0
M  V30 9 C 1.299038 -3.750000 0.000000 0
M  V30 10 O 0.000000 -1.500000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 4 5 CFG=3
M  V30 5 1 4 6
M  V30 6 1 6 7 CFG=1
M  V30 7 1 6 8
M  V30 8 1 8 9 CFG=1
M  V30 9 1 8 10
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC8 ATOMS=(2 4 6)
M  V30 MDLV30/STEREL1 ATOMS=(1 8)
M  V30 MDLV30/STERAC7 ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);

    // Parse
    const auto &groups = mol->getStereoGroups();
    REQUIRE(groups.size() == 3);

    INFO("Reference mol");
    run_checks(*mol);

    MolPickler::pickleMol(*mol, pkl);

    REQUIRE(!pkl.empty());
  }
  RWMol mol2(pkl);

  INFO("roundtripped mol");
  run_checks(mol2);
}
