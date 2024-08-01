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