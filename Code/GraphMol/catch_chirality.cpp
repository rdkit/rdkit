//
//  Copyright (C) 2020-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("bond StereoInfo", "[unittest]") {
  SECTION("basics") {
    {
      auto mol = "CC=C(C#C)C=C"_smiles;
      REQUIRE(mol);
      auto sinfo = Chirality::detail::getStereoInfo(mol->getBondWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Bond_Double);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 5);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Unspecified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::None);
    }
    {
      auto mol = "CC=NC=N"_smiles;
      REQUIRE(mol);
      auto sinfo = Chirality::detail::getStereoInfo(mol->getBondWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Bond_Double);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == Chirality::StereoInfo::NOATOM);
    }
  }
  SECTION("stereo") {
    {
      auto mol = "C/C=C(/C#C)C"_smiles;
      REQUIRE(mol);

      CHECK(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
      CHECK(mol->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      CHECK(mol->getBondWithIdx(1)->getStereoAtoms()[1] == 3);

      auto sinfo = Chirality::detail::getStereoInfo(mol->getBondWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Bond_Double);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 5);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Specified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::Bond_Trans);
    }
    {  // check an example where one of the stereo atoms isn't the first
       // neighbor
      auto mol = "C/C=C(/C)C#C"_smiles;
      REQUIRE(mol);

      CHECK(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
      CHECK(mol->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      CHECK(mol->getBondWithIdx(1)->getStereoAtoms()[1] == 4);

      auto sinfo = Chirality::detail::getStereoInfo(mol->getBondWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Bond_Double);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 4);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Specified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::Bond_Trans);
    }
    {
      auto mol = "C/C=C(\\C#C)C"_smiles;
      REQUIRE(mol);

      CHECK(mol->getBondWithIdx(1)->getStereoAtoms().size() == 2);
      CHECK(mol->getBondWithIdx(1)->getStereoAtoms()[0] == 0);
      CHECK(mol->getBondWithIdx(1)->getStereoAtoms()[1] == 3);

      auto sinfo = Chirality::detail::getStereoInfo(mol->getBondWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Bond_Double);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 5);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Specified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::Bond_Cis);
    }
    {  // any bonds
      auto mol = "CC=C(C#C)C"_smiles;
      REQUIRE(mol);

      mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);

      auto sinfo = Chirality::detail::getStereoInfo(mol->getBondWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Bond_Double);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == Chirality::StereoInfo::NOATOM);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 5);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Unknown);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::None);
    }
  }
}
TEST_CASE("isBondPotentialStereoBond", "[unittest]") {
  SECTION("basics") {
    {
      auto mol = "CC=C(C#C)C=C"_smiles;
      REQUIRE(mol);
      CHECK(
          Chirality::detail::isBondPotentialStereoBond(mol->getBondWithIdx(1)));
      CHECK(!Chirality::detail::isBondPotentialStereoBond(
          mol->getBondWithIdx(5)));
      CHECK(!Chirality::detail::isBondPotentialStereoBond(
          mol->getBondWithIdx(3)));
      CHECK(!Chirality::detail::isBondPotentialStereoBond(
          mol->getBondWithIdx(4)));
    }
    {
      auto mol = "CC=NC=N"_smiles;
      REQUIRE(mol);
      CHECK(
          Chirality::detail::isBondPotentialStereoBond(mol->getBondWithIdx(1)));
      CHECK(!Chirality::detail::isBondPotentialStereoBond(
          mol->getBondWithIdx(3)));
    }
    {
      SmilesParserParams ps;
      ps.removeHs = false;
      std::unique_ptr<ROMol> mol{SmilesToMol("[H]C=CC=C([H])[H]", ps)};
      REQUIRE(mol);
      CHECK(!Chirality::detail::isBondPotentialStereoBond(
          mol->getBondWithIdx(1)));
      CHECK(!Chirality::detail::isBondPotentialStereoBond(
          mol->getBondWithIdx(3)));
    }
  }
  SECTION("ring size") {
    {
      auto m = "C1=CCCCC1"_smiles;
      REQUIRE(m);
      CHECK(
          !Chirality::detail::isBondPotentialStereoBond(m->getBondWithIdx(0)));
    }
    {
      auto m = "C1=CCCCCC1"_smiles;
      REQUIRE(m);
      CHECK(
          !Chirality::detail::isBondPotentialStereoBond(m->getBondWithIdx(0)));
    }
    {
      auto m = "C12=C(CCCC2)CCCCCC1"_smiles;
      REQUIRE(m);
      CHECK(
          !Chirality::detail::isBondPotentialStereoBond(m->getBondWithIdx(0)));
    }
  }
}

TEST_CASE("atom StereoInfo", "[unittest]") {
  SECTION("basics") {
    {
      auto mol = "CC(F)(Cl)CNC(C)C"_smiles;
      REQUIRE(mol);
      auto sinfo = Chirality::detail::getStereoInfo(mol->getAtomWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == 2);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 4);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Unspecified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::None);

      sinfo = Chirality::detail::getStereoInfo(mol->getAtomWithIdx(6));
      CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(sinfo.centeredOn == 6);
      REQUIRE(sinfo.controllingAtoms.size() == 3);
      CHECK(sinfo.controllingAtoms[0] == 5);
      CHECK(sinfo.controllingAtoms[1] == 7);
      CHECK(sinfo.controllingAtoms[2] == 8);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Unspecified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::None);
    }

    {
      auto mol = "C[C@](F)(Cl)CNC(C)C"_smiles;
      REQUIRE(mol);
      auto sinfo = Chirality::detail::getStereoInfo(mol->getAtomWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 4);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == 2);
      CHECK(sinfo.controllingAtoms[2] == 3);
      CHECK(sinfo.controllingAtoms[3] == 4);
      CHECK(sinfo.specified == Chirality::StereoSpecified::Specified);
      CHECK(sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW);
    }

    {
      auto mol = "CN1CC1N(F)C"_smiles;
      REQUIRE(mol);
      auto sinfo = Chirality::detail::getStereoInfo(mol->getAtomWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 3);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == 2);
      CHECK(sinfo.controllingAtoms[2] == 3);
    }

    {
      auto mol = "O[As](F)C[As]C[As]"_smiles;
      REQUIRE(mol);

      auto sinfo = Chirality::detail::getStereoInfo(mol->getAtomWithIdx(1));
      CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(sinfo.centeredOn == 1);
      REQUIRE(sinfo.controllingAtoms.size() == 3);
      CHECK(sinfo.controllingAtoms[0] == 0);
      CHECK(sinfo.controllingAtoms[1] == 2);
      CHECK(sinfo.controllingAtoms[2] == 3);

      sinfo = Chirality::detail::getStereoInfo(mol->getAtomWithIdx(4));
      CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(sinfo.centeredOn == 4);
      REQUIRE(sinfo.controllingAtoms.size() == 2);
      CHECK(sinfo.controllingAtoms[0] == 3);
      CHECK(sinfo.controllingAtoms[1] == 5);
    }
  }
}
TEST_CASE("isAtomPotentialTetrahedralCenter", "[unittest]") {
  SECTION("basics") {
    {
      auto mol = "CC(F)(Cl)CNC(C)(C)C"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(0)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(4)));
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(6)));
    }
    {
      auto mol = "CN1CC1N(F)C"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(4)));
    }
    {
      auto mol = "O=S(F)CC[S+]([O-])CS=O"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(5)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(8)));
    }
    {
      auto mol = "O=[Se](F)CC[Se+]([O-])C[Se]=O"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(5)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(8)));
    }
    {
      auto mol = "OP(F)CPCP"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(4)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(6)));
    }
    {
      auto mol = "O[As](F)C[As]C[As]"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(4)));
      CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(6)));
    }
    {
      auto mol = "O[P]([O-])(=O)OC"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
          mol->getAtomWithIdx(1)));
    }
  }
}
TEST_CASE("isAtomPotentialStereoAtom", "[unittest]") {
  SECTION("basics") {
    {
      auto mol = "CC(F)(Cl)CNC(C)(C)C"_smiles;
      REQUIRE(mol);
      for (const auto atom : mol->atoms()) {
        CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(atom) ==
              Chirality::detail::isAtomPotentialStereoAtom(atom));
      }
    }
    {
      auto mol = "CN1CC1N(F)C"_smiles;
      REQUIRE(mol);
      for (const auto atom : mol->atoms()) {
        CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(atom) ==
              Chirality::detail::isAtomPotentialStereoAtom(atom));
      }
    }
    {
      auto mol = "O=S(F)CC[S+]([O-])CS=O"_smiles;
      REQUIRE(mol);
      for (const auto atom : mol->atoms()) {
        CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(atom) ==
              Chirality::detail::isAtomPotentialStereoAtom(atom));
      }
    }
  }
}

TEST_CASE("possible stereochemistry on atoms", "[chirality]") {
  SECTION("specified") {
    {
      auto mol = "CC(C)(O)[C@](Cl)(F)I"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[0].centeredOn == 4);
      std::vector<unsigned> catoms = {1, 5, 6, 7};
      CHECK(stereoInfo[0].controllingAtoms == catoms);
    }
    {
      auto mol = "C[C@@H](O)[C@H](C)[C@H](C)O"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 3);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[0].centeredOn == 1);

      CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[1].centeredOn == 3);

      CHECK(stereoInfo[2].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[2].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[2].centeredOn == 5);
    }

    {
      auto mol = "FC(F)(F)[C@@H](O)[C@H](C)[C@H](C(F)(F)F)O"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 3);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[0].centeredOn == 4);

      CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[1].centeredOn == 6);

      CHECK(stereoInfo[2].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[2].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[2].centeredOn == 8);
    }
  }
  SECTION("simple unspecified") {
    {
      auto mol = "CC(C)(O)C(Cl)(F)I"_smiles;
      REQUIRE(mol);

      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
      CHECK(stereoInfo[0].centeredOn == 4);
      std::vector<unsigned> catoms = {1, 5, 6, 7};
      CHECK(stereoInfo[0].controllingAtoms == catoms);
    }
  }
  SECTION("atoms with unknown set, real") {
    auto mol = "FC(O)C"_smiles;
    REQUIRE(mol);
    mol->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BondDir::UNKNOWN);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    CHECK(stereoInfo[0].centeredOn == 1);
  }
  SECTION("atoms with unknown set, not real") {
    auto mol = "CC(O)C"_smiles;
    REQUIRE(mol);
    mol->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BondDir::UNKNOWN);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Isotopes") {
    {
      auto mol = "O[C@H](F)[18OH]"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
      CHECK(stereoInfo[0].centeredOn == 1);
      std::vector<unsigned> catoms = {0, 2, 3};
      CHECK(stereoInfo[0].controllingAtoms == catoms);
    }
  }
}

TEST_CASE("possible stereochemistry on bonds", "[chirality]") {
  SECTION("simplest") {
    {
      auto mol = "CC=CC"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      std::vector<unsigned> catoms = {0, Chirality::StereoInfo::NOATOM, 3,
                                      Chirality::StereoInfo::NOATOM};
      CHECK(stereoInfo[0].controllingAtoms == catoms);
    }
    {
      auto mol = "CC=C(C)C"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      CHECK(stereoInfo.size() == 0);
    }
    {
      auto mol = "CC=C"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      CHECK(stereoInfo.size() == 0);
    }
    {
      auto mol = "CC(F)=C(Cl)C"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 2);
      std::vector<unsigned> catoms = {0, 2, 4, 5};
      CHECK(stereoInfo[0].controllingAtoms == catoms);
    }
    {
      auto mol = "CC=C(Cl)C"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      std::vector<unsigned> catoms = {0, Chirality::StereoInfo::NOATOM, 3, 4};
      CHECK(stereoInfo[0].controllingAtoms == catoms);
    }
  }
  SECTION("bond with unknown set, real") {
    auto mol = "CC=C(C)F"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
  }
  SECTION("bond with unknown set, not real") {
    auto mol = "CC=C(C)C"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    CHECK(stereoInfo.size() == 0);
  }
}

TEST_CASE("para-stereocenters and assignStereochemistry", "[chirality]") {
  SECTION("simplest") {
    auto mol = "CC(F)C(C)C(C)F"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 3);

    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].controllingAtoms.size() == 3);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 3);
    CHECK(stereoInfo[1].controllingAtoms.size() == 3);
    CHECK(stereoInfo[2].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[2].centeredOn == 5);
    CHECK(stereoInfo[2].controllingAtoms.size() == 3);
  }

  SECTION("including bonds") {
    // thanks to Salome Rieder for this nasty example
    auto mol = "CC=CC(C=CC)C(C)C(C=CC)C=CC"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    CHECK(stereoInfo.size() == 7);

    std::sort(stereoInfo.begin(), stereoInfo.end(),
              [](const Chirality::StereoInfo &a,
                 const Chirality::StereoInfo &b) -> bool {
                return (a.type < b.type) && (a.centeredOn < b.centeredOn) &&
                       (a.specified < b.specified) &&
                       (a.descriptor < b.descriptor) &&
                       (a.controllingAtoms < b.controllingAtoms);
              });
    REQUIRE(stereoInfo.size() == 7);

    CHECK(stereoInfo[6].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[6].centeredOn == 13);
    CHECK(stereoInfo[5].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[5].centeredOn == 10);
    CHECK(stereoInfo[4].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[4].centeredOn == 4);
    CHECK(stereoInfo[3].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[3].centeredOn == 1);

    CHECK(stereoInfo[2].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[2].centeredOn == 9);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 7);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 3);
  }

  SECTION("sugar fun") {
    auto mol = "C1(O)C(O)C(O)C(O)C(O)C1O"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 6);
    for (const auto &si : stereoInfo) {
      CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(si.centeredOn % 2 == 0);
      CHECK(si.specified == Chirality::StereoSpecified::Unspecified);
    }
  }
}

TEST_CASE("ring stereochemistry", "[chirality]") {
  SECTION("specified") {
    auto mol = "C[C@H]1CC[C@@H](C)CC1"_smiles;
    REQUIRE(mol);
    // std::cerr << "------------ 1 -------------" << std::endl;

    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 4);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
  }
  SECTION("unspecified") {
    auto mol = "CC1CCC(C)CC1"_smiles;
    REQUIRE(mol);
    // std::cerr << "------------ 2 -------------" << std::endl;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);

    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 4);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
  }
  SECTION("four ring") {
    auto mol = "C[C@H]1C[C@@H](C)C1"_smiles;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);

    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 3);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
  }
  SECTION("four ring unspecified") {
    auto mol = "CC1CC(C)C1"_smiles;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);

    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 3);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
  }
}
#if 0
// FIX: the double bond stereo in rings isn't working. This also fails with the canonicalizer, so it's not unique to this code
TEST_CASE("tricky recursive example from Dan Nealschneider", "[chirality]") {
  SECTION("adapted") {
    auto mol = "CC=C1CCC(O)CC1"_smiles;
    REQUIRE(mol);
    mol->updatePropertyCache();
    MolOps::setBondStereoFromDirections(*mol);
    std::cerr << "------------ 1 -------------" << std::endl;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 5);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[1].centeredOn == 1);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
  }
  SECTION("simplified") {
    // can't sanitize this because the current (2020.03) assignStereochemistry
    // code doesn't recognize the stereo here and removes it
    SmilesParserParams ps;
    ps.sanitize = false;
    ps.removeHs = false;
    std::unique_ptr<ROMol> mol(SmilesToMol("C/C=C1/C[C@H](O)C1", ps));
    REQUIRE(mol);
    mol->updatePropertyCache();
    MolOps::setBondStereoFromDirections(*mol);
    std::cerr << "------------ 2 -------------" << std::endl;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 4);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[1].centeredOn == 1);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
  }
  // FIX this still isn't working
  SECTION("unspecified") {
    auto mol = "CC=C1C[CH](O)C1"_smiles;
    REQUIRE(mol);
    std::cerr << "------------ 3 -------------" << std::endl;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 4);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[1].centeredOn == 1);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
  }
}
#endif

TEST_CASE("unknown stereo", "[chirality]") {
  SECTION("atoms") {
    auto mol = "CC(O)C[C@@H](O)F"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    mol->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BondDir::UNKNOWN);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 4);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
  }
  SECTION("atoms2") {
    // artificial situation: "squiggly bond" overrides the specified atomic
    // stereo
    auto mol = "C[C@H](O)C[C@@H](O)F"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    mol->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BondDir::UNKNOWN);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 4);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
  }
  SECTION("bonds") {
    {
      auto mol = "CC=CC"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getBondBetweenAtoms(1, 2));
      mol->getBondBetweenAtoms(1, 2)->setBondDir(Bond::BondDir::EITHERDOUBLE);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    }
    {
      auto mol = "CC=CC=C"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getBondBetweenAtoms(1, 2));
      mol->getBondBetweenAtoms(1, 2)->setBondDir(Bond::BondDir::EITHERDOUBLE);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    }
  }
  SECTION("bonds with squiggle bonds") {
    {  // to begin atom
      auto mol = "CC=CC"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getBondBetweenAtoms(0, 1));
      mol->getBondBetweenAtoms(0, 1)->setBondDir(Bond::BondDir::UNKNOWN);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    }
    {  // to end atom
      auto mol = "CC=CC"_smiles;
      REQUIRE(mol);
      REQUIRE(mol->getBondBetweenAtoms(2, 3));
      mol->getBondBetweenAtoms(2, 3)->setBondDir(Bond::BondDir::UNKNOWN);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unknown);
    }
  }
}

TEST_CASE("cleaning chirality", "[chirality]") {
  SECTION("atoms") {
    auto mol = "CC(O)C"_smiles;
    REQUIRE(mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    {
      // by default we don't clean up, so the chiral center survives even though
      // we don't get any results:
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      CHECK(stereoInfo.size() == 0);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    }
    {
      bool cleanIt = true;
      auto stereoInfo = Chirality::findPotentialStereo(*mol, cleanIt);
      CHECK(stereoInfo.size() == 0);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("bonds") {
    auto mol = "CC=C(C)C"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOCIS);
    {
      // by default we don't clean up, so the stereo bond survives even though
      // we don't get any results:
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      CHECK(stereoInfo.size() == 0);
      CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOCIS);
    }
    {
      bool cleanIt = true;
      auto stereoInfo = Chirality::findPotentialStereo(*mol, cleanIt);
      CHECK(stereoInfo.size() == 0);
      CHECK(mol->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREONONE);
    }
  }
}

TEST_CASE("flagPossible", "[chirality]") {
  SECTION("atoms") {
    auto mol = "CC(O)[C@H](F)O"_smiles;
    REQUIRE(mol);
    {
      // by default we do use flagPossible:
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
      CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[1].centeredOn == 3);
      CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
    }
    {
      bool cleanIt = false;
      bool flagPossible = false;
      auto stereoInfo =
          Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].centeredOn == 3);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
    }
  }
  SECTION("bonds") {
    auto mol = "CC=C/C=C/C"_smiles;
    REQUIRE(mol);
    {
      // by default we do use flagPossible
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
      CHECK(stereoInfo[1].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[1].centeredOn == 3);
      CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Specified);
    }
    {
      bool cleanIt = true;
      bool flagPossible = false;
      auto stereoInfo =
          Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
      REQUIRE(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Bond_Double);
      CHECK(stereoInfo[0].centeredOn == 3);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Specified);
    }
  }
}

TEST_CASE("cleanup after removing possible centers", "[chirality]") {
  SECTION("atoms1") {
    auto mol = "FC(Cl)(F)C(C(Cl)(F)F)I"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    CHECK(stereoInfo.empty());
  }
  SECTION("bonds1") {
    auto mol = "FC(Cl)(F)C(C(Cl)(F)F)=CF"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    CHECK(stereoInfo.empty());
  }
  SECTION("atoms2") {
    auto mol = "ClC(F)(F)C(=CC(F)C=C(C(F)(F)Cl)C(F)(F)Cl)C(Cl)(F)F"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    CHECK(stereoInfo.empty());
  }
}

TEST_CASE("findPotentialStereo problems related to #3490", "[chirality][bug]") {
  SECTION("example 1") {
    auto mol = "CC1CC(O)C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 3);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
  }
  SECTION("example 2a") {
    auto mol = "C(C(C)C1)C12CCN2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("example 2b") {
    auto mol = "CC(C1)CC12CCN2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("example 2c") {
    auto mol = "C([C@H](C)C1)[C@]12CCN2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("example 2d") {
    auto mol = "C[C@H](C1)C[C@]12CCN2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("example 3") {
    auto mol = "C(C(C)C1)C12CN(C3)CCCCC23"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 4);  // [1, 4, 6, 12]
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 4);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(stereoInfo[2].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[2].centeredOn == 6);
    CHECK(stereoInfo[2].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(stereoInfo[3].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[3].centeredOn == 12);
    CHECK(stereoInfo[3].specified == Chirality::StereoSpecified::Unspecified);
  }
}
TEST_CASE("ring stereo finding is overly aggressive", "[chirality][bug]") {
  SECTION("Finding too much 1a") {
    auto mol = "CC1CCCCC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 1b") {
    auto mol = "CC1CCC(C)CC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("Finding too much 1c") {
    auto mol = "C[C@H]1CCC(C)CC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("Finding too much 1d") {
    auto mol = "CC1(C)CCCCC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 1e") {
    auto mol = "CC1(C)CCC(C)CC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 1f") {
    auto mol = "C2CC2C1(C2CC2)CCCCC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 1g") {
    auto mol = "CC1CC2(CCC2)C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 1h") {
    auto mol = "CC1CC2(CC(C)C2)C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 3);
  }

  SECTION("Finding too much 2a") {
    auto mol = "CC1CCNCC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 2b") {
    auto mol = "CC1CCN(C)CC1"_smiles;  // 3-coordinate N is not stereogenic
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }
  SECTION("Finding too much 3a") {
    auto mol = "CC1CCC1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 0);
  }

  SECTION("Finding too much 3b") {
    auto mol = "CC1CC(C)C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
  SECTION("fused rings 1") {
    auto mol = "C1CCC2CCCCC2C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }

  SECTION("fused rings 2") {
    auto mol = "C1CC2CCCC2C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }

  SECTION("cages 1") {
    auto mol = "CC1CN2CCC1CC2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 3);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[1].centeredOn == 3);
    CHECK(stereoInfo[2].centeredOn == 6);
  }
  SECTION("cages 1b") {
    auto mol = "O1CN2CCC1CC2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].centeredOn == 2);
    CHECK(stereoInfo[1].centeredOn == 5);
  }
  SECTION("cages 2") {
    auto mol = "C1CC2(O)CCC1(C)CC2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].centeredOn == 2);
    CHECK(stereoInfo[1].centeredOn == 6);
  }
  SECTION("cages 3") {
    auto mol = "C1CC2(O)CCC1CC2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].centeredOn == 2);
    CHECK(stereoInfo[1].centeredOn == 6);
  }
  SECTION("adamantyl") {
    auto mol = "CC12CC3CC(CC(C3)C1)C2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 4);
  }
  SECTION("bug 1a") {
    // example that came up during testing
    auto mol = "C(=O)C(C(C)N2C=C2)C(=O)"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].centeredOn == 3);
  }
  SECTION("bug 1b") {
    // example that came up during testing
    auto mol = "C(=O)C(C(CC)c2ccc(Cl)cc2)C(=O)"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].centeredOn == 3);
  }

  SECTION("bug 1c") {
    // example that came up during testing
    auto mol = "O=CC(C=O)C(C)n2cccc2"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].centeredOn == 5);
  }

  SECTION("bug 1c") {
    // example that came up during testing
    auto mol = "C(=O)C(C(C)n2cccc2)C(=O)"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].centeredOn == 3);
  }

  SECTION("bug 1d") {
    // example that came up during testing
    auto mol = "C(O)C(C(C)n2cccc2)C(O)"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    REQUIRE(stereoInfo.size() == 1);
    CHECK(stereoInfo[0].centeredOn == 3);
  }
  SECTION("just a bug") {
    // example that came up during testing

    auto mol = "CC1=CN(C2OC(CNC(=O)C3c4ccccc4Sc4ccccc43)CC2)C(=O)NC1=O"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    auto stereoInfo =
        Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
    CHECK(stereoInfo.size() == 2);
  }
}

TEST_CASE(
    "github #3631: Ring stereochemistry not properly removed from N atoms",
    "[chirality][bug]") {
  SECTION("basics") {
    SmilesParserParams ps;
    ps.sanitize = false;
    ps.removeHs = false;
    std::unique_ptr<RWMol> mol{SmilesToMol("C[N@]1C[C@@](F)(Cl)C1", ps)};
    REQUIRE(mol);
    MolOps::sanitizeMol(*mol);

    CHECK(mol->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(3)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    bool cleanIt = true;
    bool flagPossible = true;
    bool force = true;
    {
      RWMol mol2(*mol);
      auto stereoInfo =
          Chirality::findPotentialStereo(mol2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 0);
    }
    {
      RWMol mol2(*mol);
      MolOps::assignStereochemistry(mol2, cleanIt, force, flagPossible);
      CHECK(mol2.getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(mol2.getAtomWithIdx(3)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("default behavior") {
    auto mol = "C[N@]1C[C@@](F)(Cl)C1"_smiles;
    REQUIRE(mol);
    auto smiles = MolToSmiles(*mol);
    CHECK(smiles == "CN1CC(F)(Cl)C1");
    bool cleanIt = true;
    bool flagPossible = true;
    bool force = true;
    CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(3)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    {
      RWMol mol2(*mol);
      auto stereoInfo =
          Chirality::findPotentialStereo(mol2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 0);
    }
    {
      RWMol mol2(*mol);
      MolOps::assignStereochemistry(mol2, cleanIt, force, flagPossible);
      CHECK(mol2.getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(mol2.getAtomWithIdx(3)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("don't overcorrect") {
    auto mol = "C[N@]1O[C@@](F)(Cl)C1"_smiles;
    REQUIRE(mol);
    bool cleanIt = true;
    bool flagPossible = true;
    bool force = true;
    {
      RWMol mol2(*mol);
      auto stereoInfo =
          Chirality::findPotentialStereo(mol2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 1);
      CHECK(stereoInfo[0].centeredOn == 3);
    }
    {
      RWMol mol2(*mol);
      MolOps::assignStereochemistry(mol2, cleanIt, force, flagPossible);
      CHECK(mol2.getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(mol2.getAtomWithIdx(3)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
}

TEST_CASE("N Chirality in rings") {
  SECTION("basics 4 coordinate") {
    {
      auto mol = "CC1CC2CC[N@@+]1(C)OC2"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(6)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(6)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {
      auto mol = "C[N@@+](F)(Cl)O"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(1)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("basics 3 coordinate") {
    {
      auto mol = "CC1CC2CC[N@@]1OC2"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(6)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(6)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {
      auto mol = "C1CC[N@]2OCCCC2C1"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(3)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(3)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("ring stereo") {
    {  // real chirality
      auto mol = "C[C@H]1CC[N@@+](C)(O)OC1"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(4)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(4)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(1)->getAtomicNum() == 6);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {  // ring stereo
      auto mol = "C[C@H]1CC[N@@+](C)(O)CC1"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(4)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(4)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(mol->getAtomWithIdx(1)->getAtomicNum() == 6);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {  // three-ring degree-three ring stereo
      auto mol = "C[C@H]1[C@@H](C)[N@]1C"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(4)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(4)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {  // CHEMBL79374
      auto mol = "Cn1ncc([C@]23CC[N@](CC2)C3)n1"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(8)->getAtomicNum() == 7);
      CHECK(mol->getAtomWithIdx(8)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {  // derived from CHEMBL79374
      auto mol = "Cn1ncc([C@]23CC[C@](CC2)C3)n1"_smiles;
      REQUIRE(mol);
      CHECK(mol->getAtomWithIdx(8)->getAtomicNum() == 6);
      CHECK(mol->getAtomWithIdx(8)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
}

TEST_CASE(
    "Github #4115: RemoveStereochemistry should also remove stereogroups") {
  SECTION("basics") {
    auto mol = "C[C@H](O)[C@@H](C)F |o1:1,3,r|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(3)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getStereoGroups().size() == 1);
    MolOps::removeStereochemistry(*mol);
    CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getAtomWithIdx(3)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol->getStereoGroups().empty());
  }
}

TEST_CASE(
    "Github #4155: Problem finding stereocenters in bridged bicyclics with "
    "4-rings") {
  SECTION("specified") {
    std::vector<std::string> smis = {
        "C[C@H]1CC[C@H](CC1)C(N)=O", "C[C@]12CC[C@](CC1)(C2)C(N)=O",
        "C[C@H]1C[C@H](C1)C(N)=O", "C[C@]12C[C@](C1)(CC2)C(N)=O"};
    for (const auto &smi : smis) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bool cleanIt = true;
      bool flagPossible = true;
      auto stereoInfo =
          Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].centeredOn == 1);
    }
  }
  SECTION("unspecified") {
    std::vector<std::string> smis = {
        "CC1CCC(CC1)C(N)=O", "CC12CCC(CC1)(C2)C(N)=O", "CC1CC(C1)C(N)=O",
        "CC12CC(C1)(CC2)C(N)=O"};
    for (const auto &smi : smis) {
      std::unique_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      bool cleanIt = true;
      bool flagPossible = true;
      auto stereoInfo =
          Chirality::findPotentialStereo(*mol, cleanIt, flagPossible);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].centeredOn == 1);
    }
  }
}

TEST_CASE("pickBondsToWedge() should avoid double bonds") {
  SECTION("simplest") {
    auto mol = "OC=C[C@H](C1CC1)C2CCC2"_smiles;
    REQUIRE(mol);
    auto wedgedBonds = pickBondsToWedge(*mol);
    REQUIRE(wedgedBonds.size() == 1);
    auto head = wedgedBonds.begin();
    CHECK(head->first == 3);
    CHECK(head->second == 3);
  }
  SECTION("simplest, specified double bond") {
    auto mol = "OC=C[C@H](C1CC1)C2CCC2"_smiles;
    REQUIRE(mol);
    mol->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    mol->getBondBetweenAtoms(1, 2)->setStereo(Bond::BondStereo::STEREOCIS);
    auto wedgedBonds = pickBondsToWedge(*mol);
    REQUIRE(wedgedBonds.size() == 1);
    auto head = wedgedBonds.begin();
    CHECK(head->first == 3);
    CHECK(head->second == 3);
  }
  SECTION("prefer unspecified bond stereo") {
    auto mol = "OC=C[C@H](C=CF)(C=CC)"_smiles;
    REQUIRE(mol);
    mol->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    mol->getBondBetweenAtoms(1, 2)->setStereo(Bond::BondStereo::STEREOCIS);
    mol->getBondBetweenAtoms(4, 5)->setStereoAtoms(3, 6);
    mol->getBondBetweenAtoms(4, 5)->setStereo(Bond::BondStereo::STEREOANY);
    auto wedgedBonds = pickBondsToWedge(*mol);
    REQUIRE(wedgedBonds.size() == 1);
    auto head = wedgedBonds.begin();
    CHECK(head->first == 6);
    CHECK(head->second == 3);
  }
}

TEST_CASE("addWavyBondsForStereoAny()") {
  SECTION("simplest") {
    auto mol = "CC=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    addWavyBondsForStereoAny(*mol);
    CHECK(mol->getBondWithIdx(0)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("don't reset flags") {
    auto mol = "CC=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    bool clearFlags = false;
    addWavyBondsForStereoAny(*mol, clearFlags);
    CHECK(mol->getBondWithIdx(0)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOANY);
  }
  SECTION("avoid double bonds") {
    auto mol = "CC=CC(CC)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(5)->setStereoAtoms(2, 7);
    mol->getBondWithIdx(5)->setStereo(Bond::BondStereo::STEREOANY);
    addWavyBondsForStereoAny(*mol);
    CHECK(mol->getBondWithIdx(6)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(5)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("avoid chiral atoms") {
    auto mol = "C[C@](F)(Cl)C(C)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(5)->setStereoAtoms(1, 7);
    mol->getBondWithIdx(5)->setStereo(Bond::BondStereo::STEREOANY);
    addWavyBondsForStereoAny(*mol);
    CHECK(mol->getBondWithIdx(4)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(5)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("prefer atoms with less neighbors") {
    auto mol = "CC(F)(Cl)C(CF)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(6)->setStereoAtoms(1, 8);
    mol->getBondWithIdx(6)->setStereo(Bond::BondStereo::STEREOANY);
    addWavyBondsForStereoAny(*mol);
    CHECK(mol->getBondWithIdx(7)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(6)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("more complex") {
    auto mol = "CC=CC(C=CO)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(6)->setStereoAtoms(2, 8);
    mol->getBondWithIdx(6)->setStereo(Bond::BondStereo::STEREOANY);
    addWavyBondsForStereoAny(*mol);
    CHECK(mol->getBondWithIdx(7)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(6)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("no solution without changing threshold") {
    auto mol = "CC=CC=CC=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOCIS);
    mol->getBondWithIdx(3)->setStereoAtoms(2, 5);
    mol->getBondWithIdx(3)->setStereo(Bond::BondStereo::STEREOANY);
    mol->getBondWithIdx(5)->setStereoAtoms(4, 7);
    mol->getBondWithIdx(5)->setStereo(Bond::BondStereo::STEREOCIS);
    addWavyBondsForStereoAny(*mol);
    // we didn't actually do anything:
    CHECK(mol->getBondWithIdx(2)->getBondDir() == Bond::BondDir::NONE);
    CHECK(mol->getBondWithIdx(3)->getStereo() == Bond::BondStereo::STEREOANY);

    bool clearDoubleBondFlags = true;
    addWavyBondsForStereoAny(*mol, clearDoubleBondFlags,
                             StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO);
    CHECK(mol->getBondWithIdx(2)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(3)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("multiple bonds to wedge") {
    auto mol = "CCC(C)=CC=C(CC)C=CC(C)=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(3)->setStereoAtoms(3, 5);
    mol->getBondWithIdx(3)->setStereo(Bond::BondStereo::STEREOCIS);
    mol->getBondWithIdx(9)->setStereoAtoms(6, 11);
    mol->getBondWithIdx(9)->setStereo(Bond::BondStereo::STEREOANY);
    mol->getBondWithIdx(5)->setStereoAtoms(4, 7);
    mol->getBondWithIdx(5)->setStereo(Bond::BondStereo::STEREOANY);
    addWavyBondsForStereoAny(*mol);
    CHECK(mol->getBondWithIdx(9)->getStereo() == Bond::BondStereo::STEREONONE);
    CHECK(mol->getBondWithIdx(5)->getStereo() == Bond::BondStereo::STEREONONE);
    CHECK(mol->getBondWithIdx(8)->getBondDir() == Bond::BondDir::UNKNOWN);
    for (const auto bond : mol->bonds()) {
      if (bond->getBondType() == Bond::BondType::SINGLE &&
          bond->getIdx() != 8) {
        CHECK(bond->getBondDir() == Bond::BondDir::NONE);
      }
    }
  }
}

TEST_CASE("Github #4215: Ring stereo being discarded in spiro systems") {
  // Note: this bug is still there when using the legacy stereochemistry
  // assignment. It's "non-trivial" to fix there and we've opted not to
  SmilesParserParams ps;
  ps.useLegacyStereo = false;
  SECTION("original failing example") {
    std::unique_ptr<RWMol> m{
        SmilesToMol("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2", ps)};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("original passing example") {
    std::unique_ptr<RWMol> m{SmilesToMol("C[C@H]1CCC2(CC1)CC[C@H](C)CC2", ps)};
    REQUIRE(m);
    // if the middle is unspecified, the two ends can't be specified
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() == Atom::CHI_UNSPECIFIED);

    {
      bool cleanIt = true;
      bool flagPossible = true;
      RWMol m2(*m);
      auto stereoInfo =
          Chirality::findPotentialStereo(m2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 3);
      for (const auto &si : stereoInfo) {
        CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
        CHECK(si.specified == Chirality::StereoSpecified::Unspecified);
        CHECK(si.descriptor == Chirality::StereoDescriptor::None);
      }
    }
    {
      bool cleanIt = true;
      bool flagPossible = false;
      RWMol m2(*m);
      auto stereoInfo =
          Chirality::findPotentialStereo(m2, cleanIt, flagPossible);
      CHECK(stereoInfo.empty());
    }
  }
  SECTION("specified chirality on spiro atom") {
    std::unique_ptr<RWMol> m{
        SmilesToMol("C[C@H]1CC[C@@]2(CC[C@H](C)CC2)CC1", ps)};
    REQUIRE(m);
    // now the middle is specified, so the two ends are as well
    CHECK(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(7)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    {
      bool cleanIt = true;
      bool flagPossible = true;
      RWMol m2(*m);
      auto stereoInfo =
          Chirality::findPotentialStereo(m2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 3);
      for (const auto &si : stereoInfo) {
        CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
        CHECK(si.specified == Chirality::StereoSpecified::Specified);
        CHECK(si.descriptor != Chirality::StereoDescriptor::None);
      }
    }
    {
      bool cleanIt = true;
      bool flagPossible = false;
      RWMol m2(*m);
      auto stereoInfo =
          Chirality::findPotentialStereo(m2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 3);
      for (const auto &si : stereoInfo) {
        CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
        CHECK(si.specified == Chirality::StereoSpecified::Specified);
        CHECK(si.descriptor != Chirality::StereoDescriptor::None);
      }
    }
  }
  SECTION("three spiro rings, unspecified spiro links") {
    std::unique_ptr<RWMol> m{
        SmilesToMol("C[C@H]1CCC2(CC1)CCC1(CC[C@H](C)CC1)CC2", ps)};
    REQUIRE(m);
    {
      bool cleanIt = true;
      bool flagPossible = true;
      RWMol m2(*m);
      auto stereoInfo =
          Chirality::findPotentialStereo(m2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 4);
      for (const auto &si : stereoInfo) {
        CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
        CHECK(si.specified == Chirality::StereoSpecified::Unspecified);
        CHECK(si.descriptor == Chirality::StereoDescriptor::None);
      }
    }
  }
  SECTION("three spiro rings, specified spiro links") {
    std::unique_ptr<RWMol> m{
        SmilesToMol("C[C@H]1CC[C@@]2(CC1)CC[C@]1(CC[C@H](C)CC1)CC2", ps)};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(12)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    {
      bool cleanIt = true;
      bool flagPossible = true;
      RWMol m2(*m);
      auto stereoInfo =
          Chirality::findPotentialStereo(m2, cleanIt, flagPossible);
      CHECK(stereoInfo.size() == 4);
      for (const auto &si : stereoInfo) {
        CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
        CHECK(si.specified == Chirality::StereoSpecified::Specified);
        CHECK(si.descriptor != Chirality::StereoDescriptor::None);
      }
    }
  }
}

TEST_CASE(
    "Github #4279: FindPotentialStereo() doesn't find *marked* ring stereo "
    "when flagPossible=False") {
  SECTION("base") {
    std::unique_ptr<RWMol> m{SmilesToMol("C[C@H]1CC[C@@H](C)CC1")};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    bool cleanIt = true;
    bool flagPossible = false;
    auto stereoInfo = Chirality::findPotentialStereo(*m, cleanIt, flagPossible);
    for (const auto &si : stereoInfo) {
      CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(si.specified == Chirality::StereoSpecified::Specified);
      CHECK(si.descriptor != Chirality::StereoDescriptor::None);
    }
    CHECK(m->getAtomWithIdx(1)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
}

TEST_CASE("StereoInfo comparisons") {
  Chirality::StereoInfo si1;
  si1.centeredOn = 3;
  CHECK(si1.type == Chirality::StereoType::Unspecified);
  si1.type = Chirality::StereoType::Atom_Tetrahedral;
  Chirality::StereoInfo si2;
  si2.centeredOn = 3;
  si2.type = Chirality::StereoType::Atom_Tetrahedral;
  CHECK(si1 == si2);
  si2.descriptor = Chirality::StereoDescriptor::Tet_CCW;
  CHECK(si1 != si2);
}

TEST_CASE("StereoGroup Testing") {
  SECTION("basics") {
    auto mol = "C[C@H](O)[C@@H](C)[C@H](F)Cl |o1:1,3,&2:5,r|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getStereoGroups().size() == 2);
    StereoGroup cp(mol->getStereoGroups()[0]);
    CHECK(cp == mol->getStereoGroups()[0]);
    CHECK(cp != mol->getStereoGroups()[1]);

    std::vector<Atom *> toRemove{mol->getAtomWithIdx(1)};
    std::vector<StereoGroup> &sgs =
        const_cast<std::vector<StereoGroup> &>(mol->getStereoGroups());
    removeGroupsWithAtoms(toRemove, sgs);
    CHECK(mol->getStereoGroups().size() == 1);
  }
}
