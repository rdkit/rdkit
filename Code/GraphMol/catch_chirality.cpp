//
//  Copyright (C) 2020-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cstdlib>

#include <catch2/catch_all.hpp>

#include <boost/noncopyable.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/test_fixtures.h>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>

using namespace RDKit;

unsigned count_wedged_bonds(const ROMol &mol) {
  unsigned nWedged = 0;
  for (const auto bond : mol.bonds()) {
    if (bond->getBondDir() != Bond::BondDir::NONE) {
      ++nWedged;
    }
  }
  return nWedged;
}

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
      CHECK(
          Chirality::detail::isBondPotentialStereoBond(mol->getBondWithIdx(3)));
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
      std::vector<unsigned> catoms = {2, 0, 4, 5};
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
  SECTION("partially specified") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> mol{SmilesToMol("C[C@H]1CCC(C)CC1", ps)};
    REQUIRE(mol);
    // std::cerr << "------------ 3 -------------" << std::endl;
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 1);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);

    CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[1].centeredOn == 4);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
  }
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

TEST_CASE("tricky recursive example from Dan Nealschneider", "[chirality]") {
  SECTION("adapted") {
    auto mol = "CC=C1CCC(O)CC1"_smiles;
    REQUIRE(mol);
    mol->updatePropertyCache();
    MolOps::setBondStereoFromDirections(*mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    REQUIRE(stereoInfo.size() == 2);
    CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
    CHECK(stereoInfo[0].centeredOn == 5);
    CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(stereoInfo[1].type == Chirality::StereoType::Bond_Double);
    CHECK(stereoInfo[1].centeredOn == 1);
    CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
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
  SECTION(
      "Removal of stereoatoms requires removing CIS/TRANS when using legacy stereo") {
    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(false);

    {
      auto mol = "N/C=C/C"_smiles;
      CHECK(mol->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREOTRANS);
      auto rwmol = dynamic_cast<RWMol *>(mol.get());
      rwmol->removeBond(0, 1);
      CHECK(mol->getBondWithIdx(0)->getStereo() ==
            Bond::BondStereo::STEREONONE);
    }
    {
      auto mol = "N/C=C/C"_smiles;
      CHECK(mol->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREOTRANS);
      auto rwmol = dynamic_cast<RWMol *>(mol.get());
      rwmol->removeBond(2, 3);
      CHECK(mol->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREONONE);
    }
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
    auto wedgedBonds = Chirality::pickBondsToWedge(*mol);
    REQUIRE(wedgedBonds.size() == 1);
    auto head = wedgedBonds.begin();
    CHECK(head->first == 3);
    CHECK(head->second->getIdx() == 3);
  }
  SECTION("simplest, specified double bond") {
    auto mol = "OC=C[C@H](C1CC1)C2CCC2"_smiles;
    REQUIRE(mol);
    mol->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    mol->getBondBetweenAtoms(1, 2)->setStereo(Bond::BondStereo::STEREOCIS);
    auto wedgedBonds = Chirality::pickBondsToWedge(*mol);
    REQUIRE(wedgedBonds.size() == 1);
    auto head = wedgedBonds.begin();
    CHECK(head->first == 3);
    CHECK(head->second->getIdx() == 3);
  }
  SECTION("prefer unspecified bond stereo") {
    auto mol = "OC=C[C@H](C=CF)(C=CC)"_smiles;
    REQUIRE(mol);
    mol->getBondBetweenAtoms(1, 2)->setStereoAtoms(0, 3);
    mol->getBondBetweenAtoms(1, 2)->setStereo(Bond::BondStereo::STEREOCIS);
    mol->getBondBetweenAtoms(4, 5)->setStereoAtoms(3, 6);
    mol->getBondBetweenAtoms(4, 5)->setStereo(Bond::BondStereo::STEREOANY);
    auto wedgedBonds = Chirality::pickBondsToWedge(*mol);
    REQUIRE(wedgedBonds.size() == 1);
    auto head = wedgedBonds.begin();
    CHECK(head->first == 6);
    CHECK(head->second->getIdx() == 3);
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
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);

  SECTION("original failing example") {
    auto m = "C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("original failing example, findPossible") {
    UseLegacyStereoPerceptionFixture inner_reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(true);
    SmilesParserParams ps2;
    ps2.sanitize = false;
    std::unique_ptr<RWMol> m{
        SmilesToMol("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2", ps2)};
    REQUIRE(m);
    MolOps::sanitizeMol(*m);
    auto stereoInfo = Chirality::findPotentialStereo(*m, true, true);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(11)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    REQUIRE(stereoInfo.size() == 4);
    for (const auto &si : stereoInfo) {
      CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
      if (si.centeredOn != 9) {
        CHECK(si.specified == Chirality::StereoSpecified::Unspecified);
        CHECK(si.descriptor == Chirality::StereoDescriptor::None);
      } else {
        CHECK(si.specified == Chirality::StereoSpecified::Specified);
        CHECK(si.descriptor != Chirality::StereoDescriptor::None);
      }
    }
  }
  SECTION("original passing example") {
    auto m = "C[C@H]1CCC2(CC1)CC[C@H](C)CC2"_smiles;
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
    auto m = "C[C@H]1CC[C@@]2(CC[C@H](C)CC2)CC1"_smiles;
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
    auto m = "C[C@H]1CCC2(CC1)CCC1(CC[C@H](C)CC1)CC2"_smiles;
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
    auto m = "C[C@H]1CC[C@@]2(CC1)CC[C@]1(CC[C@H](C)CC1)CC2"_smiles;
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
    "Github #4215: Ring stereo being discarded in spiro systems, using deprecated useLegacyStereo option") {
  // Note: this bug is still there when using the legacy stereochemistry
  // assignment. It's "non-trivial" to fix there and we've opted not to
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);
  SECTION("original failing example") {
    std::unique_ptr<RWMol> m{SmilesToMol("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2")};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("original failing example, findPossible") {
    SmilesParserParams ps2;
    ps2.sanitize = false;
    std::unique_ptr<RWMol> m{
        SmilesToMol("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2", ps2)};
    REQUIRE(m);
    MolOps::sanitizeMol(*m);
    auto stereoInfo = Chirality::findPotentialStereo(*m, true, true);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(11)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    REQUIRE(stereoInfo.size() == 4);
    for (const auto &si : stereoInfo) {
      CHECK(si.type == Chirality::StereoType::Atom_Tetrahedral);
      if (si.centeredOn != 9) {
        CHECK(si.specified == Chirality::StereoSpecified::Unspecified);
        CHECK(si.descriptor == Chirality::StereoDescriptor::None);
      } else {
        CHECK(si.specified == Chirality::StereoSpecified::Specified);
        CHECK(si.descriptor != Chirality::StereoDescriptor::None);
      }
    }
  }
  SECTION("original passing example") {
    std::unique_ptr<RWMol> m{SmilesToMol("C[C@H]1CCC2(CC1)CC[C@H](C)CC2")};
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
    std::unique_ptr<RWMol> m{SmilesToMol("C[C@H]1CC[C@@]2(CC[C@H](C)CC2)CC1")};
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
        SmilesToMol("C[C@H]1CCC2(CC1)CCC1(CC[C@H](C)CC1)CC2")};
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
        SmilesToMol("C[C@H]1CC[C@@]2(CC1)CC[C@]1(CC[C@H](C)CC1)CC2")};
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

TEST_CASE("Removing stereogroups from unspecified atoms") {
  SECTION("basics") {
    auto mol = "C[C@](O)(Cl)F |o1:1|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getStereoGroups().size() == 1);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    Chirality::cleanupStereoGroups(*mol);
    CHECK(mol->getStereoGroups().empty());
  }
  SECTION("parsing") {
    auto mol = "C[C@](C)(Cl)F |o1:1|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getStereoGroups().empty());
  }
  SECTION("partial group removal") {
    auto mol = "C[C@](C)(Cl)[C@H](F)Cl |o1:1,4|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getStereoGroups().size() == 1);
    CHECK(mol->getStereoGroups()[0].getAtoms().size() == 1);
    CHECK(mol->getStereoGroups()[0].getAtoms()[0]->getIdx() == 4);
  }
}

TEST_CASE("replaceAtom and StereoGroups") {
  SECTION("basics") {
    auto mol = "C[C@](O)(Cl)[C@H](F)Cl |o1:1,4|"_smiles;
    REQUIRE(mol);
    CHECK(mol->getStereoGroups().size() == 1);
    CHECK(mol->getStereoGroups()[0].getAtoms().size() == 2);
    CHECK(mol->getStereoGroups()[0].getAtoms()[0] == mol->getAtomWithIdx(1));

    Atom acp(*mol->getAtomWithIdx(1));
    mol->replaceAtom(1, &acp);
    CHECK(mol->getStereoGroups().size() == 1);
    CHECK(mol->getStereoGroups()[0].getAtoms().size() == 2);
    CHECK(mol->getStereoGroups()[0].getAtoms()[0] == mol->getAtomWithIdx(1));
  }
}

TEST_CASE(
    "Github #5200: FindPotentialStereo does not clean stereoflags from atoms "
    "which cannot be stereocenters") {
  auto m = "CCF"_smiles;
  REQUIRE(m);
  m->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  bool cleanIt = true;
  auto sinfo = Chirality::findPotentialStereo(*m, cleanIt);
  CHECK(sinfo.empty());
  CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
        Atom::ChiralType::CHI_UNSPECIFIED);
}

TEST_CASE(
    "Github #5196: Zero & coordinate bonds are being taken into account for "
    "chirality") {
  RDLog::LogStateSetter setter;  // disable irritating warning messages
  auto mol = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 15 18 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.136359 0.025241 -0.986870 0
M  V30 2 C 0.211183 -0.810922 0.138318 0
M  V30 3 C -0.446638 -0.713741 1.305561 0
M  V30 4 C -1.141107 0.914647 -0.916429 0
M  V30 5 R -1.628248 -0.983190 -0.411960 0
M  V30 6 H 0.392055 -0.106505 -1.920607 0
M  V30 7 H 0.974038 -1.568492 0.017171 0
M  V30 8 H -0.209921 -1.406535 2.084966 0
M  V30 9 H -1.378909 1.482059 -1.807349 0
M  V30 10 C -1.544607 0.306162 1.588191 0
M  V30 11 C -1.946856 1.186683 0.358271 0
M  V30 12 H -1.207983 0.944410 2.407927 0
M  V30 13 H -2.419549 -0.225146 1.965589 0
M  V30 14 H -3.006492 1.040978 0.144313 0
M  V30 15 H -1.830875 2.240146 0.620809 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 2 3 2
M  V30 3 2 4 1
M  V30 4 0 1 5
M  V30 5 0 2 5
M  V30 6 0 3 5
M  V30 7 0 4 5
M  V30 8 1 1 6
M  V30 9 1 2 7
M  V30 10 1 3 8
M  V30 11 1 4 9
M  V30 12 1 10 3
M  V30 13 1 11 4
M  V30 14 1 11 10
M  V30 15 1 12 10
M  V30 16 1 13 10
M  V30 17 1 14 11
M  V30 18 1 15 11
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(mol);
  SECTION("as reported") {
    MolOps::assignStereochemistryFrom3D(*mol);
    for (auto aidx : {0, 1, 2, 3}) {
      CHECK(mol->getAtomWithIdx(aidx)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("as reported - ZOBs") {
    for (auto bidx : {3, 4, 5, 6}) {
      mol->getBondWithIdx(bidx)->setBondType(Bond::BondType::ZERO);
    }
    MolOps::assignStereochemistryFrom3D(*mol);
    for (auto idx : {0, 1, 2, 3}) {
      CHECK(mol->getAtomWithIdx(idx)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("as reported - datives") {
    for (auto bidx : {3, 4, 5, 6}) {
      mol->getBondWithIdx(bidx)->setBondType(Bond::BondType::DATIVE);
    }
    MolOps::assignStereochemistryFrom3D(*mol);
    for (auto idx : {0, 1, 2, 3}) {
      CHECK(mol->getAtomWithIdx(idx)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("as reported - reversed datives") {
    // structure is bogus, but we want to test
    for (auto bidx : {3, 4, 5, 6}) {
      auto bond = mol->getBondWithIdx(bidx);
      bond->setEndAtomIdx(bond->getBeginAtomIdx());
      bond->setBeginAtomIdx(4);
      bond->setBondType(Bond::BondType::DATIVE);
    }
    MolOps::assignStereochemistryFrom3D(*mol);
    for (auto idx : {0, 1, 2, 3}) {
      CHECK(mol->getAtomWithIdx(idx)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("as reported - singles") {
    // structure is bogus, but we want to test
    for (auto bidx : {3, 4, 5, 6}) {
      mol->getBondWithIdx(bidx)->setBondType(Bond::BondType::SINGLE);
    }
    MolOps::assignStereochemistryFrom3D(*mol);
    for (auto idx : {0, 1, 2, 3}) {
      CHECK(mol->getAtomWithIdx(idx)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("assignStereochemistry") {
    auto mol = "[Fe]C(=C)O |C:1.0|"_smiles;
    REQUIRE(mol);
    for (auto bt : {Bond::BondType::DATIVE, Bond::BondType::ZERO,
                    Bond::BondType::UNSPECIFIED}) {
      mol->getAtomWithIdx(1)->setChiralTag(
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
      mol->getBondWithIdx(0)->setBondType(bt);
      bool cleanit = true;
      bool force = true;
      MolOps::assignStereochemistry(*mol, cleanit, force);
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("isAtomPotentialTetrahedralCenter() and getStereoInfo()") {
    auto mol = "[Fe]C(=C)O |C:1.0|"_smiles;
    REQUIRE(mol);
    for (auto bt : {Bond::BondType::DATIVE, Bond::BondType::ZERO,
                    Bond::BondType::UNSPECIFIED}) {
      mol->getAtomWithIdx(1)->setChiralTag(
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
      mol->getBondWithIdx(0)->setBondType(bt);
      CHECK(!Chirality::detail::isAtomPotentialStereoAtom(
          mol->getAtomWithIdx(1)));
      bool cleanit = true;
      auto sinfo = Chirality::findPotentialStereo(*mol, cleanit);
      CHECK(sinfo.empty());
      CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
}

TEST_CASE(
    "Github #5239: Precondition violation on chiral Atoms with zero order "
    "bonds") {
  RDLog::LogStateSetter setter;  // disable irritating warning messages
  auto molblock = R"CTAB(
     RDKit          3D
     
  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.446600 -0.713700 1.305600 0
M  V30 2 Fe -1.628200 -0.983200 -0.412000 0
M  V30 3 Cl -0.049300 -1.876700 2.613900 0
M  V30 4 C -1.544600 0.306200 1.588200 0
M  V30 5 F 0.673700 0.029200 0.993700 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3
M  V30 2 1 1 4 CFG=1
M  V30 3 1 1 5
M  V30 4 0 2 1
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB";
  bool sanitize = false;
  std::unique_ptr<ROMol> mol(MolBlockToMol(molblock, sanitize));
  REQUIRE(mol);
  MolOps::assignStereochemistryFrom3D(*mol);

  CHECK(mol->getAtomWithIdx(0)->getChiralTag() !=
        Atom::ChiralType::CHI_UNSPECIFIED);
}
TEST_CASE("nontetrahedral stereo from 3D", "[nontetrahedral]") {
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/test_data/nontetrahedral_3d.sdf";
  SECTION("basics") {
    SDMolSupplier suppl(pathName);
    while (!suppl.atEnd()) {
      std::unique_ptr<ROMol> m{suppl.next()};
      REQUIRE(m);
      MolOps::assignChiralTypesFrom3D(*m);
      auto ct = m->getProp<std::string>("ChiralType");
      auto cp = m->getProp<unsigned>("ChiralPermutation");
      auto atom = m->getAtomWithIdx(0);

      if (ct == "SP") {
        CHECK(atom->getChiralTag() == Atom::ChiralType::CHI_SQUAREPLANAR);
      } else if (ct == "TB") {
        CHECK(atom->getChiralTag() ==
              Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
      } else if (ct == "TH") {
        CHECK(atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL);
      } else if (ct == "OH") {
        CHECK(atom->getChiralTag() == Atom::ChiralType::CHI_OCTAHEDRAL);
      }

      CHECK(atom->getProp<unsigned>(common_properties::_chiralPermutation) ==
            cp);
    }
  }
  SECTION("disable nontetrahedral stereo") {
    AllowNontetrahedralChiralityFixture reset_non_tetrahedral_allowed;
    Chirality::setAllowNontetrahedralChirality(false);
    SDMolSupplier suppl(pathName);
    while (!suppl.atEnd()) {
      std::unique_ptr<ROMol> m{suppl.next()};
      REQUIRE(m);
      MolOps::assignChiralTypesFrom3D(*m);
      auto ct = m->getProp<std::string>("ChiralType");
      auto atom = m->getAtomWithIdx(0);

      if (ct == "TH") {
        CHECK(atom->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL);
      } else {
        CHECK(atom->getChiralTag() == Atom::ChiralType::CHI_UNSPECIFIED);
      }
    }
  }
}

TEST_CASE("assignStereochemistry shouldn't remove nontetrahedral stereo",
          "[nontetrahedral]") {
  SECTION("basics") {
    SmilesParserParams parseps;
    parseps.sanitize = false;
    parseps.removeHs = false;
    std::unique_ptr<RWMol> m{SmilesToMol("F[Pt@TB1](O)(Br)(N)Cl", parseps)};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
    bool cleanIt = true;
    bool force = true;
    MolOps::assignStereochemistry(*m, cleanIt, force);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
  }
  SECTION("standard SMILES parsing") {
    std::unique_ptr<RWMol> m{SmilesToMol("F[Pt@TB1](O)(Br)(N)Cl")};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
  }
  SECTION("SMILES parsing w/o sanitization") {
    SmilesParserParams parseps;
    // we need to skip stereo assignment
    parseps.sanitize = false;
    std::unique_ptr<RWMol> m{SmilesToMol("F[Pt@TB1](O)(Br)(N)Cl", parseps)};
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
  }
}

TEST_CASE("remove hs and non-tetrahedral stereo", "[nontetrahedral]") {
  SmilesParserParams parseps;
  parseps.sanitize = false;
  parseps.removeHs = false;
  std::vector<std::string> smiles = {"F[Pt@TB1]([H])(Br)(N)Cl",
                                     "F[Pt@TB]([H])(Br)(N)Cl"};
  for (const auto &smi : smiles) {
    std::unique_ptr<RWMol> m{SmilesToMol(smi, parseps)};
    REQUIRE(m);
    CHECK(m->getNumAtoms(6));
    {
      // the default is to not remove Hs to non-tetrahedral stereocenters
      RWMol molcp(*m);
      MolOps::removeHs(molcp);
      CHECK(molcp.getNumAtoms() == 6);
    }
    {
      // but we can enable it
      RWMol molcp(*m);
      MolOps::RemoveHsParameters ps;
      ps.removeNontetrahedralNeighbors = true;
      MolOps::removeHs(molcp, ps);
      CHECK(molcp.getNumAtoms() == 5);
    }
    {
      // but we can enable it
      RWMol molcp(*m);
      MolOps::removeAllHs(molcp);
      CHECK(molcp.getNumAtoms() == 5);
    }
  }
}

TEST_CASE("getIdealAngle", "[nontetrahedral]") {
  SECTION("TB1") {
    auto m = "S[As@TB1](F)(Cl)(Br)N"_smiles;
    REQUIRE(m);
    CHECK(Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                    m->getAtomWithIdx(0)));
    CHECK(Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                    m->getAtomWithIdx(5)));
    CHECK(!Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     m->getAtomWithIdx(2)));
    CHECK(!Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     m->getAtomWithIdx(3)));
    CHECK(!Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     m->getAtomWithIdx(4)));
    CHECK(Chirality::getTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1))
              ->getIdx() == 0);
    CHECK(Chirality::getTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1), -1)
              ->getIdx() == 5);

    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(4)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(2), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(120, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(5)),
        Catch::Matchers::WithinAbs(180, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(5), m->getAtomWithIdx(2)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(5), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(5), m->getAtomWithIdx(4)),
        Catch::Matchers::WithinAbs(90, 0.001));
  }
  SECTION("TB1 missing 1") {
    auto m = "S[As@TB1](F)(Cl)Br"_smiles;
    REQUIRE(m);

    CHECK(Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                    m->getAtomWithIdx(0)));
    CHECK(!Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     m->getAtomWithIdx(2)));
    CHECK(!Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     m->getAtomWithIdx(3)));
    CHECK(!Chirality::isTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     m->getAtomWithIdx(4)));
    CHECK(Chirality::getTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1))
              ->getIdx() == 0);
    CHECK(Chirality::getTrigonalBipyramidalAxialAtom(m->getAtomWithIdx(1),
                                                     -1) == nullptr);

    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(4)),
        Catch::Matchers::WithinAbs(90, 0.001));
    CHECK_THAT(
        Chirality::getIdealAngleBetweenLigands(
            m->getAtomWithIdx(1), m->getAtomWithIdx(2), m->getAtomWithIdx(3)),
        Catch::Matchers::WithinAbs(120, 0.001));
  }
}

TEST_CASE("getChiralPermutation", "[nontetrahedral]") {
  SECTION("TB1") {
    // clang-format off
    std::vector<std::pair<std::list<int>, unsigned int>> data = {
        {{2, 3, 4, 5, 6}, 1},
        {{2, 3, 5, 4, 6}, 2},

        {{2, 3, 4, 6, 5}, 3},
        {{2, 3, 5, 6, 4}, 4},

        {{2, 3, 6, 4, 5}, 5},
        {{2, 3, 6, 5, 4}, 6},

        {{2, 6, 3, 4, 5}, 7},
        {{2, 6, 3, 5, 4}, 8},

        {{3, 2, 4, 5, 6}, 9},
        {{3, 2, 5, 4, 6}, 11},

        {{3, 2, 4, 6, 5}, 10},
        {{3, 2, 5, 6, 4}, 12},

        {{3, 2, 6, 4, 5}, 13},
        {{3, 2, 6, 5, 4}, 14},

        {{3, 4, 2, 5, 6}, 15},
        {{3, 5, 2, 4, 6}, 20},

        {{3, 4, 2, 6, 5}, 16},
        {{3, 5, 2, 6, 4}, 19},

        {{3, 4, 5, 2, 6}, 17},
        {{3, 5, 4, 2, 6}, 18},


    };
    // clang-format on
    auto m = "CCS[As@TB1](F)(Cl)(Br)N"_smiles;
    REQUIRE(m);
    const auto atm = m->getAtomWithIdx(3);
    for (const auto &pr : data) {
      // std::cerr << "---- " << pr.second << std::endl;
      CHECK(Chirality::getChiralPermutation(atm, pr.first) == pr.second);
    }
  }

  SECTION("SP1") {
    // clang-format off
    std::vector<std::pair<std::list<int>, unsigned int>> data = {
        {{2, 3, 4, 5}, 1},
        {{2, 4, 3, 5}, 2},
        {{2, 3, 5, 4}, 3},
    };
    // clang-format on
    auto m = "CCC[Pt@SP1](F)(Cl)N"_smiles;
    REQUIRE(m);
    const auto atm = m->getAtomWithIdx(3);
    for (const auto &pr : data) {
      // std::cerr << "---- " << pr.second << std::endl;
      CHECK(Chirality::getChiralPermutation(atm, pr.first) == pr.second);
    }
  }
  SECTION("OH1") {
    // clang-format off
    std::vector<std::pair<std::list<int>, unsigned int>> data = {
        {{2, 3, 4, 5, 6, 7}, 1},
        {{2, 3, 6, 5, 4, 7}, 2},

        {{2, 3, 4, 5, 7, 6}, 3},
        {{2, 3, 6, 5, 7, 4}, 16},

        {{2, 3, 4, 7, 5, 6}, 6},
        {{2, 3, 6, 7, 5, 4}, 18},

        {{2, 3, 7, 4, 5, 6}, 19},
        {{2, 3, 7, 6, 5, 4}, 24},

        {{2, 7, 3, 4, 5, 6}, 25},
        {{2, 7, 3, 6, 5, 4}, 30},

        {{2, 3, 4, 6, 5, 7}, 4},
        {{2, 3, 6, 4, 5, 7}, 14},

        {{2, 3, 4, 6, 7, 5}, 5},
        {{2, 3, 6, 4, 7, 5}, 15},

        {{2, 3, 4, 7, 6, 5}, 7},
        {{2, 3, 6, 7, 4, 5}, 17},

        {{2, 3, 7, 4, 6, 5}, 20},
        {{2, 3, 7, 6, 4, 5}, 23},

        {{2, 7, 3, 4, 6, 5}, 26},
        {{2, 7, 3, 6, 4, 5}, 29},

        {{2, 3, 5, 6, 4, 7}, 10},
        {{2, 3, 5, 4, 6, 7}, 8},

        {{2, 3, 5, 6, 7, 4}, 11},
        {{2, 3, 5, 4, 7, 6}, 9},

        {{2, 3, 5, 7, 6, 4}, 13},
        {{2, 3, 5, 7, 4, 6}, 12},

        {{2, 3, 7, 5, 6, 4}, 22},
        {{2, 3, 7, 5, 4, 6}, 21},

        {{2, 7, 3, 5, 6, 4}, 28},
        {{2, 7, 3, 5, 4, 6}, 27},

    };
    // clang-format on
    auto m = "CCO[Co@OH1](Cl)(C)(N)(F)P"_smiles;
    REQUIRE(m);
    const auto atm = m->getAtomWithIdx(3);
    for (const auto &pr : data) {
      // std::cerr << "---- " << pr.second << std::endl;
      CHECK(Chirality::getChiralPermutation(atm, pr.first) == pr.second);
    }
  }
}

TEST_CASE("isAtomPotentialNontetrahedralCenter", "[nontetrahedral]") {
  SECTION("basics") {
    {
      auto mol = "C[S+](O)F"_smiles;
      REQUIRE(mol);
      CHECK(!Chirality::detail::isAtomPotentialNontetrahedralCenter(
          mol->getAtomWithIdx(1)));
    }
    {
      auto mol = "C[SH](O)F"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialNontetrahedralCenter(
          mol->getAtomWithIdx(1)));
    }
    {
      auto mol = "C[S@SP](O)F"_smiles;
      REQUIRE(mol);
      CHECK(Chirality::detail::isAtomPotentialNontetrahedralCenter(
          mol->getAtomWithIdx(1)));
    }
  }
}
TEST_CASE("nontetrahedral StereoInfo", "[nontetrahedral]") {
  SECTION("SP") {
    auto m = "C[Pt@SP1](F)(Cl)O"_smiles;
    REQUIRE(m);
    auto sinfo = Chirality::findPotentialStereo(*m);
    REQUIRE(sinfo.size() == 1);
    CHECK(sinfo[0].centeredOn == 1);
    CHECK(sinfo[0].type == Chirality::StereoType::Atom_SquarePlanar);
    CHECK(sinfo[0].descriptor == Chirality::StereoDescriptor::None);
    CHECK(sinfo[0].permutation == 1);
    CHECK(sinfo[0].controllingAtoms == std::vector<unsigned int>{0, 2, 3, 4});
  }
  SECTION("TB") {
    auto m = "C[Fe@TB4](F)(Cl)(O)N"_smiles;
    REQUIRE(m);
    auto sinfo = Chirality::findPotentialStereo(*m);
    REQUIRE(sinfo.size() == 1);
    CHECK(sinfo[0].centeredOn == 1);
    CHECK(sinfo[0].type == Chirality::StereoType::Atom_TrigonalBipyramidal);
    CHECK(sinfo[0].descriptor == Chirality::StereoDescriptor::None);
    CHECK(sinfo[0].permutation == 4);
    CHECK(sinfo[0].controllingAtoms ==
          std::vector<unsigned int>{0, 2, 3, 4, 5});
  }
  SECTION("TB0") {
    auto m = "C[Fe@TB](F)(Cl)(O)N"_smiles;
    REQUIRE(m);
    auto sinfo = Chirality::findPotentialStereo(*m);
    REQUIRE(sinfo.size() == 1);
    CHECK(sinfo[0].centeredOn == 1);
    CHECK(sinfo[0].specified == Chirality::StereoSpecified::Unknown);
    CHECK(sinfo[0].type == Chirality::StereoType::Atom_TrigonalBipyramidal);
    CHECK(sinfo[0].descriptor == Chirality::StereoDescriptor::None);
    CHECK(sinfo[0].permutation == 0);
    CHECK(sinfo[0].controllingAtoms ==
          std::vector<unsigned int>{0, 2, 3, 4, 5});
  }
  SECTION("perceived as possible") {
    auto m = "C[Fe](F)(Cl)(O)N"_smiles;
    REQUIRE(m);
    auto sinfo = Chirality::findPotentialStereo(*m);
    REQUIRE(sinfo.size() == 1);
    CHECK(sinfo[0].centeredOn == 1);
    CHECK(sinfo[0].specified == Chirality::StereoSpecified::Unspecified);
    CHECK(sinfo[0].type == Chirality::StereoType::Atom_TrigonalBipyramidal);
    CHECK(sinfo[0].descriptor == Chirality::StereoDescriptor::None);
    CHECK(sinfo[0].permutation == 0);
    CHECK(sinfo[0].controllingAtoms ==
          std::vector<unsigned int>{0, 2, 3, 4, 5});
  }

  SECTION("OH") {
    auto m = "C[Fe@OH9](F)(Cl)(O)(N)Br"_smiles;
    REQUIRE(m);
    auto sinfo = Chirality::findPotentialStereo(*m);
    REQUIRE(sinfo.size() == 1);
    CHECK(sinfo[0].centeredOn == 1);
    CHECK(sinfo[0].type == Chirality::StereoType::Atom_Octahedral);
    CHECK(sinfo[0].descriptor == Chirality::StereoDescriptor::None);
    CHECK(sinfo[0].permutation == 9);
    CHECK(sinfo[0].controllingAtoms ==
          std::vector<unsigned int>{0, 2, 3, 4, 5, 6});
  }
  SECTION("OH missing ligand") {
    auto m = "C[Fe@OH9](F)(Cl)(O)N"_smiles;
    REQUIRE(m);
    auto sinfo = Chirality::findPotentialStereo(*m);
    REQUIRE(sinfo.size() == 1);
    CHECK(sinfo[0].centeredOn == 1);
    CHECK(sinfo[0].type == Chirality::StereoType::Atom_Octahedral);
    CHECK(sinfo[0].descriptor == Chirality::StereoDescriptor::None);
    CHECK(sinfo[0].permutation == 9);
    CHECK(sinfo[0].controllingAtoms ==
          std::vector<unsigned int>{0, 2, 3, 4, 5});
  }
}

TEST_CASE("github #5328: assignChiralTypesFrom3D() ignores wiggly bonds") {
  SECTION("basics") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.900794 -0.086835 0.009340 0
M  V30 2 C -0.552652 0.319534 0.077502 0
M  V30 3 F -0.861497 0.413307 1.437370 0
M  V30 4 Cl -0.784572 1.925710 -0.672698 0
M  V30 5 O -1.402227 -0.583223 -0.509512 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=2
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    MolOps::assignChiralTypesFrom3D(*m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("non-tetrahedral") {
    auto m = R"CTAB(
  Mrv2108 05252216313D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.7191 0.2488 -3.5085 0
M  V30 2 As -1.0558 1.9209 -2.6345 0
M  V30 3 F -0.4636 3.422 -1.7567 0
M  V30 4 O -2.808 2.4243 -2.1757 0
M  V30 5 Cl -0.1145 2.6609 -4.5048 0
M  V30 6 Br 0.2255 0.6458 -1.079 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5 CFG=2
M  V30 5 1 2 6
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    MolOps::assignChiralTypesFrom3D(*m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
}

TEST_CASE("useLegacyStereoPerception feature flag") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;

  SECTION("original failing example") {
    Chirality::setUseLegacyStereoPerception(true);
    auto m = "C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  SECTION("use new code") {
    Chirality::setUseLegacyStereoPerception(false);
    auto m = "C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2"_smiles;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
  std::string molblock = R"CTAB(
  Mrv2108 05202206352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.5417 3.165 0 0
M  V30 2 C -5.8753 2.395 0 0
M  V30 3 C -5.8753 0.855 0 0
M  V30 4 C -4.5417 0.085 0 0 CFG=1
M  V30 5 C -3.208 0.855 0 0
M  V30 6 C -3.208 2.395 0 0
M  V30 7 C -4.5417 -1.455 0 0
M  V30 8 C -1.8743 0.085 0 0
M  V30 9 C -4.5417 6.2451 0 0 CFG=2
M  V30 10 C -5.8753 5.4751 0 0
M  V30 11 C -5.8753 3.9351 0 0
M  V30 12 C -3.208 3.9351 0 0
M  V30 13 C -3.208 5.4751 0 0
M  V30 14 C -4.5417 7.7851 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 1 6
M  V30 7 1 4 7 CFG=1
M  V30 8 1 5 8
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 12 13
M  V30 12 1 9 13
M  V30 13 1 11 1
M  V30 14 1 1 12
M  V30 15 1 9 14 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
  SECTION("original example, from mol block") {
    Chirality::setUseLegacyStereoPerception(true);
    std::unique_ptr<RWMol> m{MolBlockToMol(molblock)};
    CHECK(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  }
  SECTION("original example, from mol block, new code") {
    Chirality::setUseLegacyStereoPerception(false);
    std::unique_ptr<RWMol> m{MolBlockToMol(molblock)};
    CHECK(m->getAtomWithIdx(3)->getChiralTag() != Atom::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_UNSPECIFIED);
  }
}

TEST_CASE("wedgeMolBonds to aromatic rings") {
  auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 11 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.948889 -2.986305 0.000000 0
M  V30 2 C 2.560660 -4.435194 0.000000 0
M  V30 3 N 1.111771 -4.046965 0.000000 0
M  V30 4 C 1.500000 -2.598076 0.000000 0
M  V30 5 C 0.750000 -1.299038 0.000000 0
M  V30 6 C 1.500000 0.000000 0.000000 0
M  V30 7 C 0.750000 1.299038 0.000000 0
M  V30 8 C -0.750000 1.299038 0.000000 0
M  V30 9 C -1.500000 0.000000 0.000000 0
M  V30 10 C -0.750000 -1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5 CFG=1
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 2 9 10
M  V30 10 1 4 1
M  V30 11 1 10 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
  REQUIRE(m);
  CHECK(m->getAtomWithIdx(3)->getChiralTag() !=
        Atom::ChiralType::CHI_UNSPECIFIED);
  CHECK(m->getBondWithIdx(2)->getBondDir() == Bond::BondDir::NONE);
  CHECK(m->getBondWithIdx(3)->getBondDir() == Bond::BondDir::NONE);

  SECTION("generating mol blocks") {
    auto mb = MolToV3KMolBlock(*m);
    CHECK(mb.find("M  V30 10 1 4 1 CFG=1") == std::string::npos);
    CHECK(mb.find("M  V30 4 1 4 5 CFG=1") != std::string::npos);
  }

  SECTION("details: pickBondsWedge()") {
    // this is with aromatic bonds
    auto wedgedBonds = Chirality::pickBondsToWedge(*m);
    CHECK(wedgedBonds.at(3)->getIdx() == 3);
    RWMol cp(*m);

    // now try kekulized:
    MolOps::Kekulize(cp);
    wedgedBonds = Chirality::pickBondsToWedge(cp);
    CHECK(wedgedBonds.at(3)->getIdx() == 3);
  }
}

TEST_CASE("github 5307: AssignAtomChiralTagsFromStructure ignores Hydrogens") {
  std::string mb = R"CTAB(
     RDKit          3D
     
  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.022097 0.003215 0.016520 0
M  V30 2 H -0.669009 0.889360 -0.100909 0
M  V30 3 H -0.377788 -0.857752 -0.588296 0
M  V30 4 H 0.096421 -0.315125 1.063781 0
M  V30 5 H 0.972473 0.280302 -0.391096 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
  bool sanitize = true;
  bool removeHs = false;
  std::unique_ptr<RWMol> m{MolBlockToMol(mb, sanitize, removeHs)};
  REQUIRE(m);
  MolOps::assignChiralTypesFrom3D(*m);
  CHECK(m->getAtomWithIdx(0)->getChiralTag() !=
        Atom::ChiralType::CHI_UNSPECIFIED);
  // assignStereochemistryFrom3D() actually checks:
  MolOps::assignStereochemistryFrom3D(*m);
  CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
        Atom::ChiralType::CHI_UNSPECIFIED);
}

TEST_CASE(
    "github 5403: Incorrect perception of pseudoasymmetric centers on "
    "non-canonical molecules") {
  SECTION("basics") {
    auto mol1 = "N[C@@]12CC[C@@](CC1)(C2)C(F)F"_smiles;
    REQUIRE(mol1);

    bool clean = true;

    auto stereoInfo1 = Chirality::findPotentialStereo(*mol1, clean);
    REQUIRE(stereoInfo1.size() == 2);
    REQUIRE(stereoInfo1[0].centeredOn == 1);
    REQUIRE(stereoInfo1[1].centeredOn == 4);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol1->getAtomWithIdx(4)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);

    auto mol2 = "C1C[C@]2(C(F)F)CC[C@@]1(N)C2"_smiles;
    REQUIRE(mol2);

    auto stereoInfo2 = Chirality::findPotentialStereo(*mol2, clean);
    REQUIRE(stereoInfo2.size() == 2);
    CHECK(stereoInfo2[0].centeredOn == 2);
    CHECK(stereoInfo2[1].centeredOn == 8);
    {
      RWMol cp(*mol2);
      Chirality::findPotentialStereo(cp, true, false);
      CHECK(cp.getAtomWithIdx(2)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(cp.getAtomWithIdx(8)->getChiralTag() !=
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
  SECTION("new stereo") {
    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(false);

    auto mol1 = "N[C@@]12CC[C@@](CC1)(C2)C(F)F"_smiles;
    REQUIRE(mol1);
    CHECK(mol1->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(mol1->getAtomWithIdx(4)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
}
TEST_CASE("assignStereochemistry sets bond stereo with new stereo perception") {
  SECTION("basics") {
    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(false);
    {
      auto m = "C/C=C/C"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOTRANS);
      CHECK(m->getBondWithIdx(1)->getStereoAtoms() == std::vector<int>{0, 3});
    }
    {
      auto m = "C/C=C\\C"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOCIS);
      CHECK(m->getBondWithIdx(1)->getStereoAtoms() == std::vector<int>{0, 3});
    }
    {
      auto m = "C(/C)=C/C"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOCIS);
      CHECK(m->getBondWithIdx(1)->getStereoAtoms() == std::vector<int>{1, 3});
    }
    {
      auto m = "FC(/C)=C/C"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(2)->getStereo() == Bond::BondStereo::STEREOTRANS);
      CHECK(m->getBondWithIdx(2)->getStereoAtoms() == std::vector<int>{0, 4});
    }
    {
      auto m = "FC(/C)=C(F)/C"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(2)->getStereo() == Bond::BondStereo::STEREOCIS);
      CHECK(m->getBondWithIdx(2)->getStereoAtoms() == std::vector<int>{0, 4});
    }
  }
}

TEST_CASE("chiral duplicates") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);
  SECTION("atom basics") {
    auto mol = "C[C@](F)([C@H](F)Cl)[C@H](F)Cl"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(3)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(mol->getAtomWithIdx(6)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(mol->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("double bonds and atoms") {
    auto mol = "C/C(O)=C([C@H](F)Cl)/[C@H](F)Cl"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(4)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(mol->getAtomWithIdx(7)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(mol->getBondWithIdx(2)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("double bonds and atoms 2") {
    auto mol = "C/C(O)=C([C@H](F)Cl)/[C@@H](F)Cl"_smiles;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(4)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(mol->getAtomWithIdx(7)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    CHECK(mol->getBondWithIdx(2)->getStereo() != Bond::BondStereo::STEREONONE);
  }
  SECTION("double bonds and double bonds") {
    auto mol = "C/C(O)=C(/C=C/C)/C=C/C"_smiles;
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(4)->getStereo() == Bond::BondStereo::STEREOTRANS);
    CHECK(mol->getBondWithIdx(7)->getStereo() == Bond::BondStereo::STEREOTRANS);
    CHECK(mol->getBondWithIdx(2)->getStereo() == Bond::BondStereo::STEREONONE);
  }
  SECTION("double bonds and double bonds 2") {
    auto mol = "C/C(O)=C(/C=C/C)/C=C\\C"_smiles;
    REQUIRE(mol);
    CHECK(mol->getBondWithIdx(4)->getStereo() == Bond::BondStereo::STEREOTRANS);
    CHECK(mol->getBondWithIdx(7)->getStereo() == Bond::BondStereo::STEREOCIS);
    CHECK(mol->getBondWithIdx(2)->getStereo() != Bond::BondStereo::STEREONONE);
  }
}

TEST_CASE("more findPotential") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);
  SECTION("basics") {
    {
      auto m = "O[C@H](C)CC(C)C[C@@H](C)O"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 2);
    }
    {
      auto m = "O[C@H](C)CC(C)C[C@H](C)O"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
    {
      auto m = "O[CH](C)C[C@H](C)C[CH](C)O"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
    {
      auto m = "O[CH](C)C[CH](C)C[CH](C)O"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
  }
  SECTION("double bond impact on atoms") {
    {
      auto m = "C[CH](/C=C/C)/C=C\\C"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
    {
      auto m = "C[CH](/C=C/C)/C=C/C"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 2);
    }
    {
      auto m = "C[CH](C=CC)C=CC"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
    {
      auto m = "C[CH](C=CC)C=CC"_smiles;
      REQUIRE(m);
      m->getBondWithIdx(2)->setStereo(Bond::BondStereo::STEREOANY);
      m->getBondWithIdx(5)->setStereo(Bond::BondStereo::STEREOANY);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
  }
  SECTION("atom impact on double bonds") {
    {
      auto m = "CC=C([C@H](F)Cl)[C@@H](F)Cl"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
    {
      auto m = "CC=C([C@H](F)Cl)[C@H](F)Cl"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 2);
    }
    {
      auto m = "CC=C([CH](F)Cl)[CH](F)Cl"_smiles;
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
    {
      auto m = "CC=C([CH](F)Cl)[CH](F)Cl"_smiles;
      REQUIRE(m);
      m->getBondBetweenAtoms(2, 3)->setBondDir(Bond::UNKNOWN);
      m->getBondBetweenAtoms(2, 6)->setBondDir(Bond::UNKNOWN);
      auto si = Chirality::findPotentialStereo(*m, true, true);
      CHECK(si.size() == 3);
    }
  }
}
TEST_CASE("more findPotential and ring stereo") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);
  SECTION("simple") {
    {
      auto m = "CC1CCC(C)CC1"_smiles;
      REQUIRE(m);
      auto stereoInfo = Chirality::findPotentialStereo(*m, true, true);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].centeredOn == 1);
      CHECK(stereoInfo[0].specified == Chirality::StereoSpecified::Unspecified);

      CHECK(stereoInfo[1].type == Chirality::StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[1].centeredOn == 4);
      CHECK(stereoInfo[1].specified == Chirality::StereoSpecified::Unspecified);
    }
  }
}

TEST_CASE(
    "github 2984: RDKit misplaces stereochemistry/chirality information for "
    "small ring") {
  using namespace RDKit::Chirality;

  UseLegacyStereoPerceptionFixture reset_stereo_perception;

  SECTION("The mol with the issue") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      // With Legacy stereo, parsing the SMILES string will discard the
      // problematic stereo features, so findPotentialStereo (which uses
      // the new algorithm) will still find them, but they will be unspecified.
      // Parsing with the new algorithm will preserve the features, so they
      // can be correctly resolved.

      auto specified_status = use_legacy ? StereoSpecified::Unspecified
                                         : StereoSpecified::Specified;

      auto mol = R"SMI(CC/C=C1\C[C@H](O)C1)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].centeredOn == 5);
      CHECK(stereoInfo[0].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == specified_status);
      CHECK(stereoInfo[1].centeredOn == 2);
      CHECK(stereoInfo[1].type == StereoType::Bond_Double);
      CHECK(stereoInfo[1].specified == specified_status);
    }
  }

  // The other sections should yield the same results independently
  // of which stereo algoritm is used, new or legacy

  SECTION("Unspecified, but still potentially stereo") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      auto mol = R"SMI(CCC=C1CC(O)C1)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].centeredOn == 5);
      CHECK(stereoInfo[0].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == StereoSpecified::Unspecified);
      CHECK(stereoInfo[1].centeredOn == 2);
      CHECK(stereoInfo[1].type == StereoType::Bond_Double);
      CHECK(stereoInfo[1].specified == StereoSpecified::Unspecified);
    }
  }

  SECTION("Specified") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      auto mol = R"SMI(CC/C=C1\CC[C@H]1O)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 2);
      CHECK(stereoInfo[0].centeredOn == 6);
      CHECK(stereoInfo[0].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == StereoSpecified::Specified);
      CHECK(stereoInfo[1].centeredOn == 2);
      CHECK(stereoInfo[1].type == StereoType::Bond_Double);
      CHECK(stereoInfo[1].specified == StereoSpecified::Specified);
    }
  }

  SECTION("More than one ring (1)") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      auto mol = R"SMI(CC\C=C/1C2CCC1[C@H]2O)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 4);
      CHECK(stereoInfo[0].centeredOn == 4);
      CHECK(stereoInfo[0].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == StereoSpecified::Unspecified);
      CHECK(stereoInfo[1].centeredOn == 7);
      CHECK(stereoInfo[1].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[1].specified == StereoSpecified::Unspecified);
      CHECK(stereoInfo[2].centeredOn == 8);
      CHECK(stereoInfo[2].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[2].specified == (use_legacy
                                            ? StereoSpecified::Unspecified
                                            : StereoSpecified::Specified));
      CHECK(stereoInfo[3].centeredOn == 2);
      CHECK(stereoInfo[3].type == StereoType::Bond_Double);
      CHECK(stereoInfo[3].specified == (use_legacy
                                            ? StereoSpecified::Unspecified
                                            : StereoSpecified::Specified));
    }
  }

  SECTION("More than one ring (2)") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      // This is the same structure as previos section, but
      // the rings are defined in the opposite order
      auto mol = R"SMI(CC\C=C/1C2[C@H](O)C1CC2)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      REQUIRE(stereoInfo.size() == 4);
      CHECK(stereoInfo[0].centeredOn == 4);
      CHECK(stereoInfo[0].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[0].specified == StereoSpecified::Unspecified);

      CHECK(stereoInfo[1].centeredOn == 5);
      CHECK(stereoInfo[1].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[1].specified == (use_legacy
                                            ? StereoSpecified::Unspecified
                                            : StereoSpecified::Specified));

      CHECK(stereoInfo[2].centeredOn == 7);
      CHECK(stereoInfo[2].type == StereoType::Atom_Tetrahedral);
      CHECK(stereoInfo[2].specified == StereoSpecified::Unspecified);

      CHECK(stereoInfo[3].centeredOn == 2);
      CHECK(stereoInfo[3].type == StereoType::Bond_Double);
      CHECK(stereoInfo[3].specified == (use_legacy
                                            ? StereoSpecified::Unspecified
                                            : StereoSpecified::Specified));
    }
  }

  SECTION("Stereo not possible: symmetric opposite atom") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      auto mol = R"SMI(CCC=C1CC(O)(O)C1)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      CHECK(stereoInfo.empty());
    }
  }

  SECTION("Stereo not possible: odd-sized ring") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);

      auto mol = R"SMI(CCC=C1CCCC1)SMI"_smiles;
      REQUIRE(mol);
      auto stereoInfo = Chirality::findPotentialStereo(*mol);
      CHECK(stereoInfo.empty());
    }
  }
}

TEST_CASE("double bond stereo with new chirality perception") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);
  SECTION("chain bonds") {
    {
      auto m = "C/C=C/C"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOTRANS);
      CHECK(m->getBondWithIdx(1)->getStereoAtoms() == std::vector<int>{0, 3});
    }
  }
  SECTION("ring bonds") {
    {
      auto m = "C1/C=C/CCCCCCC1"_smiles;
      REQUIRE(m);
      CHECK(m->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOTRANS);
      CHECK(m->getBondWithIdx(1)->getStereoAtoms() == std::vector<int>{0, 3});
    }
  }
}

TEST_CASE("false positives from new stereo code") {
  SECTION("elements") {
    std::vector<std::string> examples{"P", "PC", "S", "SC", "S(F)C"};
    for (auto &smi : examples) {
      INFO(smi);
      std::unique_ptr<RWMol> m{SmilesToMol(smi)};
      REQUIRE(m);
      auto si = Chirality::findPotentialStereo(*m);
      CHECK(si.empty());
    }
  }
  SECTION("non-tetrahedral and implicit Hs") {
    std::vector<std::string> examples{
        "[SiH4]",         "[SiH3]C",      "[SH4]",     "[PH5]",
        "[PH4]C",         "[SH6]",        "[SH5]C",  //"[SiH2](C)C",
        "[PH3](C)C",      "[PH2](C)(C)C", "[SH4](C)C", "[SH3](C)(C)C",
        "[SH2](C)(C)(C)C"};
    {
      AllowNontetrahedralChiralityFixture reset_non_tetrahedral_allowed;
      Chirality::setAllowNontetrahedralChirality(false);
      for (auto &smi : examples) {
        INFO(smi);
        std::unique_ptr<RWMol> m{SmilesToMol(smi)};
        REQUIRE(m);
        auto si = Chirality::findPotentialStereo(*m);
        CHECK(si.empty());
      }
    }
    {
      AllowNontetrahedralChiralityFixture reset_non_tetrahedral_allowed;
      Chirality::setAllowNontetrahedralChirality(true);
      for (auto &smi : examples) {
        INFO(smi);
        std::unique_ptr<RWMol> m{SmilesToMol(smi)};
        REQUIRE(m);
        auto si = Chirality::findPotentialStereo(*m);
        CHECK(si.size() == 1);
      }
    }
  }
}

TEST_CASE(
    "Github #6049: Cyclobutyl group in a macrocycle triggers a stereo center") {
  SECTION("as reported") {
    auto mol = "O=S1(=O)C=CC=C2CCCCC3CC(C3)N21"_smiles;
    REQUIRE(mol);
    auto stereoInfo = Chirality::findPotentialStereo(*mol);
    for (const auto &sg : stereoInfo) {
      CHECK(sg.centeredOn != 15);
    }
  }
}

void testStereoValidationFromMol(std::string molBlock,
                                 std::string expectedSmiles, bool legacyFlag,
                                 bool canonicalFlag = false) {
  RDKit::Chirality::setUseLegacyStereoPerception(legacyFlag);

  std::unique_ptr<RWMol> mol(MolBlockToMol(molBlock, true, false, false));
  REQUIRE(mol);

  // CHECK(CIPLabeler::validateStereochem(*mol, validationFlags));

  RDKit::SmilesWriteParams smilesWriteParams;
  smilesWriteParams.doIsomericSmiles = true;
  smilesWriteParams.doKekule = false;
  smilesWriteParams.canonical = canonicalFlag;

  unsigned int flags = 0 |
                       RDKit::SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES |
                       RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_PROPS |
                       RDKit::SmilesWrite::CXSmilesFields::CX_BOND_CFG |
                       RDKit::SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO |
                       RDKit::SmilesWrite::CXSmilesFields::CX_SGROUPS |
                       RDKit::SmilesWrite::CXSmilesFields::CX_POLYMER;

  auto outSmiles = MolToCXSmiles(*mol, smilesWriteParams, flags);
  RDKit::Chirality::setUseLegacyStereoPerception(false);

  CHECK(outSmiles == expectedSmiles);
}

std::string validateStereoMolBlockSpiro = R"(
  Mrv2308 06232316112D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 14 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.2921 0.9158 0 0
M  V30 2 C -4.2921 -0.6244 0 0
M  V30 3 C -2.9583 -1.3942 0 0 CFG=1
M  V30 4 C -1.6246 -0.6244 0 0
M  V30 5 C -1.6246 0.9158 0 0
M  V30 6 C -2.9583 4.7658 0 0 CFG=2
M  V30 7 C -4.2921 3.9958 0 0
M  V30 8 C -4.2921 2.4556 0 0
M  V30 9 C -2.9583 1.6858 0 0 CFG=2
M  V30 10 C -1.6246 2.4556 0 0
M  V30 11 C -1.6246 3.9958 0 0
M  V30 12 C -2.9583 6.3058 0 0
M  V30 13 Cl -2.9583 -2.9342 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 6 7
M  V30 6 1 6 11
M  V30 7 1 7 8
M  V30 8 1 9 8 CFG=1
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 9 1
M  V30 12 1 9 5
M  V30 13 1 6 12 CFG=1
M  V30 14 1 3 13 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END

  )";

std::string validateStereoMolBlockDoubleBondNoStereo = R"(
  Mrv2308 06232316392D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8498 2.3564 0 0
M  V30 2 C 2.1835 1.5864 0 0
M  V30 3 C 2.1835 0.0464 0 0
M  V30 4 C -0.4839 1.5864 0 0
M  V30 5 O -1.8176 2.3564 0 0
M  V30 6 C -1.8176 3.8964 0 0
M  V30 7 C -3.3342 3.629 0 0
M  V30 8 O -0.4839 4.6664 0 0
M  V30 9 C 0.8498 3.8964 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 1 1 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 6 8
M  V30 8 1 8 9
M  V30 9 1 1 9
M  V30 END BOND
M  V30 END CTAB
M  END
  )";

std::string validateStereoMolBlockDoubleBondNoStereo2 = R"(
  Mrv0541 07011416342D          

 21 22  0  0  0  0            999 V2000
   -1.9814    1.4834    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9581   -2.7980    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3658   -2.0787    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8189    1.5764    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5800    0.7796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9760    2.3967    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1333   -2.7836    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9581   -1.3592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7256   -2.0690    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4882   -1.3401    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8471    2.4140    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7304   -0.6351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9007   -0.6351    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4834    0.0700    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1333   -1.3497    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3658    0.9831    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -1.9814    0.0525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7598    0.7854    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3410    0.0700    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6556    2.2396    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2722    2.7980    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
 20  1  1  0  0  0  0
  4  1  2  0  0  0  0
  5  1  1  0  0  0  0
  2  7  2  0  0  0  0
  3  2  1  0  0  0  0
  3  8  2  0  0  0  0
  6  4  1  0  0  0  0
 16  4  1  0  0  0  0
 18  5  1  0  0  0  0
 17  5  2  0  0  0  0
 21  6  2  0  0  0  0
  7  9  1  0  0  0  0
  8 15  1  0  0  0  0
  9 15  2  0  0  0  0
 10 13  1  0  0  0  0
 11 20  1  0  0  0  0
 13 12  2  0  0  0  0
 15 12  1  0  0  0  0
 14 13  1  0  0  0  0
 19 14  2  0  0  0  0
 19 18  1  0  0  0  0
 21 20  1  0  0  0  0
M  END
> <Compound Name>
VM-0411129

> <CDD Number>
CDD-2839271

$$$$

)";

std::string validateStereoError1 = R"(
  -ISIS-  -- StrEd -- 

 29 32  0  0  0  0  0  0  0  0999 V2000
   -1.2050   -4.7172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8959   -3.7660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0823   -3.5582    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7514   -4.3013    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7296   -4.0934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0385   -3.1425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3694   -2.3993    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6785   -1.4481    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0094   -0.7050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0312   -0.9129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9470   -1.1209    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3183    0.2461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2694    0.5550    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2694    1.5549    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.3183    1.8641    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0094    2.8151    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.6785    3.5582    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3694    4.5093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3912    4.7172    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2777    3.9739    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0312    3.0230    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7305    1.0550    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2694    1.0550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7694    1.9210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7695    1.9210    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2694    1.0550    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7695    0.1889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7694    0.1889    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3912   -2.6071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  2  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  9 10  1  0  0  0  0
 10 11  3  0  0  0  0
  9 12  1  0  0  0  0
 12 13  2  0  0  0  0
 13 14  1  0  0  0  0
 14 15  2  0  0  0  0
 15 16  1  0  0  0  0
 16 17  1  0  0  0  0
 17 18  1  0  0  0  0
 18 19  1  0  0  0  0
 19 20  1  0  0  0  0
 20 21  1  0  0  0  0
 16 21  1  0  0  0  0
 15 22  1  0  0  0  0
 12 22  1  0  0  0  0
 22 23  1  0  0  0  0
 23 24  1  0  0  0  0
 24 25  2  0  0  0  0
 25 26  1  0  0  0  0
 26 27  2  0  0  0  0
 27 28  1  0  0  0  0
 23 28  2  0  0  0  0
  7 29  1  0  0  0  0
  3 29  2  0  0  0  0
M  END
> <Compound Name>
Z362114294

> <CDD Number>
CDD-3164311

$$$$

)";

std::string validateStereoUniq1 = R"(
  Mrv0541 06301412152D          

 15 15  0  0  0  0            999 V2000
    1.0464   -0.3197    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0464   -1.1460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7601   -1.5571    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.4739   -1.1460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4739   -0.3197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7601    0.0914    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3309   -1.5706    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1912    0.0976    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9103   -0.3114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1853    0.9197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.9115   -1.1336    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7601    0.9135    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    4.6215   -1.5447    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    1.7601   -2.3793    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    3.1894   -1.5665    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  1  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  7  2  0  0  0  0
  3  4  1  0  0  0  0
  3 14  1  0  0  0  0
  4 15  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  5  8  1  0  0  0  0
  6 12  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  2  0  0  0  0
  9 11  2  0  0  0  0
 11 13  1  0  0  0  0
M  STY  4   1 SUP   2 SUP   3 SUP   4 SUP
M  SAL   1  1  12
M  SBL   1  1  11
M  SMT   1 R3a
M  SAP   1  1  12   6  1
M  SCL   1 CXN
M  SAL   2  1  13
M  SBL   2  1  15
M  SMT   2 R3b
M  SAP   2  1  13  11  1
M  SCL   2 CXN
M  SAL   3  1  14
M  SBL   3  1   6
M  SMT   3 R1
M  SAP   3  1  14   3  1
M  SCL   3 CXN
M  SAL   4  1  15
M  SBL   4  1   7
M  SMT   4 R2
M  SAP   4  1  15   4  1
M  SCL   4 CXN
M  END
> <Compound Name>
VM-0021367

> <CDD Number>
CDD-3832787

$$$$
)";

TEST_CASE("ValidateStereo", "[accurateCIP]") {
  SECTION("validateStereoUniqOldCanon1") {
    testStereoValidationFromMol(validateStereoUniq1,
                                "*/C=C/C(=O)C1=C(*)N(*)C(=O)NC1* |,,,|", true,
                                true);
  }

  SECTION("validateStereoUniqNewCanon1") {
    testStereoValidationFromMol(validateStereoUniq1,
                                "*/C=C/C(=O)C1=C(*)N(*)C(=O)NC1* |,,,|", false,
                                true);
  }

  SECTION("validateStereoUniqNewNoCanon1") {
    testStereoValidationFromMol(validateStereoUniq1,
                                "N1C(=O)N(*)C(*)=C(C(/C=C/*)=O)C1* |,,,|",
                                false, false);
  }

  SECTION("SprioChiralLostOldNoCanon") {
    testStereoValidationFromMol(validateStereoMolBlockSpiro,
                                "C1C[C@H](Cl)CCC12CC[C@@H](C)CC2", true, false);
  }

  SECTION("SprioChiralLostOldCanon") {
    testStereoValidationFromMol(validateStereoMolBlockSpiro,
                                "C[C@H]1CCC2(CC1)CC[C@H](Cl)CC2", true, true);
  }

  SECTION("SprioChiralNotLostNewNoCanon") {
    testStereoValidationFromMol(validateStereoMolBlockSpiro,
                                "C1C[C@H](Cl)CC[C@]12CC[C@@H](C)CC2", false,
                                false);
  }

  SECTION("SprioChiralNotLostNewCanon") {
    testStereoValidationFromMol(validateStereoMolBlockSpiro,
                                "C[C@H]1CC[C@@]2(CC1)CC[C@H](Cl)CC2", false,
                                true);
  }

  SECTION("DoubleBondMarkedStereoOldNoCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo,
                                "C1(=CC)COC(C)OC1", true, false);
  }
  SECTION("DoubleBondMarkedStereoNewNoCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo,
                                "C1(=CC)COC(C)OC1", false, false);
  }

  SECTION("DoubleBondMarkedStereoOldCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo,
                                "CC=C1COC(C)OC1", true, true);
  }

  SECTION("DoubleBondMarkedStereoNewCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo,
                                "CC=C1COC(C)OC1", false, true);
  }

  SECTION("DoubleBondStereo2OldCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo2,
                                "CC(/C=N/NC(=O)c1c(Br)cnn1C)=C\\c1ccccc1", true,
                                true);
  }
  SECTION("DoubleBondStereo2NewCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo2,
                                "CC(=C\\c1ccccc1)/C=N/NC(=O)c1c(Br)cnn1C",
                                false, true);
  }

  SECTION("DoubleBondStereo2NewNoCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo2,
                                "c1(C(=O)N/N=C/C(C)=C/c2ccccc2)c(Br)cnn1C",
                                false, false);
  }
  SECTION("DoubleBondStereo2OldNoCanon") {
    testStereoValidationFromMol(validateStereoMolBlockDoubleBondNoStereo2,
                                "c1(C(=O)N/N=C/C(C)=C/c2ccccc2)c(Br)cnn1C",
                                true, false);
  }

  SECTION("ValidateStereoError1OldCanon") {
    testStereoValidationFromMol(
        validateStereoError1,
        "COc1cccc(/C=C(\\C#N)c2nnc(N3CCOCC3)n2-c2ccccc2)c1", true, true);
  }
  SECTION("ValidateStereoError1NewCanon") {
    testStereoValidationFromMol(
        validateStereoError1,
        "COc1cccc(/C=C(\\C#N)c2nnc(N3CCOCC3)n2-c2ccccc2)c1", false, true);
  }
  SECTION("ValidateStereoError1OldNoCanon") {
    testStereoValidationFromMol(
        validateStereoError1,
        "COc1cccc(/C=C(\\C#N)c2nnc(N3CCOCC3)n2-c2ccccc2)c1", true, false);
  }
  SECTION("ValidateStereoError1NewNoCanon") {
    testStereoValidationFromMol(
        validateStereoError1,
        "COc1cccc(/C=C(\\C#N)c2nnc(N3CCOCC3)n2-c2ccccc2)c1", false, false);
  }
}

TEST_CASE(
    "assignStereochemistry should clear crossed double bonds that can't have stereo") {
  SECTION("basics") {
    auto m = "CC=C(C)C"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->setBondDir(Bond::BondDir::EITHERDOUBLE);
    bool clean = true;
    bool flag = true;
    bool force = true;
    {
      UseLegacyStereoPerceptionFixture reset_stereo_perception;
      Chirality::setUseLegacyStereoPerception(false);
      auto cp(*m);
      RDKit::MolOps::assignStereochemistry(cp, clean, force, flag);
      CHECK(cp.getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
    }
    {
      UseLegacyStereoPerceptionFixture reset_stereo_perception;
      Chirality::setUseLegacyStereoPerception(true);
      auto cp(*m);
      RDKit::MolOps::assignStereochemistry(cp, clean, force, flag);
      CHECK(cp.getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
    }
  }
  SECTION("make sure we don't mess with actual potential stereosystems") {
    auto m = "CC=C(C)[13CH3]"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->setBondDir(Bond::BondDir::EITHERDOUBLE);
    bool clean = true;
    bool flag = true;
    bool force = true;
    {
      UseLegacyStereoPerceptionFixture reset_stereo_perception;
      Chirality::setUseLegacyStereoPerception(false);
      auto cp(*m);
      RDKit::MolOps::assignStereochemistry(cp, clean, force, flag);
      // the crossed bond dir has been translated to unknown stereo:
      CHECK(cp.getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOANY);
    }
    {
      UseLegacyStereoPerceptionFixture reset_stereo_perception;
      Chirality::setUseLegacyStereoPerception(true);
      auto cp(*m);
      RDKit::MolOps::assignStereochemistry(cp, clean, force, flag);
      // the crossed bond dir has been translated to unknown stereo:
      CHECK(cp.getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOANY);
    }
  }
  SECTION("make sure stereoatoms are also cleared") {
    auto m = "CC=C(C)C"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->setBondDir(Bond::BondDir::EITHERDOUBLE);
    m->getBondWithIdx(1)->setStereoAtoms(0, 3);
    bool clean = true;
    bool flag = true;
    bool force = true;
    {
      UseLegacyStereoPerceptionFixture reset_stereo_perception;
      Chirality::setUseLegacyStereoPerception(false);
      auto cp(*m);
      RDKit::MolOps::assignStereochemistry(cp, clean, force, flag);
      CHECK(cp.getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondWithIdx(1)->getStereoAtoms().empty());
    }
    {
      UseLegacyStereoPerceptionFixture reset_stereo_perception;
      Chirality::setUseLegacyStereoPerception(true);
      auto cp(*m);
      RDKit::MolOps::assignStereochemistry(cp, clean, force, flag);
      CHECK(cp.getBondWithIdx(1)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondWithIdx(1)->getStereoAtoms().empty());
    }
  }
  SECTION("ensure we can enumerate stereo on either double bonds") {
    auto mol = R"CTAB(
  Mrv2004 11072316002D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.859 2.7821 0 0
M  V30 2 C 3.2818 2.1928 0 0
M  V30 3 C 3.8711 0.77 0 0
M  V30 4 C 3.2818 -0.6528 0 0
M  V30 5 C 1.859 -1.2421 0 0
M  V30 6 C 0.4362 -0.6528 0 0
M  V30 7 C -0.1533 0.77 0 0
M  V30 8 C 0.4362 2.1928 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 3 CFG=2
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 1 8
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    // mol->debugMol(std::cerr);
    std::string smi = MolToCXSmiles(*mol, SmilesWriteParams());
    std::unique_ptr<ROMol> f(SmilesToMol(smi));
    mol->getBondWithIdx(3)->setStereo(Bond::BondStereo::STEREOCIS);
    f->getBondWithIdx(0)->setStereo(Bond::BondStereo::STEREOCIS);
    CHECK(MolToSmiles(*mol) == "C1=C\\CCCCCC/1");
    CHECK(MolToSmiles(*f) == "C1=C\\CCCCCC/1");
    mol->getBondWithIdx(3)->setStereo(Bond::BondStereo::STEREOTRANS);
    f->getBondWithIdx(0)->setStereo(Bond::BondStereo::STEREOTRANS);
    CHECK(MolToSmiles(*mol) == "C1=C/CCCCCC/1");
    CHECK(MolToSmiles(*f) == "C1=C/CCCCCC/1");
  }
}
TEST_CASE("adding two wedges to chiral centers") {
  SECTION("basics") {
    auto mol = R"CTAB(
  Mrv2219 02112315062D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -3.6667 2.5 0 0
M  V30 2 C -2.333 3.27 0 0 CFG=1
M  V30 3 F -0.9993 2.5 0 0
M  V30 4 C -3.103 4.6037 0 0
M  V30 5 N -1.3955 4.4918 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 5);

    Chirality::BondWedgingParameters ps;
    ps.wedgeTwoBondsIfPossible = true;

    {
      RWMol cp(*mol);
      Chirality::wedgeMolBonds(cp);
      CHECK(cp.getBondBetweenAtoms(1, 3)->getBondDir() != Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 4)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 0)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 2)->getBondDir() == Bond::BondDir::NONE);
    }
    {
      RWMol cp(*mol);
      Chirality::wedgeMolBonds(cp, nullptr, &ps);
      CHECK(cp.getBondBetweenAtoms(1, 3)->getBondDir() != Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 4)->getBondDir() != Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 4)->getBondDir() !=
            cp.getBondBetweenAtoms(1, 3)->getBondDir());
      CHECK(cp.getBondBetweenAtoms(1, 0)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 2)->getBondDir() == Bond::BondDir::NONE);
    }
    {
      // Wedge a second bond after we already wedged a first one

      RWMol cp(*mol);
      Chirality::wedgeMolBonds(cp);
      CHECK(cp.getBondBetweenAtoms(1, 3)->getBondDir() != Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 4)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 0)->getBondDir() == Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 2)->getBondDir() == Bond::BondDir::NONE);

      Chirality::wedgeMolBonds(cp, nullptr, &ps);

      REQUIRE(count_wedged_bonds(cp) == 2);

      CHECK(cp.getBondBetweenAtoms(1, 3)->getBondDir() != Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 4)->getBondDir() != Bond::BondDir::NONE);
      CHECK(cp.getBondBetweenAtoms(1, 4)->getBondDir() !=
            cp.getBondBetweenAtoms(1, 3)->getBondDir());
    }
    {
      // What if the first wedged bond is not our preferred one?

      for (auto wedgedAtomIdx : {0, 2, 4}) {
        INFO("wedgedAtomIdx: " << wedgedAtomIdx);

        RWMol cp(*mol);
        REQUIRE(count_wedged_bonds(cp) == 0);

        auto bond = cp.getBondBetweenAtoms(1, wedgedAtomIdx);
        if (bond->getEndAtomIdx() == 1) {
          // One of the bonds in the input is reversed, make sure the chiral
          // atom is always at the start so that the wedge is valid!
          auto tmp = bond->getBeginAtomIdx();
          bond->setBeginAtomIdx(bond->getEndAtomIdx());
          bond->setEndAtomIdx(tmp);
        }

        // This probably disagrees with the chirality in some of the test cases,
        // but that's not relevant for this test
        bond->setBondDir(Bond::BondDir::BEGINWEDGE);

        Chirality::wedgeMolBonds(cp, nullptr, &ps);

        CHECK(count_wedged_bonds(cp) == 2);
      }
    }
  }
  SECTION(
      "more complex 1, this should only have one wedge for each of the two chiral centers") {
    std::string smi =
        "[H][C@@]12CC(=O)N1[C@@H](C(=O)O)C(C)(C)S2(=O)=O |(-2.78382,0.183015,;-1.38222,-0.351313,;-2.12923,-1.65207,;-0.828466,-2.39908,;-0.436905,-3.84707,;-0.0814577,-1.09832,;1.03095,-0.0920638,;2.49888,-0.400554,;2.96569,-1.82607,;3.50001,0.71647,;0.41769,1.27685,;1.8432,1.74365,;0.102447,2.74335,;-1.07373,1.11662,;-1.07718,2.61662,;-2.56587,1.26998,)|";
    SmilesParserParams spps;
    spps.removeHs = false;
    std::unique_ptr<RWMol> m{SmilesToMol(smi, spps)};
    REQUIRE(m);
    Chirality::BondWedgingParameters bwps;
    bwps.wedgeTwoBondsIfPossible = true;
    Chirality::wedgeMolBonds(*m, &m->getConformer(), &bwps);
    unsigned nWedged = 0;
    for (const auto bond : m->bonds()) {
      if (bond->getBondDir() != Bond::BondDir::NONE) {
        ++nWedged;
      }
    }
    CHECK(nWedged == 2);
  }
  SECTION("more complex 2, have two wedges around the chiral center") {
    std::string smi =
        "[H][C@@]12CCCN1C(=O)CN1C(=O)[C@](C)(N)O[C@]12O |(-0.888297,0.626611,;-1.19852,-0.840959,;-1.94707,-2.14084,;-3.41464,-1.83061,;-3.5731,-0.339006,;-2.20347,0.272634,;-1.74154,1.69974,;-2.74648,2.81333,;-0.274666,2.01325,;0.730277,0.899655,;2.23028,0.901335,;3.11059,2.11585,;2.6954,-0.52473,;3.44685,-1.82293,;4.06503,0.0869091,;1.48286,-1.40777,;0.26835,-0.527448,;-0.0418744,-1.99502,)|";
    SmilesParserParams spps;
    spps.removeHs = false;
    std::unique_ptr<RWMol> m{SmilesToMol(smi, spps)};
    REQUIRE(m);
    Chirality::BondWedgingParameters bwps;
    bwps.wedgeTwoBondsIfPossible = true;
    Chirality::wedgeMolBonds(*m, &m->getConformer(), &bwps);
    CHECK(m->getBondWithIdx(12)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondWithIdx(13)->getBondDir() != Bond::BondDir::NONE);
  }
  SECTION("another one") {
    auto m =
        "CC[C@@]1(O)C(=O)OCc2c1cc1-c3nc4ccccc4cc3Cn1c2=O |(-2.67178,3.55256,;-2.43493,2.07138,;-3.59925,1.12567,;-4.32841,2.43652,;-5.01681,0.635215,;-6.15033,1.61763,;-5.30084,-0.837648,;-4.16731,-1.82006,;-2.74976,-1.3296,;-2.46573,0.143259,;-1.04818,0.633713,;0.085343,-0.348697,;1.57389,-0.163687,;2.43243,1.06631,;3.92692,0.937796,;4.78546,2.1678,;6.27994,2.03928,;6.91588,0.680758,;6.05734,-0.549244,;4.56286,-0.420725,;3.70432,-1.65073,;2.20983,-1.52221,;1.11432,-2.54683,;-0.198688,-1.82156,;-1.61624,-2.31201,;-1.90027,-3.78488,)|"_smiles;
    REQUIRE(m);
    Chirality::BondWedgingParameters bwps;
    bwps.wedgeTwoBondsIfPossible = true;
    Chirality::wedgeMolBonds(*m, &m->getConformer(), &bwps);
    CHECK(m->getBondWithIdx(1)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondWithIdx(1)->getBeginAtomIdx() == 2);
    CHECK(m->getBondWithIdx(2)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondWithIdx(2)->getBeginAtomIdx() == 2);
  }
  SECTION("favor degree 1") {
    auto m =
        "[H][C@@]12CC[C@@](C)(O)[C@H](CC[C@@](C)(O)C=C)[C@@]1(C)CCCC2(C)C |(3.59567,-1.0058,;2.33379,-0.194852,;2.4456,-1.69068,;1.20608,-2.53542,;-0.145252,-1.88434,;-1.63777,-1.73471,;-0.551787,-3.3282,;-0.257063,-0.388514,;-1.60839,0.262569,;-2.84791,-0.582176,;-4.19924,0.068907,;-5.55057,0.71999,;-4.85032,-1.28242,;-3.54816,1.42024,;-4.3929,2.65976,;0.982456,0.456231,;-0.368873,1.10731,;0.870645,1.95206,;2.11016,2.7968,;3.46149,2.14572,;3.5733,0.649893,;5.02699,1.01975,;4.35205,-0.632117,)|"_smiles;
    REQUIRE(m);
    Chirality::BondWedgingParameters bwps;
    bwps.wedgeTwoBondsIfPossible = true;
    Chirality::wedgeMolBonds(*m, &m->getConformer(), &bwps);
    CHECK(m->getBondWithIdx(9)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondWithIdx(10)->getBondDir() != Bond::BondDir::NONE);
    CHECK(m->getBondWithIdx(11)->getBondDir() == Bond::BondDir::NONE);
  }
}

TEST_CASE(
    "RDKit Issue #6217: Atoms may get flagged with non-tetrahedral stereo even when it is not allowed",
    "[bug][stereo][non-tetrahedral]") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);

  AllowNontetrahedralChiralityFixture reset_non_tetrahedral_allowed;
  Chirality::setAllowNontetrahedralChirality(false);

  auto m = "CS(=O)(=O)O"_smiles;
  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);

  auto stereoInfo = Chirality::findPotentialStereo(*m);
  CHECK(stereoInfo.size() == 0);

  auto at = m->getAtomWithIdx(1);

  auto sinfo = Chirality::detail::getStereoInfo(at);
  CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);

  REQUIRE(at->getAtomicNum() == 16);
  CHECK(!at->hasProp(common_properties::_ChiralityPossible));
}

TEST_CASE(
    "RDKit Issue #6239: Tri-coordinate atom with implicit + neighbor H atom is found potentially chiral",
    "[bug][stereo]") {
  // Parametrize test to run under legacy and new stereo perception
  const auto legacy_stereo = GENERATE(true, false);
  INFO("Legacy stereo perception == " << legacy_stereo);

  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(legacy_stereo);

  auto p = SmilesParserParams();
  p.removeHs = false;

  std::unique_ptr<RWMol> m{SmilesToMol("[H]C(C)CC", p)};

  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 5);

  auto at = m->getAtomWithIdx(1);
  REQUIRE(at->getAtomicNum() == 6);
  REQUIRE(at->getDegree() == 3);

  CHECK(!at->hasProp(common_properties::_ChiralityPossible));

  CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(at));

  auto stereoInfo = Chirality::findPotentialStereo(*m);
  CHECK(stereoInfo.size() == 0);

  auto sinfo = Chirality::detail::getStereoInfo(at);
  CHECK(sinfo.type == Chirality::StereoType::Atom_Tetrahedral);
}

TEST_CASE("double bonded N with H should be stereogenic", "[bug][stereo]") {
  SECTION("assign stereo") {
    // Parametrize test to run under legacy and new stereo perception
    const auto legacy_stereo = GENERATE(true, false);
    INFO("Legacy stereo perception == " << legacy_stereo);

    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(legacy_stereo);
    auto m = "[H]/N=C/F"_smiles;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(1)->getStereo() != Bond::BondStereo::STEREONONE);
  }
  SECTION("find potential stereo") {
    auto m = "[H]/N=C/F"_smiles;
    REQUIRE(m);
    CHECK(Chirality::detail::isBondPotentialStereoBond(m->getBondWithIdx(1)));
    bool cleanIt = false;
    bool flagPossible = true;
    auto si = Chirality::findPotentialStereo(*m, cleanIt, flagPossible);
    CHECK(si.size() == 1);
  }
}

TEST_CASE("Issue in GitHub #6473", "[bug][stereo]") {
  constexpr const char *mb = R"CTAB(
     RDKit          2D

  6  5  0  0  0  0  0  0  0  0999 V2000
    2.0443    0.2759    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7038    2.5963    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.3828    2.5961    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.0444    1.8228    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6359    1.8229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3827   -0.4985    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0
  4  2  1  0
  4  3  2  0
  2  5  1  0
  6  1  1  0
M  END)CTAB";

  UseLegacyStereoPerceptionFixture reset_stereo_perception;

  auto use_legacy_stereo = GENERATE(true, false);
  CAPTURE(use_legacy_stereo);
  Chirality::setUseLegacyStereoPerception(use_legacy_stereo);

  std::unique_ptr<ROMol> mol(MolBlockToMol(mb));
  REQUIRE(mol);

  // Removal of the stereogenic H will cause loss of stereo information
  // on the imine double bond, and although the bond is still detected as
  // potentially stereo, it will be reverted to "unspecified"
  auto bond = mol->getBondWithIdx(2);
  REQUIRE(bond->getBondType() == Bond::BondType::DOUBLE);
  CHECK(Chirality::detail::isBondPotentialStereoBond(bond));
  CHECK(bond->getStereo() == Bond::BondStereo::STEREONONE);

  if (!use_legacy_stereo) {
    auto sinfo = Chirality::detail::getStereoInfo(bond);
    REQUIRE(sinfo.type == Chirality::StereoType::Bond_Double);
    CHECK(sinfo.specified == Chirality::StereoSpecified::Unspecified);
    REQUIRE(sinfo.controllingAtoms.size() == 4);
    CHECK(sinfo.controllingAtoms[0] == 0);
    CHECK(sinfo.controllingAtoms[1] == 1);
    CHECK(sinfo.controllingAtoms[2] == Chirality::StereoInfo::NOATOM);
    CHECK(sinfo.controllingAtoms[3] == Chirality::StereoInfo::NOATOM);
  }
}

TEST_CASE("GitHub Issue #6640", "[bug][stereo]") {
  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(false);

  auto p = SmilesParserParams();
  p.sanitize = false;
  p.removeHs = false;
  std::string smiles{"NC1=NC(=N)N=C(N)C1C"};
  std::unique_ptr<RWMol> mol(SmilesToMol(smiles, p));
  REQUIRE(mol);

  MolOps::removeHs(*mol);

  auto cleanIt = true;
  auto force = true;
  auto flagPossibleStereoCenters = true;
  MolOps::assignStereochemistry(*mol, cleanIt, force,
                                flagPossibleStereoCenters);
}

TEST_CASE("zero bond-length chirality cases") {
  SECTION("basics") {
    {
      auto m = R"CTAB(derived from CHEMBL3183068
  Mrv2211 07202306222D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 10 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N 3.231 0 0 0
M  V30 2 C 4.0137 -1.3378 0 0
M  V30 3 C 2.6757 -2.1004 0 0
M  V30 4 C 2.6757 -3.6389 0 0
M  V30 5 C 1.3378 -1.3378 0 0
M  V30 6 C 5.4672 -0.99 0 0
M  V30 7 C 1.3378 -4.4216 0 0
M  V30 8 C 0 -2.1205 0 0
M  V30 9 C 0 -3.659 0 0
M  V30 10 F 4.0137 -1.3378 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 2 6
M  V30 4 1 2 10
M  V30 5 2 3 4
M  V30 6 1 3 5
M  V30 7 1 4 7
M  V30 8 2 5 8
M  V30 9 2 7 9
M  V30 10 1 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
    {
      auto m = R"CTAB(derived from CHEMBL3183068
  Mrv2211 07202306222D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 10 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N 3.231 0 0 0
M  V30 2 C 4.0137 -1.3378 0 0
M  V30 3 C 2.6757 -2.1004 0 0
M  V30 4 C 2.6757 -3.6389 0 0
M  V30 5 C 1.3378 -1.3378 0 0
M  V30 6 C 4.0137 -1.3378 0 0
M  V30 7 C 1.3378 -4.4216 0 0
M  V30 8 C 0 -2.1205 0 0
M  V30 9 C 0 -3.659 0 0
M  V30 10 F 4.5135 -2.6451 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 2 6
M  V30 4 1 2 10
M  V30 5 2 3 4
M  V30 6 1 3 5
M  V30 7 1 4 7
M  V30 8 2 5 8
M  V30 9 2 7 9
M  V30 10 1 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
    }
  }
}

TEST_CASE("t-shaped chirality cases") {
  SECTION("ChEMBL example") {
    {
      auto m = R"CTAB(CHEMBL3183068
     RDKit          2D

 11 11  0  0  1  0  0  0  0  0999 V2000
    1.7309    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1502   -0.7167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4334   -1.1252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4334   -1.9494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7167   -0.7167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8669   -0.2974    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7167   -2.3687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.1360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.9602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5836   -0.7060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5694   -1.4334    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  2  6  1  0
  2 11  1  0
  3  4  2  0
  3  5  1  0
  4  7  1  0
  5  8  2  0
  6 10  1  0
  7  9  2  0
  8  9  1  0
M  END
> <chembl_id>
CHEMBL3183068

> <chembl_pref_name>
None
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
    {
      auto m = R"CTAB(CHEMBL3183068 (with an H removed)
     RDKit          2D

 10 10  0  0  1  0  0  0  0  0999 V2000
    1.7309    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1502   -0.7167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4334   -1.1252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4334   -1.9494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7167   -0.7167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8669   -0.2974    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7167   -2.3687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.1360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.9602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5836   -0.7060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  2  6  1  0
  3  4  2  0
  3  5  1  0
  4  7  1  0
  5  8  2  0
  6 10  1  0
  7  9  2  0
  8  9  1  0
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
  }

  SECTION("four-coordinate") {
    {
      auto m = R"CTAB(
  Mrv2211 07202306492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N -3.3332 0.9919 0 0
M  V30 2 C -2.5555 -0.3373 0 0
M  V30 3 O -3.885 -1.095 0 0
M  V30 4 C -1.2263 0.4404 0 0
M  V30 5 F -1.7854 -1.6709 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
    {
      auto m = R"CTAB(
  Mrv2211 07202306492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N -3.3332 0.9919 0 0
M  V30 2 C -2.5555 -0.3373 0 0
M  V30 3 O -3.885 -1.095 0 0
M  V30 4 C -1.2263 0.4404 0 0
M  V30 5 F -1.7854 -1.6709 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
    {
      auto m = R"CTAB(
  Mrv2211 07202306492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 1
M  V30 BEGIN ATOM
M  V30 1 N -3.3332 0.9919 0 0
M  V30 2 C -2.5555 -0.3373 0 0
M  V30 3 O -3.885 -1.095 0 0
M  V30 4 C -1.2263 0.4404 0 0
M  V30 5 F -1.7854 -1.6709 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    }
  }
}

TEST_CASE("almost linear, degree 4") {
  SECTION("from chembl 1") {
    auto m = R"CTAB(CHEMBL3680147
     RDKit          2D

 15 16  0  0  1  0  0  0  0  0999 V2000
    3.8912   -4.9570    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8933   -3.7570    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5951   -3.0039    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5973   -1.5031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2073   -4.4014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2131   -3.2886    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.4656   -1.9881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9978   -2.2972    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.9557   -0.8928    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
 10  5  1  0
  2 11  1  0
 11 12  1  0
 12 13  1  0
 13 14  2  0
 14  2  1  0
 13 15  1  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
  SECTION("from chembl 2") {
    auto m = R"CTAB(CHEMBL76346
     RDKit          2D

 16 17  0  0  1  0  0  0  0  0999 V2000
    3.4292   -0.4250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.7167   -0.0125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7167    0.8125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.4292   -1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1500   -0.0125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7667   -1.7292    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1500    0.8208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4292    1.2333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0250   -2.5167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8500   -2.5167    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0042   -0.4250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1042   -1.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4250    2.0583    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2250   -2.7292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2375   -3.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0042   -3.5250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  3  2  1  0
  4  1  1  0
  5  1  1  0
  6  4  1  0
  7  5  2  0
  8  7  1  0
  9  6  1  0
 10 12  1  0
 11  2  2  0
 12  4  1  0
 13  8  1  0
  9 14  1  6
 15  9  1  0
 16 14  1  0
  3  8  2  0
 10  9  1  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(8)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
  SECTION("from chembl 3") {
    auto m = R"CTAB(CHEMBL3577363
     RDKit          2D

 15 16  0  0  0  0  0  0  0  0999 V2000
   -2.5989    1.5003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3863    2.3198    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8514    3.7459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3514    3.7442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8133    2.3171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3383   -1.3500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6369    0.8981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1471    4.7175    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  1  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  2  0
 11  6  1  0
  8 12  1  0
  9  1  1  0
  6 13  1  0
  1 14  1  6
  3 15  2  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
  SECTION("from chembl 4") {
    auto m =
        R"CTAB(derived from CHEMBL2373651. This was wrong in the RDKit implementation
  Mrv2211 07212313282D
            
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 1
M  V30 BEGIN ATOM
M  V30 1 O -7.4486 0.4751 0 0
M  V30 2 C -6.1148 -0.2949 0 0
M  V30 3 C -6.1148 1.2451 0 0
M  V30 4 C -4.7811 -1.0649 0 0
M  V30 5 C -4.7811 2.0151 0 0
M  V30 6 H -6.1148 1.2451 0 0
M  V30 7 H -6.1148 -0.2949 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 3 1
M  V30 3 1 2 3
M  V30 4 1 2 4
M  V30 5 1 3 5
M  V30 6 1 3 6 CFG=1
M  V30 7 1 2 7 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }
  SECTION("cut from CHEMBL4578507") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 10.694361 -35.424753 0.000000 0
M  V30 2 C 9.463880 -34.700253 0.000000 0
M  V30 3 C 8.217865 -35.405903 0.000000 0
M  V30 4 C 8.207044 -36.840940 0.000000 0
M  V30 5 N 9.449393 -37.570327 0.000000 0
M  V30 6 C 6.965044 -36.121153 0.000000 0
M  V30 7 O 8.215247 -33.985178 0.000000 0
M  V30 8 H 6.885805 -37.382701 0.000000 0
M  V30 9 H 9.169786 -33.305883 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 3
M  V30 4 1 4 5
M  V30 5 1 4 6
M  V30 6 1 2 7
M  V30 7 1 4 8 CFG=3
M  V30 8 1 2 9 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(3)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }

  SECTION("derived from CHEMBL4578507") {
    auto m = R"CTAB(CHEMBL4578507
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 17 19 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 6.828300 -21.526500 0.000000 0
M  V30 2 C 7.542200 -21.113800 0.000000 0
M  V30 3 C 7.537900 -20.291100 0.000000 0
M  V30 4 C 6.829700 -19.884700 0.000000 0
M  V30 5 C 6.121000 -21.114600 0.000000 0
M  V30 6 C 6.127300 -20.296500 0.000000 0
M  V30 7 C 5.422300 -19.881400 0.000000 0
M  V30 8 C 4.708400 -20.285700 0.000000 0
M  V30 9 C 4.702200 -21.107900 0.000000 0
M  V30 10 N 5.414000 -21.525800 0.000000 0
M  V30 11 C 3.990600 -20.695500 0.000000 0
M  V30 12 O 4.706900 -19.471700 0.000000 0
M  V30 13 C 3.988500 -19.876800 0.000000 0
M  V30 14 O 3.284900 -21.104500 0.000000 0
M  V30 15 F 8.242300 -19.879800 0.000000 0
M  V30 16 H 3.945200 -21.418300 0.000000 0 MASS=2
M  V30 17 H 5.253800 -19.082500 0.000000 0 MASS=2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 5 1
M  V30 2 1 1 2
M  V30 3 2 2 3
M  V30 4 1 3 4
M  V30 5 2 4 6
M  V30 6 1 5 6
M  V30 7 1 5 10
M  V30 8 1 6 7
M  V30 9 1 7 8
M  V30 10 1 9 8
M  V30 11 1 9 10
M  V30 12 1 9 11
M  V30 13 1 7 12
M  V30 14 1 11 13
M  V30 15 1 13 12
M  V30 16 1 11 14 CFG=3
M  V30 17 1 3 15
M  V30 18 1 9 16 CFG=3
M  V30 19 1 7 17 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(6)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(8)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(10)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }

  SECTION("derived from CHEMBL123021") {
    auto m = R"CTAB(CHEMBL123021
     RDKit          2D

 16 17  0  0  1  0  0  0  0  0999 V2000
   -2.8208   -0.5792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8208    0.2458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5375   -0.1542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4083   -0.5792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1125   -0.9917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4000    0.2458    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1125    0.6583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6833    0.6583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1083   -1.8167    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1125    1.4833    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6833    1.4833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6958   -0.9917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0250   -0.5875    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0375    1.8958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3161    1.1145    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3141    0.2829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  1  3  1  0
  4  5  1  0
  5  1  1  0
  6  4  2  0
  7  2  1  0
  8  6  1  0
  9  5  2  0
 10  7  2  0
 11  8  2  0
 12  4  1  0
 13 12  1  0
 14 11  1  0
  2  3  1  0
  7  6  1  0
  2 15  1  6
  1 16  1  6
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }

  SECTION("derived from CHEMBL85809") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 7.185783 -5.677981 0.000000 0
M  V30 2 C 7.207427 -7.106480 0.000000 0
M  V30 3 O 9.032790 -6.219079 0.000000 0
M  V30 4 O 10.172645 -8.527707 0.000000 0
M  V30 5 C 10.886895 -7.106480 0.000000 0
M  V30 6 C 7.900033 -8.527707 0.000000 0
M  V30 7 C 10.865251 -5.677981 0.000000 0
M  V30 8 H 7.185783 -8.534979 0.000000 0 MASS=2
M  V30 9 H 10.865251 -8.534979 0.000000 0 MASS=2
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=1
M  V30 2 1 3 2
M  V30 3 1 4 6
M  V30 4 1 5 3 
M  V30 5 1 6 2
M  V30 6 1 5 7 CFG=1
M  V30 7 1 2 8 CFG=3
M  V30 8 1 5 9 CFG=3
M  V30 9 1 4 5
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(2 2 5)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    CHECK(m->getAtomWithIdx(4)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }

  SECTION("derived from CHEMBL2333552") {
    auto m = R"CTAB(blah
     RDKit          2D

 18 20  0  0  0  0  0  0  0  0999 V2000
   35.6738   -9.2984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   34.9598   -8.8850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   34.9588   -9.7100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   32.1164   -7.2274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   32.1164   -8.0541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   32.8299   -8.4634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   32.8299   -6.8099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   33.5434   -7.2274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   33.5399   -8.0541    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   34.9720   -7.2335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   34.2572   -6.8149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   31.3992   -6.8162    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   32.8308   -9.2901    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   33.5319   -8.8767    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   34.9684   -8.0602    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   34.2507   -8.4688    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   35.7656   -8.2712    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   34.2458   -7.6408    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  3  2  1  0
  4  5  1  0
  4  7  1  0
  5  6  1  0
  6  9  1  0
  8  7  2  0
  8  9  1  0
  8 11  1  0
  9 16  1  0
 15 10  1  0
 10 11  1  0
  4 12  2  0
  6 13  1  6
  9 14  1  6
 16 15  1  0
  2 16  1  0
 15  2  1  0
 15 17  1  6
 16 18  1  6
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(5)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(m->getAtomWithIdx(8)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    CHECK(m->getAtomWithIdx(14)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(m->getAtomWithIdx(15)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }
  SECTION("CHEMBL94022") {
    auto m = R"CTAB(CHEMBL94022
     RDKit          2D

 16 18  0  0  1  0  0  0  0  0999 V2000
    2.0917   -2.6875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2667   -2.7042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6917   -3.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6750   -3.2667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4792   -4.0667    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.5000   -3.2042    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5500   -2.2917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1792   -4.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8125   -3.9667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1708   -2.7042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5500   -1.4667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8833   -2.2917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1708   -1.0542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8833   -1.4667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6792   -3.2917    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.7917   -2.2542    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  3  1  1  0
  1  4  1  1
  5  4  2  0
  6  4  1  0
  7  2  1  0
  8  5  1  0
  9  6  1  0
 10  7  2  0
 11  7  1  0
 12 10  1  0
 13 11  2  0
 14 13  1  0
  2 15  1  1
  1 16  1  6
  3  2  1  0
  8  9  1  0
 14 12  2  0
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  }
  SECTION("overlapping neighbors") {
    auto m = R"CTAB(derived from CHEMBL3752539
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.237143 -0.714286 0.000000 0
M  V30 2 C -1.237143 0.714286 0.000000 0
M  V30 3 C -2.474381 1.428571 0.000000 0
M  V30 4 C -3.711524 0.714286 0.000000 0
M  V30 5 C -2.474381 -1.428571 0.000000 0
M  V30 6 C -3.711524 -0.714286 0.000000 0
M  V30 7 O -2.474381 -2.857143 0.000000 0
M  V30 8 H -4.536286 -0.238095 0.000000 0
M  V30 9 H -1.649524 -1.904762 0.000000 0
M  V30 10 C -4.948762 -1.428571 0.000000 0
M  V30 11 C -2.474381 -1.428571 0.000000 0
M  V30 12 C -3.711524 0.714286 0.000000 0
M  V30 13 O -3.711524 -0.714286 0.000000 0
M  V30 14 H -2.474381 2.380952 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 5
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 5 1 4 6
M  V30 6 1 5 6
M  V30 7 1 5 7
M  V30 8 1 6 10
M  V30 9 1 6 8 CFG=1
M  V30 10 1 5 9 CFG=3
M  V30 11 1 1 11 CFG=3
M  V30 12 1 3 12
M  V30 13 1 12 13
M  V30 14 1 11 13
M  V30 15 1 3 14 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 6)
M  V30 MDLV30/STERAC1 ATOMS=(3 1 3 5)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getAtomWithIdx(2)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("bond atoms overlapping central atom at the end of wedge bonds") {
    auto m = R"CTAB(CHEMBL3612237
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
   -0.6828   -1.6239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6828    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6090    0.5905    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9007    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9007   -1.6239    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6090   -2.3805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9007    1.3471    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6828    1.3471    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6090    2.0668    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6319    1.4849    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5072    2.6784    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7280    0.9964    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6828    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9007    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  6  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  4  7  1  0
  2  8  1  0
  8  9  1  0
  7  9  1  0
  3 10  1  1
 10 11  1  0
 10 12  2  0
  2 13  1  1
  4 14  1  1
M  END)CTAB"_ctab;
    REQUIRE(m);
    // if the bond is wedged, then we should have chirality even if the bonded
    // atom overlaps the central atom
    CHECK(m->getAtomWithIdx(1)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(3)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
    CHECK(m->getAtomWithIdx(2)->getChiralTag() !=
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
}

TEST_CASE("github #6931: atom maps influencing chirality perception") {
  SECTION("basics") {
    auto m = "[CH3:1]C([CH3:2])(O)F"_smiles;
    REQUIRE(m);
    bool cleanIt = true;
    bool force = true;
    bool flagPossibleStereoCenters = true;
    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(false);
    MolOps::assignStereochemistry(*m, cleanIt, force,
                                  flagPossibleStereoCenters);
    CHECK(
        !m->getAtomWithIdx(1)->hasProp(common_properties::_ChiralityPossible));
  }
}

TEST_CASE(
    "Github Issue #6981: Parsing a Mol leaks the \"_needsDetectBondStereo\" property",
    "[bug][stereo]") {
  // Parametrize test to run under legacy and new stereo perception
  const auto legacy_stereo = GENERATE(true, false);
  INFO("Legacy stereo perception == " << legacy_stereo);

  UseLegacyStereoPerceptionFixture reset_stereo_perception;
  Chirality::setUseLegacyStereoPerception(legacy_stereo);

  auto m = R"CTAB(
  Mrv2311 12122315472D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -9.2083 1.8333 0 0
M  V30 2 C -8.0639 0.8029 0 0
M  V30 3 C -6.5239 0.8029 0 0
M  V30 4 C -5.7539 -0.5308 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;

  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 4);

  CHECK(m->hasProp("_needsDetectBondStereo") == false);
}

TEST_CASE(
    "Github Issue #7076: new stereo code not properly handling crossed double bonds") {
  // Parametrize test to run under legacy and new stereo perception
  SECTION("second part") {
    const auto legacy_stereo = GENERATE(true, false);
    INFO("Legacy stereo perception == " << legacy_stereo);

    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(legacy_stereo);

    auto m = R"CTAB(
  Mrv2211 01252410552D          

 10  9  0  0  0  0            999 V2000
    0.0000   -1.4364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7108   -3.4884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -3.8988    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8208   -1.4364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.2572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7108   -2.6676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4217   -2.2572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -4.7196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2439   -5.4378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7108   -5.1300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0  0  0  0
  1  5  2  0  0  0  0
  2  3  1  0  0  0  0
  2  6  2  3  0  0  0
  3  8  2  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
M  END)CTAB"_ctab;

    REQUIRE(m);
    CHECK(m->getBondWithIdx(3)->getStereo() == Bond::BondStereo::STEREOANY);
  }
  SECTION("original") {
    std::string ctab = R"CTAB(7630532
     RDKit          2D

 37 41  0  0  0  0  0  0  0  0999 V2000
    0.0000   -1.7500    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -4.2500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -4.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7500   -1.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000   -1.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2475   -0.8825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4975   -2.6175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2475   -2.6175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4975   -0.8825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8660   -3.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2527   -0.8825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5027   -2.6175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2527   -2.6175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5027   -0.8825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    2.0104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7604   -1.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0104   -1.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7321   -2.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -5.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5155   -6.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8660   -6.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3801   -6.1225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8631   -7.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5126   -7.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7306   -5.7475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2507   -6.6251    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7337   -7.7526    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3832   -8.1276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6012   -6.2501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2566   -7.6302    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6071   -7.2552    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  2  0
  2  3  1  0
  2 14  2  3
  3 25  2  0
  4  8  2  0
  4 11  1  0
  5  9  2  0
  5 12  1  0
  6 10  2  0
  6 13  1  0
  7 14  1  0
  8 15  1  0
  9 16  1  0
 10 17  1  0
 11 18  2  0
 12 19  2  0
 13 20  2  0
 14 24  1  0
 15 21  2  0
 16 22  2  0
 17 23  2  0
 18 21  1  0
 19 22  1  0
 20 23  1  0
 25 26  1  0
 25 27  1  0
 26 28  2  0
 26 30  1  0
 27 29  2  0
 27 31  1  0
 28 32  1  0
 29 33  1  0
 30 34  2  0
 31 35  2  0
 32 36  2  0
 33 37  2  0
 34 36  1  0
 35 37  1  0
M  END)CTAB";
    const auto legacy_stereo = GENERATE(true, false);
    INFO("Legacy stereo perception == " << legacy_stereo);

    UseLegacyStereoPerceptionFixture reset_stereo_perception;
    Chirality::setUseLegacyStereoPerception(legacy_stereo);
    {
      // normal file parsing
      auto m = std::unique_ptr<RWMol>(MolBlockToMol(ctab));
      REQUIRE(m);
      CHECK(m->getBondWithIdx(5)->getStereo() == Bond::BondStereo::STEREOANY);
      CHECK(m->getBondWithIdx(6)->getStereo() == Bond::BondStereo::STEREONONE);
    }
    {
      // no sanitization during parsing
      bool sanitize = false;
      bool removeHs = false;
      auto m = std::unique_ptr<RWMol>(MolBlockToMol(ctab, sanitize, removeHs));
      REQUIRE(m);
      MolOps::sanitizeMol(*m);
      bool cleanIt = true;
      bool force = true;
      MolOps::assignStereochemistry(*m, cleanIt, force);
      CHECK(m->getBondWithIdx(5)->getStereo() == Bond::BondStereo::STEREOANY);
      CHECK(m->getBondWithIdx(6)->getStereo() == Bond::BondStereo::STEREONONE);
    }
  }
}

TEST_CASE(
    "github #3369: support new CIP code and StereoGroups in addStereoAnnotations()") {
  auto m1 =
      "C[C@@H]1N[C@H](C)[C@@H]([C@H](C)[C@@H]1C)C1[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]1C/C=C/C |a:5,o1:1,8,o2:14,16,&1:18,&2:3,6,r|"_smiles;
  REQUIRE(m1);
  SECTION("defaults") {
    ROMol m2(*m1);
    Chirality::addStereoAnnotations(m2);

    std::string txt;
    CHECK(m2.getAtomWithIdx(5)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "abs (S)");
    CHECK(m2.getAtomWithIdx(3)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "and2");
  }
  SECTION("double bonds") {
    ROMol m2(*m1);
    REQUIRE(m2.getBondBetweenAtoms(20, 21));
    m2.getBondBetweenAtoms(20, 21)->setStereo(Bond::BondStereo::STEREOTRANS);
    // initially no label is assigned since we have TRANS
    Chirality::addStereoAnnotations(m2);
    CHECK(
        !m2.getBondBetweenAtoms(20, 21)->hasProp(common_properties::bondNote));

    CIPLabeler::assignCIPLabels(m2);
    std::string txt;
    CHECK(m2.getBondBetweenAtoms(20, 21)->getPropIfPresent(
        common_properties::_CIPCode, txt));
    CHECK(txt == "E");
    Chirality::addStereoAnnotations(m2);
    CHECK(m2.getBondBetweenAtoms(20, 21)->getPropIfPresent(
        common_properties::bondNote, txt));
    CHECK(txt == "(E)");
  }

  SECTION("custom labels") {
    ROMol m2(*m1);
    CIPLabeler::assignCIPLabels(m2);

    std::string absLabel = "abs [{cip}]";
    std::string orLabel = "o{id} ({cip})";
    std::string andLabel = "&{id} ({cip})";
    std::string cipLabel = "[{cip}]";
    std::string bondLabel = "[{cip}]";
    Chirality::addStereoAnnotations(m2, absLabel, orLabel, andLabel, cipLabel,
                                    bondLabel);

    std::string txt;

    CHECK(m2.getAtomWithIdx(5)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "abs [S]");

    CHECK(m2.getAtomWithIdx(3)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "&2 (R)");

    CHECK(m2.getAtomWithIdx(1)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "o1 (S)");

    CHECK(m2.getAtomWithIdx(11)->getPropIfPresent(common_properties::atomNote,
                                                  txt));
    CHECK(txt == "[R]");

    CHECK(m2.getBondBetweenAtoms(20, 21)->getPropIfPresent(
        common_properties::bondNote, txt));
    CHECK(txt == "[E]");
  }

  SECTION("empty labels") {
    ROMol m2(*m1);
    CIPLabeler::assignCIPLabels(m2);

    std::string absLabel = "";
    std::string orLabel = "o{id} ({cip})";
    std::string andLabel = "&{id} ({cip})";
    std::string cipLabel = "[{cip}]";
    std::string bondLabel = "";
    Chirality::addStereoAnnotations(m2, absLabel, orLabel, andLabel, cipLabel,
                                    bondLabel);

    std::string txt;

    CHECK(!m2.getAtomWithIdx(5)->getPropIfPresent(common_properties::atomNote,
                                                  txt));

    CHECK(m2.getAtomWithIdx(3)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "&2 (R)");

    CHECK(m2.getAtomWithIdx(1)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "o1 (S)");

    CHECK(m2.getAtomWithIdx(11)->getPropIfPresent(common_properties::atomNote,
                                                  txt));
    CHECK(txt == "[R]");

    CHECK(!m2.getBondBetweenAtoms(20, 21)->getPropIfPresent(
        common_properties::bondNote, txt));
  }
}

TEST_CASE("do not wedge bonds to attachment points") {
  SECTION("basics") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0 ATTCHPT=1
M  V30 2 F -1.299038 0.750000 0.000000 0
M  V30 3 Cl -0.000000 -1.500000 0.000000 0
M  V30 4 O 1.299038 0.750000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 4);
    MolOps::expandAttachmentPoints(*m);
    CHECK(m->getNumAtoms() == 5);
    CHECK(m->getAtomWithIdx(4)->getTotalValence() == 1);
    Chirality::wedgeMolBonds(*m);
    CHECK(m->getBondBetweenAtoms(0, 4)->getBondDir() == Bond::BondDir::NONE);
  }
  SECTION("cage") {
    auto m = R"CTAB(
  Mrv2305 03052406362D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 9 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.125 2.5817 0 0 CFG=1
M  V30 2 C -4.4587 1.8117 0 0
M  V30 3 C -4.4587 0.2716 0 0
M  V30 4 C -3.125 -0.4984 0 0
M  V30 5 C -1.7913 0.2716 0 0
M  V30 6 N -1.7913 1.8117 0 0
M  V30 7 O -2.5357 1.1589 0 0
M  V30 8 * -3.125 4.1217 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 1 7
M  V30 7 1 7 4
M  V30 8 1 1 8
M  V30 9 1 1 6 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    m->getAtomWithIdx(7)->setProp(common_properties::_fromAttachPoint, 1);
    Chirality::wedgeMolBonds(*m);
    CHECK(m->getBondBetweenAtoms(0, 7)->getBondDir() == Bond::BondDir::NONE);
  }
}

TEST_CASE(
    "github #7203: Remove invalid stereo groups on MolOps::assignStereochemistry") {
  SECTION("single-atom groups") {
    UseLegacyStereoPerceptionFixture reset_stereo_perception;

    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);
      INFO(use_legacy);
      auto m = "C[C@H](N)[C@@H](C)C |o1:1,o2:3|"_smiles;
      REQUIRE(m);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
      CHECK(m->getAtomWithIdx(3)->getChiralTag() ==
            Atom::ChiralType::CHI_UNSPECIFIED);
      CHECK(m->getStereoGroups().size() == 1);
    }
  }
  SECTION("two-atom groups") {
    for (auto use_legacy : {false, true}) {
      Chirality::setUseLegacyStereoPerception(use_legacy);
      INFO(use_legacy);
      {  // no removal necessary
        auto m = "C[C@H](N)[C@@H](F)C |o1:1,3|"_smiles;
        REQUIRE(m);
        CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
              Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
        CHECK(m->getAtomWithIdx(3)->getChiralTag() ==
              Atom::ChiralType::CHI_TETRAHEDRAL_CW);
        CHECK(m->getStereoGroups().size() == 1);
        CHECK(m->getStereoGroups()[0].getAtoms().size() == 2);
      }
      {  // removal
        auto m = "C[C@H](N)[C@@H](C)C |o1:1,3|"_smiles;
        REQUIRE(m);
        CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
              Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
        CHECK(m->getAtomWithIdx(3)->getChiralTag() ==
              Atom::ChiralType::CHI_UNSPECIFIED);
        CHECK(m->getStereoGroups().size() == 1);
        CHECK(m->getStereoGroups()[0].getAtoms().size() == 1);
      }
    }
  }
}

TEST_CASE(
    "GitHub Issue #7346: Trigonal Pyramid Carbon may or not have a parity depending on atom ordering",
    "[bug]") {
  // First atom not in the same plane as the rest
  const auto m0 = R"CTAB(
                    3D

  5  4  0  0  1  0            999 V2000
   -0.2052    0.0804   -0.7454 F   0  0  0  0  0  0
   -0.5696    0.2231   -2.0688 C   0  0  0  0  0  0
   -2.2748   -0.0294   -1.1593 Br  0  0  0  0  0  0
    0.2659   -1.3655   -2.0764 Cl  0  0  0  0  0  0
   -0.1867    1.1526   -2.4961 D   0  0  0  0  0  0
  1  2  1  0  0  0
  2  3  1  0  0  0
  2  4  1  0  0  0
  2  5  1  0  0  0
M  END)CTAB"_ctab;

  // All atoms in the same plane except the rest
  const auto m1 = R"CTAB(
                    3D

  5  4  0  0  1  0            999 V2000
   -0.1867    1.1526   -2.4961 D   0  0  0  0  0  0
   -0.5696    0.2231   -2.0688 C   0  0  0  0  0  0
   -2.2748   -0.0294   -1.1593 Br  0  0  0  0  0  0
    0.2659   -1.3655   -2.0764 Cl  0  0  0  0  0  0
   -0.2052    0.0804   -0.7454 F   0  0  0  0  0  0
  1  2  1  0  0  0
  2  3  1  0  0  0
  2  4  1  0  0  0
  2  5  1  0  0  0
M  END)CTAB"_ctab;

  REQUIRE(m0);
  REQUIRE(m1);

  const auto cAt0 = m0->getAtomWithIdx(1);
  const auto cAt1 = m1->getAtomWithIdx(1);
  REQUIRE(cAt0->getAtomicNum() == 6);
  REQUIRE(cAt1->getAtomicNum() == 6);

  // Two atoms changes positions in the molblocks, but their coordinates
  // didn't change, so they must have opposite parities in order to
  // preserve the absolute chirality
  CHECK(cAt0->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  CHECK(cAt1->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW);

  CIPLabeler::assignCIPLabels(*m0);
  CIPLabeler::assignCIPLabels(*m1);

  CHECK(cAt0->getProp<std::string>(common_properties::_CIPCode) ==
        cAt1->getProp<std::string>(common_properties::_CIPCode));
}

TEST_CASE(
    "Github #7300: incorrect chiral carbon perception for some phosphates") {
  SECTION("basics") {
    auto m = "NC(CP(=O)(O)O)CP(=O)(O)O"_smiles;
    REQUIRE(m);
    auto sis = Chirality::findPotentialStereo(*m);
    CHECK(sis.empty());
  }
  SECTION("make sure Ps can still yield chiral centers") {
    auto m = "NC(CP(=O)(O)[O-])CP(=O)(O)[O-]"_smiles;
    REQUIRE(m);
    auto sis = Chirality::findPotentialStereo(*m);
    CHECK(sis.size() == 3);
  }
}

TEST_CASE("github 7371") {
  SECTION("basics") {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 20 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -0.024300 -1.051200 0.000000 0
M  V30 2 C 0.690200 -1.463700 0.000000 0
M  V30 3 O 0.690200 -2.288700 0.000000 0
M  V30 4 N 1.404500 -1.051200 0.000000 0
M  V30 5 C 2.119100 -1.463700 0.000000 0
M  V30 6 C 2.119100 -2.288700 0.000000 0
M  V30 7 C 1.404500 -2.701200 0.000000 0
M  V30 8 C 2.833500 -2.701200 0.000000 0
M  V30 9 C 3.547900 -2.288700 0.000000 0
M  V30 10 N 3.547900 -1.463700 0.000000 0
M  V30 11 C 2.833600 -1.051200 0.000000 0
M  V30 12 C 3.363800 -0.419100 0.000000 0
M  V30 13 C 4.176300 -0.562300 0.000000 0
M  V30 14 C 3.081600 0.356100 0.000000 0
M  V30 15 C 1.404500 -0.226000 0.000000 0
M  V30 16 N 2.155300 0.191900 0.000000 0
M  V30 17 C 2.161900 1.051300 0.000000 0
M  V30 18 C 1.417800 1.480900 0.000000 0
M  V30 19 C 0.676800 1.045500 0.000000 0
M  V30 20 C 0.690200 0.186400 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 4 2
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 1 6 8
M  V30 8 2 8 9
M  V30 9 1 9 10
M  V30 10 2 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 12 14
M  V30 14 1 4 15
M  V30 15 2 15 16
M  V30 16 1 16 17
M  V30 17 2 17 18
M  V30 18 1 18 19
M  V30 19 2 19 20
M  V30 20 1 5 11 CFG=1
M  V30 21 1 20 15
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(3)->getStereo() == Bond::BondStereo::STEREOATROPCW);
    // after reading in, there's no bond wedging:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
    }
    Chirality::wedgeMolBonds(*m, &m->getConformer());
    // now there is:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      if (bnd->getIdx() != 19) {
        CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
      } else {
        CHECK(bnd->getBondDir() == Bond::BondDir::BEGINWEDGE);
      }
    }
  }
  SECTION("3d") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 20 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.382755 -0.645531 -2.970173 0
M  V30 2 C 1.191902 -0.080895 -1.951313 0
M  V30 3 O 2.258061 0.469395 -2.326427 0
M  V30 4 N 0.799020 -0.142115 -0.557889 0
M  V30 5 C -0.412487 -0.780398 -0.245786 0
M  V30 6 C -0.494802 -2.127175 0.018641 0
M  V30 7 C 0.666160 -3.010952 -0.001627 0
M  V30 8 C -1.739073 -2.694470 0.315567 0
M  V30 9 C -2.888797 -1.963937 0.355103 0
M  V30 10 N -2.786157 -0.648937 0.093334 0
M  V30 11 C -1.586505 -0.059517 -0.200695 0
M  V30 12 C -1.679580 1.394735 -0.453406 0
M  V30 13 C -1.637015 2.096392 0.888082 0
M  V30 14 C -3.022840 1.757751 -1.085608 0
M  V30 15 C 1.696224 0.471157 0.348870 0
M  V30 16 N 1.608949 1.758263 0.730152 0
M  V30 17 C 2.457188 2.346328 1.581612 0
M  V30 18 C 3.501956 1.626784 2.127981 0
M  V30 19 C 3.644663 0.298166 1.774403 0
M  V30 20 C 2.724579 -0.259689 0.878536 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 4 2
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 1 6 8
M  V30 8 2 8 9
M  V30 9 1 9 10
M  V30 10 2 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 12 14
M  V30 14 1 4 15
M  V30 15 2 15 16
M  V30 16 1 16 17
M  V30 17 2 17 18
M  V30 18 1 18 19
M  V30 19 2 19 20
M  V30 20 1 5 11 CFG=1
M  V30 21 1 20 15
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getBondWithIdx(3)->getStereo() == Bond::BondStereo::STEREOATROPCW);
    // after reading in, there's no bond wedging:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
    }
    Chirality::wedgeMolBonds(*m, &m->getConformer());
    // now there is:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      if (bnd->getIdx() != 4) {
        CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
      } else {
        CHECK(bnd->getBondDir() == Bond::BondDir::BEGINWEDGE);
      }
    }
  }
  SECTION("favor larger rings") {
    auto m = R"CTAB(
  Mrv2401 04262410272D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 14 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.5796 -3.2152 0 0
M  V30 2 N 2.1037 -4.6798 0 0
M  V30 3 N 0.5637 -4.6798 0 0
M  V30 4 C 0.0878 -3.2152 0 0
M  V30 5 C -1.3768 -2.7393 0 0
M  V30 6 C 1.3337 -2.31 0 0
M  V30 7 C 1.3337 -0.77 0 0
M  V30 8 C 0 0 0 0
M  V30 9 Cl -1.3337 -0.77 0 0
M  V30 10 C 0 1.54 0 0
M  V30 11 C 1.3337 2.31 0 0
M  V30 12 C 2.6674 1.54 0 0
M  V30 13 C 2.6674 -0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 1 6 4 CFG=1
M  V30 6 2 1 6
M  V30 7 1 6 7
M  V30 8 2 7 8
M  V30 9 1 8 9
M  V30 10 1 8 10
M  V30 11 2 10 11
M  V30 12 1 11 12
M  V30 13 2 12 13
M  V30 14 1 7 13
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getBondBetweenAtoms(5, 6)->getStereo() ==
          Bond::BondStereo::STEREOATROPCW);
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
    }
    Chirality::wedgeMolBonds(*m, &m->getConformer());
    // now there is:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      if (bnd->getIdx() != 13) {
        CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
      } else {
        CHECK(bnd->getBondDir() == Bond::BondDir::BEGINWEDGE);
      }
    }
  }
  SECTION("favor larger rings 3D") {
    auto m = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 14 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.472629 -1.845227 -0.264711 0
M  V30 2 N 2.829375 -1.724790 -0.334942 0
M  V30 3 N 3.127980 -0.437811 -0.267564 0
M  V30 4 C 1.973684 0.284699 -0.153598 0
M  V30 5 C 1.905135 1.771086 -0.051971 0
M  V30 6 C 0.912735 -0.591792 -0.150090 0
M  V30 7 C -0.495016 -0.243797 -0.046092 0
M  V30 8 C -1.102614 -0.154763 1.189160 0
M  V30 9 Cl -0.084572 -0.483042 2.603743 0
M  V30 10 C -2.434893 0.170808 1.361855 0
M  V30 11 C -3.202022 0.420663 0.247454 0
M  V30 12 C -2.626093 0.340227 -1.000486 0
M  V30 13 C -1.279199 0.010190 -1.157641 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 1 6 4 CFG=3
M  V30 6 2 1 6
M  V30 7 1 6 7
M  V30 8 2 7 8
M  V30 9 1 8 9
M  V30 10 1 8 10
M  V30 11 2 10 11
M  V30 12 1 11 12
M  V30 13 2 12 13
M  V30 14 1 7 13
M  V30 END BOND
M  V30 END CTAB
M  END

)CTAB"_ctab;
    REQUIRE(m);
    CHECK(m->getBondBetweenAtoms(5, 6)->getStereo() ==
          Bond::BondStereo::STEREOATROPCW);
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
    }
    Chirality::wedgeMolBonds(*m, &m->getConformer());
    // now there is:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      if (bnd->getIdx() != 7) {
        CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
      } else {
        CHECK(bnd->getBondDir() == Bond::BondDir::BEGINWEDGE);
      }
    }
  }

  SECTION("avoid macrocycles") {
    auto m =
        "C1=C(C2=CCCCCCCCCCC2)CCC1 |(-1.71835,-1.22732,;-1.61655,-0.232318,;-0.751953,0.270082,;-0.754753,1.27008,;0.109847,1.77248,;0.977247,1.27488,;1.84185,1.77728,;2.70925,1.27968,;2.71205,0.279682,;1.84745,-0.222718,;1.85025,-1.22272,;0.985647,-1.72512,;0.118247,-1.22752,;0.115447,-0.227518,;-2.53135,0.171882,;-3.19835,-0.573118,;-2.69595,-1.43772,)|"_smiles;
    REQUIRE(m);
    m->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOATROPCW);
    Chirality::wedgeMolBonds(*m, &m->getConformer());
    // now there is:
    for (const auto bnd : m->bonds()) {
      INFO(bnd->getIdx());
      if (bnd->getIdx() != 13) {
        CHECK(bnd->getBondDir() == Bond::BondDir::NONE);
      } else {
        CHECK(bnd->getBondDir() == Bond::BondDir::BEGINWEDGE);
      }
    }
  }
}

TEST_CASE(
    "github #7434: planar amide nitrogen incorrectly flagged as _ChiralityPossible") {
  SECTION("as reported") {
    auto m1 = "O=C1CCCCC[C@@H]2CN1CCO2"_smiles;
    REQUIRE(m1);
    // new stereo
    CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
        m1->getAtomWithIdx(7)));
    CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
        m1->getAtomWithIdx(9)));

    // force old stereo
    UseLegacyStereoPerceptionFixture resetStereoPerception{true};
    bool cleanIt = true;
    bool flagPossible = true;
    auto si = Chirality::findPotentialStereo(*m1, cleanIt, flagPossible);
    CHECK(si.size() == 1);
  }
  SECTION("details") {
    auto m1 = "C1CCCCC[C@@H]2CN1CCO2"_smiles;
    REQUIRE(m1);
    CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
        m1->getAtomWithIdx(6)));
    CHECK(Chirality::detail::isAtomPotentialTetrahedralCenter(
        m1->getAtomWithIdx(8)));

    // force old stereo
    UseLegacyStereoPerceptionFixture resetStereoPerception{true};
    bool cleanIt = true;
    bool flagPossible = true;
    {
      auto si = Chirality::findPotentialStereo(*m1, cleanIt, flagPossible);
      CHECK(si.size() == 2);
    }
    m1->getAtomWithIdx(8)->setHybridization(Atom::HybridizationType::SP2);
    CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
        m1->getAtomWithIdx(8)));
    {
      auto si = Chirality::findPotentialStereo(*m1, cleanIt, flagPossible);
      CHECK(si.size() == 1);
    }
    m1->getAtomWithIdx(8)->setHybridization(Atom::HybridizationType::SP3);
    m1->getBondBetweenAtoms(0, 8)->setIsConjugated(true);
    CHECK(!Chirality::detail::isAtomPotentialTetrahedralCenter(
        m1->getAtomWithIdx(8)));
    {
      auto si = Chirality::findPotentialStereo(*m1, cleanIt, flagPossible);
      CHECK(si.size() == 1);
    }
  }
}

TEST_CASE("github #7438: expose code for simplified stereo labels") {
  SECTION("basics") {
    {
      auto m = "C[C@H](N)[C@@H](C)F |o1:1,3|"_smiles;
      REQUIRE(m);
      {
        ROMol m2(*m);
        Chirality::simplifyEnhancedStereo(m2);
        std::string note;
        CHECK(m2.getPropIfPresent(common_properties::molNote, note));
        CHECK(note == "OR enantiomer");
        CHECK(m2.getStereoGroups().empty());
      }
      {
        ROMol m2(*m);
        bool removeAffectedStereoGroups = false;
        Chirality::simplifyEnhancedStereo(m2, removeAffectedStereoGroups);
        std::string note;
        CHECK(m2.getPropIfPresent(common_properties::molNote, note));
        CHECK(note == "OR enantiomer");
        CHECK(m2.getStereoGroups().size() == 1);
      }
    }
    {
      auto m = "C[C@H](N)[C@@H](C)F |&1:1,3|"_smiles;
      REQUIRE(m);
      {
        ROMol m2(*m);
        Chirality::simplifyEnhancedStereo(m2);
        std::string note;
        CHECK(m2.getPropIfPresent(common_properties::molNote, note));
        CHECK(note == "AND enantiomer");
        CHECK(m2.getStereoGroups().empty());
      }
      {
        ROMol m2(*m);
        bool removeAffectedStereoGroups = false;
        Chirality::simplifyEnhancedStereo(m2, removeAffectedStereoGroups);
        std::string note;
        CHECK(m2.getPropIfPresent(common_properties::molNote, note));
        CHECK(note == "AND enantiomer");
        CHECK(m2.getStereoGroups().size() == 1);
      }
    }
  }

  SECTION("incomplete") {
    auto m = "C[C@H](N)[C@@H](C)F |o1:1|"_smiles;
    REQUIRE(m);
    {
      ROMol m2(*m);
      Chirality::simplifyEnhancedStereo(m2);
      CHECK(!m2.hasProp(common_properties::molNote));
      CHECK(m2.getStereoGroups().size() == 1);
    }
  }
}