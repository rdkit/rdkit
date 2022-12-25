#include "RDGeneral/test.h"
#include "catch.hpp"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/DetermineBonds/DetermineBonds.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/AddHs.cpp>
#include <iostream>
#include <fstream>
#include <GraphMol/Resonance.h>

using namespace RDKit;

TEST_CASE("Determine Connectivity") {
  SECTION("Van der Waals") {
    unsigned int numTests = 39;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName =
          rdbase + "/Code/GraphMol/DetermineBonds/test_data/connectivity/" +
          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);

      determineConnectivity(*mol, false);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();

      REQUIRE(orig->getNumAtoms() == numAtoms);
      for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
          const auto origBond = orig->getBondBetweenAtoms(i, j);
          const auto molBond = mol->getBondBetweenAtoms(i, j);
          if (origBond) {
            CHECK(molBond);
          } else {
            CHECK(!molBond);
          }
        }
      }
    }
  }  // SECTION

  SECTION("Hueckel") {
    unsigned int numTests = 39;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName =
          rdbase + "/Code/GraphMol/DetermineBonds/test_data/connectivity/" +
          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);
      int charge = MolOps::getFormalCharge(*orig);

      determineConnectivity(*mol, true, charge);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();

      REQUIRE(orig->getNumAtoms() == numAtoms);
      for (unsigned int i = 0; i < numAtoms; i++) {
        for (unsigned int j = i + 1; j < numAtoms; j++) {
          const auto origBond = orig->getBondBetweenAtoms(i, j);
          const auto molBond = mol->getBondBetweenAtoms(i, j);
          if (origBond) {
            CHECK(molBond);
          } else {
            CHECK(!molBond);
          }
        }
      }
    }
  }  // SECTION

  SECTION("DetermineBondOrdering using charged fragments") {
    unsigned int numTests = 38;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName =
          rdbase +
          "/Code/GraphMol/DetermineBonds/test_data/charged_fragments/" +
          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);
      SmilesWriteParams params = {false, false, true, false, false, false, -1};
      std::string canonSmiles = MolToSmiles(*orig, params);
      int charge = MolOps::getFormalCharge(*orig);

      determineBonds(*mol, false, charge);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();
      REQUIRE(orig->getNumAtoms() == numAtoms);

      ResonanceMolSupplier resMolSuppl(
          *mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                    ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
      bool valid = false;
      for (unsigned int i = 0; i < resMolSuppl.length(); i++) {
        std::unique_ptr<ROMol> firstResMol(resMolSuppl[i]);
        std::unique_ptr<RWMol> resMol(new RWMol(*firstResMol));
        MolOps::setAromaticity(*resMol);

        std::string molSmiles = MolToSmiles(*resMol, params);
        if (molSmiles == canonSmiles) {
          CHECK(true);
          valid = true;
          break;
        }
      }
      if (!valid) {
        CHECK(false);
      }
    }
  }  // SECTION

  SECTION("DetermineBondOrdering using radicals") {
    unsigned int numTests = 10;
    for (unsigned int i = 0; i < numTests; i++) {
      std::string rdbase = getenv("RDBASE");
      std::string fName = rdbase +
                          "/Code/GraphMol/DetermineBonds/test_data/radicals/" +
                          "test" + std::to_string(i) + ".xyz";
      std::unique_ptr<RWMol> mol(XYZFileToMol(fName));
      REQUIRE(mol);
      std::string smiles = mol->getProp<std::string>("_FileComments");
      std::unique_ptr<RWMol> orig(SmilesToMol(smiles));
      REQUIRE(orig);
      SmilesWriteParams params = {false, false, true, false, false, false, -1};
      std::string canonSmiles = MolToSmiles(*orig, params);
      int charge = MolOps::getFormalCharge(*orig);

      determineBonds(*mol, false, charge, 1.3, false);
      MolOps::removeAllHs(*mol, false);

      auto numAtoms = mol->getNumAtoms();
      REQUIRE(orig->getNumAtoms() == numAtoms);

      ResonanceMolSupplier resMolSuppl(
          *mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                    ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
      bool valid = false;
      for (unsigned int i = 0; i < resMolSuppl.length(); i++) {
        std::unique_ptr<ROMol> firstResMol(resMolSuppl[i]);
        std::unique_ptr<RWMol> resMol(new RWMol(*firstResMol));
        MolOps::setAromaticity(*resMol);

        std::string molSmiles = MolToSmiles(*resMol, params);
        if (molSmiles == canonSmiles) {
          CHECK(true);
          valid = true;
          break;
        }
      }
      if (!valid) {
        CHECK(false);
      }
    }
  }  // SECTION
}

TEST_CASE("Github #5894: extra stereocenters assigned") {
  SECTION("as reported") {
    std::string xyz = R"XYZ(24

C	  -2.321570	  -2.154410	   0.000001
N	  -1.361840	  -1.072025	   0.000000
C	  -0.022813	  -1.380036	   0.000000
O	   0.540584	  -2.460253	  -0.000001
C	   0.939066	  -0.223023	   0.000000
C	   0.424348	   1.032682	  -0.000001
N	   1.421181	   1.964870	   0.000000
C	   2.538290	   1.259598	   0.000003
N	   2.287404	  -0.088227	   0.000001
C	   3.264961	  -1.144161	  -0.000001
N	  -0.930470	   1.289618	  -0.000003
C	  -1.392807	   2.667466	  -0.000002
C	  -1.856929	   0.240323	   0.000000
O	  -3.074555	   0.448852	   0.000002
H	  -2.955990	  -2.068408	  -0.888068
H	  -2.955969	  -2.068422	   0.888087
H	  -1.835402	  -3.133634	  -0.000012
H	   3.539988	   1.670786	   0.000006
H	   3.124842	  -1.749165	   0.899229
H	   4.268839	  -0.711765	  -0.000003
H	   3.124838	  -1.749164	  -0.899232
H	  -1.013179	   3.171571	  -0.894505
H	  -2.483808	   2.734391	  -0.000015
H	  -1.013199	   3.171563	   0.894514
)XYZ";
    std::unique_ptr<RWMol> m(XYZBlockToMol(xyz));
    REQUIRE(m);
    determineBonds(*m);
    CHECK(m->getAtomWithIdx(0)->getChiralTag() ==
          Atom::ChiralType::CHI_UNSPECIFIED);
  }
  SECTION("make sure chirality does still get assigned when appropriate") {
    std::string xyz = R"XYZ(8

C     -0.817192   -0.078725   -0.028772
C      0.680862    0.048830    0.076987
F      0.990227    0.019664    1.437282
Cl     1.147716    1.625471   -0.563047
Br     1.617608   -1.365187   -0.810838
H     -1.246798    0.864941    0.386221
H     -1.184702   -0.165336   -1.061440
H     -1.187721   -0.949657    0.563607
)XYZ";
    std::unique_ptr<RWMol> m(XYZBlockToMol(xyz));
    REQUIRE(m);
    determineBonds(*m);
    CHECK(m->getAtomWithIdx(1)->getChiralTag() ==
          Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
}