#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
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

      bool useHueckel = false;
      int charge = 0;
      double factor = 1.3;
      bool useVdw = true;
      determineConnectivity(*mol, useHueckel, charge, factor, useVdw);
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

  SECTION("connect the dots") {
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
      bool useHueckel = false;
      int charge = 0;
      double factor = 1.3;
      bool useVdw = false;
      determineConnectivity(*mol, useHueckel, charge, factor, useVdw);
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
#ifdef RDK_BUILD_YAEHMOP_SUPPORT
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
#endif
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

TEST_CASE("Github #5902: memory-hungry iterating") {
  SECTION("as reported") {
    std::string xyz = R"XYZ(102

C	 -16.385236	   2.590667	  -1.663038
O	 -15.382835	   2.646285	  -0.960157
C	 -17.486770	   1.590496	  -1.480412
C	 -17.426479	   0.848393	  -0.171161
O	 -17.870759	   1.354322	   0.863406
C	 -16.784469	  -0.515016	  -0.145056
C	 -15.293614	  -0.428148	  -0.317466
O	 -14.787174	  -0.463307	  -1.447143
C	 -14.460481	  -0.183298	   0.922568
C	 -13.219203	  -1.029666	   0.959157
O	 -13.240340	  -2.109801	   1.554415
C	 -12.037401	  -0.539162	   0.163235
C	 -10.675409	  -0.899627	   0.728926
O	 -10.353810	  -2.052168	   1.041750
C	  -9.678800	   0.216332	   0.848768
C	  -9.147864	   0.649135	  -0.485909
O	  -9.912272	   1.071583	  -1.359994
C	  -7.690467	   0.449234	  -0.773698
C	  -6.741448	   1.277715	   0.057015
O	  -7.004096	   1.629209	   1.199682
C	  -5.435743	   1.642402	  -0.599522
C	  -4.255239	   1.805063	   0.325633
O	  -3.375176	   2.631013	   0.050431
C	  -4.126808	   0.890066	   1.513706
C	  -2.695932	   0.672215	   1.954409
O	  -2.307667	   0.973239	   3.085947
C	  -1.766786	  -0.029114	   0.995972
C	  -1.859422	  -1.516153	   1.175568
O	  -2.957466	  -2.072109	   1.268212
C	  -0.569966	  -2.301235	   1.262650
C	   0.221471	  -2.266086	  -0.018571
O	  -0.211035	  -1.709162	  -1.027850
C	   1.555492	  -2.971945	   0.014884
C	   2.524452	  -2.659473	  -1.104628
O	   3.264932	  -3.547770	  -1.518099
C	   2.620413	  -1.237289	  -1.608953
C	   3.910719	  -0.925729	  -2.341491
O	   3.900504	  -0.522439	  -3.508184
C	   5.209118	  -1.007197	  -1.576049
C	   5.573147	   0.323512	  -0.961128
O	   4.700920	   1.082029	  -0.511564
C	   7.034888	   0.723585	  -0.905299
C	   7.889661	  -0.178933	  -0.039347
O	   7.444628	  -1.248384	   0.388334
C	   9.281043	   0.321054	   0.264631
C	  10.122198	  -0.671110	   1.023787
O	  10.081243	  -0.746355	   2.253367
C	  11.026182	  -1.578678	   0.231610
C	  12.195542	  -0.837792	  -0.351694
O	  12.271585	  -0.660090	  -1.574015
C	  13.266671	  -0.338377	   0.593379
C	  14.060350	   0.796605	  -0.011446
O	  13.666866	   1.962103	   0.067237
C	  15.329669	   0.450233	  -0.751762
C	  16.556123	   0.572388	   0.121941
O	  16.513756	   0.263215	   1.321202
C	  17.847471	   0.947579	  -0.556078
C	  18.654332	   1.924866	   0.240100
O	  18.228761	   3.055525	   0.471293
C	  20.010064	   1.488749	   0.720576
H	 -16.494708	   3.261706	  -2.536900
H	 -18.442907	   2.131691	  -1.518739
H	 -17.488497	   0.888823	  -2.321515
H	 -17.029335	  -1.012326	   0.800140
H	 -17.211375	  -1.134412	  -0.938223
H	 -15.042339	  -0.403273	   1.821906
H	 -14.187361	   0.876763	   0.981376
H	 -12.095859	   0.552470	   0.070144
H	 -12.093974	  -0.980080	  -0.838512
H	  -8.865857	  -0.137720	   1.483178
H	 -10.135800	   1.076876	   1.340128
H	  -7.529359	   0.701228	  -1.829927
H	  -7.449930	  -0.612058	  -0.654954
H	  -5.175830	   0.885782	  -1.343615
H	  -5.608620	   2.601535	  -1.122150
H	  -4.674700	   1.328318	   2.342078
H	  -4.557640	  -0.078949	   1.259616
H	  -0.749695	   0.315160	   1.189385
H	  -2.007106	   0.220191	  -0.039496
H	  -0.812981	  -3.345089	   1.505133
H	   0.035142	  -1.881843	   2.078326
H	   1.351982	  -4.048535	  -0.031261
H	   2.057067	  -2.740493	   0.959212
H	   2.544216	  -0.539871	  -0.757467
H	   1.790939	  -1.072278	  -2.303407
H	   5.987490	  -1.330007	  -2.267575
H	   5.141386	  -1.749617	  -0.788286
H	   7.089611	   1.745875	  -0.514724
H	   7.442524	   0.749886	  -1.927161
H	   9.218256	   1.229224	   0.873587
H	   9.787357	   0.591511	  -0.674079
H	  10.427055	  -2.026073	  -0.569973
H	  11.386156	  -2.408966	   0.846235
H	  13.933970	  -1.165689	   0.850429
H	  12.808573	   0.029000	   1.515898
H	  15.271198	  -0.580045	  -1.127757
H	  15.413536	   1.114025	  -1.621868
H	  18.422549	   0.036295	  -0.723464
H	  17.663459	   1.399062	  -1.532264
H	  19.934486	   0.489546	   1.157926
H	  20.693645	   1.489300	  -0.120543
H	  20.371657	   2.180532	   1.492305
)XYZ";
    std::unique_ptr<RWMol> m(XYZBlockToMol(xyz));
    REQUIRE(m);
    determineBonds(*m);
    CHECK(m->getNumAtoms() == 102);
  }
}

TEST_CASE("Github #6121: Single Atom in DetermineBonds") {
  SECTION("as reported") {
    std::string xyz = R"XYZ(1

H   0.0         0.0           0.0
)XYZ";
    std::unique_ptr<RWMol> m(XYZBlockToMol(xyz));
    REQUIRE(m);
    determineBonds(*m);
  }
}

TEST_CASE("github 6961: P-H bonds not found in phosphine") {
  SECTION("as reported") {
    std::string xyz = R"XYZ(4
xyz file
P 9.9999767321286015 9.9999428968490651 9.9298216136095618 
H 8.8082284983002523 9.9999330478847508 10.7216030817151875 
H 10.5974890657086007 11.0338788274478361 10.7168666854072114 
H 10.5976057038625981 8.9661452278177478 10.7170086192680003)XYZ";
    std::unique_ptr<RWMol> m(XYZBlockToMol(xyz));
    REQUIRE(m);
    determineBonds(*m);
    CHECK(m->getNumBonds() == 3);
    CHECK(m->getBondBetweenAtoms(0, 1));
    CHECK(m->getBondBetweenAtoms(0, 2));
    CHECK(m->getBondBetweenAtoms(0, 3));
  }
}

TEST_CASE(
    "github #7299: DetermineBondOrders() does not assign single bonds correctly") {
  SECTION("as reported") {
    RWMol m;
    m.addAtom(new Atom(6));
    m.addAtom(new Atom(8));
    m.addAtom(new Atom(8));
    m.addAtom(new Atom(8));
    m.addAtom(new Atom(1));
    m.addAtom(new Atom(1));
    m.addBond(0, 1, Bond::UNSPECIFIED);
    m.addBond(0, 2, Bond::UNSPECIFIED);
    m.addBond(0, 3, Bond::UNSPECIFIED);
    m.addBond(2, 4, Bond::UNSPECIFIED);
    m.addBond(3, 5, Bond::UNSPECIFIED);
    int charge = 0;
    bool allowChargedFragments = true;
    bool embedChiral = false;
    determineBondOrders(m, charge, allowChargedFragments, embedChiral);
    for (auto bnd : m.bonds()) {
      if (bnd->getIdx()) {
        CHECK(bnd->getBondType() == Bond::SINGLE);
      } else {
        CHECK(bnd->getBondType() == Bond::DOUBLE);
      }
    }
  }
}

TEST_CASE(
    "github #7331: DetermineBondOrders() makes incorrect assumptions about valence") {
  SECTION("as reported") {
    // do not anything here that needs implicit Hs
    std::vector<std::string> smiles = {
        "O=NO[Cl+][O-]",
        "[O-][I+3]([O-])([O-])[O-]",
        "[O-][I+2]([O-])[O-]",
        "F[P-](F)(F)(F)(F)F",
        "F[C+](F)F",
        "F[C-](F)F",
        "F[N+](F)(F)F",
        "F[N-]F",
        "F[Cl+]F",
        "F[Br+]F",
        "O=[Cl+]",
    };
    for (const auto &smi : smiles) {
      INFO(smi);
      auto m = v2::SmilesParse::MolFromSmiles(smi);
      REQUIRE(m);
      int charge = 0;
      for (auto atom : m->atoms()) {
        charge += atom->getFormalCharge();
      }
      bool allowChargedFragments = true;
      bool embedChiral = false;
      RWMol m2(*m);
      determineBondOrders(m2, charge, allowChargedFragments, embedChiral);
      for (auto bnd : m2.bonds()) {
        CHECK(bnd->getBondType() ==
              m->getBondWithIdx(bnd->getIdx())->getBondType());
      }
    }
  }
}