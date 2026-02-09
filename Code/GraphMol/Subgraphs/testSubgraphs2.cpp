//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>

#include <algorithm>

using namespace std;
using namespace RDKit;

std::string canon(const std::string &smiles) {
  unique_ptr<ROMol> m(SmilesToMol(smiles));
  const bool useStereo = true;
  return MolToSmiles(*m, useStereo);
}

TEST_CASE("pathToSubmol") {
  std::string smiles = "CC1CC1";
  RWMol *mol = SmilesToMol(smiles);
  REQUIRE(mol);

  PATH_LIST sgs;
  sgs = findAllSubgraphsOfLengthN(*mol, 3, false, 0);
  REQUIRE(sgs.size() == 3);
  for (const auto &tmp : sgs) {
    REQUIRE(tmp.size() == 3);
    CHECK(tmp[0] == 0);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, tmp, false);
    smiles = MolToSmiles(*frag, true, false, 0, false);
    if (tmp[1] == 1) {
      if (tmp[2] == 2) {
        CHECK(smiles == "CCCC");
      } else if (tmp[2] == 3) {
        CHECK(smiles == "CC(C)C");
      } else {
        FAIL();
      }
    } else if (tmp[1] == 3) {
      if (tmp[2] == 2) {
        CHECK(smiles == "CCCC");
      } else if (tmp[2] == 1) {
        CHECK(smiles == "CC(C)C");
      } else {
        FAIL();
      }
    } else {
      FAIL();
    }
    delete frag;
  }
  delete mol;
}

TEST_CASE("Atom Environments") {
  {
    std::string smiles = "CC1CC1";
    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 1, 0);
    REQUIRE(pth.size() == 1);
    REQUIRE(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 2, 0);
    REQUIRE(pth.size() == 3);
    REQUIRE(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 3, 0);
    REQUIRE(pth.size() == 4);
    REQUIRE(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 4, 0);
    REQUIRE(pth.size() == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 1, 1);
    REQUIRE(pth.size() == 3);

    pth = findAtomEnvironmentOfRadiusN(*mol, 2, 1);
    REQUIRE(pth.size() == 4);

    pth = findAtomEnvironmentOfRadiusN(*mol, 3, 1);
    REQUIRE(pth.size() == 0);

    delete mol;
  }

  {
    std::string smiles = "CC1CC1";
    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);
    ROMol *mH = MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mH, 1, 0);
    REQUIRE(pth.size() == 1);
    REQUIRE(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mH, 1, 0, true);
    REQUIRE(pth.size() == 4);

    delete mol;
    delete mH;
  }

  {
    std::string smiles = "O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))";
    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);
    smiles = MolToSmiles(*mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 9);
    REQUIRE(pth.size() == 8);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, false);
    smiles = MolToSmiles(*frag, true, false, 0, false);
    REQUIRE(canon(smiles) == canon("C(C(C(O)C)C(C)C)C"));
    delete frag;
    delete mol;
  }
}

TEST_CASE("Atom Environments (Extension)") {
  {
    std::string smiles = "C=NC";
    unsigned int rootedAtAtom = 2;

    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);
    ROMol *mH = MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth =
        findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, false, true);
    REQUIRE(pth.size() == 2);
    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, true, true);
    REQUIRE(pth.size() == 5);

    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, false, true);
    REQUIRE(pth.size() == 0);
    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, true, true);
    REQUIRE(pth.size() == 7);
    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, true, false);
    REQUIRE(pth.size() == 7);
    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, false, false);
    REQUIRE(pth.size() == 2);

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, true);
    REQUIRE(pth.size() == 0);
    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, true);
    REQUIRE(pth.size() == 0);
    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, false);
    REQUIRE(pth.size() == 7);
    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, false);
    REQUIRE(pth.size() == 2);

    delete mol;
    delete mH;
  }

  {
    std::string smiles = "C=NO";  // C(0)=N(1)O(2)
    unsigned int rootedAtAtom = 2;
    std::unordered_map<unsigned int, unsigned int> cAtomMap;

    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);
    ROMol *mH = MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, false,
                                                 true, &cAtomMap);
    REQUIRE(cAtomMap.size() == 3);
    REQUIRE(cAtomMap[rootedAtAtom] == 0);
    REQUIRE(cAtomMap[1] == 1);
    REQUIRE(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, false, false,
                                       &cAtomMap);
    REQUIRE(cAtomMap.size() == 3);
    REQUIRE(cAtomMap[rootedAtAtom] == 0);
    REQUIRE(cAtomMap[1] == 1);
    REQUIRE(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, true, true,
                                       &cAtomMap);
    REQUIRE(cAtomMap.size() == 4);
    REQUIRE(cAtomMap[rootedAtAtom] == 0);
    REQUIRE(cAtomMap[1] == 1);
    REQUIRE(cAtomMap[5] == 1);
    REQUIRE(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, true, false,
                                       &cAtomMap);
    REQUIRE(cAtomMap.size() == 4);
    REQUIRE(cAtomMap[rootedAtAtom] == 0);
    REQUIRE(cAtomMap[1] == 1);
    REQUIRE(cAtomMap[5] == 1);
    REQUIRE(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, true,
                                       &cAtomMap);
    REQUIRE(pth.size() == 0);
    REQUIRE(cAtomMap.size() == 0);

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, true,
                                       &cAtomMap);
    REQUIRE(pth.size() == 0);
    REQUIRE(cAtomMap.size() == 0);

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, false,
                                       &cAtomMap);
    REQUIRE(pth.size() == 2);
    REQUIRE(cAtomMap.size() == 3);
    REQUIRE(cAtomMap[rootedAtAtom] == 0);
    REQUIRE(cAtomMap[1] == 1);
    REQUIRE(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, false,
                                       &cAtomMap);
    REQUIRE(pth.size() == 5);
    REQUIRE(cAtomMap.size() == 6);
    REQUIRE(cAtomMap[rootedAtAtom] == 0);
    REQUIRE(cAtomMap[1] == 1);
    REQUIRE(cAtomMap[5] == 1);
    REQUIRE(cAtomMap[0] == 2);
    REQUIRE(cAtomMap[3] == 3);
    REQUIRE(cAtomMap[4] == 3);
    cAtomMap.clear();

    delete mol;
    delete mH;
  }

  {
    std::string smiles = "c1cc[nH]c1";
    unsigned int rootedAtAtom = 1;
    std::unordered_map<unsigned int, unsigned int> cAtomMap;

    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);

    // This test worked on ring system to guarantee that the atom ID
    // is marked correctly with the search radius
    PATH_TYPE pth;
    unsigned int size;
    for (size = 2; size < 4; size++) {
      pth = findAtomEnvironmentOfRadiusN(*mol, size, rootedAtAtom, false, true,
                                         &cAtomMap);
      REQUIRE(cAtomMap.size() == 5);
      REQUIRE(cAtomMap[rootedAtAtom] == 0);
      REQUIRE(cAtomMap[0] == 1);
      REQUIRE(cAtomMap[2] == 1);
      REQUIRE(cAtomMap[3] == 2);
      REQUIRE(cAtomMap[4] == 2);
      cAtomMap.clear();
    }

    for (size = 4; size < 6; size++) {
      pth = findAtomEnvironmentOfRadiusN(*mol, size, rootedAtAtom, false, false,
                                         &cAtomMap);
      REQUIRE(cAtomMap.size() == 5);
      REQUIRE(cAtomMap[rootedAtAtom] == 0);
      REQUIRE(cAtomMap[0] == 1);
      REQUIRE(cAtomMap[2] == 1);
      REQUIRE(cAtomMap[3] == 2);
      REQUIRE(cAtomMap[4] == 2);
      cAtomMap.clear();
    }
    delete mol;
  }

  {
    std::string smiles = "c1cc[nH]c1";
    unsigned int rootedAtAtom = 1;
    std::unordered_map<unsigned int, unsigned int> cAtomMap;

    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);

    PATH_TYPE pth1 = findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom);
    PATH_TYPE pth2 = findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom, false);
    PATH_TYPE pth3 = findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom, true);
    PATH_TYPE pth4 =
        findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom, false, false);
    PATH_TYPE pth5 =
        findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom, false, true);
    PATH_TYPE pth6 = findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom, false,
                                                  false, &cAtomMap);
    cAtomMap.clear();
    PATH_TYPE pth7 = findAtomEnvironmentOfRadiusN(*mol, 2, rootedAtAtom, false,
                                                  true, &cAtomMap);
    cAtomMap.clear();
    delete mol;
  }
}

TEST_CASE("Github Issue103: stereochemistry and pathToSubmol") {
  {
    std::string smiles = "O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 12);
    REQUIRE(pth.size() == 5);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, false);
    smiles = MolToSmiles(*frag, true);
    REQUIRE(canon(smiles) == canon("C=CC(O)C=C"));
    delete frag;
    delete mol;
  }
  {
    std::string smiles = "O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 12);
    REQUIRE(pth.size() == 5);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, false);
    smiles = MolToSmarts(*frag);
    REQUIRE(smiles == "[#6]=[#6]-[#6@H](-[#8])-[#6]=[#6]");
    delete frag;
    delete mol;
  }
  {
    std::string smiles = "O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol = SmilesToMol(smiles);
    REQUIRE(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 12);
    REQUIRE(pth.size() == 5);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, true);
    smiles = MolToSmarts(*frag);
    REQUIRE(smiles == "[#6]=[#6]-[#6@](-[#8])-[#6]=[#6]");
    delete frag;
    delete mol;
  }
}

TEST_CASE(
    "Github Issue103: more stereochemistry and pathToSubmol (path needs to be in sorted order)") {
  std::string smiles = "I[C@](F)(Br)O";
  std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
  std::vector<int> path = {0, 3, 2, 1};
  const bool useQuery = false;
  std::unique_ptr<ROMol> mol2(Subgraphs::pathToSubmol(*mol, path, useQuery));
  REQUIRE(MolToSmiles(*mol2) == MolToSmiles(*mol));
}
