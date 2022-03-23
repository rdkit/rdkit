//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>

#include <iostream>
#include <algorithm>

using namespace std;
using namespace RDKit;

std::string canon(std::string smiles) {
  unique_ptr<ROMol> m(SmilesToMol(smiles));
  const bool useStereo = true;
  return MolToSmiles(*m, useStereo);
}

void test1() {
  std::cout << "-----------------------\n Test1: pathToSubmol" << std::endl;
  {
    std::string smiles = "CC1CC1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_LIST sgs;
    sgs = findAllSubgraphsOfLengthN(*mol, 3, false, 0);
    TEST_ASSERT(sgs.size() == 3);
    for (const auto &tmp : sgs) {
      TEST_ASSERT(tmp[0] == 0);
      TEST_ASSERT(tmp.size() == 3);
      ROMol *frag = Subgraphs::pathToSubmol(*mol, tmp, false);
      smiles = MolToSmiles(*frag, true, false, 0, false);
      if (tmp[1] == 1) {
        if (tmp[2] == 2) {
          TEST_ASSERT(smiles == "CCCC");
        } else if (tmp[2] == 3) {
          TEST_ASSERT(smiles == "CC(C)C");
        } else {
          TEST_ASSERT(0);
        }
      } else if (tmp[1] == 3) {
        if (tmp[2] == 2) {
          TEST_ASSERT(smiles == "CCCC");
        } else if (tmp[2] == 1) {
          TEST_ASSERT(smiles == "CC(C)C");
        } else {
          TEST_ASSERT(0);
        }
      } else {
        TEST_ASSERT(0);
      }
      delete frag;
    }
    delete mol;
  }
  std::cout << "Finished" << std::endl;
}

void test2() {
  std::cout << "-----------------------\n Test2: Atom Environments"
            << std::endl;
  {
    std::string smiles = "CC1CC1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 1, 0);
    TEST_ASSERT(pth.size() == 1);
    TEST_ASSERT(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 2, 0);
    TEST_ASSERT(pth.size() == 3);
    TEST_ASSERT(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 3, 0);
    TEST_ASSERT(pth.size() == 4);
    TEST_ASSERT(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 4, 0);
    TEST_ASSERT(pth.size() == 0);

    pth = findAtomEnvironmentOfRadiusN(*mol, 1, 1);
    TEST_ASSERT(pth.size() == 3);

    pth = findAtomEnvironmentOfRadiusN(*mol, 2, 1);
    TEST_ASSERT(pth.size() == 4);

    pth = findAtomEnvironmentOfRadiusN(*mol, 3, 1);
    TEST_ASSERT(pth.size() == 0);

    delete mol;
  }

  {
    std::string smiles = "CC1CC1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    ROMol *mH = MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mH, 1, 0);
    TEST_ASSERT(pth.size() == 1);
    TEST_ASSERT(pth[0] == 0);

    pth = findAtomEnvironmentOfRadiusN(*mH, 1, 0, true);
    TEST_ASSERT(pth.size() == 4);

    delete mol;
    delete mH;
  }

  {
    std::string smiles = "O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    smiles = MolToSmiles(*mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 9);
    TEST_ASSERT(pth.size() == 8);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, false);
    smiles = MolToSmiles(*frag, true, false, 0, false);
    TEST_ASSERT(canon(smiles) == canon("C(C(C(O)C)C(C)C)C"));
    delete frag;
    delete mol;
  }

  std::cout << "Finished" << std::endl;
}

void test3() {
  std::cout << "-----------------------" << std::endl;
  std::cout << "Test 3: Atom Environments (Extension)" << std::endl;

  {
    std::string smiles = "C=NC";
    unsigned int rootedAtAtom = 2;

    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    ROMol *mH = MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth =
        findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, false, true);
    TEST_ASSERT(pth.size() == 2);
    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, true, true);
    TEST_ASSERT(pth.size() == 5);

    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, false, true);
    TEST_ASSERT(pth.size() == 0);
    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, true, true);
    TEST_ASSERT(pth.size() == 7);
    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, true, false);
    TEST_ASSERT(pth.size() == 7);
    pth = findAtomEnvironmentOfRadiusN(*mH, 3, rootedAtAtom, false, false);
    TEST_ASSERT(pth.size() == 2);

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, true);
    TEST_ASSERT(pth.size() == 0);
    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, true);
    TEST_ASSERT(pth.size() == 0);
    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, false);
    TEST_ASSERT(pth.size() == 7);
    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, false);
    TEST_ASSERT(pth.size() == 2);

    delete mol;
    delete mH;
  }

  {
    std::string smiles = "C=NO";  // C(0)=N(1)O(2)
    unsigned int rootedAtAtom = 2;
    std::unordered_map<unsigned int, unsigned int> cAtomMap;

    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);
    ROMol *mH = MolOps::addHs(static_cast<const ROMol &>(*mol));

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, false,
                                                 true, &cAtomMap);
    TEST_ASSERT(cAtomMap.size() == 3);
    TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
    TEST_ASSERT(cAtomMap[1] == 1);
    TEST_ASSERT(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, false, false,
                                       &cAtomMap);
    TEST_ASSERT(cAtomMap.size() == 3);
    TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
    TEST_ASSERT(cAtomMap[1] == 1);
    TEST_ASSERT(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, true, true,
                                       &cAtomMap);
    TEST_ASSERT(cAtomMap.size() == 4);
    TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
    TEST_ASSERT(cAtomMap[1] == 1);
    TEST_ASSERT(cAtomMap[5] == 1);
    TEST_ASSERT(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 2, rootedAtAtom, true, false,
                                       &cAtomMap);
    TEST_ASSERT(cAtomMap.size() == 4);
    TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
    TEST_ASSERT(cAtomMap[1] == 1);
    TEST_ASSERT(cAtomMap[5] == 1);
    TEST_ASSERT(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, true,
                                       &cAtomMap);
    TEST_ASSERT(pth.size() == 0);
    TEST_ASSERT(cAtomMap.size() == 0);

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, true,
                                       &cAtomMap);
    TEST_ASSERT(pth.size() == 0);
    TEST_ASSERT(cAtomMap.size() == 0);

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, false, false,
                                       &cAtomMap);
    TEST_ASSERT(pth.size() == 2);
    TEST_ASSERT(cAtomMap.size() == 3);
    TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
    TEST_ASSERT(cAtomMap[1] == 1);
    TEST_ASSERT(cAtomMap[0] == 2);
    cAtomMap.clear();

    pth = findAtomEnvironmentOfRadiusN(*mH, 4, rootedAtAtom, true, false,
                                       &cAtomMap);
    TEST_ASSERT(pth.size() == 5);
    TEST_ASSERT(cAtomMap.size() == 6);
    TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
    TEST_ASSERT(cAtomMap[1] == 1);
    TEST_ASSERT(cAtomMap[5] == 1);
    TEST_ASSERT(cAtomMap[0] == 2);
    TEST_ASSERT(cAtomMap[3] == 3);
    TEST_ASSERT(cAtomMap[4] == 3);
    cAtomMap.clear();

    delete mol;
    delete mH;
  }

  {
    std::string smiles = "c1cc[nH]c1";
    unsigned int rootedAtAtom = 1;
    std::unordered_map<unsigned int, unsigned int> cAtomMap;

    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    // This test worked on ring system to guarantee that the atom ID
    // is marked correctly with the search radius
    PATH_TYPE pth;
    unsigned int size;
    for (size = 2; size < 4; size++) {
      pth = findAtomEnvironmentOfRadiusN(*mol, size, rootedAtAtom, false, true,
                                         &cAtomMap);
      TEST_ASSERT(cAtomMap.size() == 5);
      TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
      TEST_ASSERT(cAtomMap[0] == 1);
      TEST_ASSERT(cAtomMap[2] == 1);
      TEST_ASSERT(cAtomMap[3] == 2);
      TEST_ASSERT(cAtomMap[4] == 2);
      cAtomMap.clear();
    }

    for (size = 4; size < 6; size++) {
      pth = findAtomEnvironmentOfRadiusN(*mol, size, rootedAtAtom, false, false,
                                         &cAtomMap);
      TEST_ASSERT(cAtomMap.size() == 5);
      TEST_ASSERT(cAtomMap[rootedAtAtom] == 0);
      TEST_ASSERT(cAtomMap[0] == 1);
      TEST_ASSERT(cAtomMap[2] == 1);
      TEST_ASSERT(cAtomMap[3] == 2);
      TEST_ASSERT(cAtomMap[4] == 2);
      cAtomMap.clear();
    }
    delete mol;
  }

  {
    std::string smiles = "c1cc[nH]c1";
    unsigned int rootedAtAtom = 1;
    std::unordered_map<unsigned int, unsigned int> cAtomMap;

    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

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
  std::cout << "Finished" << std::endl;
}

void testGithubIssue103() {
  std::cout << "-----------------------\n Testing github Issue103: "
               "stereochemistry and pathToSubmol"
            << std::endl;
  {
    std::string smiles = "O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 12);
    TEST_ASSERT(pth.size() == 5);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, false);
    smiles = MolToSmiles(*frag, true);
    TEST_ASSERT(canon(smiles) == canon("C=CC(O)C=C"));
    delete frag;
    delete mol;
  }
  {
    std::string smiles = "O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 12);
    TEST_ASSERT(pth.size() == 5);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, false);
    smiles = MolToSmarts(*frag);
    TEST_ASSERT(smiles == "[#6]=[#6]-[#6@H](-[#8])-[#6]=[#6]");
    delete frag;
    delete mol;
  }
  {
    std::string smiles = "O=C(O)C(=O)C[C@@]1(C(=O)O)C=C[C@H](O)C=C1";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    PATH_TYPE pth = findAtomEnvironmentOfRadiusN(*mol, 2, 12);
    TEST_ASSERT(pth.size() == 5);
    ROMol *frag = Subgraphs::pathToSubmol(*mol, pth, true);
    smiles = MolToSmarts(*frag);
    TEST_ASSERT(smiles == "[#6]=[#6]-[#6@](-[#8])-[#6]=[#6]");
    delete frag;
    delete mol;
  }

  std::cout << "Finished" << std::endl;
}

void testGithubIssue2647() {
  std::cout << "-----------------------\n Testing github Issue103: "
               "more stereochemistry and pathToSubmol (path needs to be in "
               "sorted order)"
            << std::endl;
  std::string smiles = "I[C@](F)(Br)O";
  std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
  std::vector<int> path = {0, 3, 2, 1};
  const bool useQuery = false;
  std::unique_ptr<ROMol> mol2(Subgraphs::pathToSubmol(*mol, path, useQuery));
  TEST_ASSERT(MolToSmiles(*mol2) == MolToSmiles(*mol));
  std::cout << "Finished" << std::endl;
}

// -------------------------------------------------------------------
int main() {
  test1();
  test2();
  test3();
  testGithubIssue103();
  testGithubIssue2647();
  return 0;
}
