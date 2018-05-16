// $Id$
//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/test.h>
#include <list>
#include <vector>
#include <string>
#include <stdio.h>
#include <ctype.h>
#include "../RDKitBase.h"
#include "../SmilesParse/SmilesParse.h"
#include "MolHash.h"

using namespace RDKit::MolHash;
namespace RDKit {

void test1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test1 DEFAULT ARGUMENTS"
                       << std::endl;
  std::cout << "Hash size = " << 8 * sizeof(HashCodeType) << " bits.\n";
  const char* smi[] = {
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1",
      "CN(C)c1ccc(CC(=O)NCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1",
      "CN(C)c1ccc(CC(NCCCCCC(NO)=O)=O)cc1",
      "CC(C)Cc1ccc(C(C)C(=O)NC23CC4CC(C2)CC(C3)C4)cc1",
      "c1cc([N+]([O-])=O)ccc1CC(=O)NC1CCCCCC1",
      "CC1(C)NC(C)(C)CC(NC(=O)Cc2ccccc2)C1",
      "CC(C)CC(NC(CNC(CNC(C(Cc1c2ccccc2[nH]c1)NC(C(Cc1cnc[nH]1)NC(CNC(C(C(C)O)"
      "NC(C(C(C)(C)S)NC(C(Cc1ccccc1)NC(C(CCCNC(=N)N)NC(C(N)CCC(N)=O)=O)=O)=O)="
      "O)=O)=O)=O)=O)=O)=O)C(NC(C(=O)N1CCCC1C(=O)NC(CS)C(NC(CC(N)=O)C(NCC(=O)"
      "N1CCCC1C(O)=O)=O)=O)Cc1ccc(O)cc1)=O",  // CHEMBL527084
      "CCCCCCC1C23C4=c5c6c7c8c5c5c9c%10c%11c%12c(c%108)c8c7c7c%10c%13c%14c%15c%"
      "16c%17c%18c%19c%20c(c%21c%22c%23c(c9C(C25C[N+]1(C)C)C1c2c3c3c5c9c2-c(c%"
      "231)c(c%22%19)C%18C9C1(C5=C%13C(C43)c6%10)C%14%17C[N+](C)(C)C1CCCCCC)c%"
      "21%11)c%12c(c%16%20)c8c%157",  // CHEMBL439119
      "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)"
      "NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)="
      "O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O",  // CHEMBL439258
      "NCCCCC(NC(CN)=O)C(NCC(NC(CC(C)C)C(=O)N1Cc2ccccc2CC1C(N1CC2CCCCC2C1C(NCC("
      "NC(CC(C)C)C(=O)N1Cc2ccccc2CC1C(=O)N1CC2CCCCC2C1C(NCC(NC(C("
      "N1Cc2ccccc2CC1C(N1CC2CCCCC2C1C(NCC(NC(CC(C)C)C(=O)N1Cc2ccccc2CC1C(=O)"
      "N1CC2CCCCC2C1C(NCC(NC(C(N1Cc2c(cccc2)CC1C(N1CC2CCCCC2C1C(NCC(NC(CC(C)C)"
      "C(=O)N1Cc2ccccc2CC1C(N1CC2CCCCC2C1C(NCC(NC(CC(C)C)C(NC(C(N)=O)CCCNC(=N)"
      "N)=O)=O)=O)=O)=O)=O)=O)=O)CCCCN)=O)=O)=O)=O)=O)=O)CCCCN)=O)=O)=O)=O)=O)="
      "O)=O",  // CHEMBL441746
      "CC(C)CC(NC(C(C(C)C)NC(C(N)CCC(O)=O)=O)=O)C(NC(C(O)C(=O)NC(CC(O)=O)C(NCC("
      "=O)NC(CCC(O)=O)C(NC(C(O)=O)Cc1ccccc1)=O)=O)Cc1ccccc1)=O",  // CHEMBL384606
      "CCC(C)C1C(=O)N2CCCC2C(=O)NC2CSSCC3NC(=O)C(C(C)C)NC(=O)C(CCCCN)NC(=O)C("
      "CC(N)=O)NC(=O)C(CCCCN)NC(=O)C4CSSCC(C(=O)NC(C(C)C)C(=O)NC(Cc5ccccc5)C(="
      "O)N1)NC(=O)C(CO)NC(=O)C(CCC(O)=O)NC(=O)CNC(=O)C(NC(=O)C1CCCN1C(=O)C(C(C)"
      "CC)NC(=O)CNC(=O)C(CC(N)=O)NC(=O)C(CCCNC(=N)N)NC(=O)C(Cc1ccc(O)cc1)NC3=O)"
      "CSSCC(NC(=O)CNC(=O)C(C)NC(=O)C(C(C)C)NC(=O)C(C(C)O)NC(=O)C(C(C)O)NC(=O)"
      "C(CC(C)C)NC2=O)C(=O)NC(CO)C(=O)N4",  // CHEMBL526869
  };
  for (auto& i : smi) {
    ROMOL_SPTR mol = ROMOL_SPTR(SmilesToMol(i));
    std::vector<unsigned> atomsToUse;
    std::vector<unsigned> bondsToUse;
    std::vector<boost::uint32_t> atomCodes(mol->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol->getNumBonds());

    unsigned n;
    n = mol->getNumAtoms();
    atomsToUse.resize(n);
    for (unsigned i = 0; i < n; i++) atomsToUse[i] = i;
    n = mol->getNumBonds();
    bondsToUse.resize(n);
    for (unsigned i = 0; i < n; i++) bondsToUse[i] = i;

    n = mol->getNumAtoms();
    for (unsigned i = 0; i < n; i++)
      atomCodes[i] =
          1;  // + mol->getAtomWithIdx(i)->getAtomicNum(); //res0 != res1,2,3
    n = mol->getNumBonds();
    for (unsigned i = 0; i < n; i++) bondCodes[i] = 1;

    fillAtomBondCodes(*mol, CF_NO_LABELS, &atomCodes, &bondCodes);

    HashCodeType res0 = generateMoleculeHashCode(*mol);
    HashCodeType res1 = generateMoleculeHashCode(*mol, &atomsToUse, nullptr,
                                                 &atomCodes, &bondCodes);
    HashCodeType res2 = generateMoleculeHashCode(*mol, nullptr, &bondsToUse,
                                                 &atomCodes, &bondCodes);
    HashCodeType res3 = generateMoleculeHashCode(*mol, &atomsToUse, &bondsToUse,
                                                 &atomCodes, &bondCodes);

    std::cout << res0 << " = " << encode(&res0, sizeof(res0)) << std::endl;
    std::cout << res1 << " = " << encode(&res1, sizeof(res1)) << std::endl;
    std::cout << res2 << " = " << encode(&res2, sizeof(res2)) << std::endl;
    std::cout << res3 << " = " << encode(&res3, sizeof(res3)) << std::endl
              << std::endl;

    //            bool passed = 0 != res0 && res0 == res1 && res0 == res2 &&
    //            res0 == res3;
    //            TEST_ASSERT(passed);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test2 CHIRALITY == ATOM"
                       << std::endl;
  std::cout << "Hash size = " << 8 * sizeof(HashCodeType) << " bits.\n";
  const char* smi[] = {
      // equal non-chiral hash
      "C[C@H](F)Cl",   "C[C@@H](F)Cl", "CC(F)Cl",
      "[13CH3]C(F)Cl", "C[C@H](Cl)F",  "C[C@@H](Cl)F",
  };

  std::vector<HashCodeType> HashNonChiral;

  for (auto& i : smi) {
    ROMOL_SPTR mol = ROMOL_SPTR(SmilesToMol(i));
    std::vector<boost::uint32_t> atomCodes(mol->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol->getNumBonds());

    fillAtomBondCodes(*mol, CF_ELEMENT | CF_CHARGE /*|CF_VALENCE*/
                                | CF_ATOM_AROMATIC,
                      &atomCodes, &bondCodes);

    //            fillAtomBondCodes(*mol, CF_ATOM_ALL &(~(CF_BOND_CHIRALITY |
    //            CF_ATOM_CHIRALITY | CF_ISOTOPE)), &atomCodes, &bondCodes);
    HashCodeType res = generateMoleculeHashCode(*mol, nullptr, nullptr,
                                                &atomCodes, &bondCodes);
    HashNonChiral.push_back(res);
    std::cout << res << " = " << encode(&res, sizeof(res)) << " | " << i
              << std::endl;
  }
  bool passed = true;
  for (size_t i = 0; i < HashNonChiral.size(); i++)
    for (size_t j = 0; j < HashNonChiral.size(); j++)
      if (i != j && HashNonChiral[i] != HashNonChiral[j]) passed = false;
  TEST_ASSERT(passed);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test21() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test21 CHIRALITY == BOND"
                       << std::endl;
  std::cout << "Hash size = " << 8 * sizeof(HashCodeType) << " bits.\n";
  const char* smi[] = {
      // equal non-chiral BOND hash
      "C/C=C/C", "CC=CC", "C/C=C\\C",
  };

  std::vector<HashCodeType> HashNonChiral;

  for (auto& i : smi) {
    ROMOL_SPTR mol = ROMOL_SPTR(SmilesToMol(i));
    std::vector<boost::uint32_t> atomCodes(mol->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol->getNumBonds());

    fillAtomBondCodes(*mol, CF_BOND_ALL & (~(CF_BOND_CHIRALITY)), &atomCodes,
                      &bondCodes);
    HashCodeType res = generateMoleculeHashCode(*mol, nullptr, nullptr,
                                                &atomCodes, &bondCodes);
    HashNonChiral.push_back(res);
    std::cout << res << " = " << encode(&res, sizeof(res)) << " | " << i
              << std::endl;
  }
  bool passed = true;
  for (size_t i = 0; i < HashNonChiral.size(); i++)
    for (size_t j = 0; j < HashNonChiral.size(); j++)
      if (i != j && HashNonChiral[i] != HashNonChiral[j]) passed = false;
  TEST_ASSERT(passed);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test3() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test3 CHIRALITY DIFF" << std::endl;
  const char* smi[] = {
      // different chiral hash
      "C[C@H](F)Cl", "C[C@@H](F)Cl", "CC(F)Cl", "[13CH3]C(F)Cl",

      "C[C@H]1CC[C@H](C)CC1", "C[C@H]1CC[C@@H](C)CC1", "CC1CCC(C)CC1",
  };

  std::vector<HashCodeType> HashChiral;

  for (auto& i : smi) {
    ROMOL_SPTR mol = ROMOL_SPTR(SmilesToMol(i));
    std::vector<boost::uint32_t> atomCodes(mol->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol->getNumBonds());

    fillAtomBondCodes(*mol, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType resC = generateMoleculeHashCode(*mol, nullptr, nullptr,
                                                 &atomCodes, &bondCodes);
    HashChiral.push_back(resC);

    std::cout << resC << " = " << encode(&resC, sizeof(resC)) << "  " << i
              << std::endl;
  }

  bool passed = true;
  for (size_t i = 0; i < HashChiral.size(); i++)
    for (size_t j = 0; j < HashChiral.size(); j++)
      if (i != j && HashChiral[i] == HashChiral[j]) passed = false;
  TEST_ASSERT(passed);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test3a() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test3a CHIRALITY EQUAL" << std::endl;
  {
    const char* smi[] = {
        "C[C@H](F)Cl", "C[C@@H](Cl)F",
    };
    ROMOL_SPTR mol1 = ROMOL_SPTR(SmilesToMol(smi[0]));
    ROMOL_SPTR mol2 = ROMOL_SPTR(SmilesToMol(smi[1]));

    std::vector<boost::uint32_t> atomCodes(mol1->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol2->getNumBonds());

    fillAtomBondCodes(*mol1, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash1 = generateMoleculeHashCode(*mol1, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    fillAtomBondCodes(*mol2, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash2 = generateMoleculeHashCode(*mol2, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    std::cout << hash1 << " " << hash2 << std::endl;
    TEST_ASSERT(hash1 == hash2);
  }
  {
    const char* smi[] = {
        "C[C@@H](F)Cl", "C[C@H](Cl)F",
    };
    ROMOL_SPTR mol1 = ROMOL_SPTR(SmilesToMol(smi[0]));
    ROMOL_SPTR mol2 = ROMOL_SPTR(SmilesToMol(smi[1]));

    std::vector<boost::uint32_t> atomCodes(mol1->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol2->getNumBonds());

    fillAtomBondCodes(*mol1, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash1 = generateMoleculeHashCode(*mol1, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    fillAtomBondCodes(*mol2, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash2 = generateMoleculeHashCode(*mol2, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    std::cout << hash1 << " " << hash2 << std::endl;
    TEST_ASSERT(hash1 == hash2);
  }

  {
    const char* smi[] = {
        "C/C=C/Cl", "Cl/C=C/C",
    };
    ROMOL_SPTR mol1 = ROMOL_SPTR(SmilesToMol(smi[0]));
    ROMOL_SPTR mol2 = ROMOL_SPTR(SmilesToMol(smi[1]));

    std::vector<boost::uint32_t> atomCodes(mol1->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol2->getNumBonds());

    fillAtomBondCodes(*mol1, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash1 = generateMoleculeHashCode(*mol1, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    fillAtomBondCodes(*mol2, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash2 = generateMoleculeHashCode(*mol2, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    std::cout << hash1 << " " << hash2 << std::endl;
    TEST_ASSERT(hash1 == hash2);
  }
  {
    const char* smi[] = {
        "C/C=C/Cl", "C/C=C\\Cl",
    };
    ROMOL_SPTR mol1 = ROMOL_SPTR(SmilesToMol(smi[0]));
    ROMOL_SPTR mol2 = ROMOL_SPTR(SmilesToMol(smi[1]));

    std::vector<boost::uint32_t> atomCodes(mol1->getNumAtoms());
    std::vector<boost::uint32_t> bondCodes(mol2->getNumBonds());

    fillAtomBondCodes(*mol1, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash1 = generateMoleculeHashCode(*mol1, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    fillAtomBondCodes(*mol2, CF_BOND_CHIRALITY | CF_ATOM_CHIRALITY | CF_ISOTOPE,
                      &atomCodes, &bondCodes);
    HashCodeType hash2 = generateMoleculeHashCode(*mol2, nullptr, nullptr,
                                                  &atomCodes, &bondCodes);
    std::cout << hash1 << " " << hash2 << std::endl;
    TEST_ASSERT(hash1 != hash2);
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test4() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test4 STRING" << std::endl;
  const char* smi[] = {
      // different chiral hash and equal non-chiral hash
      "C[C@H](F)Cl", "C[C@@H](F)Cl", "CC(F)Cl", "[13CH3]C(F)Cl",
      // different chiral hash
      "C[C@H]1CC[C@H](C)CC1", "C[C@H]1CC[C@@H](C)CC1", "CC1CCC(C)CC1",
  };

  for (auto& i : smi) {
    ROMOL_SPTR mol = ROMOL_SPTR(SmilesToMol(i));
    std::cout << generateMoleculeHashSet(*mol, nullptr, nullptr) << "  " << i
              << std::endl;
  }
  TEST_ASSERT(true);  // there is no any exseption
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test5() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing MolHash test5 " << std::endl;
  const char* smi[] = {// different chiral hash and equal non-chiral hash
                       // groups of 3
                       "C[CH](F)Cl", "C[C@H](F)Cl", "C[C@@H](F)Cl",
                       //
                       "c1cc(C[CH](F)Cl)cnc1", "c1cc(C[C@H](F)Cl)cnc1",
                       "c1cc(C[C@@H](F)Cl)cnc1"};

  for (size_t i = 0; i < sizeof(smi) / sizeof(smi[0]); i += 3) {
    ROMOL_SPTR mol1 = ROMOL_SPTR(SmilesToMol(smi[i]));
    TEST_ASSERT(mol1);
    ROMOL_SPTR mol2 = ROMOL_SPTR(SmilesToMol(smi[i + 1]));
    TEST_ASSERT(mol2);
    ROMOL_SPTR mol3 = ROMOL_SPTR(SmilesToMol(smi[i + 2]));
    TEST_ASSERT(mol3);
    {
      std::string hash1 = generateMoleculeHashSet(*mol1);
      std::string hash2 = generateMoleculeHashSet(*mol2);
      std::string hash3 = generateMoleculeHashSet(*mol3);
      TEST_ASSERT(hash1 != hash2);
      TEST_ASSERT(hash1 != hash3);
      TEST_ASSERT(hash3 != hash2);
    }
    // {
    //   std::string hash1=generateMoleculeHashSet(*mol1);
    //   std::string hash2=generateMoleculeHashSet(*mol2);
    //   std::cout << hash1 <<"  "<< smi[i] << std::endl;
    //   std::cout << hash2 <<"  "<< smi[i+1] << std::endl;
    //   TEST_ASSERT(hash1!=hash2);
    // }
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void doUnitTest() {
  std::cout << "Hash size = " << 8 * sizeof(HashCodeType) << " bits.\n";

  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test1();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test2();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test21();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test3();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test3a();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test4();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  test5();
}

//=============================================================================
// investigation test case for computing of a probability of the hash code
// collisions
//=============================================================================

std::string getSmilesOnly(const char* smiles, std::string* id = nullptr) {
  const char* sp = strchr(smiles, ' ');
  unsigned n = (sp ? sp - smiles + 1 : strlen(smiles));
  if (id) *id = std::string(smiles + n);
  return std::string(smiles, n);
}

HashCodeType computeHash(const ROMol& mol, CodeFlags flags) {
  std::vector<boost::uint32_t> atomCodes;
  std::vector<boost::uint32_t> bondCodes;

  fillAtomBondCodes(mol, flags, &atomCodes, &bondCodes);

  std::vector<unsigned> atomsToUse;
  std::vector<unsigned> bondsToUse;

  unsigned n = mol.getNumAtoms();
  for (unsigned i = 0; i < n; i++) {
    const Atom* atom = mol.getAtomWithIdx(i);
    if (1) atomsToUse.push_back(atom->getIdx());
  }

  n = mol.getNumBonds();
  for (unsigned i = 0; i < n; i++) {
    const Bond* bond = mol.getBondWithIdx(i);
    if (1) bondsToUse.push_back(bond->getIdx());
  }

  return generateMoleculeHashCode(mol, &atomsToUse, &bondsToUse, &atomCodes,
                                  &bondCodes);
}

//   {num atoms, num bonds} - {formula hash}
// - {non-chiral atom hashes} - {non-chiral bond hashes}
// - {chiral information}

#pragma pack(push, 1)
struct HashResult {
  unsigned Line;  // molecule Id [1, ...)
  HashCodeType Hash;
  //        unsigned   ChiralInfo;
  HashResult(unsigned id = 0)
      : Line(id),
        Hash(0)  //, ChiralInfo(0)
  {}
};
#pragma pack(pop)

bool lessHashResult_ALL(const HashResult& r, const HashResult& l) {
  return r.Hash < l.Hash;
}

void analyzeResults(std::list<HashResult>& res) {
  std::cerr << "\nANALYZING " << res.size() << " Results ...\n";
  std::cout << "Collisions found:\n";
  //        std::sort(res.begin(), res.end(), lessHashResult_ALL);
  unsigned rn = 0, cn = 0;
  for (auto r0 = res.begin(); r0 != res.end(); r0++) {
    std::cerr << "Result: " << ++rn << "\r";
    if (0 == r0->Line)  // collision has been already found
      continue;
    unsigned hashCollision = 0;
    std::vector<unsigned> cl;
    // use binary search of collision in sorted list to improve performance
    //........
    auto r1 = r0;
    for (auto r = ++r1; r != res.end(); r++) {
      if (0 == r->Line)  // collision has been already found
        continue;
      if (r->Hash == r0->Hash)  // collision found
      {
        ++hashCollision;
        cl.push_back(r->Line);
        //                    std::cout<<r0->Id<<"="<<r->Line<<"\n"; // TEMP
        //                    TEST
        r->Line =
            0;  // mark as already processed collision to exclude duplicates
      }
    }
    if (0 != hashCollision)  // collision found
    {
      cn += hashCollision;
      std::cout << "mol line " << r0->Line << ": " << hashCollision
                << " collisions with: ";
      for (unsigned int i : cl) std::cout << i << " ";
      std::cout << "lines.\n";
    }
  }
  std::cout << "Total: " << cn << " hash collisions found in " << res.size()
            << " molecules.\n";
}

void testFileSMILES(const char* file, HashCodeType bitMask) {
  unsigned line = 0;
  std::list<HashResult> res;
  std::cout << "FILE: " << file << "\n";

  FILE* f = fopen(file, "rt");
  if (!f) {
    perror("Could not OPEN smi file");
    return;
  }
  char smiles[4096];
  while (fgets(smiles, sizeof(smiles), f) && line <= 1000999) {
    for (size_t i = strlen(smiles) - 1; i > 0 && smiles[i] < ' '; i--)
      smiles[i] = '\0';  // remove LF
    std::string id;
    std::cerr << "\rLine: " << ++line << " ";
    if ('#' != smiles[0] && ' ' != smiles[0] &&
        '/' != smiles[0]                    // commented to skip
        && nullptr == strchr(smiles, '.'))  // skip ions
    {
      ROMOL_SPTR mol;
      try {
        mol = ROMOL_SPTR(SmilesToMol(getSmilesOnly(smiles, &id)));
      } catch (...)  // internal RDKit error: Invar::Invariant& ex
      {
        std::cerr << " RDKit error: " << smiles << "/n";
        continue;
      }
      res.push_back(HashResult(line));
      HashResult& r = res.back();
      //                r.ChiralInfo = 0;//mol-();
      r.Hash = computeHash(*mol, CF_ALL) & bitMask;
    } else
      std::cerr << " skipped: " << smiles << "/n";
  }
  fclose(f);
  std::cout << "\nDONE. " << res.size() << " molecules processed.\n";
  analyzeResults(res);
  std::cout << "Test COMPLETED.\n";
}

void checkCollisions(const char* file, boost::uint32_t bits = 0) {
  HashCodeType bitMask = 0;
  if (0 == bits || 8 * sizeof(HashCodeType) < bits)
    bits = 8 * sizeof(HashCodeType);
  for (unsigned i = 0; i < bits; i++) bitMask |= 1ULL << i;
  std::cout << "Hash size = " << bits << " bits. Mask = " << bitMask << "\n";

  if (0 == strcmp(file + strlen(file) - 4, ".smi"))
    testFileSMILES(file, bitMask);
  else
    std::cout << "UNKNOWN File Extention.\n";
}

}  // RDKit

int main(int argc, char* argv[]) {
  RDKit::doUnitTest();

  if (2 == argc)
    RDKit::checkCollisions(argv[1]);
  else if (3 == argc && isdigit(*argv[2]))
    RDKit::checkCollisions(argv[1], atoi(argv[2]));
  else if (1 != argc)
    std::cout << "UNKNOWN Argument.\n";
  return 0;
}
