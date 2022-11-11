//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

std::string smis[] = {
    "C[C@@H]1CCC[C@H](C)[C@H]1C",
    "N[C@@]1(C[C@H]([18F])C1)C(=O)O",
    "CC(C)CCCC[C@@H]1C[C@H](/C=C/[C@]2(C)CC[C@H](O)CC2)[C@@H](O)[C@H]1O",
    "COC(=O)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(\\C)/C=C/C(=O)[O-]",
    "[O:1]=[C:2]([CH2:3][C:4]1=[CH:5][CH:6]=[CH:7][CH:8]=[CH:9]1)[NH2:10]",
    "Cl[C@H]1[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H]1Cl",
    "O=S(=O)(NC[C@H]1CC[C@H](CNCc2ccc3ccccc3c2)CC1)c1ccc2ccccc2c1",
    "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)-c1cccc(-"
    "c2ccc(cc2)C[n+]2ccn(c2)Cc2ccc(cc2)C3=O)c1C(=O)O.[Br-].[Br-]",
    "CCCCCCC1C23C4=c5c6c7c8c9c%10c%11c%12c%13c%14c%15c%16c%17c%18c%19c%20c%21c%"
    "22c%23c(c5c5c6c6c8c%11c8c%11c%12c%15c%12c(c%20%16)c%21c%15c%23c5c(c68)c%"
    "15c%12%11)C2(C[N+]1(C)C)C%22C%19c1c-%18c2c5c(c13)C4C7C9=C5C1(C2C%17%14)C("
    "CCCCCC)[N+](C)(C)CC%10%131.[I-].[I-]",
    "C12C3C4C5C6C7C8C1C1C9C5C5C%10C2C2C%11C%12C%13C3C3C7C%10C7C4C%11C1C3C(C5C8%"
    "12)C(C62)C7C9%13"};

void testAtomPairFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test atom-pair fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp;
  std::uint32_t c1, c2, c3;

  FingerprintGenerator<std::uint32_t> *atomPairGenerator =
      AtomPair::getAtomPairGenerator<std::uint32_t>();

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 1);

  delete mol;
  delete fp;
  mol = SmilesToMol("CC=O.Cl");
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 3);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 1);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 1);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c2, c3, 1, false)) == 1);

  delete mol;
  delete fp;
  delete atomPairGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testAtomPairArgs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "Test atom-pair fp generator with different arguments" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp;
  std::uint32_t c1, c2, c3;

  FingerprintGenerator<std::uint32_t> *atomPairGenerator =
      AtomPair::getAtomPairGenerator<std::uint32_t>(2);

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 1);
  TEST_ASSERT(fp->getNonzeroElements().size() == 1);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 0);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 1);

  delete fp;
  delete atomPairGenerator;

  atomPairGenerator = AtomPair::getAtomPairGenerator<std::uint32_t>(1, 1);
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 2);
  TEST_ASSERT(fp->getNonzeroElements().size() == 1);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 0);

  delete fp;
  delete atomPairGenerator;

  atomPairGenerator =
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, 30, true);
  atomPairGenerator->getOptions()->df_countSimulation = false;
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, true);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, true);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, true);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, true)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, true)) == 1);

  delete fp;

  atomPairGenerator->getOptions()->df_includeChirality = false;
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 1);

  delete mol;
  delete fp;
  delete atomPairGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testAtomPairOld() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test compatibility between atom pair "
                           "implementation for FingerprintGenerator"
                           " and old atom pairs implementation"
                        << std::endl;
  {
    ROMol *mol;
    SparseIntVect<std::int32_t> *fp1, *fp2;
    SparseIntVect<std::uint32_t> *fpu;

    FingerprintGenerator<std::uint32_t> *atomPairGenerator =
        AtomPair::getAtomPairGenerator<std::uint32_t>();

    for (const auto &sm : smis) {
      mol = SmilesToMol(sm);
      fp1 = AtomPairs::getAtomPairFingerprint(*mol);
      fpu = atomPairGenerator->getSparseCountFingerprint(*mol);
      fp2 = new SparseIntVect<std::int32_t>(fpu->size());
      std::map<std::uint32_t, int> nz = fpu->getNonzeroElements();
      for (auto &it : nz) {
        fp2->setVal(static_cast<std::int32_t>(it.first), it.second);
      }

      TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
      TEST_ASSERT(*fp1 == *fp2);

      delete mol;
      delete fp1;
      delete fp2;
      delete fpu;
    }

    delete atomPairGenerator;

    BOOST_LOG(rdErrorLog) << "  done" << std::endl;
  }
}

// todo this test needs to be updated since the fingerprint size logic is
// changed, count simulation no longer makes fingerprints larger
void testAtomPairNonSparseBitvector() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test consistency between different output types "
                           "for folded atom pair fingerprint"
                        << std::endl;
  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp1;
    ExplicitBitVect *fp2;

    FingerprintGenerator<std::uint32_t> *atomPairGenerator =
        AtomPair::getAtomPairGenerator<std::uint32_t>();
    std::vector<std::uint32_t> defaultCountBounds = {1, 2, 4, 8};

    mol = SmilesToMol("CCC");
    fp1 = atomPairGenerator->getCountFingerprint(*mol);
    fp2 = atomPairGenerator->getFingerprint(*mol);

    std::map<std::uint32_t, int> nz = fp1->getNonzeroElements();
    for (const auto &it : nz) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); ++i) {
        bool isSet = static_cast<bool>(
            fp2->getBit(it.first * defaultCountBounds.size() + i));
        TEST_ASSERT(isSet ==
                    (it.second >= static_cast<long>(defaultCountBounds[i])));
      }
    }

    delete mol;
    delete fp1;
    delete fp2;

    mol = SmilesToMol("CC=O.Cl");
    fp1 = atomPairGenerator->getCountFingerprint(*mol);
    fp2 = atomPairGenerator->getFingerprint(*mol);

    nz = fp1->getNonzeroElements();
    for (const auto &it : nz) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); ++i) {
        bool isSet = static_cast<bool>(
            fp2->getBit(it.first * defaultCountBounds.size() + i));
        TEST_ASSERT(isSet ==
                    (it.second >= static_cast<long>(defaultCountBounds[i])));
      }
    }

    delete mol;
    delete fp1;
    delete fp2;
    delete atomPairGenerator;

    BOOST_LOG(rdErrorLog) << "  done" << std::endl;
  }
}

void testAtomPairOutput() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test atom-pair additional output" << std::endl;
  auto mol = "CCC"_smiles;
  TEST_ASSERT(mol);
  AdditionalOutput additionalOutput;

  auto atomPairGenerator = AtomPair::getAtomPairGenerator<std::uint32_t>();
  additionalOutput.allocateAtomToBits();
  additionalOutput.allocateAtomCounts();

  auto fp = atomPairGenerator->getSparseCountFingerprint(*mol, nullptr, nullptr,
                                                         -1, &additionalOutput);
  auto c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  auto c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  auto c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  auto v = additionalOutput.atomToBits->at(0);
  TEST_ASSERT(std::find(v.begin(), v.end(),
                        (AtomPair::getAtomPairCode(c1, c2, 1, false))) !=
              v.end());
  TEST_ASSERT(std::find(v.begin(), v.end(),
                        (AtomPair::getAtomPairCode(c1, c3, 2, false))) !=
              v.end());
  v = additionalOutput.atomToBits->at(1);
  TEST_ASSERT(std::find(v.begin(), v.end(),
                        (AtomPair::getAtomPairCode(c1, c2, 1, false))) !=
              v.end());
  TEST_ASSERT(std::find(v.begin(), v.end(),
                        (AtomPair::getAtomPairCode(c2, c3, 1, false))) !=
              v.end());
  v = additionalOutput.atomToBits->at(2);
  TEST_ASSERT(std::find(v.begin(), v.end(),
                        (AtomPair::getAtomPairCode(c1, c3, 2, false))) !=
              v.end());
  TEST_ASSERT(std::find(v.begin(), v.end(),
                        (AtomPair::getAtomPairCode(c2, c3, 1, false))) !=
              v.end());
  TEST_ASSERT(additionalOutput.atomCounts->at(0) == 2);
  TEST_ASSERT(additionalOutput.atomCounts->at(1) == 2);
  TEST_ASSERT(additionalOutput.atomCounts->at(2) == 2);

  delete fp;
  delete atomPairGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

// Old version test name testHashedAtomPairs
void testFoldedAtomPairs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Folded Atom Pairs." << std::endl;

  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>(1, 30, false, true,
                                                      nullptr, true, 1024);
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getCountFingerprint(*mol);
    SparseIntVect<std::uint32_t> *fp2;
    fp2 = generator->getCountFingerprint(*mol);
    TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
    TEST_ASSERT(*fp1 == *fp2);

    delete mol;
    delete fp2;
    mol = SmilesToMol("c1ccccn1");
    fp2 = generator->getCountFingerprint(*mol);
    RANGE_CHECK(0.0, DiceSimilarity(*fp1, *fp2), 1.0);

    delete mol;
    delete fp1;
    delete fp2;
    delete generator;
  }

  {
    FingerprintGenerator<std::uint32_t> *generator1 =
        AtomPair::getAtomPairGenerator<std::uint32_t>();
    FingerprintGenerator<std::uint32_t> *generator2 =
        AtomPair::getAtomPairGenerator<std::uint32_t>(1, 2);
    FingerprintGenerator<std::uint32_t> *generator3 =
        AtomPair::getAtomPairGenerator<std::uint32_t>(1, 3);
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator1->getCountFingerprint(*mol);
    SparseIntVect<std::uint32_t> *fp2;
    fp2 = generator3->getCountFingerprint(*mol);
    TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
    TEST_ASSERT(*fp1 == *fp2);

    delete fp2;
    fp2 = generator2->getCountFingerprint(*mol);
    RANGE_CHECK(0.0, DiceSimilarity(*fp1, *fp2), 1.0);

    delete mol;
    delete fp1;
    delete fp2;
    delete generator1;
    delete generator2;
    delete generator3;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRootedAtomPairs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Rooted Atom Pairs." << std::endl;

  FingerprintGenerator<std::uint32_t> *generator =
      AtomPair::getAtomPairGenerator<std::uint32_t>();
  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp1, *fp2;
  std::vector<std::uint32_t> roots;

  mol = SmilesToMol("OCCCCC");
  fp1 = generator->getSparseCountFingerprint(*mol);
  SparseIntVect<std::uint32_t>::StorageType nz1 = fp1->getNonzeroElements();
  TEST_ASSERT(nz1.size() > 0);

  roots.push_back(0);
  fp2 = generator->getSparseCountFingerprint(*mol, &roots);
  SparseIntVect<std::uint32_t>::StorageType nz2 = fp2->getNonzeroElements();
  TEST_ASSERT(nz2.size() > 0);
  TEST_ASSERT(nz2.size() < nz1.size());

  for (SparseIntVect<std::uint32_t>::StorageType::const_iterator bIt =
           nz2.begin();
       bIt != nz2.end(); ++bIt) {
    TEST_ASSERT(bIt->second <= fp2->getVal(bIt->first));
  }

  delete mol;
  delete fp1;
  delete fp2;
  delete generator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIgnoreAtomPairs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test ignoring atoms in Atom Pairs."
                        << std::endl;

  FingerprintGenerator<std::uint32_t> *generator =
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, 5, false, true, nullptr,
                                                    true, 4096);

  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp1, *fp2;
    std::vector<std::uint32_t> roots;

    mol = SmilesToMol("OCCCCC");
    fp1 = generator->getSparseCountFingerprint(*mol);
    SparseIntVect<std::uint32_t>::StorageType nz1 = fp1->getNonzeroElements();
    TEST_ASSERT(nz1.size() > 0);

    roots.push_back(0);
    fp2 = generator->getSparseCountFingerprint(*mol, nullptr, &roots);
    SparseIntVect<std::uint32_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() == nz1.size() - 5);

    for (SparseIntVect<std::uint32_t>::StorageType::const_iterator bIt =
             nz2.begin();
         bIt != nz2.end(); ++bIt) {
      TEST_ASSERT(bIt->second <= fp2->getVal(bIt->first));
    }

    delete mol;
    delete fp1;
    delete fp2;
  }
  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp2;
    std::vector<std::uint32_t> roots;

    mol = SmilesToMol("OCCCCC");
    roots.push_back(0);
    fp2 = generator->getSparseCountFingerprint(*mol, &roots, &roots);
    SparseIntVect<std::uint32_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() == 0);

    delete mol;
    delete fp2;
  }

  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp1, *fp2;
    std::vector<std::uint32_t> roots;

    mol = SmilesToMol("OCCCCC");
    fp1 = generator->getCountFingerprint(*mol);
    SparseIntVect<std::uint32_t>::StorageType nz1 = fp1->getNonzeroElements();
    TEST_ASSERT(nz1.size() > 0);

    roots.push_back(0);
    fp2 = generator->getCountFingerprint(*mol, nullptr, &roots);
    SparseIntVect<std::uint32_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() < nz1.size());

    for (SparseIntVect<std::uint32_t>::StorageType::const_iterator bIt =
             nz2.begin();
         bIt != nz2.end(); ++bIt) {
      TEST_ASSERT(bIt->second <= fp2->getVal(bIt->first));
    }

    delete mol;
    delete fp1;
    delete fp2;
  }
  delete generator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testChiralPairs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Atom Pairs including info about chirality."
                        << std::endl;

  FingerprintGenerator<std::uint32_t> *generator =
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, 5, false, true, nullptr,
                                                    true, 4096);
  FingerprintGenerator<std::uint32_t> *generatorChirality =
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, 5, true, true, nullptr,
                                                    true, 4096);
  ROMol *m1, *m2, *m3;

  m1 = SmilesToMol("CC[CH](F)Cl");
  TEST_ASSERT(m1);
  m2 = SmilesToMol("CC[C@H](F)Cl");
  TEST_ASSERT(m1);
  m3 = SmilesToMol("CC[C@@H](F)Cl");
  TEST_ASSERT(m1);

  {
    SparseIntVect<std::uint32_t> *fp1, *fp2, *fp3;
    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 10);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 10);
    fp2 = generator->getSparseCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 10);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 10);
    fp3 = generator->getSparseCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 10);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 10);

    TEST_ASSERT((*fp1) == (*fp2));
    TEST_ASSERT((*fp1) == (*fp3));
    TEST_ASSERT((*fp2) == (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;

    fp1 = generatorChirality->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 10);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 10);
    fp2 = generatorChirality->getSparseCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 10);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 10);
    fp3 = generatorChirality->getSparseCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 10);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 10);

    TEST_ASSERT((*fp1) != (*fp2));
    TEST_ASSERT((*fp1) != (*fp3));
    TEST_ASSERT((*fp2) != (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;
  }

  {
    SparseIntVect<std::uint32_t> *fp1, *fp2, *fp3;
    fp1 = generator->getCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 10);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 10);
    fp2 = generator->getCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 10);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 10);
    fp3 = generator->getCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 10);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 10);

    TEST_ASSERT((*fp1) == (*fp2));
    TEST_ASSERT((*fp1) == (*fp3));
    TEST_ASSERT((*fp2) == (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;

    fp1 = generatorChirality->getCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 10);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 10);
    fp2 = generatorChirality->getCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 10);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 10);
    fp3 = generatorChirality->getCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 10);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 10);

    TEST_ASSERT((*fp1) != (*fp2));
    TEST_ASSERT((*fp1) != (*fp3));
    TEST_ASSERT((*fp2) != (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;
  }

  delete m1;
  delete m2;
  delete m3;
  delete generator;
  delete generatorChirality;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFPOld() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test morgan fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp, *fpOld;
  int radius = 3;

  FingerprintGenerator<std::uint32_t> *morganGenerator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(radius);

  for (const auto &sm : smis) {
    mol = SmilesToMol(sm);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    fpOld = MorganFingerprints::getFingerprint(*mol, radius);
    TEST_ASSERT(DiceSimilarity(*fp, *fpOld) == 1.0);
    TEST_ASSERT(*fp == *fpOld);

    delete mol;
    delete fp;
    delete fpOld;
  }

  delete morganGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Morgan Fingerprints." << std::endl;

  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("CCCCC");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    fp = morganGenerator->getCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 5);
    delete fp;
    fp = morganGenerator->getCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 5);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 7);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 7);
    delete fp;
    delete morganGenerator;

    delete mol;
  }
  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("O=C(O)CC1CC1");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 6);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = MorganFingerprints::getFingerprint(*mol, 1);
    TEST_ASSERT(fp->getNonzeroElements().size() == 12);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 16);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 17);
    delete fp;
    delete morganGenerator;

    delete mol;
  }

  {
    // test that the results aren't order dependent, i.e. that we're
    // "canonicalizing" the fps correctly
    ROMol *mol, *mol2;
    SparseIntVect<std::uint32_t> *fp, *fp2;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("O=C(O)CC1CC1");
    mol2 = SmilesToMol("OC(=O)CC1CC1");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    fp2 = morganGenerator->getSparseCountFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 6);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 6);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    fp2 = morganGenerator->getSparseCountFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 12);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 12);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    fp2 = morganGenerator->getSparseCountFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 16);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 16);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    fp2 = morganGenerator->getSparseCountFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 17);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 17);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    delete mol;
    delete mol2;
  }

  {
    // symmetry test:
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    mol = SmilesToMol("OCCCCO");
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 7);
    SparseIntVect<std::uint32_t>::StorageType::const_iterator iter;
    for (iter = fp->getNonzeroElements().begin();
         iter != fp->getNonzeroElements().end(); ++iter) {
      TEST_ASSERT(iter->second == 2 || iter->second == 4);
    }

    delete fp;
    delete mol;
    delete morganGenerator;
  }

  {
    // chirality test:
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("CC(F)(Cl)C(F)(Cl)C");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 8);
    delete mol;
    delete fp;
    delete morganGenerator;

    mol = SmilesToMol("CC(F)(Cl)[C@](F)(Cl)C");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    delete fp;
    delete morganGenerator;
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 8);
    delete fp;
    delete morganGenerator;

    morganGenerator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(0, true, true);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    delete fp;
    delete morganGenerator;

    morganGenerator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(1, true, true);
    fp = morganGenerator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 9);
    delete fp;
    delete morganGenerator;
    delete mol;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFPFromAtoms() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Morgan Fingerprints using fromAtoms argument." << std::endl;

  FingerprintGenerator<std::uint32_t> *radius0Generator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
  FingerprintGenerator<std::uint32_t> *radius1Generator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
  FingerprintGenerator<std::uint32_t> *radius2Generator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
  FingerprintGenerator<std::uint32_t> *radius3Generator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    std::vector<std::uint32_t> atoms;
    atoms.push_back(0);

    mol = SmilesToMol("CCCCC");
    fp = radius0Generator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    fp = radius0Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 1);
    delete fp;

    fp = radius1Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    // tests issue 3415636
    fp = radius2Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 3);
    delete fp;

    delete mol;
  }

  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    std::vector<std::uint32_t> atoms;

    mol = SmilesToMol("CCCCC");
    fp = radius0Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 0);
    delete fp;

    fp = radius1Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 0);
    delete fp;

    delete mol;
  }

  {
    ROMol *mol;
    SparseIntVect<std::uint32_t> *fp;
    std::vector<std::uint32_t> atoms;
    atoms.push_back(0);

    mol = SmilesToMol("C(CC)CO");

    fp = radius0Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 1);
    delete fp;

    fp = radius1Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    fp = radius2Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 3);
    delete fp;

    // tests issue 3415636
    fp = radius3Generator->getSparseCountFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 3);
    delete fp;

    delete mol;
  }

  delete radius0Generator;
  delete radius1Generator;
  delete radius2Generator;
  delete radius3Generator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFPBitVect() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Morgan Fingerprints as bit vects."
                        << std::endl;

  {
    ROMol *mol;
    ExplicitBitVect *fp;
    FingerprintGenerator<std::uint32_t> *radius0Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(0, false);
    FingerprintGenerator<std::uint32_t> *radius1Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(1, false);
    FingerprintGenerator<std::uint32_t> *radius2Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2, false);
    FingerprintGenerator<std::uint32_t> *radius3Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(3, false);

    mol = SmilesToMol("CCCCC");
    fp = radius0Generator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 2);
    delete fp;
    fp = radius1Generator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 5);
    delete fp;
    fp = radius2Generator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 7);
    delete fp;
    fp = radius3Generator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 7);
    delete fp;

    delete mol;
    delete radius0Generator;
    delete radius1Generator;
    delete radius2Generator;
    delete radius3Generator;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFPFeatureInvs() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Morgan Fingerprints with feature invariants." << std::endl;

  {
    ROMol *mol;
    mol = SmilesToMol("Cc1ccccc1");
    TEST_ASSERT(mol);
    auto *invGen = new MorganFingerprint::MorganFeatureAtomInvGenerator();
    std::vector<std::uint32_t> *invars = invGen->getAtomInvariants(*mol);
    TEST_ASSERT((*invars)[0] == 0);
    TEST_ASSERT((*invars)[1] != 0);
    TEST_ASSERT((*invars)[1] == (*invars)[2]);
    TEST_ASSERT((*invars)[1] == (*invars)[3]);
    TEST_ASSERT((*invars)[1] == (*invars)[4]);
    TEST_ASSERT((*invars)[1] == (*invars)[5]);
    TEST_ASSERT((*invars)[1] == (*invars)[6]);
    delete mol;
    delete invGen;
    delete invars;
  }
  {
    ROMol *mol;
    mol = SmilesToMol("FCCCl");
    TEST_ASSERT(mol);
    auto *invGen = new MorganFingerprint::MorganFeatureAtomInvGenerator();
    std::vector<std::uint32_t> *invars = invGen->getAtomInvariants(*mol);
    TEST_ASSERT((*invars)[1] == (*invars)[2]);
    TEST_ASSERT((*invars)[1] == 0);
    TEST_ASSERT((*invars)[0] == (*invars)[3]);
    TEST_ASSERT((*invars)[0] != 0);
    delete mol;
    delete invGen;
    delete invars;
  }

  {
    ROMol *mol;
    mol = SmilesToMol("Cc1ncccc1O");
    TEST_ASSERT(mol);
    std::vector<const ROMol *> patterns(2);

    RWMol *p;
    p = SmartsToMol("[A]");
    patterns[0] = static_cast<const ROMol *>(p);
    p = SmartsToMol("[a]");
    patterns[1] = static_cast<const ROMol *>(p);
    auto *invGen =
        new MorganFingerprint::MorganFeatureAtomInvGenerator(&patterns);

    std::vector<std::uint32_t> *invars = invGen->getAtomInvariants(*mol);
    TEST_ASSERT((*invars)[0] != 0);
    TEST_ASSERT((*invars)[1] != 0);
    TEST_ASSERT((*invars)[0] != (*invars)[1]);
    TEST_ASSERT((*invars)[1] == (*invars)[2]);
    TEST_ASSERT((*invars)[0] == (*invars)[7]);
    delete mol;
    delete invGen;
    delete invars;
    delete patterns[0];
    delete patterns[1];
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFPOptions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test additional Morgan fingerprints options."
                        << std::endl;

  {
    ROMol *m1, *m2;
    ExplicitBitVect *fp1, *fp2;
    FingerprintGenerator<std::uint32_t> *generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    FingerprintGenerator<std::uint32_t> *generatorNoBoundType =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2, true, false,
                                                             false);
    std::vector<std::uint32_t> invars(3);
    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;

    m1 = SmilesToMol("CCC");
    TEST_ASSERT(m1);
    m2 = SmilesToMol("CC=C");
    TEST_ASSERT(m2);

    fp1 =
        generator->getFingerprint(*m1, nullptr, nullptr, -1, nullptr, &invars);
    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;
    fp2 =
        generator->getFingerprint(*m2, nullptr, nullptr, -1, nullptr, &invars);
    TEST_ASSERT((*fp1) != (*fp2));
    delete fp1;
    delete fp2;

    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;
    fp1 = generatorNoBoundType->getFingerprint(*m1, nullptr, nullptr, -1,
                                               nullptr, &invars);
    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;
    fp2 = generatorNoBoundType->getFingerprint(*m2, nullptr, nullptr, -1,
                                               nullptr, &invars);
    TEST_ASSERT((*fp1) == (*fp2));
    delete fp1;
    delete fp2;

    delete m1;
    delete m2;
    delete generator;
    delete generatorNoBoundType;
  }

  {
    ROMol *m1, *m2, *m3;
    ExplicitBitVect *fp1, *fp2, *fp3;
    FingerprintGenerator<std::uint32_t> *generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    FingerprintGenerator<std::uint32_t> *generatorChirality =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2, true, true);

    m1 = SmilesToMol("C[C@H](F)Cl");
    TEST_ASSERT(m1);
    m2 = SmilesToMol("C[C@@H](F)Cl");
    TEST_ASSERT(m2);
    m3 = SmilesToMol("CC(F)Cl");
    TEST_ASSERT(m3);

    fp1 = generator->getFingerprint(*m1);
    fp2 = generator->getFingerprint(*m2);
    fp3 = generator->getFingerprint(*m3);
    TEST_ASSERT((*fp1) == (*fp2));
    TEST_ASSERT((*fp1) == (*fp3));
    TEST_ASSERT((*fp2) == (*fp3));
    delete fp1;
    delete fp2;
    delete fp3;

    fp1 = generatorChirality->getFingerprint(*m1);
    fp2 = generatorChirality->getFingerprint(*m2);
    fp3 = generatorChirality->getFingerprint(*m3);
    TEST_ASSERT((*fp1) != (*fp2));
    TEST_ASSERT((*fp1) != (*fp3));
    TEST_ASSERT((*fp2) != (*fp3));
    delete fp1;
    delete fp2;
    delete fp3;

    delete m1;
    delete m2;
    delete m3;
    delete generator;
    delete generatorChirality;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue695() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test GitHub Issue 695: Include cis/trans "
                           "stereochemistry when useChirality=true with the "
                           "morgan fingerprints"
                        << std::endl;

  FingerprintGenerator<std::uint32_t> *generator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
  FingerprintGenerator<std::uint32_t> *generatorChirality =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(1, true, true);
  {
    ROMol *m1 = SmilesToMol("CC=CC");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp;
    SparseIntVect<std::uint32_t>::StorageType::const_iterator iter;

    fp = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    iter = fp->getNonzeroElements().find(736731344);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246703798);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246728737);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(3545353036);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    delete fp;

    fp = generatorChirality->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    iter = fp->getNonzeroElements().find(736731344);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246703798);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246728737);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(3545353036);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);

    delete fp;

    delete m1;
  }
  {
    ROMol *m1 = SmilesToMol("C/C=C/C");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp;
    SparseIntVect<std::uint32_t>::StorageType::const_iterator iter;

    fp = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    iter = fp->getNonzeroElements().find(736731344);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246703798);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246728737);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(3545353036);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    delete fp;

    fp = generatorChirality->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);

    iter = fp->getNonzeroElements().find(736735794);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246703798);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246728737);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(3545353036);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    delete fp;

    delete m1;
  }
  {
    ROMol *m1 = SmilesToMol("C/C=C\\C");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp;
    SparseIntVect<std::uint32_t>::StorageType::const_iterator iter;

    fp = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    iter = fp->getNonzeroElements().find(736731344);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246703798);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246728737);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(3545353036);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    delete fp;

    fp = generatorChirality->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    iter = fp->getNonzeroElements().find(736735858);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246703798);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(2246728737);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    iter = fp->getNonzeroElements().find(3545353036);
    TEST_ASSERT(iter != fp->getNonzeroElements().end() && iter->second == 2);
    delete fp;

    delete m1;
  }
  delete generator;
  delete generatorChirality;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue874() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Single atoms setting radius 1 bits in Morgan fingerprints  "
      << std::endl;
  {
    std::string smiles = "Cl";
    ROMol *m1 = SmilesToMol(smiles);
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumAtoms() == 1);
    FingerprintGenerator<std::uint32_t> *radius0Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    FingerprintGenerator<std::uint32_t> *radius1Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(1);

    SparseIntVect<std::uint32_t> *fp;
    fp = radius0Generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp->getNonzeroElements().size() == 1);
    delete fp;
    fp = radius1Generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp->getNonzeroElements().size() == 1);
    delete fp;

    delete m1;
    delete radius0Generator;
    delete radius1Generator;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testInvariantGenerators() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test different invariant generator combinations"
                        << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp;
  int radius = 3;

  AtomInvariantsGenerator *atomInvariantsGenerator =
      new MorganFingerprint::MorganAtomInvGenerator();
  BondInvariantsGenerator *bondInvariantsGenerator =
      new MorganFingerprint::MorganBondInvGenerator();

  FingerprintGenerator<std::uint32_t> *morganGenerator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(
          radius, true, false, true, false, atomInvariantsGenerator,
          bondInvariantsGenerator);

  mol = SmilesToMol("CCCCC");
  fp = morganGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getNonzeroElements().size() == 7);

  delete mol;
  delete fp;
  delete morganGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator = new MorganFingerprint::MorganAtomInvGenerator();
  bondInvariantsGenerator = new MorganFingerprint::MorganBondInvGenerator();

  morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(
      radius, true, false, true, false, atomInvariantsGenerator,
      bondInvariantsGenerator);

  mol = SmilesToMol("CCCCC");
  fp = morganGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getNonzeroElements().size() == 7);

  delete mol;
  delete fp;
  delete morganGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator =
      new MorganFingerprint::MorganFeatureAtomInvGenerator();
  bondInvariantsGenerator = new MorganFingerprint::MorganBondInvGenerator();

  morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(
      radius, true, false, true, false, atomInvariantsGenerator,
      bondInvariantsGenerator);

  mol = SmilesToMol("CCCCC");
  fp = morganGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getNonzeroElements().size() == 5);

  delete mol;
  delete fp;
  delete morganGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator = new AtomPair::AtomPairAtomInvGenerator();
  bondInvariantsGenerator = new MorganFingerprint::MorganBondInvGenerator();

  morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(
      radius, true, false, true, false, atomInvariantsGenerator,
      bondInvariantsGenerator);

  mol = SmilesToMol("CCCCC");
  fp = morganGenerator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getNonzeroElements().size() == 7);

  delete mol;
  delete fp;
  delete morganGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator = new MorganFingerprint::MorganAtomInvGenerator();
  bondInvariantsGenerator = nullptr;

  FingerprintGenerator<std::uint32_t> *atomPairGenerator =
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, radius, false, true,
                                                    atomInvariantsGenerator,
                                                    bondInvariantsGenerator);

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);

  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  delete mol;
  delete fp;
  delete atomPairGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator =
      new MorganFingerprint::MorganFeatureAtomInvGenerator();
  bondInvariantsGenerator = nullptr;

  atomPairGenerator = AtomPair::getAtomPairGenerator<std::uint32_t>(
      1, radius, false, true, atomInvariantsGenerator, bondInvariantsGenerator);

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getSparseCountFingerprint(*mol);

  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  delete mol;
  delete fp;
  delete atomPairGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testCustomInvariants() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test invariant overriding" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp1, *fp2;
  int radius = 3;

  AtomInvariantsGenerator *atomInvariantsGenerator =
      new AtomPair::AtomPairAtomInvGenerator();

  FingerprintGenerator<std::uint32_t> *morganGenerator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(
          radius, true, false, true, false, atomInvariantsGenerator);

  FingerprintGenerator<std::uint32_t> *defaultMorganGenerator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(radius);

  mol = SmilesToMol("CCCCC");
  std::vector<std::uint32_t> *customInvariants =
      atomInvariantsGenerator->getAtomInvariants(*mol);
  fp1 = morganGenerator->getSparseCountFingerprint(*mol);
  fp2 = defaultMorganGenerator->getSparseCountFingerprint(
      *mol, nullptr, nullptr, -1, nullptr, customInvariants);
  TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
  TEST_ASSERT(*fp1 == *fp2);

  delete mol;
  delete fp1;
  delete fp2;
  delete morganGenerator;
  delete defaultMorganGenerator;
  delete atomInvariantsGenerator;
  delete customInvariants;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRDKitFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test RDKit fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint64_t> *fp, *fpTemp, *fpOld;

  // I won't lie: having to do this makes my head hurt, but fixing it to create
  // a ctor that takes a Parameters object is more effort than I can devote at
  // the moment
  unsigned int minPath = 1;
  unsigned int maxPath = 7;
  bool useHs = true;
  bool branchedPaths = true;
  bool useBondOrder = true;
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
  bool countSimulation = false;
  const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  std::uint32_t fpSize = 2048;
  std::uint32_t numBitsPerFeature = 1;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGenerator(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>(
          minPath, maxPath, useHs, branchedPaths, useBondOrder,
          atomInvariantsGenerator, countSimulation, countBounds, fpSize,
          numBitsPerFeature));

  for (const auto &sm : smis) {
    mol = SmilesToMol(sm);
    fp = fpGenerator->getSparseCountFingerprint(*mol);
    fpTemp = getUnfoldedRDKFingerprintMol(*mol);

    // Old and new versions produce different length results, but data should be
    // the same
    std::map<std::uint64_t, int> nz = fpTemp->getNonzeroElements();
    fpOld = new SparseIntVect<std::uint64_t>(fp->getLength());
    for (auto &it : nz) {
      fpOld->setVal(it.first, it.second);
    }

    TEST_ASSERT(DiceSimilarity(*fp, *fpOld) == 1.0);
    TEST_ASSERT(*fp == *fpOld);

    delete mol;
    delete fp;
    delete fpOld;
    delete fpTemp;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRDKFPUnfolded() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test unfolded version of RDKFP   " << std::endl;
  {
    FingerprintGenerator<std::uint64_t> *generator =
        RDKitFP::getRDKitFPGenerator<std::uint64_t>();
    ROMol *m1 = SmilesToMol("c1ccccc1N");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint64_t> *fp1;
    SparseIntVect<std::uint64_t>::StorageType::const_iterator iter;

    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 38);
    iter = fp1->getNonzeroElements().find(374073638);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 6);
    iter = fp1->getNonzeroElements().find(464351883);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 2);
    iter = fp1->getNonzeroElements().find(1949583554);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 6);
    iter = fp1->getNonzeroElements().find(4105342207);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 1);
    iter = fp1->getNonzeroElements().find(794080973);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 1);
    iter = fp1->getNonzeroElements().find(3826517238);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 2);

    delete m1;
    delete fp1;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        RDKitFP::getRDKitFPGenerator<std::uint64_t>();
    ROMol *m1 = SmilesToMol("Cl");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint64_t> *fp1;
    SparseIntVect<std::uint64_t>::StorageType::const_iterator iter;

    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 0);

    delete m1;
    delete fp1;
    delete generator;
  }
  // additional outputs are not fully functional yet
  /*{
    FingerprintGenerator<std::uint64_t> *generator =
        RDKitFP::getRDKitFPGenerator<std::uint64_t>();
    ROMol *m1 = SmilesToMol("CCCO");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint64_t> *fp1;
    SparseIntVect<std::uint64_t>::StorageType::const_iterator iter;
    std::map<std::uint64_t, std::vector<std::vector<int>>> bitInfo;
    std::map<std::uint64_t, std::vector<std::vector<int>>>::const_iterator
        iter2;

    fp1 = getUnfoldedRDKFingerprintMol(*m1, 1, 7, true, true, true, nullptr,
                                       nullptr, nullptr, &bitInfo);
    TEST_ASSERT(fp1);

    TEST_ASSERT(fp1->getNonzeroElements().size() == 5);
    iter = fp1->getNonzeroElements().find(1524090560);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 1);
    iter = fp1->getNonzeroElements().find(1940446997);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 1);
    iter = fp1->getNonzeroElements().find(3977409745);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 1);
    iter = fp1->getNonzeroElements().find(4274652475);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 1);
    iter = fp1->getNonzeroElements().find(4275705116);
    TEST_ASSERT(iter != fp1->getNonzeroElements().end() && iter->second == 2);

    iter2 = bitInfo.find(4275705116);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0][0] == 0);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second.size() == 2 &&
                iter2->second[1][0] == 1);
    iter2 = bitInfo.find(4274652475);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0][0] == 2);
    iter2 = bitInfo.find(3977409745);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0].size() == 3 &&
                iter2->second[0][0] == 0);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0].size() == 3 &&
                iter2->second[0][1] == 1);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0].size() == 3 &&
                iter2->second[0][2] == 2);
    iter2 = bitInfo.find(1524090560);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0].size() == 2 &&
                iter2->second[0][0] == 1);
    TEST_ASSERT(iter2 != bitInfo.end() && iter2->second[0].size() == 2 &&
                iter2->second[0][1] == 2);

    delete m1;
  }*/
}

void testTopologicalTorsionFPOld() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test topological torsion fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::int64_t> *fpSigned;
  SparseIntVect<std::uint64_t> *fp, *fpOld;

  FingerprintGenerator<std::uint64_t> *fpGenerator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();

  for (const auto &sm : smis) {
    mol = SmilesToMol(sm);
    fp = fpGenerator->getSparseCountFingerprint(*mol);
    fpSigned = AtomPairs::getTopologicalTorsionFingerprint(*mol);

    fpOld = new SparseIntVect<std::uint64_t>(fp->getLength());
    std::map<std::int64_t, int> nz = fpSigned->getNonzeroElements();
    for (auto &it : nz) {
      fpOld->setVal(static_cast<std::uint64_t>(it.first), it.second);
    }

    TEST_ASSERT(DiceSimilarity(*fp, *fpOld) == 1.0);
    TEST_ASSERT(*fp == *fpOld);

    delete mol;
    delete fp;
    delete fpOld;
    delete fpSigned;
  }

  delete fpGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testTorsions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Topological Torsions." << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint64_t> *fp;
  FingerprintGenerator<std::uint64_t> *generator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();

  mol = SmilesToMol("CCCC");
  fp = generator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 1);
  TEST_ASSERT(fp->getNonzeroElements().size() == 1);

  delete mol;
  delete fp;
  mol = SmilesToMol("CCCCO.Cl");
  fp = generator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 2);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  delete fp;
  delete generator;
  generator = TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
      false, 3);
  fp = generator->getSparseCountFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 3);

  delete mol;
  delete fp;
  delete generator;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

// previous version test name: testHashedTorsions
void testFoldedTorsions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Folded torsions." << std::endl;

  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getCountFingerprint(*mol);
    SparseIntVect<std::uint32_t> *fp2;
    fp2 = generator->getCountFingerprint(*mol);
    TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
    TEST_ASSERT(*fp1 == *fp2);

    delete mol;
    delete fp2;
    mol = SmilesToMol("c1ccccn1");
    fp2 = generator->getCountFingerprint(*mol);
    RANGE_CHECK(0.0, DiceSimilarity(*fp1, *fp2), 1.0);

    delete mol;
    delete fp1;
    delete fp2;
    delete generator;
  }

  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(false,
                                                                          6);
    ROMol *mol;
    mol = SmilesToMol("c1ccccc1");
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getCountFingerprint(*mol);
    SparseIntVect<std::uint32_t> *fp2;
    fp2 = generator->getCountFingerprint(*mol);
    TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
    TEST_ASSERT(*fp1 == *fp2);

    delete mol;
    delete fp2;
    mol = SmilesToMol("c1ccccn1");
    fp2 = generator->getCountFingerprint(*mol);
    RANGE_CHECK(0.0, DiceSimilarity(*fp1, *fp2), 1.0);

    delete mol;
    delete fp1;
    delete fp2;
    delete generator;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testBulkTorsions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Bulk Topological Torsions." << std::endl;

  FingerprintGenerator<std::uint64_t> *generator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
  std::string fName = getenv("RDBASE");
  fName += "/Projects/DbCLI/testData/pubchem.200.sdf";
  SDMolSupplier suppl(fName);
  while (!suppl.atEnd()) {
    ROMol *mol = suppl.next();
    SparseIntVect<std::uint64_t> *fp;
    fp = generator->getSparseCountFingerprint(*mol);
    TEST_ASSERT(fp->getTotalVal() > 1);
    delete mol;
    delete fp;
  }
  delete generator;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRootedTorsions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Rooted Topological Torsions." << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint64_t> *fp1, *fp2;
  std::vector<std::uint32_t> roots;
  FingerprintGenerator<std::uint64_t> *generator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();

  mol = SmilesToMol("OCCCC");
  roots.push_back(0);

  fp1 = generator->getSparseCountFingerprint(*mol);
  SparseIntVect<std::uint64_t>::StorageType nz1 = fp1->getNonzeroElements();
  TEST_ASSERT(nz1.size() > 0);

  fp2 = generator->getSparseCountFingerprint(*mol, &roots);
  SparseIntVect<std::uint64_t>::StorageType nz2 = fp2->getNonzeroElements();
  TEST_ASSERT(nz2.size() > 0);
  TEST_ASSERT(nz2.size() < nz1.size());

  for (SparseIntVect<std::uint64_t>::StorageType::const_iterator bIt =
           nz2.begin();
       bIt != nz2.end(); ++bIt) {
    TEST_ASSERT(bIt->second <= fp2->getVal(bIt->first));
  }

  delete mol;
  delete fp1;
  delete fp2;
  delete generator;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testIgnoreTorsions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test ignoring atoms in Topological Torsions."
                        << std::endl;

  FingerprintGenerator<std::uint64_t> *generator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
  {
    ROMol *mol;
    SparseIntVect<std::uint64_t> *fp1, *fp2;
    std::vector<std::uint32_t> roots;

    mol = SmilesToMol("OCCCC");
    roots.push_back(0);

    fp1 = generator->getSparseCountFingerprint(*mol);
    SparseIntVect<std::uint64_t>::StorageType nz1 = fp1->getNonzeroElements();
    TEST_ASSERT(nz1.size() == 2);

    fp2 = generator->getSparseCountFingerprint(*mol, &roots);
    SparseIntVect<std::uint64_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() == 1);

    for (SparseIntVect<std::uint64_t>::StorageType::const_iterator bIt =
             nz2.begin();
         bIt != nz2.end(); ++bIt) {
      TEST_ASSERT(bIt->second <= fp2->getVal(bIt->first));
    }

    delete mol;
    delete fp1;
    delete fp2;
  }
  {
    ROMol *mol;
    SparseIntVect<std::uint64_t> *fp2;
    std::vector<std::uint32_t> roots;

    mol = SmilesToMol("OCCCC");
    roots.push_back(1);

    fp2 = generator->getSparseCountFingerprint(*mol, nullptr, &roots);
    SparseIntVect<std::uint64_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() == 0);

    delete mol;
    delete fp2;
  }

  {
    ROMol *mol;
    SparseIntVect<std::uint64_t> *fp2;
    std::vector<std::uint32_t> roots;

    mol = SmilesToMol("OCCCC");
    roots.push_back(0);

    fp2 = generator->getSparseCountFingerprint(*mol, &roots, &roots);
    SparseIntVect<std::uint64_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() == 0);

    delete mol;
    delete fp2;
  }
  delete generator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testPairsAndTorsionsOptions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing atom pair and torsions options" << std::endl;
  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>();
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    SparseIntVect<std::uint32_t> *fp1 =
        generator->getSparseCountFingerprint(*m1);
    SparseIntVect<std::uint32_t> *fp2 =
        generator->getSparseCountFingerprint(*m2);
    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 = generator->getSparseCountFingerprint(*m1, nullptr, nullptr, -1,
                                               nullptr, &invars);
    fp2 = generator->getSparseCountFingerprint(*m2, nullptr, nullptr, -1,
                                               nullptr, &invars);

    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>(1, 5, false, true,
                                                      nullptr, true, 1024);
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    SparseIntVect<std::uint32_t> *fp1, *fp2;
    fp1 = generator->getCountFingerprint(*m1);
    fp2 = generator->getCountFingerprint(*m2);
    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 = generator->getCountFingerprint(*m1, nullptr, nullptr, -1, nullptr,
                                         &invars);
    fp2 = generator->getCountFingerprint(*m2, nullptr, nullptr, -1, nullptr,
                                         &invars);

    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>(1, 5, false, true,
                                                      nullptr, true, 1024);
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    ExplicitBitVect *fp1, *fp2;
    fp1 = generator->getFingerprint(*m1);
    fp2 = generator->getFingerprint(*m2);
    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 =
        generator->getFingerprint(*m1, nullptr, nullptr, -1, nullptr, &invars);
    fp2 =
        generator->getFingerprint(*m2, nullptr, nullptr, -1, nullptr, &invars);

    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    SparseIntVect<std::uint64_t> *fp1, *fp2;
    fp1 = generator->getSparseCountFingerprint(*m1);
    fp2 = generator->getSparseCountFingerprint(*m2);
    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;
    UINT_VECT invars(6, 1);
    fp1 = generator->getSparseCountFingerprint(*m1, nullptr, nullptr, -1,
                                               nullptr, &invars);
    fp2 = generator->getSparseCountFingerprint(*m2, nullptr, nullptr, -1,
                                               nullptr, &invars);

    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }

  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    SparseIntVect<std::uint32_t> *fp1, *fp2;
    fp1 = generator->getCountFingerprint(*m1);
    fp2 = generator->getCountFingerprint(*m2);
    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 = generator->getCountFingerprint(*m1, nullptr, nullptr, -1, nullptr,
                                         &invars);
    fp2 = generator->getCountFingerprint(*m2, nullptr, nullptr, -1, nullptr,
                                         &invars);

    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    ExplicitBitVect *fp1, *fp2;
    fp1 = generator->getFingerprint(*m1);
    fp2 = generator->getFingerprint(*m2);
    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 =
        generator->getFingerprint(*m1, nullptr, nullptr, -1, nullptr, &invars);
    fp2 =
        generator->getFingerprint(*m2, nullptr, nullptr, -1, nullptr, &invars);

    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }

  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
            false, 4, nullptr, true, 1024, {1, 2, 4, 8});
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    SparseIntVect<std::uint32_t> *fp1 = generator->getCountFingerprint(*m1);
    SparseIntVect<std::uint32_t> *fp2 = generator->getCountFingerprint(*m2);

    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 = generator->getCountFingerprint(*m1, nullptr, nullptr, -1, nullptr,
                                         &invars);
    fp2 = generator->getCountFingerprint(*m2, nullptr, nullptr, -1, nullptr,
                                         &invars);
    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
            false, 4, nullptr, true, 1024, {1, 2, 4, 8});
    std::string smi = "C1=CC=CC=C1";
    RWMol *m1 = SmilesToMol(smi);
    TEST_ASSERT(m1);
    smi = "C1=CC=CC=N1";
    RWMol *m2 = SmilesToMol(smi);
    TEST_ASSERT(m2);
    ExplicitBitVect *fp1 = generator->getFingerprint(*m1);
    ExplicitBitVect *fp2 = generator->getFingerprint(*m2);

    TEST_ASSERT(*fp1 != *fp2);
    delete fp1;
    delete fp2;

    UINT_VECT invars(6, 1);
    fp1 =
        generator->getFingerprint(*m1, nullptr, nullptr, -1, nullptr, &invars);
    fp2 =
        generator->getFingerprint(*m2, nullptr, nullptr, -1, nullptr, &invars);
    TEST_ASSERT(*fp1 == *fp2);
    delete m1;
    delete m2;
    delete fp1;
    delete fp2;
    delete generator;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testChiralTorsions() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test Topological Torsions including info about chirality."
      << std::endl;
  // Folded size 4096 is causing a collision for this test hence the unusual
  // 4000 size
  FingerprintGenerator<std::uint64_t> *generator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
          false, 4, nullptr, true, 4096, {1, 2, 4, 8});
  FingerprintGenerator<std::uint64_t> *generatorChirality =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
          true, 4, nullptr, true, 4096, {1, 2, 4, 8});

  ROMol *m1, *m2, *m3;

  m1 = SmilesToMol("CC[CH](F)Cl");
  TEST_ASSERT(m1);
  m2 = SmilesToMol("CC[C@H](F)Cl");
  TEST_ASSERT(m1);
  m3 = SmilesToMol("CC[C@@H](F)Cl");
  TEST_ASSERT(m1);

  {
    SparseIntVect<std::uint64_t> *fp1, *fp2, *fp3;
    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 2);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 2);
    fp2 = generator->getSparseCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 2);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 2);
    fp3 = generator->getSparseCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 2);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 2);

    TEST_ASSERT((*fp1) == (*fp2));
    TEST_ASSERT((*fp1) == (*fp3));
    TEST_ASSERT((*fp2) == (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;

    fp1 = generatorChirality->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 2);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 2);
    fp2 = generatorChirality->getSparseCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 2);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 2);
    fp3 = generatorChirality->getSparseCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 2);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 2);

    TEST_ASSERT((*fp1) != (*fp2));
    TEST_ASSERT((*fp1) != (*fp3));
    TEST_ASSERT((*fp2) != (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;
  }

  {
    SparseIntVect<std::uint32_t> *fp1, *fp2, *fp3;
    fp1 = generator->getCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 2);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 2);
    fp2 = generator->getCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 2);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 2);
    fp3 = generator->getCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 2);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 2);

    TEST_ASSERT((*fp1) == (*fp2));
    TEST_ASSERT((*fp1) == (*fp3));
    TEST_ASSERT((*fp2) == (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;

    fp1 = generatorChirality->getCountFingerprint(*m1);
    TEST_ASSERT(fp1->getTotalVal() == 2);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 2);
    fp2 = generatorChirality->getCountFingerprint(*m2);
    TEST_ASSERT(fp2->getTotalVal() == 2);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 2);
    fp3 = generatorChirality->getCountFingerprint(*m3);
    TEST_ASSERT(fp3->getTotalVal() == 2);
    TEST_ASSERT(fp3->getNonzeroElements().size() == 2);

    TEST_ASSERT((*fp1) != (*fp2));
    TEST_ASSERT((*fp1) != (*fp3));
    TEST_ASSERT((*fp2) != (*fp3));

    delete fp1;
    delete fp2;
    delete fp3;
  }

  delete m1;
  delete m2;
  delete m3;
  delete generator;
  delete generatorChirality;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue25() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test GitHub Issue 25: fingerprint backwards compatibility."
      << std::endl;

  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    ROMol *m1 = SmilesToMol("CCCCO");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint64_t> *fp1;
    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    TEST_ASSERT(fp1->getTotalVal() == 2);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 2);
    TEST_ASSERT((*fp1)[4437590048LL] == 1);
    TEST_ASSERT((*fp1)[12893306913LL] == 1);
    delete fp1;
    delete m1;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(
            false, 4, nullptr, true, 1000, {1, 2, 4, 8});
    ROMol *m1 = SmilesToMol("CCCCO");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    TEST_ASSERT(fp1->getTotalVal() == 2);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 2);
    TEST_ASSERT((*fp1)[24] == 1);
    TEST_ASSERT((*fp1)[288] == 1);
    delete fp1;
    delete m1;
    delete generator;
  }

  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>();
    ROMol *m1 = SmilesToMol("CCO");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    TEST_ASSERT(fp1->getTotalVal() == 3);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 3);
    TEST_ASSERT((*fp1)[558113] == 1);
    TEST_ASSERT((*fp1)[1590306] == 1);
    TEST_ASSERT((*fp1)[1590337] == 1);
    delete fp1;
    delete m1;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>();
    ROMol *m1 = SmilesToMol("CCO");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    TEST_ASSERT(fp1->getTotalVal() == 3);
    TEST_ASSERT(fp1->getNonzeroElements().size() == 3);
    TEST_ASSERT((*fp1)[1375] == 1);
    TEST_ASSERT((*fp1)[1423] == 1);
    TEST_ASSERT((*fp1)[1503] == 1);
    delete fp1;
    delete m1;
    delete generator;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue334() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test GitHub Issue 334: explicit Hs in SMILES "
                           "modifies atom pair (and topological torsion) FP."
                        << std::endl;

  {
    FingerprintGenerator<std::uint32_t> *generator =
        AtomPair::getAtomPairGenerator<std::uint32_t>();
    ROMol *m1 = SmilesToMol("N#C");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint32_t> *fp1;
    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    delete m1;

    m1 = SmilesToMol("N#[CH]");
    SparseIntVect<std::uint32_t> *fp2;
    fp2 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp2);
    delete m1;

    TEST_ASSERT(fp1->getTotalVal() == fp2->getTotalVal());
    TEST_ASSERT(fp1->getNonzeroElements().size() ==
                fp2->getNonzeroElements().size());
    TEST_ASSERT(*fp1 == *fp2);
    delete fp1;
    delete fp2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    ROMol *m1 = SmilesToMol("N#C");
    TEST_ASSERT(m1);
    SparseIntVect<std::uint64_t> *fp1;
    fp1 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp1);
    delete m1;

    m1 = SmilesToMol("N#[CH]");
    SparseIntVect<std::uint64_t> *fp2;
    fp2 = generator->getSparseCountFingerprint(*m1);
    TEST_ASSERT(fp2);
    delete m1;

    TEST_ASSERT(fp1->getTotalVal() == fp2->getTotalVal());
    TEST_ASSERT(fp1->getNonzeroElements().size() ==
                fp2->getNonzeroElements().size());
    TEST_ASSERT(*fp1 == *fp2);
    delete fp1;
    delete fp2;
    delete generator;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue811() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test GitHub Issue 811: rooted atom fingerprint "
                           "non identical for the same molecules #811"
                        << std::endl;
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(false,
                                                                          7);
    std::string dirName = getenv("RDBASE");
    dirName += "/Code/GraphMol/Fingerprints/testData/";

    ROMol *m1 = MolFileToMol(dirName + "github811a.mol");
    TEST_ASSERT(m1);
    ROMol *m2 = MolFileToMol(dirName + "github811b.mol");
    TEST_ASSERT(m2);

    SparseIntVect<std::uint64_t> *fp1, *fp2;

    std::vector<std::uint32_t> roots;
    roots.push_back(1);

    fp1 = generator->getSparseCountFingerprint(*m1, &roots);
    SparseIntVect<std::uint64_t>::StorageType nz1 = fp1->getNonzeroElements();
    TEST_ASSERT(nz1.size() == 4);
    fp2 = generator->getSparseCountFingerprint(*m1, &roots);
    SparseIntVect<std::uint64_t>::StorageType nz2 = fp2->getNonzeroElements();
    TEST_ASSERT(nz2.size() == 4);

    TEST_ASSERT(*fp1 == *fp2);

    delete fp1;
    delete fp2;
    delete m1;
    delete m2;
    delete generator;
  }
  {
    FingerprintGenerator<std::uint64_t> *generator =
        TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();
    ROMol *m1 = SmilesToMol("C1CC1");
    TEST_ASSERT(m1);

    SparseIntVect<std::uint64_t> *fp1;

    fp1 = generator->getSparseCountFingerprint(*m1);
    SparseIntVect<std::uint64_t>::StorageType nz1 = fp1->getNonzeroElements();
    TEST_ASSERT(nz1.size() == 1);
    delete fp1;
    delete m1;
    delete generator;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testBulkFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog)
      << "    Test convenience and bulk fp calculation functions" << std::endl;

  std::vector<std::pair<FingerprintGenerator<std::uint64_t> *, FPType>>
      testPairs;

  std::vector<const ROMol *> molVect;

  for (const auto &sm : smis) {
    molVect.push_back(SmilesToMol(sm));
  }

  testPairs.emplace_back(AtomPair::getAtomPairGenerator<std::uint64_t>(),
                         FPType::AtomPairFP);

  testPairs.emplace_back(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2),
      FPType::MorganFP);

  testPairs.emplace_back(RDKitFP::getRDKitFPGenerator<std::uint64_t>(),
                         FPType::RDKitFP);

  testPairs.emplace_back(
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>(),
      FPType::TopologicalTorsionFP);

  for (const auto &it : testPairs) {
    std::vector<SparseIntVect<std::uint64_t> *> *results =
        getSparseCountFPBulk(molVect, it.second);

    std::vector<SparseIntVect<std::uint64_t> *> compareRes;

    for (const auto &m : molVect) {
      compareRes.push_back(it.first->getSparseCountFingerprint(*m));
    }

    for (unsigned long i = 0; i < results->size(); ++i) {
      TEST_ASSERT(*((*results)[i]) == *compareRes[i]);

      delete (*results)[i];
      delete compareRes[i];
    }

    delete results;
  }

  for (auto &&m : molVect) {
    delete m;
  }
  for (auto &&t : testPairs) {
    delete t.first;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
#if 1
  testAtomPairFP();
  testAtomPairArgs();
  testAtomPairOld();
  // testAtomPairNonSparseBitvector();
  testAtomPairOutput();
  testFoldedAtomPairs();
  testRootedAtomPairs();
  testIgnoreAtomPairs();
  testChiralPairs();
  testMorganFP();
  testMorganFPOld();
  testMorganFPFromAtoms();
  testMorganFPBitVect();
  testMorganFPFeatureInvs();
  testMorganFPOptions();
  testGitHubIssue695();
  testGitHubIssue874();
  testInvariantGenerators();
  testCustomInvariants();
  testRDKitFP();
  testRDKFPUnfolded();
  testTopologicalTorsionFPOld();
  testTorsions();
  testFoldedTorsions();
  testBulkTorsions();
  testRootedTorsions();
  testIgnoreTorsions();
  testPairsAndTorsionsOptions();
  testChiralTorsions();
#endif
  testGitHubIssue25();
  testGitHubIssue334();
  testGitHubIssue811();
  testBulkFP();

  return 0;
}
