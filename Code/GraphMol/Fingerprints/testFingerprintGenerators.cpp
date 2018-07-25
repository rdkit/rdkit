

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.cpp>
#include <GraphMol/Fingerprints/MorganGenerator.cpp>
#include <GraphMol/Fingerprints/RDKitFPGenerator.cpp>

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

void testAtomCodes() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test atom codes AtomPairGenerator" << std::endl;

  ROMol *mol;
  std::uint32_t tgt;
  std::uint32_t c1, c2, c3;

  mol = SmilesToMol("C=C");
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false) ==
              AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false));
  tgt = 1 | (1 | 1 << AtomPair::numPiBits) << AtomPair::numBranchBits;
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false) == tgt);
  tgt = 1 << AtomPair::numBranchBits |
        1 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(0), 1, false) == tgt);

  delete mol;
  mol = SmilesToMol("C#CO");
  tgt = 1 | 2 << AtomPair::numBranchBits |
        1 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false) == tgt);
  tgt = 2 | 2 << AtomPair::numBranchBits |
        1 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false) == tgt);
  tgt = 1 | 0 << AtomPair::numBranchBits |
        3 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false) == tgt);

  delete mol;
  mol = SmilesToMol("CC(O)C(O)(O)C");
  tgt = 1 | 0 << AtomPair::numBranchBits |
        1 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(1), 2, false) == tgt);
  tgt = 2 | 0 << AtomPair::numBranchBits |
        1 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(3), 2, false) == tgt);

  delete mol;
  mol = SmilesToMol("C=CC(=O)O");
  tgt = 0 | 0 << AtomPair::numBranchBits |
        3 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(4), 1, false) == tgt);
  tgt = 3 | 1 << AtomPair::numBranchBits |
        1 << (AtomPair::numBranchBits + AtomPair::numPiBits);
  TEST_ASSERT(AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false) == tgt);

  delete mol;

  mol = SmilesToMol("CCCCC");
  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  tgt = 1 | (std::min(c1, c2) | std::max(c1, c2) << AtomPair::codeSize)
                << AtomPair::numPathBits;
  TEST_ASSERT(AtomPair::getAtomPairCode(c1, c2, 1, false) == tgt);
  TEST_ASSERT(AtomPair::getAtomPairCode(c2, c1, 1, false) == tgt);
  tgt = 2 | (std::min(c1, c3) | std::max(c1, c3) << AtomPair::codeSize)
                << AtomPair::numPathBits;
  TEST_ASSERT(AtomPair::getAtomPairCode(c1, c3, 2, false) == tgt);
  TEST_ASSERT(AtomPair::getAtomPairCode(c3, c1, 2, false) == tgt);

  delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testAtomPairFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test atom-pair fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp;
  std::uint32_t c1, c2, c3;

  FingerprintGenerator<std::uint32_t> *atomPairGenerator =
      AtomPair::getAtomPairGenerator();

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getFingerprint(*mol);
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
  fp = atomPairGenerator->getFingerprint(*mol);
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
      AtomPair::getAtomPairGenerator(2);

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 1);
  TEST_ASSERT(fp->getNonzeroElements().size() == 1);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 0);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 1);

  delete fp;
  delete atomPairGenerator;

  atomPairGenerator = AtomPair::getAtomPairGenerator(1, 1);
  fp = atomPairGenerator->getFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 2);
  TEST_ASSERT(fp->getNonzeroElements().size() == 1);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, false)) == 0);

  delete fp;
  delete atomPairGenerator;

  atomPairGenerator = AtomPair::getAtomPairGenerator(1, 30, true);
  fp = atomPairGenerator->getFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, true);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, true);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, true);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, true)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, true)) == 1);

  atomPairGenerator = AtomPair::getAtomPairGenerator(1, 30, false, true, false);
  fp = atomPairGenerator->getFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 2);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, false)) == 1);
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
    SparseIntVect<boost::int32_t> *fp1, *fp2;
    SparseIntVect<boost::uint32_t> *fpu;

    FingerprintGenerator<std::uint32_t> *atomPairGenerator =
        AtomPair::getAtomPairGenerator();

    BOOST_FOREACH (std::string sm, smis) {
      mol = SmilesToMol(sm);
      fp1 = AtomPairs::getAtomPairFingerprint(*mol);
      fpu = atomPairGenerator->getFingerprint(*mol);
      fp2 = new SparseIntVect<boost::int32_t>(fpu->size());
      std::map<boost::uint32_t, int> nz = fpu->getNonzeroElements();
      for (std::map<boost::uint32_t, int>::iterator it = nz.begin();
           it != nz.end(); it++) {
        fp2->setVal(static_cast<boost::int32_t>(it->first), it->second);
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

void testAtomPairBitvector() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test consistency between different output types "
                           "for unfolded atom pair fingerprint"
                        << std::endl;
  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp1;
    SparseBitVect *fp2;

    FingerprintGenerator<std::uint32_t> *atomPairGenerator =
        AtomPair::getAtomPairGenerator();
    std::vector<std::uint32_t> defaultCountBounds = {1, 2, 4, 8};

    mol = SmilesToMol("CCC");
    fp1 = atomPairGenerator->getFingerprint(*mol);
    fp2 = atomPairGenerator->getFingerprintAsBitVect(*mol);

    std::map<boost::uint32_t, int> nz = fp1->getNonzeroElements();
    for (std::map<boost::uint32_t, int>::iterator it = nz.begin();
         it != nz.end(); it++) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); i++) {
        bool isSet = fp2->getBit(it->first * defaultCountBounds.size() + i);
        TEST_ASSERT(isSet == (it->second >= defaultCountBounds[i]));
      }
    }

    delete mol;
    delete fp1;
    delete fp2;

    mol = SmilesToMol("CC=O.Cl");
    fp1 = atomPairGenerator->getFingerprint(*mol);
    fp2 = atomPairGenerator->getFingerprintAsBitVect(*mol);

    nz = fp1->getNonzeroElements();
    for (std::map<boost::uint32_t, int>::iterator it = nz.begin();
         it != nz.end(); it++) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); i++) {
        bool isSet = fp2->getBit(it->first * defaultCountBounds.size() + i);
        TEST_ASSERT(isSet == (it->second >= defaultCountBounds[i]));
      }
    }

    delete mol;
    delete fp1;
    delete fp2;
    delete atomPairGenerator;

    BOOST_LOG(rdErrorLog) << "  done" << std::endl;
  }
}

void testAtomPairFoldedBitvector() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test consistency between different output types "
                           "for folded atom pair fingerprint"
                        << std::endl;
  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp1;
    ExplicitBitVect *fp2;

    FingerprintGenerator<std::uint32_t> *atomPairGenerator =
        AtomPair::getAtomPairGenerator();
    std::vector<std::uint32_t> defaultCountBounds = {1, 2, 4, 8};

    mol = SmilesToMol("CCC");
    fp1 = atomPairGenerator->getFoldedFingerprint(*mol);
    fp2 = atomPairGenerator->getFoldedFingerprintAsBitVect(*mol);

    std::map<boost::uint32_t, int> nz = fp1->getNonzeroElements();
    for (std::map<boost::uint32_t, int>::iterator it = nz.begin();
         it != nz.end(); it++) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); i++) {
        bool isSet = fp2->getBit(it->first * defaultCountBounds.size() + i);
        TEST_ASSERT(isSet == (it->second >= defaultCountBounds[i]));
      }
    }

    delete mol;
    delete fp1;
    delete fp2;

    mol = SmilesToMol("CC=O.Cl");
    fp1 = atomPairGenerator->getFoldedFingerprint(*mol);
    fp2 = atomPairGenerator->getFoldedFingerprintAsBitVect(*mol);

    nz = fp1->getNonzeroElements();
    for (std::map<boost::uint32_t, int>::iterator it = nz.begin();
         it != nz.end(); it++) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); i++) {
        bool isSet = fp2->getBit(it->first * defaultCountBounds.size() + i);
        TEST_ASSERT(isSet == (it->second >= defaultCountBounds[i]));
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

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp;
  std::uint32_t c1, c2, c3;
  AdditionalOutput additionalOutput = {nullptr, nullptr, nullptr, nullptr};
  std::vector<std::uint64_t> v;

  FingerprintGenerator<std::uint32_t> *atomPairGenerator =
      AtomPair::getAtomPairGenerator();
  mol = SmilesToMol("CCC");
  additionalOutput.atomToBits =
      new std::vector<std::vector<std::uint64_t>>(mol->getNumAtoms());
  fp = atomPairGenerator->getFingerprint(*mol, nullptr, nullptr, -1,
                                         &additionalOutput);
  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  v = additionalOutput.atomToBits->at(0);
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

  delete mol;
  delete fp;
  delete additionalOutput.atomToBits;
  delete atomPairGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMorganFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test morgan fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint32_t> *fp, *fpOld;
  int radius = 3;

  FingerprintGenerator<std::uint32_t> *morganGenerator =
      MorganFingerprint::getMorganGenerator<std::uint32_t>(radius);

  BOOST_FOREACH (std::string sm, smis) {
    mol = SmilesToMol(sm);
    fp = morganGenerator->getFingerprint(*mol);
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
  fp = morganGenerator->getFingerprint(*mol);
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
  fp = morganGenerator->getFingerprint(*mol);
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
  fp = morganGenerator->getFingerprint(*mol);
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
  fp = morganGenerator->getFingerprint(*mol);
  TEST_ASSERT(fp->getNonzeroElements().size() == 7);

  delete mol;
  delete fp;
  delete morganGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator = new MorganFingerprint::MorganAtomInvGenerator();
  bondInvariantsGenerator = nullptr;

  FingerprintGenerator<std::uint32_t> *atomPairGenerator =
      AtomPair::getAtomPairGenerator(1, radius, false, true, true,
                                     atomInvariantsGenerator,
                                     bondInvariantsGenerator);

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getFingerprint(*mol);

  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  delete mol;
  delete fp;
  delete atomPairGenerator;
  delete atomInvariantsGenerator;
  delete bondInvariantsGenerator;

  atomInvariantsGenerator =
      new MorganFingerprint::MorganFeatureAtomInvGenerator();
  bondInvariantsGenerator = nullptr;

  atomPairGenerator = AtomPair::getAtomPairGenerator(
      1, radius, false, true, true, atomInvariantsGenerator,
      bondInvariantsGenerator);

  mol = SmilesToMol("CCC");
  fp = atomPairGenerator->getFingerprint(*mol);

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
  fp1 = morganGenerator->getFingerprint(*mol);
  fp2 = defaultMorganGenerator->getFingerprint(*mol, nullptr, nullptr, -1,
                                               nullptr, customInvariants);
  TEST_ASSERT(DiceSimilarity(*fp1, *fp2) == 1.0);
  TEST_ASSERT(*fp1 == *fp2);

  delete mol;
  delete fp1, fp2;
  delete morganGenerator, defaultMorganGenerator;
  delete atomInvariantsGenerator;
  delete customInvariants;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testRDKitFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test RDKit fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::uint64_t> *fp, *fpTemp, *fpOld;

  FingerprintGenerator<std::uint64_t> *fpGenerator =
      RDKitFP::getRDKitFPGenerator<std::uint64_t>();

  BOOST_FOREACH (std::string sm, smis) {
    mol = SmilesToMol(sm);
    fp = fpGenerator->getFingerprint(*mol);
    fpTemp = getUnfoldedRDKFingerprintMol(*mol);

    // Old and new versions produce different length results, but data should be
    // the same
    std::map<std::uint64_t, int> nz = fpTemp->getNonzeroElements();
    fpOld = new SparseIntVect<std::uint64_t>(fp->getLength());
    for (auto it = nz.begin(); it != nz.end(); it++) {
      fpOld->setVal(it->first, it->second);
    }

    TEST_ASSERT(DiceSimilarity(*fp, *fpOld) == 1.0);
    TEST_ASSERT(*fp == *fpOld);

    delete mol;
    delete fp;
    delete fpOld;
    delete fpTemp;
  }

  delete fpGenerator;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  testAtomCodes();
  testAtomPairFP();
  testAtomPairArgs();
  testAtomPairOld();
  testAtomPairBitvector();
  testAtomPairFoldedBitvector();
  testAtomPairOutput();
  testMorganFP();
  testInvariantGenerators();
  testCustomInvariants();
  testRDKitFP();

  return 0;
}
