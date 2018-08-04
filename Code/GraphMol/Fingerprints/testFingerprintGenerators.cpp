

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>

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
      AtomPair::getAtomPairGenerator<std::uint32_t>(2);

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

  atomPairGenerator = AtomPair::getAtomPairGenerator<std::uint32_t>(1, 1);
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

  atomPairGenerator =
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, 30, true);
  fp = atomPairGenerator->getFingerprint(*mol);
  TEST_ASSERT(fp->getTotalVal() == 3);
  TEST_ASSERT(fp->getNonzeroElements().size() == 2);

  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, true);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, true);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, true);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c2, 1, true)) == 2);
  TEST_ASSERT(fp->getVal(AtomPair::getAtomPairCode(c1, c3, 2, true)) == 1);

  atomPairGenerator = AtomPair::getAtomPairGenerator<std::uint32_t>(
      1, 30, false, true, nullptr, false);
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
        AtomPair::getAtomPairGenerator<std::uint32_t>();

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
        AtomPair::getAtomPairGenerator<std::uint32_t>();
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
        AtomPair::getAtomPairGenerator<std::uint32_t>();
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
      AtomPair::getAtomPairGenerator<std::uint32_t>();
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

void testMorganFPOld() {
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

void testMorganFP() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Morgan Fingerprints." << std::endl;

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("CCCCC");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    fp = morganGenerator->getFoldedFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = MorganFingerprints::getFingerprint(*mol, 1);
    TEST_ASSERT(fp->getNonzeroElements().size() == 5);
    delete fp;
    fp = morganGenerator->getFoldedFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 5);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 7);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 7);
    delete fp;
    delete morganGenerator;

    delete mol;
  }
  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("O=C(O)CC1CC1");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 6);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = MorganFingerprints::getFingerprint(*mol, 1);
    TEST_ASSERT(fp->getNonzeroElements().size() == 12);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 16);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 17);
    delete fp;
    delete morganGenerator;

    delete mol;
  }

  {
    // test that the results aren't order dependent, i.e. that we're
    // "canonicalizing" the fps correctly
    ROMol *mol, *mol2;
    SparseIntVect<boost::uint32_t> *fp, *fp2;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("O=C(O)CC1CC1");
    mol2 = SmilesToMol("OC(=O)CC1CC1");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getFingerprint(*mol);
    fp2 = morganGenerator->getFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 6);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 6);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getFingerprint(*mol);
    fp2 = morganGenerator->getFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 12);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 12);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    fp = morganGenerator->getFingerprint(*mol);
    fp2 = morganGenerator->getFingerprint(*mol2);
    TEST_ASSERT(fp->getNonzeroElements().size() == 16);
    TEST_ASSERT(fp2->getNonzeroElements().size() == 16);
    TEST_ASSERT(*fp == *fp2);
    delete fp;
    delete fp2;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(3);
    fp = morganGenerator->getFingerprint(*mol);
    fp2 = morganGenerator->getFingerprint(*mol2);
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
    SparseIntVect<boost::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    mol = SmilesToMol("OCCCCO");
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 7);
    SparseIntVect<boost::uint32_t>::StorageType::const_iterator iter;
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
    SparseIntVect<boost::uint32_t> *fp;
    FingerprintGenerator<std::uint32_t> *morganGenerator;

    mol = SmilesToMol("CC(F)(Cl)C(F)(Cl)C");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    delete fp;
    delete morganGenerator;

    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 8);
    delete mol;
    delete fp;
    delete morganGenerator;

    mol = SmilesToMol("CC(F)(Cl)[C@](F)(Cl)C");
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    delete fp;
    delete morganGenerator;
    morganGenerator = MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 8);
    delete fp;
    delete morganGenerator;

    morganGenerator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(0, true, true);
    fp = morganGenerator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 4);
    delete fp;
    delete morganGenerator;

    morganGenerator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(1, true, true);
    fp = morganGenerator->getFingerprint(*mol);
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
    SparseIntVect<boost::uint32_t> *fp;
    std::vector<boost::uint32_t> atoms;
    atoms.push_back(0);

    mol = SmilesToMol("CCCCC");
    fp = radius0Generator->getFingerprint(*mol);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    fp = radius0Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 1);
    delete fp;

    fp = radius1Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    // tests issue 3415636
    fp = radius2Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 3);
    delete fp;

    delete mol;
  }

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    std::vector<boost::uint32_t> atoms;

    mol = SmilesToMol("CCCCC");
    fp = radius0Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 0);
    delete fp;

    fp = radius1Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 0);
    delete fp;

    delete mol;
  }

  {
    ROMol *mol;
    SparseIntVect<boost::uint32_t> *fp;
    std::vector<boost::uint32_t> atoms;
    atoms.push_back(0);

    mol = SmilesToMol("C(CC)CO");

    fp = radius0Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 1);
    delete fp;

    fp = radius1Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 2);
    delete fp;

    fp = radius2Generator->getFingerprint(*mol, &atoms);
    TEST_ASSERT(fp->getNonzeroElements().size() == 3);
    delete fp;

    // tests issue 3415636
    fp = radius3Generator->getFingerprint(*mol, &atoms);
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
        MorganFingerprint::getMorganGenerator<std::uint32_t>(0);
    FingerprintGenerator<std::uint32_t> *radius1Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(1);
    FingerprintGenerator<std::uint32_t> *radius2Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2);
    FingerprintGenerator<std::uint32_t> *radius3Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(3);

    mol = SmilesToMol("CCCCC");
    fp = radius0Generator->getFoldedFingerprintAsBitVect(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 4);
    delete fp;
    fp = radius1Generator->getFoldedFingerprintAsBitVect(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 9);
    delete fp;
    fp = radius2Generator->getFoldedFingerprintAsBitVect(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 12);
    delete fp;
    fp = radius3Generator->getFoldedFingerprintAsBitVect(*mol);
    TEST_ASSERT(fp->getNumOnBits() == 12);
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
    MorganFingerprint::MorganFeatureAtomInvGenerator *invGen =
        new MorganFingerprint::MorganFeatureAtomInvGenerator();
    std::vector<boost::uint32_t> *invars = invGen->getAtomInvariants(*mol);
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
    MorganFingerprint::MorganFeatureAtomInvGenerator *invGen =
        new MorganFingerprint::MorganFeatureAtomInvGenerator();
    std::vector<boost::uint32_t> *invars = invGen->getAtomInvariants(*mol);
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
    MorganFingerprint::MorganFeatureAtomInvGenerator *invGen =
        new MorganFingerprint::MorganFeatureAtomInvGenerator(&patterns);

    std::vector<boost::uint32_t> *invars = invGen->getAtomInvariants(*mol);
    TEST_ASSERT((*invars)[0] != 0);
    TEST_ASSERT((*invars)[1] != 0);
    TEST_ASSERT((*invars)[0] != (*invars)[1]);
    TEST_ASSERT((*invars)[1] == (*invars)[2]);
    TEST_ASSERT((*invars)[0] == (*invars)[7]);
    delete mol;
    delete invGen;
    delete invars;
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
    std::vector<boost::uint32_t> invars(3);
    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;

    m1 = SmilesToMol("CCC");
    TEST_ASSERT(m1);
    m2 = SmilesToMol("CC=C");
    TEST_ASSERT(m2);

    fp1 = generator->getFoldedFingerprintAsBitVect(*m1, nullptr, nullptr, -1,
                                                   nullptr, &invars);
    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;
    fp2 = generator->getFoldedFingerprintAsBitVect(*m2, nullptr, nullptr, -1,
                                                   nullptr, &invars);
    TEST_ASSERT((*fp1) != (*fp2));
    delete fp1;
    delete fp2;

    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;
    fp1 = generatorNoBoundType->getFoldedFingerprintAsBitVect(
        *m1, nullptr, nullptr, -1, nullptr, &invars);
    invars[0] = 1;
    invars[1] = 1;
    invars[2] = 1;
    fp2 = generatorNoBoundType->getFoldedFingerprintAsBitVect(
        *m2, nullptr, nullptr, -1, nullptr, &invars);
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

    fp1 = generator->getFoldedFingerprintAsBitVect(*m1);
    fp2 = generator->getFoldedFingerprintAsBitVect(*m2);
    fp3 = generator->getFoldedFingerprintAsBitVect(*m3);
    TEST_ASSERT((*fp1) == (*fp2));
    TEST_ASSERT((*fp1) == (*fp3));
    TEST_ASSERT((*fp2) == (*fp3));
    delete fp1;
    delete fp2;
    delete fp3;

    fp1 = generatorChirality->getFoldedFingerprintAsBitVect(*m1);
    fp2 = generatorChirality->getFoldedFingerprintAsBitVect(*m2);
    fp3 = generatorChirality->getFoldedFingerprintAsBitVect(*m3);
    TEST_ASSERT((*fp1) != (*fp2));
    TEST_ASSERT((*fp1) != (*fp3));
    TEST_ASSERT((*fp2) != (*fp3));
    delete fp1;
    delete fp2;
    delete fp3;

    delete m1;
    delete m2;
    delete generator;
    delete generatorChirality;
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
      AtomPair::getAtomPairGenerator<std::uint32_t>(1, radius, false, true,
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

  atomPairGenerator = AtomPair::getAtomPairGenerator<std::uint32_t>(
      1, radius, false, true, atomInvariantsGenerator, bondInvariantsGenerator);

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

void testTopologicalTorsionFPOld() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test topological torsion fp generator" << std::endl;

  ROMol *mol;
  SparseIntVect<std::int64_t> *fpSigned;
  SparseIntVect<std::uint64_t> *fp, *fpOld;

  FingerprintGenerator<std::uint64_t> *fpGenerator =
      TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>();

  BOOST_FOREACH (std::string sm, smis) {
    mol = SmilesToMol(sm);
    fp = fpGenerator->getFingerprint(*mol);
    fpSigned = AtomPairs::getTopologicalTorsionFingerprint(*mol);

    fpOld = new SparseIntVect<boost::uint64_t>(fp->getLength());
    std::map<boost::int64_t, int> nz = fpSigned->getNonzeroElements();
    for (std::map<boost::int64_t, int>::iterator it = nz.begin();
         it != nz.end(); it++) {
      fpOld->setVal(static_cast<boost::uint64_t>(it->first), it->second);
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

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  testAtomPairFP();
  testAtomPairArgs();
  testAtomPairOld();
  testAtomPairBitvector();
  testAtomPairFoldedBitvector();
  testAtomPairOutput();
  testMorganFP();
  testMorganFPOld();
  testMorganFPFromAtoms();
  testMorganFPBitVect();
  testMorganFPFeatureInvs();
  testMorganFPOptions();
  testInvariantGenerators();
  testCustomInvariants();
  testRDKitFP();
  testTopologicalTorsionFPOld();

  return 0;
}
