

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

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

  FingerprintGenerator *atomPairGenerator = AtomPair::getAtomPairGenerator();

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

    FingerprintGenerator *atomPairGenerator = AtomPair::getAtomPairGenerator();

    mol = SmilesToMol("CCC");
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

    mol = SmilesToMol("CC=O.Cl");
    fp1 = AtomPairs::getAtomPairFingerprint(*mol);
    fpu = atomPairGenerator->getFingerprint(*mol);
    fp2 = new SparseIntVect<boost::int32_t>(fpu->size());
    nz = fpu->getNonzeroElements();
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

    FingerprintGenerator *atomPairGenerator = AtomPair::getAtomPairGenerator();
    std::vector<std::uint32_t> defaultCountBounds = {1, 2, 4, 8};

    mol = SmilesToMol("CCC");
    fp1 = atomPairGenerator->getFingerprint(*mol);
    fp2 = atomPairGenerator->getFingerprintAsBitVect(*mol);

    std::map<boost::uint32_t, int> nz = fp1->getNonzeroElements();
    for (std::map<boost::uint32_t, int>::iterator it = nz.begin();
         it != nz.end(); it++) {
      for (unsigned int i = 0; i < defaultCountBounds.size(); i++) {
        if (it->second >= defaultCountBounds[i]) {
          bool isSet = fp2->getBit(it->first * defaultCountBounds.size() + i);
          TEST_ASSERT(isSet);
        }
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
        if (it->second >= defaultCountBounds[i]) {
          bool isSet = fp2->getBit(it->first * defaultCountBounds.size() + i);
          TEST_ASSERT(isSet);
        }
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

  FingerprintGenerator *atomPairGenerator = AtomPair::getAtomPairGenerator();
  mol = SmilesToMol("CCC");
  additionalOutput.atomToBits =
      new std::vector<std::vector<std::uint64_t>>(mol->getNumAtoms());
  fp = atomPairGenerator->getFingerprint(*mol, nullptr, nullptr, -1,
                                         &additionalOutput);
  c1 = AtomPair::getAtomCode(mol->getAtomWithIdx(0), 0, false);
  c2 = AtomPair::getAtomCode(mol->getAtomWithIdx(1), 0, false);
  c3 = AtomPair::getAtomCode(mol->getAtomWithIdx(2), 0, false);
  std::cout << additionalOutput.atomToBits->at(0).empty() << std::endl;
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
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  testAtomCodes();
  testAtomPairFP();
  testAtomPairOld();
  testAtomPairBitvector();
  testAtomPairOutput();

  return 0;
}