

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDBoost/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

void testAtomPairFPGenerator() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "Test compatibility between atom pair "
                           "implementation for FingerprintGenerator"
                           " and old atom pairs implementation"
                        << std::endl;
  {
    ROMol *mol;
    SparseIntVect<boost::int32_t> *fp1, *fp2;
    SparseIntVect<boost::uint32_t> *fpu;

    FingerprintGenerator atomPairGenerator = AtomPair::getAtomPairGenerator();

    mol = SmilesToMol("CCC");
    fp1 = AtomPairs::getAtomPairFingerprint(*mol);
    fpu = atomPairGenerator.getFingerprint(*mol);
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
    fpu = atomPairGenerator.getFingerprint(*mol);
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

    atomPairGenerator.cleanUpResources();
  }
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  testAtomPairFPGenerator();

  return 0;
}