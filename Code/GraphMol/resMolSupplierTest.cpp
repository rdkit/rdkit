#include <RDGeneral/test.h>
#include <iostream>
#include <string>
#include <map>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Resonance.h>

using namespace RDKit;

void addFormalChargeIndices(const ROMol *mol,
                            std::map<unsigned int, int> &fcMap) {
  for (ROMol::ConstAtomIterator it = mol->beginAtoms(); it != mol->endAtoms();
       ++it) {
    int fc = (*it)->getFormalCharge();
    unsigned int idx = (*it)->getIdx();
    if (fc) {
      if (fcMap.find(idx) == fcMap.end()) {
        fcMap[idx] = fc;
      } else {
        fcMap[idx] += fc;
      }
    }
  }
}

int getTotalFormalCharge(const ROMol *mol) {
  int totalFormalCharge = 0;
  for (ROMol::ConstAtomIterator it = mol->beginAtoms(); it != mol->endAtoms();
       ++it) {
    totalFormalCharge += (*it)->getFormalCharge();
  }
  return totalFormalCharge;
}

void cmpFormalChargeBondOrder(const ROMol *mol1, const ROMol *mol2) {
  TEST_ASSERT(mol1->getNumAtoms() == mol2->getNumAtoms());
  TEST_ASSERT(mol1->getNumBonds() == mol2->getNumBonds());
  for (unsigned int i = 0; i < mol1->getNumAtoms(); ++i) {
    TEST_ASSERT(mol1->getAtomWithIdx(i)->getFormalCharge() ==
                mol2->getAtomWithIdx(i)->getFormalCharge());
  }
  for (unsigned int i = 0; i < mol1->getNumBonds(); ++i) {
    TEST_ASSERT(mol1->getBondWithIdx(i)->getBondType() ==
                mol2->getBondWithIdx(i)->getBondType());
  }
}

void testBaseFunctionality() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testBaseFunctionality" << std::endl;

  ResonanceMolSupplier *resMolSuppl;
  RWMol *mol;
  mol = SmilesToMol("CC");
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT((resMolSuppl->getAtomConjGrpIdx(0) == -1) &&
              (resMolSuppl->getAtomConjGrpIdx(1) == -1))
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->getNumConjGrps() == 0);
  TEST_ASSERT(resMolSuppl->length() == 1);
  TEST_ASSERT(resMolSuppl->getNumConjGrps() == 0);
  delete resMolSuppl;
  delete mol;

  mol = SmilesToMol("NC(=[NH2+])c1ccc(cc1)C(=O)[O-]");
  int totalFormalCharge = getTotalFormalCharge(mol);
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(!resMolSuppl->getIsEnumerated());
  TEST_ASSERT(resMolSuppl->length() == 4);
  TEST_ASSERT(resMolSuppl->getIsEnumerated());
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(!resMolSuppl->getIsEnumerated());
  resMolSuppl->enumerate();
  TEST_ASSERT(resMolSuppl->getIsEnumerated());
  auto *smol0 = (*resMolSuppl)[0];
  auto *smol1 = (*resMolSuppl)[1];
  TEST_ASSERT((smol0->getBondBetweenAtoms(0, 1)->getBondType() !=
               smol1->getBondBetweenAtoms(0, 1)->getBondType()) ||
              (smol0->getBondBetweenAtoms(9, 10)->getBondType() !=
               smol1->getBondBetweenAtoms(9, 10)->getBondType()));
  delete smol0;
  delete smol1;
  delete resMolSuppl;

  resMolSuppl =
      new ResonanceMolSupplier((ROMol &)*mol, ResonanceMolSupplier::KEKULE_ALL);
  TEST_ASSERT(resMolSuppl->length() == 8);
  std::set<Bond::BondType> bondTypeSet;
  smol0 = (*resMolSuppl)[0];
  smol1 = (*resMolSuppl)[1];
  // check that we actually have two alternate Kekule structures
  bondTypeSet.insert(smol0->getBondBetweenAtoms(3, 4)->getBondType());
  bondTypeSet.insert(smol1->getBondBetweenAtoms(3, 4)->getBondType());
  TEST_ASSERT(bondTypeSet.size() == 2);
  bondTypeSet.clear();
  delete smol0;
  delete smol1;
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                         ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
  TEST_ASSERT(resMolSuppl->length() == 32);
  for (unsigned int i = 0; i < resMolSuppl->length(); ++i) {
    ROMol *resMol = (*resMolSuppl)[i];
    TEST_ASSERT(getTotalFormalCharge(resMol) == totalFormalCharge);
    delete resMol;
  }
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    TEST_ASSERT(getTotalFormalCharge(resMol) == totalFormalCharge);
    delete resMol;
  }
  resMolSuppl->reset();
  smol0 = (*resMolSuppl)[0];
  smol1 = resMolSuppl->next();
  cmpFormalChargeBondOrder(smol0, smol1);
  delete smol0;
  delete smol1;

  resMolSuppl->moveTo(12);
  smol0 = (*resMolSuppl)[12];
  smol1 = resMolSuppl->next();
  cmpFormalChargeBondOrder(smol0, smol1);
  delete resMolSuppl;
  delete smol0;
  delete smol1;

  resMolSuppl =
      new ResonanceMolSupplier((ROMol &)*mol,
                               ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS |
                                   ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                                   ResonanceMolSupplier::UNCONSTRAINED_ANIONS,
                               10);
  TEST_ASSERT(resMolSuppl->length() == 10);
  delete resMolSuppl;

  resMolSuppl =
      new ResonanceMolSupplier((ROMol &)*mol,
                               ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS |
                                   ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                                   ResonanceMolSupplier::UNCONSTRAINED_ANIONS,
                               0);
  TEST_ASSERT(resMolSuppl->length() == 0);
  delete resMolSuppl;
  delete mol;

  mol = SmilesToMol("CC(C)C(C(=O)OC(C#N)c1cccc(Oc2ccccc2)c1)c3ccc(OC(F)F)cc3");
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
                                         ResonanceMolSupplier::KEKULE_ALL, 3);
  TEST_ASSERT(!resMolSuppl->getIsEnumerated());
  TEST_ASSERT(resMolSuppl->length() == 3);
  TEST_ASSERT(resMolSuppl->getIsEnumerated());
  delete resMolSuppl;
  delete mol;
}

void testBenzylCation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testBenzylCation" << std::endl;
  RWMol *mol = SmilesToMol("[CH2+]c1ccccc1");
  ResonanceMolSupplier *resMolSuppl;
  std::map<unsigned int, int> fcMap;
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 4);
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    addFormalChargeIndices(resMol, fcMap);
    delete resMol;
  }
  unsigned int indices[] = {0, 2, 4, 6};
  for (unsigned int &idx : indices) {
    TEST_ASSERT(fcMap.find(idx) != fcMap.end());
    TEST_ASSERT(fcMap[idx] == 1);
  }
  delete resMolSuppl;
  delete mol;
}

void testBenzylAnion() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testBenzylAnion" << std::endl;
  RWMol *mol = SmilesToMol("[CH2-]c1ccccc1");
  ResonanceMolSupplier *resMolSuppl;
  std::map<unsigned int, int> fcMap;
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 4);
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    addFormalChargeIndices(resMol, fcMap);
    delete resMol;
  }
  unsigned int indices[] = {0, 2, 4, 6};
  for (unsigned int &idx : indices) {
    TEST_ASSERT(fcMap.find(idx) != fcMap.end());
    TEST_ASSERT(fcMap[idx] == -1);
  }
  delete resMolSuppl;
  delete mol;
}

void testButadiene() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testButadiene" << std::endl;
  RWMol *mol = SmilesToMol("C=CC=C");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 1);
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS |
                         ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
  TEST_ASSERT(resMolSuppl->length() == 3);
  std::map<unsigned int, int> fcMap;
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    addFormalChargeIndices(resMol, fcMap);
    delete resMol;
  }
  unsigned int indices[] = {0, 3};
  for (unsigned int &idx : indices) {
    TEST_ASSERT(fcMap.find(idx) != fcMap.end());
    TEST_ASSERT(fcMap[idx] == 0);
  }
  delete resMolSuppl;
  delete mol;
}

void testZwitterion() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testZwitterion" << std::endl;
  RWMol *mol = SmilesToMol("NC=CC=CC=O");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 1);
  ROMol *mol0 = (*resMolSuppl)[0];
  TEST_ASSERT(mol0->getAtomWithIdx(0)->getFormalCharge() == 0);
  TEST_ASSERT(mol0->getAtomWithIdx(6)->getFormalCharge() == 0);
  delete mol0;
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 2);
  mol0 = (*resMolSuppl)[0];
  ROMol *mol1 = (*resMolSuppl)[1];
  TEST_ASSERT(mol1->getAtomWithIdx(0)->getFormalCharge() == 1);
  TEST_ASSERT(mol1->getAtomWithIdx(6)->getFormalCharge() == -1);
  delete mol0;
  delete mol1;
  delete resMolSuppl;
  delete mol;
}

void testChargeMigration() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testChargeMigration" << std::endl;
  RWMol *mol = SmilesToMol("C=CC=CC=C[CH2+]");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 4);
  std::map<unsigned int, int> fcMap;
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    addFormalChargeIndices(resMol, fcMap);
    delete resMol;
  }
  unsigned int indices[] = {0, 2, 4, 6};
  for (unsigned int &idx : indices) {
    TEST_ASSERT(fcMap.find(idx) != fcMap.end());
    TEST_ASSERT(fcMap[idx] == 1);
  }
  delete resMolSuppl;
  delete mol;
}

void testChargeSeparation1() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testChargeSeparation1" << std::endl;
  RWMol *mol = SmilesToMol("[NH2+]=C1C=CC(C=C1)=CN=CN");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 3);
  std::map<unsigned int, int> fcMap;
  for (unsigned int i = 0; i < 3; ++i) {
    ROMol *resMol = (*resMolSuppl)[i];
    addFormalChargeIndices(resMol, fcMap);
    delete resMol;
  }
  unsigned int indices[] = {0, 8, 10};
  for (unsigned int &idx : indices) {
    TEST_ASSERT(fcMap.find(idx) != fcMap.end());
    TEST_ASSERT(fcMap[idx] == 1);
  }
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 4);
  resMolSuppl->moveTo(3);
  ROMol *resMol = resMolSuppl->next();
  TEST_ASSERT(resMol->getAtomWithIdx(0)->getFormalCharge() == 1);
  TEST_ASSERT(resMol->getAtomWithIdx(10)->getFormalCharge() == 1);
  TEST_ASSERT(resMol->getAtomWithIdx(8)->getFormalCharge() == -1);
  delete resMol;
  delete resMolSuppl;
  delete mol;
}

void testChargeSeparation2() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testChargeSeparation2" << std::endl;
  RWMol *mol = SmilesToMol("NC(C(=O)[O-])=CC=CC=O");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 2);
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 4);
  // less charge separation in the first 2,
  // more charge separation in the last 6
  // the last 4 feature carbanions
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
  TEST_ASSERT(resMolSuppl->length() == 8);
  std::map<unsigned int, int> fcMap;
  for (unsigned int i = 0; i < 2; ++i) {
    ROMol *resMol = (*resMolSuppl)[i];
    addFormalChargeIndices(resMol, fcMap);
    TEST_ASSERT(fcMap.size() == 1);
    fcMap.clear();
    delete resMol;
  }
  for (unsigned int i = 2; i < 8; ++i) {
    ROMol *resMol = (*resMolSuppl)[i];
    addFormalChargeIndices(resMol, fcMap);
    TEST_ASSERT(fcMap.size() == 3);
    fcMap.clear();
    delete resMol;
  }
  for (unsigned int i = 0; i < 4; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 6) && ((*it)->getFormalCharge() == -1)) {
        haveCarbanion = true;
      }
    }
    TEST_ASSERT(!haveCarbanion);
    delete resMol;
  }
  for (unsigned int i = 4; i < 8; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 6) && ((*it)->getFormalCharge() == -1)) {
        haveCarbanion = true;
      }
    }
    TEST_ASSERT(haveCarbanion);
    delete resMol;
  }
  delete resMolSuppl;
  delete mol;
}

void testMultipleConjGroups1() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testMultipleConjGroups1" << std::endl;
  RWMol *mol = SmilesToMol("NC(C(=O)[O-])=CC=C(CCC(=O)[O-])C=O");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 8);
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
  TEST_ASSERT(resMolSuppl->length() == 16);
  for (unsigned int i = 0; i < 8; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 6) && ((*it)->getFormalCharge() == -1)) {
        haveCarbanion = true;
      }
    }
    TEST_ASSERT(!haveCarbanion);
    delete resMol;
  }
  for (unsigned int i = 8; i < 16; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 6) && ((*it)->getFormalCharge() == -1)) {
        haveCarbanion = true;
      }
    }
    TEST_ASSERT(haveCarbanion);
    delete resMol;
  }
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::UNCONSTRAINED_ANIONS, 8);
  TEST_ASSERT(resMolSuppl->length() == 8);
  for (unsigned int i = 0; i < 8; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 6) && ((*it)->getFormalCharge() == -1)) {
        haveCarbanion = true;
      }
    }
    TEST_ASSERT(!haveCarbanion);
    delete resMol;
  }
  delete resMolSuppl;
  delete mol;
}

void testMultipleConjGroups2() {
  // if breadth-first resonance structure enumeration works properly
  // all N(2-) should be at the end of the supplier
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testMultipleConjGroups2" << std::endl;
  RWMol *mol = SmilesToMol("[N-]=[N+]=NCN=[N+]=[N-]");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 9);
  bool haveFirstDoubleNeg = false;
  for (unsigned int i = 0; i < resMolSuppl->length(); ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveDoubleNeg = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 7) && ((*it)->getFormalCharge() == -2)) {
        haveDoubleNeg = true;
      }
    }
    if (!haveFirstDoubleNeg && haveDoubleNeg) {
      haveFirstDoubleNeg = true;
    }
    TEST_ASSERT(!haveFirstDoubleNeg || (haveFirstDoubleNeg && haveDoubleNeg));
    delete resMol;
  }
  delete resMolSuppl;
  delete mol;
}

void testDimethylMalonate() {
  // carbanion should be last
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testDimethylMalonate" << std::endl;
  RWMol *mol = SmilesToMol("COC(=O)[CH-]C(=O)OC");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 3);
  for (unsigned int i = 0; i < 3; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 6) && ((*it)->getFormalCharge() == -1)) {
        haveCarbanion = true;
      }
    }
    TEST_ASSERT(((i < 2) && !haveCarbanion) || ((i == 2) && haveCarbanion));
    delete resMol;
  }
  delete resMolSuppl;
  delete mol;
}

void testMethylAcetate() {
  // cation on oxygen should appear only when UNCONSTRAINED_CATIONS is set
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testMethylAcetate" << std::endl;
  RWMol *mol = SmilesToMol("CC(=O)OC");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 1);
  const ROMol *resMol = (*resMolSuppl)[0];
  bool haveCation = false;
  for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
       it != resMol->endAtoms(); ++it) {
    if (((*it)->getAtomicNum() == 8) && ((*it)->getFormalCharge() == 1)) {
      haveCation = true;
    }
  }
  delete resMol;
  TEST_ASSERT(!haveCation);
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier(
      (ROMol &)*mol, ResonanceMolSupplier::UNCONSTRAINED_CATIONS);
  TEST_ASSERT(resMolSuppl->length() == 2);
  for (unsigned int i = 0; i < 2; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCation = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
         it != resMol->endAtoms(); ++it) {
      if (((*it)->getAtomicNum() == 8) && ((*it)->getFormalCharge() == 1)) {
        haveCation = true;
      }
    }
    TEST_ASSERT(((i == 0) && !haveCation) || ((i == 1) && haveCation));
    delete resMol;
  }
  delete resMolSuppl;
  delete mol;
}

void testPyranium2_6dicarboxylate() {
  // when there is no alternative cation on elements right of N
  // should be accepted even though the total formal charge is not
  // positive
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testPyranium2_6dicarboxylate" << std::endl;
  RWMol *mol = SmilesToMol("[o+]1c(C(=O)[O-])cccc1C(=O)[O-]");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 4);
  for (unsigned int i = 0; i < 4; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    TEST_ASSERT(resMol->getAtomWithIdx(0)->getFormalCharge() == 1);
    delete resMol;
  }
  delete resMolSuppl;
  delete mol;
}

void testSubstructMatchAcetate() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testSubstructMatchAcetate" << std::endl;
  RWMol *mol = SmilesToMol("CC(=O)[O-]");
  RWMol *query = SmartsToMol("C(=O)[O-]");
  ResonanceMolSupplier *resMolSuppl;
  unsigned int n;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  std::vector<MatchVectType> matchVect;
  n = SubstructMatch(*mol, *query, matchVect, false, true, false, false);
  TEST_ASSERT(n == 1);
  matchVect.clear();
  n = SubstructMatch(*resMolSuppl, *query, matchVect, false, true, false,
                     false);
  TEST_ASSERT(n == 2);
  delete mol;
  delete query;
  delete resMolSuppl;
}

void testSubstructMatchDMAP() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testSubstructMatchDMAP" << std::endl;
  RWMol *mol = SmilesToMol("C(C)Nc1cc[nH+]cc1");
  RWMol *query = SmartsToMol("[#7+]");
  ResonanceMolSupplier *resMolSuppl;
  unsigned int n;
  MatchVectType p;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  std::vector<MatchVectType> matchVect;
  n = SubstructMatch(*mol, *query, matchVect, false, true, false, false);
  TEST_ASSERT((n == 1) && (matchVect.size() == 1));
  p = matchVect[0];
  TEST_ASSERT((p.size() == 1) && (p[0].second == 6));
  matchVect.clear();
  n = SubstructMatch(*resMolSuppl, *query, matchVect, false, true, false,
                     false);
  TEST_ASSERT((n == 2) && (matchVect.size() == 2));
  std::vector<unsigned int> v;
  p = matchVect[0];
  TEST_ASSERT(p.size() == 1);
  v.push_back(p[0].second);
  p = matchVect[1];
  TEST_ASSERT(p.size() == 1);
  v.push_back(p[0].second);
  std::sort(v.begin(), v.end());
  TEST_ASSERT(v[0] == 2);
  TEST_ASSERT(v[1] == 6);
  delete mol;
  delete query;
  delete resMolSuppl;
}

void setResidueFormalCharge(RWMol *mol, std::vector<RWMol *> &res, int fc) {
  for (std::vector<RWMol *>::const_iterator it = res.begin(); it != res.end();
       ++it) {
    std::vector<MatchVectType> matchVect;
    SubstructMatch(*mol, *(*it), matchVect);
    for (std::vector<MatchVectType>::const_iterator it = matchVect.begin();
         it != matchVect.end(); ++it) {
      mol->getAtomWithIdx((*it).back().second)->setFormalCharge(fc);
    }
  }
}

void getBtVectVect(ResonanceMolSupplier *resMolSuppl,
                   std::vector<std::vector<unsigned int>> &btVect2) {
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    std::vector<unsigned int> bt;
    bt.reserve(resMol->getNumBonds());
    for (ROMol::BondIterator bi = resMol->beginBonds();
         bi != resMol->endBonds(); ++bi) {
      bt.push_back(static_cast<unsigned int>((*bi)->getBondTypeAsDouble()));
    }
    btVect2.push_back(bt);
    delete resMol;
  }
  for (unsigned int i = 0; i < btVect2.size(); ++i) {
    bool same = true;
    for (unsigned int j = 0; same && (j < btVect2[i].size()); ++j) {
      if (!i) {
        continue;
      }
      if (same) {
        same = (btVect2[i][j] == btVect2[i - 1][j]);
      }
    }
    if (i) {
      TEST_ASSERT(!same);
    }
  }
}

void testCrambin() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testCrambin" << std::endl;
  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/FileParsers/test_data/1CRN.pdb";
  RWMol *crambin = PDBFileToMol(pathName);
  TEST_ASSERT(crambin);
  RWMol *query;
  unsigned int n;
  std::vector<RWMol *> res;
  // protonate NH2
  query = SmartsToMol("[Nh2][Ch;Ch2]");
  TEST_ASSERT(query);
  res.push_back(query);
  // protonate Arg
  query = SmartsToMol("[Nh][C]([Nh2])=[Nh]");
  TEST_ASSERT(query);
  res.push_back(query);
  setResidueFormalCharge(crambin, res, 1);
  for (std::vector<RWMol *>::const_iterator it = res.begin(); it != res.end();
       ++it) {
    delete *it;
  }
  res.clear();
  // deprotonate COOH
  query = SmartsToMol("C(=O)[Oh]");
  TEST_ASSERT(query);
  res.push_back(query);
  setResidueFormalCharge(crambin, res, -1);
  for (std::vector<RWMol *>::const_iterator it = res.begin(); it != res.end();
       ++it) {
    delete *it;
  }
  auto *resMolSupplST = new ResonanceMolSupplier((ROMol &)*crambin);
  TEST_ASSERT(resMolSupplST);
  // crambin has 2 Arg (3 resonance structures each); 1 Asp, 1 Glu
  // and 1 terminal COO- (2 resonance structures each)
  // so possible resonance structures are 3^2 * 2^3 = 72
  TEST_ASSERT(resMolSupplST->length() == 72);
  TEST_ASSERT(resMolSupplST->getNumConjGrps() == 56);
  RWMol *carboxylateQuery = SmartsToMol("C(=O)[O-]");
  RWMol *guanidiniumQuery = SmartsToMol("NC(=[NH2+])N");
  std::vector<MatchVectType> matchVect;
  n = SubstructMatch(*crambin, *carboxylateQuery, matchVect, false, true, false,
                     false, 1000);
  TEST_ASSERT(n == 3);
  n = SubstructMatch(*crambin, *carboxylateQuery, matchVect, true, true, false,
                     false, 1000);
  TEST_ASSERT(n == 3);
  n = SubstructMatch(*crambin, *guanidiniumQuery, matchVect, false, true, false,
                     false, 1000);
  TEST_ASSERT(n == 0);
  n = SubstructMatch(*crambin, *guanidiniumQuery, matchVect, true, true, false,
                     false, 1000);
  TEST_ASSERT(n == 0);
  n = SubstructMatch(*resMolSupplST, *carboxylateQuery, matchVect, false, true,
                     false, false, 1000);
  TEST_ASSERT(n == 6);
  n = SubstructMatch(*resMolSupplST, *carboxylateQuery, matchVect, true, true,
                     false, false, 1000);
  TEST_ASSERT(n == 3);
  n = SubstructMatch(*resMolSupplST, *guanidiniumQuery, matchVect, false, true,
                     false, false, 1000);
  TEST_ASSERT(n == 8);
  n = SubstructMatch(*resMolSupplST, *guanidiniumQuery, matchVect, true, true,
                     false, false, 1000);
  TEST_ASSERT(n == 2);
#ifdef RDK_TEST_MULTITHREADED
  std::vector<std::vector<unsigned int>> btVect2ST;
  getBtVectVect(resMolSupplST, btVect2ST);
  auto *resMolSupplMT = new ResonanceMolSupplier((ROMol &)*crambin, 0, 1000);
  TEST_ASSERT(resMolSupplMT);
  resMolSupplMT->setNumThreads(0);
  TEST_ASSERT(resMolSupplST->length() == resMolSupplMT->length());
  std::vector<std::vector<unsigned int>> btVect2MT;
  getBtVectVect(resMolSupplMT, btVect2MT);
  TEST_ASSERT(btVect2ST.size() == btVect2MT.size());
  for (unsigned int i = 0; i < btVect2ST.size(); ++i) {
    for (unsigned int j = 0; j < btVect2ST[i].size(); ++j) {
      TEST_ASSERT(btVect2ST[i][j] == btVect2MT[i][j]);
    }
  }
  ResonanceMolSupplier *ptr[2] = {resMolSupplST, resMolSupplMT};
  for (auto &i : ptr) {
    n = SubstructMatch(*i, *carboxylateQuery, matchVect, false, true, false,
                       false, 1000, 0);
    TEST_ASSERT(n == 6);
    n = SubstructMatch(*i, *carboxylateQuery, matchVect, true, true, false,
                       false, 1000, 0);
    TEST_ASSERT(n == 3);
    n = SubstructMatch(*i, *guanidiniumQuery, matchVect, false, true, false,
                       false, 1000, 0);
    TEST_ASSERT(n == 8);
    n = SubstructMatch(*i, *guanidiniumQuery, matchVect, true, true, false,
                       false, 1000, 0);
    TEST_ASSERT(n == 2);
  }
  delete resMolSupplMT;
#endif
  delete resMolSupplST;
  delete carboxylateQuery;
  delete guanidiniumQuery;
  delete crambin;
}

void testConjGrpPerception() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testConjGrpPerception" << std::endl;
  RWMol *mol1 = MolBlockToMol(R"SDF(
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
    3.7539   -1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4317   -0.5660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1571   -1.3568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1651   -0.6484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4397   -1.4393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3921   -2.9385    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7619   -0.7309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8095    0.7684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1316    1.4768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5349    1.5592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2127    0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0619    1.6417    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3841    0.9333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6587    1.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  1  0
  5  7  4  0
  7  8  4  0
  8  9  1  0
  8 10  4  0
 10 11  4  0
 11 12  4  0
 12 13  4  0
 13 14  1  0
 13  2  4  0
 11  4  4  0
M  END
)SDF");
  RWMol *mol2 = MolBlockToMol(R"SDF(
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
    1.0619   -1.6417    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2127   -0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5349   -1.5592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8095   -0.7684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7619    0.7309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4397    1.4393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1651    0.6484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1571    1.3568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4317    0.5660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7539    1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3841   -0.9333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6587   -1.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1316   -1.4768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3921    2.9385    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  2  3  4  0
  2  7  4  0
  7  8  4  0
  8  9  4  0
  9 10  1  0
  9 11  4  0
 11 12  1  0
 11  1  4  0
  6  7  4  0
  4 13  1  0
  6 14  1  0
M  END
)SDF");
  auto *resMolSuppl1 = new ResonanceMolSupplier(
      static_cast<ROMol &>(*mol1), ResonanceMolSupplier::KEKULE_ALL);
  TEST_ASSERT(resMolSuppl1->length() == 3);
  auto *resMolSuppl2 = new ResonanceMolSupplier(
      static_cast<ROMol &>(*mol2), ResonanceMolSupplier::KEKULE_ALL);
  TEST_ASSERT(resMolSuppl2->length() == 3);
  delete resMolSuppl1;
  delete resMolSuppl2;
  delete mol1;
  delete mol2;
}

void testGitHub1166() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testGitHub1166" << std::endl;
  RWMol *mol = SmilesToMol("NC(=[NH2+])c1ccc(cc1)C(=O)[O-]");
  auto *resMolSuppl = new ResonanceMolSupplier(
      *mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION |
                ResonanceMolSupplier::KEKULE_ALL);
  TEST_ASSERT(resMolSuppl->length() == 8);
  // check that formal charges on odd indices are in the same position
  // as on even indices
  for (unsigned int i = 0; i < resMolSuppl->length(); i += 2) {
    auto *smol0 = (*resMolSuppl)[i];
    auto *smol1 = (*resMolSuppl)[i + 1];
    TEST_ASSERT(smol0->getNumAtoms() == smol1->getNumAtoms());
    for (unsigned int atomIdx = 0; atomIdx < smol0->getNumAtoms(); ++atomIdx) {
      TEST_ASSERT(smol0->getAtomWithIdx(atomIdx)->getFormalCharge() ==
                  smol1->getAtomWithIdx(atomIdx)->getFormalCharge());
    }
    // check that bond orders are alternate on aromatic bonds between
    // structures on odd indices and structures on even indices
    TEST_ASSERT(smol0->getNumBonds() == smol1->getNumBonds());
    for (unsigned int bondIdx = 0; bondIdx < smol0->getNumBonds(); ++bondIdx) {
      TEST_ASSERT(
          (!smol0->getBondWithIdx(bondIdx)->getIsAromatic() &&
           !smol1->getBondWithIdx(bondIdx)->getIsAromatic() &&
           (smol0->getBondWithIdx(bondIdx)->getBondType() ==
            smol1->getBondWithIdx(bondIdx)->getBondType())) ||
          (smol0->getBondWithIdx(bondIdx)->getIsAromatic() &&
           smol1->getBondWithIdx(bondIdx)->getIsAromatic() &&
           (static_cast<int>(
                smol0->getBondWithIdx(bondIdx)->getBondTypeAsDouble() +
                smol1->getBondWithIdx(bondIdx)->getBondTypeAsDouble()) == 3)));
    }
    delete smol0;
    delete smol1;
  }
  delete resMolSuppl;
  delete mol;
}

void testGitHub3048() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testGitHub3048" << std::endl;
  RWMol *mol = SmilesToMol("C1CN3N(C1)c2ccccc2N=C3N");
  auto *resMolSuppl =
      new ResonanceMolSupplier(*mol, ResonanceMolSupplier::KEKULE_ALL);
  // This caused a segfault due to a null ptr being accessed (#3048)
  TEST_ASSERT(resMolSuppl->length() == 2);
  delete resMolSuppl;
  delete mol;
}

void testGitHub2597() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testGitHub2597" << std::endl;
  {
    class MyCallBack : public ResonanceMolSupplierCallback {
      bool operator()() override {
        TEST_ASSERT(getNumConjGrps() == 1);
        return (getNumStructures(0) < 12);
      }
    };
    class MyCallBack2 : public ResonanceMolSupplierCallback {
      bool operator()() override {
        TEST_ASSERT(getNumConjGrps() == 1);
        return (getNumDiverseStructures(0) < 8);
      }
    };
    RWMol *mol = SmilesToMol(
        "ClC1=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=NC(NC2=CC=CC3=C2C(=O)C2=CC="
        "CC=C2C3=O)=N1");
    auto *resMolSuppl = new ResonanceMolSupplier(*mol);
    TEST_ASSERT(resMolSuppl->length() == 1);
    delete resMolSuppl;
    resMolSuppl =
        new ResonanceMolSupplier(*mol, ResonanceMolSupplier::KEKULE_ALL);
    TEST_ASSERT(resMolSuppl->length() == 32);
    delete resMolSuppl;
    resMolSuppl = new ResonanceMolSupplier(
        *mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION, 10);
    TEST_ASSERT(!resMolSuppl->getProgressCallback());
    TEST_ASSERT(resMolSuppl->length() == 10);
    TEST_ASSERT(!resMolSuppl->wasCanceled());
    delete resMolSuppl;
    resMolSuppl = new ResonanceMolSupplier(
        *mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
    MyCallBack *callback = new MyCallBack();
    resMolSuppl->setProgressCallback(callback);
    TEST_ASSERT(resMolSuppl->getProgressCallback() == callback);
    TEST_ASSERT(resMolSuppl->length() == 12);
    TEST_ASSERT(resMolSuppl->wasCanceled());
    delete resMolSuppl;
    resMolSuppl = new ResonanceMolSupplier(
        *mol, ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
    resMolSuppl->setProgressCallback(new MyCallBack2());
    TEST_ASSERT(resMolSuppl->getProgressCallback());
    resMolSuppl->setProgressCallback(nullptr);
    TEST_ASSERT(!resMolSuppl->getProgressCallback());
    resMolSuppl->setProgressCallback(new MyCallBack2());

    TEST_ASSERT(resMolSuppl->length() == 9);
    TEST_ASSERT(resMolSuppl->wasCanceled());
    delete resMolSuppl;
    delete mol;
  }
}

void testGitHub3349() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testGitHub3349" << std::endl;
  RWMol *mol = SmilesToMol("CC(=O)[O-]->[*]");
  auto *resMolSuppl =
      new ResonanceMolSupplier(*mol, ResonanceMolSupplier::KEKULE_ALL);
  // This erroneously returned a single resonance structure
  // as dative and zero-order bonds were not accounted for (#3349)
  TEST_ASSERT(resMolSuppl->length() == 2);
  delete resMolSuppl;
  delete mol;
}

void testGitHub5406() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
                       << "testGitHub5406 and 4884" << std::endl;
  for (auto smiles_string : {"CC=[N+]=[N-]", "O=[N+][O-]"}) {
    std::unique_ptr<RWMol> mol(SmilesToMol(smiles_string));
    MolOps::addHs(*mol);
    std::unique_ptr<ResonanceMolSupplier> suppl =
        std::make_unique<ResonanceMolSupplier>((ROMol &)*mol);
    TEST_ASSERT(suppl->length() == 0);
  }
}

int main() {
  RDLog::InitLogs();
#if 1
  testBaseFunctionality();
  testBenzylCation();
  testBenzylAnion();
  testButadiene();
  testZwitterion();
  testChargeMigration();
  testChargeSeparation1();
  testChargeSeparation2();
  testMultipleConjGroups1();
  testMultipleConjGroups2();
  testDimethylMalonate();
  testMethylAcetate();
  testPyranium2_6dicarboxylate();
  testSubstructMatchAcetate();
  testSubstructMatchDMAP();
  testCrambin();
  testGitHub1166();
  testConjGrpPerception();
  testGitHub3048();
  testGitHub2597();
  testGitHub3349();
  testGitHub5406();
#endif
  return 0;
}
