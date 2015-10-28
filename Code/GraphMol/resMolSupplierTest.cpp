#include <iostream>
#include <string>
#include <map>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Resonance.h>


using namespace RDKit;


void addFormalChargeIndices(const ROMol *mol,
  std::map<unsigned int, int> &fcMap) {
  for (ROMol::ConstAtomIterator it = mol->beginAtoms();
    it != mol->endAtoms(); ++it) {
    int fc = (*it)->getFormalCharge();
    unsigned int idx = (*it)->getIdx();
    if (fc) {
      if (fcMap.find(idx) == fcMap.end())
        fcMap[idx] = fc;
      else
        fcMap[idx] += fc;
    }
  }
}

int getTotalFormalCharge(const ROMol *mol) {
  int totalFormalCharge = 0;
  for (ROMol::ConstAtomIterator it = mol->beginAtoms();
    it != mol->endAtoms(); ++it)
    totalFormalCharge += (*it)->getFormalCharge();
  return totalFormalCharge;
}

void cmpFormalChargeBondOrder(const ROMol *mol1, const ROMol *mol2) {
  TEST_ASSERT(mol1->getNumAtoms() == mol2->getNumAtoms());
  TEST_ASSERT(mol1->getNumBonds() == mol2->getNumBonds());
  for (unsigned int i = 0; i < mol1->getNumAtoms(); ++i)
    TEST_ASSERT(mol1->getAtomWithIdx(i)->getFormalCharge()
      == mol2->getAtomWithIdx(i)->getFormalCharge());
  for (unsigned int i = 0; i < mol1->getNumBonds(); ++i)
    TEST_ASSERT(mol1->getBondWithIdx(i)->getBondType()
      == mol2->getBondWithIdx(i)->getBondType());
}

void testBaseFunctionality() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
    << "testBaseFunctionality" << std::endl;
  RWMol *mol = SmilesToMol("NC(=[NH2+])c1ccc(cc1)C(=O)[O-]");
  int totalFormalCharge = getTotalFormalCharge(mol);
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 4);
  TEST_ASSERT(((*resMolSuppl)[0]->getBondBetweenAtoms(0, 1)->getBondType()
    != (*resMolSuppl)[1]->getBondBetweenAtoms(0, 1)->getBondType())
   || ((*resMolSuppl)[0]->getBondBetweenAtoms(9, 10)->getBondType()
    != (*resMolSuppl)[1]->getBondBetweenAtoms(9, 10)->getBondType()));
  delete resMolSuppl;
  
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::KEKULE_ALL);
  TEST_ASSERT(resMolSuppl->length() == 8);
  std::map<Bond::BondType, bool> bondTypeMap;
  // check that we actually have two alternate Kekule structures
  bondTypeMap[(*resMolSuppl)[0]->getBondBetweenAtoms
    (3, 4)->getBondType()] = true;
  bondTypeMap[(*resMolSuppl)[1]->getBondBetweenAtoms
    (3, 4)->getBondType()] = true;
  TEST_ASSERT(bondTypeMap.size() == 2);
  bondTypeMap.clear();
  delete resMolSuppl;
  
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS
    | ResonanceMolSupplier::UNCONSTRAINED_CATIONS
    | ResonanceMolSupplier::UNCONSTRAINED_ANIONS);
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
  cmpFormalChargeBondOrder((*resMolSuppl)[0], resMolSuppl->next());
  resMolSuppl->moveTo(12);
  cmpFormalChargeBondOrder((*resMolSuppl)[12], resMolSuppl->next());
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS
    | ResonanceMolSupplier::UNCONSTRAINED_CATIONS
    | ResonanceMolSupplier::UNCONSTRAINED_ANIONS, 10);
  TEST_ASSERT(resMolSuppl->length() == 10);
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS
    | ResonanceMolSupplier::UNCONSTRAINED_CATIONS
    | ResonanceMolSupplier::UNCONSTRAINED_ANIONS, 0);
  TEST_ASSERT(resMolSuppl->length() == 0);
  delete resMolSuppl;
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
  unsigned int indices[] = { 0, 2, 4, 6 };
  for (unsigned int i = 0;
    i < sizeof(indices) / sizeof(unsigned int); ++i) {
    TEST_ASSERT(fcMap.find(indices[i]) != fcMap.end());
    TEST_ASSERT(fcMap[indices[i]] == 1);
  }
  delete resMolSuppl;
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
  unsigned int indices[] = { 0, 2, 4, 6 };
  for (unsigned int i = 0;
    i < sizeof(indices) / sizeof(unsigned int); ++i) {
    TEST_ASSERT(fcMap.find(indices[i]) != fcMap.end());
    TEST_ASSERT(fcMap[indices[i]] == -1);
  }
  delete resMolSuppl;
}

void testButadiene() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
    << "testButadiene" << std::endl;
  RWMol *mol = SmilesToMol("C=CC=C");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 1);
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_INCOMPLETE_OCTETS
    | ResonanceMolSupplier::UNCONSTRAINED_CATIONS
    | ResonanceMolSupplier::UNCONSTRAINED_ANIONS
    | ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 3);
  std::map<unsigned int, int> fcMap;
  while (!resMolSuppl->atEnd()) {
    ROMol *resMol = resMolSuppl->next();
    addFormalChargeIndices(resMol, fcMap);
    delete resMol;
  }
  unsigned int indices[] = { 0, 3 };
  for (unsigned int i = 0;
    i < sizeof(indices) / sizeof(unsigned int); ++i) {
    TEST_ASSERT(fcMap.find(indices[i]) != fcMap.end());
    TEST_ASSERT(fcMap[indices[i]] == 0);
  }
  delete resMolSuppl;
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
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 2);
  mol0 = (*resMolSuppl)[0];
  ROMol *mol1 = (*resMolSuppl)[1];
  TEST_ASSERT(mol1->getAtomWithIdx(0)->getFormalCharge() == 1);
  TEST_ASSERT(mol1->getAtomWithIdx(6)->getFormalCharge() == -1);
  delete mol1;
  delete resMolSuppl;
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
  unsigned int indices[] = { 0, 2, 4, 6 };
  for (unsigned int i = 0;
    i < sizeof(indices) / sizeof(unsigned int); ++i) {
    TEST_ASSERT(fcMap.find(indices[i]) != fcMap.end());
    TEST_ASSERT(fcMap[indices[i]] == 1);
  }
  delete resMolSuppl;
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
  unsigned int indices[] = { 0, 8, 10 };
  for (unsigned int i = 0;
    i < sizeof(indices) / sizeof(unsigned int); ++i) {
    TEST_ASSERT(fcMap.find(indices[i]) != fcMap.end());
    TEST_ASSERT(fcMap[indices[i]] == 1);
  }
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 4);
  resMolSuppl->moveTo(3);
  ROMol *resMol = resMolSuppl->next();
  TEST_ASSERT(resMol->getAtomWithIdx(0)->getFormalCharge() == 1);
  TEST_ASSERT(resMol->getAtomWithIdx(10)->getFormalCharge() == 1);
  TEST_ASSERT(resMol->getAtomWithIdx(8)->getFormalCharge() == -1);
  delete resMol;
  delete resMolSuppl;
}
  
void testChargeSeparation2() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
    << "testChargeSeparation2" << std::endl;
  RWMol *mol = SmilesToMol("NC(C(=O)[O-])=CC=CC=O");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 2);
  delete resMolSuppl;
  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 8);
  // less charge separation in the first 2,
  // more charge separation in the last 6
  // the last 4 feature carbanions
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
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 6)
        && ((*it)->getFormalCharge() == -1))
        haveCarbanion = true;
    }
    TEST_ASSERT(!haveCarbanion);
    delete resMol;
  }
  for (unsigned int i = 4; i < 8; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 6)
        && ((*it)->getFormalCharge() == -1))
        haveCarbanion = true;
    }
    TEST_ASSERT(haveCarbanion);
    delete resMol;
  }
  delete resMolSuppl;
}
  
void testMultipleConjGroups1() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
    << "testMultipleConjGroups1" << std::endl;
  RWMol *mol = SmilesToMol("NC(C(=O)[O-])=CC=C(CCC(=O)[O-])C=O");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION);
  TEST_ASSERT(resMolSuppl->length() == 16);
  for (unsigned int i = 0; i < 8; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 6)
        && ((*it)->getFormalCharge() == -1))
        haveCarbanion = true;
    }
    TEST_ASSERT(!haveCarbanion);
    delete resMol;
  }
  for (unsigned int i = 8; i < 16; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 6)
        && ((*it)->getFormalCharge() == -1))
        haveCarbanion = true;
    }
    TEST_ASSERT(haveCarbanion);
    delete resMol;
  }
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION, 8);
  TEST_ASSERT(resMolSuppl->length() == 8);
  for (unsigned int i = 0; i < 8; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCarbanion = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 6)
        && ((*it)->getFormalCharge() == -1))
        haveCarbanion = true;
    }
    TEST_ASSERT(!haveCarbanion);
    delete resMol;
  }
  delete resMolSuppl;
}

void testMultipleConjGroups2() {
  // if breadth-first resonance structure enumeration works properly
  // all N(2-) should be at the end of the supplier
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
    << "testMultipleConjGroups2" << std::endl;
  RWMol *mol = SmilesToMol("[N-]=[N+]=NCN=[N+]=[N-]");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 9);
  bool haveFirstDoubleNeg = false;
  for (unsigned int i = 0; i < resMolSuppl->length(); ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveDoubleNeg = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 7)
        && ((*it)->getFormalCharge() == -2))
        haveDoubleNeg = true;
    }
    if (!haveFirstDoubleNeg && haveDoubleNeg)
      haveFirstDoubleNeg = true;
    TEST_ASSERT(!haveFirstDoubleNeg
      || (haveFirstDoubleNeg && haveDoubleNeg));
    delete resMol;
  }
  delete resMolSuppl;
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
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 6)
        && ((*it)->getFormalCharge() == -1))
        haveCarbanion = true;
    }
    TEST_ASSERT(((i < 2) && !haveCarbanion)
      || ((i == 2) && haveCarbanion));
    delete resMol;
  }
  delete resMolSuppl;
}

void testMethylAcetate() {
  // cation on oxygen should appear only when UNCONSTRAINED_CATIONS
  // and ALLOW_CHARGE_SEPARATION are set
  BOOST_LOG(rdInfoLog) << "-----------------------\n"
    << "testMethylAcetate" << std::endl;
  RWMol *mol = SmilesToMol("CC(=O)OC");
  ResonanceMolSupplier *resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol);
  TEST_ASSERT(resMolSuppl->length() == 1);
  const ROMol *resMol = (*resMolSuppl)[0];
  bool haveCation = false;
  for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
    it != resMol->endAtoms(); ++it)  {
    if (((*it)->getAtomicNum() == 8)
      && ((*it)->getFormalCharge() == 1))
      haveCation = true;
  }
  delete resMol;
  TEST_ASSERT(!haveCation);
  delete resMolSuppl;

  resMolSuppl = new ResonanceMolSupplier((ROMol &)*mol,
    ResonanceMolSupplier::ALLOW_CHARGE_SEPARATION
    | ResonanceMolSupplier::UNCONSTRAINED_CATIONS);
  TEST_ASSERT(resMolSuppl->length() == 2);
  for (unsigned int i = 0; i < 2; ++i) {
    const ROMol *resMol = (*resMolSuppl)[i];
    bool haveCation = false;
    for (ROMol::ConstAtomIterator it = resMol->beginAtoms();
      it != resMol->endAtoms(); ++it)  {
      if (((*it)->getAtomicNum() == 8)
        && ((*it)->getFormalCharge() == 1))
        haveCation = true;
    }
    TEST_ASSERT(((i == 0) && !haveCation)
      || ((i == 1) && haveCation));
    delete resMol;
  }
  delete resMolSuppl;
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
  n = SubstructMatch(*mol, *query, matchVect,
		false, true, false, false);
  TEST_ASSERT(n == 1);
  matchVect.clear();
  n = SubstructMatch(*resMolSuppl, *query, matchVect,
		false, true, false, false);
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
  n = SubstructMatch(*mol, *query, matchVect,
		false, true, false, false);
  TEST_ASSERT((n == 1) && (matchVect.size() == 1));
  p = matchVect[0];
  TEST_ASSERT((p.size() == 1) && (p[0].second == 6));
  matchVect.clear();
  n = SubstructMatch(*resMolSuppl, *query, matchVect,
		false, true, false, false);
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
#endif
  return 0;
}
