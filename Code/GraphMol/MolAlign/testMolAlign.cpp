//
//  Copyright (C) 2001-2018 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include "AlignMolecules.h"
#include "O3AAlignMolecules.h"
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Vector.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>

using namespace RDKit;

void test1MolAlign() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);

  double rmsd = MolAlign::alignMol(*m2, *m1);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578) || RDKit::feq(rmsd, 1.0345));

  std::string fname3 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_trans.mol";
  ROMol *m3 = MolFileToMol(fname3);
  const Conformer &conf1 = m2->getConformer(0);
  const Conformer &conf2 = m3->getConformer(0);
  unsigned int i, nat = m3->getNumAtoms();
  for (i = 0; i < nat; i++) {
    RDGeom::Point3D pt1 = conf1.getAtomPos(i);
    RDGeom::Point3D pt2 = conf2.getAtomPos(i);
    TEST_ASSERT(RDKit::feq(pt1.x, pt2.x, 0.001));
    TEST_ASSERT(RDKit::feq(pt1.y, pt2.y, 0.001));
    TEST_ASSERT(RDKit::feq(pt1.z, pt2.z, 0.001));
  }
  RDGeom::Transform3D trans;
  rmsd = MolAlign::getAlignmentTransform(*m1, *m2, trans);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578) || RDKit::feq(rmsd, 1.0345));

  // specify conformations
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578) || RDKit::feq(rmsd, 1.0345));

  // provide an atom mapping
  delete m1;
  delete m2;
  delete m3;
}

void test1GetBestRMS() {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/MolAlign/test_data/probe_mol.sdf";
  SDMolSupplier supplier(fname, true, false);
  std::unique_ptr<ROMol> prb(supplier[1]);
  std::unique_ptr<ROMol> ref(supplier[2]);
  std::unique_ptr<ROMol> prbCopy1(new ROMol(*prb));
  std::unique_ptr<ROMol> prbCopy2(new ROMol(*prb));
  std::unique_ptr<ROMol> prbCopy3(new ROMol(*prb));
  RDGeom::Transform3D bestTrans;
  MatchVectType bestMatch;

  // alignMol() would return this for the rms: 2.50561
  // But the best rms is: 2.43449
  double rmsdInPlace = MolAlign::CalcRMS(*prbCopy1, *ref);
  TEST_ASSERT(RDKit::feq(rmsdInPlace, 2.6026));
  double rmsd = MolAlign::getBestRMS(*prb, *ref);
  TEST_ASSERT(RDKit::feq(rmsd, 2.43449));
  double rmsdCopy = MolAlign::getBestAlignmentTransform(*prbCopy1, *ref,
                                                        bestTrans, bestMatch);
  TEST_ASSERT(RDKit::feq(rmsd, rmsdCopy));
  TEST_ASSERT(bestMatch.size() == ref->getNumAtoms());

  SmilesParserParams params;
  params.removeHs = false;
  ROMOL_SPTR scaffold(SmilesToMol(
      "N1C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])C([H])([H])C1([H])[H]",
      params));
  MatchVectType scaffoldMatch;
  TEST_ASSERT(SubstructMatch(*ref, *scaffold, scaffoldMatch));
  boost::dynamic_bitset<> scaffoldIndices(ref->getNumAtoms());
  for (const auto &pair : scaffoldMatch) {
    scaffoldIndices.set(pair.second);
  }
  std::vector<MatchVectType> matches;
  TEST_ASSERT(SubstructMatch(*ref, *prb, matches, false));
  std::vector<MatchVectType> matchesPruned(matches.size());
  std::transform(matches.begin(), matches.end(), matchesPruned.begin(),
                 [&scaffoldIndices](const auto &match) {
                   MatchVectType matchPruned;
                   std::copy_if(match.begin(), match.end(),
                                std::back_inserter(matchPruned),
                                [&scaffoldIndices](const auto &pair) {
                                  return scaffoldIndices.test(pair.second);
                                });
                   return matchPruned;
                 });
  rmsdInPlace = MolAlign::CalcRMS(*prbCopy2, *ref, -1, -1, matchesPruned);
  TEST_ASSERT(RDKit::feq(rmsdInPlace, 2.5672));
  rmsd = MolAlign::getBestRMS(*prb, *ref, -1, -1, matchesPruned);
  TEST_ASSERT(RDKit::feq(rmsd, 1.14329));
  rmsdCopy = MolAlign::getBestAlignmentTransform(
      *prbCopy2, *ref, bestTrans, bestMatch, -1, -1, matchesPruned);
  TEST_ASSERT(RDKit::feq(rmsd, rmsdCopy));
  TEST_ASSERT(bestMatch.size() == scaffoldMatch.size());
  RDNumeric::DoubleVector weights(scaffoldIndices.size(), 1.0);
  for (unsigned int i = 0; i < scaffoldIndices.size(); ++i) {
    if (scaffoldIndices.test(i)) {
      weights.setVal(i, 100.0);
    }
  }
  rmsdInPlace =
      MolAlign::CalcRMS(*prbCopy3, *ref, -1, -1, matches, 1000, true, &weights);
  TEST_ASSERT(RDKit::feq(rmsdInPlace, 17.7959));
  rmsd =
      MolAlign::getBestRMS(*prb, *ref, -1, -1, matches, 1000, true, &weights);
  TEST_ASSERT(RDKit::feq(rmsd, 10.9681));
  rmsdCopy = MolAlign::getBestAlignmentTransform(*prbCopy3, *ref, bestTrans,
                                                 bestMatch, -1, -1, matches,
                                                 1000, true, &weights);
  TEST_ASSERT(RDKit::feq(rmsd, rmsdCopy));
  TEST_ASSERT(bestMatch.size() == ref->getNumAtoms());
}

void test1MolWithQueryAlign() {
  // identical to test1MolAlign except we replace one atom with a QueryAtom
  // instead

  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  auto *m1 = MolFileToMol(fname1);
  auto *a1 = new QueryAtom(6);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  auto *m2 = MolFileToMol(fname2);
  auto *a2 = new QueryAtom(6);

  // we replace the same nitrogen instead with a null
  // query  28 and 19 are the "same" atoms
  m1->replaceAtom(28, a1);
  m2->replaceAtom(19, a2);

  double rmsd = MolAlign::alignMol(*m2, *m1);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578) || RDKit::feq(rmsd, 1.0345));

  std::string fname3 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_trans.mol";

  auto *m3 = MolFileToMol(fname3);
  auto *a3 = new QueryAtom(5);
  m3->replaceAtom(0, a3);

  const Conformer &conf1 = m2->getConformer(0);
  const Conformer &conf2 = m3->getConformer(0);
  unsigned int i, nat = m3->getNumAtoms();
  for (i = 0; i < nat; i++) {
    RDGeom::Point3D pt1 = conf1.getAtomPos(i);
    RDGeom::Point3D pt2 = conf2.getAtomPos(i);
    TEST_ASSERT(RDKit::feq(pt1.x, pt2.x, 0.001));
    TEST_ASSERT(RDKit::feq(pt1.y, pt2.y, 0.001));
    TEST_ASSERT(RDKit::feq(pt1.z, pt2.z, 0.001));
  }

  RDGeom::Transform3D trans;
  rmsd = MolAlign::getAlignmentTransform(*m1, *m2, trans);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578) || RDKit::feq(rmsd, 1.0345));

  // specify conformations
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578) || RDKit::feq(rmsd, 1.0345));

  // provide an atom mapping
  delete m1;
  delete m2;
  delete m3;
  delete a1;
  delete a2;
  delete a3;
}

void test2AtomMap() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
  double rmsd = MolAlign::alignMol(*m2, *m1, 0, 0, &atomMap);
  TEST_ASSERT(RDKit::feq(rmsd, 0.8525));
  delete m1;
  delete m2;
}

void test3Weights() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));

  RDNumeric::DoubleVector wts(6);
  wts.setVal(0, 1.0);
  wts.setVal(1, 1.0);
  wts.setVal(2, 1.0);
  wts.setVal(3, 1.0);
  wts.setVal(4, 1.0);
  wts.setVal(5, 2.0);
  double rmsd = MolAlign::alignMol(*m2, *m1, 0, 0, &atomMap, &wts);
  TEST_ASSERT(RDKit::feq(rmsd, 0.9513));
  delete m1;
  delete m2;
}

void testIssue241() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/Issue241.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string res;
  MolPickler::pickleMol(*m1, res);
  auto *ref = new ROMol(res);
  DGeomHelpers::EmbedMolecule(*ref, 30, 239 * 10);
  ForceFields::ForceField *ff1 = UFF::constructForceField(*ref);
  ff1->initialize();
  ff1->minimize(200, 1e-8, 1e-6);

  std::string pkl2;
  MolPickler::pickleMol(*m1, pkl2);
  auto *probe = new ROMol(pkl2);
  DGeomHelpers::EmbedMolecule(*probe, 30, 239 * 10);
  ForceFields::ForceField *ff2 = UFF::constructForceField(*probe);
  ff2->initialize();
  ff2->minimize(200, 1e-8, 1e-6);

  double rmsd = MolAlign::alignMol(*ref, *probe);

  delete ff1;
  delete ff2;
  delete m1;
  delete ref;
  delete probe;

  TEST_ASSERT(RDKit::feq(rmsd, 0.0));
}

int main() {
  std::cout << "***********************************************************\n";
  std::cout << "Testing MolAlign\n";

#if 1
  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1MolAlign \n\n";
  test1MolAlign();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1GetBestRMS \n\n";
  test1GetBestRMS();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1MolWithQueryAlign \n\n";
  test1MolWithQueryAlign();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test2AtomMap \n\n";
  test2AtomMap();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test3Weights \n\n";
  test3Weights();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testIssue241 \n\n";
  testIssue241();

#endif

  std::cout << "***********************************************************\n";
}
