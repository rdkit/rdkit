//
//  Copyright (C) 2001-2018 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/test.h>
#include "AlignMolecules.h"
#include "O3AAlignMolecules.h"
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Vector.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

using namespace RDKit;

void test1MolAlign() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);

  double rmsd = MolAlign::alignMol(*m2, *m1);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

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
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // specify conformations
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

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
  ROMol *m1 = supplier[1];
  ROMol *m2 = supplier[2];

  // alignMol() would return this for the rms: 2.50561
  // But the best rms is: 2.43449
  double rmsd = MolAlign::getBestRMS(*m1, *m2);

  TEST_ASSERT(RDKit::feq(rmsd, 2.43449));
  delete m1;
  delete m2;
}

void test1MolWithQueryAlign() {
  // identical to test1MolAlign except we replace one atom with a QueryAtom
  // instead

  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  auto *m1 = new RWMol(*MolFileToMol(fname1));
  auto *a1 = new QueryAtom(6);
  std::string fname2 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  auto *m2 = new RWMol(*MolFileToMol(fname2));
  auto *a2 = new QueryAtom(6);

  // we replace the same nitrogen instead with a null
  // query  28 and 19 are the "same" atoms
  m1->replaceAtom(28, a1);
  m2->replaceAtom(19, a2);

  double rmsd = MolAlign::alignMol(*m2, *m1);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  std::string fname3 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_trans.mol";

  auto *m3 = new RWMol(*MolFileToMol(fname3));
  m3->replaceAtom(0, new QueryAtom(5));

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
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // specify conformations
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  TEST_ASSERT(RDKit::feq(rmsd, 0.6578));

  // provide an atom mapping
  delete m1;
  delete m2;
  delete m3;
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
  TEST_ASSERT(RDKit::feq(rmsd, 0.0));
}

void testMMFFO3A() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2";
  std::string newSdf = sdf + "_MMFFO3A.sdf";
  sdf += ".sdf";
  SDMolSupplier supplier(sdf, true, false);
  int nMol = supplier.length();
  const int refNum = 48;
  // SDWriter *newMol = new SDWriter(newSdf);
  ROMol *refMol = supplier[refNum];
  MMFF::MMFFMolProperties refMP(*refMol);
  double cumScore = 0.0;
  double cumMsd = 0.0;
  for (int prbNum = 0; prbNum < nMol; ++prbNum) {
    ROMol *prbMol = supplier[prbNum];
    MMFF::MMFFMolProperties prbMP(*prbMol);
    MolAlign::O3A o3a(*prbMol, *refMol, &prbMP, &refMP);
    double rmsd = o3a.align();
    cumScore += o3a.score();
    cumMsd += rmsd * rmsd;
    // newMol->write(prbMol);
    delete prbMol;
  }
  cumMsd /= (double)nMol;
  delete refMol;
  // newMol->close();
  // std::cerr<<cumScore<<","<<sqrt(cumMsd)<<std::endl;
  TEST_ASSERT(RDKit::feq(cumScore, 6941.8, 1));
  TEST_ASSERT(RDKit::feq(sqrt(cumMsd), .345, .001));
}

void testCrippenO3A() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2";
  std::string newSdf = sdf + "_CrippenO3A.sdf";
  sdf += ".sdf";
  SDMolSupplier supplier(sdf, true, false);
  int nMol = supplier.length();
  const int refNum = 48;
  // SDWriter *newMol = new SDWriter(newSdf);
  ROMol *refMol = supplier[refNum];
  unsigned int refNAtoms = refMol->getNumAtoms();
  std::vector<double> refLogpContribs(refNAtoms);
  std::vector<double> refMRContribs(refNAtoms);
  std::vector<unsigned int> refAtomTypes(refNAtoms);
  std::vector<std::string> refAtomTypeLabels(refNAtoms);
  Descriptors::getCrippenAtomContribs(*refMol, refLogpContribs, refMRContribs,
                                      true, &refAtomTypes, &refAtomTypeLabels);
  double cumScore = 0.0;
  double cumMsd = 0.0;
  for (int prbNum = 0; prbNum < nMol; ++prbNum) {
    ROMol *prbMol = supplier[prbNum];
    unsigned int prbNAtoms = prbMol->getNumAtoms();
    std::vector<double> prbLogpContribs(prbNAtoms);
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(*prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
    MolAlign::O3A o3a(*prbMol, *refMol, &prbLogpContribs, &refLogpContribs,
                      MolAlign::O3A::CRIPPEN);
    double rmsd = o3a.align();
    cumScore += o3a.score();
    cumMsd += rmsd * rmsd;
    // newMol->write(prbMol);
    delete prbMol;
  }
  cumMsd /= (double)nMol;
  delete refMol;
  // newMol->close();
  // std::cerr<<cumScore<<","<<sqrt(cumMsd)<<std::endl;
  TEST_ASSERT(RDKit::feq(cumScore, 4918.1, 1));
  TEST_ASSERT(RDKit::feq(sqrt(cumMsd), .304, .001));
}

void testMMFFO3AMolHist() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2";
  std::string newSdf = sdf + "_MMFFO3A.sdf";
  sdf += ".sdf";
  SDMolSupplier supplier(sdf, true, false);
  int nMol = supplier.length();
  const int refNum = 48;
  // SDWriter *newMol = new SDWriter(newSdf);
  ROMol *refMol = supplier[refNum];
  MMFF::MMFFMolProperties refMP(*refMol);
  double *refDmat = MolOps::get3DDistanceMat(*refMol);
  MolAlign::MolHistogram refHist(*refMol, refDmat);
  double cumScore = 0.0;
  double cumMsd = 0.0;
  for (int prbNum = 0; prbNum < nMol; ++prbNum) {
    ROMol *prbMol = supplier[prbNum];
    MMFF::MMFFMolProperties prbMP(*prbMol);
    double *prbDmat = MolOps::get3DDistanceMat(*prbMol);
    MolAlign::MolHistogram prbHist(*prbMol, prbDmat);

    MolAlign::O3A o3a(*prbMol, *refMol, &prbMP, &refMP, MolAlign::O3A::MMFF94,
                      -1, -1, false, 50, 0, nullptr, nullptr, nullptr, &prbHist,
                      &refHist);
    double rmsd = o3a.align();
    cumScore += o3a.score();
    cumMsd += rmsd * rmsd;
    // newMol->write(prbMol);
    delete prbMol;
  }
  cumMsd /= (double)nMol;
  delete refMol;
  // newMol->close();
  // std::cerr<<cumScore<<","<<sqrt(cumMsd)<<std::endl;
  TEST_ASSERT(RDKit::feq(cumScore, 6941.8, 1));
  TEST_ASSERT(RDKit::feq(sqrt(cumMsd), .345, .001));
}

void testCrippenO3AMolHist() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2";
  std::string newSdf = sdf + "_CrippenO3A.sdf";
  sdf += ".sdf";
  SDMolSupplier supplier(sdf, true, false);
  int nMol = supplier.length();
  const int refNum = 48;
  // SDWriter *newMol = new SDWriter(newSdf);
  ROMol *refMol = supplier[refNum];
  unsigned int refNAtoms = refMol->getNumAtoms();
  std::vector<double> refLogpContribs(refNAtoms);
  std::vector<double> refMRContribs(refNAtoms);
  std::vector<unsigned int> refAtomTypes(refNAtoms);
  std::vector<std::string> refAtomTypeLabels(refNAtoms);
  Descriptors::getCrippenAtomContribs(*refMol, refLogpContribs, refMRContribs,
                                      true, &refAtomTypes, &refAtomTypeLabels);
  double *refDmat = MolOps::get3DDistanceMat(*refMol);
  MolAlign::MolHistogram refHist(*refMol, refDmat);
  double cumScore = 0.0;
  double cumMsd = 0.0;
  for (int prbNum = 0; prbNum < nMol; ++prbNum) {
    ROMol *prbMol = supplier[prbNum];
    unsigned int prbNAtoms = prbMol->getNumAtoms();
    std::vector<double> prbLogpContribs(prbNAtoms);
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(*prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
    double *prbDmat = MolOps::get3DDistanceMat(*prbMol);
    MolAlign::MolHistogram prbHist(*prbMol, prbDmat);

    MolAlign::O3A o3a(*prbMol, *refMol, &prbLogpContribs, &refLogpContribs,
                      MolAlign::O3A::CRIPPEN, -1, -1, false, 50, 0, nullptr,
                      nullptr, nullptr, &prbHist, &refHist);
    double rmsd = o3a.align();
    cumScore += o3a.score();
    cumMsd += rmsd * rmsd;
    // newMol->write(prbMol);
    delete prbMol;
  }
  cumMsd /= (double)nMol;
  delete refMol;
  // newMol->close();
  // std::cerr<<cumScore<<","<<sqrt(cumMsd)<<std::endl;
  TEST_ASSERT(RDKit::feq(cumScore, 4918.1, 1));
  TEST_ASSERT(RDKit::feq(sqrt(cumMsd), .304, .001));
}

void testMMFFO3AConstraints() {
  ROMol *m = SmilesToMol("n1ccc(cc1)-c1ccccc1");
  TEST_ASSERT(m);
  ROMol *m1 = MolOps::addHs(*m);
  delete m;
  TEST_ASSERT(m1);
  DGeomHelpers::EmbedMolecule(*m1);
  MMFF::sanitizeMMFFMol((RWMol &)(*m1));
  MMFF::MMFFMolProperties mp(*m1);
  TEST_ASSERT(mp.isValid());
  ForceFields::ForceField *field = MMFF::constructForceField(*m1, &mp);
  field->initialize();
  field->minimize();
  RWMol *patt = SmartsToMol("nccc-cccc");
  MatchVectType matchVect;
  TEST_ASSERT(SubstructMatch(*m1, (ROMol &)*patt, matchVect));
  unsigned int nIdx = matchVect[0].second;
  unsigned int cIdx = matchVect[matchVect.size() - 1].second;
  MolTransforms::setDihedralDeg(m1->getConformer(), matchVect[2].second,
                                matchVect[3].second, matchVect[4].second,
                                matchVect[5].second, 0.0);
  ROMol m2(*m1);
  MolAlign::randomTransform(m2);
  ROMol m3(m2);
  auto *o3a = new MolAlign::O3A(m2, *m1, &mp, &mp);
  TEST_ASSERT(o3a);
  o3a->align();
  delete o3a;
  double d =
      (m2.getConformer().getAtomPos(cIdx) - m1->getConformer().getAtomPos(cIdx))
          .length();
  TEST_ASSERT(feq(d, 0.0, 1));
  MatchVectType constraintMap;
  constraintMap.push_back(std::make_pair(cIdx, nIdx));
  o3a = new MolAlign::O3A(m3, *m1, &mp, &mp, MolAlign::O3A::MMFF94, -1, -1,
                          false, 50, 0, &constraintMap);
  TEST_ASSERT(o3a);
  o3a->align();
  delete o3a;
  d = (m3.getConformer().getAtomPos(cIdx) - m1->getConformer().getAtomPos(cIdx))
          .length();
  TEST_ASSERT(feq(d, 7.0, 1.0));
  delete m1;
}

void testCrippenO3AConstraints() {
  ROMol *m = SmilesToMol("n1ccc(cc1)-c1ccccc1");
  TEST_ASSERT(m);
  ROMol *m1 = MolOps::addHs(*m);
  delete m;
  TEST_ASSERT(m1);
  DGeomHelpers::EmbedMolecule(*m1);
  MMFF::sanitizeMMFFMol((RWMol &)(*m1));
  MMFF::MMFFMolProperties mp(*m1);
  TEST_ASSERT(mp.isValid());
  ForceFields::ForceField *field = MMFF::constructForceField(*m1, &mp);
  field->initialize();
  field->minimize();
  RWMol *patt = SmartsToMol("nccc-cccc");
  MatchVectType matchVect;
  TEST_ASSERT(SubstructMatch(*m1, (ROMol &)*patt, matchVect));
  unsigned int nIdx = matchVect[0].second;
  unsigned int cIdx = matchVect[matchVect.size() - 1].second;
  MolTransforms::setDihedralDeg(m1->getConformer(), matchVect[2].second,
                                matchVect[3].second, matchVect[4].second,
                                matchVect[5].second, 0.0);
  ROMol m2(*m1);
  MolAlign::randomTransform(m2);
  ROMol m3(m2);
  unsigned int prbNAtoms = m2.getNumAtoms();
  std::vector<double> prbLogpContribs(prbNAtoms);
  std::vector<double> prbMRContribs(prbNAtoms);
  std::vector<unsigned int> prbAtomTypes(prbNAtoms);
  std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
  Descriptors::getCrippenAtomContribs(m2, prbLogpContribs, prbMRContribs, true,
                                      &prbAtomTypes, &prbAtomTypeLabels);
  auto *o3a = new MolAlign::O3A(m2, *m1, &prbLogpContribs, &prbLogpContribs,
                                MolAlign::O3A::CRIPPEN);
  TEST_ASSERT(o3a);
  o3a->align();
  delete o3a;
  double d =
      (m2.getConformer().getAtomPos(cIdx) - m1->getConformer().getAtomPos(cIdx))
          .length();
  TEST_ASSERT(feq(d, 0.0, 1));
  MatchVectType constraintMap;
  constraintMap.push_back(std::make_pair(cIdx, nIdx));
  o3a = new MolAlign::O3A(m3, *m1, &prbLogpContribs, &prbLogpContribs,
                          MolAlign::O3A::CRIPPEN, -1, -1, false, 50, 0,
                          &constraintMap);
  TEST_ASSERT(o3a);
  o3a->align();
  delete o3a;
  d = (m3.getConformer().getAtomPos(cIdx) - m1->getConformer().getAtomPos(cIdx))
          .length();
  TEST_ASSERT(feq(d, 7.0, 1.0));
  delete m1;
}

void testMMFFO3AConstraintsAndLocalOnly() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2.sdf";
  SDMolSupplier supplier(sdf, true, false);
  const int refNum = 23;
  const int prbNum = 32;
  ROMol *refMol = supplier[refNum];
  ROMol *prbMol = supplier[prbNum];
  unsigned int refNAtoms = refMol->getNumAtoms();
  std::vector<double> refLogpContribs(refNAtoms);
  std::vector<double> refMRContribs(refNAtoms);
  std::vector<unsigned int> refAtomTypes(refNAtoms);
  std::vector<std::string> refAtomTypeLabels(refNAtoms);
  Descriptors::getCrippenAtomContribs(*refMol, refLogpContribs, refMRContribs,
                                      true, &refAtomTypes, &refAtomTypeLabels);
  unsigned int prbNAtoms = prbMol->getNumAtoms();
  std::vector<double> prbLogpContribs(prbNAtoms);
  std::vector<double> prbMRContribs(prbNAtoms);
  std::vector<unsigned int> prbAtomTypes(prbNAtoms);
  std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
  Descriptors::getCrippenAtomContribs(*prbMol, prbLogpContribs, prbMRContribs,
                                      true, &prbAtomTypes, &prbAtomTypeLabels);
  RWMol *patt = SmartsToMol("S");
  MatchVectType matchVect;
  TEST_ASSERT(SubstructMatch(*refMol, (ROMol &)*patt, matchVect));
  delete patt;
  unsigned int refSIdx = matchVect[0].second;
  matchVect.clear();
  patt = SmartsToMol("O");
  TEST_ASSERT(SubstructMatch(*prbMol, (ROMol &)*patt, matchVect));
  delete patt;
  unsigned int prbOIdx = matchVect[0].second;
  std::vector<double> distOS(2);
  distOS[0] = 2.7;
  distOS[1] = 0.4;
  std::vector<double> weights(2);
  weights[0] = 0.1;
  weights[1] = 100.0;
  for (unsigned int i = 0; i < 2; ++i) {
    MatchVectType constraintMap;
    constraintMap.push_back(std::make_pair(prbOIdx, refSIdx));
    RDNumeric::DoubleVector constraintWeights(1);
    constraintWeights[0] = weights[i];
    auto *o3a =
        new MolAlign::O3A(*prbMol, *refMol, &prbLogpContribs, &refLogpContribs,
                          MolAlign::O3A::CRIPPEN, -1, -1, false, 50, 0,
                          &constraintMap, &constraintWeights);
    TEST_ASSERT(o3a);
    o3a->align();
    delete o3a;
    o3a = new MolAlign::O3A(*prbMol, *refMol, &prbLogpContribs,
                            &refLogpContribs, MolAlign::O3A::CRIPPEN, -1, -1,
                            false, 50, MolAlign::O3_LOCAL_ONLY);
    TEST_ASSERT(o3a);
    o3a->align();
    delete o3a;
    double d = (prbMol->getConformer().getAtomPos(prbOIdx) -
                refMol->getConformer().getAtomPos(refSIdx))
                   .length();
    TEST_ASSERT(feq(d, distOS[i], 0.1));
  }
}

void testCrippenO3AConstraintsAndLocalOnly() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2.sdf";
  SDMolSupplier supplier(sdf, true, false);
  const int refNum = 23;
  const int prbNum = 32;
  ROMol *refMol = supplier[refNum];
  ROMol *prbMol = supplier[prbNum];
  MMFF::MMFFMolProperties refMP(*refMol);
  TEST_ASSERT(refMP.isValid());
  MMFF::MMFFMolProperties prbMP(*prbMol);
  TEST_ASSERT(prbMP.isValid());
  RWMol *patt = SmartsToMol("S");
  MatchVectType matchVect;
  TEST_ASSERT(SubstructMatch(*refMol, (ROMol &)*patt, matchVect));
  delete patt;
  unsigned int refSIdx = matchVect[0].second;
  matchVect.clear();
  patt = SmartsToMol("O");
  TEST_ASSERT(SubstructMatch(*prbMol, (ROMol &)*patt, matchVect));
  delete patt;
  unsigned int prbOIdx = matchVect[0].second;
  std::vector<double> distOS(2);
  distOS[0] = 3.2;
  distOS[1] = 0.3;
  std::vector<double> weights(2);
  weights[0] = 10.0;
  weights[1] = 100.0;
  for (unsigned int i = 0; i < 2; ++i) {
    MatchVectType constraintMap;
    constraintMap.push_back(std::make_pair(prbOIdx, refSIdx));
    RDNumeric::DoubleVector constraintWeights(1);
    constraintWeights[0] = weights[i];
    auto *o3a = new MolAlign::O3A(*prbMol, *refMol, &prbMP, &refMP,
                                  MolAlign::O3A::MMFF94, -1, -1, false, 50, 0,
                                  &constraintMap, &constraintWeights);
    TEST_ASSERT(o3a);
    o3a->align();
    delete o3a;
    o3a = new MolAlign::O3A(*prbMol, *refMol, &prbMP, &refMP,
                            MolAlign::O3A::MMFF94, -1, -1, false, 50,
                            MolAlign::O3_LOCAL_ONLY);
    TEST_ASSERT(o3a);
    o3a->align();
    delete o3a;
    double d = (prbMol->getConformer().getAtomPos(prbOIdx) -
                refMol->getConformer().getAtomPos(refSIdx))
                   .length();
    TEST_ASSERT(feq(d, distOS[i], 0.1));
  }
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock_o3a_mmff(ROMol *refMol, const std::vector<ROMol *> &mols,
                       const std::vector<double> &rmsds,
                       const std::vector<double> &scores, unsigned int count,
                       unsigned int idx) {
  for (unsigned int rep = 0; rep < 30; ++rep) {
    MMFF::MMFFMolProperties refMP(*refMol);
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) continue;
      if (!(rep % 10)) {
        BOOST_LOG(rdErrorLog) << "Rep: " << rep << " Mol:" << i << std::endl;
      }
      ROMol prbMol(*mols[i]);
      MMFF::MMFFMolProperties prbMP(prbMol);
      MolAlign::O3A o3a(prbMol, *refMol, &prbMP, &refMP);
      double rmsd = o3a.align();
      double score = o3a.score();
      TEST_ASSERT(feq(rmsd, rmsds[i]));
      TEST_ASSERT(feq(score, scores[i]));
    }
  }
}
void runblock_o3a_crippen(ROMol *refMol, const std::vector<ROMol *> &mols,
                          const std::vector<double> &rmsds,
                          const std::vector<double> &scores, unsigned int count,
                          unsigned int idx) {
  ROMol refMolCopy(*refMol);
  for (unsigned int rep = 0; rep < 30; ++rep) {
    unsigned int refNAtoms = refMolCopy.getNumAtoms();
    std::vector<double> refLogpContribs(refNAtoms);
    std::vector<double> refMRContribs(refNAtoms);
    std::vector<unsigned int> refAtomTypes(refNAtoms);
    std::vector<std::string> refAtomTypeLabels(refNAtoms);
    Descriptors::getCrippenAtomContribs(refMolCopy, refLogpContribs,
                                        refMRContribs, true, &refAtomTypes,
                                        &refAtomTypeLabels);
    for (unsigned int i = 0; i < mols.size(); ++i) {
      if (i % count != idx) continue;
      if (!(rep % 10)) {
        BOOST_LOG(rdErrorLog) << "Rep: " << rep << " Mol:" << i << std::endl;
      }
      ROMol prbMol(*mols[i]);
      unsigned int prbNAtoms = prbMol.getNumAtoms();
      std::vector<double> prbLogpContribs(prbNAtoms);
      std::vector<double> prbMRContribs(prbNAtoms);
      std::vector<unsigned int> prbAtomTypes(prbNAtoms);
      std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
      Descriptors::getCrippenAtomContribs(prbMol, prbLogpContribs,
                                          prbMRContribs, true, &prbAtomTypes,
                                          &prbAtomTypeLabels);
      MolAlign::O3A o3a(prbMol, refMolCopy, &prbLogpContribs, &refLogpContribs,
                        MolAlign::O3A::CRIPPEN);
      double rmsd = o3a.align();
      double score = o3a.score();
      TEST_ASSERT(feq(rmsd, rmsds[i]));
      TEST_ASSERT(feq(score, scores[i]));
    }
  }
}
}
#include <thread>
#include <future>
void testMMFFO3AMultiThread() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2.sdf";

  SDMolSupplier suppl(sdf, true, false);

  std::vector<ROMol *> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) continue;
    mols.push_back(mol);
  }

  std::cerr << "generating reference data" << std::endl;
  std::vector<double> rmsds(mols.size(), 0.0);
  std::vector<double> scores(mols.size(), 0.0);
  const int refNum = 48;
  ROMol *refMol = mols[refNum];
  MMFF::MMFFMolProperties refMP(*refMol);

  for (unsigned int i = 0; i < mols.size(); ++i) {
    ROMol prbMol(*mols[i]);
    MMFF::MMFFMolProperties prbMP(prbMol);
    MolAlign::O3A o3a(prbMol, *refMol, &prbMP, &refMP);
    rmsds[i] = o3a.align();
    scores[i] = o3a.score();
  }

  std::vector<std::future<void>> tg;

  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(std::async(std::launch::async, runblock_o3a_mmff, refMol,
                               mols, rmsds, scores, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }

  BOOST_FOREACH (ROMol *mol, mols) { delete mol; }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testCrippenO3AMultiThread() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2.sdf";

  SDMolSupplier suppl(sdf, true, false);

  std::vector<ROMol *> mols;
  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) continue;
    mols.push_back(mol);
  }

  std::cerr << "generating reference data" << std::endl;
  std::vector<double> rmsds(mols.size(), 0.0);
  std::vector<double> scores(mols.size(), 0.0);
  const int refNum = 48;
  ROMol *refMol = mols[refNum];
  unsigned int refNAtoms = refMol->getNumAtoms();
  std::vector<double> refLogpContribs(refNAtoms);
  std::vector<double> refMRContribs(refNAtoms);
  std::vector<unsigned int> refAtomTypes(refNAtoms);
  std::vector<std::string> refAtomTypeLabels(refNAtoms);
  Descriptors::getCrippenAtomContribs(*refMol, refLogpContribs, refMRContribs,
                                      true, &refAtomTypes, &refAtomTypeLabels);

  for (unsigned int i = 0; i < mols.size(); ++i) {
    ROMol prbMol(*mols[i]);
    unsigned int prbNAtoms = prbMol.getNumAtoms();
    std::vector<double> prbLogpContribs(prbNAtoms);
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
    MolAlign::O3A o3a(prbMol, *refMol, &prbLogpContribs, &refLogpContribs,
                      MolAlign::O3A::CRIPPEN);
    rmsds[i] = o3a.align();
    scores[i] = o3a.score();
  }

  std::vector<std::future<void>> tg;

  std::cerr << "processing" << std::endl;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(std::async(std::launch::async, runblock_o3a_crippen, refMol,
                               mols, rmsds, scores, count, i));
  }
  for (auto &fut : tg) {
    fut.get();
  }

  BOOST_FOREACH (ROMol *mol, mols) { delete mol; }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#endif

void testGetO3AForProbeConfs() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/ref_e2.sdf";

  SDMolSupplier suppl(sdf, true, false);
  ROMol *refMol = suppl[13];
  TEST_ASSERT(refMol);

  sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/probe_mol.sdf";
  SDMolSupplier psuppl(sdf, true, false);
  ROMol *prbMol = psuppl.next();
  TEST_ASSERT(prbMol);
  while (!psuppl.atEnd()) {
    ROMol *mol = psuppl.next();
    if (!mol) continue;
    auto *conf = new Conformer(mol->getConformer());
    prbMol->addConformer(conf, true);
    delete mol;
  }
  TEST_ASSERT(prbMol->getNumConformers() == 50);

  MMFF::MMFFMolProperties refMP(*refMol);
  MMFF::MMFFMolProperties prbMP(*prbMol);

  std::vector<std::pair<double, double>> oscores;
  for (unsigned int i = 0; i < prbMol->getNumConformers(); ++i) {
    MolAlign::O3A o3a(*prbMol, *refMol, &prbMP, &refMP, MolAlign::O3A::MMFF94,
                      i);
    double rmsd = o3a.align();
    double score = o3a.score();
    oscores.push_back(std::make_pair(rmsd, score));
  }

  {
    std::vector<boost::shared_ptr<MolAlign::O3A>> o3s;
    MolAlign::getO3AForProbeConfs(*prbMol, *refMol, &prbMP, &refMP, o3s);
    TEST_ASSERT(o3s.size() == prbMol->getNumConformers());
    for (unsigned int i = 0; i < prbMol->getNumConformers(); ++i) {
      TEST_ASSERT(feq(oscores[i].first, o3s[i]->align()));
      TEST_ASSERT(feq(oscores[i].second, o3s[i]->score()));
    }
  }
#ifdef RDK_TEST_MULTITHREADED
  {
    ROMol prbMol2(*prbMol);
    unsigned int nDups = 10;
    for (unsigned int j = 0; j < nDups; ++j) {
      for (unsigned int i = 0; i < prbMol->getNumConformers(); ++i) {
        prbMol2.addConformer(new Conformer(prbMol->getConformer(i)), true);
      }
    }

    std::vector<boost::shared_ptr<MolAlign::O3A>> o3s;
    MolAlign::getO3AForProbeConfs(prbMol2, *refMol, &prbMP, &refMP, o3s, 4);
    TEST_ASSERT(o3s.size() == prbMol2.getNumConformers());
    for (unsigned int i = 0; i < prbMol2.getNumConformers(); ++i) {
      TEST_ASSERT(
          feq(oscores[i % prbMol->getNumConformers()].first, o3s[i]->align()));
      TEST_ASSERT(
          feq(oscores[i % prbMol->getNumConformers()].second, o3s[i]->score()));
    }
  }

#endif
  delete refMol;
  delete prbMol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testO3AMultiThreadBug() {
  std::string rdbase = getenv("RDBASE");
  std::string sdf = rdbase + "/Code/GraphMol/MolAlign/test_data/bzr_data.sdf";

  SDMolSupplier suppl(sdf, true, false);

  std::vector<ROMol *> mols;
  while (!suppl.atEnd()) {
    ROMol *mol = suppl.next();
    if (!mol) continue;

    while (mol->getNumConformers() < 50) {
      auto *conf = new Conformer(mol->getConformer(0));
      mol->addConformer(conf, true);
    }
    mols.push_back(mol);
  }
  TEST_ASSERT(mols.size() == 10);

  auto *refMol = new ROMol(*mols[0]);
  TEST_ASSERT(refMol);

  MMFF::MMFFMolProperties refMP(*refMol);

#ifdef RDK_TEST_MULTITHREADED
  {
    for (auto &mol : mols) {
      ROMol prbMol = *mol;
      TEST_ASSERT(prbMol.getNumConformers() == 50);

      MMFF::MMFFMolProperties prbMP(prbMol);

      std::vector<std::pair<double, double>> oscores;
      for (unsigned int i = 0; i < prbMol.getNumConformers(); ++i) {
        MolAlign::O3A o3a(prbMol, *refMol, &prbMP, &refMP,
                          MolAlign::O3A::MMFF94, i);
        double rmsd = o3a.align();
        double score = o3a.score();
        oscores.push_back(std::make_pair(rmsd, score));
      }

      ROMol prbMol2 = *mol;
      std::vector<boost::shared_ptr<MolAlign::O3A>> o3s;
      MolAlign::getO3AForProbeConfs(prbMol2, *refMol, &prbMP, &refMP, o3s, 0);
      TEST_ASSERT(o3s.size() == prbMol2.getNumConformers());
      for (unsigned int i = 0; i < prbMol2.getNumConformers(); ++i) {
        TEST_ASSERT(
            feq(oscores[i % prbMol.getNumConformers()].first, o3s[i]->align()));
        TEST_ASSERT(feq(oscores[i % prbMol.getNumConformers()].second,
                        o3s[i]->score()));
      }
    }
  }

#endif
  delete refMol;
  for (auto &mol : mols) delete mol;

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
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

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testMMFFO3A \n\n";
  testMMFFO3A();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testMMFFO3A with pre-computed dmat and MolHistogram\n\n";
  testMMFFO3AMolHist();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testMMFFO3A with constraints\n\n";
  testMMFFO3AConstraints();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testMMFFO3A with variable weight constraints followed by "
               "local-only optimization\n\n";
  testMMFFO3AConstraintsAndLocalOnly();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testCrippenO3A \n\n";
  testCrippenO3A();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testCrippenO3A with pre-computed dmat and MolHistogram\n\n";
  testCrippenO3AMolHist();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testCrippenO3A with constraints\n\n";
  testCrippenO3AConstraints();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testCrippenO3A with variable weight constraints followed by "
               "local-only optimization\n\n";
  testCrippenO3AConstraintsAndLocalOnly();

#ifdef RDK_TEST_MULTITHREADED
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testMMFFO3A multithreading\n\n";
  testMMFFO3AMultiThread();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test O3A multithreading bug\n\n";
  testO3AMultiThreadBug();
#endif

#ifdef RDK_TEST_MULTITHREADED
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testCrippenO3A multithreading\n\n";
  testCrippenO3AMultiThread();
#endif

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test getO3AForProbeConfs\n\n";
  testGetO3AForProbeConfs();
#endif

#ifdef RDK_TEST_MULTITHREADED
  std::cout << "\t---------------------------------\n";
  std::cout << "\t test O3A multithreading bug\n\n";
  testO3AMultiThreadBug();
#endif

  std::cout << "***********************************************************\n";
}
