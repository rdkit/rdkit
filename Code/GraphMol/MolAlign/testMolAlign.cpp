// $Id$
//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AlignMolecules.h"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Vector.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>

using namespace RDKit;

void test1MolAlign() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  
  double rmsd = MolAlign::alignMol(*m2, *m1);
  CHECK_INVARIANT(RDKit::feq(rmsd, 0.6578), "");
  
  std::string fname3 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_trans.mol";
  ROMol *m3 = MolFileToMol(fname3);
  const Conformer &conf1 = m2->getConformer(0);
  const Conformer &conf2 = m3->getConformer(0);
  unsigned int i, nat = m3->getNumAtoms();
  for (i = 0; i < nat; i++) {
    RDGeom::Point3D pt1 = conf1.getAtomPos(i);
    RDGeom::Point3D pt2 = conf2.getAtomPos(i);
    CHECK_INVARIANT(RDKit::feq(pt1.x, pt2.x, 0.001), "");
    CHECK_INVARIANT(RDKit::feq(pt1.y, pt2.y, 0.001), "");
    CHECK_INVARIANT(RDKit::feq(pt1.z, pt2.z, 0.001), "");
  }
  RDGeom::Transform3D trans;
  rmsd = MolAlign::getAlignmentTransform(*m1, *m2, trans);
  CHECK_INVARIANT(RDKit::feq(rmsd, 0.6578), "");

  // specify conformations
  rmsd = MolAlign::alignMol(*m1, *m2, 0, 0);
  CHECK_INVARIANT(RDKit::feq(rmsd, 0.6578), "");

  // provide an atom mapping
  delete m1;
  delete m2;
  delete m3;
}

void test2AtomMap() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
  double rmsd = MolAlign::alignMol(*m2, *m1, 0, 0, &atomMap);
  CHECK_INVARIANT(RDKit::feq(rmsd, 0.8525), "");
  delete m1;
  delete m2;
  
}

void test3Weights() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string fname2 = rdbase + "/Code/GraphMol/MolAlign/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
  
  RDNumeric::DoubleVector wts(6);
  wts.setVal(0, 1.0); wts.setVal(1, 1.0);
  wts.setVal(2, 1.0); wts.setVal(3, 1.0);
  wts.setVal(4, 1.0); wts.setVal(5, 2.0);
  double rmsd = MolAlign::alignMol(*m2, *m1, 0, 0, &atomMap, &wts);
  CHECK_INVARIANT(RDKit::feq(rmsd, 0.9513), "");
  delete m1;
  delete m2;
}

void testIssue241() {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolAlign/test_data/Issue241.mol";
  ROMol *m1 = MolFileToMol(fname1);
  std::string res;
  MolPickler::pickleMol(*m1, res);
  ROMol *ref = new ROMol(res);
  DGeomHelpers::EmbedMolecule(*ref, 30, 239*10);
  ForceFields::ForceField *ff1=UFF::constructForceField(*ref);
  ff1->initialize();
  ff1->minimize(200, 1e-8, 1e-6);

  std::string pkl2;
  MolPickler::pickleMol(*m1, pkl2);
  ROMol *probe = new ROMol(pkl2);
  DGeomHelpers::EmbedMolecule(*probe, 30, 239*10);
  ForceFields::ForceField *ff2=UFF::constructForceField(*probe);
  ff2->initialize();
  ff2->minimize(200, 1e-8, 1e-6);

  double rmsd = MolAlign::alignMol(*ref, *probe);
  CHECK_INVARIANT(RDKit::feq(rmsd, 0.0), "");
}
    
int main() {
  std::cout << "***********************************************************\n";
  std::cout << "Testing MolAlign\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1MolAlign \n\n";
  test1MolAlign();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test2AtomMap \n\n";
  test2AtomMap();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test3Weights \n\n";
  test3Weights();
    
  std::cout << "\t---------------------------------\n";
  std::cout << "\t testIssue241 \n\n";
  testIssue241();
  std::cout << "***********************************************************\n";
  return(0);

}

