//  $Id$
// 
//   Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <Geometry/Transform3D.h>
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <Geometry/point.h>
#include "MolTransforms.h"

using namespace RDKit;
using namespace MolTransforms;
bool comparePts(const RDGeom::Point3D &pt1, const RDGeom::Point3D &pt2, double tol=1.0e-4) {
  RDGeom::Point3D tpt = pt1;
  tpt -= pt2;
  return (tpt.length() < tol);
}

void test1Canonicalization() {
  ROMol *mol = SmilesToMol("C", 0, 1);
  Conformer *conf = new Conformer(1);
  conf->setAtomPos(0, RDGeom::Point3D(4.0, 5.0, 6.0));
  int cid = mol->addConformer(conf, true);
  RDGeom::Point3D pt = computeCentroid(*conf);
  CHECK_INVARIANT(comparePts(pt, RDGeom::Point3D(4.0, 5.0, 6.0)), "");
  
  RDGeom::Transform3D *trans = computeCanonicalTransform(*conf);
  transformConformer(*conf, *trans);
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(0.0, 0.0, 0.0)), "");

  conf->setAtomPos(0, RDGeom::Point3D(4.0, 5.0, 6.0));
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(0.0, 0.0, 0.0)), "");

  delete mol;
  //delete conf;
  delete trans;
  // lets try two points now
  mol = SmilesToMol("CC", 0, 1);
  conf = new Conformer(2);
  conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
  cid = mol->addConformer(conf, true);
  trans = computeCanonicalTransform(*conf);
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");

  conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(0.0, 1.5, 0.0));
  trans = computeCanonicalTransform(*conf);
  canonicalizeConformer(*conf);
  
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");
  delete mol;
  delete trans;
  
  mol = SmilesToMol("CC", 0, 1);
  conf = new Conformer(2);
  conf->setAtomPos(0, RDGeom::Point3D(0.0, 0.0, 0.0));
  conf->setAtomPos(1, RDGeom::Point3D(1.5, 0.0, 0.0));
  cid = mol->addConformer(conf, true);
  trans = computeCanonicalTransform(*conf);
  transformConformer(*conf, *trans);
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");
  delete mol;
  delete trans;

  mol = SmilesToMol("C1CC1", 0, 1);
  conf = new Conformer(3);
  conf->setAtomPos(0, RDGeom::Point3D(0.58, -0.66, -0.08));
  conf->setAtomPos(1, RDGeom::Point3D(-0.88, -0.18, -0.04));
  conf->setAtomPos(2, RDGeom::Point3D(.26, 0.82, 0.14));
  cid = mol->addConformer(conf, true);
  //trans = computeCanonicalTransform(*conf);
  //transformConformer(*conf, *trans);
  canonicalizeConformer(*conf);
  CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.6418, 0.6158, 0.0)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(1), RDGeom::Point3D(-0.2029, -0.8602, 0.0)), "");
  CHECK_INVARIANT(comparePts(conf->getAtomPos(2), RDGeom::Point3D(0.8447, 0.2445, 0.0)), "");
  MolToMolFile(*mol, "junk.mol", 0);
  //CHECK_INVARIANT(comparePts(conf->getAtomPos(0), RDGeom::Point3D(-0.75, 0.0, 0.0)), "");
  //CHECK_INVARIANT(comparePts(conf->getAtomPos(1), RDGeom::Point3D(0.75, 0.0, 0.0)), "");
  delete mol;
  
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/MolTransforms/test_data/1oir.mol";
  mol = MolFileToMol(fname1);
  std::string fname2 = rdbase + "/Code/GraphMol/MolTransforms/test_data/1oir_canon.mol";
  ROMol *mol2 = MolFileToMol(fname2);
  
  Conformer &conf1 = mol->getConformer(0); 
  canonicalizeConformer(conf1);

  Conformer &conf2 = mol2->getConformer();
  unsigned int i, nats = mol->getNumAtoms();
  for (i = 0; i < nats; ++i) {
    CHECK_INVARIANT(comparePts(conf1.getAtomPos(i), conf2.getAtomPos(i)), "");
  }


  delete mol;
  delete mol2;
}

void test1(){
  std::cout << " ----------> Test1 "<< std::endl;

  std::cout << " Finished <---------- "<< std::endl;
}



int main() {
  //test1();
  std::cout << "***********************************************************\n";
  std::cout << "Testing MolTransforms\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1Canonicalization \n\n";
  test1Canonicalization();
  std::cout << "***********************************************************\n";
  return(0);
}
