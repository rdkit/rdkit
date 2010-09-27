// $Id$
// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <Geometry/UniformGrid3D.h>
#include "ShapeEncoder.h"
#include "ShapeUtils.h"
#include <GraphMol/RDKitBase.h>
//#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/FileParsers.h>
//#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
//#include <ForceField/ForceField.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <Geometry/GridUtils.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

void test1Encode() {
  RDGeom::UniformGrid3D grd(30.0, 16.0, 10.0);
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  ROMol *m = MolFileToMol(fname1);
  MolTransforms::canonicalizeMol(*m);
  
  MolShapes::EncodeShape(*m, grd, 0);
  //RDGeom::writeGridToFile(grd, "junk.grd");
  //MolToMolFile(m, "junk.mol", 0);
  CHECK_INVARIANT(grd.getOccupancyVect()->getTotalVal() == 9250, "");
}

void test2Compare() {
  RDGeom::UniformGrid3D grd(30.0, 16.0, 10.0);
  std::string rdbase = getenv("RDBASE");
  std::string fname1 = rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir.mol";
  ROMol *m = MolFileToMol(fname1);
  MolTransforms::canonicalizeMol(*m);

  ROMol *mdup = MolFileToMol(fname1);
  MolTransforms::canonicalizeMol(*mdup);

  double dist = MolShapes::tanimotoDistance(*m, *mdup);
  CHECK_INVARIANT(dist == 0.0, "");

  delete m;
  delete mdup;

  m = MolFileToMol(fname1);
  std::string fname2 = rdbase + "/Code/GraphMol/ShapeHelpers/test_data/1oir_conf.mol";
  ROMol *m2 = MolFileToMol(fname2);
  double rmsd = MolAlign::alignMol(*m, *m2);
  dist = MolShapes::tanimotoDistance(*m, *m2);
  CHECK_INVARIANT(RDKit::feq(dist, 0.2813), "");
  delete m2;
  
  m2 = MolFileToMol(fname2);
  MatchVectType atomMap;
  atomMap.push_back(std::pair<int, int>(18, 27));
  atomMap.push_back(std::pair<int, int>(13, 23));
  atomMap.push_back(std::pair<int, int>(21, 14));
  atomMap.push_back(std::pair<int, int>(24, 7));
  atomMap.push_back(std::pair<int, int>(9, 19));
  atomMap.push_back(std::pair<int, int>(16, 30));
  rmsd = MolAlign::alignMol(*m2, *m, 0, 0, &atomMap);
  dist = MolShapes::tanimotoDistance(*m, *m2);
  
  CHECK_INVARIANT(RDKit::feq(dist, 0.3244), "");
  delete m;
  delete m2;
}

void test3Methane() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/ShapeHelpers/test_data/methane.mol";
  ROMol *m = MolFileToMol(fname);
  RDGeom::Point3D dims, offSet;
  MolShapes::computeConfDimsAndOffset(m->getConformer(), dims, offSet, 0, 3.0);
  std::cout << dims << " " << offSet << "\n";
  RDGeom::UniformGrid3D grd(6.5, 6.5, 6.5);
  //dims.x, dims.y, dims.z, 0.5, DiscreteValueVect::TWOBITVALUE, &offSet);
  MolShapes::EncodeShape(*m, grd, 0);
  RDGeom::writeGridToFile(grd, "methane.grd");
}


int main() {

#if 1
  std::cout << "***********************************************************\n";
  std::cout << "Testing ShapeHelpers\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test1Encode \n\n";
  test1Encode();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t test2Compare \n\n";
  test2Compare();
  std::cout << "***********************************************************\n";
#endif
  //test3Methane();
  return 0;
}

