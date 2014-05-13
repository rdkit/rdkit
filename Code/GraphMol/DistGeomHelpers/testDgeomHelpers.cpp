// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/DistGeomUtils.h>
#include <DistGeom/TriangleSmooth.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <iostream>
#include "BoundsMatrixBuilder.h"
#include "Embedder.h"
#include <stdlib.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ROMol.h>
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <RDGeneral/FileParseException.h>
#include <ForceField/ForceField.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <math.h>
#include <RDBoost/Exceptions.h>

#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

using namespace RDKit;

void test1() {
  std::string smiString = "CC1=C(C(C)=CC=C2)C2=CC=C1 c1ccccc1C C/C=C/CC \
                           C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4 C1CCCCS1(=O)(=O) c1ccccc1 \
                           C1CCCC1 C1CCCCC1 \
                           C1CC1(C)C C12(C)CC1CC2";
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/initCoords.sdf";
  SDMolSupplier sdsup(fname);
  //SDWriter writer("foo.sdf");

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();++token){
    std::string smi= *token;
    RWMol *m = SmilesToMol(smi, 0, 1); 
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1,true,false,2,true,1,0,1e-2);
    CHECK_INVARIANT(cid >= 0, "");
    ROMol *m2 = sdsup.next();
    //BOOST_LOG(rdInfoLog) << ">>> " << smi << std::endl;
    //writer.write(*m);
    //writer.flush();

    //ROMol *m2 = NULL;
    if(m2){
      unsigned int nat = m->getNumAtoms();
    
      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
#if 0
      BOOST_LOG(rdInfoLog) << "-----------------------" << std::endl;
      BOOST_LOG(rdInfoLog) << MolToMolBlock(*m2) << std::endl;
      BOOST_LOG(rdInfoLog) << "---" << std::endl;
      BOOST_LOG(rdInfoLog) << MolToMolBlock(*m) << std::endl;
      BOOST_LOG(rdInfoLog) << "-----------------------" << std::endl;
#endif
      for (unsigned int i = 0; i < nat; i++) {
	RDGeom::Point3D pt1i = conf1.getAtomPos(i);
	RDGeom::Point3D pt2i = conf2.getAtomPos(i);
	for(unsigned int j=i+1;j<nat;j++){
	  RDGeom::Point3D pt1j = conf1.getAtomPos(j);
	  RDGeom::Point3D pt2j = conf2.getAtomPos(j);
	  double d1=(pt1j-pt1i).length();
	  double d2=(pt2j-pt2i).length();
	  if(m->getBondBetweenAtoms(i,j)){
            //BOOST_LOG(rdInfoLog) << ">1> " <<i<<","<<j<<":"<< d1 << " " << d2 << std::endl;
	    TEST_ASSERT(fabs(d1-d2)/d1<0.06);
          }else{
            //BOOST_LOG(rdInfoLog) << ">2> " <<i<<","<<j<<":"<< d1 << " " << d2 << " "<<fabs(d1-d2)/d1<<std::endl;
	    TEST_ASSERT(fabs(d1-d2)/d1<0.12);
	  }
	}
      }
    }
    delete m;
    delete m2;
  }
}

void computeDistMat(const RDGeom::PointPtrVect &origCoords, RDNumeric::DoubleSymmMatrix &distMat) {
  unsigned int N = origCoords.size();
  CHECK_INVARIANT(N == distMat.numRows(), "");
  unsigned int i, j;
  RDGeom::Point3D pti, ptj;
  double d;
  for (i = 1; i < N; i++) {
    pti = *(RDGeom::Point3D*)origCoords[i];
    for (j = 0; j < i; j++) {
      ptj = *(RDGeom::Point3D*)origCoords[j];
      ptj -= pti;
      d = ptj.length();
      distMat.setVal(i,j, d);
    }
  }
}

void computeMolDmat(ROMol &mol, RDNumeric::DoubleSymmMatrix &distMat) {
  RDGeom::PointPtrVect origCoords;
  unsigned int i, nat = mol.getNumAtoms();
  Conformer &conf = mol.getConformer(0);
  for (i = 0; i < nat; i++) {
    origCoords.push_back(&conf.getAtomPos(i));
  }
  computeDistMat(origCoords, distMat);
}

void test2() {
  // check for in ring distances, and distances containing two bonds in a ring
  std::string smi = "Cc1c(C=CC(C)N2)c2[nH]n1";
  ROMol *mol = SmilesToMol(smi, 0, 1);
  DistGeom::BoundsMatPtr bm;
  RDNumeric::DoubleSymmMatrix *dmat;
  int cid;
  unsigned int nat = mol->getNumAtoms();

#if 1
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);

  std::cerr<<"go"<<std::endl;
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);

  TEST_ASSERT( (bm->getUpperBound(0,9) - bm->getLowerBound(0,9)) < 0.13);
  TEST_ASSERT( (bm->getUpperBound(0,9) - dmat->getVal(0,9) > -0.1 )
               && (bm->getLowerBound(0,9) - dmat->getVal(0,9) < 0.10));

  TEST_ASSERT( (bm->getUpperBound(10,7) - bm->getLowerBound(10,7)) < 0.13);
  TEST_ASSERT( (bm->getUpperBound(10,7) - dmat->getVal(10,7) > -0.1 )
               && (bm->getLowerBound(10,7) - dmat->getVal(10,7) < 0.10 ));

  TEST_ASSERT( (bm->getUpperBound(2,5) - bm->getLowerBound(2,5)) < 0.20);
  TEST_ASSERT( (bm->getUpperBound(2,5) - dmat->getVal(2,5) > -0.1 )
               && (bm->getLowerBound(2,5) - dmat->getVal(2,5) < 0.10 ));

  TEST_ASSERT( (bm->getUpperBound(8,4) - bm->getLowerBound(8,4)) > 1.);
  TEST_ASSERT( (bm->getUpperBound(8,4) - bm->getLowerBound(8,4)) < 1.2);
  TEST_ASSERT( (bm->getUpperBound(8,4) - dmat->getVal(8,4) > -0.1 )
               && (bm->getLowerBound(8,4) - dmat->getVal(8,4) < 0.10));
  
  TEST_ASSERT( (bm->getUpperBound(8,6) - bm->getLowerBound(8,6)) > 1.0);
  TEST_ASSERT( (bm->getUpperBound(8,6) - bm->getLowerBound(8,6)) < 1.2);
  TEST_ASSERT( (bm->getUpperBound(8,6) - dmat->getVal(8,6) > -0.1)
               && (bm->getLowerBound(8,6) - dmat->getVal(8,6) < 0.10 ));
  
  delete mol;
  delete dmat;

  // chain of singles
  smi = "CCCC";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) > 1.0);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < 1.3);
  TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
               && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10 ));
  
  delete mol;
  delete dmat;

  smi = "C=C=C=C";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  TEST_ASSERT( bm->getUpperBound(0,3) > dmat->getVal(0,3));
  // this is kinda goofy but this linear molecule doesn't satisfy the bounds completely
  TEST_ASSERT(fabs(bm->getLowerBound(0,3) - dmat->getVal(0,3)) < 0.2);
  
  delete mol;
  delete dmat;
#endif
  std::cerr << "-------------------------------------\n\n";
  smi = "C/C=C/C";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  //std::cerr << "\n-----\n" << MolToMolBlock(mol,false,cid) << std::endl;;
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  //std::cerr << bm->getUpperBound(0,3) << " " << dmat->getVal(0,3) << " " << bm->getLowerBound(0,3) << std::endl;
  TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
               && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.1));
  
  delete mol;
  delete dmat;

  smi = "C/C=C\\C";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
               && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10));
  
  delete mol;
  delete dmat;

  smi = "CC=CC";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < 1.13);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) > 1.);
  TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
               && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10));
  
  delete mol;
  delete dmat;

  smi = "O=S-S=O";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
               && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10 ));
  
  delete mol;
  delete dmat;

#if 0
  // this next block of tests all handle the special case that led to Issue284
  smi = "COC=O";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  //TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  double x,y,z;
  x = bm->getUpperBound(0,3);
  y = bm->getLowerBound(0,3);
  z = dmat->getVal(0,3);
  //TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
  //             && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10 ));
  
  delete mol;
  delete dmat;

  smi = "C[NH]C=O";
  mol = SmilesToMol(smi, 0, 1);
  nat = mol->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
  DGeomHelpers::setTopolBounds(*mol, bm);
  cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1);
  TEST_ASSERT(cid>-1);
  dmat = new RDNumeric::DoubleSymmMatrix(nat, 0.0);
  computeMolDmat(*mol, *dmat);
  TEST_ASSERT( (bm->getUpperBound(0,3) - bm->getLowerBound(0,3)) < .13);
  //TEST_ASSERT( (bm->getUpperBound(0,3) - dmat->getVal(0,3) > -0.1)
  //             && (bm->getLowerBound(0,3) - dmat->getVal(0,3) < 0.10));
  
  delete mol;
  delete dmat;
#endif
}

void test3() {
  // check embedding based based on distances calculated from previous created
  // (good) coordinates 
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/combi_coords.sdf";
  SDMolSupplier sdsup(fname);
  
  unsigned int i, j, nat;
  bool gotCoords;
  while (!sdsup.atEnd()) {
    ROMol *mol = sdsup.next();
    std::string mname;
    mol->getProp("_Name", mname);
    RDGeom::PointPtrVect origCoords, newCoords;
    nat = mol->getNumAtoms();
    Conformer &conf = mol->getConformer(0);
    for (i = 0; i < nat; i++) {
      origCoords.push_back(&conf.getAtomPos(i));
    }
    RDNumeric::DoubleSymmMatrix distMat(nat, 0.0);
    computeDistMat(origCoords, distMat);
    
    gotCoords = DistGeom::computeInitialCoords(distMat, origCoords);
    
    CHECK_INVARIANT(gotCoords, "");
    RDNumeric::DoubleSymmMatrix distMatNew(nat, 0.0);
    computeDistMat(origCoords, distMatNew);
    
    for (i = 1; i < nat; i++) {
      for (j = 0; j < i; j++) {
        CHECK_INVARIANT(RDKit::feq(distMat.getVal(i,j), distMatNew.getVal(i,j), 0.01), "");
      }
    }
  }
}

void test4() {
  std::string smi = "c1cc(C(F)(F)F)ccc1/C=N/NC(=O)c(n2)c[n]3cc(C(F)(F)F)cc(c23)Cl";
  ROMol *m = SmilesToMol(smi, 0, 1);
  DGeomHelpers::EmbedMolecule(*m, 10, 1);//etCoords(*m, iter);
  std::string fname = "test.mol";
  MolToMolFile(*m, fname);
  delete m;
}

void test5() {
  // some real CDK2 molecules lets see how many fail
  std::string rdbase = getenv("RDBASE");
  std::string smifile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);
  
  ROMol *mol;
  int i = 0;
  int cid;
  while (1) {
    try {
      i++;
      mol = smiSup.next();
      cid = DGeomHelpers::EmbedMolecule(*mol, 10, 1); //getCoords(*mol, iter);
      TEST_ASSERT(cid>-1);
      delete mol;
    } catch (FileParseException &) {
      break;
    }
  }
}

void test6() {
  std::string smi;
  ROMol *m;
  DistGeom::BoundsMatPtr bm;

  m = SmilesToMol("CC");
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  TEST_ASSERT(bm);
  DGeomHelpers::initBoundsMat(bm,0.0,1000.0);
  TEST_ASSERT(feq(bm->getLowerBound(0,1),0.0));
  TEST_ASSERT(feq(bm->getLowerBound(1,0),0.0));
  TEST_ASSERT(feq(bm->getUpperBound(0,1),1000.0));
  TEST_ASSERT(feq(bm->getUpperBound(1,0),1000.0));
  
  DGeomHelpers::setTopolBounds(*m,bm);
  TEST_ASSERT(bm->getLowerBound(0,1)>0.0);
  TEST_ASSERT(bm->getUpperBound(0,1)<1000.0);
  TEST_ASSERT(bm->getLowerBound(0,1)<bm->getUpperBound(0,1));
  
  delete m;
}

void testIssue215() {
  std::string smi;
  ROMol *m;
  DistGeom::BoundsMatPtr bm;
  bool ok;

  m = SmilesToMol("C=C1C2CC1C2");
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  TEST_ASSERT(bm);
  DGeomHelpers::initBoundsMat(bm,0.0,1000.0);
  DGeomHelpers::setTopolBounds(*m,bm);

  // this was the specific problem:
  TEST_ASSERT(bm->getUpperBound(0,4)<100.0);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  
  delete m;
}

void _computeStats(const std::vector<double> &vals, double &mean,
                   double &stdDev) {
  unsigned int N = vals.size();
  mean = 0.0;
  double sqSum=0.0;
  stdDev = 0.0;
  std::vector<double>::const_iterator vci;
  for (vci = vals.begin(); vci != vals.end(); vci++) {
    mean += (*vci);
    sqSum += (*vci)*(*vci);
  }
  mean /= N;
  sqSum /= N;
  stdDev = sqrt(sqSum - (mean*mean));
}

  
void testTemp() {
  //ROMol *m = SmilesToMol("CN(C)c(cc1)ccc1\\C=C(/C#N)c(n2)c(C#N)c([n]3cccc3)[n]2c(c4)cccc4");
  //ROMol *m = SmilesToMol("C12C(=O)N(Cc3ccccc3)C(=O)C1C(C(=O)OC)(NC2c4ccc(Cl)c(c4)[N+]([O-])=O)Cc5c6ccccc6[nH]c5");
  //c1(n2cccc2)n(c3ccccc3)ncc1
  //ROMol *m = SmilesToMol("N#CCc1c(C#N)cn(C)n1");
  //ROMol *m = SmilesToMol("Cc1c(C#N)cn(C)n1");
  std::string rdbase = getenv("RDBASE");
  std::string smifile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);
  ROMol *m;
  unsigned int cnt = 0;
  while (!smiSup.atEnd()) {
    m = smiSup.next();
    unsigned int i;
    std::vector<double> energies;
    double energy; 
    for (i = 0; i < 20; i++) {
      int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, false);
      if (cid >=0) {
        ForceFields::ForceField *ff=UFF::constructForceField(*m, 10, cid);
        ff->initialize();
        energy = ff->calcEnergy();
        energies.push_back(energy);
        delete ff;
      }
    }
    double mean, stdDev;
    _computeStats(energies, mean, stdDev);
    std::string mname;
    cnt++;
    m->getProp("_Name", mname);
    BOOST_LOG(rdDebugLog) << cnt << "," << mname << "," << mean << "," << stdDev << "\n";
    delete m;
  }
}

void test15Dists() {
  
  ROMol *m = SmilesToMol("c1ccccc1C");
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DGeomHelpers::initBoundsMat(mat);
  DistGeom::BoundsMatPtr mmat(mat);
  DGeomHelpers::setTopolBounds(*m, mmat);
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(2,6), 4.32, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(2,6), 4.16, 0.01), "");
  delete m;
  
  m = SmilesToMol("CC1=C(C(C)=CC=C2)C2=CC=C1");
  nat = m->getNumAtoms();
  mmat.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(mmat);
  DGeomHelpers::setTopolBounds(*m, mmat);
  
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0,4), 2.31, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0,4), 2.47, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(4,11), 4.11, 0.01), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(4,11), 4.27, 0.01) , "");
  
  delete m;

  m = SmilesToMol("C/C=C/C=C/C", 0, 1);
  nat = m->getNumAtoms();
 
  mmat.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(mmat);
  DGeomHelpers::setTopolBounds(*m, mmat);
  
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(0,4), 4.1874), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(0,4), 4.924), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getLowerBound(1,5), 4.1874), "");
  CHECK_INVARIANT(RDKit::feq(mmat->getUpperBound(1,5), 4.924) , "");
 
  delete m;
  m = SmilesToMol("NCc(c1)cccc1");
  delete m;
}

void testMultipleConfs() {
  std::string smi = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC";
  ROMol *m = SmilesToMol(smi, 0, 1);
  INT_VECT cids = DGeomHelpers::EmbedMultipleConfs(*m, 10, 30, 100, true,
						   false,-1);
  INT_VECT_CI ci;
  SDWriter writer("junk.sdf");
  double energy;
 
  for (ci = cids.begin(); ci != cids.end(); ci++) {
    writer.write(*m, *ci);
    ForceFields::ForceField *ff=UFF::constructForceField(*m, 10, *ci);
    ff->initialize();
    energy = ff->calcEnergy();
    //BOOST_LOG(rdInfoLog) << energy << std::endl;
    TEST_ASSERT(energy>100.0);
    TEST_ASSERT(energy<300.0);
    delete ff;
  }

}

void testOrdering() {
  std::string smi = "CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DGeomHelpers::initBoundsMat(mat);
  DistGeom::BoundsMatPtr mmat(mat);
  DGeomHelpers::setTopolBounds(*m, mmat);
  delete m;

  std::string smi2 = "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2";
  ROMol *m2 = SmilesToMol(smi2, 0, 1);
  DistGeom::BoundsMatrix *mat2 = new DistGeom::BoundsMatrix(nat);
  DGeomHelpers::initBoundsMat(mat2);
  DistGeom::BoundsMatPtr mmat2(mat2);
  DGeomHelpers::setTopolBounds(*m2, mmat2);
  delete m2;
}

#if 1
void testIssue227() {
  std::string smi = "CCOP1(OCC)=CC(c2ccccc2)=C(c2ccc([N+]([O-])=O)cc2)N=C1c1ccccc1";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  bool ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;

  smi = "OC(=O)c1cc2cc(c1)-c1c(O)c(ccc1)-c1cc(C(O)=O)cc(c1)OCCOCCO2";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *nmat = new DistGeom::BoundsMatrix(nat);
  bm.reset(nmat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;

}
#endif

void testIssue236() {
  std::string smi = "Cn1c2n(-c3ccccc3)c(=O)c3c(nc4ccc([N+]([O-])=O)cc4c3)c2c(=O)n(C)c1=O";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  bool ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  
  smi = "Cc1cccc2c1c(C3=CCC3)c(C)cc2";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *nmat = new DistGeom::BoundsMatrix(nat);
  bm.reset(nmat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  delete m;
}

void testIssue244() {
  // the bug here was a crash during the setting of the topological
  // bounds, so just completing the calls to setTopolBounds() indicates
  // success
  std::string smi = "NC1C(=O)N2C1SCC(Cl)=C2C(O)=O";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  delete m;

  smi = "NC1C(=O)N2C1SCC(Cl)=C2C(O)=O";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  delete m;
  
  smi = "CC1(C)SC2N(C1C(O)=O)C(=O)C2N";
  m = SmilesToMol(smi, 0, 1);
  nat = m->getNumAtoms();
  bm.reset(new DistGeom::BoundsMatrix(nat));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  delete m;
  
}

void testIssue251() {
  std::string smi = "COC=O";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
  TEST_ASSERT(RDKit::feq(bm->getLowerBound(0,3), 2.67, 0.01));
  TEST_ASSERT(RDKit::feq(bm->getUpperBound(0,3), 2.79, 0.01));
  delete m;
}

void testIssue276() {
  bool ok;
  std::string smi = "CP1(C)=CC=CN=C1C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);
#if 0
  for(unsigned int i=0;i<nat;i++){
    for(unsigned int j=0;j<nat;j++){
      if(i<j) std::cout << bm->getUpperBound(i,j) << " ";
      else if(i>j) std::cout << bm->getLowerBound(i,j) << " ";
      else std::cout << "0.00000" << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif  
  
  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);
  
  delete m;
}

void testIssue284() {
  bool ok;
  std::string smi = "CNC(=O)C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  // amide bonds are cis-oid:
  TEST_ASSERT(bm->getLowerBound(0,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(0,3)<3.0);

  delete m;

  smi = "CN(C)C(=O)C";
  m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  nat = m->getNumAtoms();

  DistGeom::BoundsMatrix *mat2 = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm2(mat2);
  DGeomHelpers::initBoundsMat(bm2);
  DGeomHelpers::setTopolBounds(*m, bm2);

  ok = DistGeom::triangleSmoothBounds(bm2);
  TEST_ASSERT(ok);

  // we've got no information to tell us cis-oid vs trans-oid here, so
  // the windows are huge:
  TEST_ASSERT(bm2->getLowerBound(0,4)<3.0);
  TEST_ASSERT(bm2->getUpperBound(0,4)>3.5);
  TEST_ASSERT(bm2->getLowerBound(2,4)<3.0);
  TEST_ASSERT(bm2->getUpperBound(2,4)>3.5);
  TEST_ASSERT(bm->getLowerBound(0,3)<bm2->getLowerBound(0,4));
  TEST_ASSERT(bm->getUpperBound(0,3)<bm2->getUpperBound(0,4));
  TEST_ASSERT(bm->getLowerBound(0,3)<bm2->getLowerBound(2,4));
  TEST_ASSERT(bm->getUpperBound(0,3)<bm2->getUpperBound(2,4));

  delete m;
}

void testIssue285() {
  bool ok;
  std::string smi = "CNC(=O)C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  unsigned int tgtNumber=10;
  INT_VECT cids = DGeomHelpers::EmbedMultipleConfs(*m, tgtNumber);
  TEST_ASSERT(cids.size()==tgtNumber);

  std::vector<std::string> molBlocks;
  for(INT_VECT_CI cid=cids.begin();cid!=cids.end();++cid){
    molBlocks.push_back(MolToMolBlock(*m,true,*cid));
  }
  for(std::vector<std::string>::const_iterator mbI=molBlocks.begin();
      mbI!=molBlocks.end();++mbI){
    for(std::vector<std::string>::const_iterator mbJ=mbI+1;
	mbJ!=molBlocks.end();++mbJ){
      TEST_ASSERT((*mbI)!=(*mbJ));
    }
    //std::cerr << (*mbI) << "\n$$$$\n";
  }
  delete m;
}


void testIssue355() {
  bool ok;
  std::string smi = "CNC(=O)C";
  ROMol *m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  unsigned int nat = m->getNumAtoms();
  DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(nat);
  DistGeom::BoundsMatPtr bm(mat);
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  TEST_ASSERT(bm->getLowerBound(0,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(0,3)<3.0);

  TEST_ASSERT(bm->getUpperBound(0,4)>3.2);
  TEST_ASSERT(bm->getLowerBound(0,4)>3.2);

  ok = DistGeom::triangleSmoothBounds(bm);
  TEST_ASSERT(ok);

  delete m;

  smi="CNC(=O)NC";
  m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  TEST_ASSERT(bm->getLowerBound(0,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(0,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(0,4)>3.2);
  TEST_ASSERT(bm->getLowerBound(0,4)>3.2);
  TEST_ASSERT(bm->getLowerBound(5,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(5,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(5,1)>3.2);
  TEST_ASSERT(bm->getLowerBound(5,1)>3.2);

  delete m;

  smi="CNC(=O)Nc1ccccc1";
  m = SmilesToMol(smi, 0, 1);
  TEST_ASSERT(m);
  bm.reset(new DistGeom::BoundsMatrix(m->getNumAtoms()));
  DGeomHelpers::initBoundsMat(bm);
  DGeomHelpers::setTopolBounds(*m, bm);

  TEST_ASSERT(bm->getLowerBound(0,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(0,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(0,4)>3.2);
  TEST_ASSERT(bm->getLowerBound(0,4)>3.2);
  TEST_ASSERT(bm->getLowerBound(5,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(5,3)<3.0);
  TEST_ASSERT(bm->getUpperBound(5,1)>3.2);
  TEST_ASSERT(bm->getLowerBound(5,1)>3.2);

  delete m;
}


void testRandomCoords() {
  std::string smiString = "CC1=C(C(C)=CC=C2)C2=CC=C1 c1ccccc1C C/C=C/CC \
                           C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4 C1CCCCS1(=O)(=O) c1ccccc1 \
                           C1CCCC1 C1CCCCC1 \
                           C1CC1(C)C C12(C)CC1CC2";
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/initCoords.random.sdf";
  SDMolSupplier sdsup(fname,true,false);
  //SDWriter writer("foo.sdf");
  //SDWriter writer(fname);

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();++token){
    std::string smi= *token;
    ROMol *m = SmilesToMol(smi, 0, 1);
    RWMol *m2 = (RWMol *)MolOps::addHs(*m);
    delete m;
    m=m2;
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, 1, true, true, 2,true,1,0,1e-2);
    CHECK_INVARIANT(cid >= 0, "");
    //writer.write(*m);
    //writer.flush();
#if 1
    m2 = static_cast<RWMol *>(sdsup.next());
    //ROMol *m2 = NULL;
    if(m2){
      TEST_ASSERT(m->getNumAtoms()==m2->getNumAtoms());
      unsigned int nat = m->getNumAtoms();
    
      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
#if 0
      BOOST_LOG(rdInfoLog) << "-----------------------" << std::endl;
      BOOST_LOG(rdInfoLog) << MolToMolBlock(*m2) << std::endl;
      BOOST_LOG(rdInfoLog) << "---" << std::endl;
      BOOST_LOG(rdInfoLog) << MolToMolBlock(*m) << std::endl;
      BOOST_LOG(rdInfoLog) << "-----------------------" << std::endl;
#endif
      for (unsigned int i = 0; i < nat; i++) {
	RDGeom::Point3D pt1i = conf1.getAtomPos(i);
	RDGeom::Point3D pt2i = conf2.getAtomPos(i);
	for(unsigned int j=i+1;j<nat;j++){
	  RDGeom::Point3D pt1j = conf1.getAtomPos(j);
	  RDGeom::Point3D pt2j = conf2.getAtomPos(j);
	  double d1=(pt1j-pt1i).length();
	  double d2=(pt2j-pt2i).length();
	  if(m->getBondBetweenAtoms(i,j)){
	    TEST_ASSERT(fabs(d1-d2)/d1<0.05);
	  }else{
	    TEST_ASSERT(fabs(d1-d2)/d1<0.1);
	  }
	}
      }
    }
    delete m2;
#endif
    delete m;
  }
}

void testIssue1989539() {
  {
    std::string smi="c1ccccc1.Cl";
    ROMol *m = SmilesToMol(smi, 0, 1);
    RWMol *m2 = (RWMol *)MolOps::addHs(*m);
    delete m;
    m=m2;
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
  {
    std::string smi="[Cl-].c1ccccc1C[NH3+]";
    RWMol *m = SmilesToMol(smi, 0, 1);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
}


void testConstrainedEmbedding() {
  std::string rdbase = getenv("RDBASE");
  std::string fname = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf";
  SDMolSupplier sdsup(fname);

  ROMol *ref=sdsup.next();
  {
    ROMol *test = new ROMol(*ref);
    std::map<int,RDGeom::Point3D> coords;
    coords[0]=ref->getConformer().getAtomPos(0);
    coords[1]=ref->getConformer().getAtomPos(1);
    coords[2]=ref->getConformer().getAtomPos(2);
    coords[3]=ref->getConformer().getAtomPos(3);
    coords[4]=ref->getConformer().getAtomPos(4);

#if 1
    int cid = DGeomHelpers::EmbedMolecule(*test,30,22,true,false,2.,true,1,&coords);
    TEST_ASSERT(cid>-1);
    
    MatchVectType alignMap;
    alignMap.push_back(std::make_pair(0,0));
    alignMap.push_back(std::make_pair(1,1));
    alignMap.push_back(std::make_pair(2,2));
    alignMap.push_back(std::make_pair(3,3));
    alignMap.push_back(std::make_pair(4,4));
    double ssd=MolAlign::alignMol(*test,*ref,-1,-1,&alignMap);
    BOOST_LOG(rdInfoLog)<<"ssd: "<<ssd<<std::endl;
    TEST_ASSERT(ssd<0.1);
#endif
    delete test;
  }

  {
    ROMol *test = sdsup.next();

    std::map<int,RDGeom::Point3D> coords;
    coords[4]=ref->getConformer().getAtomPos(0);
    coords[5]=ref->getConformer().getAtomPos(1);
    coords[6]=ref->getConformer().getAtomPos(2);
    coords[7]=ref->getConformer().getAtomPos(3);
    coords[8]=ref->getConformer().getAtomPos(4);
    int cid = DGeomHelpers::EmbedMolecule(*test,30,22,true,false,2.,true,1,&coords);
    TEST_ASSERT(cid>-1);
    
    MatchVectType alignMap;
    alignMap.push_back(std::make_pair(4,0));
    alignMap.push_back(std::make_pair(5,1));
    alignMap.push_back(std::make_pair(6,2));
    alignMap.push_back(std::make_pair(7,3));
    alignMap.push_back(std::make_pair(8,4));
    double ssd=MolAlign::alignMol(*test,*ref,-1,-1,&alignMap);
    BOOST_LOG(rdInfoLog)<<"ssd: "<<ssd<<std::endl;
    TEST_ASSERT(ssd<0.1);
    delete test;
  }
}

void testIssue2091864() {
  {
    std::string smi="C1C2CC12";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    delete m;
  }
  {
    std::string smi="C1CC2C3C1C23";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
  {
    std::string smi="c1ccc2c(c1)C1C3C2C13";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
}

void testIssue2091974() {
  {
    std::string smi="CCOC(OCC)(OCC)OCC";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 =MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    delete m2;
  }
  {
    std::string smi="O=N(=O)OCC(CON(=O)=O)(CON(=O)=O)CON(=O)=O";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 =MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    delete m2;
  }
}


void testIssue2835784() {
#if 1
  {
    std::string smi="C1C=C1";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
  {
    std::string smi="C1C=C1";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 =MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m2,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m2;
  }
  {
    std::string smi="C12=CCC1C2";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
#endif
  {
    std::string smi="C12=CCC1C2";
    ROMol *m = SmilesToMol(smi);
    ROMol *m2 =MolOps::addHs(*m);
    delete m;
    int cid = DGeomHelpers::EmbedMolecule(*m2);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m2,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m2;
  }
}

void testIssue3019283() {
  {
    std::string smi="C1=C2C1C1CC21";
    RWMol *m = SmilesToMol(smi);
    int cid = DGeomHelpers::EmbedMolecule(*m);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
}

void testIssue3238580() {
  {
    std::string smi="C1CCC2=CC12";
    RWMol *m = SmilesToMol(smi);

    DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());
    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);
    
    delete m;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.1.mol";
    RWMol *m = MolFileToMol(molfile);

    DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());

    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);
    
    delete m;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.2.mol";
    RWMol *m = MolFileToMol(molfile);

    DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());
    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);
    
    delete m;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.3.mol";
    RWMol *m = MolFileToMol(molfile);

    DistGeom::BoundsMatrix *mat = new DistGeom::BoundsMatrix(m->getNumAtoms());
    DistGeom::BoundsMatPtr bm(mat);
    DGeomHelpers::initBoundsMat(bm);
    DGeomHelpers::setTopolBounds(*m, bm);
    // if we get here it indicates everything was fine.
    // do bounds smoothing just to be sure:
    bool ok = DistGeom::triangleSmoothBounds(bm);
    TEST_ASSERT(ok);
    
    delete m;
  }
}

void testIssue3483968() {
  {
    std::string rdbase = getenv("RDBASE");
    std::string molfile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/Issue3483968.mol";
    RWMol *m = MolFileToMol(molfile);
    TEST_ASSERT(m);
    
    int cid = DGeomHelpers::EmbedMolecule(*m,0,-1,true,false,2.0,true,1,0,1e-3,true);
    TEST_ASSERT(cid >= 0);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(*m,10,30,
                                                           1,true,false,2.0,true,1,-1.0,0,
                                                           1e-3,true);
    TEST_ASSERT(cids.size() == 10);
    TEST_ASSERT(std::find(cids.begin(),cids.end(),-1)==cids.end());
    delete m;
  }
}

#ifdef RDK_TEST_MULTITHREADED
namespace {
  void runblock(const std::vector<ROMol *> &mols,const std::vector<double> &energies,
                unsigned int count,unsigned int idx){
    for(unsigned int j=0;j<100;j++){
      for(unsigned int i=0;i<mols.size();++i){
        if(i%count != idx) continue;
        ROMol mol(*mols[i]);
        std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(mol,10,30,0xFEED);
        TEST_ASSERT(cids.size() == 10);
        ForceFields::ForceField *field = 0;
        try {
          field = UFF::constructForceField(mol,100.0,cids[0]);
        } catch (...) {
          field = 0;
        }
        TEST_ASSERT(field);
        field->initialize();
        double eng=field->calcEnergy();
        if(!feq(eng,energies[i])){
          std::cerr<<i<<" iter "<<j<<" "<<energies[i]<<" != "<<eng<<std::endl;
        }
        
        TEST_ASSERT(feq(eng,energies[i]));
        delete field;
      }
    }
  };
}


#include <boost/thread.hpp>  
void testMultiThread(){
  std::cerr<<"building molecules"<<std::endl;
  //std::string smi="C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4";
  std::string smi="c1ccc2c(c1)C1C3C2C13";
  std::vector<ROMol *> mols;
  for(unsigned int i=0;i<100;++i){
    RWMol *m = SmilesToMol(smi);
    mols.push_back(m);
  }


  std::cerr<<"generating reference data"<<std::endl;
  std::vector<double> energies(mols.size(),0.0);
  for(unsigned int i=0;i<mols.size();++i){
    ROMol mol(*mols[i]);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(mol,10,30,0xFEED);
    TEST_ASSERT(cids.size() == 10);
    ForceFields::ForceField *field = 0;
    try {
      field = UFF::constructForceField(mol,100.0,cids[0]);
    } catch (...) {
      field = 0;
    }
    TEST_ASSERT(field);
    field->initialize();
    double eng=field->calcEnergy();
    TEST_ASSERT(eng!=0.0);
    energies[i]=eng;
    delete field;
  }

  std::cerr<<"validating reference data"<<std::endl;
  for(unsigned int i=0;i<mols.size();++i){
    ROMol mol(*mols[i]);
    std::vector<int> cids=DGeomHelpers::EmbedMultipleConfs(mol,10,30,0xFEED);
    TEST_ASSERT(cids.size() == 10);
    ForceFields::ForceField *field = 0;
    try {
      field = UFF::constructForceField(mol,100.0,cids[0]);
    } catch (...) {
      field = 0;
    }
    TEST_ASSERT(field);
    field->initialize();
    double eng=field->calcEnergy();
    TEST_ASSERT(feq(eng,energies[i]));
    delete field;
  }

  boost::thread_group tg;

  std::cerr<<"processing"<<std::endl;
  unsigned int count=4;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock,mols,energies,count,i));
  }
  tg.join_all();

  for(unsigned int i=0;i<mols.size();++i) delete mols[i];

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread(){
}
#endif

void testGithub55() {
  {
    std::string smiles = "c1cnco1";
    RWMol *core = SmilesToMol(smiles);
    TEST_ASSERT(core);
    
    int cid = DGeomHelpers::EmbedMolecule(*core);
    TEST_ASSERT(cid >= 0);

    smiles = "o1cncc1C";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    std::map<int,RDGeom::Point3D> coords;
    coords[0]=core->getConformer().getAtomPos(4);
    coords[1]=core->getConformer().getAtomPos(3);
    coords[2]=core->getConformer().getAtomPos(2);
    coords[3]=core->getConformer().getAtomPos(1);
    coords[4]=core->getConformer().getAtomPos(0);
    cid = DGeomHelpers::EmbedMolecule(*mol,50,22,true,false,2.,true,1,&coords);
    TEST_ASSERT(cid>-1);

    delete core;
    delete mol;
  }
  {
    std::string smiles = "c1cncs1";
    RWMol *core = SmilesToMol(smiles);
    TEST_ASSERT(core);
    
    int cid = DGeomHelpers::EmbedMolecule(*core);
    TEST_ASSERT(cid >= 0);

    smiles = "s1cncc1C";
    RWMol *mol = SmilesToMol(smiles);
    TEST_ASSERT(mol);

    std::map<int,RDGeom::Point3D> coords;
    coords[0]=core->getConformer().getAtomPos(4);
    coords[1]=core->getConformer().getAtomPos(3);
    coords[2]=core->getConformer().getAtomPos(2);
    coords[3]=core->getConformer().getAtomPos(1);
    coords[4]=core->getConformer().getAtomPos(0);
    cid = DGeomHelpers::EmbedMolecule(*mol,50,22,true,false,2.,true,1,&coords);
    TEST_ASSERT(cid>-1);

    delete core;
    delete mol;
  }
}

void testGithub256() {
  {
    RWMol *mol = new RWMol();
    TEST_ASSERT(mol);

    bool ok=false;
    try{
      DGeomHelpers::EmbedMolecule(*mol);
      ok=false;
    } catch (const ValueErrorException &e) {
      ok=true;
    }
    TEST_ASSERT(ok);
    delete mol;
  }
}

int main() { 
  RDLog::InitLogs();
    
  BOOST_LOG(rdInfoLog) << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing DistGeomHelpers\n";

#if 1
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test2 \n\n";
  test2(); 
 
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test3 \n\n";
  test3();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test4 \n\n";
  test4();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test5 \n\n";
  test5();
  
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test6 \n\n";
  test6();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test15Dists \n\n";
  test15Dists();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test Issue 215 \n\n";
  testIssue215();
   
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testMultipleConfs \n\n";
  testMultipleConfs();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue227 \n\n";
  testIssue227();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue236 \n\n";
  testIssue236();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testOrdering \n\n";
  testOrdering();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue244 \n\n";
  testIssue244();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue251 \n\n";
  testIssue251();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue276 \n";
  testIssue276();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue284 \n\n";
  testIssue284();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue285 \n\n";
  testIssue285();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue355 \n\n";
  testIssue355();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testRandomCoords \n\n";
  testRandomCoords();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 1989539 \n\n";
  testIssue1989539();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 2091864 \n\n";
  testIssue2091864();


  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test constrained embedding \n\n";
  testConstrainedEmbedding();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 2091974 \n\n";
  testIssue2091974();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 2835784 \n\n";
  testIssue2835784();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 3019283 \n\n";
  testIssue3019283();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 3238580 \n\n";
  testIssue3238580();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test sf.net issue 3483968 \n\n";
  testIssue3483968();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test github issue 55 \n\n";
  testGithub55();

#endif
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test1 \n\n";
  test1();
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test multi-threading \n\n";
  testMultiThread();
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test github issue 256: handling of zero-atom molecules\n\n";
  testGithub256();





  BOOST_LOG(rdInfoLog) << "*******************************************************\n";

  return(0);
}

