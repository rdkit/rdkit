// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
#include <math.h>
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
  //SDWriter writer(fname);

  boost::char_separator<char> spaceSep(" ");
  tokenizer tokens(smiString,spaceSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();++token){
    std::string smi= *token;
    RWMol *m = SmilesToMol(smi, 0, 1); 
    int cid = DGeomHelpers::EmbedMolecule(*m, 10, true);
    CHECK_INVARIANT(cid >= 0, "");

    ROMol *m2 = sdsup.next();
    //writer.write(m);
    //ROMol *m2 = NULL;
    if(m2){
      unsigned int nat = m->getNumAtoms();
    
      const Conformer &conf1 = m->getConformer(0);
      const Conformer &conf2 = m2->getConformer(0);
      //BOOST_LOG(rdDebugLog) << MolToMolBlock(m) << std::endl;
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
    delete m;
    delete m2;
  }
}

void computeDistMat(const DistGeom::PointPtrVect &origCoords, RDNumeric::DoubleSymmMatrix &distMat) {
  unsigned int N = origCoords.size();
  CHECK_INVARIANT(N == distMat.numRows(), "");
  unsigned int i, j;
  RDGeom::Point3D pti, ptj;
  double d;
  for (i = 1; i < N; i++) {
    pti = *origCoords[i];
    for (j = 0; j < i; j++) {
      ptj = *origCoords[j];
      ptj -= pti;
      d = ptj.length();
      distMat.setVal(i,j, d);
    }
  }
}

void computeMolDmat(ROMol &mol, RDNumeric::DoubleSymmMatrix &distMat) {
  DistGeom::PointPtrVect origCoords;
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
  TEST_ASSERT(fabs(bm->getLowerBound(0,3) - dmat->getVal(0,3)) < 0.1);
  
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
    DistGeom::PointPtrVect origCoords, newCoords;
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
  DGeomHelpers::EmbedMolecule(*m, 10, true);//etCoords(*m, iter);
  std::string fname = "test.mol";
  MolToMolFile(*m, fname);
  delete m;
}

void test5() {
  // some real CDK2 molecules lets see how many fail
  std::string rdbase = getenv("RDBASE");
  std::string smifile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv";
  SmilesMolSupplier smiSup(smifile, ",", 0, 1);
  
  std::string sdfile = rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/embedDistOpti2.sdf";
  SDWriter writer(sdfile);
  ROMol *mol;
  int i = 0;
  int cid;
  while (1) {
    try {
      i++;
      mol = smiSup.next();
      std::string mname, mname2;
      cid = DGeomHelpers::EmbedMolecule(*mol, 10, true); //getCoords(*mol, iter);
      TEST_ASSERT(cid>-1);
      mol->getProp("_Name", mname);
      writer.write(*mol);
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
      int cid = DGeomHelpers::EmbedMolecule(*m, 10, true, false);
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
                                      true, 1, 1e-3, 5.0);
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
  TEST_ASSERT(RDKit::feq(bm->getLowerBound(0,3), 2.75, 0.01));
  TEST_ASSERT(RDKit::feq(bm->getUpperBound(0,3), 2.87, 0.01));
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



int main() { 
  RDLog::InitLogs();
    
  BOOST_LOG(rdInfoLog) << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing DistGeomHelpers\n";

#if 1
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test1 \n\n";
  test1();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test2 \n\n";
  test2();
  
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test3 \n\n";
  test3();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t test4 \n\n";
  test4();

  //BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  //BOOST_LOG(rdInfoLog) << "\t test5 \n\n";
  //test5();
  
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
#endif
  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t testIssue355 \n\n";
  testIssue355();

  BOOST_LOG(rdInfoLog) << "*******************************************************\n";
  return(0);
}

