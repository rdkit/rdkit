//  $Id$
//
//  Copyright (C) 2003-2004 Rational Discovery LLC
//   All Rights Reserved
//
#include "TemplEnum.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;
using namespace TemplateEnum;

#include <math.h>

bool feq(RDGeom::Point3D p1,RDGeom::Point3D p2,double tol=1e-4){
  return feq(p1.x,p2.x,tol)&&feq(p1.y,p2.y,tol)&&feq(p1.z,p2.z,tol);
}

void test1(){
  // single attachment point, small list
  std::cout << " ----------> Test1 " << std::endl;
  RWMOL_SPTR_VECT library;
  std::vector<const char *>fileNames;
  fileNames.push_back("esters.2.sdf");
  fileNames.push_back("esters.2.sdf");
  library = enumFromFiles("template.1.mol",fileNames);

  CHECK_INVARIANT(library.size()==2,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==10,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==11,"");

  library.clear();
  std::cout << " <---------- Done " << std::endl;
}


void test2(){
  // single attachment point, larger list
  std::cout << " ----------> Test2 " << std::endl;
  RWMOL_SPTR_VECT library;
  std::vector<const char *>fileNames;
  fileNames.push_back("esters.sdf");
  fileNames.push_back("esters.sdf");
  library = enumFromFiles("template.1.mol",fileNames);

  CHECK_INVARIANT(library.size()==6,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==10,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==11,"");
  CHECK_INVARIANT(library[2]->getNumAtoms()==12,"");
  CHECK_INVARIANT(library[3]->getNumAtoms()==12,"");
  CHECK_INVARIANT(library[4]->getNumAtoms()==13,"");
  CHECK_INVARIANT(library[5]->getNumAtoms()==13,"");

  library.clear();
  std::cout << " <---------- Done " << std::endl;
}

void test3(){
  // two attachment points, small list
  std::cout << " ----------> Test3 " << std::endl;
  RWMOL_SPTR_VECT library;
  std::vector<const char *>fileNames;
  fileNames.push_back("esters.2.sdf");
  fileNames.push_back("esters.2.sdf");
  library = enumFromFiles("template.mol",fileNames);

  CHECK_INVARIANT(library.size()==4,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==14,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==15,"");
  CHECK_INVARIANT(library[2]->getNumAtoms()==15,"");
  CHECK_INVARIANT(library[3]->getNumAtoms()==16,"");

  library.clear();
  std::cout << " <---------- Done " << std::endl;
}

void test4(){
  // test templates that have repeated attachment points
  std::cout << " ----------> Test4 " << std::endl;
  RWMOL_SPTR_VECT library;
  std::vector<const char *>fileNames;
  fileNames.push_back("esters.2.sdf");
  fileNames.push_back("esters.2.sdf");
  library = enumFromFiles("template.2.mol",fileNames);

  CHECK_INVARIANT(library.size()==2,"");
  library[0]->debugMol(std::cout);
  std::cout << "smi0: " << MolToSmiles(*library[0],0) << std::endl;
  std::cout << "smi1: " << MolToSmiles(*library[1],0) << std::endl;
  CHECK_INVARIANT(library[0]->getNumAtoms()==14,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==16,"");


  library.clear();
  std::cout << " <---------- Done " << std::endl;
}

void test5(){
  // test working from SMILES
  std::cout << " ----------> Test5 " << std::endl;
  RWMol *m1 = SmilesToMol("[Xa]CC([Xb])CC",0,0);
  CHECK_INVARIANT(m1,"");
  markAttachmentPoints(m1,'X');
  
  RWMOL_SPTR_VECT sidechains;
  sidechains.push_back(RWMOL_SPTR(SmilesToMol("[X]OC(=O)",0,0)));
  sidechains.push_back(RWMOL_SPTR(SmilesToMol("[X]OC(=O)C",0,0)));
  sidechains.push_back(RWMOL_SPTR(SmilesToMol("[X]OC(=O)CCC",0,0)));
  prepareSidechains(&sidechains,'X');

  VECT_RWMOL_SPTR_VECT allSideChains;
  allSideChains.push_back(sidechains);
  allSideChains.push_back(sidechains);

  RWMOL_SPTR_VECT library;
  library = enumerateLibrary(m1,allSideChains,false);

  CHECK_INVARIANT(library.size()==9,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==10,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==11,"");
  CHECK_INVARIANT(library[2]->getNumAtoms()==13,"");
  CHECK_INVARIANT(library[3]->getNumAtoms()==11,"");
  CHECK_INVARIANT(library[4]->getNumAtoms()==12,"");
  CHECK_INVARIANT(library[5]->getNumAtoms()==14,"");
  CHECK_INVARIANT(library[6]->getNumAtoms()==13,"");
  CHECK_INVARIANT(library[7]->getNumAtoms()==14,"");
  CHECK_INVARIANT(library[8]->getNumAtoms()==16,"");

  library.clear();
  std::cout << " <---------- Done " << std::endl;
}

void test6(){
  // test working from SMILES with no matches
  std::cout << " ----------> Test6 " << std::endl;
  RWMol *m1 = SmilesToMol("[Xa]CC([Xb])CC",0,0);
  CHECK_INVARIANT(m1,"");
  markAttachmentPoints(m1,'X');
  
  RWMOL_SPTR_VECT sidechains;
  sidechains.push_back(RWMOL_SPTR(SmilesToMol("OC(=O)",0,0)));
  sidechains.push_back(RWMOL_SPTR(SmilesToMol("OC(=O)C",0,0)));
  sidechains.push_back(RWMOL_SPTR(SmilesToMol("OC(=O)CCC",0,0)));
  prepareSidechains(&sidechains,'X');

  VECT_RWMOL_SPTR_VECT allSideChains;
  allSideChains.push_back(sidechains);
  allSideChains.push_back(sidechains);

  RWMOL_SPTR_VECT library;
  library = enumerateLibrary(m1,allSideChains,false);

  CHECK_INVARIANT(library.size()==1,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==6,"");

  library.clear();
  std::cout << " <---------- Done " << std::endl;
}

void test7(){
  // test transforms
  std::cout << " ----------> Test7 " << std::endl;
  RWMOL_SPTR_VECT library;
  RWMOL_SPTR mol;
  Atom *at1,*at2;
  int i;
  std::vector<const char *>fileNames;
  fileNames.push_back("Ts.1.sdf");
  library = enumFromFiles("box.1.mol",fileNames);
  CHECK_INVARIANT(library.size()==1,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==8,"");
  mol = library[0];
  at1=mol->getAtomWithIdx(0);
  at2=mol->getAtomWithIdx(4);
  CHECK_INVARIANT(feq(mol->getConformer().getAtomPos(at1->getIdx()).x,mol->getConformer().getAtomPos(at2->getIdx()).x),"");
  CHECK_INVARIANT(mol->getConformer().getAtomPos(at1->getIdx()).y-mol->getConformer().getAtomPos(at2->getIdx()).y==-1.0,"");
  library.clear();

  // try another orientation of the sidechain molecule:
  fileNames.clear();
  fileNames.push_back("Ts.4.sdf");
  library = enumFromFiles("box.1.mol",fileNames);
  CHECK_INVARIANT(library.size()==1,"");
  CHECK_INVARIANT(library[0]->getNumAtoms()==8,"");
  mol = library[0];
  std::cout << MolToMolBlock(*mol);
  std::cout << std::endl;

  at1=mol->getAtomWithIdx(0);
  at2=mol->getAtomWithIdx(7);
  TEST_ASSERT(feq(mol->getConformer().getAtomPos(7),RDGeom::Point3D(-.5,2.5,0.0)));
  library.clear();

  
  // now use an SD file that has the same molecule in different
  // orientations as sidechains:
  fileNames.clear();
  fileNames.push_back("Ts.sdf");
  library = enumFromFiles("box.1.mol",fileNames);
  CHECK_INVARIANT(library.size()==4,"");

  CHECK_INVARIANT(library[0]->getNumAtoms()==8,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==8,"");
  CHECK_INVARIANT(library[2]->getNumAtoms()==8,"");
  CHECK_INVARIANT(library[3]->getNumAtoms()==8,"");

  for(i=0;i<library.size();i++){
    std::cout << "------ Mol: " << i << "------" << std::endl;
    mol = library[i];
    std::cout << MolToMolBlock(*mol);
    std::cout << std::endl;
  }
  for(i=0;i<library[0]->getNumAtoms();i++){
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[1]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[1]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[1]->getConformer().getAtomPos(i).z),"");

    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[2]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[2]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[2]->getConformer().getAtomPos(i).z),"");

    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i),library[3]->getConformer().getAtomPos(i)),"");
  }
  library.clear();

  // move the attachment point on the template.  This should
  // make no difference.
  library = enumFromFiles("box.1a.mol",fileNames);
  CHECK_INVARIANT(library.size()==4,"");

  CHECK_INVARIANT(library[0]->getNumAtoms()==8,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==8,"");
  CHECK_INVARIANT(library[2]->getNumAtoms()==8,"");
  CHECK_INVARIANT(library[3]->getNumAtoms()==8,"");

  for(i=0;i<library[0]->getNumAtoms();i++){
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[1]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[1]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[1]->getConformer().getAtomPos(i).z),"");

    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[2]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[2]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[2]->getConformer().getAtomPos(i).z),"");

    //std::cout << i << "\t" << library[0]->getConformer().getAtomPos(i) << std::endl;
    //std::cout << "\t" << library[3]->getConformer().getAtomPos(i) << std::endl;
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[3]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[3]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[3]->getConformer().getAtomPos(i).z),"");
  }
  library.clear();

  // move the attachment point on the template.  This should
  // make no difference.
  library = enumFromFiles("box.mol",fileNames);
  CHECK_INVARIANT(library.size()==4,"");

  CHECK_INVARIANT(library[0]->getNumAtoms()==20,"");
  CHECK_INVARIANT(library[1]->getNumAtoms()==20,"");
  CHECK_INVARIANT(library[2]->getNumAtoms()==20,"");
  CHECK_INVARIANT(library[3]->getNumAtoms()==20,"");

  for(i=0;i<library[0]->getNumAtoms();i++){
    at1 = library[0]->getAtomWithIdx(i);

    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[1]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[1]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[1]->getConformer().getAtomPos(i).z),"");

    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[2]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[2]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[2]->getConformer().getAtomPos(i).z),"");

    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).x,library[3]->getConformer().getAtomPos(i).x),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).y,library[3]->getConformer().getAtomPos(i).y),"");
    CHECK_INVARIANT(feq(library[0]->getConformer().getAtomPos(i).z,library[3]->getConformer().getAtomPos(i).z),"");
  }
  library.clear();


  
  std::cout << " <---------- Done " << std::endl;
}

#if 0
void testCoords(){
  std::cout << " ----------> Test Coords " << std::endl;
  RWMol *m1 = SmilesToMol("[Xa]C",0,0);
  CHECK_INVARIANT(m1,"");
  m1->getConformer().setAtomPos(1,Point3D(0,0,0));
  m1->getAtomWithIdx(0)->setPos(0.0,1.5,0.0);

  RWMol *m2 = SmilesToMol("[X]C(C)(C)C",0,0);
  TEST_ASSERT(m2);
  m2->getAtomWithIdx(0)->setPos(3.0,1,1);
  m2->getAtomWithIdx(1)->setPos(1.5,1,1);
  m2->getAtomWithIdx(2)->setPos(0.0,1,1);
  m2->getAtomWithIdx(3)->setPos(1.5,2.5,1);
  m2->getAtomWithIdx(4)->setPos(1.5,1,2.5);

  orientSidechain(m1,m2,0,0);
  std::cout << "Atoms:" << std::endl;
  for(int i=0;i<m2->getNumAtoms();i++){
    std::cout << i << " " << m2->getAtomWithIdx(i)->getPos() << std::endl;
  }

  TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getPos(),
		  m2->getAtomWithIdx(1)->getPos()));

  TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getPos(),
		  RDGeom::Point3D(0,0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(1)->getPos(),
		  RDGeom::Point3D(0,1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(2)->getPos(),
		  RDGeom::Point3D(0,3.0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(3)->getPos(),
		  RDGeom::Point3D(1.5,1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(4)->getPos(),
		  RDGeom::Point3D(0,1.5,1.5)));
  
  // ---------------------------------------------------
  std::cout << "-*-*-*-*-*-*-*-*-*-*-*-*-" << std::endl;
  m1->getAtomWithIdx(1)->setPos(0,0,0);
  m1->getAtomWithIdx(0)->setPos(0.0,-1.5,0.0);
  m2->getAtomWithIdx(0)->setPos(3.0,1,1);
  m2->getAtomWithIdx(1)->setPos(1.5,1,1);
  m2->getAtomWithIdx(2)->setPos(0.0,1,1);
  m2->getAtomWithIdx(3)->setPos(1.5,2.5,1);
  m2->getAtomWithIdx(4)->setPos(1.5,1,2.5);

  orientSidechain(m1,m2,0,0);
  std::cout << "Atoms:" << std::endl;
  for(int i=0;i<m2->getNumAtoms();i++){
    std::cout << i << " " << m2->getAtomWithIdx(i)->getPos() << std::endl;
  }
  TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getPos(),
		  m2->getAtomWithIdx(1)->getPos()));

  TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getPos(),
		  RDGeom::Point3D(0,0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(1)->getPos(),
		  RDGeom::Point3D(0,-1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(2)->getPos(),
		  RDGeom::Point3D(0,-3.0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(3)->getPos(),
		  RDGeom::Point3D(-1.5,-1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(4)->getPos(),
		  RDGeom::Point3D(0,-1.5,1.5)));

  // ---------------------------------------------------
  std::cout << "-*-*-*-*-*-*-*-*-*-*-*-*-" << std::endl;
  m1->getAtomWithIdx(1)->setPos(1,0,0);
  m1->getAtomWithIdx(0)->setPos(1.0,-1.5,0);
  m2->getAtomWithIdx(0)->setPos(3.0,1,1);
  m2->getAtomWithIdx(1)->setPos(1.5,1,1);
  m2->getAtomWithIdx(2)->setPos(0.0,1,1);
  m2->getAtomWithIdx(3)->setPos(1.5,2.5,1);
  m2->getAtomWithIdx(4)->setPos(1.5,1,2.5);

  orientSidechain(m1,m2,0,0);
  std::cout << "Atoms:" << std::endl;
  for(int i=0;i<m2->getNumAtoms();i++){
    std::cout << i << " " << m2->getAtomWithIdx(i)->getPos() << std::endl;
  }
  TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getPos(),
		  m2->getAtomWithIdx(1)->getPos()));

  TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getPos(),
		  RDGeom::Point3D(1,0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(1)->getPos(),
		  RDGeom::Point3D(1,-1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(2)->getPos(),
		  RDGeom::Point3D(1,-3.0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(3)->getPos(),
		  RDGeom::Point3D(-0.5,-1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(4)->getPos(),
		  RDGeom::Point3D(1,-1.5,1.5)));

  // ---------------------------------------------------
  std::cout << "-*-*-*-*-*-*-*-*-*-*-*-*-" << std::endl;
  m1->getAtomWithIdx(1)->setPos(3.0,0,0);
  m1->getAtomWithIdx(0)->setPos(1.5,0.0,0.0);
  m2->getAtomWithIdx(0)->setPos(3.0,1,1);
  m2->getAtomWithIdx(1)->setPos(1.5,1,1);
  m2->getAtomWithIdx(2)->setPos(0.0,1,1);
  m2->getAtomWithIdx(3)->setPos(1.5,2.5,1);
  m2->getAtomWithIdx(4)->setPos(1.5,1,2.5);

  orientSidechain(m1,m2,0,0);
  std::cout << "Atoms:" << std::endl;
  for(int i=0;i<m2->getNumAtoms();i++){
    std::cout << i << " " << m2->getAtomWithIdx(i)->getPos() << std::endl;
  }
  TEST_ASSERT(feq(m1->getAtomWithIdx(0)->getPos(),
		  m2->getAtomWithIdx(1)->getPos()));

  TEST_ASSERT(feq(m2->getAtomWithIdx(0)->getPos(),
		  RDGeom::Point3D(3,0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(1)->getPos(),
		  RDGeom::Point3D(1.5,0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(2)->getPos(),
		  RDGeom::Point3D(0,0,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(3)->getPos(),
		  RDGeom::Point3D(1.5,1.5,0)));
  TEST_ASSERT(feq(m2->getAtomWithIdx(4)->getPos(),
		  RDGeom::Point3D(1.5,0,1.5)));


  delete m1;
  delete m2;
  
  std::cout << " <---------- Done " << std::endl;
}

#endif
int main(){
#if 1
  test1();
  test2();
  test3();
  //test4();
  test5();
  test6();
#endif
  test7();
  //testCoords();
}
