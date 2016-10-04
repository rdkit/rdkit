//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifdef _MSC_VER
// disable warnings about getenv in visual C++
#define _CRT_SECURE_NO_WARNINGS
#endif


#include <iostream>
#include <fstream>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/Descriptors/MolDescriptors3D.h>

using namespace RDKit;
using namespace RDKit::Descriptors;

bool compare(const std::string &inm,double ref,double val,double tol=1e-3){
  if(fabs(ref-val)>.001){
    std::cerr<<"value mismatch: "<<inm<<" "<<ref<<" "<<val<<std::endl;
  }
  return fabs(ref-val)<tol;
}

void testPMI1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic PMI tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";

  RDKit::SDMolSupplier reader(sdfName,true,false);
  std::string fName = pathName+"/Code/GraphMol/Descriptors/test_data/PMI_egfr.out";
  std::ifstream instrm(fName.c_str());
  int nDone=0;
  while(!reader.atEnd()){
    RDKit::ROMol *m=reader.next();
    TEST_ASSERT(m);
    RDKit::ROMol mcpy(*m);
    std::string nm;
    m->getProp("_Name",nm);
    std::string inm;
    instrm>>inm;
    TEST_ASSERT(inm==nm);
    double val;
    double pmi1_m,pmi2_m,pmi3_m,pmi1_nom,pmi2_nom,pmi3_nom;
    instrm>>pmi1_m;
    instrm>>pmi2_m;
    instrm>>pmi3_m;
    instrm>>pmi1_nom;
    instrm>>pmi2_nom;
    instrm>>pmi3_nom;

    val = RDKit::Descriptors::PMI1(*m);
    TEST_ASSERT(compare(inm,pmi1_m,val));
    val = RDKit::Descriptors::PMI2(*m);
    TEST_ASSERT(compare(inm,pmi2_m,val));
    val = RDKit::Descriptors::PMI3(*m);
    TEST_ASSERT(compare(inm,pmi3_m,val));

    val = RDKit::Descriptors::PMI1(*m,-1,false);
    TEST_ASSERT(compare(inm,pmi1_nom,val));
    val = RDKit::Descriptors::PMI2(*m,-1,false);
    TEST_ASSERT(compare(inm,pmi2_nom,val));
    val = RDKit::Descriptors::PMI3(*m,-1,false);
    TEST_ASSERT(compare(inm,pmi3_nom,val));

    // now try doing it in the reverse order to make sure caching doesn't
    // screw up.
    val = RDKit::Descriptors::PMI1(mcpy,-1,false);
    TEST_ASSERT(compare(inm,pmi1_nom,val));
    val = RDKit::Descriptors::PMI2(mcpy,-1,false);
    TEST_ASSERT(compare(inm,pmi2_nom,val));
    val = RDKit::Descriptors::PMI3(mcpy,-1,false);
    TEST_ASSERT(compare(inm,pmi3_nom,val));
    val = RDKit::Descriptors::PMI1(mcpy);
    TEST_ASSERT(compare(inm,pmi1_m,val));
    val = RDKit::Descriptors::PMI2(mcpy);
    TEST_ASSERT(compare(inm,pmi2_m,val));
    val = RDKit::Descriptors::PMI3(mcpy);
    TEST_ASSERT(compare(inm,pmi3_m,val));


    delete m;
    ++nDone;
}
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testPMIEdges(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    PMI edge cases." << std::endl;

  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/linear.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::PMI1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI2(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI3(*m);
    TEST_ASSERT(val>=1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/linear_2atom.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::PMI1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI2(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI3(*m);
    TEST_ASSERT(val>=1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/planar.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::PMI1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI2(*m);
    TEST_ASSERT(val>=1e-4);
    val = RDKit::Descriptors::PMI3(*m);
    TEST_ASSERT(val>=1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/planar_3atom.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::PMI1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI2(*m);
    TEST_ASSERT(val>=1e-4);
    val = RDKit::Descriptors::PMI3(*m);
    TEST_ASSERT(val>=1e-4);

    delete m;
  }

  {
    RDKit::RWMol m;
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addConformer(new RDKit::Conformer(m.getNumAtoms()));
    double val = RDKit::Descriptors::PMI1(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI2(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::PMI3(m);
    TEST_ASSERT(fabs(val)<1e-4);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}


void testNPR1(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic NPR tests." << std::endl;

  std::string pathName = getenv("RDBASE");
  std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf";
  RDKit::SDMolSupplier reader(sdfName,true,false);

  int nDone=0;
  while(!reader.atEnd()){
    RDKit::ROMol *m=reader.next();
    TEST_ASSERT(m);
    RDKit::ROMol mcpy(*m);
    std::string nm;
    m->getProp("_Name",nm);

    double val;
    double pmi1_m,pmi2_m,pmi3_m,pmi1_nom,pmi2_nom,pmi3_nom;
    pmi1_m = RDKit::Descriptors::PMI1(*m);
    pmi2_m = RDKit::Descriptors::PMI2(*m);
    pmi3_m = RDKit::Descriptors::PMI3(*m);
    pmi1_nom = RDKit::Descriptors::PMI1(*m,-1,false);
    pmi2_nom = RDKit::Descriptors::PMI2(*m,-1,false);
    pmi3_nom = RDKit::Descriptors::PMI3(*m,-1,false);

    val = RDKit::Descriptors::NPR1(*m);
    compare(nm,pmi1_m/pmi3_m,val);
    val = RDKit::Descriptors::NPR2(*m);
    compare(nm,pmi2_m/pmi3_m,val);

    delete m;
    ++nDone;
  }
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testNPREdges(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    NPR edge cases." << std::endl;

  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/linear.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::NPR1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::NPR2(*m);
    TEST_ASSERT(fabs(val)<1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/linear_2atom.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::NPR1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::NPR2(*m);
    TEST_ASSERT(fabs(val)<1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/planar.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::NPR1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::NPR2(*m);
    TEST_ASSERT(val>=1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/planar_3atom.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::NPR1(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::NPR2(*m);
    TEST_ASSERT(val>=1e-4);

    delete m;
  }
  {
    RDKit::RWMol m;
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addConformer(new RDKit::Conformer(m.getNumAtoms()));
    double val = RDKit::Descriptors::NPR1(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::NPR2(m);
    TEST_ASSERT(fabs(val)<1e-4);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test3DEdges(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    3D descriptor edge cases." << std::endl;

  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/linear.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::radiusOfGyration(*m);
    TEST_ASSERT(fabs(val)<1e-2);
    val = RDKit::Descriptors::inertialShapeFactor(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::eccentricity(*m);
    TEST_ASSERT(fabs(1.0-val)<1e-4);
    val = RDKit::Descriptors::asphericity(*m);
    TEST_ASSERT(fabs(1.0-val)<1e-4);
    val = RDKit::Descriptors::spherocityIndex(*m);
    TEST_ASSERT(fabs(val)<1e-4);

    delete m;
  }
  {
    std::string pathName = getenv("RDBASE");
    std::string sdfName = pathName+"/Code/GraphMol/Descriptors/test_data/planar.mol";

    RDKit::ROMol *m=MolFileToMol(sdfName);
    TEST_ASSERT(m);
    double val;

    val = RDKit::Descriptors::radiusOfGyration(*m);
    TEST_ASSERT(fabs(val)>1e-2);
    val = RDKit::Descriptors::inertialShapeFactor(*m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::eccentricity(*m);
    TEST_ASSERT(fabs(1.0-val)<1e-4);
    val = RDKit::Descriptors::asphericity(*m);
    TEST_ASSERT(fabs(0.5-val)<1e-4);
    val = RDKit::Descriptors::spherocityIndex(*m);
    TEST_ASSERT(fabs(val)<1e-4);

    delete m;
  }
  { // octahedron
    RDKit::RWMol m;
    m.addAtom(new RDKit::Atom(1));
    m.addAtom(new RDKit::Atom(1));
    m.addAtom(new RDKit::Atom(1));
    m.addAtom(new RDKit::Atom(1));
    m.addAtom(new RDKit::Atom(1));
    m.addAtom(new RDKit::Atom(1));
    m.addConformer(new RDKit::Conformer(m.getNumAtoms()));
    m.getConformer().setAtomPos(0,RDGeom::Point3D(1,0,0));
    m.getConformer().setAtomPos(1,RDGeom::Point3D(-1,0,0));
    m.getConformer().setAtomPos(2,RDGeom::Point3D(0,1,0));
    m.getConformer().setAtomPos(3,RDGeom::Point3D(0,-1,0));
    m.getConformer().setAtomPos(4,RDGeom::Point3D(0,0,1));
    m.getConformer().setAtomPos(5,RDGeom::Point3D(0,0,-1));
    double val;
    val = RDKit::Descriptors::radiusOfGyration(m);
    TEST_ASSERT(fabs(val)>0.1);
    val = RDKit::Descriptors::inertialShapeFactor(m);
    TEST_ASSERT(fabs(val)>1);
    val = RDKit::Descriptors::eccentricity(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::asphericity(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::spherocityIndex(m);
    TEST_ASSERT(fabs(1.-val)<1e-4);
  }

  {
    RDKit::RWMol m;
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addAtom(new RDKit::Atom(6));
    m.addConformer(new RDKit::Conformer(m.getNumAtoms()));
    double val;
    val = RDKit::Descriptors::radiusOfGyration(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::inertialShapeFactor(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::eccentricity(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::asphericity(m);
    TEST_ASSERT(fabs(val)<1e-4);
    val = RDKit::Descriptors::spherocityIndex(m);
    TEST_ASSERT(fabs(val)<1e-4);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}



//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
  test3DEdges();
  testPMIEdges();
  testNPREdges();
  testPMI1();
  testNPR1();

}
