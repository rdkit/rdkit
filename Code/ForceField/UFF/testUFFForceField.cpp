// $Id$
//
// Copyright (C)  2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <iomanip>
#include <math.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <Geometry/point.h>

#include <ForceField/ForceField.h>
#include <ForceField/UFF/Params.h>
#include <ForceField/UFF/BondStretch.h>
#include <ForceField/UFF/AngleBend.h>
#include <ForceField/UFF/Nonbonded.h>
#include <ForceField/UFF/TorsionAngle.h>
#include <ForceField/UFF/DistanceConstraint.h>
#include <ForceField/UFF/AngleConstraint.h>
#include <ForceField/UFF/TorsionConstraint.h>
#include <ForceField/UFF/PositionConstraint.h>

#include <GraphMol/Atom.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

using namespace RDGeom;

void test1(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for force field basics." << std::endl;
  
  ForceFields::ForceField ff;
  TEST_ASSERT(ff.dimension() ==3 );

  Point3D p1(0,0,0),p2(1,0,0),p3(2,0,0),p4(0,1,0);
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);

#if 0
  Point3D f1,f2,f3,f4;
  RDGeom::PointPtrVect &fs=ff.forces();
  fs.push_back(&f1);
  fs.push_back(&f2);
  fs.push_back(&f3);
  fs.push_back(&f4);
#endif
  TEST_ASSERT(ff.positions().size()==4);
  //TEST_ASSERT(ff.forces().size()==4);

  ff.initialize();

  TEST_ASSERT(RDKit::feq(ff.distance(0,1),1.0));
  TEST_ASSERT(RDKit::feq(ff.distance(1,0),1.0));
  TEST_ASSERT(RDKit::feq(ff.distance(0,0),0.0));
  TEST_ASSERT(RDKit::feq(ff.distance(0,2),2.0));
  TEST_ASSERT(RDKit::feq(ff.distance(2,0),2.0));
  TEST_ASSERT(RDKit::feq(ff.distance(0,3),1.0));
  TEST_ASSERT(RDKit::feq(ff.distance(3,0),1.0));
  TEST_ASSERT(RDKit::feq(ff.distance(3,3),0.0));
  TEST_ASSERT(RDKit::feq(ff.distance(1,2),1.0));
  TEST_ASSERT(RDKit::feq(ff.distance(2,1),1.0));
  
  std::cerr << "  done" << std::endl;
}


void testUFF1(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for basics of UFF bond-stretch terms." << std::endl;
  
  ForceFields::UFF::AtomicParams p1,p2;
  double restLen,forceConstant;

  // sp3 carbon:
  p1.r1 = .757;
  p1.Z1 = 1.912;
  p1.GMP_Xi = 5.343;
  
  // sp3 - sp3: checks basics
  restLen=ForceFields::UFF::Utils::calcBondRestLength(1.0,&p1,&p1);
  TEST_ASSERT(RDKit::feq(restLen,1.514));
  
  forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p1,&p1);
  TEST_ASSERT(RDKit::feq(forceConstant,699.5918));
  
  // sp2 carbon:
  p2.r1 = .732;
  p2.Z1 = 1.912;
  p2.GMP_Xi = 5.343;
  // sp2 - sp2: checks rBO
  restLen=ForceFields::UFF::Utils::calcBondRestLength(2.0,&p2,&p2);
  TEST_ASSERT(RDKit::feq(restLen,1.32883,1e-5));

  forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p2,&p2);
  TEST_ASSERT(RDKit::feq(forceConstant,1034.69,1e-2));
  
  // sp3 nitrogen:
  p2.r1 = .700;
  p2.Z1 = 2.544;
  p2.GMP_Xi = 6.899;

  // Csp3 - Nsp3: checks rEN
  restLen=ForceFields::UFF::Utils::calcBondRestLength(1.0,&p1,&p2);
  TEST_ASSERT(RDKit::feq(restLen,1.451071,1e-5));

  forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p1,&p2);
  TEST_ASSERT(RDKit::feq(forceConstant,1057.27,1e-2));


  // amide bond: check we can reproduce values from the UFF paper:
  // C_R:
  p1.r1 = .729;
  p1.Z1 = 1.912;
  p1.GMP_Xi = 5.343;
  // N_R:
  p2.r1 = .699;
  p2.Z1 = 2.544;
  p2.GMP_Xi = 6.899;

  restLen=ForceFields::UFF::Utils::calcBondRestLength(ForceFields::UFF::Params::amideBondOrder,&p1,&p2);
  TEST_ASSERT(RDKit::feq(restLen,1.357,1e-3)); // NOTE: the paper has 1.366


  forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p1,&p2);
  TEST_ASSERT(RDKit::feq(forceConstant,1293.,1)); // NOTE: the paper has 1293


  std::cerr << "  done" << std::endl;
}


void testUFF2(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for UFF bond-stretch terms." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1(0,0,0),p2(1.514,0,0);
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);

  ForceFields::UFF::AtomicParams param1;
  // sp3 carbon:
  param1.r1 = .757;
  param1.Z1 = 1.912;
  param1.GMP_Xi = 5.343;

  // C_3 - C_3, r0=1.514, k01=699.5918
  ForceFields::ForceFieldContrib *bs;
  bs = new ForceFields::UFF::BondStretchContrib(&ff,0,1,1,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(bs));
  ff.initialize();

  double *p,*g;
  p = new double[6];
  g = new double[6];
  for(int i=0;i<6;i++){
    p[i] = 0.0;
    g[i] = 0.0;
  }

  double E;
  // edge case: zero bond length:
  E=bs->getEnergy(p);
  TEST_ASSERT(E>0.0);
  bs->getGrad(p,g);
  for(int i=0;i<6;i++){
    TEST_ASSERT(fabs(g[i])>0.0);
  }

  p[0] = 0;
  p[3] = 1.514;
  for(int i=0;i<6;i++){
    g[i] = 0.0;
  }
  ff.initialize();
  E=bs->getEnergy(p);
  TEST_ASSERT(RDKit::feq(E,0.0));
  bs->getGrad(p,g);
  for(int i=0;i<6;i++){
    TEST_ASSERT(RDKit::feq(g[i],0.0));
  }

  (*ff.positions()[1])[0] = 1.814;
  p[3] = 1.814;
  ff.initialize();
  E=bs->getEnergy(p);
  TEST_ASSERT(RDKit::feq(E,31.4816));
  bs->getGrad(p,g);
  TEST_ASSERT(RDKit::feq(g[0],-209.8775));
  TEST_ASSERT(RDKit::feq(g[3],209.8775));
  TEST_ASSERT(RDKit::feq(g[1],0.0));
  TEST_ASSERT(RDKit::feq(g[2],0.0));
  TEST_ASSERT(RDKit::feq(g[4],0.0));
  TEST_ASSERT(RDKit::feq(g[5],0.0));

  // try a different axis:
  for(int i=0;i<6;i++){
    g[i] = 0.0;
    p[i] = 0.0;
  }
  ff.initialize();
  (*ff.positions()[1])[0] = 0.0;
  (*ff.positions()[1])[2] = 1.814;
  p[5] = 1.814;
  E=bs->getEnergy(p);
  TEST_ASSERT(RDKit::feq(E,31.4816));
  bs->getGrad(p,g);
  TEST_ASSERT(RDKit::feq(g[2],-209.8775));
  TEST_ASSERT(RDKit::feq(g[5],209.8775));
  TEST_ASSERT(RDKit::feq(g[0],0.0));
  TEST_ASSERT(RDKit::feq(g[1],0.0));
  TEST_ASSERT(RDKit::feq(g[3],0.0));
  TEST_ASSERT(RDKit::feq(g[4],0.0));

  // try a bit of minimization
  RDGeom::Point3D d;
  ff.initialize();
  (*ff.positions()[1])[2] = 0.0;
  (*ff.positions()[1])[0] = 1.814;
  ff.minimize(10,1e-8);
  d=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[1];
  TEST_ASSERT(RDKit::feq(d.length(),1.514,1e-3));
  
  // minimize in "3D"
  ff.initialize();
  (*ff.positions()[1])[2] = 1.1;
  (*ff.positions()[1])[1] = 0.9;
  (*ff.positions()[1])[0] = 1.00;
  ff.minimize(10,1e-8);
  d=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[1];
  TEST_ASSERT(RDKit::feq(d.length(),1.514,1e-3));
  
  

  delete [] p;
  delete [] g;
  std::cerr << "  done" << std::endl;
}


void testUFF3(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for basics of UFF angle terms." << std::endl;
  
  ForceFields::UFF::AtomicParams p1,p2,p3;
  double restLen,forceConstant;

  // sp3 carbon:
  p3.r1 = .757;
  p3.Z1 = 1.912;
  p3.GMP_Xi = 5.343;
  p3.theta0 = 109.47*M_PI/180.0;
  
  // sp3 - sp3: checks basics
  restLen=ForceFields::UFF::Utils::calcBondRestLength(1.0,&p3,&p3);
  TEST_ASSERT(RDKit::feq(restLen,1.514));
  
  // C_3 - C_3 - C_3
  forceConstant=ForceFields::UFF::Utils::calcAngleForceConstant(p3.theta0,1,1,&p3,&p3,&p3);
  //TEST_ASSERT(RDKit::feq(forceConstant,699.5918));
  
  // amide bond bend:
  // C_R - N_R - C_3
  // C_R:
  p1.r1 = .729;
  p1.Z1 = 1.912;
  p1.GMP_Xi = 5.343;
  // N_R:
  p2.r1 = .699;
  p2.Z1 = 2.544;
  p2.GMP_Xi = 6.899;
  p2.theta0 = 120.0*M_PI/180.;
  restLen=ForceFields::UFF::Utils::calcBondRestLength(ForceFields::UFF::Params::amideBondOrder,&p1,&p2);
  TEST_ASSERT(RDKit::feq(restLen,1.357,1e-3));
  restLen=ForceFields::UFF::Utils::calcBondRestLength(1.0,&p2,&p3);
  TEST_ASSERT(RDKit::feq(restLen,1.450,1e-3));

  forceConstant=ForceFields::UFF::Utils::calcAngleForceConstant(p2.theta0,ForceFields::UFF::Params::amideBondOrder,1,
						   &p1,&p2,&p3);
  TEST_ASSERT(RDKit::feq(forceConstant,211.0,1e-1)); //  paper has 105.5


  std::cerr << "  done" << std::endl;
}

void testUFF4(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for UFF angle-bend terms." << std::endl;


  ForceFields::ForceField ff;
  Point3D p1(1.514,0,0),p2(0,0,0),p3(0.1,1.5,0);
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);

  ForceFields::UFF::AtomicParams param1;
  // sp3 carbon:
  param1.r1 = .757;
  param1.Z1 = 1.912;
  param1.GMP_Xi = 5.343;
  // cheat to get the angle to 90 so that testing is easier:
  param1.theta0 = 90.0*M_PI/180.;
  
  // C_3 - C_3, r0=1.514, k01=699.5918
  ForceFields::ForceFieldContrib *contrib;
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,1,1,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,2,1,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,2,1,1,&param1,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));


  RDGeom::Point3D d,v1,v2;
  double theta;
#if 1
  // ------- ------- ------- ------- ------- ------- -------
  // try a bit of minimization
  ff.initialize();
  ff.minimize(10,1e-8,1e-8);

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[1]-*(RDGeom::Point3D*)ff.positions()[2];
  theta = v1.angleTo(v2);

  TEST_ASSERT(RDKit::feq(v1.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(theta,90*M_PI/180.,1e-4));

  // ------- ------- ------- ------- ------- ------- -------
  // more complicated atomic coords:
  p1.x=1.3;
  p1.y=0.1;
  p1.z=0.1;
  p2.x=-0.1;
  p2.y=0.05;
  p2.z=-0.05;
  p3.x=0.1;
  p3.y=1.5;
  p3.z=0.05;
  ff.initialize();
  ff.minimize(10,1e-8,1e-8);

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[1]-*(RDGeom::Point3D*)ff.positions()[2];
  theta = v1.angleTo(v2);

  TEST_ASSERT(RDKit::feq(v1.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(theta,90*M_PI/180.,1e-4));
    
  // ------- ------- ------- ------- ------- ------- -------
  // try for the tetrahedral angle instead of 90:
  param1.theta0 = 109.47*M_PI/180.;
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,2,1,1,&param1,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  
  p1.x=1.3;
  p1.y=0.1;
  p1.z=0.1;
  p2.x=-0.1;
  p2.y=0.05;
  p2.z=-0.05;
  p3.x=0.1;
  p3.y=1.5;
  p3.z=0.05;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);
  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[2]-*(RDGeom::Point3D*)ff.positions()[1];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(v1.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));

#endif  

  // ------- ------- ------- ------- ------- ------- -------
  //
  // Do a series of "special cases" (i.e. test the functional forms
  // for linear, trigonal planar, square planar and octahedral)
  //
  // ------- ------- ------- ------- ------- ------- -------

  // ------- ------- ------- ------- ------- ------- -------
  // test a linear molecule:
  param1.theta0 = M_PI;
  //ff.contribs().pop_back();
  //ff.contribs().pop_back();
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,2,1,1,&param1,&param1,&param1,2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  
  p1.x=1.3;
  p1.y=0.1;
  p1.z=0.0;
  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;
  p3.x=-1.3;
  p3.y=0.1;
  p3.z=0.00;
  ff.initialize();
  ff.minimize(100,1e-8,1e-8);

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[2]-*(RDGeom::Point3D*)ff.positions()[1];
  theta = v1.angleTo(v2);

  TEST_ASSERT(RDKit::feq(v1.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),1.514,1e-3));
  std::cerr << "theta = " << theta << "; theta0 = " << param1.theta0 << std::endl;
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));


  // ------- ------- ------- ------- ------- ------- -------
  // test n=3:
  param1.theta0 = 120.*M_PI/180.0;
  //ff.contribs().pop_back();
  //ff.contribs().pop_back();
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,2,1,1,&param1,&param1,&param1,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  
  p1.x=1.3;
  p1.y=0.1;
  p1.z=0.0;
  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;
  p3.x=-.3;
  p3.y=-1.3;
  p3.z=0.00;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[2]-*(RDGeom::Point3D*)ff.positions()[1];
  theta = v1.angleTo(v2);

  TEST_ASSERT(RDKit::feq(v1.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),1.514,1e-3));
  std::cerr << "theta = " << std::fixed << std::setprecision(6) << theta
    << ", param1.theta0 = " << std::fixed << std::setprecision(6) << param1.theta0
    << std::endl;
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));


  // ------- ------- ------- ------- ------- ------- -------
  // test n=4:
  param1.theta0 = M_PI/2.0;
  //ff.contribs().pop_back();
  //ff.contribs().pop_back();
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,2,1,1,&param1,&param1,&param1,4);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  
  p1.x=1.3;
  p1.y=0.1;
  p1.z=0.0;
  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;
  p3.x=-.3;
  p3.y=-1.3;
  p3.z=0.00;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[2]-*(RDGeom::Point3D*)ff.positions()[1];
  theta = v1.angleTo(v2);

  TEST_ASSERT(RDKit::feq(v1.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),1.514,1e-3));
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));


  
#if 0
  std::cerr << " " << *ff.positions()[0] << std::endl;
  std::cerr << " " << *ff.positions()[1] << std::endl;
  std::cerr << " " << *ff.positions()[2] << std::endl;

  std::cerr << "v1: " << v1 << std::endl;
  std::cerr << "v2: " << v2 << std::endl;
  std::cerr << "FINAL: " << v1.angleTo(v2) << " " << v1.signedAngleTo(v2) << std::endl;
#endif

  std::cerr << "  done" << std::endl;
}


void testUFF5(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test Simple UFF molecule optimizations." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1,p2,p3,p4,p5,p6;
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);
  ps.push_back(&p5);
  ps.push_back(&p6);

  ForceFields::UFF::AtomicParams param1,param2;
  // sp2 carbon:
  param1.r1 = .732;
  param1.Z1 = 1.912;
  param1.GMP_Xi = 5.343;
  param1.theta0 = 120.*M_PI/180.;
  
  // H_1:
  param2.r1 = 0.354;
  param2.Z1 = 0.712;
  param2.GMP_Xi = 4.528;

  ForceFields::ForceFieldContrib *contrib;

  // build ethylene:
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,1,2,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,2,1,&param1,&param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,3,1,&param1,&param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,4,1,&param1,&param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,5,1,&param1,&param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,0,2,2,1,&param1,&param1,&param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,0,3,2,1,&param1,&param1,&param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,2,0,3,1,1,&param2,&param1,&param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,4,2,1,&param1,&param1,&param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,5,2,1,&param1,&param1,&param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,4,1,5,1,1,&param2,&param1,&param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  // dodge the fact that we're not using torsions yet by putting
  // everything in the z=0 plane:
  p1.x=-0.58;
  p1.y=-0.33;
  p1.z=0.0;

  p2.x=0.58;
  p2.y=0.33;
  p2.z=0.0;

  p3.x=-0.61;
  p3.y=-1.43;
  p3.z=0.0;

  p4.x=-1.54;
  p4.y=0.20;
  p4.z=0.0;

  p5.x=0.61;
  p5.y=1.43;
  p5.z=0.0;

  p6.x=1.54;
  p6.y=-0.20;
  p6.z=0.0;

  RDGeom::Point3D d,v1,v2;
  double theta;
  // ------- ------- ------- ------- ------- ------- -------
  // try a bit of minimization
  ff.initialize();
  ff.minimize(10,1e-8,1e-8);

  double CCDblBondLen=ForceFields::UFF::Utils::calcBondRestLength(2,&param1,&param1);
  double CHBondLen=ForceFields::UFF::Utils::calcBondRestLength(1,&param1,&param2);

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[2];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(v1.length(),CCDblBondLen,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),CHBondLen,1e-3));
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));
  v2=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[3];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(v2.length(),CHBondLen,1e-3));
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));

  v1=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[2];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(theta,param1.theta0,1e-4));
    
  std::cerr << "  done" << std::endl;
}

void testUFF6(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for UFF nonbonded terms." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1(0,0,0),p2(0.0,0,0);
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);

  ForceFields::UFF::AtomicParams param1;
  // sp3 carbon:
  param1.r1 = .757;
  param1.Z1 = 1.912;
  param1.GMP_Xi = 5.343;
  param1.x1 = 3.851;
  param1.D1 = 0.105;

  ff.initialize();
  ForceFields::ForceFieldContrib *contrib;
  contrib = new ForceFields::UFF::vdWContrib(&ff,0,1,&param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  // try a bit of minimization
  RDGeom::Point3D d;
  ff.initialize();

  // edge case: our energy at zero length should be zero:
  double E;
  E=ff.calcEnergy();
  TEST_ASSERT(RDKit::feq(E,0.0));
  
  (*ff.positions()[0])[0] = 0.0;
  (*ff.positions()[1])[0] = 4.0;
  ff.minimize(10,1e-8,1e-8);
  d=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[1];
  TEST_ASSERT(RDKit::feq(d.length(),3.851,1e-3));
  
  // minimize in "3D"
  ff.initialize();
  (*ff.positions()[0])[0] = 0.0;
  (*ff.positions()[0])[1] = 0.0;
  (*ff.positions()[0])[2] = 0.0;
  (*ff.positions()[1])[2] = 3.1;
  (*ff.positions()[1])[1] = 0.9;
  (*ff.positions()[1])[0] = 1.00;
  ff.minimize(10,1e-8,1e-8);
  d=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  TEST_ASSERT(RDKit::feq(d.length(),3.851,1e-3));
  
  std::cerr << "  done" << std::endl;
}

void testUFF7(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test UFF torsional terms." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1,p2,p3,p4;
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);

  ForceFields::UFF::AtomicParams param1,param2;
  // sp3 carbon:
  param1.r1 = .757;
  param1.Z1 = 1.912;
  param1.GMP_Xi = 5.343;
  param1.x1 = 3.851;
  param1.D1 = 0.105;
  param1.V1 = 2.119;
  param1.U1 = 2.0;
  
  // H_1:
  param2.r1 = 0.354;
  param2.Z1 = 0.712;
  param2.GMP_Xi = 4.528;

  RDGeom::Point3D d,v1,v2;
  double cosPhi;

  ForceFields::ForceFieldContrib *contrib;
  // ------- ------- ------- ------- ------- ------- -------
  // Basic SP3 - SP3
  // ------- ------- ------- ------- ------- ------- -------
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,3,1,
					 6,6,
					 RDKit::Atom::SP3,RDKit::Atom::SP3,
					 &param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
#if 1
  p1.x=0;
  p1.y=1.5;
  p1.z=0;

  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;

  p3.x=1.5;
  p3.y=0.0;
  p3.z=0.0;

  p4.x=1.5;
  p4.y=0.0;
  p4.z=1.5;

  ff.initialize();
  ff.minimize(10,1e-8,1e-8);
  cosPhi = ForceFields::UFF::Utils::calculateCosTorsion(*(RDGeom::Point3D*)ff.positions()[0],
							*(RDGeom::Point3D*)ff.positions()[1],
							*(RDGeom::Point3D*)ff.positions()[2],
							*(RDGeom::Point3D*)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi,0.5,1e-4));

  // ------- ------- ------- ------- ------- ------- -------
  // Basic SP2 - SP2
  // ------- ------- ------- ------- ------- ------- -------
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,3,1,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP2,
					 &param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  p1.x=0;
  p1.y=1.5;
  p1.z=0.1;

  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;

  p3.x=1.5;
  p3.y=0.0;
  p3.z=0.0;

  p4.x=1.5;
  p4.y=0.2;
  p4.z=1.5;

  ff.initialize();
  ff.minimize(10,1e-8,1e-8);
  cosPhi = ForceFields::UFF::Utils::calculateCosTorsion(*(RDGeom::Point3D*)ff.positions()[0],
							*(RDGeom::Point3D*)ff.positions()[1],
							*(RDGeom::Point3D*)ff.positions()[2],
							*(RDGeom::Point3D*)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi,1.0,1e-4));
  
  // ------- ------- ------- ------- ------- ------- -------
  // Basic SP2 - SP3
  // ------- ------- ------- ------- ------- ------- -------
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,3,1,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 &param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  p1.x=0;
  p1.y=1.5;
  p1.z=0.1;

  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;

  p3.x=1.5;
  p3.y=0.0;
  p3.z=0.0;

  p4.x=1.5;
  p4.y=0.2;
  p4.z=1.5;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);
  cosPhi = ForceFields::UFF::Utils::calculateCosTorsion(*(RDGeom::Point3D*)ff.positions()[0],
							*(RDGeom::Point3D*)ff.positions()[1],
							*(RDGeom::Point3D*)ff.positions()[2],
							*(RDGeom::Point3D*)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi,0.5,1e-4));

  // ------- ------- ------- ------- ------- ------- -------
  // special case for group 6 - group 6 bonds:
  // ------- ------- ------- ------- ------- ------- -------
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,3,1,
					 8,8,
					 RDKit::Atom::SP3,RDKit::Atom::SP3,
					 &param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  p1.x=0;
  p1.y=1.5;
  p1.z=0.1;

  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;

  p3.x=1.5;
  p3.y=0.0;
  p3.z=0.0;

  p4.x=1.5;
  p4.y=0.2;
  p4.z=1.5;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);
  cosPhi = ForceFields::UFF::Utils::calculateCosTorsion(*(RDGeom::Point3D*)ff.positions()[0],
							*(RDGeom::Point3D*)ff.positions()[1],
							*(RDGeom::Point3D*)ff.positions()[2],
							*(RDGeom::Point3D*)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi,0.0,1e-4));

  // ------- ------- ------- ------- ------- ------- -------
  // special case for SP3 group 6 - SP2 other group
  // ------- ------- ------- ------- ------- ------- -------
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,3,1,
					 8,6,
					 RDKit::Atom::SP3,RDKit::Atom::SP2,
					 &param1,&param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  p1.x=0;
  p1.y=1.5;
  p1.z=0.1;

  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;

  p3.x=1.5;
  p3.y=0.0;
  p3.z=0.0;

  p4.x=1.5;
  p4.y=0.2;
  p4.z=1.5;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);
  cosPhi = ForceFields::UFF::Utils::calculateCosTorsion(*(RDGeom::Point3D*)ff.positions()[0],
							*(RDGeom::Point3D*)ff.positions()[1],
							*(RDGeom::Point3D*)ff.positions()[2],
							*(RDGeom::Point3D*)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi,0.0,1e-4));

#endif

  // ------- ------- ------- ------- ------- ------- -------
  // special case for (SP2 -) SP2 - SP3
  // ------- ------- ------- ------- ------- ------- -------
  ff.contribs().pop_back();
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,3,1,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 &param1,&param1,true);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  p1.x=0;
  p1.y=1.5;
  p1.z=0.1;

  p2.x=0.0;
  p2.y=0.0;
  p2.z=0.0;

  p3.x=1.5;
  p3.y=0.0;
  p3.z=0.0;

  p4.x=1.5;
  p4.y=0.2;
  p4.z=1.5;

  ff.initialize();
  ff.minimize(100,1e-8,1e-8);
  cosPhi = ForceFields::UFF::Utils::calculateCosTorsion(*(RDGeom::Point3D*)ff.positions()[0],
							*(RDGeom::Point3D*)ff.positions()[1],
							*(RDGeom::Point3D*)ff.positions()[2],
							*(RDGeom::Point3D*)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi,0.5,1e-4));
  
  std::cerr << "  done" << std::endl;
}


void testUFFParams(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test UFF Parameter objects" << std::endl;

  ForceFields::UFF::ParamCollection *params=ForceFields::UFF::ParamCollection::getParams();
  TEST_ASSERT(params);

  const ForceFields::UFF::AtomicParams *ptr;
  ptr=(*params)("C_3");
  TEST_ASSERT(ptr);
  TEST_ASSERT(RDKit::feq(ptr->r1,0.757));
  TEST_ASSERT(RDKit::feq(ptr->theta0,109.47*M_PI/180.));
  TEST_ASSERT(RDKit::feq(ptr->x1,3.851));
  TEST_ASSERT(RDKit::feq(ptr->D1,0.105));
  TEST_ASSERT(RDKit::feq(ptr->zeta,12.73));
  TEST_ASSERT(RDKit::feq(ptr->Z1,1.912));
  TEST_ASSERT(RDKit::feq(ptr->V1,2.119));
  TEST_ASSERT(RDKit::feq(ptr->GMP_Xi,5.343));
  TEST_ASSERT(RDKit::feq(ptr->GMP_Hardness,5.063));
  TEST_ASSERT(RDKit::feq(ptr->GMP_Radius,0.759));
  
  ptr=(*params)("N_3");
  TEST_ASSERT(ptr);

  ptr=(*params)("C_5");
  TEST_ASSERT(!ptr);
}


void testUFF8(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test Simple UFF molecule optimization, part 2." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1,p2,p3,p4,p5,p6;
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);
  ps.push_back(&p5);
  ps.push_back(&p6);

  ForceFields::UFF::ParamCollection *params=ForceFields::UFF::ParamCollection::getParams();
  const ForceFields::UFF::AtomicParams *param1,*param2;

  // C_2 (sp2 carbon):
  param1 = (*params)("C_2");
  TEST_ASSERT(param1);
  // H_:
  param2 = (*params)("H_");
  TEST_ASSERT(param2);
  ForceFields::ForceFieldContrib *contrib;

  // build ethylene:
  // BONDS:
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,1,2,param1,param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,2,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,3,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,4,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,5,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  // ANGLES:
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,0,2,2,1,param1,param1,param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,0,3,2,1,param1,param1,param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,2,0,3,1,1,param2,param1,param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,4,2,1,param1,param1,param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,5,2,1,param1,param1,param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,4,1,5,1,1,param2,param1,param2,3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  // DIHEDRALS:
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,2,0,1,4,2,
					 6,6,
					 RDKit::Atom::SP3,RDKit::Atom::SP3,
					 param1,param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,2,0,1,5,2,
					 6,6,
					 RDKit::Atom::SP3,RDKit::Atom::SP3,
					 param1,param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,3,0,1,4,2,
					 6,6,
					 RDKit::Atom::SP3,RDKit::Atom::SP3,
					 param1,param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,3,0,1,5,2,
					 6,6,
					 RDKit::Atom::SP3,RDKit::Atom::SP3,
					 param1,param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  
  p1.x=-0.58;
  p1.y=-0.33;
  p1.z=0.1;

  p2.x=0.58;
  p2.y=0.33;
  p2.z=0.1;

  p3.x=-0.61;
  p3.y=-1.43;
  p3.z=0.0;

  p4.x=-1.54;
  p4.y=0.20;
  p4.z=0.0;

  p5.x=0.61;
  p5.y=1.43;
  p5.z=0.0;

  p6.x=1.54;
  p6.y=-0.20;
  p6.z=0.0;

  RDGeom::Point3D d,v1,v2;
  double theta;
  // ------- ------- ------- ------- ------- ------- -------
  // try a bit of minimization
  ff.initialize();
  ff.minimize(100,1e-8,1e-8);

  double CCDblBondLen=ForceFields::UFF::Utils::calcBondRestLength(2,param1,param1);
  double CHBondLen=ForceFields::UFF::Utils::calcBondRestLength(1,param1,param2);

  v1=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[1];
  v2=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[2];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(v1.length(),CCDblBondLen,1e-3));
  TEST_ASSERT(RDKit::feq(v2.length(),CHBondLen,1e-3));
  TEST_ASSERT(RDKit::feq(theta,param1->theta0,1e-4));
  v2=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[3];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(v2.length(),CHBondLen,1e-3));
  TEST_ASSERT(RDKit::feq(theta,param1->theta0,1e-4));

  v1=*(RDGeom::Point3D*)ff.positions()[0] - *(RDGeom::Point3D*)ff.positions()[2];
  theta = v1.angleTo(v2);
  TEST_ASSERT(RDKit::feq(theta,param1->theta0,1e-4));
    
  std::cerr << "  done" << std::endl;
}

void testUFFTorsionConflict(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test UFF Torsion Conflicts." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1,p2,p3,p4,p5,p6,p7;
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);
  ps.push_back(&p5);
  ps.push_back(&p6);
  ps.push_back(&p7);

  ForceFields::UFF::ParamCollection *params=ForceFields::UFF::ParamCollection::getParams();
  const ForceFields::UFF::AtomicParams *param1,*param2,*param3;

  // C_2 (sp2 carbon):
  param1 = (*params)("C_2");
  TEST_ASSERT(param1);
  // H_:
  param2 = (*params)("H_");
  TEST_ASSERT(param2);
  // C_3 (sp3 carbon):
  param3 = (*params)("C_3");
  TEST_ASSERT(param3);

  ForceFields::ForceFieldContrib *contrib;

  // BONDS:
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,0,1,2,param1,param1);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,2,1,param1,param3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,1,3,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,2,4,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,2,5,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::BondStretchContrib(&ff,2,6,1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

#if 1
  // ANGLES:
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,2,2.0,1.0,param1,param1,param3);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,0,1,3,2.0,1.0,param1,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,2,4,1.0,1.0,param1,param3,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,2,5,1.0,1.0,param1,param3,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,1,2,6,1.0,1.0,param1,param3,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::AngleBendContrib(&ff,2,1,3,1.0,1.0,param3,param1,param2);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
#endif

  // DIHEDRALS:
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,4,1.0,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 param1,param3,true);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,3,1,2,4,1.0,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 param1,param3,false);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,5,1.0,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 param1,param3,true);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,3,1,2,5,1.0,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 param1,param3,false);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,0,1,2,6,1.0,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 param1,param3,true);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));
  contrib = new ForceFields::UFF::TorsionAngleContrib(&ff,3,1,2,6,1.0,
					 6,6,
					 RDKit::Atom::SP2,RDKit::Atom::SP3,
					 param1,param3,false);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));


  
  p1.x=0.5411;
  p1.y=-0.7741;
  p1.z=0.0902;

  p2.x=-0.5622;
  p2.y=-0.0368;
  p2.z=0.1202;

  p3.x=-0.5101;
  p3.y=1.4485;
  p3.z=0.0816;

  p4.x=-1.5285;
  p4.y=-0.5341;
  p4.z=0.1892;

  p5.x=0.5097;
  p5.y=1.8065;
  p5.z=0.1988;

  p6.x=-1.1436;
  p6.y=1.8781;
  p6.z=0.8983;

  p7.x=-0.9145;
  p7.y=1.8185;
  p7.z=-0.8983;

  RDGeom::Point3D d,v1,v2;
  // ------- ------- ------- ------- ------- ------- -------
  // try a bit of minimization
  ff.initialize();
  ff.minimize(100,1e-8,1e-8);
    
#if 1
  std::cerr.setf(std::ios_base::fixed,std::ios_base::floatfield);
  std::cerr.precision(4);
  std::cerr << "C " << *ff.positions()[0] << std::endl;
  std::cerr << "C " << *ff.positions()[1] << std::endl;
  std::cerr << "C " << *ff.positions()[2] << std::endl;
  std::cerr << "H " << *ff.positions()[3] << std::endl;
  std::cerr << "H " << *ff.positions()[4] << std::endl;
  std::cerr << "O " << *ff.positions()[5] << std::endl;
  std::cerr << "F " << *ff.positions()[6] << std::endl;
#endif
  std::cerr << "  done" << std::endl;
}

void testUFFDistanceConstraints(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for UFF distance constraint terms." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1(0,0,0),p2(1.514,0,0);
  RDGeom::PointPtrVect &ps=ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);

  double *p,*g;
  p = new double[6];
  g = new double[6];
  for(int i=0;i<6;i++){
    p[i] = 0.0;
    g[i] = 0.0;
  }
  p[0] = 0;
  p[3] = 1.40;

  ff.initialize();

  // C_3 - C_3, r0=1.514, k01=699.5918
  ForceFields::ForceFieldContrib *bs;
  bs = new ForceFields::UFF::DistanceConstraintContrib(&ff,0,1,1.35,1.55,1000.0);
  ff.contribs().push_back(ForceFields::ContribPtr(bs));
  double E;
  E=bs->getEnergy(p);
  TEST_ASSERT(RDKit::feq(E,0.0));
  bs->getGrad(p,g);
  for(int i=0;i<6;i++){
    TEST_ASSERT(RDKit::feq(g[i],0.0));
  }

  ff.initialize();
  (*ff.positions()[1])[0] = 1.20;
  p[3] = 1.20;
  E=bs->getEnergy(p);
  TEST_ASSERT(RDKit::feq(E,11.25));
  bs->getGrad(p,g);
  TEST_ASSERT(RDKit::feq(g[0],150.0));
  TEST_ASSERT(RDKit::feq(g[3],-150.0));
  TEST_ASSERT(RDKit::feq(g[1],0.0));
  TEST_ASSERT(RDKit::feq(g[2],0.0));
  TEST_ASSERT(RDKit::feq(g[4],0.0));
  TEST_ASSERT(RDKit::feq(g[5],0.0));

  // try a bit of minimization
  RDGeom::Point3D d;
  ff.initialize();
  (*ff.positions()[1])[2] = 0.0;
  (*ff.positions()[1])[0] = 1.20;
  ff.minimize(10,1e-8);
  d=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  TEST_ASSERT(d.length()>=1.35)
  TEST_ASSERT(d.length()<=1.55)
  
  ff.initialize();
  (*ff.positions()[1])[2] = 0.0;
  (*ff.positions()[1])[0] = 1.70;
  ff.minimize(10,1e-8);
  d=*(RDGeom::Point3D*)ff.positions()[0]-*(RDGeom::Point3D*)ff.positions()[1];
  TEST_ASSERT(d.length()>=1.35)
  TEST_ASSERT(d.length()<=1.55)
  
  delete [] p;
  delete [] g;
  std::cerr << "  done" << std::endl;
}

void testUFFAllConstraints(){
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << "Unit tests for all UFF constraint terms." << std::endl;

  std::string molBlock =
    "butane\n"
    "     RDKit          3D\n"
    "butane\n"
    " 17 16  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.4280    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.7913   -0.2660    0.9927 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.9040    1.3004   -0.3485 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.5407    2.0271    0.3782 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.5407    1.5664   -1.3411 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.3320    1.3004   -0.3485 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6953    1.5162   -1.3532 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.8080    0.0192    0.0649 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.4447   -0.7431   -0.6243 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.4447   -0.1966    1.0697 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    4.8980    0.0192    0.0649 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    3.6954    2.0627    0.3408 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.7913   -0.7267   -0.7267 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3633    0.7267    0.7267 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3633   -0.9926    0.2660 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "   -0.3633    0.2660   -0.9926 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0  0  0  0\n"
    "  1 15  1  0  0  0  0\n"
    "  1 16  1  0  0  0  0\n"
    "  1 17  1  0  0  0  0\n"
    "  2  3  1  0  0  0  0\n"
    "  2  4  1  0  0  0  0\n"
    "  2 14  1  0  0  0  0\n"
    "  4  5  1  0  0  0  0\n"
    "  4  6  1  0  0  0  0\n"
    "  4  7  1  0  0  0  0\n"
    "  7  8  1  0  0  0  0\n"
    "  7  9  1  0  0  0  0\n"
    "  7 13  1  0  0  0  0\n"
    "  9 10  1  0  0  0  0\n"
    "  9 11  1  0  0  0  0\n"
    "  9 12  1  0  0  0  0\n"
    "M  END\n";
  RDKit::RWMol *mol;
  ForceFields::ForceField *field;
  
  // distance constraints
  ForceFields::UFF::DistanceConstraintContrib *dc;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  field = RDKit::UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  dc = new ForceFields::UFF::DistanceConstraintContrib(field, 1, 3, 2.0, 2.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(dc));
  field->minimize();
  TEST_ASSERT(MolTransforms::getBondLength(mol->getConformer(), 1, 3) > 1.99);
  delete field;
  field = RDKit::UFF::constructForceField(*mol);
  field->initialize();
  dc = new ForceFields::UFF::DistanceConstraintContrib(field, 1, 3, true, -0.2, 0.2, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(dc));
  field->minimize();
  TEST_ASSERT(MolTransforms::getBondLength(mol->getConformer(), 1, 3) > 1.79);
  delete field;
  delete mol;
  
  // angle constraints
  ForceFields::UFF::AngleConstraintContrib *ac;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  field = RDKit::UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  ac = new ForceFields::UFF::AngleConstraintContrib(field, 1, 3, 6, 90.0, 90.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getAngleDeg(mol->getConformer(), 1, 3, 6) == 90);
  delete field;
  field = RDKit::UFF::constructForceField(*mol);
  field->initialize();
  ac = new ForceFields::UFF::AngleConstraintContrib(field, 1, 3, 6, true, -10.0, 10.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getAngleDeg(mol->getConformer(), 1, 3, 6) == 100);
  delete field;
  field = RDKit::UFF::constructForceField(*mol);
  field->initialize();
  ac = new ForceFields::UFF::AngleConstraintContrib(field, 1, 3, 6, false, -10.0, 10.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getAngleDeg(mol->getConformer(), 1, 3, 6) == 10);
  delete field;
  delete mol;
  
  // torsion constraints
  ForceFields::UFF::TorsionConstraintContrib *tc;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  MolTransforms::setDihedralDeg(mol->getConformer(), 1, 3, 6, 8, 50.0);
  field = RDKit::UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  tc = new ForceFields::UFF::TorsionConstraintContrib(field, 1, 3, 6, 8, 10.0, 20.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(tc));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getDihedralDeg(mol->getConformer(), 1, 3, 6, 8) == 20);
  delete field;
  MolTransforms::setDihedralDeg(mol->getConformer(), 1, 3, 6, 8, -30.0);
  field = RDKit::UFF::constructForceField(*mol);
  field->initialize();
  tc = new ForceFields::UFF::TorsionConstraintContrib(field, 1, 3, 6, 8, true, -10.0, -8.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(tc));
  field->minimize();
  TEST_ASSERT((int)MolTransforms::getDihedralDeg(mol->getConformer(), 1, 3, 6, 8) == -40);
  delete field;
  MolTransforms::setDihedralDeg(mol->getConformer(), 1, 3, 6, 8, 30.0);
  field = RDKit::UFF::constructForceField(*mol);
  field->initialize();
  tc = new ForceFields::UFF::TorsionConstraintContrib(field, 1, 3, 6, 8, false, -10.0, -8.0, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(tc));
  field->minimize(500);
  TEST_ASSERT((int)MolTransforms::getDihedralDeg(mol->getConformer(), 1, 3, 6, 8) == -10);
  delete field;
  delete mol;
  
  // position constraints
  ForceFields::UFF::PositionConstraintContrib *pc;
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  field = RDKit::UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  RDGeom::Point3D p = mol->getConformer().getAtomPos(1);
  pc = new ForceFields::UFF::PositionConstraintContrib(field, 1, 0.3, 1.0e5);
  field->contribs().push_back(ForceFields::ContribPtr(pc));
  field->minimize();
  RDGeom::Point3D q = mol->getConformer().getAtomPos(1);
  TEST_ASSERT((p - q).length() < 0.3);
  delete field;
  delete mol;
  
  // fixed atoms
  mol = RDKit::MolBlockToMol(molBlock, true, false);
  TEST_ASSERT(mol);
  field = RDKit::UFF::constructForceField(*mol);
  TEST_ASSERT(field);
  field->initialize();
  RDGeom::Point3D fp = mol->getConformer().getAtomPos(1);
  field->fixedPoints().push_back(1);
  field->minimize();
  RDGeom::Point3D fq = mol->getConformer().getAtomPos(1);
  TEST_ASSERT((fp - fq).length() < 0.01);
  delete field;
  delete mol;

  std::cerr << "  done" << std::endl;
}

int main(){
#if 1
  test1();
  testUFF1();
  testUFF2();
  testUFF3();
  testUFF4();
  testUFF5();
  testUFF6();
  testUFF7();
  testUFFParams();
  testUFF8();
  testUFFTorsionConflict();
  
#endif
  testUFFDistanceConstraints();
  testUFFAllConstraints();
}
