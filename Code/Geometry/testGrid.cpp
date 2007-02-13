//  $Id$
// 
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//


#include "UniformGrid3D.h"
#include <DataStructs/DiscreteValueVect.h>
#include "point.h"
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include "GridUtils.h"

using namespace RDGeom;
using namespace RDKit;

void testUniformGrid1() {
  UniformGrid3D grd(6.0, 5.0, 4.0);
  CHECK_INVARIANT(grd.getSize() == 960, "");
  CHECK_INVARIANT(RDKit::feq(grd.getSpacing(), .5), "");
  CHECK_INVARIANT(grd.getNumX()==12, "");
  CHECK_INVARIANT(grd.getNumY()==10, "");
  CHECK_INVARIANT(grd.getNumZ()==8, "");
  
  grd.setSphereOccupancy(Point3D(0.0, 0.0, 0.0), 1.5, 0.25);
  CHECK_INVARIANT(grd.getOccupancyVect()->getTotalVal() == 523, "" ); 
  //writeGridToFile(grd, "junk.grd");
}

void testUniformGrid2() {
  // test tanimoto distance 
  UniformGrid3D grd(10.0, 10.0, 10.0);
  grd.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);
  writeGridToFile(grd, "junk.grd");
  double dist = tanimotoDistance(grd, grd);
  CHECK_INVARIANT(RDKit::feq(dist, 0.0), "");

  UniformGrid3D grd2(10.0, 10.0, 10.0);
  grd2.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd2.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  grd2.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd, grd2);
  CHECK_INVARIANT(RDKit::feq(dist, 0.25), "");
  dist = protrudeDistance(grd, grd2);
  CHECK_INVARIANT(RDKit::feq(dist, 0.25), "");
  dist = protrudeDistance(grd2, grd);
  CHECK_INVARIANT(RDKit::feq(dist, 0.0), "");

  UniformGrid3D grd3(10.0, 10.0, 10.0);
  grd3.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd3.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd, grd3);
  CHECK_INVARIANT(RDKit::feq(dist, 0.5), "");
  dist = protrudeDistance(grd, grd3);
  CHECK_INVARIANT(RDKit::feq(dist, 0.5), "");
  dist = protrudeDistance(grd3, grd);
  CHECK_INVARIANT(RDKit::feq(dist, 0.0), "");
  
  UniformGrid3D grd4(10.0, 10.0, 10.0);
  grd4.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);
  grd4.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd3, grd4);
  CHECK_INVARIANT(RDKit::feq(dist, 1.0), "");


  UniformGrid3D grd5(10.0, 10.0, 10.0);
  grd5.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  dist = tanimotoDistance(grd, grd5);
  CHECK_INVARIANT(RDKit::feq(dist, 0.75), "");
  dist = protrudeDistance(grd, grd5);
  CHECK_INVARIANT(RDKit::feq(dist, 0.75), "");
  dist = protrudeDistance(grd5, grd);
  CHECK_INVARIANT(RDKit::feq(dist, 0.00), "");
  
 
}

void testUniformGridPickling() {
  // test tanimoto distance 
  UniformGrid3D grd(10.0, 10.0, 10.0);
  grd.setSphereOccupancy(Point3D(-2.0, -2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(-2.0, 2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(2.0, -2.0, 0.0), 1.5, 0.25);
  grd.setSphereOccupancy(Point3D(2.0, 2.0, 0.0), 1.5, 0.25);
  UniformGrid3D grd2(grd.toString());
  double dist = tanimotoDistance(grd, grd2);
  CHECK_INVARIANT(RDKit::feq(dist, 0.0), "");

}

int main() {
  std::cout << "***********************************************************\n";
  std::cout << "Testing Grid\n";

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testUniformGrid1 \n\n";
  testUniformGrid1();

  std::cout << "\t---------------------------------\n";
  std::cout << "\t testUniformGrid2 \n\n";
  testUniformGrid2();


  std::cout << "\t---------------------------------\n";
  std::cout << "\t testUniformGridPickling \n\n";
  testUniformGridPickling();
  return 0;
}
