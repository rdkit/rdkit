//
// Copyright (C)  2015 Sereina Riniker
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <RDGeneral/utils.h>
#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/TorsionAngle.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleM6.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <ForceField/ForceField.h>
#include <iostream>
#include <string>
#include <math.h>

using namespace RDGeom;
using namespace RDKit;

void testTorsionAngleM6() {
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test CrystalFF torsional term." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1, p2, p3, p4;
  RDGeom::PointPtrVect &ps = ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);

  ForceFields::CrystalFF::TorsionAngleContribM6 *contrib;
  // ------- ------- ------- ------- ------- ------- -------
  // Basic SP3 - SP3
  // ------- ------- ------- ------- ------- ------- -------

  // [!#1:1][CX4H2:2]!@[CX4H2:3][!#1:4] 1 0.0 1 0.0 1 4.0 1 0.0 1 0.0 1 0.0
  std::vector<int> signs(6, 1);
  std::vector<double> v(6, 0.0);
  v[2] = 4.0;

  contrib = new ForceFields::CrystalFF::TorsionAngleContribM6(&ff, 0, 1, 2, 3,
                                                              v, signs);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  p1.x = 0;
  p1.y = 1.5;
  p1.z = 0;

  p2.x = 0.0;
  p2.y = 0.0;
  p2.z = 0.0;

  p3.x = 1.5;
  p3.y = 0.0;
  p3.z = 0.0;

  p4.x = 1.5;
  p4.y = 0.0;
  p4.z = 1.5;

  ff.initialize();
  ff.minimize(10, 1e-8, 1e-8);
  double cosPhi = ForceFields::MMFF::Utils::calcTorsionCosPhi(
      *(RDGeom::Point3D *)ff.positions()[0],
      *(RDGeom::Point3D *)ff.positions()[1],
      *(RDGeom::Point3D *)ff.positions()[2],
      *(RDGeom::Point3D *)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi, 0.5, 1e-4));

  // ------- ------- ------- ------- ------- ------- -------
  // Basic SP2 - SP2
  // ------- ------- ------- ------- ------- ------- -------

  signs[1] = -1;
  v[2] = 0.0;
  v[1] = 7.0;

  ff.contribs().pop_back();
  contrib = new ForceFields::CrystalFF::TorsionAngleContribM6(&ff, 0, 1, 2, 3,
                                                              v, signs);
  ff.contribs().push_back(ForceFields::ContribPtr(contrib));

  p1.x = 0;
  p1.y = 1.5;
  p1.z = 0.1;

  p2.x = 0.0;
  p2.y = 0.0;
  p2.z = 0.0;

  p3.x = 1.5;
  p3.y = 0.0;
  p3.z = 0.0;

  p4.x = 1.5;
  p4.y = 0.2;
  p4.z = 1.5;

  ff.initialize();
  ff.minimize(10, 1e-8, 1e-8);
  cosPhi = ForceFields::MMFF::Utils::calcTorsionCosPhi(
      *(RDGeom::Point3D *)ff.positions()[0],
      *(RDGeom::Point3D *)ff.positions()[1],
      *(RDGeom::Point3D *)ff.positions()[2],
      *(RDGeom::Point3D *)ff.positions()[3]);
  TEST_ASSERT(RDKit::feq(cosPhi, 1.0, 1e-4));
}

void testTorsionPrefs() {
  ROMol *mol;
  mol = SmilesToMol("CCCC");
  TEST_ASSERT(mol);

  ForceFields::CrystalFF::CrystalFFDetails details;
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false, 1,
                                                  false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 1);
  TEST_ASSERT(details.expTorsionAngles.size() == 1);
  TEST_ASSERT(details.expTorsionAtoms[0][0] == 0);
  TEST_ASSERT(details.expTorsionAtoms[0][3] == 3);
  TEST_ASSERT(details.expTorsionAngles[0].first.size() == 6);
  TEST_ASSERT(details.expTorsionAngles[0].second.size() == 6);

  delete mol;
  mol = SmilesToMol("CCCCC");
  TEST_ASSERT(mol);

  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false, 1,
                                                  false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 2);
  TEST_ASSERT(details.expTorsionAngles.size() == 2);
  delete mol;
}

int main() {
  testTorsionAngleM6();
  testTorsionPrefs();

  return 0;
}
