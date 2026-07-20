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
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/PlanarityContribs.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleContribs.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>

#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/TorsionAngle.h>
#include <ForceField/FiniteDifference.h>
#include <ForceField/ForceField.h>

#include <string>
#include <cmath>
#include <iostream> 


using namespace RDGeom;
using namespace RDKit;

void testTorsionAngleContribs() {
  std::cerr << "-------------------------------------" << std::endl;
  std::cerr << " Test CrystalFF torsional term." << std::endl;

  ForceFields::ForceField ff;
  Point3D p1, p2, p3, p4;
  RDGeom::PointPtrVect &ps = ff.positions();
  ps.push_back(&p1);
  ps.push_back(&p2);
  ps.push_back(&p3);
  ps.push_back(&p4);

  // ------- ------- ------- ------- ------- ------- -------
  // Basic SP3 - SP3
  // ------- ------- ------- ------- ------- ------- -------

  // [!#1:1][CX4H2:2]!@[CX4H2:3][!#1:4] 1 0.0 1 0.0 1 4.0 1 0.0 1 0.0 1 0.0
  std::vector<int> signs(6, 1);
  std::vector<double> v(6, 0.0);
  v[2] = 4.0;

  auto contrib = new ForceFields::CrystalFF::TorsionAngleContribs(&ff);
  contrib->addContrib(0, 1, 2, 3, v, signs);
  ff.contribs().emplace_back(contrib);

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
  contrib = new ForceFields::CrystalFF::TorsionAngleContribs(&ff);
  contrib->addContrib(0, 1, 2, 3, v, signs);
  ff.contribs().emplace_back(contrib);

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
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false,
                                                  false, false, 1, false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 1);
  TEST_ASSERT(details.expTorsionAngles.size() == 1);
  TEST_ASSERT(details.expTorsionAtoms[0][0] == 0);
  TEST_ASSERT(details.expTorsionAtoms[0][3] == 3);
  TEST_ASSERT(details.expTorsionAngles[0].first.size() == 6);
  TEST_ASSERT(details.expTorsionAngles[0].second.size() == 6);

  std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                         const ForceFields::CrystalFF::ExpTorsionAngle *>>
      torsionBonds;
  ForceFields::CrystalFF::getExperimentalTorsions(
      *mol, details, torsionBonds, true, false, false, false, 2, false);
  TEST_ASSERT(torsionBonds.size() == 1);
  TEST_ASSERT(std::get<0>(torsionBonds[0]) == 1);
  TEST_ASSERT(std::get<2>(torsionBonds[0])->smarts ==
              "[!#1:1][CX4H2:2]!@;-[CX4H2:3][!#1:4]");
  TEST_ASSERT(std::get<2>(torsionBonds[0])->torsionIdx == 229);

  delete mol;
  mol = SmilesToMol("CCCCC");
  TEST_ASSERT(mol);

  ForceFields::CrystalFF::getExperimentalTorsions(
      *mol, details, torsionBonds, true, false, false, false, 2, false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 2);
  TEST_ASSERT(details.expTorsionAngles.size() == 2);
  TEST_ASSERT(torsionBonds.size() == 2);
  TEST_ASSERT(std::get<0>(torsionBonds[0]) == 1);
  TEST_ASSERT(std::get<0>(torsionBonds[1]) == 2);
  delete mol;
}

void testTorsionPrefsSmallRings() {
  ROMol *mol;
  ForceFields::CrystalFF::CrystalFFDetails details;
  mol = SmilesToMol("C1COCC1");
  TEST_ASSERT(mol);

  // small ring torsion turned off
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false,
                                                  false, false, 1, false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 0);
  TEST_ASSERT(details.expTorsionAngles.size() == 0);

  // small ring torsion turned on
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, true,
                                                  false, false, 1, false);
  TEST_ASSERT(details.expTorsionAtoms.size() > 0);
  TEST_ASSERT(details.expTorsionAngles.size() > 0);
  delete mol;
}

void testTorsionPrefsBridgedSmallRings() {
  // test bridged system does not trigger any small ring torsion
  ROMol *mol;
  ForceFields::CrystalFF::CrystalFFDetails details;

  mol = SmilesToMol("O[C@H]1C[C@H]2CC[C@]1(C)C2(C)C");
  TEST_ASSERT(mol);

  // small ring torsion turned off
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, true,
                                                  false, false, 1, false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 0);
  TEST_ASSERT(details.expTorsionAngles.size() == 0);
  delete mol;
}

void testTorsionPrefsMacrocycles() {
  ROMol *mol;
  ForceFields::CrystalFF::CrystalFFDetails details;
  mol = SmilesToMol("C1COCCCCCCC1");
  TEST_ASSERT(mol);

  // macrocycle ring torsion turned off
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false,
                                                  false, false, 1, false);
  TEST_ASSERT(details.expTorsionAtoms.size() == 0);
  TEST_ASSERT(details.expTorsionAngles.size() == 0);

  // macrocycle ring torsion turned on
  ForceFields::CrystalFF::getExperimentalTorsions(*mol, details, true, false,
                                                  true, false, 1, false);
  TEST_ASSERT(details.expTorsionAtoms.size() > 0);
  TEST_ASSERT(details.expTorsionAngles.size() > 0);
  delete mol;
}

void testPlanarityContribs() {
  std::cerr << "--------------------------------------------" << std::endl;
  std::cerr << "Finite difference Test of Planarity Contribs" << std::endl;
  std::string mb = R"MOL(  // Starting from a DG benzene Conformer
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.0386    0.0113    0.1053 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1602   -0.1902   -0.0318 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7173   -0.8338   -0.0364 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4043    1.0127   -0.0371 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
M  END
)MOL";
  auto ps = v2::FileParsers::MolFileParserParams{.removeHs = false};
  auto mol = v2::FileParsers::MolFromMolBlock(mb, ps);

  TEST_ASSERT(mol);
  auto ff = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
  ff->initialize();
  TEST_ASSERT(ff);

  auto c = std::make_unique<ForceFields::CrystalFF::PlanarityContribs>(ff.get());
  c->addContrib(1,0,2,3,1.0);
  ff->contribs().push_back(std::move(c));

  TEST_ASSERT(ff->calcEnergy() > 1.0);
  const double delta = ForceFields::calcFiniteDifference(*ff);
  TEST_ASSERT(delta < 1e-6);
}

int main() {
  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog)
      << "********************************************************\n";
  BOOST_LOG(rdInfoLog) << "Testing Crystal ForceField\n";

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t SMARTS parsing\n";
  testTorsionAngleContribs();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Seeing if non-ring torsions are applied\n";
  testTorsionPrefs();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Seeing if small ring torsions are applied\n";
  testTorsionPrefsSmallRings();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Seeing if macrocycle ring torsions are applied\n";
  testTorsionPrefsMacrocycles();

  BOOST_LOG(rdInfoLog) << "\t---------------------------------\n";
  BOOST_LOG(rdInfoLog) << "\t Finite difference test of PlanarityContribs\n";
  testPlanarityContribs();

  return 0;
}
