
//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cmath>
#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ForceFieldHelpers/FFConvenience.h>

#include <ForceField/AngleConstraints.h>
#include <ForceField/DistanceConstraints.h>
#include <ForceField/UFF/Params.h>

using namespace RDKit;
namespace {

double get_dist(const ROMol &mol, unsigned int idx1, unsigned int idx2) {
  return (mol.getConformer().getAtomPos(idx1) -
          mol.getConformer().getAtomPos(idx2))
      .length();
}

double get_angle(const ROMol &mol, unsigned int idx1, unsigned int idx2,
                 unsigned int idx3) {
  const auto &pos1 = mol.getConformer().getAtomPos(idx1);
  const auto &pos2 = mol.getConformer().getAtomPos(idx2);
  const auto &pos3 = mol.getConformer().getAtomPos(idx3);
  const RDGeom::Point3D r[2] = {pos1 - pos2, pos3 - pos2};
  const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
                               std::max(1.0e-5, r[1].lengthSq())};
  double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
  cosTheta = std::clamp(cosTheta, -1.0, 1.0);
  return ForceFields::UFF::RAD2DEG * acos(cosTheta);
}
}  // namespace

TEST_CASE("Test DistanceConstraintContribs") {
  auto mol =
      "CCCO |"
      "(-1.28533,-0.0567758,0.434662;-0.175447,0.695786,-0.299881;"
      "0.918409,-0.342619,-0.555572;1.30936,-0.801512,0.71705)|"_smiles;
  REQUIRE(mol);
  SECTION("absolute distance minimization") {
    auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
    REQUIRE(forceField);
    forceField->initialize();
    auto contribs = std::make_unique<ForceFields::DistanceConstraintContribs>(
        forceField.get());
    contribs->addContrib(0, 1, 3.0, 3.0, 1);
    contribs->addContrib(0, 2, 4.0, 4.0, 1);
    forceField->contribs().push_back(std::move(contribs));
    CHECK(forceField->minimize() == 0);
    CHECK(feq(get_dist(*mol, 0, 1), 3.0));
    CHECK(feq(get_dist(*mol, 0, 2), 4.0));
  }
  SECTION("relative distance minimization") {
    auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
    REQUIRE(forceField);
    forceField->initialize();
    auto contribs = std::make_unique<ForceFields::DistanceConstraintContribs>(
        forceField.get());
    contribs->addContrib(0, 1, true, 1.0, 1.0, 1);
    contribs->addContrib(0, 2, false, 4.0, 4.0, 1);
    forceField->contribs().push_back(std::move(contribs));
    auto before = get_dist(*mol, 0, 1);
    CHECK(forceField->minimize() == 0);
    CHECK(feq(get_dist(*mol, 0, 1), 1.0 + before));
    CHECK(feq(get_dist(*mol, 0, 2), 4.0));
  }
}

TEST_CASE("Test AngleConstraintContribs") {
  auto mol =
      "CCCO |"
      "(-1.28533,-0.0567758,0.434662;-0.175447,0.695786,-0.299881;"
      "0.918409,-0.342619,-0.555572;1.30936,-0.801512,0.71705)|"_smiles;
  REQUIRE(mol);
  SECTION("absolute angle minimization") {
    auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
    REQUIRE(forceField);
    forceField->initialize();
    auto contribs = std::make_unique<ForceFields::AngleConstraintContribs>(
        forceField.get());
    contribs->addContrib(0, 1, 2, 120.0, 120.0, 1);
    contribs->addContrib(1, 2, 3, 160.0, 160.0, 1);
    forceField->contribs().push_back(std::move(contribs));
    CHECK(forceField->minimize() == 0);
    CHECK(feq(get_angle(*mol, 0, 1, 2), 120.0));
    CHECK(feq(get_angle(*mol, 1, 2, 3), 160.0));
  }
  SECTION("relative angle minimization") {
    auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
    REQUIRE(forceField);
    forceField->initialize();
    auto contribs = std::make_unique<ForceFields::AngleConstraintContribs>(
        forceField.get());
    contribs->addContrib(0, 1, 2, true, 5.0, 5.0, 1);
    contribs->addContrib(1, 2, 3, false, 160.0, 160.0, 1);
    forceField->contribs().push_back(std::move(contribs));
    auto before = get_angle(*mol, 0, 1, 2);
    CHECK(forceField->minimize() == 0);
    CHECK(feq(get_angle(*mol, 0, 1, 2), before + 5.0));
    CHECK(feq(get_angle(*mol, 1, 2, 3), 160.0));
  }
}

TEST_CASE("Force field copies preserve fixed points") {
  auto mol = "CCC |(0,0,0;2,0,0;4,0,0)|"_smiles;
  REQUIRE(mol);
  auto forceField = ForceFieldsHelper::createEmptyForceFieldForMol(*mol);
  forceField->initialize();
  auto contribs = std::make_unique<ForceFields::DistanceConstraintContribs>(
      forceField.get());
  contribs->addContrib(0, 1, 1.0, 1.0, 10.0);
  forceField->contribs().push_back(std::move(contribs));
  forceField->fixedPoints() = {0, 2};

  ForceFields::ForceField copied(*forceField);
  CHECK(copied.dimension() == forceField->dimension());
  CHECK(copied.numPoints() == forceField->numPoints());
  CHECK(copied.positions().empty());
  CHECK(copied.fixedPoints() == forceField->fixedPoints());
  REQUIRE(copied.contribs().size() == forceField->contribs().size());
  CHECK(copied.contribs()[0].get() != forceField->contribs()[0].get());

  // The fixed-point list belongs to the copy, while positions are rebound to
  // a separate conformer before initialization.
  forceField->fixedPoints().clear();
  CHECK(copied.fixedPoints() == ForceFields::INT_VECT{0, 2});
  ROMol other(*mol);
  for (auto &position : other.getConformer().getPositions()) {
    copied.positions().push_back(&position);
  }
  copied.initialize();
  CHECK(copied.calcEnergy() == forceField->calcEnergy());
  CHECK(copied.minimize() == 0);
  CHECK(other.getConformer().getAtomPos(0).x == 0.0);
  CHECK(other.getConformer().getAtomPos(2).x == 4.0);
  CHECK_THAT(other.getConformer().getAtomPos(1).x,
             Catch::Matchers::WithinAbs(1.0, 1e-3));
  CHECK(mol->getConformer().getAtomPos(1).x == 2.0);
}
