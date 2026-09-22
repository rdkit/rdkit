//
// Copyright (C) 2026 ETH Zurich and other RDKit contributors.
// Author: Niels Maeder
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cmath>
#include <numbers>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>
#include <GraphMol/MolOps.h>

#include "GaussianTorsionAngleContribs.h"
#include "TorsionPreferences.h"

using namespace RDKit;
constexpr double PI = std::numbers::pi;

inline double _numeric_dEdPhi(const std::vector<double> &heights,
                              const std::vector<double> &positions,
                              const std::vector<double> &widths, double phi,
                              double h = 1e-6) {
  const double f_plus =
      ForceFields::CrystalFF::getEnergy(heights, positions, widths, phi + h);
  const double f_minus =
      ForceFields::CrystalFF::getEnergy(heights, positions, widths, phi - h);
  return (f_plus - f_minus) / (2.0 * h);
}

TEST_CASE("GaussianTorsionContribsBasics") {
  auto match = Catch::Matchers::WithinAbs;

  SECTION("1 Peak at 0") {
    const std::vector<double> heights{1.0};
    const std::vector<double> positions{0.0};
    const std::vector<double> widths{1.0};

    SECTION("Energy") {
      const double energy_0 =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, 0.0);
      const double energy_pi =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, PI);
      const double energy_pi2 =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, PI / 2);
      const double energy_3pi2 = ForceFields::CrystalFF::getEnergy(
          heights, positions, widths, PI * 3 / 2);
      const double energy_2pi =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, 2 * PI);

      // Energies should be non-infinite
      CHECK(std::isfinite(energy_0));
      CHECK(std::isfinite(energy_pi));
      CHECK(std::isfinite(energy_pi2));
      CHECK(std::isfinite(energy_3pi2));
      CHECK(std::isfinite(energy_2pi));

      // Energies should be the same at mirrored place
      CHECK(energy_0 == energy_2pi);
      CHECK(energy_3pi2 == energy_pi2);

      // Also check that they are different when they should!
      CHECK(energy_0 != energy_pi2);

      // Check ranking is correct
      CHECK(energy_0 < energy_pi2);
      CHECK(energy_pi2 < energy_pi);
    }

    SECTION("Gradient") {
      const double grad_0 =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, 0.0);
      const double grad_pi =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, PI);
      const double grad_pi2 =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, PI / 2);
      const double grad_3pi2 = ForceFields::CrystalFF::getdEdPhi(
          heights, positions, widths, PI * 3 / 2);
      const double grad_2pi =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, 2 * PI);

      // Gradients should be finite
      CHECK(std::isfinite(grad_0));
      CHECK(std::isfinite(grad_pi));
      CHECK(std::isfinite(grad_pi2));
      CHECK(std::isfinite(grad_3pi2));
      CHECK(std::isfinite(grad_2pi));

      // Gradients should be the same at mirrored place * -1
      CHECK_THAT(grad_0, match(0.0, 1e-20));
      CHECK_THAT(grad_0, match(grad_2pi, 1e-20));

      CHECK(grad_pi2 > 0.0);
      CHECK(grad_3pi2 < 0.0);
      CHECK(grad_3pi2 == -grad_pi2);
    }
    SECTION("FiniteDiff") {
      for (std::size_t i = 0; i < 360; ++i) {
        const double phi = 2 * i / 360 * PI;
        const double grad_pi =
            ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, phi);
        const double grad_pi_num =
            _numeric_dEdPhi(heights, positions, widths, phi);
        CHECK_THAT(grad_pi, match(grad_pi_num, 1e-10));
      }
    }
  }

  SECTION("3 Peaks") {
    const std::vector<double> heights{1.0, 0.3, .5};
    const std::vector<double> positions{0.0, 1.0, PI};
    const std::vector<double> widths{1.0, .1, .25};

    SECTION("Energy") {
      const double energy_0 =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, 0.0);
      const double energy_pi =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, PI);
      const double energy_pi2 =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, PI / 2);
      const double energy_3pi2 = ForceFields::CrystalFF::getEnergy(
          heights, positions, widths, PI * 3 / 2);
      const double energy_2pi =
          ForceFields::CrystalFF::getEnergy(heights, positions, widths, 2 * PI);

      // Energies should be non-infinite
      CHECK(std::isfinite(energy_0));
      CHECK(std::isfinite(energy_pi));
      CHECK(std::isfinite(energy_pi2));
      CHECK(std::isfinite(energy_3pi2));
      CHECK(std::isfinite(energy_2pi));

      // Energies should be the same at mirrored place
      CHECK(energy_0 == energy_2pi);
      CHECK(energy_3pi2 == energy_pi2);

      // Check ranking is correct pi/2 has highest energy in this potential
      CHECK(energy_0 < energy_pi2);
      CHECK(energy_pi2 > energy_pi);
    }

    SECTION("Gradient") {
      const double grad_0 =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, 0.0);
      const double grad_pi =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, PI);
      const double grad_1 =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, PI);
      const double grad_pi2 =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, PI / 2);
      const double grad_3pi2 = ForceFields::CrystalFF::getdEdPhi(
          heights, positions, widths, PI * 3 / 2);
      const double grad_2pi =
          ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, 2 * PI);

      // Gradients should be finite
      CHECK(std::isfinite(grad_0));
      CHECK(std::isfinite(grad_pi));
      CHECK(std::isfinite(grad_pi2));
      CHECK(std::isfinite(grad_3pi2));
      CHECK(std::isfinite(grad_2pi));

      // Gradients should be the same at mirrored place * -1
      CHECK_THAT(grad_0, match(0.0, 1e-20));
      CHECK_THAT(grad_1, match(0.0, 1e-20));
      CHECK_THAT(grad_pi, match(0.0, 1e-20));

      CHECK_THAT(grad_0, match(grad_2pi, 1e-20));

      CHECK(grad_pi2 > 0.0);
      CHECK(grad_3pi2 < 0.0);
      CHECK(grad_3pi2 == -grad_pi2);
    }
    SECTION("FiniteDiff") {
      for (std::size_t i = 0; i < 360; ++i) {
        const double phi = 2 * i / 360 * PI;
        const double grad_pi =
            ForceFields::CrystalFF::getdEdPhi(heights, positions, widths, phi);
        const double grad_pi_num =
            _numeric_dEdPhi(heights, positions, widths, phi);
        CHECK_THAT(grad_pi, match(grad_pi_num, 1e-10));
      }
    }
  }
}

namespace {

constexpr std::size_t LOOKUP_SIZE = 180;
using GaussianDetails = ForceFields::CrystalFF::CrystalFFDetails<
    ForceFields::CrystalFF::GaussianExp_T>;

}  // namespace

TEST_CASE("GaussianTorsionContribsLookupTable") {
  SECTION("Basic") {
    const std::vector<double> heights{1.0};
    const std::vector<double> positions{0.0};
    const std::vector<double> widths{1.0};

    GaussianDetails details;
    // arbitrary index, should afterwards be mapped to 0
    details.torsionIdx = {123};
    details.expTorsionAngles.reserve(1);
    details.expTorsionAngles.emplace_back(heights, positions, widths, 1.0);

    ForceFields::CrystalFF::populateRefTable(details);

    REQUIRE(details.phiToEnergy.size() == 1);
    REQUIRE(details.phiToGrad.size() == 1);

    REQUIRE(details.phiToEnergy[0].size() == LOOKUP_SIZE);
    REQUIRE(details.phiToGrad[0].size() == LOOKUP_SIZE);

    auto &energies = details.phiToEnergy[0];
    auto &gradients = details.phiToGrad[0];

    CHECK(energies[0] < energies[179]);

    CHECK_THAT(gradients[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  }

  const bool useExpTorsions = true;
  const bool useSmallRingTorsions = false;
  const bool useMacrocycleTorsions = false;
  const bool useBasicKnowledge = true;
  const unsigned int version = 4;
  std::vector<
      std::tuple<unsigned int, std::vector<unsigned int>,
                 const ForceFields::CrystalFF::GaussianExpTorsionAngle *>>
      torsionBonds;
  SECTION("simple molecule") {
    auto mol = "CCCC"_smiles;
    MolOps::addHs(*mol);
    GaussianDetails details;

    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol, details, torsionBonds, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version);

    ForceFields::CrystalFF::populateRefTable(details);

    REQUIRE(details.phiToEnergy.size() == 1);
    REQUIRE(details.phiToGrad.size() == 1);

    REQUIRE(details.phiToEnergy[0].size() == LOOKUP_SIZE);
    REQUIRE(details.phiToGrad[0].size() == LOOKUP_SIZE);

    auto &energies = details.phiToEnergy[0];
    auto &gradients = details.phiToGrad[0];

    // Pattern for CCCC has a peak at 180 degrees
    CHECK(energies[0] > energies[179]);
    CHECK_THAT(gradients[179], Catch::Matchers::WithinAbs(0.0, 1e-10));
  }
  SECTION("molecule w/ aromat") {
    auto mol = "c1ccccc1CC"_smiles;
    MolOps::addHs(*mol);
    GaussianDetails details;

    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol, details, torsionBonds, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version);

    ForceFields::CrystalFF::populateRefTable(details);

    // aromat and 1 exp torsion
    REQUIRE(details.phiToEnergy.size() == 2);
    REQUIRE(details.phiToGrad.size() == 2);

    // aromatic should come last
    auto &energies = details.phiToEnergy[1];
    auto &gradients = details.phiToGrad[1];

    // minimum at 0!
    for (std::size_t i = 1; i < 180; ++i) {
      CHECK(energies[0] < energies[i]);
    }
    CHECK_THAT(gradients[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  }
  SECTION("molecule w/ amide") {
    auto mol = "CCC(=O)NC"_smiles;
    MolOps::addHs(*mol);
    GaussianDetails details;

    // needed to set cis / trans in amides...
    DistGeom::BoundsMatPtr mmat;
    mmat.reset(new DistGeom::BoundsMatrix(mol->getNumAtoms()));
    DGeomHelpers::initBoundsMat(mmat);
    auto params = DGeomHelpers::KDG;
    DGeomHelpers::setTopolBounds(*mol, mmat, params, false, true, true, true,
                                 &details.path14Configs);

    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol, details, torsionBonds, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version);

    ForceFields::CrystalFF::populateRefTable(details);

    // 2 for the amide and 2 exp torsion
    REQUIRE(details.phiToEnergy.size() == 4);
    REQUIRE(details.phiToGrad.size() == 4);

    // after reordering trans should be second last
    auto &energies_t = details.phiToEnergy[2];
    auto &gradients_t = details.phiToGrad[2];

    // and_cis last
    auto &energies_c = details.phiToEnergy[3];
    auto &gradients_c = details.phiToGrad[3];

    // minimum at 180!
    for (std::size_t i = 0; i < 179; ++i) {
      CHECK(energies_t[179] < energies_t[i]);
    }
    for (std::size_t i = 1; i < 180; ++i) {
      CHECK(energies_c[0] < energies_c[i]);
    }
    CHECK_THAT(gradients_t[179], Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(gradients_c[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  }

  SECTION("molecule w/ amide and aromat") {
    auto mol = "c1ccccc1C(=O)NC"_smiles;
    MolOps::addHs(*mol);
    GaussianDetails details;

    // needed to set cis / trans in amides...
    DistGeom::BoundsMatPtr mmat;
    mmat.reset(new DistGeom::BoundsMatrix(mol->getNumAtoms()));
    DGeomHelpers::initBoundsMat(mmat);
    auto params = DGeomHelpers::KDG;
    DGeomHelpers::setTopolBounds(*mol, mmat, params, false, true, true, true,
                                 &details.path14Configs);

    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol, details, torsionBonds, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version);

    ForceFields::CrystalFF::populateRefTable(details);

    // 2 for the amide and 2 exp torsion 1 aromat
    REQUIRE(details.phiToEnergy.size() == 5);
    REQUIRE(details.phiToGrad.size() == 5);

    // after reordering trans should be second last
    auto &energies_t = details.phiToEnergy[2];
    auto &gradients_t = details.phiToGrad[2];

    // and cis second to last
    auto &energies_c = details.phiToEnergy[3];
    auto &gradients_c = details.phiToGrad[3];

    // and aromat last
    auto &energies_a = details.phiToEnergy[4];
    auto &gradients_a = details.phiToGrad[4];


    // minimum at 180!
    for (std::size_t i = 0; i < 179; ++i) {
      CHECK(energies_t[179] < energies_t[i]);
    }
    for (std::size_t i = 1; i < 180; ++i) {
      CHECK(energies_c[0] < energies_c[i]);
    }
    for (std::size_t i = 1; i < 180; ++i) {
      CHECK(energies_a[0] < energies_c[i]);
    }
    CHECK_THAT(gradients_t[179], Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(gradients_c[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK_THAT(gradients_a[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  }
}