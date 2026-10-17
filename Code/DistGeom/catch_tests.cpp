//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>
#include "ZMatrix.h"
#include "ZMatrixUtils.h"

#ifdef RDK_TEST_MULTITHREADED
#include <csignal>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

TEST_CASE("ZMatrix baiscs") {
  DistGeom::ZMatrix zmat(5);
  zmat.addElement(4);
  zmat.addElement(1, {1.0}, {4});
  zmat.addElement(3, {1.2}, {1}, {M_PI}, {4});
  zmat.addElement(0, {0.2}, {3}, {2 * M_PI}, {1},
                  {DistGeom::TorsionRange{0.0, 2.0 * M_PI}}, {4});
  zmat.addElement(2, {0.2}, {3}, {2 * M_PI}, {1}, {}, {}, {{M_PI / 2.0, 0}});

  SECTION("Order") {
    CHECK(zmat[0].atomIdx == 4u);
    CHECK(zmat[4].torsionDependence.has_value());
  }
  SECTION("Invert torsion") {
    zmat.invertTorsionDependence(2);
    CHECK_THAT(zmat[4].torsionDependence->offset,
               Catch::Matchers::WithinAbs(-M_PI / 2.0, 0e-4));
  }
}

TEST_CASE("Torsion Candidates") {
  DistGeom::TorsionCandidates cand1 = DistGeom::TorsionRange(1.0, 2.0);
  DistGeom::TorsionCandidates cand2 = DistGeom::TorsionRange(1.9, 2.1);
  DistGeom::TorsionCandidates cand3 = DistGeom::TorsionValues({0.0, 1.3});
  DistGeom::TorsionCandidates cand4 = DistGeom::TorsionValues({-1.3, 0.0, 4.0});

  SECTION("Less") {
    CHECK(DistGeom::less(cand2, cand1));
    CHECK(DistGeom::less(cand3, cand2));
    CHECK(DistGeom::less(cand3, cand4));
    CHECK_FALSE(DistGeom::less(cand1, cand3));
  }

  SECTION("Eqaul") {
    CHECK(DistGeom::equal(DistGeom::TorsionValues{-M_PI},
                          DistGeom::TorsionValues{M_PI}));
    CHECK_FALSE(DistGeom::equal(DistGeom::TorsionValues{0.0},
                                DistGeom::TorsionValues{M_PI}));
  }

  SECTION("Contains") {
    DistGeom::TorsionCandidates fullRange =
        DistGeom::TorsionRange(0.0, 2.0 * M_PI);
    CHECK(DistGeom::contains(fullRange, 0.0));
    CHECK(DistGeom::contains(fullRange, 2.0 * M_PI));
    CHECK(DistGeom::contains(fullRange, 1.0));
    DistGeom::TorsionCandidates range2 = DistGeom::TorsionRange(2.0, 4.0);
    CHECK(DistGeom::contains(range2, 2.0));
    CHECK(DistGeom::contains(range2, 4.0));
    CHECK(DistGeom::contains(range2, M_PI));
    CHECK_FALSE(DistGeom::contains(range2, 0.0));
    CHECK_FALSE(DistGeom::contains(range2, 5.0));
    DistGeom::TorsionCandidates rangeAcross = DistGeom::TorsionRange(-4.0, 1.0);
    CHECK(DistGeom::contains(rangeAcross, 1.0));
    CHECK(DistGeom::contains(rangeAcross, -4.0));
    CHECK(DistGeom::contains(rangeAcross, 0));
    CHECK_FALSE(DistGeom::contains(rangeAcross, 2.0));

    CHECK(DistGeom::contains(cand3, 0.0));
    CHECK(DistGeom::contains(cand4, -1.3));
    CHECK(DistGeom::contains(cand4, 4.0));

    DistGeom::TorsionCandidates cand5 = DistGeom::TorsionValues({2.0 * M_PI});
    CHECK(DistGeom::contains(cand5, 2.0 * M_PI));
  }

  SECTION("Merge") {
    auto merge12 = DistGeom::merge(cand1, cand2);
    CHECK_NOTHROW(std::get<DistGeom::TorsionRange>(merge12));
    auto merge12range = std::get<DistGeom::TorsionRange>(merge12);
    CHECK_THAT(merge12range.lower, Catch::Matchers::WithinAbs(1.9, 1e-6));
    CHECK_THAT(merge12range.upper, Catch::Matchers::WithinAbs(2.0, 1e-6));

    auto merge13 = DistGeom::merge(cand1, cand3);
    CHECK_NOTHROW(std::get<DistGeom::TorsionValues>(merge13));
    auto merge13vals = std::get<DistGeom::TorsionValues>(merge13);
    CHECK(merge13vals.size() == 1);
    CHECK_THAT(merge13vals.firstValue(), Catch::Matchers::WithinAbs(1.3, 1e-6));

    auto merge43 = DistGeom::merge(cand4, cand3);
    CHECK_NOTHROW(std::get<DistGeom::TorsionValues>(merge43));
    auto merge43vals = std::get<DistGeom::TorsionValues>(merge43);

    CHECK(merge43vals.size() == 1);
    CHECK_THAT(merge43vals.firstValue(), Catch::Matchers::WithinAbs(0.0, 1e-6));
  }

  SECTION("Sample") {
    int seed = 42;
    RDKit::getRandomGenerator(seed);
    auto rng = RDKit::getDoubleRandomSource();
    auto s1 = DistGeom::sample(cand1, rng);
    CHECK(s1 >= 1.0 - 1e-6);
    CHECK(s1 <= 2.0 + 1e-6);
  }
}