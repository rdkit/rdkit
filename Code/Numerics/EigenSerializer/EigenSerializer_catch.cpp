//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"

#include <RDGeneral/Invariant.h>
#include <Numerics/EigenSerializer/EigenSerializer.h>

TEST_CASE("Eigen Matrix Serialization Test") {
  SECTION("Test 1") {
    Eigen::ArrayXXd J(10, 5);
    J.setRandom();
    RDNumeric::EigenSerializer::serialize(J, "matrix.bin");

    Eigen::ArrayXXd JCopy;
    RDNumeric::EigenSerializer::deSerialize(JCopy, "matrix.bin");

    CHECK(J.isApprox(JCopy));
  }
}