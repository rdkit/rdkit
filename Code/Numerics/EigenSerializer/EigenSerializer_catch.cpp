//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <Numerics/EigenSerializer/EigenSerializer.h>
#include <RDGeneral/Invariant.h>

#include "RDGeneral/test.h"
#include "catch.hpp"

TEST_CASE("Eigen Matrix Serialization Test") {
  SECTION("Test serialize and deserialize") {
    // make random data
    Eigen::ArrayXXd J(10, 5);
    J.setRandom();
    // write the data
    RDNumeric::EigenSerializer::serialize(J, "matrix.bin");
    // read the data
    Eigen::ArrayXXd JCopy;
    RDNumeric::EigenSerializer::deserialize(JCopy, "matrix.bin");
    // check that we read the same as we wrote
    CHECK(J.isApprox(JCopy));
  }
  SECTION("Test serializeAll and deserializeAll") {
    // make random data
    std::vector<Eigen::ArrayXXd> matrices;
    Eigen::ArrayXXd J(10, 5), K(10, 5);
    J.setRandom();
    K.setRandom();
    matrices.push_back(J);
    matrices.push_back(K);
    std::vector<std::string> labels = {"label J", "label K"};
    // write the data
    RDNumeric::EigenSerializer::serializeAll(matrices, labels, "multi.bin");
    // read the data
    std::vector<Eigen::ArrayXXd> read_matrices;
    std::vector<std::string> read_labels;
    RDNumeric::EigenSerializer::deserializeAll(read_matrices, read_labels,
                                            "multi.bin");
    // check that we read the same as we wrote
    for (size_t i = 0; i < matrices.size(); i++) {
      CHECK(read_matrices[i].isApprox(matrices[i]));
      CHECK(read_labels[i] == labels[i]);
    }
  }
}