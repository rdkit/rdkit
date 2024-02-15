//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>

using namespace RDKit;
using namespace RDKit::v2;

TEST_CASE("ForwardSDMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Data/NCI/first_200.props.sdf";
    {
      std::ifstream inStream(fName);
      FileParsers::ForwardSDMolSupplier sdsup(&inStream, false);
      auto mol = sdsup.next();
      REQUIRE(mol);
    }
  }
}

TEST_CASE("SDMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Data/NCI/first_200.props.sdf";
    {
      std::ifstream inStream(fName);
      FileParsers::SDMolSupplier sdsup(&inStream, false);
      auto mol = sdsup.next();
      REQUIRE(mol);
      auto mol2 = sdsup[10];
      REQUIRE(mol2);
    }
  }
}
