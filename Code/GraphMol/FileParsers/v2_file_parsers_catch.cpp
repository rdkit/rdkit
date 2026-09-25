//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>

#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

using namespace RDKit::v2::FileParsers;

TEST_CASE("expandAttachmentPoints") {
  SECTION("basics") {
    std::string mb = R"CTAB(
  Mrv2211 02082417032D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.7083 3.25 0 0
M  V30 2 O -3.3747 4.02 0 0 ATTCHPT=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB";
    auto mol{MolFromMolBlock(mb)};
    REQUIRE(mol);
    CHECK(mol->getNumAtoms() == 2);
    MolFileParserParams ps;
    ps.expandAttachmentPoints = true;
    auto mol2{MolFromMolBlock(mb, ps)};
    REQUIRE(mol2);
    CHECK(mol2->getNumAtoms() == 3);
    SmilesWriteParams sps;
    CHECK(MolToCXSmiles(*mol2, sps,
                        SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS) ==
          "*OC |$_AP1;;$|");
    // verify coords too
    CHECK(MolToCXSmiles(*mol2) ==
          "*OC |(-2.50869,4.52002,;-3.3747,4.02,;-4.7083,3.25,),$_AP1;;$|");
  }
}

TEST_CASE("empty column names in SmilesMolSupplier") {
  std::string rdbase = getenv("RDBASE");

  std::string fName =
      rdbase + "/Code/GraphMol/FileParsers/test_data/s1p_chembldoc89753.txt";
  SmilesMolSupplierParams params;
  params.delimiter = ",";
  params.smilesColumn = 9;
  params.nameColumn = 10;
  v2::FileParsers::SmilesMolSupplier suppl(fName, params);
  auto mol = suppl.next();
  REQUIRE(mol);
  CHECK(mol->hasProp("_Name"));
  CHECK(mol->hasProp("Column_0"));
}