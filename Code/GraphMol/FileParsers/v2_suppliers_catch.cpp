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
#include <boost/format.hpp>

using namespace RDKit;
using namespace RDKit::v2;

TEST_CASE("ForwardSDMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Data/NCI/first_200.props.sdf";
    {
      std::ifstream inStream(fName, std::ios::binary);
      FileParsers::ForwardSDMolSupplier sdsup(&inStream, false);
      auto mol = sdsup.next();
      REQUIRE(mol);
    }
    {  // v1
      std::ifstream inStream(fName, std::ios::binary);
      ForwardSDMolSupplier sdsup(&inStream, false);
      std::unique_ptr<ROMol> mol{sdsup.next()};
      REQUIRE(mol);
    }
  }
}

TEST_CASE("SDMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Data/NCI/first_200.props.sdf";
    {
      std::ifstream inStream(fName, std::ios::binary);
      FileParsers::SDMolSupplier sdsup(&inStream, false);
      auto mol = sdsup.next();
      REQUIRE(mol);
      auto mol2 = sdsup[10];
      REQUIRE(mol2);
    }
    {  // v1
      std::ifstream inStream(fName, std::ios::binary);
      SDMolSupplier sdsup(&inStream, false);
      std::unique_ptr<ROMol> mol{sdsup.next()};
      REQUIRE(mol);
    }
  }
}

TEST_CASE("SmilesMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Data/NCI/first_200.tpsa.csv";
    {
      std::ifstream inStream(fName, std::ios::binary);
      FileParsers::SmilesMolSupplierParams params;
      params.delimiter = ',';
      FileParsers::SmilesMolSupplier smsup(&inStream, false, params);
      auto mol = smsup.next();
      REQUIRE(mol);
      auto mol2 = smsup[10];
      REQUIRE(mol2);
    }
    {  // v1
      std::ifstream inStream(fName, std::ios::binary);
      SmilesMolSupplier smsup(&inStream, false, ",");
      std::unique_ptr<ROMol> mol{smsup.next()};
      REQUIRE(mol);
    }
  }
}

TEST_CASE("TDTMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/FileParsers/test_data/acd_few.tdt";
    {
      std::ifstream inStream(fName, std::ios::binary);
      FileParsers::TDTMolSupplier smsup(&inStream, false);
      auto mol = smsup.next();
      REQUIRE(mol);
      CHECK(smsup.length() == 10);
      auto mol2 = smsup[8];
      REQUIRE(mol2);
    }
    {  // v1
      std::ifstream inStream(fName, std::ios::binary);
      TDTMolSupplier smsup(&inStream, false, ",");
      std::unique_ptr<ROMol> mol{smsup.next()};
      REQUIRE(mol);
    }
  }
}

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
TEST_CASE("MaeMolSupplier") {
  SECTION("basics") {
    std::string fName = getenv("RDBASE");
    fName += "/Code/GraphMol/FileParsers/test_data/props_test.mae";

    FileParsers::MaeMolSupplier maesup(fName);
    auto nmol(maesup.next());
    TEST_ASSERT(nmol);
  }
}

TEST_CASE("TestParsingInvalidChiralityLabelsWithMaeMolSupplier") {
  constexpr const char *maeblock_template = R"DATA(f_m_ct {
  i_m_ct_stereo_status
  s_m_title
  s_st_Chirality_1
  :::
  1
  ""
  %s
  m_atom[5] {
    # First column is Index #
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord
    i_m_atomic_number
    i_m_formal_charge
    :::
    1 -1.000000 -0.060606 0.000000 6 0
    2 -1.000000 -1.560606 0.000000 6 0
    3 -0.250000 -2.859644 0.000000 17 0
    4 -2.500000 -1.560606 0.000000 9 0
    5 0.299038 -0.810606 0.000000 35 0
    :::
  }
  m_bond[4] {
    # First column is Index #
    i_m_from
    i_m_order
    i_m_to
    :::
    1 1 1 2
    2 2 1 3
    3 2 1 4
    4 2 1 5
    :::
  }
}
)DATA";
  auto invalid_chirality_label =
      GENERATE("12_R_1_3_4_5",     // missing chiral atom
               "12_ANR_1_3_4_15",  // missing substituent atom
               "2_S_1_3_4",        // incomplete substituent list
               "2_ANS_1_3_4_5_2"   // self bond and too many atoms
      );
  CAPTURE(invalid_chirality_label);

  auto maeblock =
      (boost::format(maeblock_template) % invalid_chirality_label).str();

  auto iss = std::make_unique<std::istringstream>(maeblock);
  constexpr bool takeOwnership = true;
  constexpr FileParsers::MaeMolSupplierParams params{.sanitize = false};

  FileParsers::MaeMolSupplier suppl(iss.release(), takeOwnership, params);
  auto mol = suppl.next();

  REQUIRE(mol);
  CHECK(mol->getNumAtoms() == 5);
  CHECK(mol->getNumBonds() == 4);
}

#endif
