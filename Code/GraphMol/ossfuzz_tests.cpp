//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <fstream>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit;
std::string pathName = std::string(getenv("RDBASE")) +
                       std::string("/Code/GraphMol/ossfuzz_testdata/");

TEST_CASE("Mol file parsing", "[ossfuzz]") {
  SECTION("Issue 24074") {
    std::string fname = pathName +
                        "clusterfuzz-testcase-minimized-mol_data_stream_to_mol_"
                        "fuzzer-5318835625000960";
    std::ifstream ins(fname, std::ifstream::in | std::ifstream::binary);
    unsigned int line;
    REQUIRE_THROWS_AS(MolDataStreamToMol(&ins, line), Invar::Invariant);
  }
}
