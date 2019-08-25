//
//  Copyright (C) 2019 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <RDGeneral/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <string>
#include <fstream>

using namespace RDKit;

TEST_CASE("basic options","[molinterchange]"){
    auto m="c1ccccc1"_smiles;
    REQUIRE(m);
    SECTION("basics1") {
        auto jsond = MolInterchange::MolToJSONData(*m);
        CHECK(jsond.find("defaults") != std::string::npos);
        CHECK(jsond.find("extensions") != std::string::npos);
    }
    SECTION("basics2") {
        MolInterchange::JSONWriteParameters ps;
        ps.useDefaults = false;
        auto jsond = MolInterchange::MolToJSONData(*m,ps);
        CHECK(jsond.find("defaults") == std::string::npos);
        CHECK(jsond.find("extensions") != std::string::npos);
    }
    SECTION("basics3") {
        MolInterchange::JSONWriteParameters ps;
        ps.includeExtensions = false;
        auto jsond = MolInterchange::MolToJSONData(*m,ps);
        CHECK(jsond.find("defaults") != std::string::npos);
        CHECK(jsond.find("extensions") == std::string::npos);
    }
    SECTION("validation_json") {
        MolInterchange::JSONWriteParameters ps;
        ps.doValidationJSON = true;
        auto jsond = MolInterchange::MolToJSONData(*m,ps);
        CHECK(jsond.find("\"explicitValence\":3") != std::string::npos);
        CHECK(jsond.find("\"bo\":100") != std::string::npos);
        CHECK(jsond.find("commonchem") == std::string::npos);
        CHECK(jsond.find("validation_JSON") != std::string::npos);
    }
}

