//
//  Copyright (C) 2019 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDGeneral/test.h"
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("Property list conversion", "[atom_list_properties]") {
  SECTION("basics: iprops") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->setProp("atom.iprop.foo1", "1   6 9");
    m->setProp("atom.iprop.foo2", "3 n/a   9");
    m->setProp("atom.iprop.foo3", "[?]  5 1 ?");
    m->setProp("atom.iprop.foo4", "[foo] 3 foo   9");
    FileParserUtils::applyMolListPropsToAtoms<std::int64_t>(*m, "atom.iprop.");
    CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
    CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
    CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo4"));
    CHECK(!m->getAtomWithIdx(1)->hasProp("foo4"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo4"));
    CHECK(m->getAtomWithIdx(1)->getProp<std::int64_t>("foo1") == 6);
    CHECK(m->getAtomWithIdx(2)->getProp<std::int64_t>("foo2") == 9);
    CHECK(m->getAtomWithIdx(1)->getProp<std::int64_t>("foo3") == 1);
  }
  SECTION("basics: dprops") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->setProp("atom.dprop.foo1", "1   6 9");
    m->setProp("atom.dprop.foo2", "3 n/a   9");
    m->setProp("atom.dprop.foo3", "[?]  5 1 ?");
    FileParserUtils::applyMolListPropsToAtoms<double>(*m, "atom.dprop.");
    CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
    CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
    CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->getProp<double>("foo1") == 6);
    CHECK(m->getAtomWithIdx(2)->getProp<double>("foo2") == 9);
    CHECK(m->getAtomWithIdx(1)->getProp<double>("foo3") == 1);
  }
  SECTION("basics: props") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->setProp("atom.prop.foo1", "1   6 9");
    m->setProp("atom.prop.foo2", "3 n/a   9");
    m->setProp("atom.prop.foo3", "[?]  5 1 ?");
    FileParserUtils::applyMolListPropsToAtoms<std::string>(*m, "atom.prop.");
    CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
    CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
    CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo1") == "6");
    CHECK(m->getAtomWithIdx(2)->getProp<std::string>("foo2") == "9");
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo3") == "1");
  }
  SECTION("basics: bprops") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->setProp("atom.bprop.foo1", "1   0 0");
    m->setProp("atom.bprop.foo2", "0 n/a   1");
    m->setProp("atom.bprop.foo3", "[?]  0 1 ?");
    FileParserUtils::applyMolListPropsToAtoms<bool>(*m, "atom.bprop.");
    CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
    CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
    CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->getProp<bool>("foo1") == false);
    CHECK(m->getAtomWithIdx(2)->getProp<bool>("foo2") == true);
    CHECK(m->getAtomWithIdx(1)->getProp<bool>("foo3") == true);
  }
}
TEST_CASE("processMolPropertyLists", "[atom_list_properties]") {
  SECTION("basics") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->setProp("atom.iprop.foo1", "1   6 9");
    m->setProp("atom.dprop.foo2", "3 n/a   9");
    m->setProp("atom.prop.foo3", "[?]  5 1 ?");
    m->setProp("atom.bprop.foo4", "1 0 0");
    FileParserUtils::processMolPropertyLists(*m);
    CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
    CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
    CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
    CHECK(m->getAtomWithIdx(0)->hasProp("foo4"));
    CHECK(m->getAtomWithIdx(1)->hasProp("foo4"));
    CHECK(m->getAtomWithIdx(2)->hasProp("foo4"));
    CHECK(m->getAtomWithIdx(1)->getProp<std::int64_t>("foo1") == 6);
    CHECK(m->getAtomWithIdx(2)->getProp<double>("foo2") == 9);
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo3") == "1");
    CHECK(m->getAtomWithIdx(1)->getProp<bool>("foo4") == false);
  }
}

TEST_CASE("basic SDF handling", "[SDF][atom_list_properties]") {
  std::string sdf = R"SDF(
     RDKit  2D

  3  3  0  0  0  0  0  0  0  0999 V2000
    0.8660    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4330    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4330   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  1  1  0
M  END
>  <atom.dprop.PartialCharge>  (1) 
0.008 -0.314 0.008

>  <atom.iprop.NumHeavyNeighbors>  (1) 
2 2 2

>  <atom.prop.AtomLabel>  (1)
C1 N2 C3

>  <atom.bprop.IsCarbon>  (1) 
1 0 1

>  <atom.prop.PartiallyMissing>  (1) 
one n/a three

>  <atom.iprop.PartiallyMissingInt>  (1) 
[?] 2 2 ?

$$$$
)SDF";
  SECTION("no processing") {
    RDKit::SDMolSupplier suppl;
    suppl.setData(sdf);
    suppl.setProcessPropertyLists(false);
    std::unique_ptr<RDKit::ROMol> m(suppl[0]);
    REQUIRE(m);
    CHECK(m->hasProp("atom.prop.AtomLabel"));
    CHECK(!m->getAtomWithIdx(0)->hasProp("AtomLabel"));
  }
  SECTION("with processing") {
    RDKit::SDMolSupplier suppl;
    suppl.setData(sdf);
    std::unique_ptr<RDKit::ROMol> m(suppl[0]);
    REQUIRE(m);
    CHECK(m->hasProp("atom.prop.AtomLabel"));
    CHECK(m->getAtomWithIdx(0)->hasProp("AtomLabel"));
  }
}

TEST_CASE("createAtomPropertyLists", "[atom_list_properties]") {
  SECTION("basics") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp<std::int64_t>("foo1", 1);
    m->getAtomWithIdx(2)->setProp<std::int64_t>("foo1", 9);
    FileParserUtils::createAtomIntPropertyList(*m, "foo1");
    REQUIRE(m->hasProp("atom.iprop.foo1"));
    CHECK(m->getProp<std::string>("atom.iprop.foo1") == "1 n/a 9");

    m->getAtomWithIdx(0)->setProp<double>("foo2", 1);
    m->getAtomWithIdx(1)->setProp<double>("foo2", 4);
    m->getAtomWithIdx(2)->setProp<double>("foo2", 9);
    FileParserUtils::createAtomDoublePropertyList(*m, "foo2");
    REQUIRE(m->hasProp("atom.dprop.foo2"));
    CHECK(m->getProp<std::string>("atom.dprop.foo2") == "1 4 9");

    m->getAtomWithIdx(0)->setProp<std::string>("foo3", "1");
    m->getAtomWithIdx(1)->setProp<std::string>("foo3", "4");
    FileParserUtils::createAtomStringPropertyList(*m, "foo3", "?");
    REQUIRE(m->hasProp("atom.prop.foo3"));
    CHECK(m->getProp<std::string>("atom.prop.foo3") == "[?] 1 4 ?");

    m->getAtomWithIdx(0)->setProp<bool>("foo4", 1);
    m->getAtomWithIdx(1)->setProp<bool>("foo4", 0);
    m->getAtomWithIdx(2)->setProp<bool>("foo4", 0);
    FileParserUtils::createAtomBoolPropertyList(*m, "foo4");
    REQUIRE(m->hasProp("atom.bprop.foo4"));
    CHECK(m->getProp<std::string>("atom.bprop.foo4") == "1 0 0");
  }
  SECTION("long lines") {
    auto m = "COC"_smiles;
    REQUIRE(m);
    m->getAtomWithIdx(0)->setProp<std::string>("foo1", std::string(80, 'a'));
    m->getAtomWithIdx(1)->setProp<std::string>("foo1", std::string(80, 'b'));
    m->getAtomWithIdx(2)->setProp<std::string>("foo1", std::string(80, 'c'));
    FileParserUtils::createAtomStringPropertyList(*m, "foo1");
    REQUIRE(m->hasProp("atom.prop.foo1"));
    std::string ps = m->getProp<std::string>("atom.prop.foo1");
    CHECK(ps.length() > 240);
    CHECK(ps.find("\n") != std::string::npos);
    for (auto &atom : m->atoms()) {
      atom->clearProp("foo1");
    }
    FileParserUtils::applyMolListPropsToAtoms<std::string>(*m, "atom.prop.");
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>("foo1") ==
          std::string(80, 'a'));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo1") ==
          std::string(80, 'b'));
    CHECK(m->getAtomWithIdx(2)->getProp<std::string>("foo1") ==
          std::string(80, 'c'));
  }
}