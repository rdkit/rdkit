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
#include "RDGeneral/test.h"
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("Property list conversion", "[atom_list_properties]") {
  SECTION("basics: iprops") {
      auto m = "COC"_smiles;
      REQUIRE(m);
      m->setProp("atom.iprop.foo1","1   6 9");
      m->setProp("atom.iprop.foo2","3 n/a   9");
      m->setProp("atom.iprop.foo3","[?]  5 1 ?");
      FileParserUtils::applyMolListPropsToAtoms<std::int64_t>(*m,"atom.iprop.");
      CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
      CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
      CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->getProp<std::int64_t>("foo1")==6);
      CHECK(m->getAtomWithIdx(2)->getProp<std::int64_t>("foo2")==9);
      CHECK(m->getAtomWithIdx(1)->getProp<std::int64_t>("foo3")==1);
  }
  SECTION("basics: dprops") {
      auto m = "COC"_smiles;
      REQUIRE(m);
      m->setProp("atom.dprop.foo1","1   6 9");
      m->setProp("atom.dprop.foo2","3 n/a   9");
      m->setProp("atom.dprop.foo3","[?]  5 1 ?");
      FileParserUtils::applyMolListPropsToAtoms<double>(*m,"atom.dprop.");
      CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
      CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
      CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->getProp<double>("foo1")==6);
      CHECK(m->getAtomWithIdx(2)->getProp<double>("foo2")==9);
      CHECK(m->getAtomWithIdx(1)->getProp<double>("foo3")==1);
  }
    SECTION("basics: props") {
      auto m = "COC"_smiles;
      REQUIRE(m);
      m->setProp("atom.prop.foo1","1   6 9");
      m->setProp("atom.prop.foo2","3 n/a   9");
      m->setProp("atom.prop.foo3","[?]  5 1 ?");
      FileParserUtils::applyMolListPropsToAtoms<std::string>(*m,"atom.prop.");
      CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
      CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
      CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo1")=="6");
      CHECK(m->getAtomWithIdx(2)->getProp<std::string>("foo2")=="9");
      CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo3")=="1");
  }
  SECTION("basics: bprops") {
      auto m = "COC"_smiles;
      REQUIRE(m);
      m->setProp("atom.bprop.foo1","1   0 0");
      m->setProp("atom.bprop.foo2","0 n/a   1");
      m->setProp("atom.bprop.foo3","[?]  0 1 ?");
      FileParserUtils::applyMolListPropsToAtoms<bool>(*m,"atom.bprop.");
      CHECK(m->getAtomWithIdx(0)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo1"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo2"));
      CHECK(!m->getAtomWithIdx(1)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(2)->hasProp("foo2"));
      CHECK(m->getAtomWithIdx(0)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->hasProp("foo3"));
      CHECK(!m->getAtomWithIdx(2)->hasProp("foo3"));
      CHECK(m->getAtomWithIdx(1)->getProp<bool>("foo1")==false);
      CHECK(m->getAtomWithIdx(2)->getProp<bool>("foo2")==true);
      CHECK(m->getAtomWithIdx(1)->getProp<bool>("foo3")==true);
  }
}
TEST_CASE("processMolPropertyLists", "[atom_list_properties]") {
  SECTION("basics") {
      auto m = "COC"_smiles;
      REQUIRE(m);
      m->setProp("atom.iprop.foo1","1   6 9");
      m->setProp("atom.dprop.foo2","3 n/a   9");
      m->setProp("atom.prop.foo3","[?]  5 1 ?");
      m->setProp("atom.bprop.foo4","1 0 0");
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
      CHECK(m->getAtomWithIdx(1)->getProp<std::int64_t>("foo1")==6);
      CHECK(m->getAtomWithIdx(2)->getProp<double>("foo2")==9);
      CHECK(m->getAtomWithIdx(1)->getProp<std::string>("foo3")=="1");
      CHECK(m->getAtomWithIdx(1)->getProp<bool>("foo4")==false);
  }
}
