//
//  Copyright (C) 2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Tests of handling generics in substructure searching
//

#include "catch.hpp"

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/Generics.h>

using namespace RDKit;

class _IsSubstructOf : public Catch::MatcherBase<const ROMol &> {
  ROMol const *m_mol;
  SubstructMatchParameters m_ps;
  unsigned m_count;

 public:
  _IsSubstructOf(const ROMol &m, unsigned count, SubstructMatchParameters ps)
      : m_mol(&m), m_ps(std::move(ps)), m_count(count) {}

  bool match(const ROMol &query) const override {
    return SubstructMatch(*m_mol, query, m_ps).size() == m_count;
  }

  std::string describe() const override {
    std::ostringstream ss;
    ss << "does not have expected number of matches to "
       << MolToCXSmiles(*m_mol);
    return ss.str();
  }
};

static _IsSubstructOf IsSubstructOf(const ROMol &m, unsigned count,
                                    const SubstructMatchParameters &ps) {
  return _IsSubstructOf(m, count, ps);
}

TEST_CASE("edge cases", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C* |$;Foo_p$|"_smarts;
    REQUIRE(query);
    query->getAtomWithIdx(2)->setProp(common_properties::_QueryAtomGenericLabel,
                                      "Foo");
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 1}, {"O=CC=C", 1}, {"O=C", 0}};
    SubstructMatchParameters ps;
    ps.useGenericMatchers = true;
    for (const auto &pr : tests) {
      SmilesParserParams smilesParms;
      smilesParms.removeHs = false;
      std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
      REQUIRE(mol);
      CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
    }
  }
}

TEST_CASE("alkyl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Alkyl", "ALK"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CCC", 1},  {"O=CC=C", 0},        {"O=CCC=C", 0},
          {"O=CCO", 0},  {"O=CC", 1},          {"O=CC[H]", 1},
          {"O=C[H]", 0}, {"O=CC.CC(C=O)C", 2}, {"O=CCC1CC1", 0}};
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("alkoxy", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Alkoxy", "AOX"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CCC", 0},  {"O=COC", 1}, {"O=COC=C", 0}, {"O=COC=C", 0},
          {"O=COCO", 0}, {"O=CO", 0},  {"O=CO[H]", 0}, {"O=COC1CC1", 0}};
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("alkenyl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Alkenyl", "AEL"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CCC", 0},  {"O=CC=C", 1},    {"O=CCC=C", 1},    {"O=CCC#C", 0},
          {"O=CC=N", 0}, {"O=CC=C[H]", 1}, {"O=CCC1CC=C1", 0}};
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("alkynyl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Alkynyl", "AYL"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      REQUIRE(query);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CCC", 0},  {"O=CC#C", 1},    {"O=CCC#C", 1},    {"O=CCC=C", 0},
          {"O=CC#N", 0}, {"O=CC#C[H]", 1}, {"O=CCC1CC#C1", 0}};
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("cycloalkyl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Carbocycloalkyl", "CAL"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CCC", 0},           {"O=CC1CC1", 1},
          {"O=CC1C(O)C1", 1},     {"O=CCC1CC1", 0},
          {"O=CC1C=C1", 0},       {"O=CC1OC1", 0},
          {"O=CC1CC(=O)C1", 1},   {"O=CC1CC2CCC1CC2", 1},
          {"O=CC1CC2CCC1NC2", 0}, {"O=CC1CCCC2CC=CCC12", 0}};
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("cycloalkenyl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Carbocycloalkenyl", "CEL"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CCC", 0},
          {"O=Cc1ccccc1", 1},
          {"O=CC1CC=CCC1", 1},
          {"O=CC1CC2CCC1CC2", 0},
          {"O=CC1CC2CCC1C=C2", 0},  // <- one of the SSSR rings doesn't match
          {"O=CC1CC2C=CC1NC2", 0},
          {"O=CC1CCCC2=C1CCCC2", 1},
          {"O=CC1CC=CC2C1C=CC1CCC=CC21", 1},
          {"O=CC1CC=CC2C1C=CC1CNC=CC21", 0}};
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("carboaryl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Carboaryl", "ARY"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1CCCCC1", 0},          {"O=Cc1ccccc1", 1},
          {"O=CC1=CC=CC2=C1CCCC2", 0}, {"O=Cc1cccnc1", 0},
          {"O=CC1=CC=CC2=C1CCNC2", 0},
      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("carbocyclic", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Carbocyclic", "CBC"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1CCCCC1", 1},    {"O=Cc1ccccc1", 1},
          {"O=Cc1ccc(O)cc1", 1}, {"O=CC1=CC=CC2=C1CCCC2", 1},
          {"O=Cc1cccnc1", 0},    {"O=CC1=CC=CC2=C1CCNC2", 0},
      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
TEST_CASE("cyclic", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Cyclic", "CYC"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1CCC1", 1},      {"O=C1CCCC1", 0},   {"O=CCC1CCC1", 0},
          {"O=CC1CNC1", 1},      {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1},
          {"O=CC1CC2CCC1CC2", 1}

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("acyclic", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Acyclic", "ACY"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1CC1", 0},
          {"O=CCCCC1CCC1", 0},
          {"O=CCCCC", 1},
          {"O=CCCCO", 1},

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("carboacyclic", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Carboacyclic", "ABC"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1CC1", 0},
          {"O=CCCCC1CCC1", 0},
          {"O=CCCCC", 1},
          {"O=CCCCO", 0},

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("heteroacyclic", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Heteroacyclic", "AHC"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1OC1", 0},
          {"O=CCCCC1CCC1", 0},
          {"O=CCCCC", 0},
          {"O=CCCCO", 1},

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("heterocyclic", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Heterocyclic", "CHC"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1OC1", 1},        {"O=CC1CC1", 0},    {"O=CCC1COC1", 0},
          {"O=CC1CCC2C1CNC2", 1}, {"O=Cc1ccccc1", 0}, {"O=Cc1cccnc1", 1},

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("heteroaryl", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"Heteroaryl", "HAR"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CC1OC1", 0},          {"O=CCC1COC1", 0},  {"O=CC1CCC2C1CNC2", 0},
          {"O=Cc1cccc2c1ccnc2", 1}, {"O=Cc1ccccc1", 0}, {"O=Cc1cccnc1", 1},

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}

TEST_CASE("no carbon ring", "[substructure][generics]") {
  SECTION("basics") {
    auto query = "O=C*"_smarts;
    REQUIRE(query);
    std::vector<std::string> labels = {"NoCarbonRing", "CXX"};
    for (auto label : labels) {
      query->getAtomWithIdx(2)->setProp(
          common_properties::_QueryAtomGenericLabel, label);
      std::vector<std::pair<std::string, unsigned>> tests = {
          {"O=CN1NC1", 0},       {"O=CN1NN1", 1},       {"O=CCN1NNN1", 0},
          {"O=CN1NNNN2N1N2", 1}, {"O=CN1NNNN2N1C2", 0}, {"O=CNNN", 0},
          {"O=C1NNN1", 0},

      };
      SubstructMatchParameters ps;
      ps.useGenericMatchers = true;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
    }
  }
}
