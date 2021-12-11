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
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/FileParsers/FileParsers.h>

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

TEST_CASE("Setting generic queries", "[substructure][generics]") {
  SECTION("cxsmiles") {
    auto m = "CO* |$;foo_p;ARY_p$|"_smiles;
    REQUIRE(m);
    auto atm = m->getAtomWithIdx(2);
    CHECK(atm->hasProp(common_properties::atomLabel));
    GenericGroups::setGenericQueriesFromProperties(*m);
    CHECK(atm->hasProp(common_properties::_QueryAtomGenericLabel));
    CHECK(atm->getProp<std::string>(
              common_properties::_QueryAtomGenericLabel) == "ARY");
    CHECK(!atm->hasProp(common_properties::atomLabel));
    // make sure we didn't remove the other atom label
    CHECK(m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    {
      auto smi = MolToCXSmiles(*m, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ATOM_LABELS);
      CHECK(smi == "*OC |$ARY_p;foo_p;$|");
      std::unique_ptr<ROMol> m2{SmilesToMol(smi)};
      REQUIRE(m2);
      GenericGroups::setGenericQueriesFromProperties(*m2);
      CHECK(m2->getAtomWithIdx(0)->hasProp(
          common_properties::_QueryAtomGenericLabel));
    }
    {
      auto ctab = MolToV3KMolBlock(*m);
      CHECK(ctab.find("ARY") != std::string::npos);
    }
  }
  SECTION("CTAB") {
    auto m = R"CTAB(
  Mrv2108 11042108512D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 3 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.0001 0.77 0 0
M  V30 2 N -1.3336 0.0001 0 0
M  V30 3 C -1.3336 -1.54 0 0
M  V30 4 C -0.0001 -2.31 0 0
M  V30 5 C 1.3336 -1.54 0 0
M  V30 6 N 1.3336 0.0001 0 0
M  V30 7 C -0.0001 2.31 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 3 4
M  V30 2 1 4 5
M  V30 3 2 1 2
M  V30 4 1 2 3
M  V30 5 1 1 6
M  V30 6 2 5 6
M  V30 7 1 1 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 -
M  V30 ATOMS=(1 7) -
M  V30 LABEL="AOX"
M  V30 2 DAT 0 ATOMS=(1 6) FIELDNAME=data -
M  V30 FIELDDISP="    1.3336    0.0001    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=value
M  V30 3 SUP 0 -
M  V30 ATOMS=(1 5) -
M  V30 LABEL="carbon"
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    auto atm = m->getAtomWithIdx(6);
    CHECK(getSubstanceGroups(*m).size() == 3);
    GenericGroups::setGenericQueriesFromProperties(*m);
    CHECK(atm->hasProp(common_properties::_QueryAtomGenericLabel));
    CHECK(atm->getProp<std::string>(
              common_properties::_QueryAtomGenericLabel) == "AOX");
    CHECK(getSubstanceGroups(*m).size() == 2);
    {
      auto smi = MolToCXSmiles(*m, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ATOM_LABELS);
      CHECK(smi == "Cc1ncccn1 |$AOX_p;;;;;;$|");
      std::unique_ptr<ROMol> m2{SmilesToMol(smi)};
      REQUIRE(m2);
      GenericGroups::setGenericQueriesFromProperties(*m2);
      CHECK(m2->getAtomWithIdx(0)->hasProp(
          common_properties::_QueryAtomGenericLabel));
    }

    {
      RWMol m2(*m);
      GenericGroups::convertGenericQueriesToSubstanceGroups(m2);
      CHECK(getSubstanceGroups(m2).size() == getSubstanceGroups(*m).size() + 1);
      CHECK(!m2.getAtomWithIdx(0)->hasProp(
          common_properties::_QueryAtomGenericLabel));
    }
    {
      auto ctab = MolToV3KMolBlock(*m);
      CHECK(ctab.find("AOX") != std::string::npos);
    }
  }
}

namespace {
bool no_match(const ROMol &mol, const std::vector<unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(ids);
  return false;
}
bool always_match(const ROMol &mol, const std::vector<unsigned int> &ids) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(ids);
  return true;
}
}  // namespace

TEST_CASE("generics are compatible with extraFinalChecks",
          "[substructure][generics]") {
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
      ps.extraFinalCheck = always_match;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, pr.second, ps));
      }
      ps.extraFinalCheck = no_match;
      for (const auto &pr : tests) {
        SmilesParserParams smilesParms;
        smilesParms.removeHs = false;
        std::unique_ptr<ROMol> mol{SmilesToMol(pr.first, smilesParms)};
        REQUIRE(mol);
        CHECK_THAT(*query, IsSubstructOf(*mol, 0, ps));
      }
    }
  }
}
