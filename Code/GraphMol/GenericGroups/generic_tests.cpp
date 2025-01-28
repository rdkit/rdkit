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

#include <catch2/catch_all.hpp>

#include <tuple>
#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MarvinParse/MarvinParser.h>

using namespace RDKit;

class _IsSubstructOf : public Catch::Matchers::MatcherBase<const ROMol &> {
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

std::string getCXSmilesQuery(const std::string groupToTest) {
  std::ostringstream out;

  out << "*C=O |$" << groupToTest << ";;$|";
  return out.str();
}

std::string getMolQuery(const std::string groupToTest,
                        bool supGroupFlag = false) {
  std::ostringstream out;

  if (!supGroupFlag) {
    // pad the groupToTest with spaces to make it 3 characters long
    std::string groupToTestPadded = groupToTest;
    while (groupToTestPadded.length() < 3) {
      groupToTestPadded += " ";
    }

    out << R"MOLQQQ(
  Ketcher  3212315322D 1   1.00000     0.00000     0

  3  2  0  0  0  0  0  0  0  0999 V2000
    9.2250   -1.6000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.0910   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.9570   -1.6000    0.0000 )MOLQQQ"
        << groupToTestPadded << R"MOLQQQ( 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0     0  0
  2  3  1  0     0  0
M  END
)MOLQQQ";

  } else {
    out << R"MOLQQQ(
  Ketcher  3202311 72D 1   1.00000     0.00000     0

  3  2  0  0  0  0  0  0  0  0999 V2000
    6.3840   61.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2500   61.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1160   61.5000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0     0  0
  2  3  1  0     0  0
M  STY  1   1 SUP
M  SLB  1   1   1
M  SAL   1  1   3
M  SBL   1  1   2
M  SMT   1 )MOLQQQ"
        << groupToTest << R"MOLQQQ(
M  END
)MOLQQQ";
  }

  return out.str();
}

std::string getMRVQuery(const std::string groupToTest) {
  std::ostringstream out;

  out << R"MRVQQQ(<cml xmlns="http://www.chemaxon.com" version="ChemAxon file format v20.20.0, generated by vunknown" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.chemaxon.com http://www.chemaxon.com/marvin/schema/mrvSchema_20_20_0.xsd">
    <MDocument>
        <MChemicalStruct>
            <molecule molID="m1">
                <atomArray>
                    <atom id="a1" elementType="O" x2="-1.3336693410482567" y2="0.3850084702795158" lonePair="2"/>
                    <atom id="a2" elementType="C" x2="0" y2="-0.3850084702795158"/>
                    <atom id="a3" elementType="R" x2="1.3336693410482567" y2="0.3850084702795158" sgroupRef="sg1"/>
                </atomArray>
                <bondArray>
                    <bond id="b1" atomRefs2="a1 a2" order="2"/>
                    <bond id="b2" atomRefs2="a2 a3" order="1"/>
                </bondArray>
                <molecule molID="m2" id="sg1" role="SuperatomSgroup" title=")MRVQQQ"
      << groupToTest << R"MRVQQQ(">
                    <atomArray>
                        <atom id="a4" elementType="*" x2="1.3336693410482567" y2="0.3850084702795158" sgroupAttachmentPoint="1"/>
                    </atomArray>
                    <AttachmentPointArray>
                        <attachmentPoint atom="a4" order="1" bond="b2"/>
                    </AttachmentPointArray>
                </molecule>
            </molecule>
        </MChemicalStruct>
        <MElectronContainer occupation="0 0" radical="0" id="o1">
            <MElectron atomRefs="m1.a1" difLoc="0.0 0.0 0.0"/>
            <MElectron atomRefs="m1.a1" difLoc="0.0 0.0 0.0"/>
        </MElectronContainer>
        <MElectronContainer occupation="0 0" radical="0" id="o2">
            <MElectron atomRefs="m1.a1" difLoc="0.0 0.0 0.0"/>
            <MElectron atomRefs="m1.a1" difLoc="0.0 0.0 0.0"/>
        </MElectronContainer>
    </MDocument>
</cml>
)MRVQQQ";

  return out.str();
}

void runTest(RWMol *query,
             std::vector<std::pair<std::string, unsigned>> &tests) {
  GenericGroups::setGenericQueriesFromProperties(*query, true, true);

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

void runCxsmilesTest(const std::string groupToTest,
                     std::vector<std::pair<std::string, unsigned>> &tests,
                     bool manualInsertOfGroupFlag = false) {
  std::unique_ptr<RWMol> query;
  if (!manualInsertOfGroupFlag) {
    std::string queryString = getCXSmilesQuery(groupToTest);
    query = std::unique_ptr<RWMol>(SmartsToMol(queryString));
  } else {
    query = std::unique_ptr<RWMol>(SmartsToMol("*C=O |$;;$|"));
    query->getAtomWithIdx(0)->setProp(RDKit::common_properties::atomLabel,
                                      groupToTest);
  }
  runTest(query.get(), tests);
}

void runMolTest(const std::string groupToTest,
                std::vector<std::pair<std::string, unsigned>> &tests,
                bool supGroupFlag) {
  std::string queryString = getMolQuery(groupToTest, supGroupFlag);
  auto query =
      std::unique_ptr<RWMol>(MolBlockToMol(queryString, false, false, false));
  runTest(query.get(), tests);
}

void runMRVTest(const std::string groupToTest,
                std::vector<std::pair<std::string, unsigned>> &tests) {
  std::string queryString = getMRVQuery(groupToTest);
  auto query = std::unique_ptr<RWMol>(MrvBlockToMol(queryString, false, false));
  runTest(query.get(), tests);
}

void runTest(const std::string groupToTest,
             std::vector<std::pair<std::string, unsigned>> &tests) {
  runCxsmilesTest(groupToTest, tests, false);
  runCxsmilesTest(groupToTest, tests, true);
  runMolTest(groupToTest, tests, true);
  if (groupToTest.size() <= 3) {
    runMolTest(groupToTest, tests, false);
  }
  runMolTest(groupToTest, tests, true);
  runMRVTest(groupToTest, tests);
}

TEST_CASE("alkyl", "[substructure][generics]") {
  SECTION("Basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 1},  {"O=CC=C", 0},        {"O=CCC=C", 0},
        {"O=CCO", 0},  {"O=CC", 1},          {"O=CC[H]", 1},
        {"O=C[H]", 0}, {"O=CC.CC(C=O)C", 2}, {"O=CCC1CC1", 0}};
    runTest("Alkyl", tests);
    runTest("ALK", tests);
  }
}
TEST_CASE("alkyl or h", "[substructure][generics]") {
  SECTION("Basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 1},  {"O=CC=C", 0},        {"O=CCC=C", 0},
        {"O=CCO", 0},  {"O=CC", 1},          {"O=CC[H]", 1},
        {"O=C[H]", 1}, {"O=CC.CC(C=O)C", 2}, {"O=CCC1CC1", 0}};
    runTest("AlkylH", tests);
    runTest("ALH", tests);
  }
}

TEST_CASE("carboacyclic", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CC1", 0}, {"O=CCCCC1CCC1", 0},
        {"O=CCCCC", 1},  {"O=CC([H])([H])([H])", 1},
        {"O=CCCCO", 0},  {"O=C[H]", 0},
    };
    runTest("Carboacyclic", tests);
    runTest("ABC", tests);
  }
}
TEST_CASE("carboacyclic or h", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CC1", 0}, {"O=CCCCC1CCC1", 0},
        {"O=CCCCC", 1},  {"O=CC([H])([H])([H])", 1},
        {"O=CCCCO", 0},  {"O=C[H]", 1},
    };
    runTest("CarboacyclicH", tests);
    runTest("ABH", tests);
  }
}
TEST_CASE("acyclic", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CC1", 0}, {"O=CCCCC1CCC1", 0}, {"O=CCCCC", 1},
        {"O=CCCCO", 1},  {"O=C[H]", 0},
    };
    runTest("Acyclic", tests);
    runTest("ACY", tests);
  }
}

TEST_CASE("acyclic or h", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CC1", 0}, {"O=CCCCC1CCC1", 0}, {"O=CCCCC", 1},
        {"O=CCCCO", 1},  {"O=C[H]", 1},
    };
    runTest("AcyclicH", tests);
    runTest("ACH", tests);
  }
}

TEST_CASE("Alkenyl", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},  {"O=CC=C", 1},    {"O=CCC=C", 1},     {"O=CCC#C", 0},
        {"O=CC=N", 0}, {"O=CC=C[H]", 1}, {"O=CCC1CC=C1", 0}, {"O=C[H]", 0}};

    runTest("Alkenyl", tests);
    runTest("AEL", tests);
  }
}
TEST_CASE("Alkenyl or h", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},  {"O=CC=C", 1},    {"O=CCC=C", 1},     {"O=CCC#C", 0},
        {"O=CC=N", 0}, {"O=CC=C[H]", 1}, {"O=CCC1CC=C1", 0}, {"O=C[H]", 1}};

    runTest("AlkenylH", tests);
    runTest("AEH", tests);
  }
}

TEST_CASE("Alkynyl", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},       {"O=CC#C", 1},     {"O=CCC#C", 1},
        {"O=CCC=C", 0},     {"O=CC#N", 0},     {"O=CC#C[H]", 1},
        {"O=CCC1CC#C1", 0}, {"O=C([H])[H]", 0}};
    runTest("Alkynyl", tests);
    runTest("AYL", tests);
  }
}
TEST_CASE("Alkynyl or h", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},       {"O=CC#C", 1},     {"O=CCC#C", 1},
        {"O=CCC=C", 0},     {"O=CC#N", 0},     {"O=CC#C[H]", 1},
        {"O=CCC1CC#C1", 0}, {"O=C([H])[H]", 2}};

    runTest("AlkynylH", tests);
    runTest("AYH", tests);
  }
}

TEST_CASE("Heteroacyclic", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1OC1", 0}, {"O=CCCCC1CCC1", 0}, {"O=CCCCC", 0},
        {"O=CCCCO", 1},  {"O=C[H]", 0},
    };
    runTest("Heteroacyclic", tests);
    runTest("AHC", tests);
  }
}
TEST_CASE("Heteroacyclic or h", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1OC1", 0}, {"O=CCCCC1CCC1", 0}, {"O=CCCCC", 0},
        {"O=CCCCO", 1},  {"O=C[H]", 1},
    };
    runTest("HeteroacyclicH", tests);
    runTest("AHH", tests);
  }
}

TEST_CASE("alkoxy", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},   {"O=COC", 1},     {"O=COC=C", 0},
        {"O=COC=C", 0}, {"O=COCO", 0},    {"O=CO", 0},
        {"O=CO[H]", 0}, {"O=COC1CC1", 0}, {"O=C([H])[H]", 0}};

    runTest("Alkoxy", tests);
    runTest("AOX", tests);
  }
}

TEST_CASE("alkoxy or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},   {"O=COC", 1},     {"O=COC=C", 0},
        {"O=COC=C", 0}, {"O=COCO", 0},    {"O=CO", 0},
        {"O=CO[H]", 0}, {"O=COC1CC1", 0}, {"O=C([H])[H]", 2}};

    runTest("AlkoxyH", tests);
    runTest("AOH", tests);
  }
}

TEST_CASE("carboaryl", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCCCC1", 0},          {"O=Cc1ccccc1", 1},
        {"O=CC1=CC=CC2=C1CCCC2", 0}, {"O=Cc1cccnc1", 0},
        {"O=CC1=CC=CC2=C1CCNC2", 0}, {"O=C([H])[H]", 0}};

    runTest("Carboaryl", tests);
    runTest("ARY", tests);
  }
}
TEST_CASE("carboaryl or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCCCC1", 0},          {"O=Cc1ccccc1", 1},
        {"O=CC1=CC=CC2=C1CCCC2", 0}, {"O=Cc1cccnc1", 0},
        {"O=CC1=CC=CC2=C1CCNC2", 0}, {"O=C([H])[H]", 2}};

    runTest("CarboarylH", tests);
    runTest("ARH", tests);
  }
}

TEST_CASE("cycloalkyl", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},           {"O=CC1CC1", 1},
        {"O=CC1C(O)C1", 1},     {"O=CCC1CC1", 0},
        {"O=CC1C=C1", 0},       {"O=CC1OC1", 0},
        {"O=CC1CC(=O)C1", 1},   {"O=CC1CC2CCC1CC2", 1},
        {"O=CC1CC2CCC1NC2", 0}, {"O=CC1CCCC2CC=CCC12", 0},
        {"O=C([H])[H]", 0}};
    runTest("Carbocycloalkyl", tests);
    runTest("CAL", tests);
  }
}

TEST_CASE("cycloalkyl or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},           {"O=CC1CC1", 1},
        {"O=CC1C(O)C1", 1},     {"O=CCC1CC1", 0},
        {"O=CC1C=C1", 0},       {"O=CC1OC1", 0},
        {"O=CC1CC(=O)C1", 1},   {"O=CC1CC2CCC1CC2", 1},
        {"O=CC1CC2CCC1NC2", 0}, {"O=CC1CCCC2CC=CCC12", 0},
        {"O=C([H])[H]", 2}};
    runTest("CarbocycloalkylH", tests);
    runTest("CAH", tests);
  }
}

TEST_CASE("carbocyclic", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCCCC1", 1},    {"O=Cc1ccccc1", 1},
        {"O=Cc1ccc(O)cc1", 1}, {"O=CC1=CC=CC2=C1CCCC2", 1},
        {"O=Cc1cccnc1", 0},    {"O=CC1=CC=CC2=C1CCNC2", 0},
        {"O=C([H])[H]", 0},
    };
    runTest("Carbocyclic", tests);
    runTest("CBC", tests);
  }
}

TEST_CASE("carbocyclic or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCCCC1", 1},    {"O=Cc1ccccc1", 1},
        {"O=Cc1ccc(O)cc1", 1}, {"O=CC1=CC=CC2=C1CCCC2", 1},
        {"O=Cc1cccnc1", 0},    {"O=CC1=CC=CC2=C1CCNC2", 0},
        {"O=C([H])[H]", 2},
    };
    runTest("CarbocyclicH", tests);
    runTest("CBH", tests);
  }
}

TEST_CASE("cycloalkenyl", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},
        {"O=Cc1ccccc1", 1},
        {"O=CC1CC=CCC1", 1},
        {"O=CC1CC2CCC1CC2", 0},
        {"O=CC1CC2CCC1C=C2", 0},  // <- one of the SSSR rings doesn't match
        {"O=CC1CC2C=CC1NC2", 0},
        {"O=CC1CCCC2=C1CCCC2", 1},
        {"O=CC1CC=CC2C1C=CC1CCC=CC21", 1},
        {"O=CC1CC=CC2C1C=CC1CNC=CC21", 0},
        {"O=C([H])[H]", 0},
    };
    runTest("Carbocycloalkenyl", tests);
    runTest("CEL", tests);
  }
}

TEST_CASE("cycloalkenyl or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCC", 0},
        {"O=Cc1ccccc1", 1},
        {"O=CC1CC=CCC1", 1},
        {"O=CC1CC2CCC1CC2", 0},
        {"O=CC1CC2CCC1C=C2", 0},  // <- one of the SSSR
                                  // rings doesn't match
        {"O=CC1CC2C=CC1NC2", 0},
        {"O=CC1CCCC2=C1CCCC2", 1},
        {"O=CC1CC=CC2C1C=CC1CCC=CC21", 1},
        {"O=CC1CC=CC2C1C=CC1CNC=CC21", 0},
        {"O=C([H])[H]", 2},
    };
    runTest("CarbocycloalkenylH", tests);
    runTest("CEH", tests);
  }
}

TEST_CASE("heterocyclic", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1OC1", 1},        {"O=CC1CC1", 0},    {"O=CCC1COC1", 0},
        {"O=CC1CCC2C1CNC2", 1}, {"O=Cc1ccccc1", 0}, {"O=Cc1cccnc1", 1},
        {"O=C([H])[H]", 0},
    };
    runTest("Heterocyclic", tests);
    runTest("CHC", tests);
  }
}
TEST_CASE("heterocyclic or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1OC1", 1},        {"O=CC1CC1", 0},    {"O=CCC1COC1", 0},
        {"O=CC1CCC2C1CNC2", 1}, {"O=Cc1ccccc1", 0}, {"O=Cc1cccnc1", 1},
        {"O=C([H])[H]", 2},
    };
    runTest("HeterocyclicH", tests);
    runTest("CHH", tests);
  }
}

TEST_CASE("no carbon ring", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CN1NC1", 0},       {"O=CN1NN1", 1},       {"O=CCN1NNN1", 0},
        {"O=CN1NNNN2N1N2", 1}, {"O=CN1NNNN2N1C2", 0}, {"O=CNNN", 0},
        {"O=C1NNN1", 0},       {"O=C([H])[H]", 0},
    };
    runTest("NoCarbonRing", tests);
    runTest("CXX", tests);
  }
}

TEST_CASE("no carbon ring or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CN1NC1", 0},       {"O=CN1NN1", 1},       {"O=CCN1NNN1", 0},
        {"O=CN1NNNN2N1N2", 1}, {"O=CN1NNNN2N1C2", 0}, {"O=CNNN", 0},
        {"O=C1NNN1", 0},       {"O=C([H])[H]", 2},
    };
    runTest("NoCarbonRingH", tests);
    runTest("CXH", tests);
  }
}

TEST_CASE("cyclic", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCC1", 1},       {"O=C1CCCC1", 0},   {"O=CCC1CCC1", 0},
        {"O=CC1CNC1", 1},       {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1},
        {"O=CC1CC2CCC1CC2", 1}, {"O=C([H])[H]", 0},
    };
    runTest("Cyclic", tests);
    runTest("CYC", tests);
  }
}

TEST_CASE("cyclic or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCC1", 1},       {"O=C1CCCC1", 0},   {"O=CCC1CCC1", 0},
        {"O=CC1CNC1", 1},       {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1},
        {"O=CC1CC2CCC1CC2", 1}, {"O=C([H])[H]", 2},
    };
    runTest("CyclicH", tests);
    runTest("CYH", tests);
  }
}

TEST_CASE("heteroaryl", "[substructure][generics]") {
  SECTION("basics") {
    std::string query = "*C=O |$HAH;;$|";
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1OC1", 0},          {"O=CCC1COC1", 0},  {"O=CC1CCC2C1CNC2", 0},
        {"O=Cc1cccc2c1ccnc2", 1}, {"O=Cc1ccccc1", 0}, {"O=Cc1cccnc1", 1},
        {"O=C([H])[H]", 0},
    };
    runTest("Heteroaryl", tests);
    runTest("HAR", tests);
  }
}

TEST_CASE("heteroaryl or H", "[substructure][generics]") {
  SECTION("basics") {
    std::string query = "*C=O |$HAH;;$|";
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1OC1", 0},          {"O=CCC1COC1", 0},  {"O=CC1CCC2C1CNC2", 0},
        {"O=Cc1cccc2c1ccnc2", 1}, {"O=Cc1ccccc1", 0}, {"O=Cc1cccnc1", 1},
        {"O=C([H])[H]", 2},
    };
    runTest("HeteroarylH", tests);
    runTest("HAH", tests);
  }
}

TEST_CASE("group", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCC1", 1},       {"O=CCCCC", 1},     {"O=CCC1CCC1", 1},
        {"O=CC1CNC1", 1},       {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1},
        {"O=CC1CC2CCC1CC2", 1}, {"O=C([H])[H]", 0},
    };
    runTest("Group", tests);
    runTest("G", tests);
  }
}

TEST_CASE("group or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCC1", 1},       {"O=CCCCC", 1},     {"O=CCC1CCC1", 1},
        {"O=CC1CNC1", 1},       {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1},
        {"O=CC1CC2CCC1CC2", 1}, {"O=C([H])[H]", 2},
    };
    runTest("GroupH", tests);
    runTest("GH", tests);
  }
}

TEST_CASE("group with ring", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CCCCC", 0},     {"O=CCC1CCC1", 1},   {"O=CC1CNC1", 1},
        {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1}, {"O=CC1CC2CCC1CC2", 1},
        {"O=C([H])[H]", 0},
    };
    runTest("Group*", tests);
    runTest("G*", tests);
  }
}

TEST_CASE("group with ring or H", "[substructure][generics]") {
  SECTION("basics") {
    std::vector<std::pair<std::string, unsigned>> tests = {
        {"O=CC1CCC1", 1},       {"O=CCCCC", 0},     {"O=CCC1CCC1", 1},
        {"O=CC1CNC1", 1},       {"O=Cc1ccccc1", 1}, {"O=CC1CC2C1C2", 1},
        {"O=CC1CC2CCC1CC2", 1}, {"O=C([H])[H]", 2},
    };
    runTest("GroupH*", tests);
    runTest("GH*", tests);
  }
}

TEST_CASE("from mol blocks using atom labels") {
  SECTION("basics") {
    auto qry = R"CTAB(
  Mrv2211 08032317212D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.5833 1.9167 0 0
M  V30 2 ARY -5.2497 2.6867 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(qry);

    SubstructMatchParameters ps;
    ps.useGenericMatchers = true;
    auto m1 = "CC"_smiles;
    REQUIRE(m1);
    auto m2 = "Cc1ccccc1"_smiles;
    REQUIRE(m2);

    CHECK(SubstructMatch(*m1, *qry, ps).empty());
    CHECK(SubstructMatch(*m2, *qry, ps).empty());

    GenericGroups::setGenericQueriesFromProperties(*qry, true);
    CHECK(SubstructMatch(*m1, *qry, ps).empty());
    CHECK(SubstructMatch(*m2, *qry, ps).size() == 1);
  }
}

void testGenericQueries() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing generic queries: "
                          "Set GenericQueriesFromProperties"
                       << std::endl;
  auto query = R"CTAB(alkquery.mol
Generated by WebMolKit

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 1 0 1
M  V30 BEGIN ATOM
M  V30 1 O 7.4214 -9.6318 0.0000 0
M  V30 2 * 8.5430 -10.2809 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 LABEL=ALK ATOMS=(1 2)
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;

  TEST_ASSERT(query);

  std::string target_smi{"c1c(OF)cccc1"};
  std::unique_ptr<RWMol> t{SmilesToMol(target_smi)};
  TEST_ASSERT(t);

  MolOps::AdjustQueryParameters params;
  auto m = new RWMol(*query);

  TEST_ASSERT(m);

  MolOps::adjustQueryProperties(*m, &params);
  SubstructMatchParameters params_match;
  params_match.useQueryQueryMatches = true;
  params_match.useGenericMatchers = true;
  params_match.maxMatches = 1;
  auto matchVect = SubstructMatch(*t, *m, params_match);
  TEST_ASSERT(matchVect.size() == 1);

  delete m;
  m = new RWMol(*query);
  GenericGroups::adjustQueryPropertiesWithGenericGroups(*m, &params);
  matchVect = SubstructMatch(*t, *m, params_match);
  TEST_ASSERT(matchVect.size() == 0);
  delete m;

  m = new RWMol(*query);
  GenericGroups::adjustQueryPropertiesWithGenericGroups(*m);
  matchVect = SubstructMatch(*t, *m, params_match);
  TEST_ASSERT(matchVect.size() == 0);
  delete m;
}
