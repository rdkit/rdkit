//
//  Copyright (C) 2020-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <utility>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

using matchCase = std::tuple<std::string, std::string, bool, bool>;

class _IsSubstructOf
    : public Catch::Matchers::MatcherBase<const std::string &> {
  ROMol const *m_query;
  std::string m_description;
  SubstructMatchParameters m_ps;

 public:
  _IsSubstructOf(const ROMol &m, std::string description)
      : m_query(&m), m_description(std::move(description)) {}

  _IsSubstructOf(const ROMol &m, std::string description,
                 SubstructMatchParameters ps)
      : m_query(&m),
        m_description(std::move(description)),
        m_ps(std::move(ps)) {}

  bool match(const std::string &smiles) const override {
    std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
    return !SubstructMatch(*mol, *m_query, m_ps).empty();
  }

  std::string describe() const override {
    std::ostringstream ss;
    ss << "is not a substructure of " << m_description;
    return ss.str();
  }
};

static _IsSubstructOf IsSubstructOf(const ROMol &m, const std::string &smarts,
                                    const SubstructMatchParameters &ps) {
  return _IsSubstructOf(m, smarts, ps);
}

static _IsSubstructOf IsSubstructOf(const ROMol &m, const std::string &smarts) {
  return _IsSubstructOf(m, smarts);
}

TEST_CASE("handling of bondStereoCare in adjustQueryProperties") {
  SECTION("fully specified") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
  }
  SECTION("fully unspecified") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREONONE);
  }
  SECTION("partially unspecified") {
    std::vector<std::string> mbs = {R"CTAB(keep
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(keep
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(remove
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2 STBOX=0
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(remove
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 
M  V30 2 C -5.6979 2.8332 0 0 STBOX=1
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB",
                                    R"CTAB(remove
  Mrv1810 01292006422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.0316 2.0632 0 0 STBOX=1
M  V30 2 C -5.6979 2.8332 0 0
M  V30 3 O -4.3642 2.0632 0 0
M  V30 4 F -8.3653 2.8332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 1 4
M  V30 3 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"};
    for (const auto &mb : mbs) {
      std::unique_ptr<RWMol> mol{MolBlockToMol(mb)};
      REQUIRE(mol);
      REQUIRE(mol->getBondBetweenAtoms(0, 1));
      CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
            Bond::BondStereo::STEREOE);
      MolOps::AdjustQueryParameters ps;
      ps.useStereoCareForBonds = true;
      MolOps::adjustQueryProperties(*mol, &ps);
      if (mol->getProp<std::string>(common_properties::_Name) == "keep") {
        CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
              Bond::BondStereo::STEREOE);
      } else {
        CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
              Bond::BondStereo::STEREONONE);
      }
    }
  }
  SECTION("V2000") {
    auto mol = R"CTAB(basic test
  Mrv1810 01292015042D          

  4  3  0  0  0  0            999 V2000
   -3.7669    1.1053    0.0000 C   0  0  0  0  1  0  0  0  0  0  0  0
   -3.0524    1.5178    0.0000 C   0  0  0  0  1  0  0  0  0  0  0  0
   -2.3380    1.1053    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4814    1.5178    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  2  2  0  0  0  0
M  END
)CTAB"_ctab;
    REQUIRE(mol);
    CHECK(mol->getAtomWithIdx(0)->hasProp(common_properties::molStereoCare));
    CHECK(mol->getAtomWithIdx(1)->hasProp(common_properties::molStereoCare));
    REQUIRE(mol->getBondBetweenAtoms(0, 1));
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    // property added by the CTAB parser:
    CHECK(mol->getBondBetweenAtoms(0, 1)->hasProp(
        common_properties::molStereoCare));
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
  }
  SECTION("molecule from SMILES") {
    auto mol = "C/C=C/C"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(2, 1));
    CHECK(mol->getBondBetweenAtoms(2, 1)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::AdjustQueryParameters ps;
    ps.useStereoCareForBonds = true;
    // since stereoCare is not set on the bond from SMILES,
    // stereochem will be removed:
    {
      RWMol molcp(*mol);
      MolOps::adjustQueryProperties(molcp, &ps);
      CHECK(molcp.getBondBetweenAtoms(2, 1)->getStereo() ==
            Bond::BondStereo::STEREONONE);
    }
    // but we can preserve it by setting the property:
    {
      RWMol molcp(*mol);
      molcp.getBondBetweenAtoms(2, 1)->setProp(common_properties::molStereoCare,
                                               1);
      MolOps::adjustQueryProperties(molcp, &ps);
      CHECK(molcp.getBondBetweenAtoms(2, 1)->getStereo() ==
            Bond::BondStereo::STEREOE);
    }
  }
}

TEST_CASE("adjustQueryParameters from JSON") {
  SECTION("basics") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.makeAtomsGeneric == false);
    CHECK(ps.makeBondsGeneric == false);
    CHECK(ps.makeBondsGenericFlags == MolOps::ADJUST_IGNORENONE);

    std::string json = R"JSON({"makeAtomsGeneric":true})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);

    CHECK(ps.makeAtomsGeneric == true);
    CHECK(ps.makeBondsGeneric == false);
    // the parsing updates the parameters, it doesn't replace them:

    json = R"JSON({"makeBondsGeneric":true,
      "makeBondsGenericFlags":"IGNOREDUMMIES|IGNORECHAINS"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);

    CHECK(ps.makeAtomsGeneric == true);
    CHECK(ps.makeBondsGeneric == true);
    CHECK(ps.makeBondsGenericFlags ==
          (MolOps::ADJUST_IGNOREDUMMIES | MolOps::ADJUST_IGNORECHAINS));
  }
  SECTION("useStereoCare") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.useStereoCareForBonds == false);

    std::string json = R"JSON({"useStereoCareForBonds":true})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.useStereoCareForBonds == true);
    json = R"JSON({"useStereoCareForBonds":false})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.useStereoCareForBonds == false);
  }
  SECTION("adjustHeavyDegreeFlags") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.adjustHeavyDegreeFlags != MolOps::ADJUST_IGNORENONE);
    std::string json = R"JSON({"adjustHeavyDegreeFlags":"IGNORENONE"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustHeavyDegreeFlags == MolOps::ADJUST_IGNORENONE);
  }
  SECTION("adjustRingCountFlags") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.adjustRingCountFlags != MolOps::ADJUST_IGNORERINGS);
    std::string json = R"JSON({"adjustRingCountFlags":"IGNORERINGS"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustRingCountFlags == MolOps::ADJUST_IGNORERINGS);
  }
  SECTION("makeAtomsGenericFlags") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.makeAtomsGenericFlags != MolOps::ADJUST_IGNOREALL);
    std::string json = R"JSON({"makeAtomsGenericFlags":"IGNOREALL"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.makeAtomsGenericFlags == MolOps::ADJUST_IGNOREALL);
  }
  SECTION("adjustRingChainFlags") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.adjustRingChainFlags != MolOps::ADJUST_IGNORENONDUMMIES);
    std::string json = R"JSON({"adjustRingChainFlags":"IGNORENONDUMMIES"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustRingChainFlags == MolOps::ADJUST_IGNORENONDUMMIES);
  }
  SECTION("bogus contents") {
    MolOps::AdjustQueryParameters ps;
    CHECK(ps.adjustDegree == true);
    CHECK(ps.adjustDegreeFlags ==
          (MolOps::ADJUST_IGNOREDUMMIES | MolOps::ADJUST_IGNORECHAINS));

    std::string json = R"JSON({"bogosity":true})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustDegree == true);

    json = R"JSON({"adjustDegree":"foo"})JSON";
    MolOps::parseAdjustQueryParametersFromJSON(ps, json);
    CHECK(ps.adjustDegree == true);

    json = R"JSON({"adjustDegreeFlags":"IGNORENONE|bogus"})JSON";
    CHECK_THROWS_AS(MolOps::parseAdjustQueryParametersFromJSON(ps, json),
                    ValueErrorException);
  }
}

TEST_CASE("MDL five-rings") {
  MolOps::AdjustQueryParameters ps =
      MolOps::AdjustQueryParameters::noAdjustments();
  ps.setMDLFiveRingAromaticity = true;
  SECTION("query details") {
    using extuple = std::tuple<std::string, std::string, std::string>;
    // clang-format off
    std::vector<extuple> examples = {
    // no queries, no change
    extuple{"adjustqueryprops_MDLfivering_1.mol","[#7H]1:[#6]:[#6]:[#6]:[#6]:1",""},
    // Q atom, no change
    extuple{"adjustqueryprops_MDLfivering_2.mol","[!#6&!#1]1:[#6]:[#6]:[#6]:[#6]:1",""},
    // A atom, this one changes
    extuple{"adjustqueryprops_MDLfivering_3.mol","[!#1]1:[#6]:[#6]:[#6]:[#6]:1","[!#1]1-,:[#6]=,:[#6]-,:[#6]=,:[#6]-,:1"},
    // NOTE that this is not technically correct according to the documentation, but if we make the bridging bond
    // aromatic then it won't match azulene in a normal RDKit molecule, which is certainly not the intent of this.
    extuple{"adjustqueryprops_MDLfivering_4.mol","[#6]12:[#6]:[#6]:[#6]:[#6]-1:[#6]:[#6]:[#6]:[#6]:[#6]:2",""},
    };
    // clang-format on
    for (auto tpl : examples) {
      if (std::get<2>(tpl).empty()) {
        std::get<2>(tpl) = std::get<1>(tpl);
      }
      auto fname = std::get<0>(tpl);
      std::string pathName = getenv("RDBASE");
      pathName += "/Code/GraphMol/test_data/";
      std::unique_ptr<RWMol> qry(MolFileToMol(pathName + fname));
      REQUIRE(qry);
      {
        RWMol cp(*qry);
        CHECK(std::get<1>(tpl) == MolToSmarts(cp));
        MolOps::adjustQueryProperties(cp, &ps);
        CHECK(std::get<2>(tpl) == MolToSmarts(cp));
      }
      {  // make sure ring-finding happens
        RWMol cp(*qry);
        cp.getRingInfo()->reset();
        CHECK(std::get<1>(tpl) == MolToSmarts(cp));
        MolOps::adjustQueryProperties(cp, &ps);
        CHECK(std::get<2>(tpl) == MolToSmarts(cp));
      }
    }
  }
  SECTION("edge cases") {
    SmilesParserParams smiles_ps;
    smiles_ps.sanitize = false;
    {
      std::unique_ptr<RWMol> qry{SmilesToMol("*1:c:c:c:c:1", smiles_ps)};
      REQUIRE(qry);
      std::unique_ptr<QueryAtom> qat(new QueryAtom(0));
      qat->setQuery(makeAAtomQuery());
      qat->setIsAromatic(true);
      qry->replaceAtom(0, qat.get());
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) == "[!#1]1-,:[#6]=,:[#6]-,:[#6]=,:[#6]-,:1");
    }
    {  // ring not fully aromatic
      std::unique_ptr<RWMol> qry{SmilesToMol("*1:c:C-c:c:1", smiles_ps)};
      REQUIRE(qry);
      std::unique_ptr<QueryAtom> qat(new QueryAtom(0));
      qat->setQuery(makeAAtomQuery());
      qat->setIsAromatic(true);
      qry->replaceAtom(0, qat.get());
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) == "[!#1]1:[#6]:[#6]-[#6]:[#6]:1");
    }
    {  // ring has additional dummy
      std::unique_ptr<RWMol> qry{SmilesToMol("*1:c:*:c:c:1", smiles_ps)};
      REQUIRE(qry);
      std::unique_ptr<QueryAtom> qat(new QueryAtom(0));
      qat->setQuery(makeAAtomQuery());
      qat->setIsAromatic(true);
      qry->replaceAtom(0, qat.get());
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) == "[!#1]1:[#6]:[#0]:[#6]:[#6]:1");
    }
    {  // query bond in ring
      std::unique_ptr<RWMol> qry{SmilesToMol("*1:c:c:c:c:1", smiles_ps)};
      REQUIRE(qry);
      std::unique_ptr<QueryAtom> qat(new QueryAtom(0));
      qat->setQuery(makeAAtomQuery());
      qat->setIsAromatic(true);
      qry->replaceAtom(0, qat.get());
      std::unique_ptr<QueryBond> qbnd(new QueryBond());
      qbnd->setBondType(Bond::BondType::SINGLE);
      qbnd->setQuery(makeBondOrderEqualsQuery(Bond::BondType::SINGLE));
      qry->replaceBond(0, qbnd.get());
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) == "[!#1]1-[#6]:[#6]:[#6]:[#6]:1");
    }
  }
}

TEST_CASE("conjugated five-rings") {
  MolOps::AdjustQueryParameters ps =
      MolOps::AdjustQueryParameters::noAdjustments();
  ps.adjustConjugatedFiveRings = true;
  SECTION("matching") {
    // clang-format off
    std::vector<matchCase> examples = {
    // 1,3 cyclopentadiene
    matchCase{"C1=CCC=C1","adjustqueryprops_fivering_1.mol",true,true},
    matchCase{"C1=CCC=C1","adjustqueryprops_fivering_2.mol",false,true},
    matchCase{"C1=CCC=C1","adjustqueryprops_fivering_3.mol",true,true},
    matchCase{"C1=CCC=C1","adjustqueryprops_fivering_4.mol",false,false},
    matchCase{"C1=CCC=C1","adjustqueryprops_fivering_5.mol",false,false},
    matchCase{"C1=CCC=C1","adjustqueryprops_fivering_6.mol",false,false},
    // pyrrole
    matchCase{"C1=CNC=C1","adjustqueryprops_fivering_1.mol",false,true},
    matchCase{"C1=CNC=C1","adjustqueryprops_fivering_2.mol",true,true},
    matchCase{"C1=CNC=C1","adjustqueryprops_fivering_3.mol",false,false},
    matchCase{"C1=CNC=C1","adjustqueryprops_fivering_4.mol",false,false},
    matchCase{"C1=CNC=C1","adjustqueryprops_fivering_5.mol",false,false},
    matchCase{"C1=CNC=C1","adjustqueryprops_fivering_6.mol",false,false},
    // thiophene
    matchCase{"C1=CSC=C1","adjustqueryprops_fivering_1.mol",false,false},
    matchCase{"C1=CSC=C1","adjustqueryprops_fivering_2.mol",true,true},
    matchCase{"C1=CSC=C1","adjustqueryprops_fivering_3.mol",false,false},
    matchCase{"C1=CSC=C1","adjustqueryprops_fivering_4.mol",true,true},
    matchCase{"C1=CSC=C1","adjustqueryprops_fivering_5.mol",false,false},
    matchCase{"C1=CSC=C1","adjustqueryprops_fivering_6.mol",true,true},
    // thiophene oxide
    matchCase{"C1=CS(=O)C=C1","adjustqueryprops_fivering_1.mol",false,false},
    matchCase{"C1=CS(=O)C=C1","adjustqueryprops_fivering_2.mol",false,true},
    matchCase{"C1=CS(=O)C=C1","adjustqueryprops_fivering_3.mol",false,false},
    matchCase{"C1=CS(=O)C=C1","adjustqueryprops_fivering_4.mol",false,true},
    matchCase{"C1=CS(=O)C=C1","adjustqueryprops_fivering_5.mol",false,false},
    matchCase{"C1=CS(=O)C=C1","adjustqueryprops_fivering_6.mol",true,true},
    // furan
    matchCase{"C1=COC=C1","adjustqueryprops_fivering_1.mol",false,true},
    matchCase{"C1=COC=C1","adjustqueryprops_fivering_2.mol",true,true},
    matchCase{"C1=COC=C1","adjustqueryprops_fivering_3.mol",false,false},
    matchCase{"C1=COC=C1","adjustqueryprops_fivering_4.mol",false,false},
    matchCase{"C1=COC=C1","adjustqueryprops_fivering_5.mol",false,false},
    matchCase{"C1=COC=C1","adjustqueryprops_fivering_6.mol",false,false},
    };
    // clang-format on
    for (const auto &tpl : examples) {
      auto fname = std::get<1>(tpl);
      std::string pathName = getenv("RDBASE");
      pathName += "/Code/GraphMol/test_data/";
      std::unique_ptr<RWMol> qry(MolFileToMol(pathName + fname));
      REQUIRE(qry);
      if (std::get<2>(tpl)) {
        CHECK_THAT(std::get<0>(tpl), IsSubstructOf(*qry, fname));
      } else {
        CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(*qry, fname));
      }
      MolOps::adjustQueryProperties(*qry, &ps);
      if (std::get<3>(tpl)) {
        CHECK_THAT(std::get<0>(tpl), IsSubstructOf(*qry, fname));
      } else {
        CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(*qry, fname));
      }
    }
  }
  SECTION("query details") {
    auto fname = "adjustqueryprops_fivering_2.mol";
    std::string pathName = getenv("RDBASE");
    pathName += "/Code/GraphMol/test_data/";
    std::unique_ptr<RWMol> qry(MolFileToMol(pathName + fname));
    REQUIRE(qry);
    auto smarts = MolToSmarts(*qry);
    CHECK(smarts == "[!#1]1:[#6]:[#6]:[#6]:[#6]:1");
    MolOps::adjustQueryProperties(*qry, &ps);
    smarts = MolToSmarts(*qry);
    CHECK(smarts == "[!#1]1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:1");
  }
  SECTION("some edge cases") {
    {
      auto qry = "C1=COCC1"_smiles;
      auto smarts = MolToSmarts(*qry);
      CHECK(smarts == "[#6]1=[#6]-[#8]-[#6]-[#6]-1");
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) == smarts);
    }
    {
      auto qry = "C1=CCC=C1"_smiles;
      auto smarts = MolToSmarts(*qry);
      CHECK(smarts == "[#6]1=[#6]-[#6]-[#6]=[#6]-1");
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) ==
            "[#6]1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:1");
    }
    {
      // conjugation (not bond order)
      auto qry = "C1=COOO1"_smiles;
      auto smarts = MolToSmarts(*qry);
      CHECK(smarts == "[#6]1=[#6]-[#8]-[#8]-[#8]-1");
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) ==
            "[#6]1-,=,:[#6]-,=,:[#8]-,=,:[#8]-,=,:[#8]-,=,:1");
    }
    {
      // conjugation (not bond order)
      auto qry = "O=C1C(=O)C(=O)C(=O)C1=O"_smiles;
      auto smarts = MolToSmarts(*qry);
      CHECK(smarts ==
            "[#8]=[#6]1-[#6](=[#8])-[#6](=[#8])-[#6](=[#8])-[#6]-1=[#8]");
      MolOps::adjustQueryProperties(*qry, &ps);
      CHECK(MolToSmarts(*qry) ==
            "[#8]=[#6]1-,=,:[#6](=[#8])-,=,:[#6](=[#8])-,=,:[#6](=[#8])-,=,:[#"
            "6]-,=,:1=[#8]");
    }
  }
  SECTION("edge cases: ring finding") {
    // test that ring finding happens:
    SmilesParserParams smiles_ps;
    smiles_ps.sanitize = false;
    std::unique_ptr<RWMol> qry{SmilesToMol("O1C=CC=C1", smiles_ps)};
    REQUIRE(qry);
    qry->updatePropertyCache();
    MolOps::setConjugation(*qry);
    qry->getRingInfo()->reset();
    CHECK(!qry->getRingInfo()->isInitialized());
    MolOps::adjustQueryProperties(*qry, &ps);
    CHECK(qry->getRingInfo()->isInitialized());
    CHECK(MolToSmarts(*qry) ==
          "[#8]1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:1");
  }
  SECTION("edge cases: ignore larger rings") {
    SmilesParserParams smiles_ps;
    smiles_ps.sanitize = false;
    std::unique_ptr<RWMol> qry{SmilesToMol("N1=CC=CC=C1", smiles_ps)};
    REQUIRE(qry);
    qry->updatePropertyCache();
    MolOps::setConjugation(*qry);
    MolOps::adjustQueryProperties(*qry, &ps);
    CHECK(MolToSmarts(*qry) == "[#7]1=[#6]-[#6]=[#6]-[#6]=[#6]-1");
  }
}

TEST_CASE("single bonds to degree-one neighbors") {
  MolOps::AdjustQueryParameters ps =
      MolOps::AdjustQueryParameters::noAdjustments();
  ps.adjustSingleBondsToDegreeOneNeighbors = true;
  SECTION("matching") {
    // clang-format off
    std::vector<matchCase> examples = {
      matchCase{"C2CCCc1c2nncc1","Cc1cnncc1",true,true},
      matchCase{"C2CCCc1c2nncc1","CCc1cnncc1",true,true},
      matchCase{"C2CCC(C)c1c2nncc1","CCc1cnncc1",true,true},
      matchCase{"c2cccc1c2nncc1","Cc1cnncc1",false,true},
      matchCase{"c2cccc1c2nncc1","CCc1cnncc1",false,false},
      matchCase{"c2ccc(C)c1c2nncc1","CCc1cnncc1",false,false},

      matchCase{"C2CCCc1[nH]ccc12","Cc1[nH]ccc1",true,true},
      matchCase{"C2CCCc1[nH]ccc12","CCc1[nH]ccc1",true,true},
      matchCase{"c2cccc1[nH]ccc12","Cc1[nH]ccc1",false,true},
      matchCase{"c2cccc1[nH]ccc12","CCc1[nH]ccc1",false,false},
    };
    // clang-format on
    for (const auto &tpl : examples) {
      auto smi = std::get<1>(tpl);
      std::unique_ptr<RWMol> qry(SmilesToMol(smi));
      REQUIRE(qry);
      if (std::get<2>(tpl)) {
        CHECK_THAT(std::get<0>(tpl), IsSubstructOf(*qry, smi));
      } else {
        CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(*qry, smi));
      }
      {
        RWMol cp(*qry);
        MolOps::adjustQueryProperties(cp, &ps);
        if (std::get<3>(tpl)) {
          CHECK_THAT(std::get<0>(tpl), IsSubstructOf(cp, smi));
        } else {
          CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(cp, smi));
        }
      }
      {  // make sure ring-finding happens
        RWMol cp(*qry);
        cp.getRingInfo()->reset();
        MolOps::adjustQueryProperties(cp, &ps);
        if (std::get<3>(tpl)) {
          CHECK_THAT(std::get<0>(tpl), IsSubstructOf(cp, smi));
        } else {
          CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(cp, smi));
        }
      }
    }
  }
}

TEST_CASE("single bonds to aromatic neighbors") {
  MolOps::AdjustQueryParameters ps =
      MolOps::AdjustQueryParameters::noAdjustments();
  ps.adjustSingleBondsBetweenAromaticAtoms = true;
  SECTION("matching") {
    // clang-format off
    std::vector<matchCase> examples = {
      matchCase{"c1ncccc1-c1cnncc1","c1ncccc1-c1cnncc1",true,true},
      matchCase{"C1=CC2=C(C=CC3=C2C=NN=C3)N=C1","c1ncccc1-c1cnncc1",false,true},
      matchCase{"C1CC2=C(C=CC=N2)C2=C1C=NN=C2","c1ncccc1-c1cnncc1",true,true},
      matchCase{"C1CC2=NN=CC3=C2C2=C(C=C3)N=CC=C12","c1ncccc1-c1cnncc1",false,true},
      // confirm that we don't modify ring bonds
      matchCase{"C1=CC2=C(C=CC3=C2C=NN=C3)N=C1","C1CC2=C(C=CC=N2)C2=C1C=NN=C2",false,false},
      matchCase{"C1CC2=C(C=CC=N2)C2=C1C=NN=C2","C1CC2=C(C=CC=N2)C2=C1C=NN=C2",true,true},
      matchCase{"C1CC2=C3C(C=CC4=NN=CC1=C34)=CC=N2","C1CC2=C(C=CC=N2)C2=C1C=NN=C2",false,true}, // was github #3325
    };
    // clang-format on
    for (const auto &tpl : examples) {
      auto smi = std::get<1>(tpl);
      std::unique_ptr<RWMol> qry(SmilesToMol(smi));
      REQUIRE(qry);
      if (std::get<2>(tpl)) {
        CHECK_THAT(std::get<0>(tpl), IsSubstructOf(*qry, smi));
      } else {
        CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(*qry, smi));
      }
      MolOps::adjustQueryProperties(*qry, &ps);
      if (std::get<3>(tpl)) {
        CHECK_THAT(std::get<0>(tpl), IsSubstructOf(*qry, smi));
      } else {
        CHECK_THAT(std::get<0>(tpl), !IsSubstructOf(*qry, smi));
      }
    }
  }
}

TEST_CASE(
    "github #3388: Information about charges and isotopes lost when calling "
    "AdjustQueryProperties") {
  MolOps::AdjustQueryParameters ps =
      MolOps::AdjustQueryParameters::noAdjustments();
  ps.adjustDegree = true;
  ps.adjustDegreeFlags = MolOps::AdjustQueryWhichFlags::ADJUST_IGNORENONE;
  SECTION("basics") {
    auto mol = "[13CH3]C[O-]"_smiles;
    REQUIRE(mol);
    MolOps::adjustQueryProperties(*mol, &ps);
    auto sma = MolToSmarts(*mol);
    CHECK(sma == "[#6&13*&D1]-[#6&D2]-[#8&-&D1]");
  }
  SECTION("root cause") {
    auto mol = "[13CH2-]C[O-]"_smiles;
    REQUIRE(mol);
    QueryAtom atm(*mol->getAtomWithIdx(0));
    auto sma = SmartsWrite::GetAtomSmarts(&atm);
    CHECK(sma == "[#6&13*&-]");
  }
  SECTION("root cause2") {
    // since we don't have a way to query for number of radical electrons in
    // SMARTS, we need to check that a different way:
    auto mol = "[CH2]C[O-]"_smiles;
    REQUIRE(mol);
    QueryAtom atm(*mol->getAtomWithIdx(0));
    auto descr = describeQuery(&atm);
    CHECK(descr.find("AtomNumRadicalElectrons 1 = val") != std::string::npos);
  }
}

TEST_CASE("makeAtomsGeneric") {
  SECTION("ignore mapped atoms") {
    MolOps::AdjustQueryParameters ps =
        MolOps::AdjustQueryParameters::noAdjustments();
    ps.makeAtomsGeneric = true;
    ps.makeAtomsGenericFlags =
        MolOps::AdjustQueryWhichFlags::ADJUST_IGNOREMAPPED;
    auto m = "C[CH3:1]"_smiles;
    REQUIRE(m);
    MolOps::adjustQueryProperties(*m, &ps);
    CHECK(MolToSmarts(*m) == "*-[#6H3:1]");
  }
}

TEST_CASE("other edges") {
  SECTION("adjustHeavyDegree") {
    MolOps::AdjustQueryParameters ps =
        MolOps::AdjustQueryParameters::noAdjustments();
    ps.adjustHeavyDegree = true;
    ps.adjustHeavyDegreeFlags = MolOps::ADJUST_IGNOREMAPPED;
    auto m = "C[CH3:1]"_smiles;
    MolOps::adjustQueryProperties(*m, &ps);
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(describeQuery(m->getAtomWithIdx(0)).find("AtomHeavyAtomDegree") !=
          std::string::npos);
    CHECK(!m->getAtomWithIdx(1)->hasQuery());
  }
  SECTION("adjustRingCount") {
    MolOps::AdjustQueryParameters ps =
        MolOps::AdjustQueryParameters::noAdjustments();
    ps.adjustRingCount = true;
    ps.adjustRingCountFlags = MolOps::ADJUST_IGNOREMAPPED;
    auto m = "C[CH3:1]"_smiles;
    MolOps::adjustQueryProperties(*m, &ps);
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(describeQuery(m->getAtomWithIdx(0)).find("AtomInNRings") !=
          std::string::npos);
    CHECK(!m->getAtomWithIdx(1)->hasQuery());
  }
  SECTION("adjustRingCount") {
    MolOps::AdjustQueryParameters ps =
        MolOps::AdjustQueryParameters::noAdjustments();
    ps.adjustRingChain = true;
    ps.adjustRingChainFlags = MolOps::ADJUST_IGNOREMAPPED;
    auto m = "C[CH3:1]"_smiles;
    MolOps::adjustQueryProperties(*m, &ps);
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(describeQuery(m->getAtomWithIdx(0)).find("AtomInRing") !=
          std::string::npos);
    CHECK(!m->getAtomWithIdx(1)->hasQuery());
  }
}

TEST_CASE(
    "Github #4789: AdjustQueryProperties() is inconsistent when "
    "adjustConjugatedFiveRings is set") {
  SECTION("degree one neighbors") {
    auto mol = "C1=CC=CC1C"_smiles;
    REQUIRE(mol);
    MolOps::AdjustQueryParameters ps =
        MolOps::AdjustQueryParameters::noAdjustments();
    ps.adjustConjugatedFiveRings = true;
    ps.adjustSingleBondsToDegreeOneNeighbors = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    CHECK(MolToSmarts(*mol) ==
          "[#6]1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:1[#6]");
  }
  SECTION("between aromatic atoms") {
    auto mol = "C1=CC=CC1C1=CC=CC1"_smiles;
    REQUIRE(mol);
    MolOps::AdjustQueryParameters ps =
        MolOps::AdjustQueryParameters::noAdjustments();
    ps.adjustConjugatedFiveRings = true;
    ps.adjustSingleBondsBetweenAromaticAtoms = true;
    MolOps::adjustQueryProperties(*mol, &ps);
    // clang-format off
    CHECK(MolToSmarts(*mol) ==
          "[#6]1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:1[#6]1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:1");
    // clang-format on
  }
}