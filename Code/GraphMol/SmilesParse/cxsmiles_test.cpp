//
//  Copyright (C) 2016-2023 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "catch.hpp"
#include <RDGeneral/test.h>
#include <string>
#include <GraphMol/RDKitBase.h>
#include "SmilesParse.h"
#include "SmilesWrite.h"
#include "SmartsWrite.h"
#include <RDGeneral/RDLog.h>
using namespace RDKit;

TEST_CASE("base functionality") {
  {  // it works when nothing is provided
    std::string smiles = "CC";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    delete m;
  }
}
TEST_CASE("reading 2D coordinates") {
  {
    std::string smiles = "CC |(0,.75,;0,-.75,)|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getNumConformers() == 1);
    CHECK(fabs(m->getConformer().getAtomPos(0).x) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(0).y - 0.75) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(0).z) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(1).x) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(1).y + 0.75) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(1).z) < 1e-4);

    delete m;
  }
  {
    std::string smiles = "CC |(,,;,,-.75)|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getNumConformers() == 1);
    CHECK(fabs(m->getConformer().getAtomPos(0).x) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(0).y) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(0).z) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(1).x) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(1).y) < 1e-4);
    CHECK(fabs(m->getConformer().getAtomPos(1).z + 0.75) < 1e-4);

    delete m;
  }
}

TEST_CASE("reading Atom Labels") {
  {
    std::string smiles = "CCC |$foo;;bar$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "foo");
    CHECK(m->getAtomWithIdx(2)->getProp<std::string>(
              common_properties::atomLabel) == "bar");
    CHECK(!m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    delete m;
  }
  {  // attachment points, example from the docs
    std::string smiles = "C[C@H](N*)C(*)=O |$;;;_AP1;;_AP2;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 7);
    CHECK(m->getAtomWithIdx(3)->getAtomicNum() == 0);
    CHECK(m->getAtomWithIdx(3)->getProp<std::string>(
              common_properties::atomLabel) == "_AP1");
    // we used to set an atom map for attachment points. This was github #3393
    CHECK(m->getAtomWithIdx(3)->getAtomMapNum() == 0);

    CHECK(m->getAtomWithIdx(5)->getAtomicNum() == 0);
    CHECK(m->getAtomWithIdx(5)->getProp<std::string>(
              common_properties::atomLabel) == "_AP2");
    // we used to set an atom map for attachment points. This was github #3393
    CHECK(m->getAtomWithIdx(5)->getAtomMapNum() == 0);

    delete m;
  }
  {  // query properties
    std::string smiles = "**C |$Q_e;QH_p;;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "Q_e");
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "QH_p");
    CHECK(!m->getAtomWithIdx(2)->hasProp(common_properties::atomLabel));
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(m->getAtomWithIdx(1)->hasQuery());
    CHECK(!m->getAtomWithIdx(2)->hasQuery());

    delete m;
  }
  {  // query properties2
    std::string smiles = "** |$;AH_p;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "AH_p");
    CHECK(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(m->getAtomWithIdx(0)->getQuery()->getDescription() ==
          "AtomAtomicNum");
    CHECK(m->getAtomWithIdx(1)->hasQuery());
    CHECK(m->getAtomWithIdx(1)->getQuery()->getDescription() == "AtomNull");

    delete m;
  }

  {  // query properties3
    std::string smiles = "** |$;XH_p;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "XH_p");
    CHECK(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(m->getAtomWithIdx(0)->getQuery()->getDescription() ==
          "AtomAtomicNum");
    CHECK(m->getAtomWithIdx(1)->hasQuery());
    CHECK(m->getAtomWithIdx(1)->getQuery()->getDescription() == "AtomOr");

    delete m;
  }
  {  // query properties3
    std::string smiles = "** |$MH_p;M_p;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "MH_p");
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "M_p");
    CHECK(m->getAtomWithIdx(0)->hasQuery());
    CHECK(m->getAtomWithIdx(0)->getQuery()->getDescription() == "AtomOr");
    CHECK(m->getAtomWithIdx(1)->hasQuery());
    CHECK(m->getAtomWithIdx(1)->getQuery()->getDescription() == "AtomOr");

    delete m;
  }
}

TEST_CASE("CXSMILES and mol name") {
  {
    std::string smiles = "CCC |$foo;;bar$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    params.parseName = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "foo");
    CHECK(m->getProp<std::string>("_CXSMILES_Data") == "|$foo;;bar$|");
    CHECK(!m->hasProp("_Name"));
    delete m;
  }
  {
    std::string smiles = "CCC |$foo;;bar$| ourname";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    params.parseName = true;

    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::atomLabel) == "foo");
    CHECK(m->getProp<std::string>("_CXSMILES_Data") == "|$foo;;bar$|");
    CHECK(m->getProp<std::string>(common_properties::_Name) == "ourname");
    delete m;
  }
}

TEST_CASE("coordinate bonds") {
  {
    std::string smiles = "[Fe]1C=C1 |C:1.0,2.2|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getBondBetweenAtoms(1, 2));
    CHECK(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    CHECK(m->getBondBetweenAtoms(0, 1));
    CHECK(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    CHECK(m->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 1);
    CHECK(m->getBondBetweenAtoms(0, 2));
    CHECK(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::DATIVE);
    CHECK(m->getBondBetweenAtoms(0, 2)->getBeginAtomIdx() == 2);
    delete m;
  }
  {
    std::string smiles = "C1[Fe]C=1 |C:0.0,2.1|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getBondBetweenAtoms(0, 2));
    CHECK(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::DOUBLE);
    CHECK(m->getBondBetweenAtoms(0, 1));
    CHECK(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    CHECK(m->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 0);
    CHECK(m->getBondBetweenAtoms(1, 2));
    CHECK(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DATIVE);
    CHECK(m->getBondBetweenAtoms(1, 2)->getBeginAtomIdx() == 2);
    delete m;
  }
}

TEST_CASE("radicals") {
  {
    std::string smiles = "[O]C[O] |^1:0,2|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    CHECK(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    CHECK(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);

    delete m;
  }
  {
    std::string smiles = "[O][C][O] |^1:0,2,^4:1|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    CHECK(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 2);
    CHECK(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);

    delete m;
  }
  {  // radicals and coordinate bonds
    std::string smiles = "[Fe]N([O])[O] |^1:2,3,C:1.0|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 4);
    CHECK(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    CHECK(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    CHECK(m->getAtomWithIdx(3)->getNumRadicalElectrons() == 1);
    CHECK(m->getBondBetweenAtoms(0, 1));
    CHECK(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    CHECK(m->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 1);

    delete m;
  }
}

TEST_CASE("atom values") {
  {  // testing atom values
    std::string smiles =
        "CCC1=CC=CC=C1 |$_AV:value 2;&#59;value1;value "
        "5&#59;6;;;;;$,c:4,6,t:2|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 8);
    CHECK(m->getAtomWithIdx(0)->hasProp(common_properties::molFileValue));
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::molFileValue) == "value 2");

    CHECK(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::molFileValue) == ";value1");

    CHECK(m->getAtomWithIdx(2)->hasProp(common_properties::molFileValue));
    CHECK(m->getAtomWithIdx(2)->getProp<std::string>(
              common_properties::molFileValue) == "value 5;6");

    delete m;
  }
}

TEST_CASE("atom properties") {
  {  // testing atom properties
    std::string smiles =
        "C1CN1 "
        "|atomProp:0.prop2.val2:0.prop1.val1:1.prop2.v2&#38;4:1.prop1.v1;2;3|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->hasProp("prop1"));
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>("prop1") == "val1");
    CHECK(m->getAtomWithIdx(0)->hasProp("prop2"));
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>("prop2") == "val2");
    CHECK(m->getAtomWithIdx(1)->hasProp("prop2"));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("prop2") == "v2&4");
    CHECK(m->getAtomWithIdx(1)->hasProp("prop1"));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("prop1") == "v1;2;3");

    delete m;
  }

  {  // testing atom properties + values
    std::string smiles =
        "C1CN1 "
        "|atomProp:0.prop2.val2:1.prop1.v1;2;3,$_AV:value 2;&#59;value1;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->hasProp("prop2"));
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>("prop2") == "val2");
    CHECK(m->getAtomWithIdx(1)->hasProp("prop1"));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>("prop1") == "v1;2;3");
    CHECK(m->getAtomWithIdx(0)->hasProp(common_properties::molFileValue));
    CHECK(m->getAtomWithIdx(0)->getProp<std::string>(
              common_properties::molFileValue) == "value 2");

    CHECK(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::molFileValue) == ";value1");

    delete m;
  }
}

TEST_CASE("Github1968: CXSMILES should be parsed before H removal") {
  {  // the original report
    std::string smiles = "[H]C* |$;;X$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "X");
    CHECK(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    delete m;
  }
  {
    std::string smiles = "C([H])* |$;Y;X$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);
    CHECK(m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    CHECK(m->getAtomWithIdx(1)->getProp<std::string>(
              common_properties::atomLabel) == "X");
    CHECK(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    delete m;
  }
}

TEST_CASE("enhanced stereo") {
  std::vector<unsigned int> atom_ref1({4, 5});
  {
    std::string smiles = "C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 8);

    auto &stereo_groups = m->getStereoGroups();

    CHECK(stereo_groups.size() == 2);

    auto stg = stereo_groups.begin();
    CHECK(stg->getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    {
      auto &atoms = stg->getAtoms();
      CHECK(atoms.size() == 1);
      CHECK(atoms[0]->getIdx() == 1);
    }
    ++stg;
    CHECK(stg->getGroupType() == StereoGroupType::STEREO_OR);
    {
      auto &atoms = stg->getAtoms();
      CHECK(atoms.size() == 2);
      CHECK(atoms[0]->getIdx() == 3);
      CHECK(atoms[1]->getIdx() == 5);
    }
    delete m;
  }
  {
    std::string smiles = "C[C@H](F)[C@H](C)[C@@H](C)Br |&1:3,5,a:1|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 8);

    auto &stereo_groups = m->getStereoGroups();

    CHECK(stereo_groups.size() == 2);

    auto stg = stereo_groups.begin();
    CHECK(stg->getGroupType() == StereoGroupType::STEREO_AND);
    {
      auto &atoms = stg->getAtoms();
      CHECK(atoms.size() == 2);
      CHECK(atoms[0]->getIdx() == 3);
      CHECK(atoms[1]->getIdx() == 5);
    }
    ++stg;
    CHECK(stg->getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    {
      auto &atoms = stg->getAtoms();
      CHECK(atoms.size() == 1);
      CHECK(atoms[0]->getIdx() == 1);
    }
    delete m;
  }
}

TEST_CASE("HTML char codes") {
  {
    std::string smiles = R"(CCCC* |$;;;;_AP1$,Sg:n:2:2&#44;6-7:ht|)";
    SmilesParserParams params;
    params.allowCXSMILES = true;

    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);

    CHECK(m->getNumAtoms() == 5);

    delete m;
  }
}

TEST_CASE("errors in CXSMILES") {
  {
    std::string smiles = R"(CC |failure)";
    SmilesParserParams params;

    ROMol *m = nullptr;
    try {
      m = SmilesToMol(smiles, params);
    } catch (const SmilesParseException &) {
    }
    CHECK(!m);

    params.strictCXSMILES = false;
    m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 2);

    delete m;
  }

  {  // sure partial parsing also works
    std::string smiles = "[O][C][O] |^1:0,2,^4:1 FAILURE|";
    SmilesParserParams params;
    params.strictCXSMILES = false;
    ROMol *m = SmilesToMol(smiles, params);
    REQUIRE(m);
    CHECK(m->getNumAtoms() == 3);
    CHECK(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    CHECK(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 2);
    CHECK(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    delete m;
  }
}

TEST_CASE("MDL queries") {
  {  // just ring bonds
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|rb:0:0,1:*,2:2|"_smiles;
    REQUIRE(m);
    auto sma = MolToSmarts(*m);
    CHECK(sma == "[#6&x0]-[#6&x2]1-[#6&x2]-[#6]-1");
  }

  {  // unsaturation
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|u:0,2|"_smiles;
    REQUIRE(m);
    auto sma = MolToSmarts(*m);
    CHECK(sma == "[#6&$(*=,:,#*)]-[#6]1-[#6&$(*=,:,#*)]-[#6]-1");
  }

  {  // substitution count
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|s:3:*,1:3|"_smiles;
    REQUIRE(m);
    auto sma = MolToSmarts(*m);
    CHECK(sma == "[#6]-[#6&d3]1-[#6]-[#6&d2]-1");
  }

  {  // everything together
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|rb:0:0,1:*,2:2,s:3:*,u:0|"_smiles;
    REQUIRE(m);
    auto sma = MolToSmarts(*m);
    CHECK(sma == "[#6&x0&$(*=,:,#*)]-[#6&x2]1-[#6&x2]-[#6&d2]-1");
  }
}

TEST_CASE("LINKNODES") {
  {
    auto m = "OC1CCC(F)C1 |LN:1:1.3.2.6|"_smiles;
    REQUIRE(m);
    std::string lns;
    CHECK(m->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    CHECK(lns == "1 3 2 2 3 2 7");
  }
  {
    auto m = "OC1CCC(F)C1 |LN:1:1.3.2.6,4:1.4.3.6|"_smiles;
    REQUIRE(m);
    std::string lns;
    CHECK(m->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    CHECK(lns == "1 3 2 2 3 2 7|1 4 2 5 4 5 7");
  }
  {  // linknodes with implicit outer atoms
    auto m = "OC1CCCC1 |LN:4:1.3|"_smiles;
    REQUIRE(m);
    std::string lns;
    CHECK(m->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    CHECK(lns == "1 3 2 5 4 5 6");
  }
  {  // linknodes with implicit outer atoms : fails because atom is
     // three-connected
    bool ok = false;
    try {
      auto m = "OC1CCC(F)C1 |LN:1:1.3|"_smiles;
    } catch (const RDKit::SmilesParseException &) {
      ok = true;
    }
    CHECK(ok);
  }
}

TEST_CASE("variable attachment bonds") {
  {
    auto m = "CO*.C1=CC=NC=C1 |m:2:3.5.4|"_smiles;
    REQUIRE(m);
    const auto bnd = m->getBondBetweenAtoms(1, 2);
    REQUIRE(bnd);
    CHECK(bnd->hasProp(common_properties::_MolFileBondAttach));
    CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondAttach) ==
          "ANY");
    CHECK(bnd->hasProp(common_properties::_MolFileBondEndPts));
    CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondEndPts) ==
          "(3 4 6 5)");
  }

  {
    auto m = "F*.Cl*.C1=CC=NC=C1 |m:1:9.8,3:4.5|"_smiles;
    REQUIRE(m);
    {
      const auto bnd = m->getBondBetweenAtoms(0, 1);
      REQUIRE(bnd);
      CHECK(bnd->hasProp(common_properties::_MolFileBondAttach));
      CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondAttach) ==
            "ANY");
      CHECK(bnd->hasProp(common_properties::_MolFileBondEndPts));
      CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondEndPts) ==
            "(2 10 9)");
    }
    {
      const auto bnd = m->getBondBetweenAtoms(2, 3);
      REQUIRE(bnd);
      CHECK(bnd->hasProp(common_properties::_MolFileBondAttach));
      CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondAttach) ==
            "ANY");
      CHECK(bnd->hasProp(common_properties::_MolFileBondEndPts));
      CHECK(bnd->getProp<std::string>(common_properties::_MolFileBondEndPts) ==
            "(2 5 6)");
    }
  }
}

TEST_CASE("github #6050: stereogroups not combined") {
  SECTION("basics") {
    auto m =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,o1:7,&2:1,&3:5,r|"_smiles;
    REQUIRE(m);
    CHECK(m->getStereoGroups().size() == 3);
  }
  SECTION("duplicate indices") {
    auto m =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,o3:7,&2:1,&3:5,r|"_smiles;
    REQUIRE(m);
    CHECK(m->getStereoGroups().size() == 3);
    CHECK(m->getStereoGroups()[0].getGroupType() ==
          StereoGroupType::STEREO_AND);
    CHECK(m->getStereoGroups()[0].getAtoms().size() == 2);
    CHECK(m->getStereoGroups()[1].getGroupType() == StereoGroupType::STEREO_OR);
    CHECK(m->getStereoGroups()[1].getAtoms().size() == 1);
  }
  SECTION("multiple repeats") {
    auto m =
        "C[C@@H](O)[C@H](C)[C@@H](C)[C@H](O)[C@@H](C)[C@H](C)[C@H](C)O |o1:7,&2:1,&3:3,5,o1:9,11,o1:13,r|"_smiles;
    REQUIRE(m);
    CHECK(m->getStereoGroups().size() == 3);
    CHECK(m->getStereoGroups()[0].getGroupType() == StereoGroupType::STEREO_OR);
    CHECK(m->getStereoGroups()[0].getAtoms().size() == 4);
  }
}