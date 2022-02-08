//
//  Copyright (C) 2016-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <string>
#include <GraphMol/RDKitBase.h>
#include "SmilesParse.h"
#include "SmilesWrite.h"
#include "SmartsWrite.h"
#include <RDGeneral/RDLog.h>
using namespace RDKit;

void testBase() {
  BOOST_LOG(rdInfoLog) << "testing base functionality" << std::endl;
  {  // it works when nothing is provided
    std::string smiles = "CC";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
void testCoords2D() {
  BOOST_LOG(rdInfoLog) << "testing reading 2D coordinates" << std::endl;
  {
    std::string smiles = "CC |(0,.75,;0,-.75,)|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumConformers() == 1);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(0).x) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(0).y - 0.75) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(0).z) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(1).x) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(1).y + 0.75) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(1).z) < 1e-4);

    delete m;
  }
  {
    std::string smiles = "CC |(,,;,,-.75)|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getNumConformers() == 1);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(0).x) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(0).y) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(0).z) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(1).x) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(1).y) < 1e-4);
    TEST_ASSERT(fabs(m->getConformer().getAtomPos(1).z + 0.75) < 1e-4);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAtomLabels() {
  BOOST_LOG(rdInfoLog) << "testing reading Atom Labels" << std::endl;
  {
    std::string smiles = "CCC |$foo;;bar$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::atomLabel) == "foo");
    TEST_ASSERT(m->getAtomWithIdx(2)->getProp<std::string>(
                    common_properties::atomLabel) == "bar");
    TEST_ASSERT(!m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    delete m;
  }
  {  // attachment points, example from the docs
    std::string smiles = "C[C@H](N*)C(*)=O |$;;;_AP1;;_AP2;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 7);
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomicNum() == 0);
    TEST_ASSERT(m->getAtomWithIdx(3)->getProp<std::string>(
                    common_properties::atomLabel) == "_AP1");
    // we used to set an atom map for attachment points. This was github #3393
    TEST_ASSERT(m->getAtomWithIdx(3)->getAtomMapNum() == 0);

    TEST_ASSERT(m->getAtomWithIdx(5)->getAtomicNum() == 0);
    TEST_ASSERT(m->getAtomWithIdx(5)->getProp<std::string>(
                    common_properties::atomLabel) == "_AP2");
    // we used to set an atom map for attachment points. This was github #3393
    TEST_ASSERT(m->getAtomWithIdx(5)->getAtomMapNum() == 0);

    delete m;
  }
  {  // query properties
    std::string smiles = "**C |$Q_e;QH_p;;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::atomLabel) == "Q_e");
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::atomLabel) == "QH_p");
    TEST_ASSERT(!m->getAtomWithIdx(2)->hasProp(common_properties::atomLabel));
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(!m->getAtomWithIdx(2)->hasQuery());

    delete m;
  }
  {  // query properties2
    std::string smiles = "** |$;AH_p;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::atomLabel) == "AH_p");
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAtomicNum");
    TEST_ASSERT(m->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(1)->getQuery()->getDescription() ==
                "AtomNull");

    delete m;
  }

  {  // query properties3
    std::string smiles = "** |$;XH_p;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::atomLabel) == "XH_p");
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(0)->getQuery()->getDescription() ==
                "AtomAtomicNum");
    TEST_ASSERT(m->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(1)->getQuery()->getDescription() == "AtomOr");

    delete m;
  }
  {  // query properties3
    std::string smiles = "** |$MH_p;M_p;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::atomLabel) == "MH_p");
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::atomLabel) == "M_p");
    TEST_ASSERT(m->getAtomWithIdx(0)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(0)->getQuery()->getDescription() == "AtomOr");
    TEST_ASSERT(m->getAtomWithIdx(1)->hasQuery());
    TEST_ASSERT(m->getAtomWithIdx(1)->getQuery()->getDescription() == "AtomOr");

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testCXSmilesAndName() {
  BOOST_LOG(rdInfoLog) << "testing CSXMILES and mol name" << std::endl;
  {
    std::string smiles = "CCC |$foo;;bar$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    params.parseName = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::atomLabel) == "foo");
    TEST_ASSERT(m->getProp<std::string>("_CXSMILES_Data") == "|$foo;;bar$|");
    TEST_ASSERT(!m->hasProp("_Name"));
    delete m;
  }
  {
    std::string smiles = "CCC |$foo;;bar$| ourname";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    params.parseName = true;

    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::atomLabel) == "foo");
    TEST_ASSERT(m->getProp<std::string>("_CXSMILES_Data") == "|$foo;;bar$|");
    TEST_ASSERT(m->getProp<std::string>(common_properties::_Name) == "ourname");
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testCoordinateBonds() {
  BOOST_LOG(rdInfoLog) << "testing coordinate bonds" << std::endl;
  {
    std::string smiles = "[Fe]1C=C1 |C:1.0,2.2|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 1);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2)->getBeginAtomIdx() == 2);
    delete m;
  }
  {
    std::string smiles = "C1[Fe]C=1 |C:0.0,2.1|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 2)->getBondType() == Bond::DOUBLE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 0);
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2));
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondBetweenAtoms(1, 2)->getBeginAtomIdx() == 2);
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testRadicals() {
  BOOST_LOG(rdInfoLog) << "testing radicals" << std::endl;
  {
    std::string smiles = "[O]C[O] |^1:0,2|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);

    delete m;
  }
  {
    std::string smiles = "[O][C][O] |^1:0,2,^4:1|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 2);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);

    delete m;
  }
  {  // radicals and coordinate bonds
    std::string smiles = "[Fe]N([O])[O] |^1:2,3,C:1.0|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 4);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 0);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    TEST_ASSERT(m->getAtomWithIdx(3)->getNumRadicalElectrons() == 1);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DATIVE);
    TEST_ASSERT(m->getBondBetweenAtoms(0, 1)->getBeginAtomIdx() == 1);

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAtomValues() {
  BOOST_LOG(rdInfoLog) << "testing atom values" << std::endl;
  {  // testing atom values
    std::string smiles =
        "CCC1=CC=CC=C1 |$_AV:value 2;&#59;value1;value "
        "5&#59;6;;;;;$,c:4,6,t:2|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::molFileValue) == "value 2");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::molFileValue) == ";value1");

    TEST_ASSERT(m->getAtomWithIdx(2)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(2)->getProp<std::string>(
                    common_properties::molFileValue) == "value 5;6");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testAtomProps() {
  BOOST_LOG(rdInfoLog) << "testing atom properties" << std::endl;
  {  // testing atom properties
    std::string smiles =
        "C1CN1 "
        "|atomProp:0.prop2.val2:0.prop1.val1:1.prop2.v2&#38;4:1.prop1.v1;2;3|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("prop1"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>("prop1") == "val1");
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("prop2"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>("prop2") == "val2");
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("prop2"));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>("prop2") == "v2&4");
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("prop1"));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>("prop1") ==
                "v1;2;3");

    delete m;
  }

  {  // testing atom properties + values
    std::string smiles =
        "C1CN1 "
        "|atomProp:0.prop2.val2:1.prop1.v1;2;3,$_AV:value 2;&#59;value1;$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp("prop2"));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>("prop2") == "val2");
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp("prop1"));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>("prop1") ==
                "v1;2;3");
    TEST_ASSERT(m->getAtomWithIdx(0)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(0)->getProp<std::string>(
                    common_properties::molFileValue) == "value 2");

    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::molFileValue));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::molFileValue) == ";value1");

    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1968() {
  BOOST_LOG(rdInfoLog)
      << "testing Github1968: CXSMILES should be parsed before H removal"
      << std::endl;
  {  // the original report
    std::string smiles = "[H]C* |$;;X$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::atomLabel) == "X");
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    delete m;
  }
  {
    std::string smiles = "C([H])* |$;Y;X$|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);
    TEST_ASSERT(m->getAtomWithIdx(1)->hasProp(common_properties::atomLabel));
    TEST_ASSERT(m->getAtomWithIdx(1)->getProp<std::string>(
                    common_properties::atomLabel) == "X");
    TEST_ASSERT(!m->getAtomWithIdx(0)->hasProp(common_properties::atomLabel));
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testEnhancedStereo() {
  BOOST_LOG(rdInfoLog) << "Testing CXSMILES Enhanced Stereo" << std::endl;

  std::vector<unsigned int> atom_ref1({4, 5});
  {
    std::string smiles = "C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    auto &stereo_groups = m->getStereoGroups();

    TEST_ASSERT(stereo_groups.size() == 2);

    auto stg = stereo_groups.begin();
    TEST_ASSERT(stg->getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    {
      auto &atoms = stg->getAtoms();
      TEST_ASSERT(atoms.size() == 1);
      TEST_ASSERT(atoms[0]->getIdx() == 1);
    }
    ++stg;
    TEST_ASSERT(stg->getGroupType() == StereoGroupType::STEREO_OR);
    {
      auto &atoms = stg->getAtoms();
      TEST_ASSERT(atoms.size() == 2);
      TEST_ASSERT(atoms[0]->getIdx() == 3);
      TEST_ASSERT(atoms[1]->getIdx() == 5);
    }
    delete m;
  }
  {
    std::string smiles = "C[C@H](F)[C@H](C)[C@@H](C)Br |&1:3,5,a:1|";
    SmilesParserParams params;
    params.allowCXSMILES = true;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 8);

    auto &stereo_groups = m->getStereoGroups();

    TEST_ASSERT(stereo_groups.size() == 2);

    auto stg = stereo_groups.begin();
    TEST_ASSERT(stg->getGroupType() == StereoGroupType::STEREO_AND);
    {
      auto &atoms = stg->getAtoms();
      TEST_ASSERT(atoms.size() == 2);
      TEST_ASSERT(atoms[0]->getIdx() == 3);
      TEST_ASSERT(atoms[1]->getIdx() == 5);
    }
    ++stg;
    TEST_ASSERT(stg->getGroupType() == StereoGroupType::STEREO_ABSOLUTE);
    {
      auto &atoms = stg->getAtoms();
      TEST_ASSERT(atoms.size() == 1);
      TEST_ASSERT(atoms[0]->getIdx() == 1);
    }
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testHTMLCharCodes() {
  BOOST_LOG(rdInfoLog) << "Testing CXSMILES with HTML char codes" << std::endl;

  {
    std::string smiles = R"(CCCC* |$;;;;_AP1$,Sg:n:2:2&#44;6-7:ht|)";
    SmilesParserParams params;
    params.allowCXSMILES = true;

    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);

    TEST_ASSERT(m->getNumAtoms() == 5);

    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testErrorsInCXSmiles() {
  BOOST_LOG(rdInfoLog) << "Testing handling of errors in CXSMILES" << std::endl;

  {
    std::string smiles = R"(CC |failure)";
    SmilesParserParams params;

    ROMol *m = nullptr;
    try {
      m = SmilesToMol(smiles, params);
    } catch (const SmilesParseException &) {
    }
    TEST_ASSERT(!m);

    params.strictCXSMILES = false;
    m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 2);

    delete m;
  }

  {  // sure partial parsing also works
    std::string smiles = "[O][C][O] |^1:0,2,^4:1 FAILURE|";
    SmilesParserParams params;
    params.strictCXSMILES = false;
    ROMol *m = SmilesToMol(smiles, params);
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumAtoms() == 3);
    TEST_ASSERT(m->getAtomWithIdx(0)->getNumRadicalElectrons() == 1);
    TEST_ASSERT(m->getAtomWithIdx(1)->getNumRadicalElectrons() == 2);
    TEST_ASSERT(m->getAtomWithIdx(2)->getNumRadicalElectrons() == 1);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testMDLQueriesInCXSmiles() {
  BOOST_LOG(rdInfoLog) << "Testing handling of MDL queries in CXSMILES"
                       << std::endl;

  {  // just ring bonds
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|rb:0:0,1:*,2:2|"_smiles;
    TEST_ASSERT(m);
    auto sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6&x0]-[#6&x2]1-[#6&x2]-[#6]-1");
  }

  {  // unsaturation
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|u:0,2|"_smiles;
    TEST_ASSERT(m);
    auto sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6&$(*=,:,#*)]-[#6]1-[#6&$(*=,:,#*)]-[#6]-1");
  }

  {  // substitution count
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|s:3:*,1:3|"_smiles;
    TEST_ASSERT(m);
    auto sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6]-[#6&d3]1-[#6]-[#6&d2]-1");
  }

  {  // everything together
    auto m =
        "[#6]-[#6]-1-[#6]-[#6]-1 "
        "|rb:0:0,1:*,2:2,s:3:*,u:0|"_smiles;
    TEST_ASSERT(m);
    auto sma = MolToSmarts(*m);
    TEST_ASSERT(sma == "[#6&x0&$(*=,:,#*)]-[#6&x2]1-[#6&x2]-[#6&d2]-1");
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testLinknodesInCXSmiles() {
  BOOST_LOG(rdInfoLog) << "Testing handling of LINKNODES in CXSMILES"
                       << std::endl;

  {
    auto m = "OC1CCC(F)C1 |LN:1:1.3.2.6|"_smiles;
    TEST_ASSERT(m);
    std::string lns;
    TEST_ASSERT(m->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    TEST_ASSERT(lns == "1 3 2 2 3 2 7")
  }
  {
    auto m = "OC1CCC(F)C1 |LN:1:1.3.2.6,4:1.4.3.6|"_smiles;
    TEST_ASSERT(m);
    std::string lns;
    TEST_ASSERT(m->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    TEST_ASSERT(lns == "1 3 2 2 3 2 7|1 4 2 5 4 5 7")
  }
  {  // linknodes with implicit outer atoms
    auto m = "OC1CCCC1 |LN:4:1.3|"_smiles;
    TEST_ASSERT(m);
    std::string lns;
    TEST_ASSERT(m->getPropIfPresent(common_properties::molFileLinkNodes, lns));
    TEST_ASSERT(lns == "1 3 2 5 4 5 6")
  }
  {  // linknodes with implicit outer atoms : fails because atom is
     // three-connected
    bool ok = false;
    try {
      auto m = "OC1CCC(F)C1 |LN:1:1.3|"_smiles;
    } catch (const RDKit::SmilesParseException &) {
      ok = true;
    }
    TEST_ASSERT(ok);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testVariableAttachmentInCXSmiles() {
  BOOST_LOG(rdInfoLog)
      << "Testing handling of variable attachment bonds in CXSMILES"
      << std::endl;

  {
    auto m = "CO*.C1=CC=NC=C1 |m:2:3.5.4|"_smiles;
    TEST_ASSERT(m);
    const auto bnd = m->getBondBetweenAtoms(1, 2);
    TEST_ASSERT(bnd);
    TEST_ASSERT(bnd->hasProp(common_properties::_MolFileBondAttach));
    TEST_ASSERT(bnd->getProp<std::string>(
                    common_properties::_MolFileBondAttach) == "ANY");
    TEST_ASSERT(bnd->hasProp(common_properties::_MolFileBondEndPts));
    TEST_ASSERT(bnd->getProp<std::string>(
                    common_properties::_MolFileBondEndPts) == "(3 4 6 5)");
  }

  {
    auto m = "F*.Cl*.C1=CC=NC=C1 |m:1:9.8,3:4.5|"_smiles;
    TEST_ASSERT(m);
    {
      const auto bnd = m->getBondBetweenAtoms(0, 1);
      TEST_ASSERT(bnd);
      TEST_ASSERT(bnd->hasProp(common_properties::_MolFileBondAttach));
      TEST_ASSERT(bnd->getProp<std::string>(
                      common_properties::_MolFileBondAttach) == "ANY");
      TEST_ASSERT(bnd->hasProp(common_properties::_MolFileBondEndPts));
      TEST_ASSERT(bnd->getProp<std::string>(
                      common_properties::_MolFileBondEndPts) == "(2 10 9)");
    }
    {
      const auto bnd = m->getBondBetweenAtoms(2, 3);
      TEST_ASSERT(bnd);
      TEST_ASSERT(bnd->hasProp(common_properties::_MolFileBondAttach));
      TEST_ASSERT(bnd->getProp<std::string>(
                      common_properties::_MolFileBondAttach) == "ANY");
      TEST_ASSERT(bnd->hasProp(common_properties::_MolFileBondEndPts));
      TEST_ASSERT(bnd->getProp<std::string>(
                      common_properties::_MolFileBondEndPts) == "(2 5 6)");
    }
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
#if 1
  testBase();
  testCoords2D();
  testAtomLabels();
  testCXSmilesAndName();
  testCoordinateBonds();
  testRadicals();
  testAtomValues();
#endif
  testAtomProps();
  testGithub1968();
  testEnhancedStereo();
  testHTMLCharCodes();
  testErrorsInCXSmiles();
  testMDLQueriesInCXSmiles();
  testLinknodesInCXSmiles();
  testVariableAttachmentInCXSmiles();
}
