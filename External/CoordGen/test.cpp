//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>

#include <fstream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/RDLog.h>

#include "coordgenlibs/sketcherMinimizer.h"
#include "CoordGen.h"

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test1: basics" << std::endl;
#if 1
  {
    ROMol* m = SmilesToMol("c1cc(CC)cnc1CC(=O)O");
    TEST_ASSERT(m);
    m->setProp("_Name", "test1");

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }
  {
    // ROMol* m = SmilesToMol("c1ccncc1");

    ROMol* m = SmilesToMol("ClC(O)(F)C");
    TEST_ASSERT(m);
    m->setProp("_Name", "test2");

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }

  {
    ROMol* m = SmilesToMol(
        "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(="
        "O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)"
        "[C@@H]([C@H](C)C/C=C/"
        "C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C");
    TEST_ASSERT(m);
    m->setProp("_Name", "cyclosporine a");

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }

  {
    // ROMol* m = SmilesToMol("c1ccncc1");

    ROMol* m = SmilesToMol("CCCNC=CNCOC=CC=CC=COC");
    TEST_ASSERT(m);
    m->setProp("_Name", "single-double");

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }
#endif

  {
    ROMol* m = SmilesToMol("O/C=C/C=C/C=C\\C=C/N");
    TEST_ASSERT(m);
    m->setProp("_Name", "cis-trans");

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void test2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test2: using templates" << std::endl;

  {
    ROMol* core = SmilesToMol("C1CONCN1");
    TEST_ASSERT(core);
    core->setProp("_Name", "core");

    CoordGen::addCoords(*core);
    TEST_ASSERT(core->getNumConformers() == 1);
    auto mb = MolToMolBlock(*core);
    std::cerr << mb << std::endl;

    ROMol* m = SmilesToMol("C1ONC(CCCCCCCCCC)NC1");
    TEST_ASSERT(m);
    m->setProp("_Name", "core+sidechain");

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;

    MatchVectType mv;
    SubstructMatch(*m, *core, mv);
    auto coreConf = core->getConformer();
    RDGeom::INT_POINT2D_MAP coordMap;
    for (unsigned int i = 0; i < mv.size(); ++i) {
      coordMap[mv[i].second] =
          RDGeom::Point2D(coreConf.getAtomPos(mv[i].first).x,
                          coreConf.getAtomPos(mv[i].first).y);
    }
    CoordGen::CoordGenParams params;
    params.coordMap = &coordMap;
    CoordGen::addCoords(*m, &params);
    TEST_ASSERT(m->getNumConformers() == 1);
    m->setProp("_Name", "templated");
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;

    delete m;
    delete core;
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  test1();
  test2();
}
