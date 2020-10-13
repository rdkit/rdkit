//
//  Copyright (C) 2017-2020 Greg Landrum
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
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

#include <RDGeneral/RDLog.h>

#include "coordgen/sketcherMinimizer.h"
#include <CoordGen/CoordGen.h>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test1: basics" << std::endl;
#if 1
  {
    ROMol* m = SmilesToMol("c1cc(CC)cnc1CC(=O)O");
    TEST_ASSERT(m);
    m->setProp("_Name", "test1");

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
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

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
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

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
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

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
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

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }

  {
    ROMol* m = SmilesToMol("C1C3CC2CC(CC1C2)C3");
    TEST_ASSERT(m);
    m->setProp("_Name", "admntn");

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
    TEST_ASSERT(m->getNumConformers() == 1);
    auto mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

namespace {
bool compareConfs(const ROMol* m, ROMol* templ, const MatchVectType& mv,
                  bool alignFirst = false, int molConfId = -1,
                  int templateConfId = -1, double postol = 1e-2,
                  double rmstol = 0.1) {
  PRECONDITION(m, "bad pointer");
  PRECONDITION(templ, "bad pointer");
  TEST_ASSERT(m->getNumAtoms() >= templ->getNumAtoms());

  if (alignFirst) {
    double rmsd =
        MolAlign::alignMol(*templ, *m, molConfId, templateConfId, &mv);
    if (rmsd > rmstol) {
      return false;
    }
  }

  const Conformer& conf1 = m->getConformer(molConfId);
  const Conformer& conf2 = templ->getConformer(templateConfId);
  for (unsigned int i = 0; i < templ->getNumAtoms(); i++) {
    TEST_ASSERT(m->getAtomWithIdx(mv[i].second)->getAtomicNum() ==
                templ->getAtomWithIdx(mv[i].first)->getAtomicNum());

    RDGeom::Point3D pt1i = conf1.getAtomPos(mv[i].second);
    RDGeom::Point3D pt2i = conf2.getAtomPos(mv[i].first);
    if ((pt1i - pt2i).length() >= postol) {
      return false;
    }
  }
  return true;
}
}  // namespace

void test2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test2: using templates" << std::endl;

  {
    ROMol* core = SmilesToMol("C1CON1");
    TEST_ASSERT(core);
    core->setProp("_Name", "core");

    TEST_ASSERT(CoordGen::addCoords(*core) == 0);
    TEST_ASSERT(core->getNumConformers() == 1);
    auto mb = MolToMolBlock(*core);
    std::cerr << mb << std::endl;

    ROMol* m = SmilesToMol("C1C(CCC)ON1");
    TEST_ASSERT(m);
    m->setProp("_Name", "core+sidechain");

    MatchVectType mv;
    SubstructMatch(*m, *core, mv);

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
    TEST_ASSERT(m->getNumConformers() == 1);
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    TEST_ASSERT(!compareConfs(m, core, mv));

    {
      auto coreConf = core->getConformer();
      RDGeom::INT_POINT2D_MAP coordMap;
      for (auto& i : mv) {
        coordMap[i.second] = RDGeom::Point2D(coreConf.getAtomPos(i.first).x,
                                             coreConf.getAtomPos(i.first).y);
      }
      CoordGen::CoordGenParams params;
      params.coordMap = coordMap;
      params.dbg_useFixed = true;
      TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
      TEST_ASSERT(m->getNumConformers() == 1);
      // m->setProp("_Name", "templated");
      // mb = MolToMolBlock(*m);
      // std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv));
    }
    {
      CoordGen::CoordGenParams params;
      params.templateMol = core;
      params.dbg_useFixed = true;
      TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
      TEST_ASSERT(m->getNumConformers() == 1);
      m->setProp("_Name", "templated");
      mb = MolToMolBlock(*m);
      std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv));
    }
    delete m;
    delete core;
  }

  {
    ROMol* core = SmilesToMol("C1CCCCCONCN1");
    TEST_ASSERT(core);
    core->setProp("_Name", "core");

    TEST_ASSERT(CoordGen::addCoords(*core) == 0);
    TEST_ASSERT(core->getNumConformers() == 1);
    auto mb = MolToMolBlock(*core);
    std::cerr << mb << std::endl;

    ROMol* m = SmilesToMol("C1CCCCONC(CC)NC1");
    TEST_ASSERT(m);
    m->setProp("_Name", "core+sidechain");

    MatchVectType mv;
    SubstructMatch(*m, *core, mv);

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
    TEST_ASSERT(m->getNumConformers() == 1);
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    TEST_ASSERT(!compareConfs(m, core, mv));

    {
      auto coreConf = core->getConformer();
      RDGeom::INT_POINT2D_MAP coordMap;
      for (auto& i : mv) {
        coordMap[i.second] = RDGeom::Point2D(coreConf.getAtomPos(i.first).x,
                                             coreConf.getAtomPos(i.first).y);
      }
      CoordGen::CoordGenParams params;
      params.coordMap = coordMap;
      params.dbg_useFixed = true;
      TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
      TEST_ASSERT(m->getNumConformers() == 1);
      // m->setProp("_Name", "templated");
      // mb = MolToMolBlock(*m);
      // std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv));
    }
    {
      CoordGen::CoordGenParams params;
      params.templateMol = core;
      params.dbg_useFixed = true;
      TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
      TEST_ASSERT(m->getNumConformers() == 1);
      m->setProp("_Name", "templated");
      mb = MolToMolBlock(*m);
      std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv));
    }

    delete m;
    delete core;
  }

  {
    ROMol* core = SmilesToMol("C1CCCCCONCN1");
    TEST_ASSERT(core);
    core->setProp("_Name", "core");

    TEST_ASSERT(CoordGen::addCoords(*core) == 0);
    TEST_ASSERT(core->getNumConformers() == 1);
    auto mb = MolToMolBlock(*core);
    std::cerr << mb << std::endl;

    ROMol* m = SmilesToMol("C1CCCCONC(CCCCCC)NC1");
    TEST_ASSERT(m);
    m->setProp("_Name", "core+sidechain");

    MatchVectType mv;
    SubstructMatch(*m, *core, mv);

    TEST_ASSERT(CoordGen::addCoords(*m) == 0);
    TEST_ASSERT(m->getNumConformers() == 1);
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    TEST_ASSERT(!compareConfs(m, core, mv));

    {
      auto coreConf = core->getConformer();
      RDGeom::INT_POINT2D_MAP coordMap;
      for (auto& i : mv) {
        coordMap[i.second] = RDGeom::Point2D(coreConf.getAtomPos(i.first).x,
                                             coreConf.getAtomPos(i.first).y);
      }

      CoordGen::CoordGenParams params;
      params.coordMap = coordMap;
      params.dbg_useFixed = true;
      TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
      TEST_ASSERT(m->getNumConformers() == 1);
      // m->setProp("_Name", "templated");
      // mb = MolToMolBlock(*m);
      // std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv));
    }
    {
      CoordGen::CoordGenParams params;
      params.templateMol = core;
      params.dbg_useFixed = true;
      TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
      TEST_ASSERT(m->getNumConformers() == 1);
      m->setProp("_Name", "templated");
      mb = MolToMolBlock(*m);
      std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv));
    }
    delete m;
    delete core;
  }

  {
    ROMol* core = SmilesToMol("C1CCCC2C1NCC2");
    TEST_ASSERT(core);
    core->setProp("_Name", "core");

    CoordGen::addCoords(*core);
    TEST_ASSERT(core->getNumConformers() == 1);
    auto mb = MolToMolBlock(*core);
    std::cerr << mb << std::endl;

    ROMol* m = SmilesToMol("C1C(CCC)CC(CC3CC3)C2C1N(C(C)C)CC2");
    TEST_ASSERT(m);
    m->setProp("_Name", "core+sidechain");

    MatchVectType mv;
    SubstructMatch(*m, *core, mv);

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    // This is a rigid core: if we provide the matching substructure,
    // and do alignment, the conformations will still match.
    TEST_ASSERT(!compareConfs(m, core, mv, false));

    {
      CoordGen::CoordGenParams params;
      params.templateMol = core;
      CoordGen::addCoords(*m, &params);
      TEST_ASSERT(m->getNumConformers() == 1);
      m->setProp("_Name", "templated");
      mb = MolToMolBlock(*m);
      std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv, true, -1, -1, 0.3));
    }
    delete m;
    delete core;
  }

  {
    ROMol* core = SmilesToMol("CC(N)CC");
    TEST_ASSERT(core);
    core->setProp("_Name", "core");

    CoordGen::addCoords(*core);
    TEST_ASSERT(core->getNumConformers() == 1);
    auto mb = MolToMolBlock(*core);
    std::cerr << mb << std::endl;

    ROMol* m = SmilesToMol("CC(N)CC(O)C");
    TEST_ASSERT(m);
    m->setProp("_Name", "core+sidechain");

    MatchVectType mv;
    SubstructMatch(*m, *core, mv);

    CoordGen::addCoords(*m);
    TEST_ASSERT(m->getNumConformers() == 1);
    mb = MolToMolBlock(*m);
    std::cerr << mb << std::endl;
    // This mol is just slightly bigger than the core, providing
    // the matching substructure and doing alignment will cause
    // the conformations to match.
    TEST_ASSERT(!compareConfs(m, core, mv, false));

    {
      CoordGen::CoordGenParams params;
      params.templateMol = core;
      CoordGen::addCoords(*m, &params);
      TEST_ASSERT(m->getNumConformers() == 1);
      m->setProp("_Name", "templated");
      mb = MolToMolBlock(*m);
      std::cerr << mb << std::endl;
      TEST_ASSERT(compareConfs(m, core, mv, true, -1, -1, 0.05));
    }
    delete m;
    delete core;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub1929() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github1929: make sure coordgen works with bogus file names"
      << std::endl;
  {
    ROMol* m = SmilesToMol("c1cc(CC)cnc1CC(=O)O");
    TEST_ASSERT(m);
    m->setProp("_Name", "test1");
    CoordGen::CoordGenParams params;
    params.templateFileDir = "I_do_not_exist";

    TEST_ASSERT(CoordGen::addCoords(*m, &params) == 0);
    TEST_ASSERT(m->getNumConformers() == 1);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testGithub3131() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "testing github3131: results from coordgen are sometimes not centered"
      << std::endl;
  {
    auto m1 =
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC="
        "C5"_smiles;
    TEST_ASSERT(m1);
    TEST_ASSERT(CoordGen::addCoords(*m1) == 0);
    TEST_ASSERT(m1->getNumConformers() == 1);
    auto center = MolTransforms::computeCentroid(m1->getConformer());
    TEST_ASSERT(feq(center.x, 0.0));
    TEST_ASSERT(feq(center.y, 0.0));
  }

  {
    auto m1 =
        "CCC1=C2N=C(C=C(N2N=C1)NCC3=C[N+](=CC=C3)[O-])N4CCCC[C@H]4CCO"_smiles;
    TEST_ASSERT(m1);
    TEST_ASSERT(CoordGen::addCoords(*m1) == 0);
    TEST_ASSERT(m1->getNumConformers() == 1);
    auto center = MolTransforms::computeCentroid(m1->getConformer());
    TEST_ASSERT(feq(center.x, 0.0));
    TEST_ASSERT(feq(center.y, 0.0));
  }

  {
    // make sure that it's not recentered if we provide a coordmap:
    auto m1 =
        "CCC1=C2N=C(C=C(N2N=C1)NCC3=C[N+](=CC=C3)[O-])N4CCCC[C@H]4CCO"_smiles;
    TEST_ASSERT(m1);
    CoordGen::CoordGenParams params;
    params.coordMap[0] = {10.0, 10.0};
    params.coordMap[1] = {11.0, 10.0};
    TEST_ASSERT(CoordGen::addCoords(*m1, &params) == 0);
    TEST_ASSERT(m1->getNumConformers() == 1);
    auto center = MolTransforms::computeCentroid(m1->getConformer());
    TEST_ASSERT(!feq(center.x, 0.0));
    TEST_ASSERT(!feq(center.y, 0.0));
  }

  {
    // make sure that it's not recentered if we provide a template:
    auto templateMol =
        "C1=C2N=C(C=C(N2N=C1)NCC3=C[N+](=CC=C3))N4CCCC[C@H]4"_smiles;
    TEST_ASSERT(templateMol);
    TEST_ASSERT(CoordGen::addCoords(*templateMol) == 0);
    TEST_ASSERT(templateMol->getNumConformers() == 1);

    auto center = MolTransforms::computeCentroid(templateMol->getConformer());
    TEST_ASSERT(feq(center.x, 0.0));
    TEST_ASSERT(feq(center.y, 0.0));

    auto m1 =
        "CCC1=C2N=C(C=C(N2N=C1)NCC3=C[N+](=CC=C3)[O-])N4CCCC[C@H]4CCO"_smiles;
    TEST_ASSERT(m1);
    CoordGen::CoordGenParams params;
    params.templateMol = templateMol.get();
    TEST_ASSERT(CoordGen::addCoords(*m1, &params) == 0);
    TEST_ASSERT(m1->getNumConformers() == 1);
    center = MolTransforms::computeCentroid(m1->getConformer());
    TEST_ASSERT(!feq(center.x, 0.0));
    TEST_ASSERT(!feq(center.y, 0.0));
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void testCoordgenMinimize() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testing coordgen minimize" << std::endl;
  {
    auto m1 =
        R"CTAB(
  Mrv2014 08052005142D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.1121 -0.3399 0 0
M  V30 2 C -3.3708 0.5474 0 0
M  V30 3 C -1.9835 0.9311 0 0
M  V30 4 C -0.7248 0.0437 0 0
M  V30 5 C 0.7926 0.3064 0 0
M  V30 6 O 1.3239 1.7518 0 0
M  V30 7 O 0.612 -1.2514 0 0
M  V30 8 C 1.3429 -0.2989 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 5 7
M  V30 7 1 7 8
M  V30 8 1 4 1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    TEST_ASSERT(m1);
    TEST_ASSERT(m1->getNumConformers() == 1);
    auto ref = R"CTAB(
     RDKit          2D

  8  8  0  0  0  0  0  0  0  0999 V2000
   -1.5738   -1.1409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7946   -0.3071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9122    0.8790    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6914    0.0452    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7517    0.2885    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2665    1.6721    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.6941   -0.8483    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1500   -0.6002    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  2  0
  5  7  1  0
  7  8  1  0
  4  1  1  0
M  END
)CTAB"_ctab;
    TEST_ASSERT(ref);
    TEST_ASSERT(ref->getNumConformers() == 1);

    CoordGen::CoordGenParams ps;
    ps.minimizeOnly = true;
    CoordGen::addCoords(*m1, &ps);
    ROMol m2(*m1);
    double rmsd = MolAlign::alignMol(m2, *ref);
    TEST_ASSERT(rmsd < 0.1);
  }
  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
#if 1
  test2();
  test1();
  testGithub1929();
  testGithub3131();
#endif
  testCoordgenMinimize();
}
