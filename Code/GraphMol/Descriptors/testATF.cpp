//  Created by Guillaume GODIN
//  Copyright (C) 2012-2020 Greg Landrum
//   @@ All Rights Reserved @@
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <chrono>  // for high_resolution_clock

#include <GraphMol/Descriptors/AtomFeat.h>

using namespace RDKit;

void test1() {
  auto m = "CO"_smiles;

  TEST_ASSERT(m);

  std::vector<double> res;

  int atomid = 0;
  RDKit::Descriptors::AtomFeatVect(*m, res, atomid);

  std::vector<double> exp{0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
                          0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 1., 0., 0.5};

  // 49 features
  TEST_ASSERT(res.size() == 49);

  for (std::size_t i = 0; i < res.size(); i++) {
    TEST_ASSERT(fabs(res[i] - exp[i]) < 0.001);
    // std::cout << res[i] << "," ;
  }
}

void test2() {
  auto m = "C1CC=C1[NH3+]"_smiles;

  TEST_ASSERT(m);

  std::vector<double> res;

  std::vector<double> molatf{
      0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0.2, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0, 0.2, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0, 0.2,
      0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0, 0,
      0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0.2, 0, 0, 1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 1, 0, 0.2};

  for (unsigned int j = 0; j < m->getNumAtoms(); j++) {
    RDKit::Descriptors::AtomFeatVect(*m, res, j);

    // 49 features
    TEST_ASSERT(res.size() == 49);

    for (std::size_t i = 0; i < res.size(); i++) {
      TEST_ASSERT(fabs(res[i] - molatf[i + j * 49]) < 0.001);
    }

    res.clear();
    res.resize(49);
  }
}

void test3() {
  auto m = "O[C@H](F)[C@H](F)O"_smiles;

  TEST_ASSERT(m);

  std::vector<double> molatf{
      0, 0, 0, 1, 0, 0, 0,        0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,        0, 0, 1, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1,        0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0.166667, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,        0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0,        0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,        0, 0, 1, 0, 0, 0, 0.166667,
      0, 0, 0, 0, 0, 1, 0,        0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,        0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1,        0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0.166667, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,        0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0,        0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,        0, 0, 1, 0, 0, 0, 0.166667,
      0, 0, 0, 0, 0, 1, 0,        0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,        0, 0, 1, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1,        0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0.166667, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0,        0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 1,        0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0,        0, 0, 1, 0, 0, 0, 0.166667};

  std::vector<double> res;

  for (unsigned int j = 0; j < m->getNumAtoms(); j++) {
    RDKit::Descriptors::AtomFeatVect(*m, res, j);

    // 49 features
    TEST_ASSERT(res.size() == 49);

    for (std::size_t i = 0; i < res.size(); i++) {
      TEST_ASSERT(fabs(res[i] - molatf[i + j * 49]) < 0.001);
    }

    res.clear();
    res.resize(49);
  }
}

void test4() {
  auto m = "O[C@H](F)[C@H](F)O"_smiles;

  TEST_ASSERT(m);

  std::vector<double> res;

  std::vector<double> molatf{0, 0, 0, 1, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 0, 1, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0.166667, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 0, 1, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0.166667, 0, 0, 0,
                             0, 0, 0, 1, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 0, 0, 0};

  for (unsigned int j = 0; j < m->getNumAtoms(); j++) {
    RDKit::Descriptors::AtomFeatVect(*m, res, j, true);

    // 52 features
    TEST_ASSERT(res.size() == 52);

    for (std::size_t i = 0; i < res.size(); i++) {
      TEST_ASSERT(fabs(res[i] - molatf[i + j * 52]) < 0.001);
    }

    res.clear();
    res.resize(52);
  }
}

void test5() {
  auto m = "O[C@@H](F)[C@@H](F)O"_smiles;

  TEST_ASSERT(m);

  std::vector<double> res;

  std::vector<double> molatf{0, 0, 0, 1, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 1, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0.166667, 0, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 1, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0, 0.166667, 0, 0, 0,
                             0, 0, 0, 1, 0, 0, 0, 0, 0, 0,        0, 0, 0,
                             0, 0, 1, 0, 0, 0, 0, 0, 0, 0,        1, 0, 0,
                             0, 1, 0, 0, 0, 0, 0, 0, 1, 0,        0, 0, 0,
                             0, 0, 0, 0, 0, 1, 0, 0, 0, 0.166667, 0, 0, 0};

  for (unsigned int j = 0; j < m->getNumAtoms(); j++) {
    RDKit::Descriptors::AtomFeatVect(*m, res, j, true);

    // 52 features
    TEST_ASSERT(res.size() == 52);

    for (std::size_t i = 0; i < res.size(); i++) {
      TEST_ASSERT(fabs(res[i] - molatf[i + j * 52]) < 0.001);
    }

    res.clear();
    res.resize(52);
  }
}

int main() {
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  test4();
  test5();
}
