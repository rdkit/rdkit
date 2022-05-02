//
//  Copyright (C) 2004-2018 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/ForceFieldHelpers/FFConvenience.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Params.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/math/special_functions/round.hpp>

using namespace RDKit;
#ifdef RDK_TEST_MULTITHREADED
namespace {
void runblock_mmff(const std::vector<ROMol *> &mols) {
  for (auto mol : mols) {
    ForceFields::ForceField *field = MMFF::constructForceField(*mol);
    TEST_ASSERT(field);
    field->initialize();
    field->minimize(1);
    delete field;
  }
}
}  // namespace
#include <thread>
#include <future>
void testMMFFMultiThread() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test MMFF multithreading" << std::endl;

  std::string pathName = getenv("RDBASE");
  pathName += "/Code/GraphMol/ForceFieldHelpers/MMFF/test_data";
  SDMolSupplier suppl(pathName + "/bulk.sdf");
  unsigned int count = 24;
  std::vector<std::vector<ROMol *>> mols;
  for (unsigned int i = 0; i < count; ++i) {
    mols.emplace_back();
  }

  while (!suppl.atEnd() && mols.size() < 100) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
      for (unsigned int i = 0; i < count; ++i) {
        if (i == 0) {
          mols[i].push_back(mol);
        } else {
          mols[i].push_back(new ROMol(*mol));
        }
      }
    } catch (...) {
      continue;
    }
  }

  std::vector<std::future<void>> tg;

  std::cerr << "processing" << std::endl;
  for (unsigned int i = 0; i < count; ++i) {
    std::cerr << " launch :" << i << std::endl;
    std::cerr.flush();
    tg.emplace_back(std::async(std::launch::async, runblock_mmff, mols[i]));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  std::cerr << "done" << std::endl;
  for (unsigned int i = 0; i < count; ++i) {
    for (auto *mol : mols[i]) {
      delete mol;
    }
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

#endif
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
  // we get a ton of warnings here about missing Hs... disable them
  boost::logging::disable_logs("rdApp.warning");

#ifdef RDK_TEST_MULTITHREADED
  testMMFFMultiThread();
#endif
}
