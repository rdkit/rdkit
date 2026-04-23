//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/MolOps.h>
#include <ForceField/FiniteDifference.h>
#include <ForceField/ForceField.h>

#include <memory>
#include <string>

using namespace RDKit;

namespace {

constexpr double FD_TOLERANCE = 1e-4;

std::unique_ptr<ForceFields::ForceField> buildFF(const std::string &ffName,
                                                 RWMol &mol) {
  if (ffName == "UFF") {
    MolOps::sanitizeMol(mol);
    auto [params, allFound] = UFF::getAtomTypes(mol);
    if (!allFound) {
      return nullptr;
    }
    return std::unique_ptr<ForceFields::ForceField>(
        UFF::constructForceField(mol));
  }
  if (ffName == "MMFF94" || ffName == "MMFF94s") {
    MMFF::sanitizeMMFFMol(mol);
    MMFF::MMFFMolProperties props(mol, ffName);
    if (!props.isValid()) {
      return nullptr;
    }
    return std::unique_ptr<ForceFields::ForceField>(
        MMFF::constructForceField(mol, &props));
  }
  if (ffName == "ETKDG") {
    MolOps::sanitizeMol(mol);
    return DGeomHelpers::testing::getETKDGForceField(mol,
                                                     DGeomHelpers::ETKDGv3);
  }
  return nullptr;
}

}  // namespace

TEST_CASE("ForceField gradient validation: MMFF suite", "[ffvalidation]") {
  auto [ffName, sdfName] = GENERATE(
      table<std::string, std::string>({{"UFF", "MMFF94_dative.sdf"},
                                       {"UFF", "MMFF94_hypervalent.sdf"},
                                       {"MMFF94", "MMFF94_dative.sdf"},
                                       {"MMFF94", "MMFF94_hypervalent.sdf"},
                                       {"MMFF94s", "MMFF94s_dative.sdf"},
                                       {"MMFF94s", "MMFF94s_hypervalent.sdf"},
                                       {"ETKDG", "MMFF94_dative.sdf"},
                                       {"ETKDG", "MMFF94_hypervalent.sdf"}}));

  const char *rdbase = getenv("RDBASE");
  REQUIRE(rdbase);
  std::string path = std::string(rdbase) + "/Code/ForceField/MMFF/test_data/";
  SDMolSupplier suppl(path + sdfName, false, false);

  for (unsigned int i = 0; i < suppl.length(); ++i) {
    std::unique_ptr<RWMol> mol(dynamic_cast<RWMol *>(suppl[i]));
    if (!mol) {
      continue;
    }
    CAPTURE(ffName, sdfName, i);

    auto ff = buildFF(ffName, *mol);
    if (!ff) {
      continue;
    }
    ff->initialize();
    // If UFF can't handle MMFF molecules, it will return a large energy.
    if (ff->calcEnergy() > 1e6) {
      continue;
    }

    double delta = ForceFields::calcFiniteDifference(*ff);
    CHECK(delta < FD_TOLERANCE);
  }
}
