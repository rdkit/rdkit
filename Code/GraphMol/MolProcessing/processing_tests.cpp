//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors

//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>
#include "RDGeneral/test.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolProcessing/MolProcessing.h>
#include <RDGeneral/RDLog.h>

using namespace RDKit;

TEST_CASE("getFingerprintsForMolsInFile") {
  std::string dirName = getenv("RDBASE");
  dirName += "/Data/NCI/";
  SECTION("SDF") {
    std::string fileName = dirName + "first_200.props.sdf";
    auto res = MolProcessing::getFingerprintsForMolsInFile<>(fileName);

    CHECK(res.size() == 200);
    for (auto i = 0u; i < res.size(); ++i) {
      INFO(i);
      CHECK(res[i]);
    }
  }
  SECTION("SMILES") {
    std::string fileName = dirName + "first_5K.smi";
    GeneralMolSupplier::SupplierOptions options;
    options.titleLine = false;
    boost::logging::disable_logs("rdApp.*");
    {
      auto res =
          MolProcessing::getFingerprintsForMolsInFile<>(fileName, options);
      CHECK(res.size() == 4999);
    }
    {
      RDKit::GeneralMolSupplier::SupplierOptions options;
      options.numWriterThreads = 1;
      options.titleLine = false;
      auto res = RDKit::MolProcessing::getFingerprintsForMolsInFile<>(fileName,
                                                                      options);
      CHECK(res.size() == 4999);
    }
  }
}