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
#include <GraphMol/MolProcessing/Processing.h>
#include <RDGeneral/RDLog.h>

TEST_CASE("getFingerprintsForMolsInFile") {
#if 1
  std::string dirName = getenv("RDBASE");
  dirName += "/Data/NCI/";
  std::string fileName = dirName + "first_200.props.sdf";
  auto res = RDKit::MolProccesing::getFingerprintsForMolsInFile<>(fileName);

  CHECK(res.first.size() == 201);
  CHECK(res.second.size() == res.first.size());
  CHECK(res.second.count() <= res.first.size());
  for (auto i = 0u; i < res.first.size(); ++i) {
    if (res.second[i]) {
      CHECK(res.first[i]);
    } else {
      CHECK(!res.first[i]);
    }
  }
#else
  std::string fileName =
      "/home/glandrum/Datasets/COD/COD_2021aug02.organic.sdf.gz";
  RDKit::GeneralMolSupplier::SupplierOptions options;
  boost::logging::disable_logs("rdApp.*");
  {
    std::cerr << "defaults!" << std::endl;
    auto res = RDKit::MolProccesing::getFingerprintsForMolsInFile<>(fileName);
  }
  {
    std::cerr << "one thread" << std::endl;
    RDKit::GeneralMolSupplier::SupplierOptions options;
    options.numWriterThreads = 1;
    auto res =
        RDKit::MolProccesing::getFingerprintsForMolsInFile<>(fileName, options);
  }
#endif
}