//
//  Copyright (c) 2024 Greg Landrum
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDGeneral/test.h"
#include <string>
#include <catch2/catch_all.hpp>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/MultithreadedMolSupplier.h>
#include <GraphMol/FileParsers/MultithreadedSDMolSupplier.h>
#include <GraphMol/FileParsers/MultithreadedSmilesMolSupplier.h>

using namespace RDKit;

TEST_CASE("multithreaded supplier destruction without reading") {
  v2::FileParsers::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = 4;
  std::string rdbase = getenv("RDBASE");
  std::string sdpath = rdbase + "/Data/NCI/first_200.props.sdf";
  auto sdsuppl = v2::FileParsers::MultithreadedSDMolSupplier(sdpath, params);
}

TEST_CASE("callbacks SDF") {
  v2::FileParsers::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = 4;
  std::string rdbase = getenv("RDBASE");
  std::string sdpath = rdbase + "/Data/NCI/first_200.props.sdf";
  auto sdsuppl = v2::FileParsers::MultithreadedSDMolSupplier(sdpath, params);
  // std::string smiPath = rdbase + "/Data/NCI/first_5K.smi";
  // auto smisuppl =
  //     v2::FileParsers::MultithreadedSmilesMolSupplier(smiPath, params);
  // {
  //   std::string smiPath = rdbase + "/Data/NCI/first_5K.smi";
  //   auto suppl =
  //       v2::FileParsers::MultithreadedSmilesMolSupplier(smiPath, params);
  // }
  SECTION("nextCallback") {
    std::map<unsigned int, unsigned int> callbackNats;
    auto callback =
        [&](RWMol &mol,
            const v2::FileParsers::MultithreadedMolSupplier &suppl) {
          callbackNats[suppl.getLastRecordId()] = mol.getNumAtoms();
        };
    auto &suppl = sdsuppl;
    suppl.setNextCallback(callback);
    unsigned int nMols = 0;
    while (!suppl.atEnd()) {
      auto mol = suppl.next();
      if (!mol) {
        continue;
      }
      auto nats = mol->getNumAtoms();
      CHECK(callbackNats[suppl.getLastRecordId()] == nats);
      ++nMols;
    }
    CHECK(nMols == callbackNats.size());
  }
  SECTION("writeCallback") {
    auto callback = [](RWMol &mol, const std::string &, unsigned int recordId) {
      MolOps::addHs(mol);
      mol.setProp("recordId", recordId);
    };
    auto &suppl = sdsuppl;
    suppl.setWriteCallback(callback);
    while (!suppl.atEnd()) {
      auto mol = suppl.next();
      if (!mol) {
        continue;
      }

      CHECK(!MolOps::needsHs(*mol));
      CHECK(mol->hasProp("recordId"));
      CHECK(mol->getProp<unsigned int>("recordId") == suppl.getLastRecordId());
    }
  }
  SECTION("readCallback") {
    auto callback = [](const std::string &sdf, unsigned int recordId) {
      auto res = sdf;
      auto pos = sdf.find("$$$$");
      if (pos == std::string::npos) {
        return res;
      }
      res.replace(pos, 4,
                  std::string(">  <recordId>\n") + std::to_string(recordId) +
                      "\n\n$$$$");
      return res;
    };
    auto &suppl = sdsuppl;

    suppl.setReadCallback(callback);
    while (!suppl.atEnd()) {
      auto mol = suppl.next();
      if (!mol) {
        continue;
      }
      CHECK(mol->hasProp("recordId"));
      CHECK(mol->getProp<unsigned int>("recordId") == suppl.getLastRecordId());
    }
  }
}

TEST_CASE("callbacks smiles") {
  v2::FileParsers::MultithreadedMolSupplier::Parameters params;
  params.numWriterThreads = 4;
  std::string rdbase = getenv("RDBASE");
  std::string smiPath = rdbase + "/Data/NCI/first_5K.smi";
  auto suppl = v2::FileParsers::MultithreadedSmilesMolSupplier(smiPath, params);
  SECTION("nextCallback") {
    std::map<unsigned int, unsigned int> callbackNats;
    auto callback =
        [&](RWMol &mol,
            const v2::FileParsers::MultithreadedMolSupplier &suppl) {
          callbackNats[suppl.getLastRecordId()] = mol.getNumAtoms();
        };
    suppl.setNextCallback(callback);
    unsigned int nMols = 0;
    while (!suppl.atEnd()) {
      auto mol = suppl.next();
      if (!mol) {
        continue;
      }
      auto nats = mol->getNumAtoms();
      CHECK(callbackNats[suppl.getLastRecordId()] == nats);
      ++nMols;
    }
    CHECK(nMols == callbackNats.size());
  }
  SECTION("writeCallback") {
    auto callback = [](RWMol &mol, const std::string &, unsigned int recordId) {
      MolOps::addHs(mol);
      mol.setProp("recordId", recordId);
    };
    suppl.setWriteCallback(callback);
    while (!suppl.atEnd()) {
      auto mol = suppl.next();
      if (!mol) {
        continue;
      }

      CHECK(!MolOps::needsHs(*mol));
      CHECK(mol->hasProp("recordId"));
      CHECK(mol->getProp<unsigned int>("recordId") == suppl.getLastRecordId());
    }
  }
}
