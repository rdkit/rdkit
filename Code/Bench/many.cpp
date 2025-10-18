
#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

using namespace RDKit;

TEST_CASE("MolPickler::pickleMol", "[pickle][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    std::string buf;
    BENCHMARK("MolPickler::pickleMol " + std::string(smiles)) {
      buf.clear();
      MolPickler::pickleMol(*mol, buf);
      return buf.size();
    };
  }
}

TEST_CASE("MolPickler::molFromPickle", "[pickle][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    std::string pickled;
    MolPickler::pickleMol(*mol, pickled);
    BENCHMARK("MolPickler::molFromPickle " + std::string(smiles)) {
      ROMol res;
      MolPickler::molFromPickle(pickled, res);
      return res;
    };
  }
}

TEST_CASE("ROMol copy constructor", "[mol][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK_ADVANCED("ROMol copy constructor")(
        Catch::Benchmark::Chronometer meter) {
      std::vector<Catch::Benchmark::storage_for<ROMol>> storage(meter.runs());
      meter.measure([&](int i) { storage[i].construct(*mol); });
    };
  }
}

TEST_CASE("ROMol destructor", "[mol][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK_ADVANCED("ROMol destructor")(
        Catch::Benchmark::Chronometer meter) {
      std::vector<Catch::Benchmark::destructable_object<ROMol>> storage(
          meter.runs());
      for (auto &o : storage) {
        o.construct(*mol);
      }
      meter.measure([&](int i) { storage[i].destruct(); });
    };
  }
}

TEST_CASE("MolOps::addHs", "[molops][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("MolOps::addHs " + std::string(smiles)) {
      auto mol_copy = RWMol(*mol);
      MolOps::addHs(mol_copy);
      return mol_copy;
    };
  }
}

TEST_CASE("MolOps::calcNumBridgeheadAtoms", "[molops][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("MolOps::calcNumBridgeheadAtoms " + std::string(smiles)) {
      return Descriptors::calcNumBridgeheadAtoms(*mol);
    };
  }
}

TEST_CASE("MolOps::FindSSR", "[molops][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("MolOps::FindSSR " + std::string(smiles)) {
      return MolOps::findSSSR(*mol);
    };
  }
}

TEST_CASE("meta bench_common::nth_random", "![meta][many]") {
  BENCHMARK_ADVANCED("meta bench_common::nth_random ")(auto meter) {
    meter.measure([&](int i) { return bench_common::nth_random(i); });
  };
}

TEST_CASE("memory pressure test", "[size][many]") {
  std::vector<ROMol> cases;
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    cases.push_back(*mol);
  }
  REQUIRE(!cases.empty());

  const size_t N = 100000;
  std::vector<ROMol> mols;
  mols.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    mols.emplace_back(cases[i % cases.size()]);
  }
  REQUIRE(mols.size() == N);

  BENCHMARK("memory pressure test", i) {
    auto a = bench_common::nth_random(i);
    auto src_idx = a % mols.size();
    auto dst_idx = (a / mols.size()) % mols.size();
    ROMol temp(mols[src_idx]);
    mols[dst_idx] = std::move(temp);
  };
}

TEST_CASE("MolOps::calcNumSpiroAtoms", "[molops][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("MolOps::calcNumSpiroAtoms " + std::string(smiles)) {
      return Descriptors::calcNumSpiroAtoms(*mol);
    };
  }
}

TEST_CASE("ROMol::getNumHeavyAtoms", "[mol][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("ROMol::getNumHeavyAtoms " + std::string(smiles)) {
      return mol->getNumHeavyAtoms();
    };
  }
}

TEST_CASE("MolOps::getMolFrags", "[molops][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("MolOps::getMolFrags " + std::string(smiles)) {
      std::vector<std::unique_ptr<ROMol>> frags;
      MolOps::getMolFrags(*mol, frags);
      return frags.size();
    };
  }
}

TEST_CASE("MorganFingerprints::getFingerprint", "[fingerprint][many]") {
  for (auto smiles : bench_common::CASES) {
    auto mol = v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    BENCHMARK("MorganFingerprints::getFingerprint " + std::string(smiles)) {
      std::unique_ptr<SparseIntVect<uint32_t>> fp(
          MorganFingerprints::getFingerprint(*mol, 20, nullptr, nullptr, false,
                                             false, false, false));
      return fp->getTotalVal();
    };
  }
}

// TODO:
// - fill in more benches
// - Substructure matching
