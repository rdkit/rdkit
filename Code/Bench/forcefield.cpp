#include <catch2/catch_all.hpp>

#include <vector>

#include "bench_common.hpp"

#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>

using namespace RDKit;

namespace {

constexpr unsigned int numConfs = 5;

DGeomHelpers::EmbedParameters makeEmbedParams() {
  auto params = DGeomHelpers::ETKDGv3;
  params.randomSeed = 42;
  params.numThreads = 1;
  // keep pathological molecules from dominating the timings
  params.maxIterations = 10;
  return params;
}

// EmbedMultipleConfs takes the parameters by non-const reference
DGeomHelpers::EmbedParameters &embedParams() {
  static auto params = makeEmbedParams();
  return params;
}

std::vector<RWMol> samplesWithHs() {
  std::vector<RWMol> mols;
  for (auto &mol : bench_common::load_samples()) {
    RWMol withHs(mol);
    MolOps::addHs(withHs);
    mols.push_back(std::move(withHs));
  }
  return mols;
}

// One single-conformer molecule per conformer ETKDG produced, so that each
// minimization in the benchmarks below starts from the same geometry.
std::vector<RWMol> startingConformers() {
  std::vector<RWMol> result;
  for (auto &mol : samplesWithHs()) {
    const auto cids =
        DGeomHelpers::EmbedMultipleConfs(mol, numConfs, embedParams());
    for (const auto cid : cids) {
      RWMol single(mol, /*quickCopy=*/false, cid);
      result.push_back(std::move(single));
    }
  }
  return result;
}

template <typename Optimize>
void benchmarkMinimization(Catch::Benchmark::Chronometer &meter,
                           const std::vector<RWMol> &starts,
                           Optimize optimize) {
  std::vector<std::vector<RWMol>> work(meter.runs(), starts);
  meter.measure([&](int run) {
    double total = 0.0;
    for (auto &mol : work[run]) {
      total += optimize(mol).second;
    }
    return total;
  });
}

}  // namespace

TEST_CASE("ETKDGv3 embedding", "[forcefield][embed]") {
  const auto mols = samplesWithHs();
  BENCHMARK_ADVANCED("ETKDGv3 EmbedMultipleConfs")
  (Catch::Benchmark::Chronometer meter) {
    std::vector<std::vector<RWMol>> work(meter.runs(), mols);
    meter.measure([&](int run) {
      size_t total = 0;
      for (auto &mol : work[run]) {
        total +=
            DGeomHelpers::EmbedMultipleConfs(mol, numConfs, embedParams())
                .size();
      }
      return total;
    });
  };
}

TEST_CASE("Force field minimization", "[forcefield][minimize]") {
  const auto starts = startingConformers();

  BENCHMARK_ADVANCED("MMFFOptimizeMolecule")
  (Catch::Benchmark::Chronometer meter) {
    benchmarkMinimization(meter, starts, [](RWMol &mol) {
      return MMFF::MMFFOptimizeMolecule(mol);
    });
  };

  BENCHMARK_ADVANCED("UFFOptimizeMolecule")
  (Catch::Benchmark::Chronometer meter) {
    benchmarkMinimization(meter, starts, [](RWMol &mol) {
      return UFF::UFFOptimizeMolecule(mol);
    });
  };
}
