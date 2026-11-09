#include <catch2/catch_all.hpp>

#include <vector>

#include "bench_common.hpp"

#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

using namespace RDKit;

namespace {

// One embedded conformer per sample molecule (heavy atoms only)
std::vector<RWMol> embeddedSamples() {
  auto params = DGeomHelpers::ETKDGv3;
  params.randomSeed = 42;
  params.numThreads = 1;
  params.maxIterations = 10;
  std::vector<RWMol> result;
  for (auto &mol : bench_common::load_samples()) {
    RWMol withHs(mol);
    MolOps::addHs(withHs);
    if (DGeomHelpers::EmbedMolecule(withHs, params) < 0) {
      continue;
    }
    MolOps::removeHs(withHs);
    result.push_back(std::move(withHs));
  }
  return result;
}

// Aligns each molecule onto the next one in the set
void benchmarkAlignment(Catch::Benchmark::Chronometer &meter,
                        const std::vector<RWMol> &mols,
                        const GaussianShape::ShapeOverlayOptions &opts) {
  std::vector<std::vector<RWMol>> work(meter.runs(), mols);
  meter.measure([&](int run) {
    double total = 0.0;
    auto &fits = work[run];
    for (size_t i = 0; i + 1 < mols.size(); ++i) {
      total += GaussianShape::AlignMolecule(
          mols[i + 1], fits[i], GaussianShape::ShapeInputOptions(),
          GaussianShape::ShapeInputOptions(), nullptr, opts)[0];
    }
    return total;
  });
}

}  // namespace

TEST_CASE("Gaussian shape alignment", "[shape]") {
  const auto mols = embeddedSamples();

  BENCHMARK_ADVANCED("AlignMolecule shape only")
  (Catch::Benchmark::Chronometer meter) {
    GaussianShape::ShapeOverlayOptions opts;
    opts.optimMode = GaussianShape::OptimMode::SHAPE_ONLY;
    benchmarkAlignment(meter, mols, opts);
  };

  BENCHMARK_ADVANCED("AlignMolecule shape plus color")
  (Catch::Benchmark::Chronometer meter) {
    GaussianShape::ShapeOverlayOptions opts;
    opts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR;
    benchmarkAlignment(meter, mols, opts);
  };
}
