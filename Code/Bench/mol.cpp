#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("ROMol copy constructor", "[mol]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("ROMol copy constructor")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<Catch::Benchmark::storage_for<ROMol>> storage(meter.runs() *
                                                              samples.size());
    meter.measure([&](int i) {
      for (size_t sample = 0; sample < samples.size(); ++sample) {
        storage[i * samples.size() + sample].construct(samples[sample]);
      }
    });
  };
}

TEST_CASE("ROMol destructor", "[mol]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("ROMol destructor")(Catch::Benchmark::Chronometer meter) {
    std::vector<Catch::Benchmark::destructable_object<ROMol>> storage(
        meter.runs() * samples.size());
    for (size_t i = 0; i < storage.size(); ++i) {
      storage[i].construct(samples[i % samples.size()]);
    }
    meter.measure([&](int i) {
      for (size_t sample = 0; sample < samples.size(); ++sample) {
        storage[i * samples.size() + sample].destruct();
      }
    });
  };
}

TEST_CASE("memory pressure test", "[mol][size]") {
  auto cases = bench_common::load_samples();
  REQUIRE(!cases.empty());

  const size_t N = 100000;
  std::vector<ROMol> mols;
  mols.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    mols.emplace_back(cases[i % cases.size()]);
  }
  REQUIRE(mols.size() == N);

  BENCHMARK("memory pressure test", i) {
    // copy from one random location to another
    auto a = bench_common::nth_random(i);
    auto src_idx = a % mols.size();
    auto dst_idx = (a / mols.size()) % mols.size();
    ROMol temp(mols[src_idx]);
    mols[dst_idx] = std::move(temp);
  };
}

TEST_CASE("ROMol::getNumHeavyAtoms", "[mol]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("ROMol::getNumHeavyAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += mol.getNumHeavyAtoms();
    }
    return sum;
  };
}
