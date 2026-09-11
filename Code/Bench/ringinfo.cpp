#include <catch2/catch_all.hpp>

#include <cstdint>
#include <utility>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/ROMol.h>

#include "bench_common.hpp"

using namespace RDKit;
using bench_common::Dataset;

namespace {

struct RingSet {
  RingInfo::VECT_INT_VECT atomRings;
  RingInfo::VECT_INT_VECT bondRings;
};

std::vector<RingSet> extract_ring_sets(const std::vector<ROMol> &mols) {
  std::vector<RingSet> result;
  result.reserve(mols.size());
  for (const auto &mol : mols) {
    const auto *info = mol.getRingInfo();
    REQUIRE(info);
    REQUIRE(info->isInitialized());
    auto &entry = result.emplace_back();
    for (const auto &ring : info->atomRings()) {
      entry.atomRings.emplace_back(ring.begin(), ring.end());
    }
    for (const auto &ring : info->bondRings()) {
      entry.bondRings.emplace_back(ring.begin(), ring.end());
    }
  }
  return result;
}

template <typename Info>
void populate_ring_info(Info &info, const RingSet &rings) {
  info.initialize(FIND_RING_TYPE_SYMM_SSSR);
  if constexpr (requires { info.addRings(rings.atomRings, rings.bondRings); }) {
    info.addRings(rings.atomRings, rings.bondRings);
  } else {
    for (size_t i = 0; i < rings.atomRings.size(); ++i) {
      info.addRing(rings.atomRings[i], rings.bondRings[i]);
    }
  }
}

void populate_incrementally(RingInfo &info, const RingSet &rings) {
  info.initialize(FIND_RING_TYPE_SYMM_SSSR);
  for (size_t i = 0; i < rings.atomRings.size(); ++i) {
    info.addRing(rings.atomRings[i], rings.bondRings[i]);
  }
}

template <typename Op>
void measure_readonly(const std::vector<ROMol> &mols,
                      Catch::Benchmark::Chronometer meter, Op op) {
  meter.measure([&](int) {
    uint64_t result = 0;
    for (const auto &mol : mols) {
      result += op(mol, *mol.getRingInfo());
    }
    return result;
  });
}

template <typename Op>
void measure_ring_perception(const std::vector<ROMol> &samples,
                             Catch::Benchmark::Chronometer meter, Op op) {
  std::vector<ROMol> work;
  work.reserve(meter.runs() * samples.size());
  for (int run = 0; run < meter.runs(); ++run) {
    for (const auto &sample : samples) {
      auto &mol = work.emplace_back(sample);
      mol.getRingInfo()->reset();
    }
  }
  meter.measure([&](int run) {
    uint64_t result = 0;
    for (size_t sample = 0; sample < samples.size(); ++sample) {
      auto &mol = work[run * samples.size() + sample];
      result += op(mol);
    }
    return result;
  });
}

#define RINGINFO_BENCHMARKS(DATASET, LABEL, TAG)                               \
  TEST_CASE("RingInfo raw operations " LABEL, "[ringinfo][raw]" TAG) {         \
    auto mols = bench_common::load_samples(DATASET);                           \
    const auto ringSets = extract_ring_sets(mols);                             \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo construct " LABEL)                            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RingInfo> storage(meter.runs() * ringSets.size());           \
      meter.measure([&](int run) {                                             \
        uint64_t result = 0;                                                   \
        for (size_t i = 0; i < ringSets.size(); ++i) {                         \
          auto &info = storage[run * ringSets.size() + i];                     \
          populate_ring_info(info, ringSets[i]);                               \
          result += info.numRings();                                           \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo incremental addRing construction " LABEL)     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RingInfo> storage(meter.runs() * ringSets.size());           \
      meter.measure([&](int run) {                                             \
        uint64_t result = 0;                                                   \
        for (size_t i = 0; i < ringSets.size(); ++i) {                         \
          auto &info = storage[run * ringSets.size() + i];                     \
          populate_incrementally(info, ringSets[i]);                           \
          result += info.numRings();                                           \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo numAtomRings " LABEL)                         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &mol, const RingInfo &ri) { \
        uint64_t result = 0;                                                   \
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {                 \
          result += ri.numAtomRings(i);                                        \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo numBondRings " LABEL)                         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &mol, const RingInfo &ri) { \
        uint64_t result = 0;                                                   \
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {                 \
          result += ri.numBondRings(i);                                        \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo atomMembers traversal " LABEL)                \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &mol, const RingInfo &ri) { \
        uint64_t result = 0;                                                   \
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {                 \
          for (const auto ring : ri.atomMembers(i)) {                          \
            result += ring + 1;                                                \
          }                                                                    \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo bondMembers traversal " LABEL)                \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &mol, const RingInfo &ri) { \
        uint64_t result = 0;                                                   \
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {                 \
          for (const auto ring : ri.bondMembers(i)) {                          \
            result += ring + 1;                                                \
          }                                                                    \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo atomRings traversal " LABEL)                  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &, const RingInfo &ri) {    \
        uint64_t result = 0;                                                   \
        for (const auto &ring : ri.atomRings()) {                              \
          for (const auto atom : ring) {                                       \
            result += atom + 1;                                                \
          }                                                                    \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo bondRings traversal " LABEL)                  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &, const RingInfo &ri) {    \
        uint64_t result = 0;                                                   \
        for (const auto &ring : ri.bondRings()) {                              \
          for (const auto bond : ring) {                                       \
            result += bond + 1;                                                \
          }                                                                    \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo ring-size queries " LABEL)                    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &mol, const RingInfo &ri) { \
        uint64_t result = 0;                                                   \
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {                 \
          const auto minimum = ri.minAtomRingSize(i);                          \
          result += minimum;                                                   \
          result += minimum && ri.isAtomInRingOfSize(i, minimum);              \
        }                                                                      \
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {                 \
          const auto minimum = ri.minBondRingSize(i);                          \
          result += minimum;                                                   \
          result += minimum && ri.isBondInRingOfSize(i, minimum);              \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo ring-size vector materialization " LABEL)     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &mol, const RingInfo &ri) { \
        uint64_t result = 0;                                                   \
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {                 \
          const auto sizes = ri.atomRingSizes(i);                              \
          result += sizes.size();                                              \
        }                                                                      \
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {                 \
          const auto sizes = ri.bondRingSizes(i);                              \
          result += sizes.size();                                              \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo same-ring queries " LABEL)                    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_readonly(mols, meter, [](const ROMol &, const RingInfo &ri) {    \
        uint64_t result = 0;                                                   \
        for (const auto &ring : ri.atomRings()) {                              \
          if (ring.size() > 1) {                                               \
            result += ri.areAtomsInSameRing(ring.front(), ring.back());        \
            result += ri.areAtomsInSameRingOfSize(ring.front(), ring.back(),   \
                                                  ring.size());                \
          }                                                                    \
        }                                                                      \
        for (const auto &ring : ri.bondRings()) {                              \
          if (ring.size() > 1) {                                               \
            result += ri.areBondsInSameRing(ring.front(), ring.back());        \
            result += ri.areBondsInSameRingOfSize(ring.front(), ring.back(),   \
                                                  ring.size());                \
          }                                                                    \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
                                                                               \
    BENCHMARK_ADVANCED("RingInfo fused initialization and traversal " LABEL)   \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RingInfo> infos;                                             \
      infos.reserve(meter.runs() * ringSets.size());                           \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &rings : ringSets) {                                   \
          populate_ring_info(infos.emplace_back(), rings);                     \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        uint64_t result = 0;                                                   \
        for (size_t i = 0; i < ringSets.size(); ++i) {                         \
          auto &info = infos[run * ringSets.size() + i];                       \
          for (unsigned int ring = 0; ring < info.numRings(); ++ring) {        \
            result += info.isRingFused(ring);                                  \
            result += info.numFusedBonds(ring);                                \
            result += info.numFusedRingNeighbors(ring);                        \
            if (ring + 1 < info.numRings()) {                                  \
              result += info.areRingsFused(ring, ring + 1);                    \
            }                                                                  \
          }                                                                    \
        }                                                                      \
        return result;                                                         \
      });                                                                      \
    };                                                                         \
  }                                                                            \
                                                                               \
  TEST_CASE("Ring perception " LABEL, "[ringinfo][perception]" TAG) {          \
    auto samples = bench_common::load_samples(DATASET, false);                 \
    BENCHMARK_ADVANCED("MolOps::findSSSR " LABEL)                              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_ring_perception(samples, meter, [](ROMol &mol) {                 \
        return static_cast<uint64_t>(MolOps::findSSSR(mol));                   \
      });                                                                      \
    };                                                                         \
    BENCHMARK_ADVANCED("MolOps::symmetrizeSSSR " LABEL)                        \
    (Catch::Benchmark::Chronometer meter) {                                    \
      measure_ring_perception(samples, meter, [](ROMol &mol) {                 \
        return static_cast<uint64_t>(MolOps::symmetrizeSSSR(mol));             \
      });                                                                      \
    };                                                                         \
  }

RINGINFO_BENCHMARKS(Dataset::Rings_2, "rings 2", "[rings_2]")
RINGINFO_BENCHMARKS(Dataset::Rings_3, "rings 3", "[rings_3]")
RINGINFO_BENCHMARKS(Dataset::Rings_4, "rings 4", "[rings_4]")
RINGINFO_BENCHMARKS(Dataset::Rings_5, "rings 5", "[rings_5]")
RINGINFO_BENCHMARKS(Dataset::Rings_6, "rings 6", "[rings_6]")

}  // namespace
