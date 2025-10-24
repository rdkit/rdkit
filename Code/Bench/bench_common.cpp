#include <catch2/catch_all.hpp>
#include <boost/random/splitmix64.hpp>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "bench_common.hpp"

namespace bench_common {

std::vector<RDKit::ROMol> load_samples() {
  std::vector<RDKit::ROMol> ret;
  for (auto smiles : SAMPLES) {
    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    ret.emplace_back(std::move(*mol));
  }
  return ret;
}

uint64_t nth_random(uint64_t n) noexcept {
  boost::random::splitmix64 gen;
  gen.seed(n);
  return gen.next();
}

}  // namespace bench_common
