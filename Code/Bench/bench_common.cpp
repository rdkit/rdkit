#include <catch2/catch_all.hpp>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "bench_common.hpp"

namespace bench_common {

std::vector<RDKit::ROMol> load_samples() {
  std::vector<RDKit::ROMol> ret;
  for (auto smiles : SAMPLES) {
    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    ret.push_back(std::move(*mol));
  }
  return ret;
}

}  // namespace bench_common
