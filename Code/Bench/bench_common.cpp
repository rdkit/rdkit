#include <catch2/catch_all.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <unordered_map>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "bench_common.hpp"

namespace bench_common {

const char *dataset_name(Dataset dataset) {
  switch (dataset) {
    case Dataset::Canonical:
      return "canonical";
    case Dataset::Rings_2:
      return "rings_2";
    case Dataset::Rings_3:
      return "rings_3";
    case Dataset::Rings_4:
      return "rings_4";
    case Dataset::Rings_5:
      return "rings_5";
    case Dataset::Rings_6:
      return "rings_6";
  }
  return "unknown";
}

namespace {

std::vector<std::string> load_dataset_file(Dataset dataset) {
  const auto *rdbase = std::getenv("RDBASE");
  if (!rdbase) {
    throw std::runtime_error(
        "RDBASE is required to locate Code/Bench/data/*.smi");
  }
  auto path = std::filesystem::path(rdbase) / "Code" / "Bench" / "data" /
              (std::string(dataset_name(dataset)) + ".smi");
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("could not open dataset: " + path.string());
  }

  std::vector<std::string> result;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    const auto tab = line.find('\t');
    const auto space = line.find(' ');
    const auto end = std::min(tab, space);
    result.push_back(end == std::string::npos ? line : line.substr(0, end));
  }
  return result;
}

const std::vector<std::string> &canonical_smiles() {
  static const std::vector<std::string> result(std::begin(SAMPLES),
                                               std::end(SAMPLES));
  return result;
}

}  // namespace

const std::vector<std::string> &dataset_smiles(Dataset dataset) {
  if (dataset == Dataset::Canonical) {
    return canonical_smiles();
  }
  static std::mutex mutex;
  static std::unordered_map<int, std::vector<std::string>> cache;
  std::scoped_lock lock(mutex);
  const auto key = static_cast<int>(dataset);
  auto found = cache.find(key);
  if (found != cache.end()) {
    return found->second;
  }
  return cache.emplace(key, load_dataset_file(dataset)).first->second;
}

std::vector<RDKit::ROMol> load_samples() {
  return load_samples(Dataset::Canonical);
}

std::vector<RDKit::ROMol> load_samples(Dataset dataset, bool sanitize) {
  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.sanitize = sanitize;
  params.removeHs = sanitize;

  std::vector<RDKit::ROMol> ret;
  const auto &smiles = dataset_smiles(dataset);
  ret.reserve(smiles.size());
  for (const auto &text : smiles) {
    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(text, params);
    REQUIRE(mol);
    ret.push_back(std::move(*mol));
  }
  return ret;
}

}  // namespace bench_common
