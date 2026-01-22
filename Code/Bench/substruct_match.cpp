#include <catch2/catch_all.hpp>
#include <string>
#include <vector>
#include <istream>
#include <filesystem>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

namespace bench_substruct_match {

std::filesystem::path relative_to_rdbase(
    const std::filesystem::path &relative) {
  char *rdbase = std::getenv("RDBASE");
  if (!rdbase) {
    throw std::runtime_error("RDBASE environment variable not set");
  }
  std::filesystem::path path(rdbase);
  path /= relative;
  return path;
}

std::vector<std::string> parse_smarts_examples(std::istream &in) {
  std::vector<std::string> examples;
  std::string line;
  while (std::getline(in, line)) {
    if (line.starts_with("#") || line.starts_with(" ") || line.empty()) {
      continue;
    }
    std::istringstream ss(line);
    std::string smarts;
    std::getline(ss, smarts, ' ');
    examples.push_back(line);
  }
  return examples;
}

std::vector<std::string> load_smarts_examples(std::filesystem::path &path) {
  std::ifstream file(path);
  return parse_smarts_examples(file);
}

std::vector<ROMol> load_smarts_queries(const std::filesystem::path &relative) {
  auto path = relative_to_rdbase(relative);
  auto examples = load_smarts_examples(path);
  std::vector<ROMol> queries;
  for (const auto &example : examples) {
    auto query = v2::SmilesParse::MolFromSmarts(example);
    REQUIRE(query);
    queries.push_back(std::move(*query));
  }
  return queries;
}

}  // namespace bench_substruct_match

TEST_CASE("ROMol::GetSubstructMatch", "[substruct_match]") {
  auto samples = bench_common::load_samples();

  auto query = v2::SmilesParse::MolFromSmarts(
      "[#6,#7]1:[#6,#7]:[#6](:[#6,#7]:[#6,#7]:[#6]:1-[#17,#35,#53,#9])-[#5](-[#8])-[#8]");
  BENCHMARK("ROMol::GetSubstructMatch") {
    auto total = 0;
    for (auto &mol : samples) {
      auto matches = SubstructMatch(mol, *query);
      total += matches.size();
    }
    return total;
  };
}

TEST_CASE("ROMol::GetSubstructMatch RLewis", "[substruct_match]") {
  auto mols = bench_common::load_samples();
  auto queries = bench_substruct_match::load_smarts_queries(
      "Data/SmartsLib/RLewis_smarts.txt");

  BENCHMARK("ROMol::GetSubstructMatch RLewis") {
    auto total = 0;
    for (const auto &mol : mols) {
      for (const auto &query : queries) {
        auto matches = SubstructMatch(mol, query);
        total += matches.size();
      }
    }
    return total;
  };
}

TEST_CASE("ROMol::GetSubstructMatch patty", "[substruct_match]") {
  auto mols = bench_common::load_samples();
  auto queries = bench_substruct_match::load_smarts_queries(
      "Data/SmartsLib/patty_rules.txt");

  BENCHMARK("ROMol::GetSubstructMatch patty") {
    auto total = 0;
    for (const auto &mol : mols) {
      for (const auto &query : queries) {
        auto matches = SubstructMatch(mol, query);
        total += matches.size();
      }
    }
    return total;
  };
}
