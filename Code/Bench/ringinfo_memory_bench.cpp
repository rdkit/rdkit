// Copyright (c) 2026, RDKit contributors
// All rights reserved.
//
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the root of the RDKit source tree.

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef __linux__
#include <unistd.h>
#endif

#include <GraphMol/RingInfo.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace {

constexpr std::size_t kDefaultCount = 100'000;

struct RingSet {
  RDKit::RingInfo::VECT_INT_VECT atoms;
  RDKit::RingInfo::VECT_INT_VECT bonds;
};

std::filesystem::path defaultInputPath() {
  const auto *rdbase = std::getenv("RDBASE");
  if (!rdbase) {
    return "Code/Bench/data/rings_4.smi";
  }
  return std::filesystem::path(rdbase) / "Code" / "Bench" / "data" /
         "rings_4.smi";
}

std::size_t parseCount(const char *text) {
  std::size_t consumed = 0;
  const auto result = std::stoull(text, &consumed);
  if (consumed != std::string(text).size() || !result) {
    throw std::invalid_argument("count must be a positive integer");
  }
  return result;
}

std::vector<std::string> readSmiles(const std::filesystem::path &path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("could not open input: " + path.string());
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
  if (result.empty()) {
    throw std::runtime_error("input contains no SMILES: " + path.string());
  }
  return result;
}

std::vector<RingSet> loadRingSets(const std::filesystem::path &path) {
  std::vector<RingSet> result;
  for (const auto &smiles : readSmiles(path)) {
    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles);
    if (!mol) {
      throw std::runtime_error("could not parse input SMILES");
    }
    const auto *info = mol->getRingInfo();
    auto &rings = result.emplace_back();
    for (const auto &ring : info->atomRings()) {
      rings.atoms.emplace_back(ring.begin(), ring.end());
    }
    for (const auto &ring : info->bondRings()) {
      rings.bonds.emplace_back(ring.begin(), ring.end());
    }
  }
  return result;
}

template <typename Info>
void populate(Info &info, const RingSet &rings) {
  info.initialize(RDKit::FIND_RING_TYPE_SYMM_SSSR);
  if constexpr (requires { info.addRings(rings.atoms, rings.bonds); }) {
    info.addRings(rings.atoms, rings.bonds);
  } else {
    for (size_t i = 0; i < rings.atoms.size(); ++i) {
      info.addRing(rings.atoms[i], rings.bonds[i]);
    }
  }
}

std::size_t residentBytes() {
#ifdef __linux__
  std::ifstream statm("/proc/self/statm");
  std::size_t totalPages = 0, residentPages = 0;
  if (!(statm >> totalPages >> residentPages)) {
    throw std::runtime_error("could not read /proc/self/statm");
  }
  return residentPages * static_cast<std::size_t>(::sysconf(_SC_PAGESIZE));
#else
  throw std::runtime_error(
      "resident memory measurement is currently implemented for Linux only");
#endif
}

std::size_t nonnegativeDifference(std::size_t after, std::size_t before) {
  return after >= before ? after - before : 0;
}

}  // namespace

int main(int argc, char *argv[]) {
  if (argc > 3) {
    std::cerr << "Usage: " << argv[0] << " [ring-smiles-file] [count]\n";
    return 2;
  }
  try {
    const auto input =
        argc >= 2 ? std::filesystem::path(argv[1]) : defaultInputPath();
    const auto count = argc == 3 ? parseCount(argv[2]) : kDefaultCount;
    const auto patterns = loadRingSets(input);

    const auto baselineRss = residentBytes();
    std::vector<RDKit::RingInfo> infos(count);
    const auto objectRss = residentBytes();
    for (size_t i = 0; i < count; ++i) {
      populate(infos[i], patterns[i % patterns.size()]);
    }
    const auto populatedRss = residentBytes();
    const auto objectDelta = nonnegativeDifference(objectRss, baselineRss);
    const auto populationDelta = nonnegativeDifference(populatedRss, objectRss);
    const auto retainedDelta = nonnegativeDifference(populatedRss, baselineRss);
    const auto objectBytesPerRingInfo =
        static_cast<double>(objectDelta) / count;
    const auto ringBytesPerRingInfo =
        static_cast<double>(populationDelta) / count;
    const auto totalBytesPerRingInfo =
        static_cast<double>(retainedDelta) / count;

    std::size_t rings = 0, members = 0;
    for (const auto &info : infos) {
      rings += info.numRings();
      for (const auto &ring : info.atomRings()) {
        members += ring.size();
      }
    }
    std::cout << "input=" << input << " count=" << count
              << " sizeof_ringinfo=" << sizeof(RDKit::RingInfo)
              << " object_bytes=" << count * sizeof(RDKit::RingInfo)
              << " baseline_rss_bytes=" << baselineRss
              << " object_rss_bytes=" << objectRss
              << " populated_rss_bytes=" << populatedRss
              << " object_rss_delta_bytes=" << objectDelta
              << " population_rss_delta_bytes=" << populationDelta
              << " retained_delta_bytes=" << retainedDelta
              << " object_bytes_per_ringinfo=" << objectBytesPerRingInfo
              << " ring_bytes_per_ringinfo=" << ringBytesPerRingInfo
              << " total_bytes_per_ringinfo=" << totalBytesPerRingInfo
              << " rings=" << rings << " atom_members=" << members << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
