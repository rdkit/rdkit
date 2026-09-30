// Copyright (c) 2026, RDKit contributors
// All rights reserved.
//
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the root of the RDKit source tree.

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace {

constexpr std::size_t kDefaultCount = 135'000;

std::string defaultInputPath() {
  const auto *home = std::getenv("HOME");
  if (!home) {
    return "data/enamine_scrambled.cxsmiles";
  }
  return std::string(home) + "/data/enamine_scrambled.cxsmiles";
}

std::size_t parseCount(const char *value) {
  try {
    const auto text = std::string(value);
    std::size_t consumed = 0;
    const auto parsed = std::stoull(text, &consumed);
    if (consumed != text.size()) {
      throw std::invalid_argument("trailing characters");
    }
    return parsed;
  } catch (const std::exception &) {
    throw std::invalid_argument(std::string("invalid molecule count: ") +
                                value);
  }
}

void printUsage(const char *program) {
  std::cerr << "Usage: " << program << " [input.cxsmiles] [count]\n"
            << "Reads complete lines in file order, parses CXSMILES, then "
               "sanitizes each molecule.\n"
            << "Defaults: " << defaultInputPath() << " " << kDefaultCount
            << "\n";
}

}  // namespace

int main(int argc, char *argv[]) {
  if (argc > 3) {
    printUsage(argv[0]);
    return 2;
  }
  const auto inputPath = argc >= 2 ? argv[1] : defaultInputPath();
  RDKit::Chirality::setUseLegacyStereoPerception(false);
  std::size_t requested = kDefaultCount;
  try {
    if (argc == 3) {
      requested = parseCount(argv[2]);
    }
  } catch (const std::invalid_argument &error) {
    std::cerr << error.what() << '\n';
    printUsage(argv[0]);
    return 2;
  }
  std::ifstream input(inputPath);
  if (!input) {
    std::cerr << "cannot open input: " << inputPath << '\n';
    return 2;
  }
  RDKit::v2::SmilesParse::SmilesParserParams parserParams;
  parserParams.sanitize = false;
  parserParams.allowCXSMILES = true;
  parserParams.strictCXSMILES = true;
  parserParams.parseName = true;
  std::size_t read = 0, parsed = 0, sanitized = 0, failures = 0, atoms = 0;
  std::string line;
  const auto start = std::chrono::steady_clock::now();
  while (read < requested && std::getline(input, line)) {
    const auto tab = line.find('\t');
    const auto smiles = line.substr(0, tab);
    if (read == 0 && smiles == "smiles") {
      continue;
    }
    ++read;
    if (smiles.empty()) {
      ++failures;
      continue;
    }
    try {
      auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles, parserParams);
      if (!mol) {
        ++failures;
        continue;
      }
      ++parsed;
      RDKit::MolOps::sanitizeMol(*mol);
      ++sanitized;
      atoms += mol->getNumAtoms();
    } catch (const std::exception &) {
      ++failures;
    }
  }
  const auto elapsed = std::chrono::steady_clock::now() - start;
  const auto seconds = std::chrono::duration<double>(elapsed).count();
  std::cout << "input=" << inputPath << " requested=" << requested
            << " read=" << read << " parsed=" << parsed
            << " sanitized=" << sanitized << " failures=" << failures
            << " atoms=" << atoms << " elapsed_s=" << seconds << '\n';
  return read == requested ? 0 : 1;
}
