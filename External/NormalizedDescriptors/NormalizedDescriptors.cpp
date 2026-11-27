//
//  Copyright (C) 2026 Brian Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "NormalizedDescriptors.h"

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <locale>
#include <sstream>

namespace RDKit {
namespace NormalizedDescriptors {

CumulativeHistogram::CumulativeHistogram(std::vector<double> edges,
                                         std::vector<double> fractions)
    : d_edges(std::move(edges)), d_fractions(std::move(fractions)) {
  if (d_edges.empty()) {
    throw ValueErrorException("CumulativeHistogram requires at least one bin");
  }
  if (d_edges.size() != d_fractions.size()) {
    throw ValueErrorException(
        "CumulativeHistogram edge and fraction sizes differ");
  }
  if (!std::is_sorted(d_edges.begin(), d_edges.end())) {
    throw ValueErrorException("CumulativeHistogram edges must be sorted");
  }
}

double CumulativeHistogram::normalize(double value) const {
  if (d_edges.empty() || std::isnan(value)) {
    return 0.0;
  }
  // the number of edges < value, as python's bisect(bins, (value,))
  auto p =
      std::lower_bound(d_edges.begin(), d_edges.end(), value) - d_edges.begin();
  if (static_cast<size_t>(p) < d_fractions.size()) {
    return d_fractions[p];
  }
  return 1.0;
}

void HistogramSet::addTable(const std::string &name,
                            CumulativeHistogram table) {
  auto it = d_index.find(name);
  if (it != d_index.end()) {
    d_tables[it->second] = std::move(table);
    return;
  }
  auto idx = d_tables.size();
  d_tables.push_back(std::move(table));
  d_names.push_back(name);
  d_index[name] = idx;
}

bool HistogramSet::hasTable(const std::string &name) const {
  return d_index.find(name) != d_index.end();
}

const CumulativeHistogram &HistogramSet::getTable(
    const std::string &name) const {
  auto it = d_index.find(name);
  if (it == d_index.end()) {
    throw KeyErrorException(name);
  }
  return d_tables[it->second];
}

const CumulativeHistogram &HistogramSet::getTable(size_t idx) const {
  if (idx >= d_tables.size()) {
    throw IndexErrorException(static_cast<int>(idx));
  }
  return d_tables[idx];
}

int HistogramSet::getTableIndex(const std::string &name) const {
  auto it = d_index.find(name);
  return it == d_index.end() ? -1 : static_cast<int>(it->second);
}

double HistogramSet::normalize(const std::string &name, double value) const {
  auto it = d_index.find(name);
  if (it == d_index.end()) {
    return 0.0;
  }
  return d_tables[it->second].normalize(value);
}

double HistogramSet::normalize(size_t idx, double value) const {
  return getTable(idx).normalize(value);
}

namespace {
// returns false at end of stream, skips blank and comment lines
bool getDataLine(std::istream &inStream, std::string &line,
                 unsigned int &lineNum) {
  while (std::getline(inStream, line)) {
    ++lineNum;
    auto start = line.find_first_not_of(" \t\r\n");
    if (start != std::string::npos && line[start] != '#') {
      return true;
    }
  }
  return false;
}

[[noreturn]] void parseError(const std::string &msg, unsigned int lineNum) {
  std::ostringstream errout;
  errout << msg << " on line " << lineNum;
  throw ValueErrorException(errout.str());
}
}  // namespace

void HistogramSet::loadFromStream(std::istream &inStream) {
  std::string line;
  unsigned int lineNum = 0;
  while (getDataLine(inStream, line, lineNum)) {
    std::istringstream ls(line);
    ls.imbue(std::locale::classic());
    std::string keyword, name, extra;
    size_t nbins;
    if (!(ls >> keyword) || keyword != "histogram") {
      parseError("expected a 'histogram' line", lineNum);
    }
    if (!(ls >> name >> nbins) || nbins == 0 || (ls >> extra)) {
      parseError("bad histogram header", lineNum);
    }
    std::vector<double> edges(nbins), fractions(nbins);
    for (size_t i = 0; i < nbins; ++i) {
      if (!getDataLine(inStream, line, lineNum)) {
        parseError("histogram for " + name + " ends early", lineNum);
      }
      std::istringstream bs(line);
      bs.imbue(std::locale::classic());
      if (!(bs >> edges[i] >> fractions[i]) || (bs >> extra)) {
        parseError("bad histogram bin for " + name, lineNum);
      }
    }
    addTable(name, CumulativeHistogram(std::move(edges), std::move(fractions)));
  }
}

void HistogramSet::loadFromFile(const std::string &fileName) {
  std::ifstream inStream(fileName);
  if (!inStream) {
    throw BadFileException("could not open histogram table file " + fileName);
  }
  loadFromStream(inStream);
}

void HistogramSet::writeToStream(std::ostream &outStream) const {
  auto oldPrecision =
      outStream.precision(std::numeric_limits<double>::max_digits10);
  for (size_t i = 0; i < d_tables.size(); ++i) {
    const auto &table = d_tables[i];
    outStream << "histogram " << d_names[i] << " " << table.getEdges().size()
              << "\n";
    for (size_t j = 0; j < table.getEdges().size(); ++j) {
      outStream << table.getEdges()[j] << " " << table.getFractions()[j]
                << "\n";
    }
  }
  outStream.precision(oldPrecision);
}

std::string getDefaultTablePath() {
  const char *rdbase = std::getenv("RDBASE");
  std::string base = rdbase ? rdbase : "";
  return base +
         "/External/NormalizedDescriptors/data/normalized_descriptor_histograms.txt";
}

std::shared_ptr<const HistogramSet> defaultTablesPtr() {
  // thread-safe one-time initialization
  static const auto tables = [] {
    auto res = std::make_shared<HistogramSet>();
    res->loadFromFile(getDefaultTablePath());
    return std::shared_ptr<const HistogramSet>(std::move(res));
  }();
  return tables;
}

const HistogramSet &getDefaultTables() { return *defaultTablesPtr(); }

}  // namespace NormalizedDescriptors
}  // namespace RDKit
