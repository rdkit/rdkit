//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file contains an implementation of synthonspace substructure search
// similar to that described in
// 'Fast Substructure Search in Combinatorial Library Spaces',
// Thomas Liphardt and Thomas Sander,
// J. Chem. Inf. Model. 2023, 63, 16, 5133–5141
// https://doi.org/10.1021/acs.jcim.3c00290
//
// The algorithm allows the substructure searching of a very large library
// of structures that is described in synthon format (such as Enamine REAL)
// without enumerating the individual structures during the search process.
//
// It is not a direct implementation of the published algorithm, as,
// for example, it uses a different fingerprint for the initial synthon
// screening.

#include <regex>

#include <GraphMol/MolPickler.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit::SynthonSpaceSearch {

const std::vector<std::shared_ptr<ROMol>> &SynthonSet::getConnectorRegions()
    const {
  return d_connectorRegions;
}

const std::unique_ptr<ExplicitBitVect> &SynthonSet::getConnRegFP() const {
  return d_connRegFP;
}

void SynthonSet::writeToDBStream(std::ostream &os) const {
  streamWrite(os, d_id);
  streamWrite(os, getConnectorRegions().size());
  for (const auto &cr : getConnectorRegions()) {
    MolPickler::pickleMol(*cr, os, PicklerOps::AllProps);
  }
  auto connRegFPstr = getConnRegFP()->toString();
  streamWrite(os, connRegFPstr);
  streamWrite(os, d_connectors.size());
  for (size_t i = 0; i < d_connectors.size(); ++i) {
    if (d_connectors[i]) {
      streamWrite(os, true);
    } else {
      streamWrite(os, false);
    }
  }
  streamWrite(os, d_synthons.size());
  for (const auto &rs : d_synthons) {
    streamWrite(os, rs.size());
    for (const auto &r : rs) {
      r->writeToDBStream(os);
    }
  }
}

void SynthonSet::readFromDBStream(std::istream &is, std::uint32_t) {
  streamRead(is, d_id, 0);
  size_t numConnRegs;
  streamRead(is, numConnRegs);
  d_connectorRegions.resize(numConnRegs);
  for (size_t i = 0; i < numConnRegs; ++i) {
    d_connectorRegions[i] = std::make_unique<ROMol>();
    MolPickler::molFromPickle(is, *d_connectorRegions[i]);
  }
  std::string pickle;
  streamRead(is, pickle, 0);
  d_connRegFP = std::make_unique<ExplicitBitVect>(pickle);
  size_t connSize;
  streamRead(is, connSize);
  d_connectors.resize(connSize);
  bool s;
  for (size_t i = 0; i < connSize; ++i) {
    streamRead(is, s);
    d_connectors[i] = s;
  }
  size_t numRS;
  streamRead(is, numRS);
  d_synthons.clear();
  for (size_t i = 0; i < numRS; ++i) {
    size_t numR;
    streamRead(is, numR);
    d_synthons.emplace_back();
    d_synthons[i].resize(numR);
    for (size_t j = 0; j < numR; ++j) {
      d_synthons[i][j] = std::make_unique<Synthon>();
      d_synthons[i][j]->readFromDBStream(is);
    }
  }
}

void SynthonSet::addSynthon(int synthonSetNum,
                            std::unique_ptr<Synthon> newSynthon) {
  if (static_cast<size_t>(synthonSetNum) >= d_synthons.size()) {
    d_synthons.resize(synthonSetNum + 1);
  }
  d_synthons[synthonSetNum].emplace_back(std::move(newSynthon));
}

void SynthonSet::assignConnectorsUsed() {
  // Find instances of "[1*]", "[1*:1]", "[2*]", "[2*:2]" etc.
  // and set d_connectors accordingly.
  static std::vector<std::regex> connRegexs;
  if (connRegexs.empty()) {
    for (size_t i = 0; i < MAX_CONNECTOR_NUM; ++i) {
      connRegexs.emplace_back(R"(\[)" + std::to_string(i + 1) + R"(\*\])");
      connRegexs.emplace_back(R"(\[)" + std::to_string(i + 1) + R"(\*\:)" +
                              std::to_string(i + 1) + R"(\])");
    }
  }
  d_connectors.resize(MAX_CONNECTOR_NUM, false);
  for (auto &reagSet : d_synthons) {
    for (auto &reag : reagSet) {
      for (size_t i = 0; i < MAX_CONNECTOR_NUM; ++i) {
        if (std::regex_search(reag->getSmiles(), connRegexs[2 * i]) ||
            std::regex_search(reag->getSmiles(), connRegexs[2 * i + 1])) {
          d_connectors.set(i);
        }
      }
    }
  }
}

const std::vector<int> &SynthonSet::getNumConnectors() const {
  return d_numConnectors;
}

void SynthonSet::buildConnectorRegions() {
  // Remove any empty reagent sets.  This most often happens if the
  // synthon number in the reaction starts from 1, not 0.  Idorsia
  // use 0, the stuff from ChemSpace uses 1.
  d_synthons.erase(
      std::remove_if(d_synthons.begin(), d_synthons.end(),
                     [](const std::vector<std::unique_ptr<Synthon>> &r)
                         -> bool { return r.empty(); }),
      d_synthons.end());

  std::set<std::string> smis;
  for (const auto &rset : d_synthons) {
    for (const auto &r : rset) {
      for (const auto &cr : r->getConnRegions()) {
        auto smi = MolToSmiles(*cr);
        if (smis.insert(smi).second) {
          d_connectorRegions.push_back(cr);
        }
      }
    }
  }
  if (!getConnectorRegions().empty()) {
    d_connRegFP.reset(PatternFingerprintMol(*getConnectorRegions().front()));
    for (size_t i = 1; i < getConnectorRegions().size(); ++i) {
      std::unique_ptr<ExplicitBitVect> fp(
          PatternFingerprintMol(*getConnectorRegions()[i]));
      *d_connRegFP |= *fp;
    }
  } else {
    d_connRegFP = std::make_unique<ExplicitBitVect>(2048);
  }
  // It should be the case that all synthons in a synthon set
  // have the same number of connections, so just do the 1st
  // one of each.
  for (const auto &synthonSet : d_synthons) {
    d_numConnectors.push_back(
        details::countConnections(synthonSet.front()->getSmiles()));
  }
}
}  // namespace RDKit::SynthonSpaceSearch