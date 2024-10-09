//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <regex>

#include <GraphMol/MolPickler.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/HyperspaceSearch/ReactionSet.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/StreamOps.h>

namespace RDKit::HyperspaceSearch {

const std::vector<std::shared_ptr<ROMol>> &ReactionSet::connectorRegions()
    const {
  if (d_connectorRegions.empty()) {
    std::set<std::string> smis;
    for (const auto &rset : d_reagents) {
      for (const auto &r : rset) {
        for (const auto &cr : r->connRegions()) {
          auto smi = MolToSmiles(*cr);
          if (smis.insert(smi).second) {
            d_connectorRegions.push_back(cr);
          }
        }
      }
    }
  }
  return d_connectorRegions;
}

const std::unique_ptr<ExplicitBitVect> &ReactionSet::connRegFP() const {
  if (!d_connRegFP) {
    if (!connectorRegions().empty()) {
      d_connRegFP.reset(PatternFingerprintMol(*connectorRegions().front()));
      for (size_t i = 1; i < connectorRegions().size(); ++i) {
        std::unique_ptr<ExplicitBitVect> fp(
            PatternFingerprintMol(*connectorRegions()[i]));
        *d_connRegFP |= *fp;
      }
    } else {
      d_connRegFP.reset(new ExplicitBitVect(2048));
    }
  }
  return d_connRegFP;
}

void ReactionSet::writeToDBStream(std::ostream &os) const {
  streamWrite(os, d_id);
  streamWrite(os, connectorRegions().size());
  for (const auto &cr : connectorRegions()) {
    MolPickler::pickleMol(*cr, os, PicklerOps::AllProps);
  }
  auto connRegFPstr = connRegFP()->toString();
  streamWrite(os, connRegFPstr);
  streamWrite(os, d_connectors.size());
  for (size_t i = 0; i < d_connectors.size(); ++i) {
    if (d_connectors[i]) {
      streamWrite(os, true);
    } else {
      streamWrite(os, false);
    }
  }
  streamWrite(os, d_reagents.size());
  for (const auto &rs : d_reagents) {
    streamWrite(os, rs.size());
    for (const auto &r : rs) {
      r->writeToDBStream(os);
    }
  }
}

void ReactionSet::readFromDBStream(std::istream &is) {
  streamRead(is, d_id, 0);
  size_t numConnRegs = 3;
  streamRead(is, numConnRegs);
  d_connectorRegions.resize(numConnRegs);
  for (size_t i = 0; i < numConnRegs; ++i) {
    d_connectorRegions[i] = std::unique_ptr<ROMol>(new ROMol);
    MolPickler::molFromPickle(is, *d_connectorRegions[i]);
  }
  std::string pickle;
  streamRead(is, pickle, 0);
  d_connRegFP.reset(new ExplicitBitVect(pickle));
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
  d_reagents.clear();
  for (size_t i = 0; i < numRS; ++i) {
    size_t numR;
    streamRead(is, numR);
    d_reagents.push_back(std::vector<std::unique_ptr<Reagent>>());
    for (size_t j = 0; j < numR; ++j) {
      d_reagents[i].emplace_back(new Reagent);
      d_reagents[i][j]->readFromDBStream(is);
    }
  }
}

void ReactionSet::addReagent(int reagentSetNum, const std::string &smiles,
                             const std::string &reagentId) {
  if (static_cast<size_t>(reagentSetNum) >= d_reagents.size()) {
    for (size_t i = d_reagents.size();
         i < static_cast<size_t>(reagentSetNum) + 1; ++i) {
      d_reagents.push_back(std::vector<std::unique_ptr<Reagent>>());
    }
  }
  d_reagents[reagentSetNum].emplace_back(new Reagent(smiles, reagentId));
}

void ReactionSet::assignConnectorsUsed() {
  const static std::vector<std::regex> connRegexs{
      std::regex(R"(\[1\*\])"), std::regex(R"(\[2\*\])"),
      std::regex(R"(\[3\*\])"), std::regex(R"(\[4\*\])")};
  d_connectors.resize(4, false);
  // Remove any empty reagent sets.  This most often happens if the
  // synthon number in the reaction starts from 1, not 0.  Idorsia
  // use 0, the stuff from ChemSpace uses 1.
  d_reagents.erase(
      remove_if(d_reagents.begin(), d_reagents.end(),
                [&](const std::vector<std::unique_ptr<Reagent>> &r) -> bool {
                  return r.empty();
                }),
      d_reagents.end());
  for (auto &reagSet : d_reagents) {
    for (auto &reag : reagSet) {
      for (size_t i = 0; i < 4; ++i) {
        if (std::regex_search(reag->smiles(), connRegexs[i])) {
          d_connectors.set(i);
        }
      }
    }
  }
}

}  // namespace RDKit::HyperspaceSearch