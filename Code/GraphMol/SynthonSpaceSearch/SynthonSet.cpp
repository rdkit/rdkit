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

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
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

void SynthonSet::readFromDBStream(std::istream &is) {
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
    for (size_t i = d_synthons.size();
         i < static_cast<size_t>(synthonSetNum) + 1; ++i) {
      d_synthons.emplace_back();
    }
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
      remove_if(d_synthons.begin(), d_synthons.end(),
                [](const std::vector<std::unique_ptr<Synthon>> &r) -> bool {
                  return r.empty();
                }),
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

namespace {
// Take the synthons and build molecules from them.  longVecNum is the number
// of the vector containing the synthon set of interest.  The other members
// of synthon are assumed to be a single molecule, and each product is
// a molecule made from the corresponding member of longVecNum and the first
// element of the other vectors.

std::vector<std::unique_ptr<ROMol>> buildSampleMolecules(
    const std::vector<std::vector<ROMol *>> &synthons, size_t longVecNum) {
  std::vector<std::unique_ptr<ROMol>> sampleMolecules;
  sampleMolecules.reserve(synthons[longVecNum].size());

  MolzipParams mzparams;
  mzparams.label = MolzipLabel::Isotope;

  for (size_t i = 0; i < synthons[longVecNum].size(); ++i) {
    auto combMol = std::make_unique<ROMol>();
    for (size_t j = 0; j < synthons.size(); ++j) {
      if (j == longVecNum) {
        combMol.reset(combineMols(*combMol, *synthons[j][i]));
      } else {
        combMol.reset(combineMols(*combMol, *synthons[j].front()));
      }
    }
    auto sampleMol = molzip(*combMol, mzparams);
    MolOps::sanitizeMol(*dynamic_cast<RWMol *>(sampleMol.get()));
    sampleMolecules.push_back(std::move(combMol));
  }

  return sampleMolecules;
}
}  // namespace

void SynthonSet::buildSynthonFingerprints(
    const std::unique_ptr<FingerprintGenerator<std::uint64_t>> &fpGenerator) {
  // It's not as straightforward as making a fingerprint for each of the
  // synthons in turn because sometimes 2 aliphatic synthons come together
  // to form an aromatic ring, for example.  Such as:
  // [1*]=CC=C[2*] and [1*]Nc1c([2*])cccc1 giving c1ccc2ncccc2c1.  Each
  // synthon is built into a product and the atoms and bonds tracked.
  // Properties of the atoms and bonds are mapped back from the products
  // onto the synthons.  This way, when a molecule is split up, the
  // fingerprints make from the fragments should be consistent with those
  // of the corresponding synthons.

  d_synthonFPs.clear();

  tagSynthonAtomsAndBonds();

  // Now build sets of sample molecules using each synthon set in turn.
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    d_synthonFPs.emplace_back();
    d_synthonFPs[i].reserve(d_synthons[i].size());
    std::vector<boost::dynamic_bitset<>> synthsToUse(d_synthons.size());
    for (size_t j = 0; j < d_synthons.size(); ++j) {
      if (j == i) {
        synthsToUse[j] = boost::dynamic_bitset<>(d_synthons[j].size(), true);
      } else {
        synthsToUse[j] = boost::dynamic_bitset<>(d_synthons[j].size(), false);
        synthsToUse[j][0] = true;
      }
    }
    auto synthons = getSynthons(synthsToUse);
    auto sampleMols = buildSampleMolecules(synthons, i);
    makeSynthonFPs(i, sampleMols, fpGenerator);
  }
}

std::vector<std::vector<ROMol *>> SynthonSet::getSynthons(
    const std::vector<boost::dynamic_bitset<>> &reqSynths) const {
  PRECONDITION(reqSynths.size() == d_synthons.size(), "wrong synthon details");
  for (size_t i = 0; i < reqSynths.size(); ++i) {
    PRECONDITION(reqSynths[i].size() == d_synthons[i].size(),
                 "wrong synthon details");
  }
  std::vector<std::vector<ROMol *>> synthons(d_synthons.size(),
                                             std::vector<ROMol *>());
  for (size_t i = 0; i < reqSynths.size(); ++i) {
    synthons[i].reserve(reqSynths[i].count());
    for (size_t j = 0; j < reqSynths[i].size(); ++j) {
      if (reqSynths[i][j]) {
        synthons[i].push_back(d_synthons[i][j]->getMol().get());
      }
    }
  }
  return synthons;
}

void SynthonSet::tagSynthonAtomsAndBonds() {
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    for (auto &syn : d_synthons[i]) {
      syn->tagAtomsAndBonds(static_cast<int>(i));
    }
  }
}

namespace {
// transfer information for bonds created by molzip
void fixSynthonAtomAndBond(const Atom *sampleMolAtom, const Bond *bond,
                           RWMol &synthCp) {
  auto synthAt = synthCp.getAtomWithIdx(sampleMolAtom->getProp<int>("idx"));
  for (auto nbor : synthCp.atomNeighbors(synthAt)) {
    if (!nbor->getAtomicNum() && nbor->getIsotope() <= MAX_CONNECTOR_NUM) {
      nbor->setIsAromatic(sampleMolAtom->getIsAromatic());
      auto synthBond =
          synthCp.getBondBetweenAtoms(synthAt->getIdx(), nbor->getIdx());
      synthBond->setIsAromatic(bond->getIsAromatic());
      synthBond->setBondType(bond->getBondType());
    }
  }
}
}  // namespace

void SynthonSet::makeSynthonFPs(
    size_t synthSetNum, const std::vector<std::unique_ptr<ROMol>> &sampleMols,
    const std::unique_ptr<FingerprintGenerator<std::uint64_t>> &fpGenerator) {
  for (size_t i = 0; i < sampleMols.size(); ++i) {
    RWMol synthCp(*d_synthons[synthSetNum][i]->getMol());
    // transfer the aromaticity of the atom in the sample molecule to the
    // corresponding atom in the synthon
    for (const auto &atom : sampleMols[i]->atoms()) {
      int molNum = atom->getProp<int>("molNum");
      if (static_cast<size_t>(molNum) == synthSetNum) {
        int atIdx = atom->getProp<int>("idx");
        synthCp.getAtomWithIdx(atIdx)->setIsAromatic(atom->getIsAromatic());
      }
    }
    // likewise for the bonds, but not all bonds will have the tags,
    // because some are formed by molzip.
    for (const auto &bond : sampleMols[i]->bonds()) {
      if (bond->hasProp("molNum")) {
        int molNum = bond->getProp<int>("molNum");
        if (static_cast<size_t>(molNum) == synthSetNum) {
          int bondIdx = bond->getProp<int>("idx");
          auto sbond = synthCp.getBondWithIdx(bondIdx);
          sbond->setIsAromatic(bond->getIsAromatic());
          sbond->setBondType(bond->getBondType());
        }
      } else {
        // it came from molzip, so in the synth, one end atom corresponds to
        // an atom in sampleMol and the other end was a dummy.
        fixSynthonAtomAndBond(bond->getBeginAtom(), bond, synthCp);
        fixSynthonAtomAndBond(bond->getEndAtom(), bond, synthCp);
      }
    }
    d_synthonFPs[synthSetNum].emplace_back(
        fpGenerator->getFingerprint(synthCp));
  }
}

}  // namespace RDKit::SynthonSpaceSearch