//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// The algorithm allows the substructure searching of a very large library
// of structures that is described in synthon format (such as Enamine REAL)
// without enumerating the individual structures during the search process.
//
// It is not a direct implementation of the published algorithm, as,
// for example, it uses a different fingerprint for the initial synthon
// screening.

#include <cmath>
#include <numeric>
#include <random>
#include <regex>
#include <unordered_map>

#include <../External/pubchem_shape/PubChemShape.hpp>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/CIPLabeler/Descriptor.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/RDThreads.h>

namespace RDKit::SynthonSpaceSearch {

const std::vector<std::shared_ptr<ROMol>> &SynthonSet::getConnectorRegions()
    const {
  return d_connectorRegions;
}
const std::vector<std::string> &SynthonSet::getConnectorRegionSmiles() const {
  return d_connRegSmis;
}
const std::vector<std::unique_ptr<ExplicitBitVect>> &SynthonSet::getConnRegFPs()
    const {
  return d_connRegFPs;
}
const std::unique_ptr<ExplicitBitVect> &SynthonSet::getAddFP() const {
  return d_addFP;
}
const std::unique_ptr<ExplicitBitVect> &SynthonSet::getSubtractFP() const {
  return d_subtractFP;
}

const std::pair<std::string, Synthon *>
SynthonSet::getSubstructureOrderedSynthon(size_t setNum,
                                          size_t synthonNum) const {
  if (setNum >= d_synthons.size() || synthonNum >= d_synthons[setNum].size()) {
    return std::make_pair("", nullptr);
  }
  return d_synthons[setNum][d_substructureSearchOrders[setNum][synthonNum]];
}

const std::pair<std::string, Synthon *>
SynthonSet::getFingerprintOrderedSynthon(size_t setNum,
                                         size_t synthonNum) const {
  if (setNum >= d_synthons.size() || synthonNum >= d_synthons[setNum].size()) {
    return std::make_pair("", nullptr);
  }
  return d_synthons[setNum][d_fingerprintSearchOrders[setNum][synthonNum]];
}

const std::pair<std::string, Synthon *> SynthonSet::getRascalOrderedSynthon(
    size_t setNum, size_t synthonNum) const {
  if (setNum >= d_synthons.size() || synthonNum >= d_synthons[setNum].size()) {
    return std::make_pair("", nullptr);
  }
  return d_synthons[setNum][d_rascalSearchOrders[setNum][synthonNum]];
}

size_t SynthonSet::getSubstructureOrderedSynthonNum(size_t setNum,
                                                    size_t synthonNum) const {
  if (setNum >= d_synthons.size() || synthonNum >= d_synthons[setNum].size()) {
    // Obviously this will be a very large number.
    return -1;
  }
  return d_substructureSearchOrders[setNum][synthonNum];
}

size_t SynthonSet::getFingerprintOrderedSynthonNum(size_t setNum,
                                                   size_t synthonNum) const {
  if (setNum >= d_synthons.size() || synthonNum >= d_synthons[setNum].size()) {
    // Obviously this will be a very large number.
    return -1;
  }
  return d_fingerprintSearchOrders[setNum][synthonNum];
}

size_t SynthonSet::getRascalOrderedSynthonNum(size_t setNum,
                                              size_t synthonNum) const {
  if (setNum >= d_synthons.size() || synthonNum >= d_synthons[setNum].size()) {
    // Obviously this will be a very large number.
    return -1;
  }
  return d_rascalSearchOrders[setNum][synthonNum];
}

namespace {
void writeBitSet(std::ostream &os, const boost::dynamic_bitset<> &bitset) {
  streamWrite(os, bitset.size());
  for (unsigned int i = 0; i < bitset.size(); ++i) {
    if (bitset[i]) {
      streamWrite(os, true);
    } else {
      streamWrite(os, false);
    }
  }
}
void writeSearchOrderSet(std::ostream &os,
                         const std::vector<std::vector<size_t>> &orderSet) {
  streamWrite(os, static_cast<std::uint64_t>(orderSet.size()));
  for (unsigned int i = 0; i < orderSet.size(); ++i) {
    streamWrite(os, static_cast<std::uint64_t>(orderSet[i].size()));
    for (const auto v : orderSet[i]) {
      streamWrite(os, static_cast<std::uint64_t>(v));
    }
  }
}
}  // namespace

void SynthonSet::writeToDBStream(std::ostream &os) const {
  streamWrite(os, d_id);
  streamWrite(os, static_cast<std::uint64_t>(getConnectorRegions().size()));
  for (const auto &cr : getConnectorRegions()) {
    MolPickler::pickleMol(*cr, os, PicklerOps::AllProps);
  }
  streamWrite(os, static_cast<std::uint64_t>(getConnRegFPs().size()));
  for (const auto &fp : getConnRegFPs()) {
    auto connRegFPstr = fp->toString();
    streamWrite(os, connRegFPstr);
  }

  writeBitSet(os, d_connectors);
  streamWrite(os, static_cast<std::uint64_t>(d_synthConnPatts.size()));
  for (const auto &scp : d_synthConnPatts) {
    writeBitSet(os, scp);
  }

  streamWrite(os, static_cast<std::uint64_t>(d_synthons.size()));
  for (const auto &rs : d_synthons) {
    streamWrite(os, static_cast<std::uint64_t>(rs.size()));
    for (const auto &[id, synthon] : rs) {
      streamWrite(os, id);
      streamWrite(os, synthon->getSmiles());
    }
  }

  if (d_addFP) {
    streamWrite(os, true);
    streamWrite(os, d_addFP->toString());
    streamWrite(os, d_subtractFP->toString());
  } else {
    streamWrite(os, false);
  }
  streamWrite(os, static_cast<std::uint64_t>(d_numConnectors.size()));
  for (const auto &c : d_numConnectors) {
    streamWrite(os, c);
  }
  writeSearchOrderSet(os, d_substructureSearchOrders);
  writeSearchOrderSet(os, d_fingerprintSearchOrders);
  writeSearchOrderSet(os, d_rascalSearchOrders);
}

namespace {
void readBitSet(std::istream &is, boost::dynamic_bitset<> &bitset) {
  size_t bsSize;
  streamRead(is, bsSize);
  bitset.resize(bsSize);
  bool s;
  for (size_t i = 0; i < bsSize; ++i) {
    streamRead(is, s);
    bitset[i] = s;
  }
}
void readSearchOrderSet(std::istream &is,
                        std::vector<std::vector<size_t>> &orderSet) {
  std::uint64_t s;
  streamRead(is, s);
  orderSet.resize(s);
  for (size_t i = 0; i < s; ++i) {
    std::uint64_t t;
    streamRead(is, t);
    orderSet[i].resize(t);
    for (size_t j = 0; j < t; ++j) {
      std::uint64_t u;
      streamRead(is, u);
      orderSet[i][j] = u;
    }
  }
}
}  // namespace

void SynthonSet::readFromDBStream(std::istream &is, const SynthonSpace &space,
                                  std::uint32_t version) {
  PRECONDITION(version >= 3000, "Binary database version no longer supported.");
  streamRead(is, d_id, 0);
  std::uint64_t numConnRegs;
  streamRead(is, numConnRegs);
  d_connectorRegions.resize(numConnRegs);
  d_connRegSmis.resize(numConnRegs);
  for (std::uint64_t i = 0; i < numConnRegs; ++i) {
    d_connectorRegions[i] = std::make_unique<ROMol>();
    MolPickler::molFromPickle(is, *d_connectorRegions[i]);
    d_connRegSmis[i] = MolToSmiles(*d_connectorRegions[i]);
  }
  std::uint64_t numConnRegFPs;
  streamRead(is, numConnRegFPs);
  for (std::uint64_t i = 0; i < numConnRegFPs; ++i) {
    std::string pickle;
    streamRead(is, pickle, 0);
    d_connRegFPs.push_back(std::make_unique<ExplicitBitVect>(pickle));
  }
  readBitSet(is, d_connectors);
  std::uint64_t numSynthConnPatts;
  streamRead(is, numSynthConnPatts);
  d_synthConnPatts.resize(numSynthConnPatts);
  for (std::uint64_t i = 0; i < numSynthConnPatts; ++i) {
    boost::dynamic_bitset<> synthConnPatt;
    readBitSet(is, synthConnPatt);
    d_synthConnPatts[i] = synthConnPatt;
  }

  std::uint64_t numRS;
  streamRead(is, numRS);
  d_synthons.clear();
  for (std::uint64_t i = 0; i < numRS; ++i) {
    std::uint64_t numR;
    streamRead(is, numR);
    d_synthons.emplace_back();
    d_synthons[i].resize(numR);
    for (std::uint64_t j = 0; j < numR; ++j) {
      std::string synthonName;
      streamRead(is, synthonName, 0);
      std::string smiles;
      streamRead(is, smiles, 0);
      // The synthons should be in the pool by now.
      if (auto synth = space.getSynthonFromPool(smiles); synth == nullptr) {
        std::cout << "smiles " << smiles << " not in pool" << std::endl;
        throw std::runtime_error("Database file " + space.getInputFileName() +
                                 " appears corrupted.");
      } else {
        d_synthons[i][j] = std::pair(synthonName, synth);
      }
    }
  }

  bool haveAddFP;
  streamRead(is, haveAddFP);
  if (haveAddFP) {
    std::string fString;
    streamRead(is, fString, 0);
    d_addFP = std::make_unique<ExplicitBitVect>(fString);
    streamRead(is, fString, 0);
    d_subtractFP = std::make_unique<ExplicitBitVect>(fString);
  }
  std::uint64_t numConns;
  streamRead(is, numConns);
  d_numConnectors.resize(numConns);
  for (unsigned int i = 0; i < numConns; ++i) {
    streamRead(is, d_numConnectors[i]);
  }

  if (version < 3020) {
    initializeSearchOrders();
  } else {
    readSearchOrderSet(is, d_substructureSearchOrders);
    readSearchOrderSet(is, d_fingerprintSearchOrders);
    readSearchOrderSet(is, d_rascalSearchOrders);
  }
}

void SynthonSet::enumerateToStream(std::ostream &os) const {
  std::vector<size_t> numSynthons;
  for (const auto &synths : d_synthons) {
    numSynthons.push_back(synths.size());
  }
  details::Stepper stepper(numSynthons);
  while (stepper.d_currState[0] != numSynthons[0]) {
    auto prod = buildProduct(stepper.d_currState);
    auto prodName = buildProductName(stepper.d_currState);
    os << MolToSmiles(*prod) << " " << prodName << std::endl;
    stepper.step();
  }
}

void SynthonSet::addSynthon(const int synthonSetNum, Synthon *newSynthon,
                            const std::string &synthonId) {
  if (std::cmp_greater_equal(synthonSetNum, d_synthons.size())) {
    d_synthons.resize(synthonSetNum + 1);
  }
  d_synthons[synthonSetNum].push_back(std::make_pair(synthonId, newSynthon));
}

std::vector<std::unique_ptr<ROMol>> SynthonSet::buildSampleMolecules(
    const std::vector<std::vector<std::unique_ptr<RWMol>>> &synthonMols,
    size_t longVecNum) const {
  std::vector<std::unique_ptr<ROMol>> sampleMolecules;
  sampleMolecules.reserve(d_synthons[longVecNum].size());

  // Find a small synthon in each of the non-longVecNum sets.
  std::vector<size_t> exSynthons(synthonMols.size(), 0);
  for (size_t i = 0; i < synthonMols.size(); ++i) {
    if (i == longVecNum) {
      continue;
    }
    unsigned int minAtoms = synthonMols[i][0]->getNumAtoms();
    for (size_t j = 1; j < synthonMols[i].size(); ++j) {
      if (synthonMols[i][j]->getNumAtoms() < minAtoms) {
        minAtoms = synthonMols[i][j]->getNumAtoms();
        exSynthons[i] = j;
      }
    }
  }

  MolzipParams mzparams;
  mzparams.label = MolzipLabel::Isotope;

  for (size_t i = 0; i < d_synthons[longVecNum].size(); ++i) {
    exSynthons[longVecNum] = i;
    auto sampleMol = buildMolecule(exSynthons);
    sampleMolecules.push_back(std::move(sampleMol));
  }
  return sampleMolecules;
}

namespace {
int findMolNumFrag(const std::vector<std::unique_ptr<ROMol>> &molFrags,
                   size_t molNum) {
  for (size_t i = 0; i < molFrags.size(); ++i) {
    for (const auto &atom : molFrags[i]->atoms()) {
      if (atom->hasProp("molNum") &&
          atom->getProp<unsigned int>("molNum") == molNum) {
        return i;
      }
    }
  }
  return -1;
}
}  // namespace

void SynthonSet::makeSynthonSearchMols() {
  // Each synthon is built into a product and then split on the
  // newly formed bonds.  The fragment from the original synthon of
  // interest is then stored as the synthon searchMol.  This way,
  // when a query molecule is split up, the fragments should be
  // consistent with those of the corresponding synthons.

  auto synthonMolCopies = copySynthons();
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  for (unsigned int i = 1; i <= MAX_CONNECTOR_NUM; ++i) {
    dummyLabels.emplace_back(i, i);
  }

  // Now build sets of sample molecules using each synthon set in turn.
  for (size_t synthSetNum = 0; synthSetNum < d_synthons.size(); ++synthSetNum) {
    auto sampleMols = buildSampleMolecules(synthonMolCopies, synthSetNum);
    for (size_t j = 0; j < sampleMols.size(); ++j) {
      std::vector<unsigned int> splitBonds;
      for (const auto &bond : sampleMols[j]->bonds()) {
        if (!bond->hasProp("molNum")) {
          splitBonds.push_back(bond->getIdx());
        }
      }

      const std::unique_ptr<ROMol> fragMol(MolFragmenter::fragmentOnBonds(
          *sampleMols[j], splitBonds, true, &dummyLabels));
      std::vector<std::unique_ptr<ROMol>> molFrags;
      MolOps::getMolFrags(*fragMol, molFrags, false);
      int fragWeWant = findMolNumFrag(molFrags, synthSetNum);
      unsigned int otf;
      sanitizeMol(*static_cast<RWMol *>(molFrags[fragWeWant].get()), otf,
                  MolOps::SANITIZE_SYMMRINGS);
      d_synthons[synthSetNum][j].second->setSearchMol(
          std::move(molFrags[fragWeWant]));
    }
  }
}

void SynthonSet::removeEmptySynthonSets() {
  d_synthons.erase(
      std::remove_if(d_synthons.begin(), d_synthons.end(),
                     [](const std::vector<std::pair<std::string, Synthon *>> &r)
                         -> bool { return r.empty(); }),
      d_synthons.end());
}

void SynthonSet::buildConnectorRegions() {
  std::set<std::string> smis;
  for (const auto &rset : d_synthons) {
    for (const auto &r : rset) {
      for (const auto &cr : r.second->getConnRegions()) {
        if (auto smi = MolToSmiles(*cr); smis.insert(smi).second) {
          d_connectorRegions.push_back(cr);
          d_connRegSmis.push_back(smi);
        }
      }
    }
  }
  if (!getConnectorRegions().empty()) {
    d_connRegFPs.reserve(d_connectorRegions.size());
    for (const auto &connReg : getConnectorRegions()) {
      std::unique_ptr<ExplicitBitVect> fp(
          PatternFingerprintMol(*connReg, PATT_FP_NUM_BITS));
      d_connRegFPs.push_back(std::move(fp));
    }
  }
  // It should be the case that all synthons in a synthon set
  // have the same number of connections, so just do the 1st
  // one of each.
  for (const auto &synthonSet : d_synthons) {
    d_numConnectors.push_back(
        details::countConnections(*synthonSet.front().second->getOrigMol()));
  }
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
  d_connectors.resize(MAX_CONNECTOR_NUM + 1, false);
  d_synthConnPatts.clear();
  for (const auto &synthSet : d_synthons) {
    // We only need to look at the first synthon in each set, as they
    // should all be the same.
    d_synthConnPatts.emplace_back();
    d_synthConnPatts.back().resize(MAX_CONNECTOR_NUM + 1, false);
    const auto &reag = synthSet.front();
    for (size_t i = 0; i < MAX_CONNECTOR_NUM; ++i) {
      if (std::regex_search(reag.second->getSmiles(), connRegexs[2 * i]) ||
          std::regex_search(reag.second->getSmiles(), connRegexs[2 * i + 1])) {
        d_connectors.set(i + 1);
        d_synthConnPatts.back().set(i + 1);
      }
    }
  }
}

const std::vector<int> &SynthonSet::getNumConnectors() const {
  return d_numConnectors;
}

std::uint64_t SynthonSet::getNumProducts() const {
  uint64_t thisSize = 1;
  for (const auto &r : d_synthons) {
    thisSize *= r.size();
  }
  return thisSize;
}

bool SynthonSet::hasFingerprints() const {
  if (d_synthons.empty()) {
    return false;
  }
  return static_cast<bool>(d_synthons.front().front().second->getFP());
}

bool SynthonSet::hasAddAndSubtractFPs() const {
  return static_cast<bool>(d_addFP);
}

void SynthonSet::buildSynthonFingerprints(
    const FingerprintGenerator<std::uint64_t> &fpGen) {
  d_addFP.reset();
  d_subtractFP.reset();

  // The synthons should have had transferProductBondsToSynthons
  // applied to them by now, so that they have a searchMol.

  for (size_t synthSetNum = 0; synthSetNum < d_synthons.size(); ++synthSetNum) {
    for (const auto &synth : d_synthons[synthSetNum]) {
      if (ControlCHandler::getGotSignal()) {
        return;
      }
      synth.second->setFP(std::unique_ptr<ExplicitBitVect>(
          fpGen.getFingerprint(*synth.second->getSearchMol())));
    }
  }
}

void SynthonSet::buildAddAndSubtractFPs(
    const FingerprintGenerator<std::uint64_t> &fpGen, unsigned int numBits) {
  PRECONDITION(hasFingerprints(), "No fingerprints for synthons.");
  d_addFP.reset();
  d_subtractFP.reset();
  std::vector<std::vector<size_t>> synthonNums(d_synthons.size());
  std::vector<size_t> numSynthons(d_synthons.size());
  std::vector<int> naddbitcounts(numBits, 0);
  std::vector<int> nsubbitcounts(numBits, 0);
  size_t totSamples = 1;
  // Sample the synthons evenly across their size ranges.
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    std::vector<std::pair<size_t, std::pair<std::string, Synthon *>>>
        sortedSynthons(d_synthons[i].size());
    for (size_t j = 0; j < d_synthons[i].size(); ++j) {
      sortedSynthons[j] = std::make_pair(j, d_synthons[i][j]);
    }
    std::sort(sortedSynthons.begin(), sortedSynthons.end(),
              [](const std::pair<size_t, std::pair<std::string, Synthon *>> &a,
                 const std::pair<size_t, std::pair<std::string, Synthon *>> &b)
                  -> bool {
                auto as = a.second.second;
                auto bs = b.second.second;
                if (as->getOrigMol()->getNumAtoms() ==
                    bs->getOrigMol()->getNumAtoms()) {
                  return a.second.first < b.second.first;
                }
                return as->getOrigMol()->getNumAtoms() <
                       bs->getOrigMol()->getNumAtoms();
              });
    size_t stride = d_synthons[i].size() / 40;
    if (!stride) {
      stride = 1;
    }
    for (size_t j = 0; j < d_synthons[i].size(); j += stride) {
      synthonNums[i].push_back(j);
    }
    numSynthons[i] = synthonNums[i].size();
    totSamples *= numSynthons[i];
  }
  details::Stepper stepper(numSynthons);
  std::vector<size_t> theseSynthNums(synthonNums.size(), 0);
  while (stepper.d_currState[0] != numSynthons[0]) {
    if (ControlCHandler::getGotSignal()) {
      return;
    }
    for (size_t i = 0; i < stepper.d_currState.size(); ++i) {
      theseSynthNums[i] = synthonNums[i][stepper.d_currState[i]];
    }
    auto prod = buildProduct(theseSynthNums);
    std::unique_ptr<ExplicitBitVect> prodFP(fpGen.getFingerprint(*prod));
    ExplicitBitVect approxFP(*d_synthons[0][theseSynthNums[0]].second->getFP());
    for (size_t j = 1; j < d_synthons.size(); ++j) {
      approxFP |= *d_synthons[j][theseSynthNums[j]].second->getFP();
    }
    // addFP is what's in the productFP and not in approxFP
    // and subtractFP is vice versa.  The former captures the bits of
    // the molecule formed by the joining the fragments, the latter
    // the bits connecting the dummy atoms.
    std::unique_ptr<ExplicitBitVect> addFP(
        new ExplicitBitVect(*prodFP & ~approxFP));
    IntVect v;
    addFP->getOnBits(v);
    for (auto i : v) {
      naddbitcounts[i]++;
    }
    std::unique_ptr<ExplicitBitVect> subtractFP(
        new ExplicitBitVect(approxFP & ~(*prodFP)));
    subtractFP->getOnBits(v);
    for (auto i : v) {
      nsubbitcounts[i]++;
    }
    stepper.step();
  }

  // This is the fraction of products that must set a bit for
  // it to be included.  Arrived at by empirical means.
  double frac = 0.75;
  d_addFP = std::make_unique<ExplicitBitVect>(numBits);
  for (size_t i = 0; i < naddbitcounts.size(); ++i) {
    if (naddbitcounts[i] > int(totSamples * frac)) {
      d_addFP->setBit(i);
    }
  }
  d_subtractFP = std::make_unique<ExplicitBitVect>(numBits);
  for (size_t i = 0; i < nsubbitcounts.size(); ++i) {
    if (nsubbitcounts[i] > int(totSamples * frac)) {
      d_subtractFP->setBit(i);
    }
  }

  // Take the complement of the subtract FP so it can be used directly
  *d_subtractFP = ~(*d_subtractFP);
}

std::string SynthonSet::buildProductName(
    const std::vector<size_t> &synthNums) const {
  std::vector<std::string> synths(synthNums.size());
  for (size_t i = 0; i < synthNums.size(); ++i) {
    synths[i] = d_synthons[i][synthNums[i]].first;
  }
  return details::buildProductName(d_id, synths);
}

std::unique_ptr<ROMol> SynthonSet::buildProduct(
    const std::vector<size_t> &synthNums) const {
  std::vector<const ROMol *> synths(synthNums.size());
  for (size_t i = 0; i < synthNums.size(); ++i) {
    synths[i] = d_synthons[i][synthNums[i]].second->getOrigMol().get();
  }
  return details::buildProduct(synths);
}

std::unique_ptr<RWMol> SynthonSet::copySynthon(size_t synthonSetNum,
                                               size_t synthonIdx) const {
  std::unique_ptr<RWMol> synthon(new RWMol(
      *d_synthons[synthonSetNum][synthonIdx].second->getOrigMol().get()));
  for (auto &atom : synthon->atoms()) {
    atom->setProp<unsigned int>("molNum", synthonSetNum);
    atom->setProp<unsigned int>("idx", atom->getIdx());
  }
  for (auto &bond : synthon->bonds()) {
    bond->setProp<unsigned int>("molNum", synthonSetNum);
    bond->setProp<unsigned int>("idx", bond->getIdx());
  }

  return synthon;
}

std::vector<std::vector<std::unique_ptr<RWMol>>> SynthonSet::copySynthons()
    const {
  // Synthons are shared, so we need to copy the molecules into a new
  // set that we can fiddle with without upsetting anything else.
  std::vector<std::vector<std::unique_ptr<RWMol>>> synthonMolCopies(
      d_synthons.size());
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    synthonMolCopies[i].reserve(d_synthons[i].size());
    for (size_t j = 0; j < d_synthons[i].size(); ++j) {
      synthonMolCopies[i].emplace_back(copySynthon(i, j));
    }
  }
  return synthonMolCopies;
}

void SynthonSet::initializeSearchOrders() {
  if (d_substructureSearchOrders.empty()) {
    d_substructureSearchOrders = orderSynthonsForSearch(
        [](const Synthon *synth1, const Synthon *synth2) -> bool {
          return synth1->getSearchMol()->getNumAtoms() <
                 synth2->getSearchMol()->getNumAtoms();
        });
  }
  if (d_rascalSearchOrders.empty()) {
    d_rascalSearchOrders = orderSynthonsForSearch(
        [](const Synthon *synth1, const Synthon *synth2) -> bool {
          return synth1->getSearchMol()->getNumAtoms() +
                     synth1->getSearchMol()->getNumBonds() <
                 synth2->getSearchMol()->getNumAtoms() +
                     synth2->getSearchMol()->getNumBonds();
        });
  }
  if (hasFingerprints() && d_fingerprintSearchOrders.empty()) {
    d_fingerprintSearchOrders = orderSynthonsForSearch(
        [](const Synthon *synth1, const Synthon *synth2) -> bool {
          return synth1->getFP()->getNumOnBits() <
                 synth2->getFP()->getNumOnBits();
        });
  }
}

std::vector<std::vector<size_t>> SynthonSet::orderSynthonsForSearch(
    const std::function<bool(const Synthon *synth1, const Synthon *synth2)>
        &cmp) {
  std::vector<std::vector<size_t>> synthonOrders(d_synthons.size());
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    synthonOrders[i].resize(d_synthons[i].size());
    std::iota(synthonOrders[i].begin(), synthonOrders[i].end(), 0);
  }

  for (size_t i = 0; i < d_synthons.size(); ++i) {
    std::vector<std::pair<const Synthon *, size_t>> synthons(
        d_synthons[i].size());
    for (size_t j = 0; j < d_synthons[i].size(); ++j) {
      synthons[j] = std::make_pair(d_synthons[i][j].second, j);
    }
    std::ranges::sort(synthons,
                      [&](const std::pair<const Synthon *, size_t> &a,
                          const std::pair<const Synthon *, size_t> &b) -> bool {
                        return cmp(a.first, b.first);
                      });
    synthonOrders[i].resize(synthons.size());
    for (size_t j = 0; j < d_synthons[i].size(); ++j) {
      synthonOrders[i][j] = synthons[j].second;
    }
  }
  return synthonOrders;
}

std::unique_ptr<SampleMolRec> SynthonSet::makeSampleMolecule(
    Synthon *synthon) const {
  size_t synthonSetNum = 0;
  size_t synthonIdx = 0;
  bool foundIt = false;
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    for (size_t j = 0; j < d_synthons[i].size(); ++j) {
      if (synthon == d_synthons[i][j].second) {
        synthonSetNum = i;
        synthonIdx = j;
        foundIt = true;
        break;
      }
    }
    if (foundIt) {
      break;
    }
  }
  auto sampleMol = std::make_unique<SampleMolRec>();
  if (!foundIt) {
    return sampleMol;
  }

  std::vector<size_t> exSynthons(d_synthons.size(), 0);
  exSynthons[synthonSetNum] = synthonIdx;
  if (!d_substructureSearchOrders.empty()) {
    for (size_t i = 0; i < d_substructureSearchOrders.size(); ++i) {
      if (i != synthonSetNum) {
        exSynthons[i] = d_substructureSearchOrders[i].front();
      }
    }
  }
  sampleMol->d_synthonSet = this;
  sampleMol->d_synthonNums = exSynthons;
  sampleMol->d_synthon = synthon;
  sampleMol->d_synthonSetNum = synthonSetNum;
  for (size_t i = 0; i < d_synthons.size(); ++i) {
    sampleMol->d_numAtoms +=
        d_synthons[i][exSynthons[i]].second->getNumHeavyAtoms();
  }

  return sampleMol;
}

std::unique_ptr<ROMol> SynthonSet::buildMolecule(
    const std::vector<size_t> &synthonNums) const {
  auto combMol = std::make_unique<ROMol>();
  std::string molName = "";
  for (size_t i = 0; i < synthonNums.size(); ++i) {
    auto synthMol = copySynthon(i, synthonNums[i]);
    combMol.reset(combineMols(*combMol, *synthMol));
    std::string sep = i ? ";" : "";
    molName += sep + d_synthons[i][synthonNums[i]].first;
  }

  MolzipParams mzparams;
  mzparams.label = MolzipLabel::Isotope;
  combMol->setProp<std::string>(common_properties::_Name, molName);
  try {
    auto sampleMol = molzip(*combMol, mzparams);
    MolOps::sanitizeMol(*dynamic_cast<RWMol *>(sampleMol.get()));
    return sampleMol;
  } catch (std::exception &e) {
    std::string msg("Error:: in reaction " + getId() +
                    " :: building molecule from synthons :");
    for (size_t i = 0; i < synthonNums.size(); ++i) {
      std::string sep = i ? " and " : " ";
      msg += sep + d_synthons[i][synthonNums[i]].first + " (" +
             d_synthons[i][synthonNums[i]].second->getSmiles() + ")";
    }
    msg += "\n" + std::string(e.what()) + "\n";
    BOOST_LOG(rdErrorLog) << msg;
    throw(e);
  }
}

void SynthonSet::assessRingFormers() {
  d_numRingFormers = 0;
  for (size_t i = 1; i < d_synthons.size(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      if ((d_synthConnPatts[i] & d_synthConnPatts[j]).count() > 1) {
        ++d_numRingFormers;
      }
    }
  }
}

}  // namespace RDKit::SynthonSpaceSearch