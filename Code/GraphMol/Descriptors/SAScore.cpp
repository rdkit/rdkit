//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>

#include "MolDescriptors.h"
#include "SAScore.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <memory>

namespace RDKit {
namespace Descriptors {

namespace {
const char *tableMagic = "RDSASCR1";
constexpr std::size_t magicLength = 8;

//! reads nElements of T from the little-endian table
template <typename T>
void readArray(std::istream &inStream, std::vector<T> &vect,
               std::size_t nElements, const std::string &filename) {
  vect.resize(nElements);
  inStream.read(reinterpret_cast<char *>(vect.data()),
                static_cast<std::streamsize>(nElements * sizeof(T)));
  if (!inStream) {
    throw ValueErrorException("fragment table " + filename +
                              " ended before the data did");
  }
  // a plain if, so that the swap is compiled everywhere and dropped by the
  // optimizer on the little-endian platforms where it is a no-op
  if (HOST_ENDIAN_ORDER != LITTLE_ENDIAN_ORDER) {
    for (auto &value : vect) {
      value = EndianSwapBytes<LITTLE_ENDIAN_ORDER, HOST_ENDIAN_ORDER>(value);
    }
  }
}
}  // namespace

SAScoreFragmentTable::SAScoreFragmentTable(const std::string &filename) {
  std::ifstream inStream(filename, std::ios_base::binary);
  if (!inStream || inStream.bad()) {
    throw BadFileException("could not open fragment table " + filename);
  }

  char magic[magicLength];
  inStream.read(magic, magicLength);
  if (!inStream || std::memcmp(magic, tableMagic, magicLength) != 0) {
    throw ValueErrorException(filename +
                              " is not an SA score fragment table (bad magic)");
  }

  std::uint32_t nScores = 0;
  std::uint32_t nKeys = 0;
  try {
    streamRead(inStream, nScores);
    streamRead(inStream, nKeys);
  } catch (const std::runtime_error &) {
    throw ValueErrorException("truncated header in fragment table " + filename);
  }

  // trust the counts only once the file is actually big enough to hold them,
  // so that a corrupt header cannot ask for an enormous allocation
  const auto headerEnd = inStream.tellg();
  inStream.seekg(0, std::ios_base::end);
  const std::uint64_t payloadBytes =
      static_cast<std::uint64_t>(inStream.tellg() - headerEnd);
  inStream.seekg(headerEnd);
  const std::uint64_t expectedBytes =
      static_cast<std::uint64_t>(nScores) * sizeof(float) +
      static_cast<std::uint64_t>(nKeys) *
          (sizeof(std::uint32_t) + sizeof(std::uint16_t));
  if (payloadBytes != expectedBytes) {
    throw ValueErrorException("fragment table " + filename + " declares " +
                              std::to_string(expectedBytes) +
                              " bytes of data but holds " +
                              std::to_string(payloadBytes));
  }

  readArray(inStream, d_scores, nScores, filename);
  readArray(inStream, d_keys, nKeys, filename);
  readArray(inStream, d_index, nKeys, filename);

  if (!std::is_sorted(d_keys.begin(), d_keys.end())) {
    throw ValueErrorException("fragment table " + filename +
                              " is not sorted by bit id");
  }
  if (std::any_of(d_index.begin(), d_index.end(),
                  [nScores](std::uint16_t i) { return i >= nScores; })) {
    throw ValueErrorException("fragment table " + filename +
                              " has an out-of-range score index");
  }
}

SAScoreFragmentTable::SAScoreFragmentTable(
    const std::vector<std::uint32_t> &bitIds,
    const std::vector<double> &contributions) {
  if (bitIds.size() != contributions.size()) {
    throw ValueErrorException(
        "got " + std::to_string(bitIds.size()) + " bit ids but " +
        std::to_string(contributions.size()) + " contributions");
  }

  // a map both sorts the ids and gives repeated ones last-one-wins semantics
  std::map<std::uint32_t, float> contribByBit;
  for (std::size_t i = 0; i < bitIds.size(); ++i) {
    contribByBit[bitIds[i]] = static_cast<float>(contributions[i]);
  }

  d_scores.reserve(contribByBit.size());
  for (const auto &entry : contribByBit) {
    d_scores.push_back(entry.second);
  }
  std::sort(d_scores.begin(), d_scores.end());
  d_scores.erase(std::unique(d_scores.begin(), d_scores.end()), d_scores.end());
  if (d_scores.size() > std::numeric_limits<std::uint16_t>::max()) {
    throw ValueErrorException("cannot index " +
                              std::to_string(d_scores.size()) +
                              " distinct contributions, the limit is 65535");
  }

  d_keys.reserve(contribByBit.size());
  d_index.reserve(contribByBit.size());
  for (const auto &[bitId, contrib] : contribByBit) {
    d_keys.push_back(bitId);
    d_index.push_back(static_cast<std::uint16_t>(
        std::lower_bound(d_scores.begin(), d_scores.end(), contrib) -
        d_scores.begin()));
  }
}

double SAScoreFragmentTable::score(std::uint32_t bitId) const {
  const auto pos = std::lower_bound(d_keys.begin(), d_keys.end(), bitId);
  if (pos == d_keys.end() || *pos != bitId) {
    return defaultScore;
  }
  return d_scores[d_index[pos - d_keys.begin()]];
}

const SAScoreFragmentTable &defaultSAScoreFragmentTable() {
  static const SAScoreFragmentTable table{[]() {
    const char *rdbase = std::getenv("RDBASE");
    if (!rdbase) {
      throw BadFileException(
          "RDBASE is not set, so the default SA score fragment table cannot be "
          "located; pass a table to calcSAScore() instead");
    }
    return std::string(rdbase) + "/Data/SA_Score/fpscores.bin";
  }()};
  return table;
}

namespace {
//! counts assigned and assignable tetrahedral centers
/*!
  Mirrors Chem.FindMolChiralCenters(mol, includeUnassigned=True), including its
  dependence on the ambient Chirality::useLegacyStereoPerception setting.
*/
unsigned int numChiralCenters(const ROMol &mol) {
  ROMol molCopy(mol);
  // the perception code only ever sets _ChiralityPossible, so stale flags have
  // to be cleared before reassigning
  for (auto atom : molCopy.atoms()) {
    atom->clearProp(common_properties::_ChiralityPossible);
  }
  constexpr bool cleanIt = true;
  constexpr bool force = true;
  constexpr bool flagPossibleStereoCenters = true;
  MolOps::assignStereochemistry(molCopy, cleanIt, force,
                                flagPossibleStereoCenters);

  unsigned int res = 0;
  for (const auto atom : molCopy.atoms()) {
    if (atom->hasProp(common_properties::_CIPCode) ||
        atom->hasProp(common_properties::_ChiralityPossible)) {
      ++res;
    }
  }
  return res;
}

const FingerprintGenerator<std::uint32_t> *sascoreFingerprintGenerator() {
  static const std::unique_ptr<FingerprintGenerator<std::uint32_t>> generator{
      MorganFingerprint::getMorganGenerator<std::uint32_t>(2)};
  return generator.get();
}
}  // namespace

double calcSAScore(const ROMol &mol, const SAScoreFragmentTable &table) {
  const unsigned int nAtoms = mol.getNumAtoms();
  if (!nAtoms) {
    throw ValueErrorException("SA score is not defined for an empty molecule");
  }

  const std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
      sascoreFingerprintGenerator()->getSparseCountFingerprint(mol)};
  const auto &nonzeroElements = fp->getNonzeroElements();

  double fragmentScore = 0.0;
  std::uint64_t nFragments = 0;
  for (const auto &[bitId, count] : nonzeroElements) {
    nFragments += count;
    fragmentScore += table.score(bitId) * count;
  }
  fragmentScore /= nFragments;

  unsigned int nMacrocycles = 0;
  if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
  }
  for (const auto &ring : mol.getRingInfo()->atomRings()) {
    if (ring.size() > 8) {
      ++nMacrocycles;
    }
  }

  const double sizePenalty = std::pow(nAtoms, 1.005) - nAtoms;
  const double stereoPenalty = std::log10(numChiralCenters(mol) + 1);
  const double spiroPenalty = std::log10(calcNumSpiroAtoms(mol) + 1);
  const double bridgePenalty = std::log10(calcNumBridgeheadAtoms(mol) + 1);
  // The paper uses log10(nMacrocycles + 1); this flat penalty gives better
  // results when two or more macrocycles are present.
  const double macrocyclePenalty = nMacrocycles > 0 ? std::log10(2.0) : 0.0;

  const double complexityScore = -sizePenalty - stereoPenalty - spiroPenalty -
                                 bridgePenalty - macrocyclePenalty;

  // not in the original publication; added to make highly symmetrical
  // molecules easier to synthesize
  const std::size_t numBits = nonzeroElements.size();
  const double symmetryScore =
      nAtoms > numBits ? std::log(static_cast<double>(nAtoms) / numBits) * 0.5
                       : 0.0;

  constexpr double rawMin = -4.0;
  constexpr double rawMax = 2.5;
  const double raw = fragmentScore + complexityScore + symmetryScore;
  double sascore = 11.0 - (raw - rawMin + 1.0) / (rawMax - rawMin) * 9.0;

  // sascorer.py subtracts 9 here, which makes log() diverge at the branch
  // boundary and slams scores just above 8 down to 1; see github #8251
  if (sascore > 8.0) {
    sascore = 8.0 + std::log(sascore - 7.0);
  }
  return std::clamp(sascore, 1.0, 10.0);
}

double calcSAScore(const ROMol &mol) {
  return calcSAScore(mol, defaultSAScoreFragmentTable());
}

}  // end of namespace Descriptors
}  // end of namespace RDKit
