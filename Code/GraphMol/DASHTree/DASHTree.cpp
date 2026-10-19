//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "DASHTree.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <GraphMol/ROMol.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>

#include "DASHTreeImpl.h"

namespace RDKit {
namespace DASH {

using namespace detail;

namespace {

//! bytes one value of \p dtype occupies, or 0 if it is not a type we know
std::size_t dtypeWidth(std::uint8_t dtype) {
  switch (dtype) {
    case DTYPE_FLOAT16:
      return 2;
    case DTYPE_FLOAT32:
    case DTYPE_INT32:
      return 4;
    case DTYPE_FLOAT64:
      return 8;
    default:
      return 0;
  }
}

//! reads a NUL-padded fixed-width name out of the property directory
std::string propertyEntryName(const char *entry) {
  // a name of the full width carries no terminator, so never look past it
  return std::string(entry,
                     std::find(entry, entry + dashPropertyNameLen, '\0'));
}

//! Bound on the embeddings one descent keeps alive. Atoms that look alike
//! multiply them: a dozen at most for drug-like molecules, about a hundred for
//! tetraphenylmethane. Embeddings past the bound are dropped, the one place
//! where the atom order could still reach the result.
constexpr std::size_t maxEmbeddings = 256;

}  // namespace

// ---------------------------------------------------------------------------
//  DASHTree::Impl -- mapping the container
// ---------------------------------------------------------------------------

DASHTree::Impl::Impl(const std::string &filename,
                     const std::vector<std::string> &properties)
    : d_filename(filename) {
  try {
    d_map.open(filename);
  } catch (const std::exception &e) {
    throw BadFileException("could not open DASH tree '" + filename +
                           "': " + e.what());
  }
  if (!d_map.is_open()) {
    throw BadFileException("could not open DASH tree '" + filename + "'");
  }
  const std::size_t size = d_map.size();
  if (size < dashHeaderSize) {
    throw BadFileException("'" + filename + "' is too short to be a DASH tree");
  }
  const char *base = d_map.data();
  std::memcpy(&d_header, base, sizeof(d_header));

  if (std::memcmp(d_header.magic, dashMagic, dashMagicLen) != 0) {
    throw BadFileException("'" + filename + "' is not a DASH tree container");
  }
  if (d_header.endianId != dashEndianId) {
    throw ValueErrorException(
        "endianness mismatch in DASH tree '" + filename +
        "': the file was written on a machine of the opposite byte order");
  }
  if (d_header.version != dashFormatVersion) {
    throw ValueErrorException(
        "'" + filename + "' is DASH container version " +
        std::to_string(d_header.version) + "; this build reads version " +
        std::to_string(dashFormatVersion) +
        (d_header.version > dashFormatVersion
             ? ". Please rebuild with a more recent version of the RDKit."
             : ". Please regenerate it with python -m rdkit.Chem.DASHTree.convert."));
  }
  if (d_header.fileSize != size) {
    throw ValueErrorException("'" + filename +
                              "' is truncated: the header says " +
                              std::to_string(d_header.fileSize) +
                              " bytes, the file is " + std::to_string(size));
  }
  if (!d_header.numBranches || !d_header.numNodes) {
    throw ValueErrorException("'" + filename + "' contains no tree");
  }

  // Every array must lie inside the file and be aligned for its type: the
  // arrays are read in place, so a bad offset would be an out-of-bounds or
  // misaligned load rather than a parse error.
  auto resolve = [&](std::uint64_t offset, std::uint64_t nBytes,
                     std::size_t alignment, const char *what) -> const char * {
    if (offset > size || nBytes > size - offset) {
      throw ValueErrorException("'" + filename + "' is corrupt: " + what +
                                " runs past the end of the file");
    }
    if (offset % alignment != 0) {
      throw ValueErrorException("'" + filename + "' is corrupt: " + what +
                                " is not " + std::to_string(alignment) +
                                "-byte aligned");
    }
    return base + offset;
  };

  const std::uint64_t nNodes = d_header.numNodes;
  d_branchRoot = reinterpret_cast<const std::uint32_t *>(
      resolve(d_header.offBranchRoot, d_header.numBranches * 4ull, 4,
              "the branch table"));
  d_node = reinterpret_cast<const NodeRecord *>(
      resolve(d_header.offNodeRecord, nNodes * sizeof(NodeRecord),
              alignof(NodeRecord), "the node records"));
  d_nodeAttn = reinterpret_cast<const float *>(
      resolve(d_header.offNodeAttn, nNodes * 4ull, 4, "the attentions"));
  resolve(d_header.offPropertyDir,
          d_header.numProperties * (std::uint64_t)dashPropertyEntrySize, 8,
          "the property directory");

  // the roots must ascend from 0: a branch owns every id up to the next root,
  // which is what lets a node be addressed by branch and offset
  for (std::uint32_t b = 0; b < d_header.numBranches; ++b) {
    if (d_branchRoot[b] >= nNodes ||
        (b ? d_branchRoot[b] <= d_branchRoot[b - 1] : d_branchRoot[b] != 0)) {
      throw ValueErrorException(
          "'" + filename +
          "' is corrupt: branch roots out of range or out of order");
    }
  }
  // the atom-feature table: one class per branch
  if (d_header.numFeatures != d_header.numBranches) {
    throw ValueErrorException(
        "'" + filename + "' is corrupt: " +
        std::to_string(d_header.numFeatures) + " atom-feature classes for " +
        std::to_string(d_header.numBranches) + " branches");
  }
  const auto *featureTable = reinterpret_cast<const AtomFeature *>(resolve(
      d_header.offFeatureTable, d_header.numFeatures * sizeof(AtomFeature),
      alignof(AtomFeature), "the atom-feature table"));
  try {
    d_features = AtomFeatureLookup(featureTable, d_header.numFeatures);
  } catch (const ValueErrorException &e) {
    throw ValueErrorException("'" + filename + "': " + e.what());
  }

  // resolve the requested property columns
  const bool wantAll = properties.empty();
  for (std::uint32_t i = 0; i < d_header.numProperties; ++i) {
    const char *entry =
        base + d_header.offPropertyDir + i * dashPropertyEntrySize;
    const std::string name = propertyEntryName(entry);
    if (!wantAll && std::find(properties.begin(), properties.end(), name) ==
                        properties.end()) {
      continue;
    }
    std::uint8_t dtype;
    std::uint64_t offset;
    std::memcpy(&dtype, entry + dashPropertyNameLen, 1);
    std::memcpy(&offset, entry + dashPropertyNameLen + 8, 8);

    const std::size_t width = dtypeWidth(dtype);
    if (!width) {
      throw ValueErrorException("'" + filename + "': property '" + name +
                                "' has unknown storage type " +
                                std::to_string(dtype));
    }
    PropertyColumn column;
    column.name = name;
    column.dtype = static_cast<PropertyDType>(dtype);
    column.data = resolve(offset, nNodes * width, width,
                          ("property '" + name + "'").c_str());
    d_properties.push_back(column);
  }
  for (const auto &want : properties) {
    if (propertyIndex(want) < 0) {
      std::string have;
      for (const auto &name : allPropertyNames()) {
        have += (have.empty() ? "" : ", ") + name;
      }
      throw ValueErrorException("'" + filename + "' has no property '" + want +
                                "'; it has: " + have);
    }
  }
}

int DASHTree::Impl::propertyIndex(const std::string &name) const {
  for (std::size_t i = 0; i < d_properties.size(); ++i) {
    if (d_properties[i].name == name) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

const PropertyColumn &DASHTree::Impl::property(const std::string &name) const {
  const int idx = propertyIndex(name);
  if (idx < 0) {
    std::string have;
    for (const auto &column : d_properties) {
      have += (have.empty() ? "" : ", ") + column.name;
    }
    throw ValueErrorException("property '" + name +
                              "' was not loaded from this DASH tree; loaded: " +
                              (have.empty() ? "(none)" : have));
  }
  return d_properties[idx];
}

std::vector<std::string> DASHTree::Impl::allPropertyNames() const {
  std::vector<std::string> names;
  const char *base = d_map.data();
  for (std::uint32_t i = 0; i < d_header.numProperties; ++i) {
    names.push_back(propertyEntryName(base + d_header.offPropertyDir +
                                      i * dashPropertyEntrySize));
  }
  return names;
}

void DASHTree::Impl::throwCorruptChildren(std::uint32_t node) const {
  throw ValueErrorException(
      "'" + d_filename + "' is corrupt: the children of node " +
      std::to_string(node) + " lie outside its branch, or not after it");
}

// ---------------------------------------------------------------------------
//  MolMatcher -- the descent
// ---------------------------------------------------------------------------

void MolMatcher::setMolecule(const ROMol &mol) {
  // The feature tuple needs the implicit-valence cache. RDKit would throw from
  // deep inside getTotalNumHs; say so here instead, and say it once rather than
  // from N threads.
  if (mol.needsUpdatePropertyCache()) {
    throw ValueErrorException(
        "DASH needs the molecule's property cache: call "
        "updatePropertyCache() (or sanitize the molecule) first");
  }
  // only sanitization sets the conjugation flags the atom classes read, and it
  // always initialises the ring info, so that stands in for "sanitized"
  if (!mol.getRingInfo()->isInitialized()) {
    BOOST_LOG(rdWarningLog) << "DASH: the molecule is not sanitized, so its "
                               "atoms may be misclassified"
                            << std::endl;
  }

  d_numAtoms = mol.getNumAtoms();
  d_feature.resize(d_numAtoms);
  d_nbrStart.assign(d_numAtoms + 1, 0);

  for (unsigned int i = 0; i < d_numAtoms; ++i) {
    const Atom *atom = mol.getAtomWithIdx(i);
    const int feature = d_tree.d_features.index(atom);
    if (feature < 0) {
      throw ValueErrorException(
          "DASH cannot assign atom " + std::to_string(i) + " (" +
          atom->getSymbol() + ", degree " + std::to_string(atom->getDegree()) +
          ", charge " + std::to_string(atom->getFormalCharge()) + ", " +
          std::to_string(atom->getTotalNumHs(true)) +
          " H): its atom type is not one of the " +
          std::to_string(d_tree.d_features.size()) + " this tree covers");
    }
    d_feature[i] = feature;
    d_nbrStart[i + 1] = atom->getDegree();
  }
  for (unsigned int i = 0; i < d_numAtoms; ++i) {
    d_nbrStart[i + 1] += d_nbrStart[i];
  }
  d_nbrAtom.resize(d_nbrStart[d_numAtoms]);
  d_nbrKey.resize(d_nbrStart[d_numAtoms]);

  // The neighbour list is walked in RDKit's own order, which only the legacy
  // tie-breaking depends on: there, of two candidates carrying the same key,
  // the earlier one is bound.
  for (unsigned int i = 0; i < d_numAtoms; ++i) {
    const Atom *atom = mol.getAtomWithIdx(i);
    std::uint32_t at = d_nbrStart[i];
    for (const auto bond : mol.atomBonds(atom)) {
      const unsigned int nbrIdx = bond->getOtherAtomIdx(i);
      d_nbrAtom[at] = nbrIdx;
      // the key of the *neighbour*, attached to atom i by the bond between
      // them; DASH's bond descriptor is 4 if conjugated, else the bond order
      d_nbrKey[at] = packPartialMatchKey(
          d_feature[nbrIdx],
          bond->getIsConjugated() ? 4 : static_cast<int>(bond->getBondType()));
      ++at;
    }
  }
}

// The next node is the first child, in tree order, that any live embedding can
// extend to. The tree's builder created a node's children in falling order of
// attention, so that is the best-attended continuation. Every candidate
// carrying the child's key is bound, each in its own copy of the embedding:
// atoms that look alike so far -- the two CH2 next to an isoindoline NH, the
// ortho carbons of a ring -- stay in play until the tree tells them apart. So
// the result does not depend on the atom order, and symmetry-equivalent atoms
// get the same value. With legacyTieBreaking only the first candidate is
// bound, as the DASH-tree python package does.
//
// The nesting matters: children outside, candidates inside. Iterating the
// candidates on the outside -- the natural shape if you were tempted to sort
// the children by key and binary-search -- picks a different node for the
// great majority of atoms.
bool MolMatcher::step(bool legacyTieBreaking) {
  const auto [first, last] = d_tree.children(d_path.back(), d_branchEnd);
  std::uint64_t filter = 0;
  for (std::size_t i = 0; i < d_numLive; ++i) {
    filter |= d_live[i].candFilter;
  }
  for (std::uint32_t child = first; child < last; ++child) {
    const std::uint16_t key = d_tree.d_node[child].key;
    // one test throws out most children before the candidate scan
    if (!((filter >> (key & 63u)) & 1u)) {
      continue;
    }
    findMatches(key, legacyTieBreaking);
    if (!d_matches.empty()) {
      d_path.push_back(child);
      advance();
      return true;
    }
  }
  return false;
}

void MolMatcher::findMatches(std::uint16_t key, bool legacyTieBreaking) {
  d_matches.clear();
  for (std::uint32_t i = 0; i < d_numLive; ++i) {
    const Embedding &embedding = d_live[i];
    bool terminalFound = false;
    for (std::uint32_t c = 0; c < embedding.candKey.size(); ++c) {
      if (embedding.candKey[c] != key) {
        continue;
      }
      // Degree-one candidates with the same key hang off the same atom by the
      // same bond, so they are interchangeable and one stands for all.
      const std::uint32_t atom = embedding.candAtom[c];
      const bool terminal = d_nbrStart[atom + 1] - d_nbrStart[atom] == 1;
      if (terminal && terminalFound) {
        continue;
      }
      terminalFound |= terminal;
      d_matches.emplace_back(i, c);
      if (legacyTieBreaking || d_matches.size() == maxEmbeddings) {
        return;
      }
    }
  }
}

void MolMatcher::advance() {
  const std::size_t numMatches = d_matches.size();
  if (numMatches == 1) {
    // the usual case: one embedding carries on, extended in place
    const std::uint32_t from = d_matches[0].first;
    if (from) {
      std::swap(d_live[0], d_live[from]);
    }
    extend(d_live[0], d_matches[0].second);
    d_numLive = 1;
    return;
  }
  if (d_next.size() < numMatches) {
    d_next.resize(numMatches);
  }
  for (std::size_t m = 0; m < numMatches; ++m) {
    const std::uint32_t from = d_matches[m].first;
    // an embedding's last match takes it over; only earlier ones copy it
    if (m + 1 < numMatches && d_matches[m + 1].first == from) {
      d_next[m] = d_live[from];
    } else {
      std::swap(d_next[m], d_live[from]);
    }
    extend(d_next[m], d_matches[m].second);
  }
  std::swap(d_live, d_next);
  d_numLive = numMatches;
}

void MolMatcher::extend(Embedding &embedding, std::uint32_t cand) const {
  embedding.atoms.push_back(embedding.candAtom[cand]);
  // The consumed candidate stays seen, but its key is retired so it cannot be
  // bound twice.
  embedding.candKey[cand] = noMatchKey;
  addCandidates(embedding,
                static_cast<std::uint32_t>(embedding.atoms.size() - 1));
}

std::uint32_t MolMatcher::match(unsigned int atomIdx,
                                const DASHParams &params) {
  PRECONDITION(atomIdx < d_numAtoms, "atom index out of range");
  const std::uint32_t branch = static_cast<std::uint32_t>(d_feature[atomIdx]);
  const std::uint32_t root = d_tree.d_branchRoot[branch];
  d_branchEnd = d_tree.branchEnd(branch);

  d_path.clear();
  d_path.push_back(root);
  d_numLive = 1;
  Embedding &start = d_live[0];
  start.atoms.clear();
  start.candKey.clear();
  start.candAtom.clear();
  start.seen.assign((d_numAtoms + 63) / 64, 0);
  start.candFilter = 0;

  unsigned int maxDepth = params.maxDepth;

  if (d_tree.d_features.feature(branch).atomicNum == 1) {
    // Hydrogens are only ever matched implicitly, through the numHs component
    // of their heavy neighbour's feature: the descent starts at the heavy atom
    // and spends one level getting there, without accumulating its attention.
    // The sole hydrogen feature class has degree one, so an H that got past
    // setMolecule has exactly one neighbour and it is not itself an H.
    const std::uint32_t nbrLo = d_nbrStart[atomIdx];
    CHECK_INVARIANT(d_nbrStart[atomIdx + 1] - nbrLo == 1,
                    "hydrogen of unexpected degree");
    const std::uint32_t heavy = d_nbrAtom[nbrLo];
    const std::uint16_t heavyKey = packMatchKey(d_feature[heavy], -1, -1);
    start.see(heavy);
    start.candKey.push_back(heavyKey);
    start.candAtom.push_back(heavy);
    start.candFilter |= 1ull << (heavyKey & 63u);
    if (!step(params.legacyTieBreaking)) {
      // no node in the tree describes this hydrogen's environment
      d_path.clear();
      return branch;
    }
    maxDepth = std::max(maxDepth, 1u) - 1;
  } else {
    start.see(atomIdx);
    start.atoms.push_back(atomIdx);
    addCandidates(start, 0);
  }

  if (maxDepth <= 1) {
    return branch;
  }

  // The stop flag holds the threshold test precomputed for the header's
  // default threshold, so with that threshold and no increment cutoff the
  // descent never reads, or pages in, the attention column.
  const bool useStopFlag =
      params.attentionThreshold == d_tree.d_header.defaultAttentionThreshold &&
      params.attentionIncrementThreshold <= 0.0;

  double cumulativeAttention = 0.0;
  for (unsigned int depth = 1; depth < maxDepth; ++depth) {
    if (!step(params.legacyTieBreaking)) {
      break;
    }
    const std::uint32_t childNode = d_path.back();

    // Both tests are strict and both come after the step has been taken, so
    // the node that crosses the threshold stays on the path and is the first
    // one the property walk reads.
    if (useStopFlag) {
      if (d_tree.d_node[childNode].flags & nodeFlagStop) {
        break;
      }
    } else {
      const double attention = d_tree.d_nodeAttn[childNode];
      cumulativeAttention += attention;
      if (cumulativeAttention > params.attentionThreshold) {
        break;
      }
      if (attention < params.attentionIncrementThreshold) {
        break;
      }
    }
  }
  return branch;
}

// Candidates are maintained incrementally rather than rebuilt at every level.
// That is order-identical to a rebuild: an atom's position in the matched
// subgraph never changes once it joins, candidates are emitted grouped by that
// position, and the newest atom always has the highest one -- so its
// neighbours belong at the end of the list either way.
void MolMatcher::addCandidates(Embedding &embedding,
                               std::uint32_t position) const {
  const std::uint32_t atomIdx = embedding.atoms[position];
  const std::uint32_t lo = d_nbrStart[atomIdx];
  const std::uint32_t hi = d_nbrStart[atomIdx + 1];
  for (std::uint32_t n = lo; n < hi; ++n) {
    const std::uint32_t nbr = d_nbrAtom[n];
    if (embedding.see(nbr)) {
      continue;
    }
    const std::uint16_t key =
        withConAtom(d_nbrKey[n], static_cast<int>(position));
    embedding.candKey.push_back(key);
    embedding.candAtom.push_back(nbr);
    embedding.candFilter |= 1ull << (key & 63u);
  }
}

double MolMatcher::pathValue(const PropertyColumn &column) const {
  for (std::size_t i = d_path.size(); i-- > 0;) {
    const double value = column.value(d_path[i]);
    if (!std::isnan(value)) {
      return value;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------------------
//  DASHTree -- the public surface
// ---------------------------------------------------------------------------

DASHTree::DASHTree(const std::string &filename,
                   const std::vector<std::string> &properties)
    : d_impl(std::make_unique<Impl>(filename, properties)) {}

DASHTree::~DASHTree() = default;
DASHTree::DASHTree(DASHTree &&) noexcept = default;
DASHTree &DASHTree::operator=(DASHTree &&) noexcept = default;

unsigned int DASHTree::getNumBranches() const {
  return d_impl->d_header.numBranches;
}

const AtomFeature &DASHTree::getAtomFeature(unsigned int index) const {
  if (index >= d_impl->d_features.size()) {
    throw IndexErrorException(static_cast<int>(index));
  }
  return d_impl->d_features.feature(index);
}

int DASHTree::getAtomFeatureIndex(const Atom *atom) const {
  return d_impl->d_features.index(atom);
}

unsigned int DASHTree::getNumNodes() const { return d_impl->d_header.numNodes; }

const std::string &DASHTree::getFilename() const { return d_impl->d_filename; }

std::vector<std::string> DASHTree::getPropertyNames() const {
  std::vector<std::string> names;
  names.reserve(d_impl->d_properties.size());
  for (const auto &column : d_impl->d_properties) {
    names.push_back(column.name);
  }
  return names;
}

std::vector<std::string> DASHTree::getAvailablePropertyNames() const {
  return d_impl->allPropertyNames();
}

void DASHTree::getAtomNodePath(const ROMol &mol, unsigned int atomIdx,
                               std::vector<std::uint32_t> &nodePath,
                               const DASHParams &params) const {
  if (atomIdx >= mol.getNumAtoms()) {
    throw IndexErrorException(static_cast<int>(atomIdx));
  }
  MolMatcher matcher(*d_impl);
  matcher.setMolecule(mol);
  const std::uint32_t branch = matcher.match(atomIdx, params);

  // Led by the branch index, so this reads like the node paths the python
  // implementation reports; the node ids are counted from the branch root.
  nodePath.clear();
  nodePath.push_back(branch);
  const std::uint32_t root = d_impl->d_branchRoot[branch];
  for (const auto node : matcher.path()) {
    nodePath.push_back(node - root);
  }
}

double DASHTree::getAtomProperty(const ROMol &mol, unsigned int atomIdx,
                                 const std::string &property,
                                 const DASHParams &params) const {
  if (atomIdx >= mol.getNumAtoms()) {
    throw IndexErrorException(static_cast<int>(atomIdx));
  }
  const PropertyColumn &column = d_impl->property(property);
  MolMatcher matcher(*d_impl);
  matcher.setMolecule(mol);
  matcher.match(atomIdx, params);
  return matcher.pathValue(column);
}

void DASHTree::getMolProperty(const ROMol &mol, const std::string &property,
                              std::vector<double> &res,
                              const DASHParams &params) const {
  const PropertyColumn &column = d_impl->property(property);
  MolMatcher matcher(*d_impl);
  matcher.setMolecule(mol);
  res.resize(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    matcher.match(i, params);
    res[i] = matcher.pathValue(column);
  }
}

// ---------------------------------------------------------------------------
//  the tree itself
// ---------------------------------------------------------------------------

DASHTreeNode DASHTree::getRoot(unsigned int branch) const {
  return getNode(branch, 0);
}

DASHTreeNode DASHTree::getNode(unsigned int branch,
                               std::uint32_t nodeId) const {
  const std::uint32_t nBranches = d_impl->d_header.numBranches;
  if (branch >= nBranches) {
    throw IndexErrorException(static_cast<int>(branch));
  }
  // branches tile the node array in order, so one owns every id up to the
  // next root
  const std::uint32_t root = d_impl->d_branchRoot[branch];
  const std::uint32_t end = d_impl->branchEnd(branch);
  if (nodeId >= end - root) {
    throw IndexErrorException(static_cast<int>(nodeId));
  }
  return DASHTreeNode(*d_impl, root + nodeId);
}

void DASHTree::getMatchedSubstructure(const ROMol &mol, unsigned int atomIdx,
                                      std::vector<unsigned int> &atoms,
                                      const DASHParams &params) const {
  if (atomIdx >= mol.getNumAtoms()) {
    throw IndexErrorException(static_cast<int>(atomIdx));
  }
  MolMatcher matcher(*d_impl);
  matcher.setMolecule(mol);
  matcher.match(atomIdx, params);
  const auto &subgraph = matcher.subgraph();
  atoms.assign(subgraph.begin(), subgraph.end());
}

// ---------------------------------------------------------------------------
//  DASHTreeNode
// ---------------------------------------------------------------------------

DASHTreeNode::DASHTreeNode(const DASHTree::Impl &impl, std::uint32_t absoluteId)
    : dp_impl(&impl), d_abs(absoluteId) {}

unsigned int DASHTreeNode::getBranch() const {
  const std::uint32_t *begin = dp_impl->d_branchRoot;
  const std::uint32_t *end = begin + dp_impl->d_header.numBranches;
  // the roots ascend, so the branch is the last root at or below this id
  return static_cast<unsigned int>(std::upper_bound(begin, end, d_abs) -
                                   begin) -
         1;
}

std::uint32_t DASHTreeNode::getId() const {
  return d_abs - dp_impl->d_branchRoot[getBranch()];
}

std::uint16_t DASHTreeNode::getKey() const {
  return dp_impl->d_node[d_abs].key;
}

int DASHTreeNode::getAtomFeatureIndex() const { return keyAtomType(getKey()); }

int DASHTreeNode::getConAtom() const { return keyConAtom(getKey()); }

int DASHTreeNode::getConType() const { return keyConType(getKey()); }

unsigned int DASHTreeNode::getNumChildren() const {
  return dp_impl->d_node[d_abs].numChildren;
}

DASHTreeNode DASHTreeNode::getChild(unsigned int i) const {
  const unsigned int branch = getBranch();
  const auto [first, last] =
      dp_impl->children(d_abs, dp_impl->branchEnd(branch));
  if (i >= last - first) {
    throw IndexErrorException(static_cast<int>(i));
  }
  return DASHTreeNode(*dp_impl, first + i);
}

float DASHTreeNode::getAttention() const { return dp_impl->d_nodeAttn[d_abs]; }

bool DASHTreeNode::isStopNode() const {
  return dp_impl->d_node[d_abs].flags & nodeFlagStop;
}

double DASHTreeNode::getValue(const std::string &property) const {
  return dp_impl->property(property).value(d_abs);
}

}  // namespace DASH
}  // namespace RDKit
