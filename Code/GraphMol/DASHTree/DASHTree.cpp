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

#include "DASHTreeImpl.h"

namespace RDKit {
namespace DASH {

using namespace detail;

namespace {

//! reads through a mapped range so later accesses do not fault
void touchRange(const void *data, std::size_t nBytes) {
  const auto *p = static_cast<const volatile unsigned char *>(data);
  unsigned char sink = 0;
  for (std::size_t i = 0; i < nBytes; i += 4096) {
    sink ^= p[i];
  }
  if (nBytes) {
    sink ^= p[nBytes - 1];
  }
  (void)sink;
}

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
  char buf[dashPropertyNameLen + 1];
  std::memcpy(buf, entry, dashPropertyNameLen);
  buf[dashPropertyNameLen] = '\0';
  return std::string(buf);
}

}  // namespace

// ---------------------------------------------------------------------------
//  DASHTree::Impl -- mapping the container
// ---------------------------------------------------------------------------

DASHTree::Impl::Impl(const std::string &filename,
                     const std::vector<std::string> &properties, bool prefetch)
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
    throw ValueErrorException("'" + filename +
                              "' is too short to be a DASH tree");
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
             : ". Please re-run tools/dash_convert.py to regenerate it."));
  }
  if (d_header.fileSize != size) {
    throw ValueErrorException("'" + filename + "' is truncated: the header says " +
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
  if (d_header.offNodeSourceId) {
    d_nodeSourceId = reinterpret_cast<const std::uint32_t *>(
        resolve(d_header.offNodeSourceId, nNodes * 4ull, 4, "the source ids"));
  }
  resolve(d_header.offPropertyDir,
          d_header.numProperties * (std::uint64_t)dashPropertyEntrySize, 8,
          "the property directory");

  if (d_header.maxChildren > 255) {
    throw ValueErrorException("'" + filename +
                              "' is corrupt: it claims a node with more than "
                              "255 children");
  }
  // A tree whose keys the reader cannot express would silently never match, so
  // refuse it rather than return wrong numbers.
  if (d_header.maxAtomType > maxKeyAtomType ||
      d_header.maxConAtom > maxKeyConAtom ||
      d_header.maxConType > maxKeyConType) {
    throw ValueErrorException(
        "'" + filename +
        "' uses match keys outside the range this build can encode");
  }
  // the roots must ascend: a branch owns every id up to the next root, which
  // is what lets a node be addressed by branch and offset
  for (std::uint32_t b = 0; b < d_header.numBranches; ++b) {
    if (d_branchRoot[b] >= nNodes ||
        (b && d_branchRoot[b] <= d_branchRoot[b - 1])) {
      throw ValueErrorException(
          "'" + filename +
          "' is corrupt: branch roots out of range or out of order");
    }
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

  if (prefetch) {
    touchRange(d_node, nNodes * sizeof(NodeRecord));
    touchRange(d_nodeAttn, nNodes * 4ull);
    for (const auto &column : d_properties) {
      touchRange(column.data, nNodes * dtypeWidth(column.dtype));
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

  d_numAtoms = mol.getNumAtoms();
  d_feature.resize(d_numAtoms);
  d_nbrStart.assign(d_numAtoms + 1, 0);
  d_seen.assign(d_numAtoms, 0);
  d_generation = 0;

  for (unsigned int i = 0; i < d_numAtoms; ++i) {
    const Atom *atom = mol.getAtomWithIdx(i);
    const int feature = atomFeatureIndex(atom);
    if (feature < 0) {
      throw ValueErrorException(
          "DASH cannot assign atom " + std::to_string(i) + " (" +
          atom->getSymbol() + ", degree " + std::to_string(atom->getDegree()) +
          ", charge " + std::to_string(atom->getFormalCharge()) + ", " +
          std::to_string(atom->getTotalNumHs(true)) +
          " H): its atom type is not one of the " +
          std::to_string(numAtomFeatures) + " the DASH trees cover");
    }
    d_feature[i] = feature;
    d_nbrStart[i + 1] = atom->getDegree();
  }
  for (unsigned int i = 0; i < d_numAtoms; ++i) {
    d_nbrStart[i + 1] += d_nbrStart[i];
  }
  d_nbrAtom.resize(d_nbrStart[d_numAtoms]);
  d_nbrKey.resize(d_nbrStart[d_numAtoms]);

  // The neighbour list is walked in RDKit's own order, which the match depends
  // on: where two candidates carry the same key, the earlier one is bound.
  for (unsigned int i = 0; i < d_numAtoms; ++i) {
    const Atom *atom = mol.getAtomWithIdx(i);
    std::uint32_t at = d_nbrStart[i];
    for (const auto nbr : mol.atomNeighbors(atom)) {
      const unsigned int nbrIdx = nbr->getIdx();
      d_nbrAtom[at] = nbrIdx;
      // the key of the *neighbour*, attached to atom i by the bond between them
      d_nbrKey[at] = packPartialMatchKey(
          d_feature[nbrIdx], dashBondType(mol.getBondBetweenAtoms(nbrIdx, i)));
      ++at;
    }
  }
}

// The nesting here IS the algorithm, and getting it the wrong way round is not
// a small error: the first child *in tree order* that matches any live
// candidate wins, and among the candidates carrying that child's key the
// earliest one is bound. Iterating candidates on the outside -- the natural
// shape if you were tempted to sort the children by key and binary-search --
// picks a different node for the great majority of atoms.
bool MolMatcher::pickChild(std::uint32_t node, std::uint32_t &childNode,
                           std::uint32_t &childAtom) {
  const NodeRecord &record = d_tree.d_node[node];
  const std::uint32_t first = record.firstChild;
  const std::uint32_t last = first + record.numChildren;
  const std::size_t nCand = d_candKey.size();
  for (std::uint32_t child = first; child < last; ++child) {
    const std::uint16_t key = d_tree.d_node[child].key;
    // one test throws out most children before the candidate scan
    if (!((d_candFilter >> (key & 63u)) & 1u)) {
      continue;
    }
    for (std::size_t c = 0; c < nCand; ++c) {
      if (d_candKey[c] == key) {
        childNode = child;
        childAtom = d_candAtom[c];
        return true;
      }
    }
  }
  return false;
}

std::uint32_t MolMatcher::match(unsigned int atomIdx,
                                const DASHParams &params) {
  PRECONDITION(atomIdx < d_numAtoms, "atom index out of range");
  const std::uint32_t branch = static_cast<std::uint32_t>(d_feature[atomIdx]);
  const std::uint32_t root = d_tree.d_branchRoot[branch];

  d_path.clear();
  d_path.push_back(root);
  d_subgraph.clear();
  d_candKey.clear();
  d_candAtom.clear();
  d_candFilter = 0;
  ++d_generation;

  // signed, so the hydrogen decrement below behaves the way the reference
  // implementation's does when a caller passes 0 or 1
  int maxDepth = static_cast<int>(params.maxDepth);

  if (getAtomFeatureTable()[branch].atomicNum == 1) {
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
    d_candKey.push_back(heavyKey);
    d_candAtom.push_back(heavy);
    d_candFilter |= 1ull << (heavyKey & 63u);

    std::uint32_t childNode, childAtom;
    if (!pickChild(root, childNode, childAtom)) {
      // no node in the tree describes this hydrogen's environment
      d_path.clear();
      return branch;
    }
    d_path.push_back(childNode);
    d_subgraph.push_back(heavy);
    d_seen[heavy] = d_generation;
    --maxDepth;

    d_candKey.clear();
    d_candAtom.clear();
    d_candFilter = 0;
  } else {
    d_subgraph.push_back(atomIdx);
    d_seen[atomIdx] = d_generation;
  }

  if (maxDepth <= 1) {
    return branch;
  }

  // Candidates are maintained incrementally rather than rebuilt at every level.
  // That is order-identical to a rebuild: an atom's position in the matched
  // subgraph never changes once it joins, candidates are emitted grouped by
  // that position, and the newest atom always has the highest one -- so its
  // neighbours belong at the end of the list either way.
  addCandidates(0);

  // The stop flag in each node record is the answer to "is the attention
  // accumulated down to here past the threshold?", precomputed by the
  // converter for the threshold recorded in the header. When the caller wants
  // that threshold and no increment cutoff -- the default, and what every
  // published use does -- the descent never reads the attentions at all, and
  // their 37.8 MB never enters memory. Any other setting takes the slow path.
  const bool useStopFlag =
      params.attentionThreshold == d_tree.d_header.defaultAttentionThreshold &&
      params.attentionIncrementThreshold <= 0.0;

  double cumulativeAttention = 0.0;
  for (int depth = 1; depth < maxDepth; ++depth) {
    std::uint32_t childNode, childAtom;
    if (!pickChild(d_path.back(), childNode, childAtom)) {
      break;
    }
    d_path.push_back(childNode);
    const std::uint32_t position =
        static_cast<std::uint32_t>(d_subgraph.size());
    d_subgraph.push_back(childAtom);
    // The consumed candidate keeps its "seen" stamp -- it is in the subgraph
    // now, so it stays ineligible -- but its key is retired so it cannot be
    // bound twice.
    retireCandidate(childAtom);
    addCandidates(position);

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

void MolMatcher::addCandidates(std::uint32_t position) {
  const std::uint32_t atomIdx = d_subgraph[position];
  const std::uint32_t lo = d_nbrStart[atomIdx];
  const std::uint32_t hi = d_nbrStart[atomIdx + 1];
  for (std::uint32_t n = lo; n < hi; ++n) {
    const std::uint32_t nbr = d_nbrAtom[n];
    if (d_seen[nbr] == d_generation) {
      continue;
    }
    d_seen[nbr] = d_generation;
    const std::uint16_t key =
        withConAtom(d_nbrKey[n], static_cast<int>(position));
    d_candKey.push_back(key);
    d_candAtom.push_back(nbr);
    d_candFilter |= 1ull << (key & 63u);
  }
}

void MolMatcher::retireCandidate(std::uint32_t atomIdx) {
  for (std::size_t c = 0; c < d_candAtom.size(); ++c) {
    if (d_candAtom[c] == atomIdx && d_candKey[c] != noMatchKey) {
      d_candKey[c] = noMatchKey;
      return;
    }
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
                   const std::vector<std::string> &properties, bool prefetch)
    : d_impl(new Impl(filename, properties, prefetch)) {}

DASHTree::~DASHTree() = default;

unsigned int DASHTree::numBranches() const {
  return d_impl->d_header.numBranches;
}

std::uint64_t DASHTree::numNodes() const { return d_impl->d_header.numNodes; }

std::uint64_t DASHTree::mappedSize() const { return d_impl->d_map.size(); }

const std::string &DASHTree::filename() const { return d_impl->d_filename; }

bool DASHTree::hasSourceNodeIds() const {
  return d_impl->d_nodeSourceId != nullptr;
}

std::vector<std::string> DASHTree::propertyNames() const {
  std::vector<std::string> names;
  names.reserve(d_impl->d_properties.size());
  for (const auto &column : d_impl->d_properties) {
    names.push_back(column.name);
  }
  return names;
}

std::vector<std::string> DASHTree::availablePropertyNames() const {
  return d_impl->allPropertyNames();
}

bool DASHTree::hasProperty(const std::string &name) const {
  return d_impl->propertyIndex(name) >= 0;
}

void DASHTree::getAtomNodePath(const ROMol &mol, unsigned int atomIdx,
                               std::vector<std::uint32_t> &nodePath,
                               const DASHParams &params,
                               bool sourceNodeIds) const {
  if (atomIdx >= mol.getNumAtoms()) {
    throw ValueErrorException("atom index " + std::to_string(atomIdx) +
                              " is out of range");
  }
  if (sourceNodeIds && !d_impl->d_nodeSourceId) {
    throw ValueErrorException(
        "'" + d_impl->d_filename +
        "' carries no source node ids; re-run tools/dash_convert.py without "
        "--no-source-ids to compare paths against the python implementation");
  }
  MolMatcher matcher(*d_impl);
  matcher.setMolecule(mol);
  const std::uint32_t branch = matcher.match(atomIdx, params);

  // Led by the branch index, so this reads like the node paths the python
  // implementation reports. The node ids themselves are the container's own
  // breadth-first numbering unless the caller asks for the source numbering.
  nodePath.clear();
  nodePath.push_back(branch);
  const std::uint32_t root = d_impl->d_branchRoot[branch];
  for (const auto node : matcher.path()) {
    nodePath.push_back(sourceNodeIds ? d_impl->d_nodeSourceId[node]
                                     : node - root);
  }
}

double DASHTree::getAtomProperty(const ROMol &mol, unsigned int atomIdx,
                                 const std::string &property,
                                 const DASHParams &params) const {
  if (atomIdx >= mol.getNumAtoms()) {
    throw ValueErrorException("atom index " + std::to_string(atomIdx) +
                              " is out of range");
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
    throw ValueErrorException("branch " + std::to_string(branch) +
                              " is out of range; the tree has " +
                              std::to_string(nBranches));
  }
  // branches tile the node array in order, so one owns every id up to the
  // next root
  const std::uint32_t root = d_impl->d_branchRoot[branch];
  const std::uint32_t end = branch + 1 < nBranches
                                ? d_impl->d_branchRoot[branch + 1]
                                : d_impl->d_header.numNodes;
  if (nodeId >= end - root) {
    throw ValueErrorException("node " + std::to_string(nodeId) +
                              " is out of range; branch " +
                              std::to_string(branch) + " has " +
                              std::to_string(end - root) + " nodes");
  }
  return DASHTreeNode(*d_impl, root + nodeId);
}

void DASHTree::getMatchedSubstructure(const ROMol &mol, unsigned int atomIdx,
                                      std::vector<unsigned int> &atoms,
                                      const DASHParams &params) const {
  if (atomIdx >= mol.getNumAtoms()) {
    throw ValueErrorException("atom index " + std::to_string(atomIdx) +
                              " is out of range");
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

int DASHTreeNode::getAtomFeatureIndex() const {
  return keyAtomType(getKey());
}

const AtomFeature &DASHTreeNode::getFeature() const {
  const int index = getAtomFeatureIndex();
  if (index >= static_cast<int>(numAtomFeatures)) {
    throw ValueErrorException("node " + std::to_string(getId()) +
                              " names atom feature " + std::to_string(index) +
                              ", which does not exist");
  }
  return getAtomFeatureTable()[index];
}

int DASHTreeNode::getConAtom() const { return keyConAtom(getKey()); }

int DASHTreeNode::getConType() const { return keyConType(getKey()); }

unsigned int DASHTreeNode::getNumChildren() const {
  return dp_impl->d_node[d_abs].numChildren;
}

DASHTreeNode DASHTreeNode::getChild(unsigned int i) const {
  const NodeRecord &record = dp_impl->d_node[d_abs];
  if (i >= record.numChildren) {
    throw ValueErrorException(
        "node " + std::to_string(getId()) + " of branch " +
        std::to_string(getBranch()) + " has " +
        std::to_string(record.numChildren) + " children, no child " +
        std::to_string(i));
  }
  return DASHTreeNode(*dp_impl, record.firstChild + i);
}

float DASHTreeNode::getAttention() const { return dp_impl->d_nodeAttn[d_abs]; }

bool DASHTreeNode::stops() const {
  return dp_impl->d_node[d_abs].flags & nodeFlagStop;
}

double DASHTreeNode::getValue(const std::string &property) const {
  return dp_impl->property(property).value(d_abs);
}

}  // namespace DASH
}  // namespace RDKit
