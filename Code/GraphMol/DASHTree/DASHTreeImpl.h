//
//  Copyright (C) 2026 Marc Lehner and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//! \file DASHTreeImpl.h
/// \brief internals of the .dash container and the subgraph matcher.
/// This header is not installed. Include it only from within the DASHTree
/// library.

#ifndef RD_DASHTREE_IMPL_H
#define RD_DASHTREE_IMPL_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/Invariant.h>

#include "MatchKey.h"
#include "DASHTree.h"

namespace RDKit {
namespace DASH {

constexpr char dashMagic[] = "DASHTREE";
constexpr std::size_t dashMagicLen = 8;
constexpr std::uint32_t dashFormatVersion = 3;
constexpr std::size_t dashHeaderSize = 128;
constexpr std::size_t dashPropertyEntrySize = 48;
constexpr std::size_t dashPropertyNameLen = 32;

//! written into the header so that the reader refuses a byte-swapped file
constexpr std::uint32_t dashEndianId = 0xDEADBEEF;

//! \brief storage type of a property column
/*!
  A column keeps whatever width its source needs and no more. The published
  trees use three of them. The MBIS charges and their deviations are float16.
  The size column, the number of training atoms that reached each node, is
  int32. The DASH-properties columns (AM1BCC, DFTD4, ...) are float64 values
  that do not survive narrowing.
*/
enum PropertyDType : std::uint8_t {
  DTYPE_FLOAT16 = 1,
  DTYPE_FLOAT32 = 2,
  DTYPE_INT32 = 3,
  DTYPE_FLOAT64 = 4
};

//! what an int32 column stores in place of NaN for a node without a value
constexpr std::int32_t missingInt32 = std::numeric_limits<std::int32_t>::min();

//! bit 0 of NodeRecord::flags, set on a node where a descent at the default
//! threshold stops
constexpr std::uint8_t nodeFlagStop = 0x01;

//! \brief one tree node
/*!
  One 8-byte load gives a descent step everything it needs about a node: how
  the node is matched, and where its children are. Nodes are numbered
  breadth-first. This makes the children a contiguous range and removes the
  child-index array that a general graph needs.
*/
struct NodeRecord {
  std::uint16_t key;         //!< packed match key, see MatchKey.h
  std::uint8_t numChildren;  //!< at most 47 in the published trees
  std::uint8_t flags;        //!< see nodeFlagStop
  std::uint32_t firstChild;  //!< node id of the first child
};

static_assert(sizeof(NodeRecord) == 8, "NodeRecord must stay 8 bytes");
static_assert(offsetof(NodeRecord, firstChild) == 4, "NodeRecord layout");

//! \brief the container header, in its exact layout on disk
struct DASHHeader {
  char magic[8];
  std::uint32_t version;
  std::uint32_t endianId;
  std::uint32_t numBranches;
  std::uint32_t numProperties;
  std::uint32_t numNodes;
  std::uint32_t maxChildren;
  //! the cumulative-attention threshold nodeFlagStop was precomputed for
  double defaultAttentionThreshold;
  std::uint64_t offBranchRoot;
  std::uint64_t offNodeRecord;
  std::uint64_t offNodeAttn;
  std::uint64_t reserved;  //!< written as 0, ignored
  std::uint64_t offPropertyDir;
  std::uint64_t fileSize;
  std::int16_t maxAtomType;
  std::int16_t maxConAtom;
  std::int16_t maxConType;
  std::uint16_t numFeatures;      //!< rows of the atom-feature table
  std::uint64_t offFeatureTable;  //!< numFeatures AtomFeature records
};

static_assert(offsetof(DASHHeader, version) == 8, "DASH header layout");
static_assert(offsetof(DASHHeader, defaultAttentionThreshold) == 32,
              "DASH header layout");
static_assert(offsetof(DASHHeader, offBranchRoot) == 40, "DASH header layout");
static_assert(offsetof(DASHHeader, fileSize) == 80, "DASH header layout");
static_assert(offsetof(DASHHeader, maxAtomType) == 88, "DASH header layout");
static_assert(offsetof(DASHHeader, numFeatures) == 94, "DASH header layout");
static_assert(offsetof(DASHHeader, offFeatureTable) == 96,
              "DASH header layout");
static_assert(sizeof(DASHHeader) <= dashHeaderSize, "DASH header too big");

//! \brief decodes one IEEE-754 binary16 value
/*!
  The published trees store their values as float16, so the container keeps them
  that way. This is lossless (float16 -> double -> float16 round-trips exactly)
  and halves the pages a lookup has to touch. It is written by hand because
  std::float16_t is C++23, and only decoding is needed.
*/
inline double decodeFloat16(std::uint16_t h) {
  const int exponent = (h >> 10) & 0x1F;
  const int mantissa = h & 0x3FF;
  double value;
  if (exponent == 0x1F) {
    value = mantissa ? std::numeric_limits<double>::quiet_NaN()
                     : std::numeric_limits<double>::infinity();
  } else if (exponent == 0) {
    value = std::ldexp(mantissa, -24);
  } else {
    value = std::ldexp(mantissa | 0x400, exponent - 25);
  }
  return (h & 0x8000) ? -value : value;
}

//! one resolved property column
struct PropertyColumn {
  const void *data = nullptr;
  PropertyDType dtype = DTYPE_FLOAT32;

  //! the value stored for a node, as a double
  double value(std::uint32_t node) const {
    switch (dtype) {
      case DTYPE_FLOAT16:
        return decodeFloat16(static_cast<const std::uint16_t *>(data)[node]);
      case DTYPE_FLOAT32:
        return static_cast<const float *>(data)[node];
      case DTYPE_INT32: {
        const std::int32_t value =
            static_cast<const std::int32_t *>(data)[node];
        return value == missingInt32 ? std::numeric_limits<double>::quiet_NaN()
                                     : value;
      }
      case DTYPE_FLOAT64:
        return static_cast<const double *>(data)[node];
    }
    return std::numeric_limits<double>::quiet_NaN();
  }
};

//! \brief classifies atoms against the feature table of one tree
/*!
  A lookup is O(1) through a direct-lookup table over the 14 bits the five
  components fit into (Z <= 63, degree <= 7, charge -1..2, numHs <= 3). The
  table is built once, when the tree is mapped. The constructor refuses a table
  with a tuple outside that budget.
*/
class AtomFeatureLookup {
 public:
  AtomFeatureLookup() = default;
  //! \throws ValueErrorException if a tuple cannot be encoded, is a hydrogen of
  //! a degree other than one, or occurs twice
  AtomFeatureLookup(const AtomFeature *features, unsigned int numFeatures);

  unsigned int size() const { return d_numFeatures; }
  //! \pre index < size()
  const AtomFeature &feature(unsigned int index) const {
    return d_features[index];
  }
  //! branch index of an atom in a molecule, or -1
  int index(const Atom *atom) const;

 private:
  const AtomFeature *d_features = nullptr;
  unsigned int d_numFeatures = 0;
  std::vector<std::int16_t> d_index;
};

//! \brief the mapped container
class DASHTree::Impl {
 public:
  Impl(const std::string &filename, const std::vector<std::string> &properties);

  //! index of a resolved column, or -1
  int propertyIndex(const std::string &name) const;
  //! throws a descriptive ValueErrorException if \p name was not resolved
  const PropertyColumn &property(const std::string &name) const;
  //! names of every column the file holds, resolved or not
  const std::vector<std::string> &allPropertyNames() const {
    return d_allPropertyNames;
  }

  //! one past the last node id \p branch owns
  std::uint32_t branchEnd(std::uint32_t branch) const {
    return branch + 1 < d_header.numBranches ? d_branchRoot[branch + 1]
                                             : d_header.numNodes;
  }
  //! \brief the child id range of \p node, which must lie after it and before
  //! \p end, the end of its branch
  /*!
    The check runs at use, not at load. A check at load would read every
    record, but the tree is mapped so that only the visited records are read.
    \throws ValueErrorException if the file is corrupt
  */
  std::pair<std::uint32_t, std::uint32_t> children(std::uint32_t node,
                                                   std::uint32_t end) const {
    const NodeRecord &record = d_node[node];
    const std::uint64_t last =
        static_cast<std::uint64_t>(record.firstChild) + record.numChildren;
    if (record.numChildren && (record.firstChild <= node || last > end)) {
      throwCorruptChildren(node);
    }
    return {record.firstChild, static_cast<std::uint32_t>(last)};
  }
  [[noreturn]] void throwCorruptChildren(std::uint32_t node) const;

  std::string d_filename;
  boost::iostreams::mapped_file_source d_map;
  DASHHeader d_header;

  const std::uint32_t *d_branchRoot = nullptr;  //!< numBranches
  const NodeRecord *d_node = nullptr;           //!< numNodes
  const float *d_nodeAttn = nullptr;            //!< numNodes
  //! the tree's atom-feature classes, read from its table
  AtomFeatureLookup d_features;

  std::vector<PropertyColumn> d_properties;
  std::vector<std::string> d_propertyNames;     //!< of d_properties
  std::vector<std::string> d_allPropertyNames;  //!< of every column in the file
};

//! \brief matches the atoms of one molecule against a tree
/*!
  The matcher holds the per-molecule tables and the scratch the descent needs.
  It computes the neighbour features of a molecule once, not once per atom. A
  thread that processes a batch reuses the buffers.

  One matcher belongs to one thread. It borrows the tree, which is immutable.
*/
class MolMatcher {
 public:
  explicit MolMatcher(const DASHTree::Impl &tree) : d_tree(tree), d_live(1) {}

  //! \brief prepares for \p mol: feature classes, neighbour lists, partial keys
  /*!
    \throws ValueErrorException if the molecule's property cache is stale, if
            an atom has implicit hydrogens, or if an atom is in none of the
            tree's atom-feature classes. The reference python implementation
            fails at the same point.
  */
  void setMolecule(const ROMol &mol);

  //! \brief descends the tree for one atom
  /*!
    Fills path() with node ids, root first. Returns the branch index. A hydrogen
    whose heavy neighbour has no matching child yields an empty path.
  */
  std::uint32_t match(unsigned int atomIdx, const DASHParams &params);

  //! node ids of the last match, root first
  const std::vector<std::uint32_t> &path() const { return d_path; }
  //! molecule atoms of the last match, in the order the descent added them.
  //! Where atoms stayed tied to the end of the descent, this is one of the
  //! equivalent matches.
  const std::vector<std::uint32_t> &subgraph() const { return d_live[0].atoms; }

  //! \brief the deepest value on the last match's path, NaN if there is none
  /*!
    Walks the path in reverse, so that a sparse column can use the value of a
    shallower, more general substructure.
  */
  double pathValue(const PropertyColumn &column) const;

 private:
  //! one way to lay the path matched so far onto the molecule
  struct Embedding {
    std::vector<std::uint32_t> atoms;     //!< matched atoms, insertion order
    std::vector<std::uint16_t> candKey;   //!< noMatchKey once consumed
    std::vector<std::uint32_t> candAtom;  //!< atoms that could join next
    //! one bit per molecule atom: matched, or already a candidate
    std::vector<std::uint64_t> seen;
    //! bloom filter over the live candidate keys. It rejects most children in
    //! one test.
    std::uint64_t candFilter = 0;

    //! marks \p atom seen and returns whether it already was
    bool see(std::uint32_t atom) {
      const std::uint64_t bit = 1ull << (atom & 63u);
      const bool was = seen[atom >> 6] & bit;
      seen[atom >> 6] |= bit;
      return was;
    }
  };

  //! descends one level. Returns false if no child of the last node matches.
  bool step(bool legacyTieBreaking);
  //! fills d_matches with the live candidates that carry \p key
  void findMatches(std::uint16_t key, bool legacyTieBreaking);
  //! makes one embedding of the next level per entry of d_matches
  void advance();
  //! binds candidate \p cand of \p embedding
  void extend(Embedding &embedding, std::uint32_t cand) const;
  //! appends the not-yet-seen neighbours of the atom matched at \p position
  void addCandidates(Embedding &embedding, std::uint32_t position) const;

  const DASHTree::Impl &d_tree;
  unsigned int d_numAtoms = 0;
  std::uint32_t d_branchEnd = 0;  //!< of the current branch

  // per-molecule tables
  std::vector<std::int32_t> d_feature;    //!< branch index per atom
  std::vector<std::uint32_t> d_nbrStart;  //!< numAtoms + 1, CSR rows
  std::vector<std::uint32_t> d_nbrAtom;   //!< 2 * numBonds
  std::vector<std::uint16_t> d_nbrKey;    //!< 2 * numBonds, partial match keys

  // per-atom scratch. The live embeddings are the first d_numLive of d_live,
  // and the next level's embeddings are built in d_next. The unused ones keep
  // their buffers for reuse.
  std::vector<std::uint32_t> d_path;
  std::vector<Embedding> d_live;
  std::vector<Embedding> d_next;
  std::size_t d_numLive = 1;
  //! (embedding, candidate) pairs that take the next node
  std::vector<std::pair<std::uint32_t, std::uint32_t>> d_matches;
};

}  // namespace DASH
}  // namespace RDKit

#endif
