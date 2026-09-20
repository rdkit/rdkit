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
/// Not installed; include only from within the DASHTree library.

#ifndef RD_DASHTREE_IMPL_H
#define RD_DASHTREE_IMPL_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/Invariant.h>

#include "AtomFeatures.h"
#include "DASHTree.h"

namespace RDKit {
namespace DASH {
namespace detail {

const char *const dashMagic = "DASHTREE";
const std::size_t dashMagicLen = 8;
const std::size_t dashHeaderSize = 128;
const std::size_t dashPropertyEntrySize = 48;
const std::size_t dashPropertyNameLen = 32;

//! written into the header so a byte-swapped read is caught rather than acted on
const std::uint32_t dashEndianId = 0xDEADBEEF;

//! \brief storage type of a property column
/*!
  A column keeps whatever width its source needs and no more. The published
  trees use all four: the MBIS charges and their deviations are float16, the
  subtree sizes int32, and the DASH-properties columns (AM1BCC, DFTD4, ...)
  genuine float64 that does not survive narrowing.
*/
enum PropertyDType : std::uint8_t {
  DTYPE_FLOAT16 = 1,
  DTYPE_FLOAT32 = 2,
  DTYPE_INT32 = 3,
  DTYPE_FLOAT64 = 4
};

//! bit 0 of NodeRecord::flags: the descent stops here at the default threshold
const std::uint8_t nodeFlagStop = 0x01;

//! \brief one tree node
/*!
  Everything a descent step needs about a node, in one 8-byte load: how it is
  matched, and where its children are. Nodes are numbered breadth-first, which
  is what makes the children a contiguous range and so removes the child-index
  array a general graph would need.
*/
struct NodeRecord {
  std::uint16_t key;           //!< packed match key, see AtomFeatures.h
  std::uint8_t numChildren;    //!< at most 47 in the published trees
  std::uint8_t flags;          //!< see nodeFlagStop
  std::uint32_t firstChild;    //!< node id of the first child
};

static_assert(sizeof(NodeRecord) == 8, "NodeRecord must stay 8 bytes");
static_assert(offsetof(NodeRecord, firstChild) == 4, "NodeRecord layout");

//! \brief the container header, laid out exactly as it is on disk
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
  std::uint64_t offNodeSourceId;  //!< 0 when the file carries no source ids
  std::uint64_t offPropertyDir;
  std::uint64_t fileSize;
  std::int16_t maxAtomType;
  std::int16_t maxConAtom;
  std::int16_t maxConType;
};

static_assert(offsetof(DASHHeader, version) == 8, "DASH header layout");
static_assert(offsetof(DASHHeader, defaultAttentionThreshold) == 32,
              "DASH header layout");
static_assert(offsetof(DASHHeader, offBranchRoot) == 40, "DASH header layout");
static_assert(offsetof(DASHHeader, fileSize) == 80, "DASH header layout");
static_assert(offsetof(DASHHeader, maxAtomType) == 88, "DASH header layout");
static_assert(sizeof(DASHHeader) <= dashHeaderSize, "DASH header too big");

//! \brief decodes one IEEE-754 binary16 value
/*!
  The published trees store their values as float16, so the container keeps them
  that way: it is lossless (float16 -> float32 -> float16 round-trips exactly)
  and halves the pages a lookup has to touch. Written out by hand because
  std::float16_t is C++23 and this is the whole of what is needed.
*/
inline float decodeFloat16(std::uint16_t h) {
  const std::uint32_t sign = static_cast<std::uint32_t>(h & 0x8000u) << 16;
  const std::uint32_t exponent = (h >> 10) & 0x1Fu;
  const std::uint32_t mantissa = h & 0x3FFu;
  std::uint32_t bits;
  if (exponent == 0) {
    if (mantissa == 0) {
      bits = sign;  // +-0
    } else {
      // subnormal: renormalise into a float32 normal
      std::uint32_t e = 0;
      std::uint32_t m = mantissa;
      while (!(m & 0x400u)) {
        m <<= 1;
        ++e;
      }
      bits = sign | ((127 - 15 - e + 1) << 23) | ((m & 0x3FFu) << 13);
    }
  } else if (exponent == 0x1Fu) {
    bits = sign | 0x7F800000u | (mantissa << 13);  // inf / NaN
  } else {
    bits = sign | ((exponent + 127 - 15) << 23) | (mantissa << 13);
  }
  float out;
  std::memcpy(&out, &bits, sizeof(out));
  return out;
}

//! one resolved property column
struct PropertyColumn {
  std::string name;
  const void *data = nullptr;
  PropertyDType dtype = DTYPE_FLOAT32;

  //! the value stored for a node, as a double
  double value(std::uint32_t node) const {
    switch (dtype) {
      case DTYPE_FLOAT16:
        return decodeFloat16(static_cast<const std::uint16_t *>(data)[node]);
      case DTYPE_FLOAT32:
        return static_cast<const float *>(data)[node];
      case DTYPE_INT32:
        return static_cast<const std::int32_t *>(data)[node];
      case DTYPE_FLOAT64:
        return static_cast<const double *>(data)[node];
    }
    return std::numeric_limits<double>::quiet_NaN();
  }
};

}  // namespace detail

//! \brief the mapped container
class DASHTree::Impl {
 public:
  Impl(const std::string &filename, const std::vector<std::string> &properties,
       bool prefetch);

  //! index of a resolved column, or -1
  int propertyIndex(const std::string &name) const;
  //! throws a descriptive ValueErrorException if \p name was not resolved
  const detail::PropertyColumn &property(const std::string &name) const;
  //! names of every column the file holds, resolved or not
  std::vector<std::string> allPropertyNames() const;

  std::string d_filename;
  boost::iostreams::mapped_file_source d_map;
  detail::DASHHeader d_header;

  const std::uint32_t *d_branchRoot = nullptr;      //!< numBranches
  const detail::NodeRecord *d_node = nullptr;       //!< numNodes
  const float *d_nodeAttn = nullptr;                //!< numNodes
  const std::uint32_t *d_nodeSourceId = nullptr;    //!< numNodes, may be null

  std::vector<detail::PropertyColumn> d_properties;
};

//! \brief matches the atoms of one molecule against a tree
/*!
  Holds the per-molecule tables and the scratch the descent needs, so a
  molecule's neighbour features are computed once rather than once per atom, and
  a thread working through a batch reuses its buffers.

  One matcher belongs to one thread. It borrows the tree, which is immutable.
*/
class MolMatcher {
 public:
  explicit MolMatcher(const DASHTree::Impl &tree) : d_tree(tree) {}

  //! \brief prepares for \p mol: feature classes, neighbour lists, partial keys
  /*!
    \throws ValueErrorException if the molecule's property cache is stale, or if
            any atom's feature tuple is not a DASH class -- the same point at
            which the reference python implementation gives up
  */
  void setMolecule(const ROMol &mol);

  //! \brief descends the tree for one atom
  /*!
    Fills path() with node ids, root first. Returns the branch index. A hydrogen
    whose heavy neighbour has no matching child yields an empty path, which is
    the one case where no value can be read.
  */
  std::uint32_t match(unsigned int atomIdx, const DASHParams &params);

  //! node ids of the last match, root first
  const std::vector<std::uint32_t> &path() const { return d_path; }
  //! molecule atoms of the last match, in the order the descent added them
  const std::vector<std::uint32_t> &subgraph() const { return d_subgraph; }

  //! \brief the deepest value on the last match's path, NaN if there is none
  /*!
    Walks the path in reverse so a sparsely populated column falls back to a
    shallower, more general substructure.
  */
  double pathValue(const detail::PropertyColumn &column) const;

  unsigned int numAtoms() const { return d_numAtoms; }

 private:
  //! first child of \p node whose key matches a live candidate
  bool pickChild(std::uint32_t node, std::uint32_t &childNode,
                 std::uint32_t &childAtom);
  //! appends the not-yet-seen neighbours of the subgraph atom at \p position
  void addCandidates(std::uint32_t position);
  //! retires the candidate for \p atomIdx, which has just joined the subgraph
  void retireCandidate(std::uint32_t atomIdx);

  const DASHTree::Impl &d_tree;
  unsigned int d_numAtoms = 0;

  // per-molecule tables
  std::vector<std::int32_t> d_feature;    //!< branch index per atom
  std::vector<std::uint32_t> d_nbrStart;  //!< numAtoms + 1, CSR rows
  std::vector<std::uint32_t> d_nbrAtom;   //!< 2 * numBonds
  std::vector<std::uint16_t> d_nbrKey;    //!< 2 * numBonds, partial match keys

  // per-atom scratch
  std::vector<std::uint32_t> d_path;
  std::vector<std::uint32_t> d_subgraph;  //!< matched atoms, insertion order
  std::vector<std::uint16_t> d_candKey;   //!< noMatchKey once consumed
  std::vector<std::uint32_t> d_candAtom;
  std::vector<std::uint32_t> d_seen;      //!< generation stamps, O(1) to clear
  std::uint32_t d_generation = 0;
  //! bloom filter over the live candidate keys; rejects most children in one test
  std::uint64_t d_candFilter = 0;
};

}  // namespace DASH
}  // namespace RDKit

#endif
