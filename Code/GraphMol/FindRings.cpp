//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/rdmol_throw.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Rings.h>
#include <GraphMol/RDMol.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>

#include <RDGeneral/utils.h>
#include <bit>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>

using RINGINVAR = boost::dynamic_bitset<>;
// TODO: RINGINVAR_SET is never removed from, so replace it with a vector of
// bit set data and a set of integers indexing into it. To compare and insert
// if not present, add to the end of the vector, pass that new index into insert
// and check return. If not inserted, resize the vector back down.
//using RINGINVAR_SET = std::set<RINGINVAR>;
//using RINGINVAR_VECT = std::vector<RINGINVAR>;

namespace RingUtils {
const size_t MAX_BFSQ_SIZE = 200000;  // arbitrary huge value

using namespace RDKit;

RINGINVAR computeRingInvariant(INT_VECT ring, unsigned int numAtoms) {
  boost::dynamic_bitset<> res(numAtoms);
  for (auto idx : ring) {
    res.set(idx);
  }
  return res;
}
void computeRingInvariant(RINGINVAR &res, const UINT_VECT &ringAtoms,
                          uint32_t begin, uint32_t end, unsigned int numAtoms) {
  if (res.size() != numAtoms) {
    res.resize(numAtoms);
  }
  res.reset();
  for (uint32_t i = begin; i < end; ++i) {
    res.set(ringAtoms[i]);
  }
}
void computeRingInvariant(uint64_t *res, const atomindex_t *ringAtoms,
                          uint32_t begin, uint32_t end) {
  for (uint32_t i = begin; i < end; ++i) {
    const atomindex_t atom = ringAtoms[i];
    const uint32_t index = (atom >> 6);
    const uint64_t mask = (uint64_t(1) << (atom & 0x3F));
    res[index] |= mask;
  }
}

struct RingInvarSet {
  std::vector<uint64_t> bits;
  size_t elementsPerInvar;

  RingInvarSet(uint32_t numBits) : elementsPerInvar((numBits + 63) / 64) {}

  size_t tempAddOne() {
    const size_t oldSize = bits.size();
    bits.resize(oldSize + elementsPerInvar, 0);
    return oldSize;
  }
  // This does a linear search. If this becomes a bottleneck for large numbers of rings,
  // have a switch over when the number of rings passes a threshold.
  std::pair<size_t, bool> find(const uint64_t *bitsToFind, size_t internalEnd) const {
    const uint64_t *bitsData = bits.data();
    const uint64_t *const newBits = bitsData + internalEnd;
    size_t index = 0;
    for (; bitsData < newBits; bitsData += elementsPerInvar, ++index) {
      bool equal = true;
      for (size_t offset = 0; offset < elementsPerInvar; ++offset) {
        if (bitsData[offset] != bitsToFind[offset]) {
          equal = false;
          break;
        }
      }
      if (equal) {
        return std::make_pair(index, false);
      }
    }
    return std::make_pair(index, true);
  }

  std::pair<size_t, bool> insert(const UINT_VECT &ringAtoms, uint32_t begin, uint32_t end) {
    const size_t oldSize = tempAddOne();
    uint64_t *bitsData = bits.data();
    uint64_t *const newBits = bitsData + oldSize;
    computeRingInvariant(newBits, ringAtoms.data(), begin, end);
    auto ret = find(newBits, oldSize);
    if (!ret.second) {
      bits.resize(oldSize);
    }
    return ret;
  }

  bool less(const size_t aIdx, const size_t bIdx) const {
    const uint64_t *a = bits.data() + aIdx * elementsPerInvar;
    const uint64_t *b = bits.data() + bIdx * elementsPerInvar;
    for (size_t i = elementsPerInvar; i > 0; ) {
      --i;
      if (a[i] != b[i]) {
        return (a[i] < b[i]);
      }
    }
    return false;
  }
};

void convertToBonds(const atomindex_t *atoms,
                    const uint32_t ringSize, uint32_t *bonds,
                    const RDMol &mol) {
  for (unsigned int i = 0; i < (ringSize - 1); i++) {
    const uint32_t bnd = mol.getBondIndexBetweenAtoms(atoms[i], atoms[i + 1]);
    if (bnd == std::numeric_limits<std::uint32_t>::max()) {
      throw ValueErrorException("expected bond not found");
    }
    bonds[i] = bnd;
  }
  // bond from last to first atom
  const uint32_t bnd =
      mol.getBondIndexBetweenAtoms(atoms[ringSize - 1], atoms[0]);
  if (bnd == std::numeric_limits<std::uint32_t>::max()) {
    throw ValueErrorException("expected bond not found");
  }

  bonds[ringSize - 1] = bnd;
}

void convertToBonds(const std::vector<uint32_t> &begins,
                    const std::vector<atomindex_t> &atoms,
                    std::vector<uint32_t> &bonds, const RDMol &mol) {
  bonds.resize(atoms.size());
  for (uint32_t ringIdx = 0, numRings = begins.size() - 1;
      ringIdx < numRings; ++ringIdx) {
    uint32_t begin = begins[ringIdx];
    uint32_t end = begins[ringIdx + 1];
    convertToBonds(atoms.data() + begin, end - begin, bonds.data() + begin, mol);
  }
}

void convertToBonds(const INT_VECT &ring, INT_VECT &bondRing,
                    const ROMol &mol) {
  bondRing.resize(ring.size());
  static_assert(sizeof(ring[0]) == sizeof(uint32_t) &&
                sizeof(bondRing[0]) == sizeof(uint32_t));
  RingUtils::convertToBonds(
      reinterpret_cast<const uint32_t *>(ring.data()), ring.size(),
      reinterpret_cast<uint32_t *>(bondRing.data()), mol.asRDMol());
}

void convertToBonds(const VECT_INT_VECT &res, VECT_INT_VECT &brings,
                    const ROMol &mol) {
  for (const auto &ring : res) {
    INT_VECT bring;
    convertToBonds(ring, bring, mol);
    brings.push_back(bring);
  }
}

}  // end of namespace RingUtils

namespace FindRings {
using namespace RDKit;

// An optimization to create a memory workspace that gets reused
class BFSWorkspace {
 public:
  void smallestRingsBfs(const RDMol &mol, atomindex_t root,
                        UINT_VECT &ringBegins, UINT_VECT &ringAtoms,
                        const detail::BitSetWrapper &activeBonds,
                        UINT_VECT *forbidden = nullptr);

 private:
  INT_VECT d_parents;
  std::vector<unsigned int> d_depths;
  INT_VECT d_done;

  // This is a ring buffer queue
  std::vector<atomindex_t> d_bfsq;

  UINT_VECT d_ring;
};

void trimBonds(const uint32_t cand, const RDMol &tMol,
               std::vector<uint32_t> *changed, INT_VECT &atomDegrees,
               detail::BitSetWrapper &activeBonds);
void storeRingInfo(const ROMol &mol, const INT_VECT &ring) {
  INT_VECT bondIndices;
  RingUtils::convertToBonds(ring, bondIndices, mol);
  mol.getRingInfo()->addRing(ring, bondIndices);
}
void storeRingsInfo(const ROMol &mol, const VECT_INT_VECT &rings) {
  for (const auto &ring : rings) {
    storeRingInfo(mol, ring);
  }
}

void invertRingData(const std::vector<uint32_t> &ringBegins,
                    const std::vector<uint32_t> &ringContents,
                    std::vector<uint32_t> &membershipBegins,
                    std::vector<uint32_t> &memberships, uint32_t contentLimit) {
  membershipBegins.clear();
  membershipBegins.resize(contentLimit + 1, 0);

  // First, count, but write the results one place later
  for (uint32_t v : ringContents) {
    ++membershipBegins[v + 1];
  }

  // Partial sum, but write result another place later
  uint32_t sum = 0;
  for (uint32_t i = 1, n = membershipBegins.size(); i < n; ++i) {
    const uint32_t v = membershipBegins[i];
    membershipBegins[i] = sum;
    sum += v;
  }

  // Fill in memberships while incrementing membershipBegins
  memberships.resize(ringContents.size());
  const uint32_t numRings = ringBegins.size() - 1;
  for (uint32_t ringIdx = 0; ringIdx < numRings; ++ringIdx) {
    uint32_t begin = ringBegins[ringIdx];
    uint32_t end = ringBegins[ringIdx + 1];
    for (; begin != end; ++begin) {
      const uint32_t v = ringContents[begin];
      uint32_t &index = membershipBegins[v + 1];
      memberships[index] = ringIdx;
      ++index;
    }
  }
}

void buildRingInfoFromAtoms(const RDMol &mol, RingInfoCache &info) {
  // Fill in bondsInRings first
  RingUtils::convertToBonds(info.ringBegins, info.atomsInRings, info.bondsInRings, mol);

  // Fill in memberships next
  invertRingData(info.ringBegins, info.atomsInRings, info.atomMembershipBegins,
                 info.atomMemberships, mol.getNumAtoms());
  invertRingData(info.ringBegins, info.bondsInRings, info.bondMembershipBegins,
                 info.bondMemberships, mol.getNumBonds());

  // Lastly, fill in fusing data
  info.initFusedInfoFromBondMemberships();
}

void markUselessD2s(unsigned int root, const ROMol &tMol,
                    boost::dynamic_bitset<> &forb, const INT_VECT &atomDegrees,
                    const boost::dynamic_bitset<> &activeBonds) {
  // recursive function to mark any degree 2 nodes that are already represented
  // by root for the purpose of finding smallest rings.
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = tMol.getAtomBonds(tMol.getAtomWithIdx(root));
  while (beg != end) {
    const Bond *bond = tMol[*beg];
    ++beg;
    if (!activeBonds[bond->getIdx()]) {
      continue;
    }
    unsigned int oIdx = bond->getOtherAtomIdx(root);
    if (!forb[oIdx] && atomDegrees[oIdx] == 2) {
      forb[oIdx] = 1;
      markUselessD2s(oIdx, tMol, forb, atomDegrees, activeBonds);
    }
  }
}

void pickD2Nodes(const ROMol &tMol, INT_VECT &d2nodes, const INT_VECT &currFrag,
                 const INT_VECT &atomDegrees,
                 const boost::dynamic_bitset<> &activeBonds) {
  d2nodes.resize(0);

  // forb contains all d2 nodes, not just the ones we want to keep
  boost::dynamic_bitset<> forb(tMol.getNumAtoms());
  while (1) {
    int root = -1;
    for (int axci : currFrag) {
      if (atomDegrees[axci] == 2 && !forb[axci]) {
        root = axci;
        d2nodes.push_back(axci);
        forb[axci] = 1;
        break;
      }
    }
    if (root == -1) {
      break;
    } else {
      markUselessD2s(root, tMol, forb, atomDegrees, activeBonds);
    }
  }
}

void markUselessD2s(uint32_t rootAtomIndex, const RDMol &tMol,
                    detail::BitSetWrapper &forb, const INT_VECT &atomDegrees,
                    const detail::BitSetWrapper &activeBonds) {
  // recursive function to mark any degree 2 nodes that are already represented
  // by root for the purpose of finding smallest rings.
  auto [beg, end] = tMol.getAtomBonds(rootAtomIndex);
  for (; beg != end; ++beg) {
    const uint32_t bondIndex = *beg;
    if (!activeBonds[bondIndex]) {
      continue;
    }
    const BondData &bond = tMol.getBond(bondIndex);
    const uint32_t oIdx = bond.getOtherAtomIdx(rootAtomIndex);
    if (!forb[oIdx] && atomDegrees[oIdx] == 2) {
      forb.set(oIdx);
      markUselessD2s(oIdx, tMol, forb, atomDegrees, activeBonds);
    }
  }
}

void pickD2Nodes(const RDMol &tMol, UINT_VECT &d2nodes,
                 const uint32_t *currFragBegin, const uint32_t *currFragEnd,
                 const INT_VECT &atomDegrees,
                 const detail::BitSetWrapper &activeBonds,
                 detail::BitSetWrapper &forb) {
  d2nodes.resize(0);

  // forb contains all d2 nodes, not just the ones we want to keep
  forb.reset();
  for (auto it = currFragBegin; it != currFragEnd; ++it) {
    const uint32_t axci = *it;
    if (atomDegrees[axci] == 2 && !forb[axci]) {
      d2nodes.push_back(axci);
      forb.set(axci);
      // TODO: De-recurse markUselessD2s. There are at most 2 paths to traverse
      // along from axci.
      markUselessD2s(axci, tMol, forb, atomDegrees, activeBonds);
    }
  }
}

void findSSSRforDupCands(
    const RDMol &mol, UINT_VECT &resBegins, UINT_VECT &resAtoms,
    RingUtils::RingInvarSet &invars,
    const std::vector<std::pair<uint32_t, uint32_t>> &dupMap,
    const std::vector<uint32_t> &dupD2CandBegins,
    const std::vector<uint32_t> &dupD2Cands,
    const std::vector<std::pair<uint32_t, uint32_t>> &order, INT_VECT &atomDegrees,
    const detail::BitSetWrapper &activeBonds,
    // The rest are only to reuse memory allocations from the caller
    BFSWorkspace &bfs_workspace, UINT_VECT &sringBegins, UINT_VECT &sringAtoms,
    UINT_VECT &nringBegins, UINT_VECT &nringAtoms) {
  INT_VECT atomDegreesCopy;
  std::vector<uint64_t> activeBondsCopyStorage;
  detail::BitSetWrapper activeBondsCopy(activeBondsCopyStorage, activeBonds.size());
  // NOTE: The iteration order has changed now that dupD2Cands is a vector
  // instead of a map ordered on the corresponding ring invariant.
  // If this compatibility break is an issue, sort dupD2Cands by the ring invar
  // with the least significant bits being first (most significant bits last)
  const size_t numDupRings = dupD2CandBegins.size() - 1;
  for (size_t iterIdx = 0; iterIdx < numDupRings; ++iterIdx) {
    const size_t ringIdx = (order.size() > 0) ? order[iterIdx].second : 0;
    uint32_t dupCandsBegin = dupD2CandBegins[ringIdx];
    uint32_t dupCandsEnd = dupD2CandBegins[ringIdx + 1];
    PRECONDITION(dupCandsEnd >= dupCandsBegin + 2,
                 "Caller to findSSSRforDupCands included dup too small");

    // we have duplicate candidates.
    auto minSiz = static_cast<uint32_t>(std::numeric_limits<int>::max());
    nringBegins.clear();
    nringAtoms.clear();
    for (uint32_t dupCandIdx = dupCandsBegin; dupCandIdx != dupCandsEnd;
         ++dupCandIdx) {
      uint32_t dupCand = dupD2Cands[dupCandIdx];
      // now break bonds for all the d2 nodes for that give the same rings as
      // with (*dupi) and recompute smallest ring with (*dupi)
      atomDegreesCopy = atomDegrees;
      std::copy(activeBonds.data(), activeBonds.data() + activeBonds.dataSize(),
                activeBondsCopy.data());
      auto dmciBegin = std::lower_bound(
          dupMap.begin(), dupMap.end(), dupCand,
          [](const std::pair<uint32_t, uint32_t> &a, const uint32_t b) -> bool {
            return a.first < b;
          });
      auto dmciEnd = std::upper_bound(
          dupMap.begin(), dupMap.end(), dupCand,
          [](const uint32_t a, const std::pair<uint32_t, uint32_t> &b) -> bool {
            return a < b.first;
          });
      CHECK_INVARIANT(dmciEnd - dmciBegin >= 1, "duplicate could not be found");
      for (; dmciBegin != dmciEnd; ++dmciBegin) {
        uint32_t dni = dmciBegin->second;
        // atomDegreesCopy is read and written by trimBonds, so may still be
        // necessary
        trimBonds(dni, mol, nullptr, atomDegreesCopy, activeBondsCopy);
      }

      // now find the smallest ring/s around (*dupi)
      sringBegins.clear();
      sringAtoms.clear();
      bfs_workspace.smallestRingsBfs(mol, dupCand, sringBegins, sringAtoms,
                                     activeBondsCopy);

      // TODO: This could be done in-place in sringBegins and sringAtoms,
      // because smallestRingsBfs appends to the vectors, but for now,
      // stick with nringBegins and nringAtoms being separate.
      for (uint32_t ringIndex = 0, numRings = sringBegins.size();
           ringIndex < numRings; ++ringIndex) {
        const uint32_t begin = sringBegins[ringIndex];
        const uint32_t end = (ringIndex + 1 < numRings)
                                 ? sringBegins[ringIndex + 1]
                                 : sringAtoms.size();
        const uint32_t size = end - begin;

        // Only keep the smallest rings
        if (size > minSiz) {
          continue;
        } else if (size < minSiz) {
          minSiz = size;
          nringBegins.clear();
          nringAtoms.clear();
        }
        nringBegins.push_back(nringAtoms.size());
        nringAtoms.insert(nringAtoms.end(), sringAtoms.begin() + begin,
                          sringAtoms.begin() + end);
      }
    }
    for (uint32_t ringIndex = 0, numRings = nringBegins.size();
         ringIndex < numRings; ++ringIndex) {
      const uint32_t begin = nringBegins[ringIndex];
      const uint32_t end = (ringIndex + 1 < numRings)
                               ? nringBegins[ringIndex + 1]
                               : nringAtoms.size();
      PRECONDITION(end - begin == minSiz,
                   "The loop above should only have kept the smallest rings");
      if (invars.insert(nringAtoms, begin, end).second) {
        resAtoms.insert(resAtoms.end(), nringAtoms.begin() + begin,
                        nringAtoms.begin() + end);
        resBegins.push_back(resAtoms.size());
      }
    }  // end of loop over new rings found
  }    // end of loop over all set of duplicate candidates
}

void sortRingsBySize(std::vector<uint32_t> &resBegins,
                     std::vector<atomindex_t> &resAtoms,
                     size_t resBeginsFragBegin,
                     size_t resAtomsFragBegin,
                     std::vector<uint32_t> &tempBegins,
                     std::vector<atomindex_t> &tempAtoms,
                     std::vector<uint32_t> &order) {
  const uint32_t numRings = resBegins.size() - 1 - resBeginsFragBegin;

  // sort on size by sorting ring indices and then reordering
  order.resize(numRings);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&resBegins, resBeginsFragBegin](uint32_t a, uint32_t b) -> bool {
              const uint32_t aSize = resBegins[resBeginsFragBegin + a + 1] -
                                     resBegins[resBeginsFragBegin + a];
              const uint32_t bSize = resBegins[resBeginsFragBegin + b + 1] -
                                     resBegins[resBeginsFragBegin + b];
              return aSize < bSize;
            });

  std::vector<uint32_t> &reorderedBegins = tempBegins;
  std::vector<atomindex_t> &reorderedAtoms = tempAtoms;
  reorderedBegins.resize(numRings + 1);
  reorderedAtoms.resize(resAtoms.size() - resAtomsFragBegin);
  uint32_t outputAtomsIdx = 0;
  for (uint32_t ringIdx = 0; ringIdx < numRings; ++ringIdx) {
    reorderedBegins[ringIdx] = outputAtomsIdx + resAtomsFragBegin;
    const uint32_t sourceRingIdx = order[ringIdx];
    uint32_t begin = resBegins[resBeginsFragBegin + sourceRingIdx];
    const uint32_t end = resBegins[resBeginsFragBegin + sourceRingIdx + 1];
    for (; begin != end; ++begin, ++outputAtomsIdx) {
      reorderedAtoms[outputAtomsIdx] = resAtoms[begin];
    }
  }
  reorderedBegins[numRings] = resAtoms.size();
  std::copy(reorderedBegins.begin(), reorderedBegins.end(),
            resBegins.begin() + resBeginsFragBegin);
  std::copy(reorderedAtoms.begin(), reorderedAtoms.end(),
            resAtoms.begin() + resAtomsFragBegin);
}

void removeExtraRings(std::vector<uint32_t> &resBegins,
                      std::vector<atomindex_t> &resAtoms,
                      size_t resBeginsFragBegin,
                      size_t resAtomsFragBegin,
                      const RDMol &mol,
                      std::vector<atomindex_t> *extraRingAtoms,
                      std::vector<uint32_t> *extraRingBegins) {
  PRECONDITION(resAtomsFragBegin == resBegins[resBeginsFragBegin],
               "Mismatch between resBegins and resAtomsFragBegin");
  const uint32_t numRings = resBegins.size() - 1 - resBeginsFragBegin;

  // TODO: Move these vectors to the caller
  std::vector<uint32_t> tempBegins;
  std::vector<atomindex_t> tempAtoms;
  std::vector<uint32_t> tempOrder;
  sortRingsBySize(resBegins, resAtoms, resBeginsFragBegin, resAtomsFragBegin,
                  tempBegins, tempAtoms, tempOrder);

  // Use a vector of uint64_t, to avoid separately allocated bitsets
  std::vector<uint64_t> bitBrings;
  size_t bondBitElements = (mol.getNumBonds() + 63) / 64;
  bitBrings.resize(bondBitElements * numRings, 0);
  for (uint32_t ringIdx = 0, offset = 0; ringIdx < numRings;
       ++ringIdx, offset += bondBitElements) {
    const uint32_t begin = resBegins[ringIdx + resBeginsFragBegin];
    const uint32_t end = resBegins[ringIdx + resBeginsFragBegin + 1];
    for (uint32_t current = begin; current + 1 < end; ++current) {
      // change the rings from atom IDs to bondIds
      const uint32_t b = mol.getBondIndexBetweenAtoms(resAtoms[current],
                                                      resAtoms[current + 1]);
      bitBrings[offset + (b >> 6)] |= (uint64_t(1) << (b & 0x3F));
    }
    // Last bond of the ring, back to the beginning
    const uint32_t b =
        mol.getBondIndexBetweenAtoms(resAtoms[end - 1], resAtoms[begin]);
    bitBrings[offset + (b >> 6)] |= (uint64_t(1) << (b & 0x3F));
  }

  boost::dynamic_bitset<> availRings(numRings);
  availRings.set();
  boost::dynamic_bitset<> keepRings(numRings);
  std::vector<uint64_t> munionStorage;
  detail::BitSetWrapper munion(munionStorage, mol.getNumBonds(), false);
  uint64_t *const munionData = munion.data();

  // optimization - don't reallocate a new one each loop
  //boost::dynamic_bitset<> workspace(mol.getNumBonds());
  boost::dynamic_bitset<> consider(numRings);

  uint32_t numKeptRings = 0;

  for (unsigned int i = 0; i < numRings; ++i) {
    // skip this ring if we've already seen all of its bonds
    const uint64_t *bitBringsi = bitBrings.data() + (i * bondBitElements);
    bool isSubset = true;
    for (size_t idx = 0; idx < bondBitElements; ++idx) {
      isSubset &= ((bitBringsi[idx] & ~munionData[idx]) == 0);
    }
    if (isSubset) {
      availRings.set(i, 0);
    }
    if (!availRings[i]) {
      continue;
    }

    for (size_t idx = 0; idx < bondBitElements; ++idx) {
      munionData[idx] |= bitBringsi[idx];
    }
    keepRings.set(i);
    ++numKeptRings;

    // from this ring we consider all others that are still available and the
    // same size
    consider.reset();
    for (unsigned int j = i + 1; j < numRings; ++j) {
      // std::cerr<<"        "<<j<<" "<<(resBegins[j + 1] - resBegins[j])<<" -
      // "<<(resBegins[i + 1] - resBegins[i]).size()<<"  >"<<availRings[j]<<std::endl;
      if (availRings[j] && ((resBegins[resBeginsFragBegin + j + 1] -
                             resBegins[resBeginsFragBegin + j]) ==
                            (resBegins[resBeginsFragBegin + i + 1] -
                             resBegins[resBeginsFragBegin + i]))) {
        consider.set(j);
      }
    }

    while (consider.any()) {
      unsigned int bestJ = i + 1;
      int bestOverlap = -1;
      // loop over the available other rings in consideration and pick the one
      // that has the most overlapping bonds with what we've done so far.
      // this is the fix to github #526
      for (unsigned int j = i + 1; j < numRings; ++j) {
        if ((resBegins[resBeginsFragBegin + j + 1] -
             resBegins[resBeginsFragBegin + j]) !=
            (resBegins[resBeginsFragBegin + i + 1] -
             resBegins[resBeginsFragBegin + i])) {
          break;
        }
        if (!consider[j] || !availRings[j]) {
          continue;
        }
        // Count the overlap in bits
        const uint64_t *bitBringsj = bitBrings.data() + (j * bondBitElements);
        int overlap = 0;
        for (size_t idx = 0; idx < bondBitElements; ++idx) {
          overlap += std::popcount(bitBringsj[idx] & munionData[idx]);
        }
        if (overlap > bestOverlap) {
          bestOverlap = overlap;
          bestJ = j;
        }
      }
      consider.set(bestJ, 0);
      const uint64_t *bitBringsBestj =
          bitBrings.data() + (bestJ * bondBitElements);
      bool isSubset = true;
      for (size_t idx = 0; idx < bondBitElements; ++idx) {
        isSubset &= ((bitBringsBestj[idx] & ~munionData[idx]) == 0);
      }
      if (isSubset) {
        availRings.set(bestJ, 0);
      } else {
        keepRings.set(bestJ);
        ++numKeptRings;
        availRings.set(bestJ, 0);
        for (size_t idx = 0; idx < bondBitElements; ++idx) {
          munionData[idx] |= bitBringsBestj[idx];
        }
      }
    }
  }

  if (numKeptRings == numRings) {
    return;
  }
  const uint32_t numExtraRings = numRings - numKeptRings;

  // remove the extra rings from res and store them on the molecule in case we
  // wish symmetrize the SSSRs later
  uint32_t resDestRing = 0;
  uint32_t resDestAtom = resAtomsFragBegin;
  for (unsigned int i = 0; i < numRings; i++) {
    const bool keepRing = keepRings[i];
    const uint32_t resSourceAtomBegin = resBegins[resBeginsFragBegin + i];
    const uint32_t resSourceAtomEnd = resBegins[resBeginsFragBegin + i + 1];
    const uint32_t ringSize = resSourceAtomEnd - resSourceAtomBegin;
    if (!keepRing) {
      // add extra rings (there could already be some from previous fragments)
      if (extraRingAtoms != nullptr && extraRingBegins != nullptr) {
        if (extraRingBegins->size() == 0) {
          extraRingBegins->reserve(numExtraRings + 1);
          extraRingBegins->push_back(0);
        }
        extraRingAtoms->insert(extraRingAtoms->end(),
                               resAtoms.begin() + resSourceAtomBegin,
                               resAtoms.begin() + resSourceAtomEnd);
        extraRingBegins->push_back(extraRingAtoms->size());
      }
    } else {
      if (resDestRing != i) {
        std::copy(resAtoms.begin() + resSourceAtomBegin,
                  resAtoms.begin() + resSourceAtomEnd,
                  resAtoms.begin() + resDestAtom);
      }
      resDestAtom += ringSize;
      resBegins[resBeginsFragBegin + resDestRing + 1] = resDestAtom;
      ++resDestRing;
    }
  }
  resAtoms.resize(resDestAtom);
  resBegins.resize(resBeginsFragBegin + resDestRing + 1);
}

void findRingsD2nodes(const RDMol &tMol, UINT_VECT &resBegins,
                      UINT_VECT &resAtoms, RingUtils::RingInvarSet &invars,
                      const UINT_VECT &d2nodes,
                      INT_VECT &atomDegrees,
                      detail::BitSetWrapper &activeBonds,
                      detail::BitSetWrapper &ringBonds,
                      detail::BitSetWrapper &ringAtoms,
                      BFSWorkspace &bfs_workspace, UINT_VECT &sringBegins,
                      UINT_VECT &sringAtoms, UINT_VECT &nringBegins,
                      UINT_VECT &nringAtoms) {
  // place to record any duplicate rings discovered from the current d2 nodes
  struct FirstCand {
    atomindex_t cand;
    uint32_t count;
    uint32_t lastIdx;
  };
  std::vector<FirstCand> dupD2FirstCands;
  struct OtherCand {
    uint32_t ringIdx;
    atomindex_t cand;
    uint32_t prevIdx;
  };
  std::vector<OtherCand> dupD2OtherCands;
  uint32_t numRingsWithDups = 0;
  uint32_t maxDupPairs = 0;

  // TODO: These could be moved to the caller
  UINT_VECT changed;

  // here is an example of molecule where the this scheme of finding other node
  // that result in duplicates is necessary : C12=CON=C1C(C4)CC3CC2CC4C3
  // It would help to draw this molecule, and number the atoms but here is what
  // happen
  //  - there are 6 d2 node - 1, 6, 7, 9, 11, 13
  //  - both 6 and 7 find the same ring (5,6,12,13,8,7) but we do not find the 7
  //  membered ring
  //    (5,7,8,9,10,0,4)
  //  - similarly 9 and 11 find a duplicate ring (9,10,11,12,13)
  //  - when we move to 13 both the above duplicate rings are found
  //  - so we will keep track for each d2 all the other node that resulted in
  //  duplicate rings
  //  - the bonds to these nodes will be broken and we attempt to find a new
  //  ring, for e.g. by breaking
  //    bonds to 7 and 13, we will find a 7 membered ring with 6 (this is done
  //    in findSSSRforDupCands)
  for (auto &cand : d2nodes) {
    // std::cerr<<"    smallest rings bfs: "<<cand<<std::endl;
    sringBegins.clear();
    sringAtoms.clear();
    // we have to find all non duplicate possible smallest rings for each node
    bfs_workspace.smallestRingsBfs(tMol, cand, sringBegins, sringAtoms, activeBonds);
    for (uint32_t i = 0, n = sringBegins.size(); i < n; ++i) {
      const uint32_t begin = sringBegins[i];
      const uint32_t end = (i + 1 < n) ? sringBegins[i + 1] : sringAtoms.size();
      auto pair = invars.insert(sringAtoms, begin, end);
      if (pair.second) {
        resAtoms.insert(resAtoms.end(), sringAtoms.begin() + begin,
                        sringAtoms.begin() + end);
        resBegins.push_back(resAtoms.size());

        for (unsigned int i = begin; i < end - 1; ++i) {
          unsigned int bIdx =
              tMol.getBondIndexBetweenAtoms(sringAtoms[i], sringAtoms[i + 1]);
          ringBonds.set(bIdx);
          ringAtoms.set(sringAtoms[i]);
        }
        ringBonds.set(tMol.getBondIndexBetweenAtoms(sringAtoms[begin],
                                                    sringAtoms[end - 1]));
        ringAtoms.set(sringAtoms[end - 1]);
      }

      //nodeInvars[cand].push_back(invr);
      // check if this ring is duplicate with something else
      if (pair.first >= dupD2FirstCands.size()) {
        dupD2FirstCands.resize(
            std::max(pair.first + 1, 2 * dupD2FirstCands.size()),
            FirstCand{0, 0, 0});
      }
      auto &firstCandAndCount = dupD2FirstCands[pair.first];
      if (firstCandAndCount.count == 0) {
        firstCandAndCount.cand = cand;
      } else {
        numRingsWithDups += (firstCandAndCount.count == 1);
        maxDupPairs += firstCandAndCount.count;
        dupD2OtherCands.push_back(OtherCand{uint32_t(pair.first), cand, firstCandAndCount.lastIdx});
        firstCandAndCount.lastIdx = dupD2OtherCands.size() - 1;
      }
      ++firstCandAndCount.count;
    }

    // We don't want to trim the bonds connecting cand here - this can disrupt
    // a second small ring. Here is an example SC(C3C1CC(C3)CC(C2S)(O)C1)2S
    // by trimming the bond connecting to atom #4, we lose the smallest ring
    // that contains atom #7. Issue 134

    // But if there were no rings found, trimming isn't dangerous, and can
    // save wasted time for long chains.
    if (sringBegins.empty()) {
      changed.clear();
      changed.push_back(cand);
      do {
        uint32_t cand = changed[0];
        std::pop_heap(changed.begin(), changed.end());
        changed.resize(changed.size() - 1);
        trimBonds(cand, tMol, &changed, atomDegrees, activeBonds);
      } while (!changed.empty());
    }
  }

  if (numRingsWithDups == 0) {
    return;
  }
  // Reorganize data about duplicates
  const uint32_t totalDups = dupD2OtherCands.size() + numRingsWithDups;
  std::vector<atomindex_t> dupD2Cands;
  dupD2Cands.reserve(totalDups);
  std::vector<uint32_t> dup2CandBegins;
  dup2CandBegins.reserve(numRingsWithDups + 1);
  std::vector<std::pair<uint32_t, uint32_t>> order;
  if (numRingsWithDups > 1) {
    order.reserve(numRingsWithDups);
  }
  std::vector<std::pair<uint32_t, uint32_t>> dupMap;
  dupMap.reserve(maxDupPairs);
  for (uint32_t ringIdx = 0, numRings = resBegins.size() - 1; ringIdx < numRings;
       ++ringIdx) {
    auto &firstCand = dupD2FirstCands[ringIdx];
    if (firstCand.count < 2) {
      continue;
    }
    if (numRingsWithDups > 1) {
      order.push_back(std::make_pair(ringIdx, dup2CandBegins.size()));
    }
    dup2CandBegins.push_back(dupD2Cands.size());

    const size_t beginIdx = dupD2Cands.size();
    dupD2Cands.insert(dupD2Cands.end(), firstCand.count, firstCand.cand);
    size_t idx = dupD2Cands.size() - 1;
    size_t otherCandsIdx = firstCand.lastIdx;
    for (; idx > beginIdx; --idx) {
      dupD2Cands[idx] = dupD2OtherCands[otherCandsIdx].cand;
      otherCandsIdx = dupD2OtherCands[otherCandsIdx].prevIdx;
    }

    // ok we discovered this ring via another node before
    // add that node as duplicate to this node and vice versa
    for (size_t i = beginIdx + 1, end = dupD2Cands.size(); i < end; ++i) {
      for (size_t j = beginIdx; j < i; ++j) {
        uint32_t a = dupD2Cands[i];
        uint32_t b = dupD2Cands[j];
        if (a != b) {
          dupMap.push_back(std::make_pair(a, b));
          dupMap.push_back(std::make_pair(b, a));
        }
      }
    }
  }
  dup2CandBegins.push_back(dupD2Cands.size());

  // If this sort becomes a bottleneck, it's possible to count, preallocate
  // dupMap, and write to dupMap in a separate pass. That would also avoid
  // the need for std::lower_bound and std::upper_bound in findSSSRforDupCands
  std::sort(dupMap.begin(), dupMap.end());

  // Find the iteration order to match earlier behaviour, for compatibility
  std::sort(order.begin(), order.end(),
            [&invars](const std::pair<uint32_t, uint32_t> &a,
                      const std::pair<uint32_t, uint32_t> &b) -> bool {
    return invars.less(a.first, b.first);
  });

  // now deal with any d2 nodes that resulted in duplicate rings before trimming
  // their bonds.
  // it is possible that one of these nodes is involved a different small ring,
  // that is not found  because the first nodes has not be trimmed. Here is an
  // example molecule:
  // CC1=CC=C(C=C1)S(=O)(=O)O[CH]2[CH]3CO[CH](O3)[CH]4OC(C)(C)O[CH]24
  findSSSRforDupCands(tMol, resBegins, resAtoms, invars, dupMap, dup2CandBegins,
                      dupD2Cands, order, atomDegrees, activeBonds, bfs_workspace,
                      sringBegins, sringAtoms, nringBegins, nringAtoms);
}

void findRingsD3Node(const RDMol &tMol, UINT_VECT &resBegins,
                     UINT_VECT &resAtoms, RingUtils::RingInvarSet &invars,
                     int cand, const detail::BitSetWrapper &activeBonds,
                     BFSWorkspace &bfs_workspace, UINT_VECT &sringBegins,
                     UINT_VECT &sringAtoms, UINT_VECT &tringBegins,
                     UINT_VECT &tringAtoms) {
  // this is brutal - we have no degree 2 nodes - find the first possible degree
  // 3 node

  // We've got a degree three node. The goal of what follows is to find the
  // three rings in which it's involved, push those onto our results, and
  // then remove the node from consideration.  This will create a bunch of
  // degree
  // 2 nodes, which we can then chew off the next time around the loop.

  // this part is a bit different from the Figueras algorithm
  // here we try to find all the rings the rings that have a potential for
  // contributing to
  // SSSR - i.e. we try to find 3 rings for this node.
  // - each bond (that contributes to the degree 3 ) is allowed to participate
  // in exactly
  //    two of these rings.
  // - also any rings that are included in already found rings are ignored

  // ASSUME: every connection from a degree three node at this point is a
  //         ring bond
  // REVIEW: Is this valid?

  // TODO: Move these to the caller
  UINT_VECT forb;

  // first find all smallest possible rings
  sringBegins.clear();
  sringAtoms.clear();
  bfs_workspace.smallestRingsBfs(tMol, cand, sringBegins, sringAtoms,
                                 activeBonds);
  const uint32_t nsmall = sringBegins.size();

  for (uint32_t i = 0, n = nsmall; i < n; ++i) {
    const uint32_t begin = sringBegins[i];
    const uint32_t end = (i + 1 < n) ? sringBegins[i + 1] : sringAtoms.size();
    if (invars.insert(sringAtoms, begin, end).second) {
      resAtoms.insert(resAtoms.end(), sringAtoms.begin() + begin,
                      sringAtoms.begin() + end);
      resBegins.push_back(resAtoms.size());
    }
  }

  // if already found >3 rings we are done with this degree 3 node
  // if we found less than 3 we have to find other potential ring/s
  if (nsmall < 3) {
    auto [beg, end] = tMol.getAtomBonds(cand);
    auto [begAtoms, endAtoms] = tMol.getAtomNeighbors(cand);
    while (beg != end && !activeBonds[*beg]) {
      ++beg;
      ++begAtoms;
    }
    CHECK_INVARIANT(beg != end, "neighbor not found");
    int n1 = *begAtoms;

    ++beg;
    ++begAtoms;
    while (beg != end && !activeBonds[*beg]) {
      ++beg;
      ++begAtoms;
    }
    CHECK_INVARIANT(beg != end, "neighbor not found");
    int n2 = *begAtoms;

    ++beg;
    ++begAtoms;
    while (beg != end && !activeBonds[*beg]) {
      ++beg;
      ++begAtoms;
    }
    CHECK_INVARIANT(beg != end, "neighbor not found");
    int n3 = *begAtoms;

    if (nsmall == 2) {
      // we found two rings find the third one
      // first find the neighbor that is common to the two ring we found so far
      int f = -1;

      const uint32_t sringBegin0 = sringBegins[0];
      const uint32_t sringBegin1 = sringBegins[1];
      const uint32_t sringEnd0 = sringBegin1;
      const uint32_t sringEnd1 = sringAtoms.size();

      const auto begin0 = sringAtoms.begin() + sringBegin0;
      const auto begin1 = sringAtoms.begin() + sringBegin1;
      const auto end0 = sringAtoms.begin() + sringEnd0;
      const auto end1 = sringAtoms.begin() + sringEnd1;

      if ((std::find(begin0, end0, n1) != end0) &&
          (std::find(begin1, end1, n1) != end1)) {
        f = n1;
      } else if ((std::find(begin0, end0, n2) != end0) &&
                 (std::find(begin1, end1, n2) != end1)) {
        f = n2;
      } else if ((std::find(begin0, end0, n3) != end0) &&
                 (std::find(begin1, end1, n3) != end1)) {
        f = n3;
      }
      CHECK_INVARIANT(f >= 0, "third ring not found");

      // now find the smallest possible ring that does not contain f
      tringBegins.clear();
      tringAtoms.clear();
      forb.clear();
      forb.push_back(f);
      bfs_workspace.smallestRingsBfs(tMol, cand, tringBegins, tringAtoms,
                                     activeBonds, &forb);
      for (uint32_t i = 0, n = tringBegins.size(); i < n; ++i) {
        const uint32_t begin = tringBegins[i];
        const uint32_t end =
            (i + 1 < n) ? tringBegins[i + 1] : tringAtoms.size();

        if (invars.insert(tringAtoms, begin, end).second) {
          resAtoms.insert(resAtoms.end(), tringAtoms.begin() + begin,
                          tringAtoms.begin() + end);
          resBegins.push_back(resAtoms.size());
        }
      }
    }  // doing degree 3 node  - end of 2 smallest rings found for cand
    else if (nsmall == 1) {
      // we found 1 ring - we need to find two more that involve the 3rd
      // neighbor
      int f1 = -1, f2 = -1;

      const uint32_t sringBegin0 = sringBegins[0];
      const uint32_t sringEnd0 = sringAtoms.size();

      const auto begin0 = sringAtoms.begin() + sringBegin0;
      const auto end0 = sringAtoms.begin() + sringEnd0;

      // Which of our three neighbors are in the small ring?
      //   these are f1 and f2
      if (std::find(begin0, end0, n1) == end0) {
        f1 = n2, f2 = n3;
      } else if (std::find(begin0, end0, n2) == end0) {
        f1 = n1;
        f2 = n3;
      } else if (std::find(begin0, end0, n3) == end0) {
        f1 = n1;
        f2 = n2;
      }
      CHECK_INVARIANT(f1 >= 0, "rings not found");
      CHECK_INVARIANT(f2 >= 0, "rings not found");

      // now find two rings that include cand, one of these rings should include
      // f1
      // and the other should include f2

      // first ring with f1 and no f2
      tringBegins.clear();
      tringAtoms.clear();
      forb.clear();
      forb.push_back(f2);
      bfs_workspace.smallestRingsBfs(tMol, cand, tringBegins, tringAtoms,
                                     activeBonds, &forb);
      for (uint32_t i = 0, n = tringBegins.size(); i < n; ++i) {
        const uint32_t begin = tringBegins[i];
        const uint32_t end =
            (i + 1 < n) ? tringBegins[i + 1] : tringAtoms.size();

        if (invars.insert(tringAtoms, begin, end).second) {
          resAtoms.insert(resAtoms.end(), tringAtoms.begin() + begin,
                          tringAtoms.begin() + end);
          resBegins.push_back(resAtoms.size());
        }
      }

      // next the ring with f2 and no f1
      tringBegins.clear();
      tringAtoms.clear();
      forb.clear();
      forb.push_back(f1);
      bfs_workspace.smallestRingsBfs(tMol, cand, tringBegins, tringAtoms,
                                     activeBonds, &forb);
      for (uint32_t i = 0, n = tringBegins.size(); i < n; ++i) {
        const uint32_t begin = tringBegins[i];
        const uint32_t end =
            (i + 1 < n) ? tringBegins[i + 1] : tringAtoms.size();

        if (invars.insert(tringAtoms, begin, end).second) {
          resAtoms.insert(resAtoms.end(), tringAtoms.begin() + begin,
                          tringAtoms.begin() + end);
          resBegins.push_back(resAtoms.size());
        }
      }
    }  // doing node of degree 3 - end of found only 1 smallest ring
  }  // end of found less than 3 smallest ring for the degree 3 node
}

/******************************************************************************
 * SUMMARY:
 *  remove the bond in the molecule that connect to the specified atom
 *
 * ARGUMENTS:
 *  cand - the node(atom) of interest
 *  tMol - molecule of interest
 *  changed - list of the atoms that are effected the bond removal
 *             this may be accumulated over multiple calls to trimBonds
 *             it basically forms a list of atom that need to be searched for
 *             the next round of pruning
 *
 ******************************************************************************/
void trimBonds(const uint32_t cand, const RDMol &tMol, std::vector<uint32_t> *changed,
               INT_VECT &atomDegrees, detail::BitSetWrapper &activeBonds) {
  auto [beg, end] = tMol.getAtomBonds(cand);
  for (; beg != end; ++beg) {
    const uint32_t bondIndex = *beg;
    if (!activeBonds.get(bondIndex)) {
      continue;
    }
    const BondData &bond = tMol.getBond(bondIndex);
    const uint32_t oIdx = bond.getOtherAtomIdx(cand);
    if (changed != nullptr && atomDegrees[oIdx] <= 2) {
      changed->push_back(oIdx);
      std::push_heap(changed->begin(), changed->end(), std::greater());
    }
    activeBonds.reset(bondIndex);
    atomDegrees[oIdx] -= 1;
    atomDegrees[cand] -= 1;
  }
}

/*******************************************************************************
 * SUMMARY:
 *  this again is a modified version of the BFS algorithm in Figueras paper to
 *  find the smallest ring with a specified root atom.
 *    JCICS, Vol. 30, No. 5, 1996, 986-991
 *  The following are changes from the original algorithm
 *   - find all smallest rings around a node not just one
 *   - once can provided a list of node IDs that should not be include in the
 *     discovered rings
 *
 * ARGUMENTS:
 *  mol - molecule of interest
 *  root - Atom ID of the node of interest
 *  ringBegins & ringAtoms - list of rings into which the results are entered
 *                           Unlike in other situations, ringBegins will only
 *                           have as many elements as rings, with no extra
 *                           final index for the end.
 *  forbidden - list of atoms ID that should be avoided
 *
 * RETURNS:
 *  number of smallest rings found
 ***********************************************************************************/
void BFSWorkspace::smallestRingsBfs(const RDMol &mol, atomindex_t root,
                                        UINT_VECT &ringBegins,
                                        UINT_VECT &ringAtoms,
                                        const detail::BitSetWrapper &activeBonds,
                                        UINT_VECT *forbidden) {
  // this function finds the smallest ring with the given root atom.
  // if multiple smallest rings are found all of them are returned
  // if any atoms are specified in the forbidden list, those atoms are avoided.

  // FIX: this should be number of atoms in the fragment (if it's required at
  // all, see below)
  const int WHITE = 0, GRAY = 1, BLACK = 2;
  d_done.assign(mol.getNumAtoms(), WHITE);

  if (forbidden) {
    for (auto i : *forbidden) {
      d_done[i] = BLACK;
    }
  }

  d_parents.assign(mol.getNumAtoms(), -1);
  d_depths.assign(mol.getNumAtoms(), 0);

  // Don't bother overwriting the queue contents. Just treat it as empty
  if (d_bfsq.size() == 0) {
    // Start with a reasonable initial size to reduce repeated reallocation
    d_bfsq.resize(std::min(8u, mol.getNumAtoms()));
  }
  d_bfsq[0] = root;
  uint32_t queueBegin = 0;
  uint32_t queueSize = 1;

  unsigned int curSize = UINT_MAX;
  while (queueSize > 0) {
    if (queueSize >= RingUtils::MAX_BFSQ_SIZE) {
      std::string msg =
          "Maximum BFS search size exceeded.\nThis is likely due to a highly "
          "symmetric fused ring system.";
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw ValueErrorException(msg);
    }

    const int curr = d_bfsq[queueBegin];
    ++queueBegin;
    --queueSize;
    if (queueBegin == d_bfsq.size() || queueSize == 0) {
      queueBegin = 0;
    }
    d_done[curr] = BLACK;

    const unsigned int depth = d_depths[curr] + 1;
    if (depth > curSize) {
      // depth is the shortest cycle I _could_ find this round.
      break;
    }

    auto [beg, end] = mol.getAtomBonds(curr);
    auto [otherAtoms, otherAtomsEnd] = mol.getAtomNeighbors(curr);
    while (beg != end) {
      const uint32_t bondIndex = *beg;
      const uint32_t nbrIdx = *otherAtoms;
      ++beg;
      ++otherAtoms;
      if (!activeBonds[bondIndex]) {
        continue;
      }
      if (d_done[nbrIdx] == BLACK || d_parents[curr] == int(nbrIdx)) {
        continue;
      }
      if (d_done[nbrIdx] == WHITE) {
        // we have never been to this node before through via any path
        d_parents[nbrIdx] = curr;
        d_done[nbrIdx] = GRAY;
        d_depths[nbrIdx] = depth;
        if (queueSize == d_bfsq.size()) {
          // Reallocate with more space
          std::vector<atomindex_t> newQueue;
          const size_t newCapacity = 2 * queueSize;
          newQueue.reserve(newCapacity);

          // Copy all data over. Remember the queue is full.
            newQueue.insert(newQueue.end(), d_bfsq.begin() + queueBegin,
                            d_bfsq.end());
          if (queueBegin != 0) {
            newQueue.insert(newQueue.end(), d_bfsq.begin(),
                            d_bfsq.begin() + queueBegin);
          }
          // The queue is only queueSize in size, but the underlying vector
          // needs to have its size set to the queue capacity
          newQueue.resize(newCapacity);
          d_bfsq.swap(newQueue);
          queueBegin = 0;
        }
        uint32_t queueNext = queueBegin + queueSize;
        if (queueNext >= d_bfsq.size()) {
          queueNext -= d_bfsq.size();
        }
        d_bfsq[queueNext] = nbrIdx;
        ++queueSize;
      } else {
        // we have been here via a different path
        // there is a potential for ring closure here
        // stitch together the two paths

        d_ring.clear();
        d_ring.push_back(nbrIdx);
        // forwards path
        int parent = d_parents[nbrIdx];
        while (parent != -1 && parent != int(root)) {
          d_ring.push_back(parent);
          parent = d_parents[parent];
        }

        // backwards path
        const size_t backwardBegin = d_ring.size();
        d_ring.push_back(curr);
        parent = d_parents[curr];
        while (parent != -1) {
          // Is the least common ancestor not the root?
          // There shouldn't be any duplicates within a single path
          if (std::find(d_ring.begin(), d_ring.begin() + backwardBegin,
                        parent) != d_ring.begin() + backwardBegin) {
            d_ring.clear();
            break;
          }
          d_ring.push_back(parent);
          parent = d_parents[parent];
        }

        // Found a new small ring including the root.
        if (d_ring.size() > 1) {
          if (d_ring.size() <= curSize) {
            curSize = rdcast<unsigned int>(d_ring.size());
            ringBegins.push_back(ringAtoms.size());
            if (ringAtoms.size() + d_ring.size() < ringAtoms.capacity()) {
              ringAtoms.reserve(std::max(d_ring.size(), 2 * ringAtoms.size()));
            }
            for (size_t i = d_ring.size() - 1; i >= backwardBegin; --i) {
              ringAtoms.push_back(d_ring[i]);
            }
            ringAtoms.insert(ringAtoms.end(), d_ring.begin(),
                             d_ring.begin() + backwardBegin);

          } else {
            // we are done with the smallest rings
            return;
          }
        }
      }
    }  // end of loop over neighbors of current atom
  }  // moving to the next node

  // if we are here we should have found everything around the node
}

bool _atomSearchBFS(const RDMol &tMol, uint32_t startAtomIdx,
                    uint32_t endAtomIdx, const detail::BitSetWrapper &ringAtoms,
                    UINT_VECT &res, RingUtils::RingInvarSet &invars) {
  res.clear();

  // TODO: Move these to the caller
  struct LinkPair {
    size_t previousIdx;
    uint32_t currentValue;
  };
  std::vector<LinkPair> listQueue;
  size_t queueFront = 0;
  //std::deque<UINT_VECT> bfsq;

  //UINT_VECT tv;
  //tv.push_back(startAtomIdx);
  //bfsq.push_back(tv);
  // Linking to the current index isn't valid, so this represents the end of the lists.
  listQueue.push_back(LinkPair{size_t(0), startAtomIdx});
  while (queueFront < listQueue.size()) {
    if (listQueue.size() - queueFront >= RingUtils::MAX_BFSQ_SIZE) {
      std::string msg =
          "Maximum BFS search size exceeded.\nThis is likely due to a highly "
          "symmetric fused ring system.";
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw ValueErrorException(msg);
    }
    //tv = bfsq.front();
    //bfsq.pop_front();
    const size_t currQueueIdx = queueFront;
    ++queueFront;

    unsigned int currAtomIdx = listQueue[currQueueIdx].currentValue;
    auto [nbrIdx, endNbrs] = tMol.getAtomNeighbors(currAtomIdx);
    while (nbrIdx != endNbrs) {
      if (*nbrIdx == endAtomIdx) {
        if (currAtomIdx != startAtomIdx) {
          //UINT_VECT nv(tv);
          //nv.push_back(rdcast<unsigned int>(*nbrIdx));

          // make sure the ring we just found isn't already in our set
          // of rings (this was an extension of sf.net issue 249)

          // Add nbrIdx and all previous atoms
          const size_t oldInvarsSize = invars.tempAddOne();
          uint64_t *newBits = invars.bits.data() + oldInvarsSize;
          newBits[(*nbrIdx) >> 6] |= (uint64_t(1) << ((*nbrIdx) & 0x3F));
          newBits[currAtomIdx >> 6] |= (uint64_t(1) << (currAtomIdx & 0x3F));
          size_t listSize = 2;
          size_t iterQueueIdx = currQueueIdx;
          while (iterQueueIdx != 0) {
            iterQueueIdx = listQueue[iterQueueIdx].previousIdx;
            uint32_t value = listQueue[iterQueueIdx].currentValue;
            newBits[value >> 6] |= (uint64_t(1) << (value & 0x3F));
            ++listSize;
          }
          bool found = invars.find(newBits, oldInvarsSize).second;
          invars.bits.resize(oldInvarsSize);
          if (found) {
            // we're done!
            res.resize(listSize);
            // Fill them in the opposite order
            res[listSize - 1] = *nbrIdx;
            res[listSize - 2] = currAtomIdx;
            listSize -= 2;
            size_t iterQueueIdx = currQueueIdx;
            while (iterQueueIdx != 0) {
              iterQueueIdx = listQueue[iterQueueIdx].previousIdx;
              --listSize;
              res[listSize] = listQueue[iterQueueIdx].currentValue;
            }
            return true;
          }
        } else {
          // ignore this one
        }
      } else if (ringAtoms[*nbrIdx]) {
        // If nbrIdx isn't already in this path, add it
        bool found = (currAtomIdx == *nbrIdx);
        size_t iterQueueIdx = currQueueIdx;
        while (!found && iterQueueIdx != 0) {
          iterQueueIdx = listQueue[iterQueueIdx].previousIdx;
          found = (listQueue[iterQueueIdx].currentValue == *nbrIdx);
        }

        if (!found) {
          //} else if(ringAtoms[*nbrIdx]){
          // UINT_VECT nv(tv);
          // nv.push_back(rdcast<unsigned int>(*nbrIdx));
          listQueue.push_back(LinkPair{currQueueIdx, *nbrIdx});
          // bfsq.push_back(nv);
        }
      } 

      ++nbrIdx;
    }
  }
  return false;
}

bool findRingConnectingAtoms(const RDMol &tMol, const uint32_t bondIndex,
                             UINT_VECT &resBegins, UINT_VECT &resAtoms,
                             RingUtils::RingInvarSet &invars,
                             detail::BitSetWrapper &ringBonds,
                             detail::BitSetWrapper &ringAtoms) {
  PRECONDITION(bondIndex < tMol.getNumBonds(), "bad bond");
  PRECONDITION(!ringBonds[bondIndex], "not a ring bond");
  const BondData &bond = tMol.getBond(bondIndex);
  PRECONDITION(ringAtoms[bond.getBeginAtomIdx()], "not a ring atom");
  PRECONDITION(ringAtoms[bond.getEndAtomIdx()], "not a ring atom");

  // TODO: Move this to the caller
  UINT_VECT nring;
  if (_atomSearchBFS(tMol, bond.getBeginAtomIdx(), bond.getEndAtomIdx(),
                     ringAtoms, nring, invars)) {
    // TODO: This already gets checked inside _atomSearchBFS, so should always
    // insert, and it might as well be inserted there.
    if (invars.insert(nring, 0, nring.size()).second) {
      resAtoms.insert(resAtoms.end(), nring.begin(), nring.end());
      resBegins.push_back(resAtoms.size());
      for (unsigned int i = 0; i < nring.size() - 1; ++i) {
        unsigned int bIdx =
            tMol.getBondIndexBetweenAtoms(nring[i], nring[i + 1]);
        ringBonds.set(bIdx);
        ringAtoms.set(nring[i]);
      }
      ringBonds.set(
          tMol.getBondIndexBetweenAtoms(nring[0], nring[nring.size() - 1]));
      ringAtoms.set(nring[nring.size() - 1]);
    }
  } else {
    return false;
  }
  return true;
}

}  // namespace FindRings

namespace RDKit {
namespace MolOps {
int findSSSR(const ROMol &mol, VECT_INT_VECT *res, bool includeDativeBonds,
             bool includeHydrogenBonds) {
  if (!res) {
    VECT_INT_VECT rings;
    return findSSSR(mol, rings, includeDativeBonds, includeHydrogenBonds);
  } else {
    return findSSSR(mol, (*res), includeDativeBonds, includeHydrogenBonds);
  }
}

int findSSSR(const ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
             bool includeHydrogenBonds) {
  // The original code wrote to the RingInfo, even though the mol is const.
  std::vector<atomindex_t> extraRingAtoms;
  std::vector<uint32_t> extraRingBegins;
  // NOTE: The const_cast must be before calling getRingInfo(), else the
  // compatibility data sync status for the ring info won't be set correctly.
  int retValue = findSSSR(
      mol.asRDMol(), const_cast<ROMol &>(mol).asRDMol().getRingInfo(),
      includeDativeBonds, includeHydrogenBonds, &extraRingAtoms,
      &extraRingBegins);

  // Fill in res from the ring info
  const RingInfoCache &ringInfo = mol.asRDMol().getRingInfo();
  const size_t numRings = ringInfo.numRings();
  const uint32_t *data = ringInfo.atomsInRings.data();
  // Ensure the reinterpret_cast is safe, at least for atom indices
  static_assert(sizeof(int) == sizeof(uint32_t));
  const int *intData = reinterpret_cast<const int *>(data);
  res.clear();
  res.reserve(numRings);
  for (size_t ringIndex = 0; ringIndex < numRings; ++ringIndex) {
    uint32_t beginIndex = ringInfo.ringBegins[ringIndex];
    uint32_t endIndex = ringInfo.ringBegins[ringIndex + 1];
    res.emplace_back(intData + beginIndex, intData + endIndex);
  }

  if (extraRingAtoms.size() != 0) {
    // FIXME: Add extraRings property
  }
  return retValue;
}
int findSSSR(const RDMol &mol, RingInfoCache &ringInfo, bool includeDativeBonds,
             bool includeHydrogenBonds,
             std::vector<atomindex_t> *extraRingAtoms,
             std::vector<uint32_t> *extraRingBegins) {
  if (ringInfo.isInitialized()) {
    ringInfo.reset();
  }
  ringInfo.isInit = true;
  ringInfo.findRingType = FIND_RING_TYPE_SSSR;

  auto &resAtoms = ringInfo.atomsInRings;
  auto &resBegins = ringInfo.ringBegins;
  if (resBegins.size() == 0) {
    resBegins.push_back(0);
  }

  unsigned int nats = mol.getNumAtoms();
  //boost::dynamic_bitset<> activeAtoms(nats);
  //activeAtoms.set();
  uint32_t nbnds = mol.getNumBonds();
  std::vector<uint64_t> activeBondsStorage;
  detail::BitSetWrapper activeBonds(activeBondsStorage, nbnds, true);

  // Zero-order bonds are not candidates for rings, and dative bonds and
  // hydrogen bonds may also be out
  for (uint32_t bond_idx = 0; bond_idx < nbnds; ++bond_idx) {
    const BondEnums::BondType bondType = mol.getBond(bond_idx).getBondType();
    if (bondType == BondEnums::BondType::ZERO ||
        (!includeDativeBonds && isDative(bondType)) ||
        (!includeHydrogenBonds && bondType == BondEnums::BondType::HYDROGEN)) {
      activeBonds.reset(bond_idx);
    }
  }

  std::vector<uint64_t> ringBondsStorage;
  detail::BitSetWrapper ringBonds(ringBondsStorage, mol.getNumBonds());
  std::vector<uint64_t> ringAtomsStorage;
  detail::BitSetWrapper ringAtoms(ringAtomsStorage, mol.getNumAtoms());

  INT_VECT atomDegrees(nats);
  INT_VECT atomDegreesWithZeroOrderBonds(nats);
  for (unsigned int i = 0; i < nats; ++i) {
    int deg = mol.getAtomDegree(i);
    atomDegrees[i] = deg;
    atomDegreesWithZeroOrderBonds[i] = deg;
    for (auto [begin, end] = mol.getAtomBonds(i); begin != end; ++begin) {
      const BondEnums::BondType bondType = mol.getBond(*begin).getBondType();
      if (bondType == BondEnums::BondType::ZERO ||
          (!includeDativeBonds && isDative(bondType)) ||
          (!includeHydrogenBonds && bondType == BondEnums::BondType::HYDROGEN)) {
        atomDegrees[i]--;
      }
    }
  }
  if (extraRingAtoms != nullptr && extraRingBegins != nullptr) {
    extraRingAtoms->resize(0);
    extraRingBegins->resize(0);
  }

  // find the number of fragments in the molecule - we will loop over them
  std::vector<uint32_t> atomFrags;
  unsigned int nfrags = getMolFrags(mol, atomFrags);

  // Invert to find the atoms in each fragment
  // Start by counting the number of atoms in each fragment
  std::vector<uint32_t> fragStarts;
  fragStarts.resize(nfrags+1, 0);
  for (uint32_t atomFrag : atomFrags) {
    ++fragStarts[atomFrag+1];
  }
  // Find where each fragment starts
  std::partial_sum(fragStarts.begin() + 1, fragStarts.end(),
                   fragStarts.begin() + 1);

  // Fill in the atoms in each fragment
  std::vector<atomindex_t> frags;
  frags.resize(nats);
  for (atomindex_t atom_idx = 0; atom_idx < nats; ++atom_idx) {
    const uint32_t frag = atomFrags[atom_idx];
    const uint32_t fragNext = fragStarts[frag];
    ++fragStarts[frag];
    frags[fragNext] = atom_idx;
  }
  // Update starts
  std::copy_backward(fragStarts.begin(), fragStarts.end() - 1, fragStarts.end());
  fragStarts[0] = 0;

  // the following is the list of atoms that are useful in the next round of
  // trimming basically atoms that become degree 0 or 1 because of bond
  // removals initialized with atoms of degrees 0 and 1
  // It was previously an ordered set, but only the first entry was visited
  // at any time, so it's now a min heap.
  std::vector<atomindex_t> changed;

  // TODO: RINGINVAR_SET is never removed from, so replace it with a vector of
  // bit set data and a set of integers indexing into it. To compare and insert
  // if not present, add to the end of the vector, pass that new index into
  // insert and check return. If not inserted, resize the vector back down.
  //RINGINVAR_SET invars;
  RingUtils::RingInvarSet invars(nats);

  // Outside the loop to avoid repeated reallocation
  std::vector<uint64_t> doneAtsStorage;
  detail::BitSetWrapper doneAts(doneAtsStorage, nats);
  std::vector<uint64_t> forbStorage;
  detail::BitSetWrapper forb(forbStorage, nats);
  std::vector<uint32_t> possibleBonds;
  std::vector<uint64_t> deadBondsStorage;
  detail::BitSetWrapper deadBonds(deadBondsStorage, mol.getNumBonds());

  UINT_VECT d2nodes;
  FindRings::BFSWorkspace bfs_workspace;
  UINT_VECT sringBegins;
  UINT_VECT sringAtoms;
  UINT_VECT tringBegins;
  UINT_VECT tringAtoms;

  // loop over the fragments in a molecule
  for (unsigned int fi = 0; fi < nfrags; ++fi) {
    // TODO: Is this necessary, or is changed always empty at this point?
    changed.clear();

    const size_t resAtomsFragBegin = resAtoms.size();
    const size_t resBeginsFragBegin = resBegins.size() - 1;

    const size_t curFragSize = fragStarts[fi + 1] - fragStarts[fi];
    const atomindex_t *const curFragBegin = frags.data() + fragStarts[fi];
    const atomindex_t *const curFragEnd = frags.data() + fragStarts[fi + 1];

    if (curFragSize < 3) {
      continue;
    }

    int bndcnt_with_zero_order_bonds = 0;
    unsigned int nbnds = 0;
    for (auto *current = curFragBegin; current != curFragEnd; ++current) {
      auto atom_idx = *current;
      bndcnt_with_zero_order_bonds += atomDegreesWithZeroOrderBonds[atom_idx];

      int deg = atomDegrees[atom_idx];

      nbnds += deg;
      if (deg < 2) {
        changed.push_back(atom_idx);
      }
    }
    std::make_heap(changed.begin(), changed.end(), std::greater());

    // check to see if this fragment can even have a possible ring
    CHECK_INVARIANT(bndcnt_with_zero_order_bonds % 2 == 0,
                    "fragment graph has a dangling degree");
    bndcnt_with_zero_order_bonds = bndcnt_with_zero_order_bonds / 2;
    int num_possible_rings = bndcnt_with_zero_order_bonds - curFragSize + 1;
    if (num_possible_rings < 1) {
      continue;
    }

    CHECK_INVARIANT(nbnds % 2 == 0,
                    "fragment graph problem when including zero-order bonds");
    nbnds = nbnds / 2;

    doneAts.reset();
    unsigned int nAtomsDone = 0;
    while (nAtomsDone < curFragSize) {
      // std::cerr<<" ndone: "<<nAtomsDone<<std::endl;
      // std::cerr<<" activeBonds: "<<activeBonds<<std::endl;
      // std::cerr<<"  done: ";
      // trim all bonds that connect to degree 0 and 1 atoms
      while (changed.size() > 0) {
        int cand = changed[0];
        std::pop_heap(changed.begin(), changed.end());
        changed.resize(changed.size()-1);
        if (!doneAts[cand]) {
          // std::cerr<<cand<<" ";
          doneAts.set(cand);
          ++nAtomsDone;
          FindRings::trimBonds(cand, mol, &changed, atomDegrees, activeBonds);
        }
      }
      // std::cerr<<std::endl;
      // std::cerr<<"activeBonds2: "<<activeBonds<<std::endl;

      // all atoms left in the fragment should at least have a degree >= 2
      // collect all the degree two nodes;
      FindRings::pickD2Nodes(mol, d2nodes, curFragBegin, curFragEnd,
                             atomDegrees, activeBonds, forb);
      if (d2nodes.size() > 0) {  // deal with the current degree two nodes
        // place to record any duplicate rings discovered from the current d2
        // nodes
        FindRings::findRingsD2nodes(mol, resBegins, resAtoms, invars, d2nodes,
                                    atomDegrees, activeBonds, ringBonds,
                                    ringAtoms, bfs_workspace, sringBegins,
                                    sringAtoms, tringBegins, tringAtoms);
        // trim after we have dealt with all the current d2 nodes,
        for (auto d2i = d2nodes.begin(); d2i != d2nodes.end(); d2i++) {
          doneAts.set(*d2i);
          ++nAtomsDone;
          FindRings::trimBonds((*d2i), mol, &changed, atomDegrees, activeBonds);
        }
      }  // end of degree two nodes
      else if (nAtomsDone < curFragSize) {  // now deal with higher degree nodes
        // this is brutal - we have no degree 2 nodes - find the first
        // possible degree 3 node
        int cand = -1;
        for (auto aidi = curFragBegin; aidi != curFragEnd;
             aidi++) {
          unsigned int deg = atomDegrees[*aidi];
          if (deg == 3) {
            cand = (*aidi);
            break;
          }
        }

        // if we did not find a degree 3 node we are done
        // REVIEW:
        if (cand == -1) {
          break;
        }
        FindRings::findRingsD3Node(mol, resBegins, resAtoms, invars, cand,
                                   activeBonds, bfs_workspace, sringBegins,
                                   sringAtoms, tringBegins, tringAtoms);
        doneAts.set(cand);
        ++nAtomsDone;
        FindRings::trimBonds(cand, mol, &changed, atomDegrees, activeBonds);
      }  // done with degree 3 node
    }  // done finding rings in this fragment

    // calculate the cyclomatic number for the fragment:
    int nexpt = rdcast<int>((nbnds - curFragSize + 1));
    int ssiz = rdcast<int>(resBegins.size() - 1 - resBeginsFragBegin);

    // first check that we got at least the number of expected rings
    // std::cerr<<"EXPT: "<<ssiz<<" "<<nexpt<<std::endl;
    if (ssiz < nexpt) {
      // Issue 3514824: in certain highly fused ring systems, the algorithm
      // above would miss rings.
      // for this fix to apply we have to have at least one non-ring bond
      // that terminates in ring atoms. Find those bonds:
      PRECONDITION(possibleBonds.size() == 0,
                   "Bug with possibleBonds in findSSSR");
      for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
        if (!ringBonds[i]) {
          const BondData &bond = mol.getBond(i);
          if (ringAtoms[bond.getBeginAtomIdx()] &&
              ringAtoms[bond.getEndAtomIdx()]) {
            possibleBonds.push_back(i);
            break;
          }
        }
      }
      deadBonds.reset();
      while (possibleBonds.size()) {
        bool ringFound = FindRings::findRingConnectingAtoms(
            mol, possibleBonds[0], resBegins, resAtoms, invars, ringBonds, ringAtoms);
        if (!ringFound) {
          deadBonds.set(possibleBonds[0], true);
        }
        possibleBonds.clear();
        // check if we need to repeat the process:
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
          if (!ringBonds[i]) {
            const BondData &bond = mol.getBond(i);
            if (!deadBonds[i] &&
                ringAtoms[bond.getBeginAtomIdx()] &&
                ringAtoms[bond.getEndAtomIdx()]) {
              possibleBonds.push_back(i);
              break;
            }
          }
        }
      }
      ssiz = rdcast<int>(resBegins.size() - 1 - resBeginsFragBegin);
      if (ssiz < nexpt) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: could not find number of expected rings. Switching to "
               "an approximate ring finding algorithm."
            << std::endl;
        ringInfo.reset();
        fastFindRings(mol, ringInfo);
        return rdcast<int>(ringInfo.numRings());
      }
    }
    // if we have more than expected we need to do some cleanup
    // otherwise do som clean up work
    // std::cerr<<"  check: "<<ssiz<<" "<<nexpt<<std::endl;
    if (ssiz > nexpt) {
      FindRings::removeExtraRings(resBegins, resAtoms, resBeginsFragBegin,
                                  resAtomsFragBegin, mol, extraRingAtoms,
                                  extraRingBegins);
    }
  }  // done with all fragments

  FindRings::buildRingInfoFromAtoms(mol, ringInfo);

  // update the ring memberships of atoms and bonds in the molecule:
  // store the SSSR rings on the molecule as a property
  // we will ignore any existing SSSRs on the molecule - simply overwrite
  return rdcast<int>(ringInfo.numRings());
}

int symmetrizeSSSR(ROMol &mol, bool includeDativeBonds,
                   bool includeHydrogenBonds) {
  VECT_INT_VECT tmp;
  return symmetrizeSSSR(mol, tmp, includeDativeBonds, includeHydrogenBonds);
};

int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
                   bool includeHydrogenBonds) {
  res.clear();

  int retValue =
      symmetrizeSSSR(mol.asRDMol(), includeDativeBonds, includeHydrogenBonds);

  // Fill in res from the ring info
  const RingInfoCache &ringInfo = mol.asRDMol().getRingInfo();
  const size_t numRings = ringInfo.numRings();
  const uint32_t *data = ringInfo.atomsInRings.data();
  // Ensure the reinterpret_cast is safe, at least for atom indices
  static_assert(sizeof(int) == sizeof(uint32_t));
  const int *intData = reinterpret_cast<const int *>(data);
  res.reserve(numRings);
  for (size_t ringIndex = 0; ringIndex < numRings; ++ringIndex) {
    uint32_t beginIndex = ringInfo.ringBegins[ringIndex];
    uint32_t endIndex = ringInfo.ringBegins[ringIndex + 1];
    res.emplace_back(intData + beginIndex, intData + endIndex);
  }
  return retValue;
}

int symmetrizeSSSR(RDMol &mol, bool includeDativeBonds,
                   bool includeHydrogenBonds) {
  // FIX: need to set flag here the symmetrization has been done in order to
  // avoid repeating this work

  std::vector<atomindex_t> extraRingAtoms;
  std::vector<uint32_t> extraRingBegins;
  findSSSR(mol, mol.getRingInfo(), includeDativeBonds, includeHydrogenBonds,
           &extraRingAtoms, &extraRingBegins);

  // reinit as SYMM_SSSR
  mol.getRingInfo().isInit = true;
  mol.getRingInfo().findRingType = FIND_RING_TYPE_SYMM_SSSR;

  // now check if there are any extra rings on the molecule
  if (extraRingAtoms.size() == 0) {
    // no extra rings nothing to be done
    return rdcast<int>(mol.getRingInfo().numRings());
  }

  // get the rings as bond ids
  const auto &ringBegins = mol.getRingInfo().ringBegins;
  const auto &ringBonds = mol.getRingInfo().bondsInRings;
  const uint32_t numRings = ringBegins.size() - 1;

  //
  // For each "extra" ring, figure out if it could replace a single
  // ring in the SSSR. A ring could be swapped out if:
  //
  // * They are the same size
  // * The replacement doesn't remove any bonds from the union of the bonds
  //   in the SSSR.
  //
  // The latter can be checked by determining if the SSSR ring is the unique
  // provider of any ring bond. If it is, the replacement ring must also
  // provide that bond.
  //
  // May miss extra rings that would need to swap two (or three...) rings
  // to be included.

  // counts of each bond can be determined from bondMembershipBegins
  const auto &bondMembershipBegins = mol.getRingInfo().bondMembershipBegins;

  const uint32_t numExtraRings = extraRingBegins.size() - 1;
  uint32_t maxExtraRingSize = 0;
  for (uint32_t extraRingIdx = 0; extraRingIdx < numExtraRings;
       ++extraRingIdx) {
    const uint32_t extraRingBegin = extraRingBegins[extraRingIdx];
    const uint32_t extraRingEnd = extraRingBegins[extraRingIdx + 1];
    const uint32_t extraRingSize = extraRingEnd - extraRingBegin;
    maxExtraRingSize = std::max(maxExtraRingSize, extraRingSize);
  }
  bool addedRings = false;
  std::vector<uint32_t> extraRing;
  extraRing.reserve(maxExtraRingSize);
  for (uint32_t extraRingIdx = 0; extraRingIdx < numExtraRings; ++extraRingIdx) {
    const uint32_t extraRingBegin = extraRingBegins[extraRingIdx];
    const uint32_t extraRingEnd = extraRingBegins[extraRingIdx + 1];
    const uint32_t extraRingSize = extraRingEnd - extraRingBegin;
    extraRing.resize(extraRingSize);
    RingUtils::convertToBonds(extraRingAtoms.data() + extraRingBegin,
                              extraRingSize, extraRing.data(),
                              mol);
    for (uint32_t ringIdx = 0; ringIdx < numRings; ++ringIdx) {
      uint32_t ringBegin = ringBegins[ringIdx];
      const uint32_t ringEnd = ringBegins[ringIdx + 1];
      const uint32_t ringSize = ringEnd - ringBegin;
      if (ringSize != extraRingSize) {
        continue;
      }

      // If `ring` is the only provider of some bond, extraRing must also
      // provide that bond.
      bool shareBond = false;
      bool replacesAllUniqueBonds = true;
      for (; ringBegin != ringEnd; ++ringBegin) {
        const uint32_t bondID = ringBonds[ringBegin];
        const uint32_t bondCount =
            bondMembershipBegins[bondID + 1] - bondMembershipBegins[bondID];
        if (bondCount == 1 || !shareBond) {
          auto position = find(extraRing.begin(), extraRing.end(), bondID);
          if (position != extraRing.end()) {
            shareBond = true;
          } else if (bondCount == 1) {
            // 1 means `ring` is the only ring in the SSSR to provide this
            // bond, and extraRing did not provide it (so extraRing is not an
            // acceptable substitution in the SSSR for ring)
            replacesAllUniqueBonds = false;
          }
        }
      }

      if (shareBond && replacesAllUniqueBonds) {
        auto &ringAtoms = mol.getRingInfo().atomsInRings;
        ringAtoms.insert(ringAtoms.end(),
                         extraRingAtoms.begin() + extraRingBegin,
                         extraRingAtoms.begin() + extraRingEnd);
        mol.getRingInfo().ringBegins.push_back(ringAtoms.size());
        // TODO: Avoid recomputing the bond info in buildRingInfoFromAtoms
        addedRings = true;
        break;
      }
    }
  }

  if (addedRings) {
    FindRings::buildRingInfoFromAtoms(mol, mol.getRingInfo());
  }

  return rdcast<int>(mol.getRingInfo().numRings());
}

namespace {
void _DFS(const ROMol &mol, const Atom *atom, INT_VECT &atomColors,
          std::vector<const Atom *> &traversalOrder, VECT_INT_VECT &res,
          const Atom *fromAtom = nullptr) {
  // std::cerr<<"  dfs: "<<atom->getIdx()<<" from
  // "<<(fromAtom?fromAtom->getIdx():-1)<<std::endl;
  PRECONDITION(atom, "bad atom");
  PRECONDITION(atomColors[atom->getIdx()] == 0, "bad color");
  atomColors[atom->getIdx()] = 1;
  traversalOrder.push_back(atom);

  for (const auto nbr : mol.atomNeighbors(atom)) {
    unsigned int nbrIdx = nbr->getIdx();
    // std::cerr<<"   "<<atom->getIdx()<<"       consider: "<<nbrIdx<<"
    // "<<atomColors[nbrIdx]<<std::endl;
    if (atomColors[nbrIdx] == 0) {
      if (nbr->getDegree() < 2) {
        atomColors[nbr->getIdx()] = 2;
      } else {
        _DFS(mol, nbr, atomColors, traversalOrder, res, atom);
      }
    } else if (atomColors[nbrIdx] == 1) {
      if (fromAtom && nbrIdx != fromAtom->getIdx()) {
        INT_VECT cycle;
        auto lastElem =
            std::find(traversalOrder.rbegin(), traversalOrder.rend(), atom);
        for (auto rIt = lastElem;  // traversalOrder.rbegin();
             rIt != traversalOrder.rend() && (*rIt)->getIdx() != nbrIdx;
             ++rIt) {
          cycle.push_back((*rIt)->getIdx());
        }
        cycle.push_back(nbrIdx);
        res.push_back(cycle);
        // std::cerr<<"    cycle from "<<atom->getIdx()<<" :";
        // std::copy(cycle.begin(),cycle.end(),std::ostream_iterator<int>(std::cerr,"
        // "));
        // std::cerr<<std::endl;
      }
    }
  }
  atomColors[atom->getIdx()] = 2;
  traversalOrder.pop_back();
  // std::cerr<<"  done "<<atom->getIdx()<<std::endl;
}
}  // end of anonymous namespace
void fastFindRings(const ROMol &mol) {
  // FIXME: Replace this with call to non-compatibility function, using const_cast
  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }

  mol.getRingInfo()->initialize(FIND_RING_TYPE_FAST);

  VECT_INT_VECT res;
  res.resize(0);

  unsigned int nats = mol.getNumAtoms();

  INT_VECT atomColors(nats, 0);

  for (unsigned int i = 0; i < nats; ++i) {
    if (atomColors[i]) {
      continue;
    }
    if (mol.getAtomWithIdx(i)->getDegree() < 2) {
      atomColors[i] = 2;
      continue;
    }
    std::vector<const Atom *> traversalOrder;
    _DFS(mol, mol.getAtomWithIdx(i), atomColors, traversalOrder, res);
  }

  FindRings::storeRingsInfo(mol, res);
}
void fastFindRings(const RDMol &mol, RingInfoCache& rings) {
  std::vector<uint32_t> &ringBegins = rings.ringBegins;
  std::vector<uint32_t> &atomsInRings = rings.atomsInRings;
  std::vector<uint32_t> &bondsInRings = rings.bondsInRings;
  ringBegins.resize(0);
  atomsInRings.resize(0);
  bondsInRings.resize(0);

  uint32_t numAtoms = mol.getNumAtoms();
  uint32_t numBonds = mol.getNumBonds();

  // atomVisited:
  // 0 = unvisited
  // 1 = atom in traversalStack
  // 2 = atom fully visited
  INT_VECT &atomVisited = rings.tempPerAtomInts;
  atomVisited.resize(0);
  atomVisited.resize(numAtoms, 0);

  std::vector<uint32_t> &atomRingCounts = rings.atomMembershipBegins;
  std::vector<uint32_t> &bondRingCounts = rings.bondMembershipBegins;
  atomRingCounts.clear();
  bondRingCounts.clear();
  atomRingCounts.resize(numAtoms, 0);
  bondRingCounts.resize(numBonds, 0);
  //struct RingStackEntry {
  //  uint32_t atomIndex;
  //  uint32_t atomNeighborIndex;
  //};
  // first is atomIndex, second is atomNeighborIndex
  std::vector<std::pair<uint32_t, uint32_t>> &traversalStack =
      rings.tempPerAtomPairs;

  // First, mark all atoms with degree less than 2 as unvisitable.
  for (uint32_t i = 0; i < numAtoms; ++i) {
    if (mol.getAtomDegree(i) < 2) {
      atomVisited[i] = 2;
    }
  }

  for (uint32_t i = 0; i < numAtoms; ++i) {
    if (atomVisited[i] != 0) {
      continue;
    }
    atomVisited[i] = 1;
    traversalStack.resize(1);
    traversalStack[0] = std::make_pair(i, 0);
    while (traversalStack.size() != 0) {
      auto [atomIndex, atomNeighborIndex] = traversalStack.back();
      uint32_t edgeIndex = mol.getAtomBondStarts()[atomIndex] + atomNeighborIndex;
      const uint32_t endEdgeIndex = mol.getAtomBondStarts()[atomIndex + 1];
      uint32_t nbrAtomIndex = mol.getOtherAtomIndices()[edgeIndex];
      //uint32_t bondIndex = mol.getBondDataIndices()[edgeIndex];

      if (atomVisited[nbrAtomIndex] == 0) {
        // Unvisited atom, so recurse
        atomVisited[nbrAtomIndex] = 1;
        traversalStack.push_back(std::make_pair(nbrAtomIndex, 0));
        continue;
      }

      if (atomVisited[nbrAtomIndex] == 1 && traversalStack.size() >= 3 &&
          nbrAtomIndex != traversalStack[traversalStack.size() - 2].first) {
        // Atom already in the traversal stack, but not the previous 2,
        // so we have a cycle.
        assert(traversalStack.back().first == atomIndex);
        ringBegins.push_back(atomsInRings.size());
        atomsInRings.push_back(atomIndex);
        // FIXME: Double-check the order of the bonds in the rings matches the original code!!!
        uint32_t bondIndex = mol.getBondDataIndices()[edgeIndex];
        bondsInRings.push_back(bondIndex);
        ++atomRingCounts[atomIndex];
        ++bondRingCounts[bondIndex];
        size_t i; 
        for (i = traversalStack.size() - 2;
             traversalStack[i].first != nbrAtomIndex; --i) {
          auto [ringAtomIndex, ringNeighborIndex] = traversalStack[i];
          atomsInRings.push_back(ringAtomIndex);
          uint32_t ringEdgeIndex =
              mol.getAtomBondStarts()[ringAtomIndex] + ringNeighborIndex;
          uint32_t ringBondIndex = mol.getBondDataIndices()[ringEdgeIndex];
          bondsInRings.push_back(ringBondIndex);
          ++atomRingCounts[ringAtomIndex];
          ++bondRingCounts[ringBondIndex];
          // This shouldn't loop forever, but we can double-check in a debug
          // build.
          assert(i > 1 || traversalStack[0].first == nbrAtomIndex);
        }
        atomsInRings.push_back(nbrAtomIndex);
        uint32_t ringEdgeIndex =
            mol.getAtomBondStarts()[nbrAtomIndex] + traversalStack[i].second;
        uint32_t ringBondIndex = mol.getBondDataIndices()[ringEdgeIndex];
        bondsInRings.push_back(ringBondIndex);
        ++atomRingCounts[nbrAtomIndex];
        ++bondRingCounts[ringBondIndex];
      }

      // Pick the next available edge, possibly up the stack.
      ++edgeIndex;
      if (edgeIndex != endEdgeIndex) {
        ++traversalStack.back().second;
      } else {
        while (true) {
          // TODO: Check if the next line is correct or not
          atomVisited[traversalStack.back().first] = 2;
          traversalStack.pop_back();
          if (traversalStack.size() == 0) {
            break;
          }
          std::pair<uint32_t, uint32_t> &pair = traversalStack.back();
          auto [atomIndex, atomNeighborIndex] = pair; 
          const uint32_t numNeighbors = mol.getAtomDegree(atomIndex);
          ++atomNeighborIndex;
          if (atomNeighborIndex != numNeighbors) {
            // Found an edge
            pair.second = atomNeighborIndex;
            break;
          }
        }
      }
    }
    atomVisited[i] = 2;
  }
  const size_t numRings = ringBegins.size();
  ringBegins.push_back(atomsInRings.size());

  // Record the rings that each atom and bond are in, using
  // atomRingCounts and bondRingCounts to preallocate,
  // and then iterate through, filling the final arrays in.
  assert(atomsInRings.size() == bondsInRings.size());
  uint32_t sum = 0;
  for (uint32_t i = 0; i < numAtoms; ++i) {
    uint32_t newSum = sum + atomRingCounts[i];
    atomRingCounts[i] = sum;
    sum = newSum;
  }
  assert(sum == atomsInRings.size());
  sum = 0;
  for (uint32_t i = 0; i < numBonds; ++i) {
    uint32_t newSum = sum + bondRingCounts[i];
    bondRingCounts[i] = sum;
    sum = newSum;
  }
  assert(sum == bondsInRings.size());

  // Fill in the atom-to-ring and bond-to-ring mappings
  std::vector<uint32_t> &atomMemberships = rings.atomMemberships;
  std::vector<uint32_t> &bondMemberships = rings.bondMemberships;
  atomMemberships.resize(atomsInRings.size());
  bondMemberships.resize(bondsInRings.size());
  for (size_t ring = 0; ring < numRings; ++ring) {
    uint32_t ringBegin = ringBegins[ring];
    uint32_t ringEnd = ringBegins[ring+1];
    for (; ringBegin < ringEnd; ++ringBegin) {
      uint32_t atom = atomsInRings[ringBegin];
      uint32_t destIndex = atomRingCounts[atom];
      atomMemberships[destIndex] = ring;
      atomRingCounts[atom] = destIndex + 1;

      uint32_t bond = bondsInRings[ringBegin];
      destIndex = bondRingCounts[bond];
      bondMemberships[destIndex] = ring;
      bondRingCounts[bond] = destIndex + 1;
    }
  }
  // atomRingCounts and bondRingCounts effectively got shifted over by
  // the incrementing, so add back the 0 at the beginning.
  atomRingCounts.push_back(0);
  memmove(atomRingCounts.data() + 1, atomRingCounts.data(),
          sizeof(uint32_t) * (atomRingCounts.size() - 1));
  atomRingCounts[0] = 0;

  bondRingCounts.push_back(0);
  memmove(bondRingCounts.data() + 1, bondRingCounts.data(),
          sizeof(uint32_t) * (bondRingCounts.size() - 1));
  bondRingCounts[0] = 0;

  rings.initFusedInfoFromBondMemberships();

  rings.isInit = true;
  rings.findRingType = FIND_RING_TYPE_FAST;
}

#ifdef RDK_USE_URF
void findRingFamilies(const ROMol &mol) {
  if (mol.getRingInfo()->isInitialized()) {
    // return if we've done this before
    if (mol.getRingInfo()->areRingFamiliesInitialized()) {
      return;
    }
  } else {
    mol.getRingInfo()->initialize();
  }

  RDL_graph *graph = RDL_initNewGraph(mol.getNumAtoms());
  for (ROMol::ConstBondIterator cbi = mol.beginBonds(); cbi != mol.endBonds();
       ++cbi) {
    RDL_addUEdge(graph, (*cbi)->getBeginAtomIdx(), (*cbi)->getEndAtomIdx());
  }
  RDL_data *urfdata = RDL_calculate(graph);
  if (urfdata == nullptr) {
    RDL_deleteGraph(graph);
    mol.getRingInfo()->dp_urfData.reset();
    throw ValueErrorException("Cannot get URFs");
  }
  mol.getRingInfo()->dp_urfData.reset(urfdata, &RDL_deleteData);
  for (unsigned int i = 0; i < RDL_getNofURF(urfdata); ++i) {
    RDL_node *nodes = nullptr;
    unsigned nNodes = RDL_getNodesForURF(urfdata, i, &nodes);
    if (nNodes == RDL_INVALID_RESULT) {
      free(nodes);
      throw ValueErrorException("Cannot get URF nodes");
    }
    RDL_edge *edges = nullptr;
    unsigned nEdges = RDL_getEdgesForURF(urfdata, i, &edges);
    if (nEdges == RDL_INVALID_RESULT) {
      free(nodes);
      free(edges);
      throw ValueErrorException("Cannot get URF edges");
    }
    INT_VECT nvect(nNodes), evect(nEdges);
    for (unsigned int ridx = 0; ridx < nNodes; ++ridx) {
      nvect[ridx] = nodes[ridx];
    }
    for (unsigned int ridx = 0; ridx < nEdges; ++ridx) {
      unsigned int bidx = edges[ridx][0];
      unsigned int eidx = edges[ridx][1];
      evect[ridx] = mol.getBondBetweenAtoms(bidx, eidx)->getIdx();
    }
    mol.getRingInfo()->addRingFamily(nvect, evect);
    free(nodes);
    free(edges);
  }
}
#else
void findRingFamilies(const ROMol &mol) {
  BOOST_LOG(rdErrorLog)
      << "This version of the RDKit was built without URF support" << std::endl;
}
#endif
}  // namespace MolOps

}  // namespace RDKit
