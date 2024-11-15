//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Rings.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>

#include <RDGeneral/utils.h>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>

using RINGINVAR = boost::dynamic_bitset<>;
using RINGINVAR_SET = std::set<RINGINVAR>;
using RINGINVAR_VECT = std::vector<RINGINVAR>;

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

void convertToBonds(const INT_VECT &ring, INT_VECT &bondRing,
                    const ROMol &mol) {
  const auto rsiz = rdcast<unsigned int>(ring.size());
  bondRing.resize(rsiz);
  for (unsigned int i = 0; i < (rsiz - 1); i++) {
    const Bond *bnd = mol.getBondBetweenAtoms(ring[i], ring[i + 1]);
    if (!bnd) {
      throw ValueErrorException("expected bond not found");
    }
    bondRing[i] = bnd->getIdx();
  }
  // bond from last to first atom
  const Bond *bnd = mol.getBondBetweenAtoms(ring[rsiz - 1], ring[0]);
  if (!bnd) {
    throw ValueErrorException("expected bond not found");
  }

  bondRing[rsiz - 1] = bnd->getIdx();
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
  int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings,
                       boost::dynamic_bitset<> &activeBonds,
                       INT_VECT *forbidden = nullptr);

 private:
  INT_VECT d_parents;
  std::vector<unsigned int> d_depths;
  INT_VECT d_done;
};

void trimBonds(unsigned int cand, const ROMol &tMol, INT_SET &changed,
               INT_VECT &atomDegrees, boost::dynamic_bitset<> &activeBonds);
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

using RINGINVAR_INT_VECT_MAP = std::map<RINGINVAR, INT_VECT>;

void findSSSRforDupCands(const ROMol &mol, VECT_INT_VECT &res,
                         RINGINVAR_SET &invars, const INT_INT_VECT_MAP dupMap,
                         const RINGINVAR_INT_VECT_MAP &dupD2Cands,
                         INT_VECT &atomDegrees,
                         boost::dynamic_bitset<> activeBonds) {
  BFSWorkspace bfs_workspace;
  for (const auto &dupD2Cand : dupD2Cands) {
    const INT_VECT &dupCands = dupD2Cand.second;
    if (dupCands.size() > 1) {
      // we have duplicate candidates.
      VECT_INT_VECT nrings;
      auto minSiz = static_cast<unsigned int>(MAX_INT);
      for (int dupCand : dupCands) {
        // now break bonds for all the d2 nodes for that give the same rings as
        // with (*dupi) and recompute smallest ring with (*dupi)
        INT_VECT atomDegreesCopy = atomDegrees;
        boost::dynamic_bitset<> activeBondsCopy = activeBonds;
        INT_SET changed;
        auto dmci = dupMap.find(dupCand);
        CHECK_INVARIANT(dmci != dupMap.end(), "duplicate could not be found");
        for (int dni : dmci->second) {
          trimBonds(dni, mol, changed, atomDegreesCopy, activeBondsCopy);
        }

        // now find the smallest ring/s around (*dupi)
        VECT_INT_VECT srings;
        bfs_workspace.smallestRingsBfs(mol, dupCand, srings, activeBondsCopy);
        for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end();
             ++sri) {
          if (sri->size() < minSiz) {
            minSiz = rdcast<unsigned int>(sri->size());
          }
          nrings.push_back((*sri));
        }
      }
      for (const auto &nring : nrings) {
        if (nring.size() == minSiz) {
          auto invr = RingUtils::computeRingInvariant(nring, mol.getNumAtoms());
#if 0
          std::cerr << "    res: " << invr << " | ";
          std::copy(nring.begin(), nring.end(),
                    std::ostream_iterator<int>(std::cerr, " "));
          std::cerr << std::endl;
#endif
          if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
          }
        }
      }  // end of loop over new rings found
    }  // end if (dupCand.size() > 1)
  }  // end of loop over all set of duplicate candidates
}

auto compRingSize = [](const auto &v1, const auto &v2) {
  return v1.size() < v2.size();
};

void removeExtraRings(VECT_INT_VECT &res, unsigned int, const ROMol &mol) {
  // sort on size
  std::sort(res.begin(), res.end(), compRingSize);

#if 0
        std::cerr<<"\n\nSORTED\n";
        for(VECT_INT_VECT::const_iterator iter=res.begin();
            iter!=res.end();++iter){
          std::cerr<<iter-res.begin()<<": ";
          std::copy(iter->begin(),iter->end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;
        }
#endif

  // change the rings from atom IDs to bondIds
  VECT_INT_VECT brings;
  RingUtils::convertToBonds(res, brings, mol);
  std::vector<boost::dynamic_bitset<>> bitBrings;
  bitBrings.reserve(brings.size());
  for (VECT_INT_VECT_CI vivi = brings.begin(); vivi != brings.end(); ++vivi) {
    boost::dynamic_bitset<> lring(mol.getNumBonds());
    for (int ivi : *vivi) {
      lring.set(ivi);
    }
    bitBrings.push_back(lring);
  }

  boost::dynamic_bitset<> availRings(res.size());
  availRings.set();
  boost::dynamic_bitset<> keepRings(res.size());
  boost::dynamic_bitset<> munion(mol.getNumBonds());

  // optimization - don't reallocate a new one each loop
  boost::dynamic_bitset<> workspace(mol.getNumBonds());

  for (unsigned int i = 0; i < res.size(); ++i) {
    // skip this ring if we've already seen all of its bonds
    if (bitBrings[i].is_subset_of(munion)) {
      availRings.set(i, 0);
    }
    if (!availRings[i]) {
      continue;
    }

    munion |= bitBrings[i];
    keepRings.set(i);

    // from this ring we consider all others that are still available and the
    // same size
    boost::dynamic_bitset<> consider(res.size());
    for (unsigned int j = i + 1; j < res.size(); ++j) {
      // std::cerr<<"        "<<j<<" "<<brings[j].size()<<" -
      // "<<brings[i].size()<<"  >"<<availRings[j]<<std::endl;
      if (availRings[j] && (brings[j].size() == brings[i].size())) {
        consider.set(j);
      }
    }

    while (consider.any()) {
      unsigned int bestJ = i + 1;
      int bestOverlap = -1;
      // loop over the available other rings in consideration and pick the one
      // that has the most overlapping bonds with what we've done so far.
      // this is the fix to github #526
      for (unsigned int j = i + 1;
           j < res.size() && brings[j].size() == brings[i].size(); ++j) {
        if (!consider[j] || !availRings[j]) {
          continue;
        }
        workspace = bitBrings[j];
        workspace &= munion;
        int overlap = rdcast<int>(workspace.count());
        if (overlap > bestOverlap) {
          bestOverlap = overlap;
          bestJ = j;
        }
      }
      consider.set(bestJ, 0);
      if (bitBrings[bestJ].is_subset_of(munion)) {
        availRings.set(bestJ, 0);
      } else {
        keepRings.set(bestJ);
        availRings.set(bestJ, 0);
        munion |= bitBrings[bestJ];
      }
    }
  }
  // remove the extra rings from res and store them on the molecule in case we
  // wish symmetrize the SSSRs later
  VECT_INT_VECT extras;
  VECT_INT_VECT temp = res;
  res.resize(0);
  for (unsigned int i = 0; i < temp.size(); i++) {
    if (keepRings[i]) {
      res.push_back(temp[i]);
    } else {
      extras.push_back(temp[i]);
    }
  }

  mol.setProp(common_properties::extraRings, extras, true);
}

void findRingsD2nodes(const ROMol &tMol, VECT_INT_VECT &res,
                      RINGINVAR_SET &invars, const INT_VECT &d2nodes,
                      INT_VECT &atomDegrees,
                      boost::dynamic_bitset<> &activeBonds,
                      boost::dynamic_bitset<> &ringBonds,
                      boost::dynamic_bitset<> &ringAtoms) {
  // place to record any duplicate rings discovered from the current d2 nodes
  RINGINVAR_INT_VECT_MAP dupD2Cands;
  INT_SET changed;

  INT_INT_VECT_MAP dupMap;
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
  std::map<int, RINGINVAR_VECT> nodeInvars;
  std::map<int, RINGINVAR_VECT>::const_iterator nici;
  BFSWorkspace bfs_workspace;
  for (auto &cand : d2nodes) {
    // std::cerr<<"    smallest rings bfs: "<<cand<<std::endl;
    VECT_INT_VECT srings;
    // we have to find all non duplicate possible smallest rings for each node
    bfs_workspace.smallestRingsBfs(tMol, cand, srings, activeBonds);
    for (const auto &nring : srings) {
      auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
      if (invars.find(invr) == invars.end()) {
        res.push_back(nring);
        invars.insert(invr);
        for (unsigned int i = 0; i < nring.size() - 1; ++i) {
          unsigned int bIdx =
              tMol.getBondBetweenAtoms(nring[i], nring[i + 1])->getIdx();
          ringBonds.set(bIdx);
          ringAtoms.set(nring[i]);
        }
        ringBonds.set(
            tMol.getBondBetweenAtoms(nring[0], nring[nring.size() - 1])
                ->getIdx());
        ringAtoms.set(nring[nring.size() - 1]);
#if 0
          std::cerr<<"    res: "<<invr<<" | ";
          std::copy(nring.begin(),nring.end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;
#endif
      }

      nodeInvars[cand].push_back(invr);
      // check if this ring is duplicate with something else
      for (auto &nici : nodeInvars) {
        if (nici.first != cand) {
          if (std::find(nici.second.begin(), nici.second.end(), invr) !=
              nici.second.end()) {
            // ok we discovered this ring via another node before
            // add that node as duplicate to this node and vice versa
            dupMap[cand].push_back(nici.first);
            dupMap[nici.first].push_back(cand);
          }
        }
      }
      dupD2Cands[invr].push_back(cand);
    }

    // We don't want to trim the bonds connecting cand here - this can disrupt
    // a second small ring. Here is an example SC(C3C1CC(C3)CC(C2S)(O)C1)2S
    // by trimming the bond connecting to atom #4, we lose the smallest ring
    // that contains atom #7. Issue 134

    // But if there were no rings found, trimming isn't dangerous, and can
    // save wasted time for long chains.
    if (srings.empty()) {
      changed = {cand};
      while (!changed.empty()) {
        int cand = *(changed.begin());
        changed.erase(changed.begin());
        trimBonds(cand, tMol, changed, atomDegrees, activeBonds);
      }
    }
  }

  // now deal with any d2 nodes that resulted in duplicate rings before trimming
  // their bonds.
  // it is possible that one of these nodes is involved a different small ring,
  // that is not found  because the first nodes has not be trimmed. Here is an
  // example molecule:
  // CC1=CC=C(C=C1)S(=O)(=O)O[CH]2[CH]3CO[CH](O3)[CH]4OC(C)(C)O[CH]24
  findSSSRforDupCands(tMol, res, invars, dupMap, dupD2Cands, atomDegrees,
                      activeBonds);
}

void findRingsD3Node(const ROMol &tMol, VECT_INT_VECT &res,
                     RINGINVAR_SET &invars, int cand, INT_VECT &,
                     boost::dynamic_bitset<> activeBonds) {
  // this is brutal - we have no degree 2 nodes - find the first possible degree
  // 3 node
  int nsmall;

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

  // first find all smallest possible rings
  VECT_INT_VECT srings;
  BFSWorkspace bfs_workspace;
  nsmall = bfs_workspace.smallestRingsBfs(tMol, cand, srings, activeBonds);

  for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end(); ++sri) {
    const INT_VECT &nring = (*sri);
    auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
    if (invars.find(invr) == invars.end()) {
      res.push_back(nring);
      invars.insert(invr);
    }
  }

  // if already found >3 rings we are done with this degree 3 node
  // if we found less than 3 we have to find other potential ring/s
  if (nsmall < 3) {
    int n1 = -1, n2 = -1, n3 = -1;

    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = tMol.getAtomBonds(tMol.getAtomWithIdx(cand));
    while (beg != end && !activeBonds[tMol[*beg]->getIdx()]) {
      ++beg;
    }
    CHECK_INVARIANT(beg != end, "neighbor not found");
    n1 = tMol[*beg]->getOtherAtomIdx(cand);

    ++beg;
    while (beg != end && !activeBonds[tMol[*beg]->getIdx()]) {
      ++beg;
    }
    CHECK_INVARIANT(beg != end, "neighbor not found");
    n2 = tMol[*beg]->getOtherAtomIdx(cand);

    ++beg;
    while (beg != end && !activeBonds[tMol[*beg]->getIdx()]) {
      ++beg;
    }
    CHECK_INVARIANT(beg != end, "neighbor not found");
    n3 = tMol[*beg]->getOtherAtomIdx(cand);

    if (nsmall == 2) {
      // we found two rings find the third one
      // first find the neighbor that is common to the two ring we found so far
      int f = -1;

      if ((std::find(srings[0].begin(), srings[0].end(), n1) !=
           srings[0].end()) &&
          (std::find(srings[1].begin(), srings[1].end(), n1) !=
           srings[1].end())) {
        f = n1;
      } else if ((std::find(srings[0].begin(), srings[0].end(), n2) !=
                  srings[0].end()) &&
                 (std::find(srings[1].begin(), srings[1].end(), n2) !=
                  srings[1].end())) {
        f = n2;
      } else if ((std::find(srings[0].begin(), srings[0].end(), n3) !=
                  srings[0].end()) &&
                 (std::find(srings[1].begin(), srings[1].end(), n3) !=
                  srings[1].end())) {
        f = n3;
      }
      CHECK_INVARIANT(f >= 0, "third ring not found");

      // now find the smallest possible ring that does not contain f
      VECT_INT_VECT trings;
      INT_VECT forb;
      forb.push_back(f);
      bfs_workspace.smallestRingsBfs(tMol, cand, trings, activeBonds, &forb);
      for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
        const INT_VECT &nring = (*sri);
        auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());

        if (invars.find(invr) == invars.end()) {
          res.push_back(nring);
          invars.insert(invr);
        }
      }
    }  // doing degree 3 node  - end of 2 smallest rings found for cand
    if (nsmall == 1) {
      // we found 1 ring - we need to find two more that involve the 3rd
      // neighbor
      int f1 = -1, f2 = -1;
      // Which of our three neighbors are in the small ring?
      //   these are f1 and f2
      if (std::find(srings[0].begin(), srings[0].end(), n1) ==
          srings[0].end()) {
        f1 = n2, f2 = n3;
      } else if (std::find(srings[0].begin(), srings[0].end(), n2) ==
                 srings[0].end()) {
        f1 = n1;
        f2 = n3;
      } else if (std::find(srings[0].begin(), srings[0].end(), n3) ==
                 srings[0].end()) {
        f1 = n1;
        f2 = n2;
      }
      CHECK_INVARIANT(f1 >= 0, "rings not found");
      CHECK_INVARIANT(f2 >= 0, "rings not found");

      // now find two rings that include cand, one of these rings should include
      // f1
      // and the other should include f2

      // first ring with f1 and no f2
      VECT_INT_VECT trings;
      INT_VECT forb;
      forb.push_back(f2);
      bfs_workspace.smallestRingsBfs(tMol, cand, trings, activeBonds, &forb);
      for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
        const INT_VECT &nring = (*sri);
        auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
        if (invars.find(invr) == invars.end()) {
          res.push_back(nring);
          invars.insert(invr);
        }
      }

      // next the ring with f2 and no f1
      trings.clear();
      forb.clear();
      forb.push_back(f1);
      bfs_workspace.smallestRingsBfs(tMol, cand, trings, activeBonds, &forb);
      for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
        const INT_VECT &nring = (*sri);
        auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
        if (invars.find(invr) == invars.end()) {
          res.push_back(nring);
          invars.insert(invr);
        }
      }
    }  // doing node of degree 3 - end of found only 1 smallest ring
  }  // end of found less than 3 smallest ring for the degree 3 node
}

int greatestComFac(long curfac, long nfac) {
  long small;
  long large;
  long rem;

  // Determine which of the numbers is the larger, and which is the smaller
  large = (curfac > nfac) ? curfac : nfac;
  small = (curfac < nfac) ? curfac : nfac;

  // Keep looping until no remainder, as this means it is a factor of both
  while (small != 0) {
    // Set the larger var to the smaller, and set the smaller to the remainder
    // of (large / small)
    rem = (large % small);
    large = small;
    small = rem;
  }

  // By here nLarge will hold the largest common factor, so just return it
  return large;
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
void trimBonds(unsigned int cand, const ROMol &tMol, INT_SET &changed,
               INT_VECT &atomDegrees, boost::dynamic_bitset<> &activeBonds) {
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = tMol.getAtomBonds(tMol.getAtomWithIdx(cand));
  while (beg != end) {
    const Bond *bond = tMol[*beg];
    ++beg;
    if (!activeBonds[bond->getIdx()]) {
      continue;
    }
    unsigned int oIdx = bond->getOtherAtomIdx(cand);
    if (atomDegrees[oIdx] <= 2) {
      changed.insert(oIdx);
    }
    activeBonds[bond->getIdx()] = 0;
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
 *  rings - list of rings into which the results are entered
 *  forbidden - list of atoms ID that should be avoided
 *
 * RETURNS:
 *  number of smallest rings found
 ***********************************************************************************/
int BFSWorkspace::smallestRingsBfs(const ROMol &mol, int root,
                                   VECT_INT_VECT &rings,
                                   boost::dynamic_bitset<> &activeBonds,
                                   INT_VECT *forbidden) {
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

  std::deque<int> bfsq;
  bfsq.push_back(root);

  INT_VECT ring;

  unsigned int curSize = UINT_MAX;
  while (bfsq.size() > 0) {
    if (bfsq.size() >= RingUtils::MAX_BFSQ_SIZE) {
      std::string msg =
          "Maximum BFS search size exceeded.\nThis is likely due to a highly "
          "symmetric fused ring system.";
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw ValueErrorException(msg);
    }

    const int curr = bfsq.front();
    bfsq.pop_front();
    d_done[curr] = BLACK;

    const unsigned int depth = d_depths[curr] + 1;
    if (depth > curSize) {
      // depth is the shortest cycle I _could_ find this round.
      break;
    }

    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(mol.getAtomWithIdx(curr));
    while (beg != end) {
      const Bond *bond = mol[*beg];
      ++beg;
      if (!activeBonds[bond->getIdx()]) {
        continue;
      }
      int nbrIdx = bond->getOtherAtomIdx(curr);
      if (d_done[nbrIdx] == BLACK || d_parents[curr] == nbrIdx) {
        continue;
      }
      if (d_done[nbrIdx] == WHITE) {
        // we have never been to this node before through via any path
        d_parents[nbrIdx] = curr;
        d_done[nbrIdx] = GRAY;
        d_depths[nbrIdx] = depth;
        bfsq.push_back(nbrIdx);
      } else {
        // we have been here via a different path
        // there is a potential for ring closure here
        // stitch together the two paths

        ring = {nbrIdx};
        // forwards path
        int parent = d_parents[nbrIdx];
        while (parent != -1 && parent != root) {
          ring.push_back(parent);
          parent = d_parents[parent];
        }

        // backwards path
        ring.insert(ring.begin(), curr);
        parent = d_parents[curr];
        while (parent != -1) {
          // Is the least common ancestor not the root?
          if (std::find(ring.begin(), ring.end(), parent) != ring.end()) {
            ring.clear();
            break;
          }
          ring.insert(ring.begin(), parent);
          parent = d_parents[parent];
        }

        // Found a new small ring including the root.
        if (ring.size() > 1) {
          if (ring.size() <= curSize) {
            curSize = rdcast<unsigned int>(ring.size());
            rings.push_back(ring);
          } else {
            // we are done with the smallest rings
            return rdcast<unsigned int>(rings.size());
          }
        }
      }
    }  // end of loop over neighbors of current atom
  }  // moving to the next node

  // if we are here we should have found everything around the node
  return rdcast<unsigned int>(rings.size());
}

bool _atomSearchBFS(const ROMol &tMol, unsigned int startAtomIdx,
                    unsigned int endAtomIdx, boost::dynamic_bitset<> &ringAtoms,
                    INT_VECT &res, RINGINVAR_SET &invars) {
  res.clear();
  std::deque<INT_VECT> bfsq;

  INT_VECT tv;
  tv.push_back(startAtomIdx);
  bfsq.push_back(tv);
  while (!bfsq.empty()) {
    if (bfsq.size() >= RingUtils::MAX_BFSQ_SIZE) {
      std::string msg =
          "Maximum BFS search size exceeded.\nThis is likely due to a highly "
          "symmetric fused ring system.";
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw ValueErrorException(msg);
    }
    tv = bfsq.front();
    bfsq.pop_front();

    unsigned int currAtomIdx = tv.back();
    ROMol::ADJ_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) =
        tMol.getAtomNeighbors(tMol.getAtomWithIdx(currAtomIdx));
    while (nbrIdx != endNbrs) {
      if (*nbrIdx == endAtomIdx) {
        if (currAtomIdx != startAtomIdx) {
          INT_VECT nv(tv);

          nv.push_back(rdcast<unsigned int>(*nbrIdx));
          // make sure the ring we just found isn't already in our set
          // of rings (this was an extension of sf.net issue 249)
          auto invr = RingUtils::computeRingInvariant(nv, tMol.getNumAtoms());
          if (invars.find(invr) == invars.end()) {
            // we're done!
            res.resize(nv.size());
            std::copy(nv.begin(), nv.end(), res.begin());
            return true;
          }
        } else {
          // ignore this one
        }
      } else if (ringAtoms[*nbrIdx] &&
                 std::find(tv.begin(), tv.end(), *nbrIdx) == tv.end()) {
        //} else if(ringAtoms[*nbrIdx]){
        INT_VECT nv(tv);
        nv.push_back(rdcast<unsigned int>(*nbrIdx));

        bfsq.push_back(nv);
      }
      ++nbrIdx;
    }
  }
  return false;
}

bool findRingConnectingAtoms(const ROMol &tMol, const Bond *bond,
                             VECT_INT_VECT &res, RINGINVAR_SET &invars,
                             boost::dynamic_bitset<> &ringBonds,
                             boost::dynamic_bitset<> &ringAtoms) {
  PRECONDITION(bond, "bad bond");
  PRECONDITION(!ringBonds[bond->getIdx()], "not a ring bond");
  PRECONDITION(ringAtoms[bond->getBeginAtomIdx()], "not a ring atom");
  PRECONDITION(ringAtoms[bond->getEndAtomIdx()], "not a ring atom");

  INT_VECT nring;
  if (_atomSearchBFS(tMol, bond->getBeginAtomIdx(), bond->getEndAtomIdx(),
                     ringAtoms, nring, invars)) {
    auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
    if (invars.find(invr) == invars.end()) {
      res.push_back(nring);
      invars.insert(invr);
#if 0
        std::cerr<<"    local: "<<invr<<" | ";
        std::copy(nring.begin(),nring.end(),std::ostream_iterator<int>(std::cerr," "));
        std::cerr<<std::endl;
#endif
      for (unsigned int i = 0; i < nring.size() - 1; ++i) {
        unsigned int bIdx =
            tMol.getBondBetweenAtoms(nring[i], nring[i + 1])->getIdx();
        ringBonds.set(bIdx);
        ringAtoms.set(nring[i]);
      }
      ringBonds.set(tMol.getBondBetweenAtoms(nring[0], nring[nring.size() - 1])
                        ->getIdx());
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
int findSSSR(const ROMol &mol, VECT_INT_VECT *res, bool includeDativeBonds) {
  if (!res) {
    VECT_INT_VECT rings;
    return findSSSR(mol, rings, includeDativeBonds);
  } else {
    return findSSSR(mol, (*res), includeDativeBonds);
  }
}

int findSSSR(const ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds) {
  res.resize(0);
  if (mol.getRingInfo()->isInitialized()) {
    mol.getRingInfo()->reset();
  }
  mol.getRingInfo()->initialize(FIND_RING_TYPE_SSSR);
  RINGINVAR_SET invars;

  unsigned int nats = mol.getNumAtoms();
  boost::dynamic_bitset<> activeAtoms(nats);
  activeAtoms.set();
  int nbnds = mol.getNumBonds();
  boost::dynamic_bitset<> activeBonds(nbnds);
  activeBonds.set();

  // Zero-order bonds are not candidates for rings, and dative bonds may also be
  // out
  ROMol::EDGE_ITER firstB, lastB;
  boost::tie(firstB, lastB) = mol.getEdges();
  while (firstB != lastB) {
    const Bond *bond = mol[*firstB];
    if (bond->getBondType() == Bond::ZERO ||
        (!includeDativeBonds && isDative(*bond))) {
      activeBonds[bond->getIdx()] = 0;
    }
    ++firstB;
  }

  boost::dynamic_bitset<> ringBonds(nbnds);
  boost::dynamic_bitset<> ringAtoms(nats);

  INT_VECT atomDegrees(nats);
  INT_VECT atomDegreesWithZeroOrderBonds(nats);
  for (unsigned int i = 0; i < nats; ++i) {
    const Atom *atom = mol.getAtomWithIdx(i);
    int deg = atom->getDegree();
    atomDegrees[i] = deg;
    atomDegreesWithZeroOrderBonds[i] = deg;
    for (const auto bond : mol.atomBonds(atom)) {
      if (bond->getBondType() == Bond::ZERO ||
          (!includeDativeBonds && isDative(*bond))) {
        atomDegrees[i]--;
      }
    }
  }

  // find the number of fragments in the molecule - we will loop over them
  VECT_INT_VECT frags;
  INT_VECT curFrag;
  unsigned int nfrags = getMolFrags(mol, frags);
  // loop over the fragments in a molecule
  for (unsigned int fi = 0; fi < nfrags; ++fi) {
    VECT_INT_VECT fragRes;
    curFrag = frags[fi];

    if (curFrag.size() < 3) {
      continue;
    }

    // the following is the list of atoms that are useful in the next round of
    // trimming basically atoms that become degree 0 or 1 because of bond
    // removals initialized with atoms of degrees 0 and 1
    INT_SET changed;
    int bndcnt_with_zero_order_bonds = 0;
    unsigned int nbnds = 0;
    for (auto atom_idx : curFrag) {
      bndcnt_with_zero_order_bonds += atomDegreesWithZeroOrderBonds[atom_idx];

      int deg = atomDegrees[atom_idx];

      nbnds += deg;
      if (deg < 2) {
        changed.insert(atom_idx);
      }
    }

    // check to see if this fragment can even have a possible ring
    CHECK_INVARIANT(bndcnt_with_zero_order_bonds % 2 == 0,
                    "fragment graph has a dangling degree");
    bndcnt_with_zero_order_bonds = bndcnt_with_zero_order_bonds / 2;
    int num_possible_rings = bndcnt_with_zero_order_bonds - curFrag.size() + 1;
    if (num_possible_rings < 1) {
      continue;
    }

    CHECK_INVARIANT(nbnds % 2 == 0,
                    "fragment graph problem when including zero-order bonds");
    nbnds = nbnds / 2;

    boost::dynamic_bitset<> doneAts(nats);
    unsigned int nAtomsDone = 0;
    while (nAtomsDone < curFrag.size()) {
      // std::cerr<<" ndone: "<<nAtomsDone<<std::endl;
      // std::cerr<<" activeBonds: "<<activeBonds<<std::endl;
      // std::cerr<<"  done: ";
      // trim all bonds that connect to degree 0 and 1 atoms
      while (changed.size() > 0) {
        int cand = *(changed.begin());
        changed.erase(changed.begin());
        if (!doneAts[cand]) {
          // std::cerr<<cand<<" ";
          doneAts.set(cand);
          ++nAtomsDone;
          FindRings::trimBonds(cand, mol, changed, atomDegrees, activeBonds);
        }
      }
      // std::cerr<<std::endl;
      // std::cerr<<"activeBonds2: "<<activeBonds<<std::endl;

      // all atoms left in the fragment should at least have a degree >= 2
      // collect all the degree two nodes;
      INT_VECT d2nodes;
      FindRings::pickD2Nodes(mol, d2nodes, curFrag, atomDegrees, activeBonds);
#if 0
          std::cerr<<"d2nodes: ";
          std::copy(d2nodes.begin(),d2nodes.end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;
#endif
      if (d2nodes.size() > 0) {  // deal with the current degree two nodes
        // place to record any duplicate rings discovered from the current d2
        // nodes
        FindRings::findRingsD2nodes(mol, fragRes, invars, d2nodes, atomDegrees,
                                    activeBonds, ringBonds, ringAtoms);
#if 0
            std::cerr<<"  d2nodes post: ";
            std::copy(d2nodes.begin(),d2nodes.end(),std::ostream_iterator<int>(std::cerr," "));
            std::cerr<<std::endl;
            std::cerr<<"  ring bonds: "<<ringBonds<<std::endl;
#endif
        INT_VECT_CI d2i;
        // trim after we have dealt with all the current d2 nodes,
        for (d2i = d2nodes.begin(); d2i != d2nodes.end(); d2i++) {
          doneAts.set(*d2i);
          ++nAtomsDone;
          FindRings::trimBonds((*d2i), mol, changed, atomDegrees, activeBonds);
        }
      }  // end of degree two nodes
      else if (nAtomsDone <
               curFrag.size()) {  // now deal with higher degree nodes
        // this is brutal - we have no degree 2 nodes - find the first
        // possible degree 3 node
        int cand = -1;
        for (INT_VECT_CI aidi = curFrag.begin(); aidi != curFrag.end();
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
        FindRings::findRingsD3Node(mol, fragRes, invars, cand, atomDegrees,
                                   activeBonds);
        doneAts.set(cand);
        ++nAtomsDone;
        FindRings::trimBonds(cand, mol, changed, atomDegrees, activeBonds);
      }  // done with degree 3 node
    }  // done finding rings in this fragment

#if 0
        std::cerr<<"\n\nFOUND:\n";
        for(VECT_INT_VECT::const_iterator iter=fragRes.begin();
            iter!=fragRes.end();++iter){
          std::cerr<<iter-fragRes.begin()<<": ";
          std::copy(iter->begin(),iter->end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;
        }
#endif
    // calculate the cyclomatic number for the fragment:
    int nexpt = rdcast<int>((nbnds - curFrag.size() + 1));
    int ssiz = rdcast<int>(fragRes.size());

    // first check that we got at least the number of expected rings
    // std::cerr<<"EXPT: "<<ssiz<<" "<<nexpt<<std::endl;
    if (ssiz < nexpt) {
      // Issue 3514824: in certain highly fused ring systems, the algorithm
      // above would miss rings.
      // for this fix to apply we have to have at least one non-ring bond
      // that terminates in ring atoms. Find those bonds:
      std::vector<const Bond *> possibleBonds;
      for (unsigned int i = 0; i < nbnds; ++i) {
        if (!ringBonds[i]) {
          const Bond *bnd = mol.getBondWithIdx(i);
          if (ringAtoms[bnd->getBeginAtomIdx()] &&
              ringAtoms[bnd->getEndAtomIdx()]) {
            possibleBonds.push_back(bnd);
            break;
          }
        }
      }
      boost::dynamic_bitset<> deadBonds(mol.getNumBonds());
      while (possibleBonds.size()) {
        bool ringFound = FindRings::findRingConnectingAtoms(
            mol, possibleBonds[0], fragRes, invars, ringBonds, ringAtoms);
        if (!ringFound) {
          deadBonds.set(possibleBonds[0]->getIdx(), 1);
        }
        possibleBonds.clear();
        // check if we need to repeat the process:
        for (unsigned int i = 0; i < nbnds; ++i) {
          if (!ringBonds[i]) {
            const Bond *bnd = mol.getBondWithIdx(i);
            if (!deadBonds[bnd->getIdx()] &&
                ringAtoms[bnd->getBeginAtomIdx()] &&
                ringAtoms[bnd->getEndAtomIdx()]) {
              possibleBonds.push_back(bnd);
              break;
            }
          }
        }
      }
      ssiz = rdcast<int>(fragRes.size());
      if (ssiz < nexpt) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: could not find number of expected rings. Switching to "
               "an approximate ring finding algorithm."
            << std::endl;
        mol.getRingInfo()->reset();
        fastFindRings(mol);
        res.clear();
        res = mol.getRingInfo()->atomRings();
        return rdcast<int>(res.size());
      }
    }
    // if we have more than expected we need to do some cleanup
    // otherwise do som clean up work
    // std::cerr<<"  check: "<<ssiz<<" "<<nexpt<<std::endl;
    if (ssiz > nexpt) {
      FindRings::removeExtraRings(fragRes, nexpt, mol);
    }

#if 0
        std::cerr<<"\n\nKEEPING:\n";
        for(VECT_INT_VECT::const_iterator iter=fragRes.begin();
            iter!=fragRes.end();++iter){
          std::copy(iter->begin(),iter->end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;
        }
#endif

    res.reserve(res.size() + fragRes.size());
    for (VECT_INT_VECT::const_iterator iter = fragRes.begin();
         iter != fragRes.end(); ++iter) {
      res.push_back(*iter);
    }
  }  // done with all fragments

  FindRings::storeRingsInfo(mol, res);

  // update the ring memberships of atoms and bonds in the molecule:
  // store the SSSR rings on the molecule as a property
  // we will ignore any existing SSSRs on the molecule - simply overwrite
  return rdcast<int>(res.size());
}

int symmetrizeSSSR(ROMol &mol, bool includeDativeBonds) {
  VECT_INT_VECT tmp;
  return symmetrizeSSSR(mol, tmp, includeDativeBonds);
};

int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds) {
  res.clear();
  VECT_INT_VECT sssrs;

  // FIX: need to set flag here the symmetrization has been done in order to
  // avoid repeating this work
  findSSSR(mol, sssrs, includeDativeBonds);

  // reinit as SYMM_SSSR
  mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);

  res.reserve(sssrs.size());
  for (const auto &r : sssrs) {
    res.emplace_back(r);
  }

  // now check if there are any extra rings on the molecule
  if (!mol.hasProp(common_properties::extraRings)) {
    // no extra rings nothing to be done
    return rdcast<int>(res.size());
  }
  const VECT_INT_VECT &extras =
      mol.getProp<VECT_INT_VECT>(common_properties::extraRings);

  // convert the rings to bond ids
  VECT_INT_VECT bondsssrs;
  RingUtils::convertToBonds(sssrs, bondsssrs, mol);

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

  // counts of each bond
  std::vector<int> bondCounts(mol.getNumBonds(), 0);
  for (const auto &r : bondsssrs) {
    for (const auto &b : r) {
      bondCounts[b] += 1;
    }
  }

  INT_VECT extraRing;
  for (auto &extraAtomRing : extras) {
    RingUtils::convertToBonds(extraAtomRing, extraRing, mol);
    for (auto &ring : bondsssrs) {
      if (ring.size() != extraRing.size()) {
        continue;
      }

      // If `ring` is the only provider of some bond, extraRing must also
      // provide that bond.
      bool shareBond = false;
      bool replacesAllUniqueBonds = true;
      for (auto &bondID : ring) {
        const int bondCount = bondCounts[bondID];
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
        res.push_back(extraAtomRing);
        FindRings::storeRingInfo(mol, extraAtomRing);
        break;
      }
    }
  }

  if (mol.hasProp(common_properties::extraRings)) {
    mol.clearProp(common_properties::extraRings);
  }
  return rdcast<int>(res.size());
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
