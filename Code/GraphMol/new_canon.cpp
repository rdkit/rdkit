//
//  Copyright (C) 2014 Greg Landrum
//  Adapted from pseudo-code from Roger Sayle
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "new_canon.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Atropisomers.h>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <cassert>
// #define VERBOSE_CANON 1

namespace RDKit {
namespace Canon {

namespace {
void flipIfNeeded(Bond::BondStereo &st1,
                  const canon_atom *const *controllingAtoms) {
  CHECK_INVARIANT(controllingAtoms[0], "missing controlling atom");
  CHECK_INVARIANT(controllingAtoms[2], "missing controlling atom");
  bool flip = false;
  if (controllingAtoms[1] &&
      controllingAtoms[1]->index > controllingAtoms[0]->index) {
    flip = !flip;
  }
  if (controllingAtoms[3] &&
      controllingAtoms[3]->index > controllingAtoms[2]->index) {
    flip = !flip;
  }
  if (flip) {
    if (st1 == Bond::BondStereo::STEREOCIS) {
      st1 = Bond::BondStereo::STEREOTRANS;
    } else if (st1 == Bond::BondStereo::STEREOTRANS) {
      st1 = Bond::BondStereo::STEREOCIS;
    }
  }
}
}  // namespace

int bondholder::compareStereo(const bondholder &o) const {
  auto st1 = stype;
  auto st2 = o.stype;
  if (st1 == Bond::BondStereo::STEREONONE) {
    if (st2 == Bond::BondStereo::STEREONONE) {
      return 0;
    } else {
      return -1;
    }
  }
  if (st2 == Bond::BondStereo::STEREONONE) {
    return 1;
  }
  if (st1 == Bond::BondStereo::STEREOANY) {
    if (st2 == Bond::BondStereo::STEREOANY) {
      return 0;
    } else {
      return -1;
    }
  }
  if (st2 == Bond::BondStereo::STEREOANY) {
    return 1;
  }
  // we have some kind of specified stereo on both bonds, work is required

  // if both have absolute stereo labels we can compare them directly
  if ((st1 == Bond::BondStereo::STEREOE || st1 == Bond::BondStereo::STEREOZ) &&
      (st2 == Bond::BondStereo::STEREOE || st2 == Bond::BondStereo::STEREOZ)) {
    if (st1 < st2) {
      return -1;
    } else if (st1 > st2) {
      return 1;
    }
    return 0;
  }

  // check to see if we need to flip the controlling atoms due to atom ranks
  flipIfNeeded(st1, controllingAtoms);
  flipIfNeeded(st2, o.controllingAtoms);
  if (st1 < st2) {
    return -1;
  } else if (st1 > st2) {
    return 1;
  }
  return 0;
}

void CreateSinglePartition(unsigned int nAtoms, int *order, int *count,
                           canon_atom *atoms) {
  PRECONDITION(order, "bad pointer");
  PRECONDITION(count, "bad pointer");
  PRECONDITION(atoms, "bad pointer");

  for (unsigned int i = 0; i < nAtoms; i++) {
    atoms[i].index = 0;
    order[i] = i;
    count[i] = 0;
  }
  count[0] = nAtoms;
}

void ActivatePartitions(unsigned int nAtoms, int *order, int *count,
                        int &activeset, int *next, int *changed) {
  PRECONDITION(order, "bad pointer");
  PRECONDITION(count, "bad pointer");
  PRECONDITION(next, "bad pointer");
  PRECONDITION(changed, "bad pointer");
  unsigned int i, j;
  activeset = -1;
  for (i = 0; i < nAtoms; i++) {
    next[i] = -2;
  }

  i = 0;
  do {
    j = order[i];
    if (count[j] > 1) {
      next[j] = activeset;
      activeset = j;
      i += count[j];
    } else {
      i++;
    }
  } while (i < nAtoms);

  for (i = 0; i < nAtoms; i++) {
    j = order[i];
    int flag = 1;
    // #define SKIP_NODE_CHANGED_OPTIMIZATION 0
    // #ifndef SKIP_NODE_CHANGED_OPTIMIZATION
    //         if(count[j]){
    //           std::cout << "j " << j << std::endl;
    //           flag=(next[j]!=-2);
    //         }
    // #endif
    changed[j] = flag;
  }
}

void compareRingAtomsConcerningNumNeighbors(Canon::canon_atom *atoms,
                                            unsigned int nAtoms,
                                            const ROMol &mol) {
  PRECONDITION(atoms, "bad pointer");
  RingInfo *ringInfo = mol.getRingInfo();

  auto visited = std::make_unique<char[]>(nAtoms);
  auto lastLevelNbrs = std::make_unique<char[]>(nAtoms);
  auto currentLevelNbrs = std::make_unique<char[]>(nAtoms);
  auto revisitedNeighbors = std::make_unique<int[]>(nAtoms);
  for (unsigned idx = 0; idx < nAtoms; ++idx) {
    const Canon::canon_atom &a = atoms[idx];
    if (!ringInfo->isInitialized() ||
        ringInfo->numAtomRings(a.atom->getIdx()) < 1) {
      continue;
    }
    std::deque<int> neighbors;
    neighbors.push_back(idx);
    unsigned currentRNIdx = 0;
    atoms[idx].neighborNum.reserve(1000);
    atoms[idx].revistedNeighbors.assign(1000, 0);

    memset(visited.get(), 0, nAtoms * sizeof(char));
    memset(lastLevelNbrs.get(), 0, nAtoms * sizeof(char));
    memset(currentLevelNbrs.get(), 0, nAtoms * sizeof(char));
    memset(revisitedNeighbors.get(), 0, nAtoms * sizeof(int));
    std::vector<int> nextLevelNbrs;
    while (!neighbors.empty()) {
      unsigned int numLevelNbrs = 0;
      nextLevelNbrs.resize(0);
      while (!neighbors.empty()) {
        int nidx = neighbors.front();
        neighbors.pop_front();
        const Canon::canon_atom &atom = atoms[nidx];
        if (!ringInfo->isInitialized() ||
            ringInfo->numAtomRings(atom.atom->getIdx()) < 1) {
          continue;
        }
        lastLevelNbrs[nidx] = 1;
        visited[nidx] = 1;
        for (unsigned int j = 0; j < atom.degree; j++) {
          int iidx = atom.nbrIds[j];
          if (!visited[iidx]) {
            currentLevelNbrs[iidx] = 1;
            numLevelNbrs++;
            visited[iidx] = 1;
            nextLevelNbrs.push_back(iidx);
          }
        }
      }
      for (unsigned i = 0; i < nAtoms; ++i) {
        if (currentLevelNbrs[i]) {
          const Canon::canon_atom &natom = atoms[i];
          for (unsigned int k = 0; k < natom.degree; k++) {
            int jidx = natom.nbrIds[k];
            if (currentLevelNbrs[jidx] || lastLevelNbrs[jidx]) {
              revisitedNeighbors[jidx] += 1;
            }
          }
        }
      }
      memset(lastLevelNbrs.get(), 0, nAtoms * sizeof(char));
      for (unsigned i = 0; i < nAtoms; ++i) {
        if (currentLevelNbrs[i]) {
          lastLevelNbrs[i] = 1;
        }
      }
      memset(currentLevelNbrs.get(), 0, nAtoms * sizeof(char));
      std::vector<int> tmp;
      tmp.reserve(30);
      for (unsigned i = 0; i < nAtoms; ++i) {
        if (revisitedNeighbors[i] > 0) {
          tmp.push_back(revisitedNeighbors[i]);
        }
      }
      std::sort(tmp.begin(), tmp.end());
      tmp.push_back(-1);
      for (int i : tmp) {
        if (currentRNIdx >= atoms[idx].revistedNeighbors.size()) {
          atoms[idx].revistedNeighbors.resize(
              atoms[idx].revistedNeighbors.size() + 1000);
        }
        atoms[idx].revistedNeighbors[currentRNIdx] = i;
        currentRNIdx++;
      }
      memset(revisitedNeighbors.get(), 0, nAtoms * sizeof(int));

      atoms[idx].neighborNum.push_back(numLevelNbrs);
      atoms[idx].neighborNum.push_back(-1);

      neighbors.insert(neighbors.end(), nextLevelNbrs.begin(),
                       nextLevelNbrs.end());
    }
    atoms[idx].revistedNeighbors.resize(currentRNIdx);
  }
}

namespace detail {
template <typename T>
void rankWithFunctor(T &ftor, bool breakTies, int *order, bool useSpecial,
                     bool useChirality,
                     const boost::dynamic_bitset<> *atomsInPlay,
                     const boost::dynamic_bitset<> *bondsInPlay) {
  PRECONDITION(order, "bad pointer");
  const ROMol &mol = *ftor.dp_mol;
  canon_atom *atoms = ftor.dp_atoms;
  unsigned int nAts = mol.getNumAtoms();

  //  auto order = std::make_unique<int[]>(mol.getNumAtoms());

  auto count = std::make_unique<int[]>(nAts);
  auto next = std::make_unique<int[]>(nAts);
  auto changed = std::make_unique<int[]>(nAts);
  memset(changed.get(), 1, nAts * sizeof(int));
  auto touched = std::make_unique<char[]>(nAts);
  memset(touched.get(), 0, nAts * sizeof(char));
  int activeset;
  CreateSinglePartition(nAts, order, count.get(), atoms);
// ActivatePartitions(nAts,order,count,activeset,next,changed);
// RefinePartitions(mol,atoms,ftor,false,order,count,activeset,next,changed,touched);
#ifdef VERBOSE_CANON
  std::cerr << "1--------" << std::endl;
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    std::cerr << order[i] + 1 << " "
              << " index: " << atoms[order[i]].index
              << " count: " << count[order[i]] << std::endl;
  }
#endif
  ftor.df_useNbrs = true;
  ActivatePartitions(nAts, order, count.get(), activeset, next.get(),
                     changed.get());
#ifdef VERBOSE_CANON
  std::cerr << "1a--------" << std::endl;
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    std::cerr << order[i] + 1 << " "
              << " index: " << atoms[order[i]].index
              << " count: " << count[order[i]] << std::endl;
  }
#endif
  RefinePartitions(mol, atoms, ftor, true, order, count.get(), activeset,
                   next.get(), changed.get(), touched.get());
#ifdef VERBOSE_CANON
  std::cerr << "2--------" << std::endl;
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    std::cerr << order[i] + 1 << " "
              << " index: " << atoms[order[i]].index
              << " count: " << count[order[i]] << std::endl;
  }
#endif
  bool ties = false;
  for (unsigned i = 0; i < nAts; ++i) {
    if (!count[i]) {
      ties = true;
    }
  }
  if (useChirality && ties) {
    SpecialChiralityAtomCompareFunctor scftor(atoms, mol, atomsInPlay,
                                              bondsInPlay);
    ActivatePartitions(nAts, order, count.get(), activeset, next.get(),
                       changed.get());
    RefinePartitions(mol, atoms, scftor, true, order, count.get(), activeset,
                     next.get(), changed.get(), touched.get());
#ifdef VERBOSE_CANON
    std::cerr << "2a--------" << std::endl;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      std::cerr << order[i] + 1 << " "
                << " index: " << atoms[order[i]].index
                << " count: " << count[order[i]] << std::endl;
    }
#endif
  }
  ties = false;
  unsigned symRingAtoms = 0;
  unsigned ringAtoms = 0;
  bool branchingRingAtom = false;
  RingInfo *ringInfo = mol.getRingInfo();
  for (unsigned i = 0; i < nAts; ++i) {
    if (ringInfo->isInitialized() && ringInfo->numAtomRings(order[i])) {
      if (count[order[i]] > 2) {
        symRingAtoms += count[order[i]];
      }
      ringAtoms++;
      if (ringInfo->isInitialized() && ringInfo->numAtomRings(order[i]) > 1 &&
          count[order[i]] > 1) {
        branchingRingAtom = true;
      }
    }
    if (!count[i]) {
      ties = true;
    }
  }
  //      std::cout << " " << ringAtoms << " "  << symRingAtoms << std::endl;
  if (useSpecial && ties && ringAtoms > 0 &&
      static_cast<float>(symRingAtoms) / ringAtoms > 0.5 && branchingRingAtom) {
    SpecialSymmetryAtomCompareFunctor sftor(atoms, mol, atomsInPlay,
                                            bondsInPlay);
    compareRingAtomsConcerningNumNeighbors(atoms, nAts, mol);
    ActivatePartitions(nAts, order, count.get(), activeset, next.get(),
                       changed.get());
    RefinePartitions(mol, atoms, sftor, true, order, count.get(), activeset,
                     next.get(), changed.get(), touched.get());
#ifdef VERBOSE_CANON
    std::cerr << "2b--------" << std::endl;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      std::cerr << order[i] + 1 << " "
                << " index: " << atoms[order[i]].index
                << " count: " << count[order[i]] << std::endl;
    }
#endif
  }
  if (breakTies) {
    BreakTies(mol, atoms, ftor, true, order, count.get(), activeset, next.get(),
              changed.get(), touched.get());
#ifdef VERBOSE_CANON
    std::cerr << "3--------" << std::endl;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      std::cerr << order[i] + 1 << " "
                << " index: " << atoms[order[i]].index
                << " count: " << count[order[i]] << std::endl;
    }
#endif
  }
}
}  // namespace detail
namespace {
bool hasRingNbr(const ROMol &mol, const Atom *at) {
  PRECONDITION(at, "bad pointer");
  for (const auto nbr : mol.atomNeighbors(at)) {
    if ((nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
         nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
        nbr->hasProp(common_properties::_ringStereoAtoms)) {
      return true;
    }
  }
  return false;
}

void getNbrs(const ROMol &mol, const Atom *at, int *ids) {
  PRECONDITION(at, "bad pointer");
  PRECONDITION(ids, "bad pointer");
  ROMol::ADJ_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomNeighbors(at);
  unsigned int idx = 0;
  while (beg != end) {
    ids[idx++] = static_cast<int>(*beg++);
  }
}

bondholder makeBondHolder(const Bond *bond, unsigned int otherIdx,
                          bool includeChirality,
                          const std::vector<Canon::canon_atom> &atoms) {
  PRECONDITION(bond, "bad pointer");
  Bond::BondStereo stereo = Bond::STEREONONE;
  if (includeChirality) {
    stereo = bond->getStereo();
  }
  Bond::BondType bt =
      bond->getIsAromatic() ? Bond::AROMATIC : bond->getBondType();
  bondholder res(bt, stereo, otherIdx, 0, bond->getIdx());
  if (includeChirality) {
    res.stype = bond->getStereo();
    if (res.stype == Bond::BondStereo::STEREOCIS ||
        res.stype == Bond::BondStereo::STEREOTRANS) {
      res.controllingAtoms[0] = &atoms[bond->getStereoAtoms()[0]];
      res.controllingAtoms[2] = &atoms[bond->getStereoAtoms()[1]];
      if (bond->getBeginAtom()->getDegree() > 2) {
        for (const auto nbr :
             bond->getOwningMol().atomNeighbors(bond->getBeginAtom())) {
          if (nbr->getIdx() != bond->getEndAtomIdx() &&
              nbr->getIdx() !=
                  static_cast<unsigned int>(bond->getStereoAtoms()[0])) {
            res.controllingAtoms[1] = &atoms[nbr->getIdx()];
          }
        }
      }
      if (bond->getEndAtom()->getDegree() > 2) {
        for (const auto nbr :
             bond->getOwningMol().atomNeighbors(bond->getEndAtom())) {
          if (nbr->getIdx() != bond->getBeginAtomIdx() &&
              nbr->getIdx() !=
                  static_cast<unsigned int>(bond->getStereoAtoms()[1])) {
            res.controllingAtoms[3] = &atoms[nbr->getIdx()];
          }
        }
      }
    }

    if (res.stype == Bond::BondStereo::STEREOATROPCCW ||
        res.stype == Bond::BondStereo::STEREOATROPCW) {
      Atropisomers::AtropAtomAndBondVec atropAtomAndBondVecs[2];
      CHECK_INVARIANT(Atropisomers::getAtropisomerAtomsAndBonds(
                          bond, atropAtomAndBondVecs, bond->getOwningMol()),
                      "Could not find atropisomer controlling atoms")

      res.controllingAtoms[0] =
          &atoms[atropAtomAndBondVecs[0]
                     .second[0]
                     ->getOtherAtom(atropAtomAndBondVecs[0].first)
                     ->getIdx()];
      res.controllingAtoms[2] =
          &atoms[atropAtomAndBondVecs[1]
                     .second[0]
                     ->getOtherAtom(atropAtomAndBondVecs[1].first)
                     ->getIdx()];
      if (atropAtomAndBondVecs[0].second.size() > 1) {
        res.controllingAtoms[1] =
            &atoms[atropAtomAndBondVecs[0]
                       .second[1]
                       ->getOtherAtom(atropAtomAndBondVecs[0].first)
                       ->getIdx()];
      }
      if (atropAtomAndBondVecs[1].second.size() > 1) {
        res.controllingAtoms[3] =
            &atoms[atropAtomAndBondVecs[1]
                       .second[1]
                       ->getOtherAtom(atropAtomAndBondVecs[1].first)
                       ->getIdx()];
      }
    }
  }
  return res;
}
void getBonds(const ROMol &mol, const Atom *at, std::vector<bondholder> &nbrs,
              bool includeChirality,
              const std::vector<Canon::canon_atom> &atoms) {
  PRECONDITION(at, "bad pointer");
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(at);
  while (beg != end) {
    const Bond *bond = (mol)[*beg];
    ++beg;
    nbrs.push_back(makeBondHolder(bond, bond->getOtherAtomIdx(at->getIdx()),
                                  includeChirality, atoms));
  }
  std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
}

void getChiralBonds(const ROMol &mol, const Atom *at,
                    std::vector<bondholder> &nbrs) {
  PRECONDITION(at, "bad pointer");
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(at);
  while (beg != end) {
    const Bond *bond = (mol)[*beg];
    ++beg;
    unsigned int nbrIdx = bond->getOtherAtomIdx(at->getIdx());
    const Atom *nbr = mol.getAtomWithIdx(nbrIdx);
    unsigned int degreeNbr = nbr->getDegree();
    unsigned int nReps = 1;
    unsigned int stereo = 0;
    // FIX: Since the user can set the stereoatoms, the use of STEREOCIS and
    // STEREOTRANS here isn't actually canonical. In order for that to work we
    // would need to be able to canonicalize the CIS/TRANS assignment, which
    // is not currently being done
    switch (bond->getStereo()) {
      case Bond::STEREOZ:
      case Bond::STEREOCIS:
        stereo = 1;
        break;
      case Bond::STEREOE:
      case Bond::STEREOTRANS:
        stereo = 2;
        break;
      default:
        stereo = 0;
    }
    if (bond->getBondType() == Bond::DOUBLE && nbr->getAtomicNum() == 15 &&
        (degreeNbr == 4 || degreeNbr == 3)) {
      // a special case for chiral phosphorous compounds
      // (this was leading to incorrect assignment of
      // R/S labels ):
      nReps = 1;
      // general justification of this is:
      // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
      // Multiple-Bond Unsaturation and Aromaticity". It contains several
      // conventions of which convention (b) is the one applying here:
      // "(b) Contributions by d orbitals to bonds of quadriligant atoms are
      // neglected."
      // FIX: this applies to more than just P
    } else {
      nReps =
          static_cast<unsigned int>(floor(2. * bond->getBondTypeAsDouble()));
    }
    unsigned int symclass =
        nbr->getAtomicNum() * ATNUM_CLASS_OFFSET + nbrIdx + 1;
    bondholder bh(
        bondholder(Bond::SINGLE, stereo, nbrIdx, symclass, bond->getIdx()));
    auto iPos = std::lower_bound(nbrs.begin(), nbrs.end(), bh);
    nbrs.insert(iPos, nReps, bh);
  }
  std::reverse(nbrs.begin(), nbrs.end());

  if (!at->needsUpdatePropertyCache()) {
    for (unsigned int ii = 0; ii < at->getTotalNumHs(); ++ii) {
      nbrs.emplace_back(Bond::SINGLE, Bond::STEREONONE, ATNUM_CLASS_OFFSET,
                        ATNUM_CLASS_OFFSET, 0);
      nbrs.emplace_back(Bond::SINGLE, Bond::STEREONONE, ATNUM_CLASS_OFFSET,
                        ATNUM_CLASS_OFFSET, 0);
    }
  }
}

void basicInitCanonAtom(const ROMol &mol, Canon::canon_atom &atom,
                        const int &idx) {
  atom.atom = mol.getAtomWithIdx(idx);
  atom.index = idx;
  atom.p_symbol = nullptr;
  atom.degree = atom.atom->getDegree();
  atom.nbrIds = std::make_unique<int[]>(atom.degree);
  getNbrs(mol, atom.atom, atom.nbrIds.get());
}

void advancedInitCanonAtom(const ROMol &mol, Canon::canon_atom &atom,
                           const int &) {
  atom.totalNumHs = atom.atom->getTotalNumHs();
  atom.isRingStereoAtom =
      (atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
       atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
      atom.atom->hasProp(common_properties::_ringStereoAtoms);
  atom.hasRingNbr = hasRingNbr(mol, atom.atom);
}
}  // end anonymous namespace

void initCanonAtoms(const ROMol &mol, std::vector<Canon::canon_atom> &atoms,
                    bool includeChirality, bool includeStereoGroups) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    basicInitCanonAtom(mol, atoms[i], i);
    advancedInitCanonAtom(mol, atoms[i], i);
    atoms[i].bonds.reserve(atoms[i].degree);
    getBonds(mol, atoms[i].atom, atoms[i].bonds, includeChirality, atoms);
  }
  if (includeChirality && includeStereoGroups) {
    unsigned int sgidx = 1;
    for (auto &sg : mol.getStereoGroups()) {
      for (auto atom : sg.getAtoms()) {
        atoms[atom->getIdx()].whichStereoGroup = sgidx;
        atoms[atom->getIdx()].typeOfStereoGroup = sg.getGroupType();
      }
      ++sgidx;
    }
  }
}
namespace detail {

void initFragmentCanonAtoms(const ROMol &mol,
                            std::vector<Canon::canon_atom> &atoms,
                            bool includeChirality,
                            const std::vector<std::string> *atomSymbols,
                            const std::vector<std::string> *bondSymbols,
                            const boost::dynamic_bitset<> &atomsInPlay,
                            const boost::dynamic_bitset<> &bondsInPlay,
                            bool needsInit) {
  needsInit = true;
  PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
               "bad atom symbols");
  PRECONDITION(!bondSymbols || bondSymbols->size() == mol.getNumBonds(),
               "bad bond symbols");
  // start by initializing the atoms
  for (const auto atom : mol.atoms()) {
    auto i = atom->getIdx();
    auto &atomsi = atoms[i];
    atomsi.atom = atom;
    atomsi.index = i;
    // we don't care about overall degree, so we start that at zero, and
    // then count the degree in the fragment itself below.
    atomsi.degree = 0;
    if (atomsInPlay[i]) {
      if (atomSymbols) {
        atomsi.p_symbol = &(*atomSymbols)[i];
      } else {
        atomsi.p_symbol = nullptr;
      }
      if (needsInit) {
        atomsi.nbrIds = std::make_unique<int[]>(atom->getDegree());
        advancedInitCanonAtom(mol, atomsi, i);
        atomsi.bonds.reserve(4);
      }
    }
  }

  // now deal with the bonds in the fragment.
  // these set the atomic degrees and update the neighbor lists
  if (needsInit) {
    for (const auto bond : mol.bonds()) {
      if (!bondsInPlay[bond->getIdx()] ||
          !atomsInPlay[bond->getBeginAtomIdx()] ||
          !atomsInPlay[bond->getEndAtomIdx()]) {
        continue;
      }
      Canon::canon_atom &begAt = atoms[bond->getBeginAtomIdx()];
      Canon::canon_atom &endAt = atoms[bond->getEndAtomIdx()];
      begAt.nbrIds[begAt.degree++] = bond->getEndAtomIdx();
      endAt.nbrIds[endAt.degree++] = bond->getBeginAtomIdx();
      begAt.bonds.push_back(
          makeBondHolder(bond, bond->getEndAtomIdx(), includeChirality, atoms));
      endAt.bonds.push_back(makeBondHolder(bond, bond->getBeginAtomIdx(),
                                           includeChirality, atoms));
      if (bondSymbols) {
        begAt.bonds.back().p_symbol = &(*bondSymbols)[bond->getIdx()];
        endAt.bonds.back().p_symbol = &(*bondSymbols)[bond->getIdx()];
      }
    }
  } else {
    if (bondSymbols) {
      for (auto &atom : atoms) {
        for (auto &bond : atom.bonds) {
          bond.p_symbol = &(*bondSymbols)[bond.bondIdx];
        }
      }
    }
  }

  // and now we can do the last bit for each atom
  for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    if (!atomsInPlay[i]) {
      continue;
    }
    auto &atomsi = atoms[i];
    if (needsInit) {
      // this is the fix for github #1567: we let the atom's degree
      // in the original molecule influence its degree in the fragment
      atomsi.totalNumHs += (mol.getAtomWithIdx(i)->getDegree() - atomsi.degree);
    }

    // and sort our list of neighboring bonds
    std::sort(atomsi.bonds.begin(), atomsi.bonds.end(), bondholder::greater);
  }
}

void initChiralCanonAtoms(const ROMol &mol,
                          std::vector<Canon::canon_atom> &atoms) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    basicInitCanonAtom(mol, atoms[i], i);
    getChiralBonds(mol, atoms[i].atom, atoms[i].bonds);
  }
}

}  // namespace detail
void updateAtomNeighborIndex(canon_atom *atoms, std::vector<bondholder> &nbrs) {
  PRECONDITION(atoms, "bad pointer");
  for (auto &nbr : nbrs) {
    unsigned nbrIdx = nbr.nbrIdx;
    unsigned newSymClass = atoms[nbrIdx].index;
    nbr.nbrSymClass = newSymClass;
  }
  std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
}

// This routine calculates the number of swaps that would be required to
// determine what the smiles chirality value would be for a given chiral atom
// given that the atom is visited first from the atom of interest.
// This is used to determine which of two atoms has priority based on the
// neighbor's chirality
//
// If the chiral neighbor has two equivlent (at least so far) neighbors that are
// not the atom of interest, it cannot be used to determine the priority of the
// atom of interest.  For this reason, we keep track of the number of neighbors
// that have the same priority so far.  If any two are the same, we do NOT use
// that neighbor to determine the priority of the atom of interest.

void updateAtomNeighborNumSwaps(
    canon_atom *atoms, std::vector<bondholder> &nbrs, unsigned int atomIdx,
    std::vector<std::pair<unsigned int, unsigned int>> &result) {
  bool isRingAtom = queryIsAtomInRing(atoms[atomIdx].atom);
  for (auto &nbr : nbrs) {
    unsigned nbrIdx = nbr.nbrIdx;

    std::list<unsigned int> neighborsSeen;
    bool tooManySimilarNbrs = false;
    if (isRingAtom && atoms[nbrIdx].atom->getChiralTag() != 0) {
      std::vector<int> ref, probe;
      for (unsigned i = 0; i < atoms[nbrIdx].degree; ++i) {
        auto nbrNbrId =
            atoms[nbrIdx].nbrIds[i];  // id of the neighbor's neighbor
        ref.push_back(nbrNbrId);
        if ((int)atomIdx != nbrNbrId) {
          if ((std::find(neighborsSeen.begin(), neighborsSeen.end(),
                         atoms[nbrNbrId].index) != neighborsSeen.end())) {
            tooManySimilarNbrs = true;
          } else {
            neighborsSeen.push_back(atoms[nbrNbrId].index);
          }
        }
      }

      probe.push_back(atomIdx);
      for (auto &bond : atoms[nbrIdx].bonds) {
        if (bond.nbrIdx != atomIdx) {
          probe.push_back(bond.nbrIdx);
        }
      }

      if (tooManySimilarNbrs) {
        result.emplace_back(nbr.nbrSymClass, 0);
      } else {
        int nSwaps = static_cast<int>(countSwapsToInterconvert(ref, probe));
        if (atoms[nbrIdx].atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
          if (nSwaps % 2) {
            result.emplace_back(nbr.nbrSymClass, 2);
          } else {
            result.emplace_back(nbr.nbrSymClass, 1);
          }
        } else if (atoms[nbrIdx].atom->getChiralTag() ==
                   Atom::CHI_TETRAHEDRAL_CCW) {
          if (nSwaps % 2) {
            result.emplace_back(nbr.nbrSymClass, 1);
          } else {
            result.emplace_back(nbr.nbrSymClass, 2);
          }
        }
      }
    } else {
      result.emplace_back(nbr.nbrSymClass, 0);
    }
  }
  sort(result.begin(), result.end());
}

void rankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res,
                  bool breakTies, bool includeChirality, bool includeIsotopes,
                  bool includeAtomMaps, bool includeChiralPresence,
                  bool includeStereoGroups) {
  if (!mol.getNumAtoms()) {
    return;
  }

  bool clearRings = false;
  if (!mol.getRingInfo()->isFindFastOrBetter()) {
    MolOps::fastFindRings(mol);
    clearRings = true;
  }
  res.resize(mol.getNumAtoms());

  std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
  initCanonAtoms(mol, atoms, includeChirality, includeStereoGroups);
  AtomCompareFunctor ftor(&atoms.front(), mol);
  ftor.df_useIsotopes = includeIsotopes;
  ftor.df_useChirality = includeChirality;
  ftor.df_useChiralityRings = includeChirality;
  ftor.df_useAtomMaps = includeAtomMaps;
  ftor.df_useChiralPresence = includeChiralPresence;

  auto order = std::make_unique<int[]>(mol.getNumAtoms());
  detail::rankWithFunctor(ftor, breakTies, order.get(), true, includeChirality);

  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res[order[i]] = atoms[order[i]].index;
  }

  if (clearRings) {
    mol.getRingInfo()->reset();
  }
}  // end of rankMolAtoms()

void rankFragmentAtoms(const ROMol &mol, std::vector<unsigned int> &res,
                       const boost::dynamic_bitset<> &atomsInPlay,
                       const boost::dynamic_bitset<> &bondsInPlay,
                       const std::vector<std::string> *atomSymbols,
                       const std::vector<std::string> *bondSymbols,
                       bool breakTies, bool includeChirality,
                       bool includeIsotopes, bool includeAtomMaps,
                       bool includeChiralPresence) {
  PRECONDITION(atomsInPlay.size() == mol.getNumAtoms(), "bad atomsInPlay size");
  PRECONDITION(bondsInPlay.size() == mol.getNumBonds(), "bad bondsInPlay size");
  PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
               "bad atomSymbols size");
  PRECONDITION(!bondSymbols || bondSymbols->size() == mol.getNumBonds(),
               "bad bondSymbols size");

  if (!mol.getNumAtoms()) {
    return;
  }

  bool clearRings = false;
  if (!mol.getRingInfo()->isFindFastOrBetter()) {
    MolOps::fastFindRings(mol);
    clearRings = true;
  }
  res.resize(mol.getNumAtoms());

  std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
  detail::initFragmentCanonAtoms(mol, atoms, includeChirality, atomSymbols,
                                 bondSymbols, atomsInPlay, bondsInPlay, true);

  AtomCompareFunctor ftor(&atoms.front(), mol, &atomsInPlay, &bondsInPlay);
  ftor.df_useIsotopes = includeIsotopes;
  ftor.df_useChirality = includeChirality;
  ftor.df_useAtomMaps = includeAtomMaps;
  ftor.df_useChiralityRings = includeChirality;
  ftor.df_useChiralPresence = includeChiralPresence;

  auto order = std::make_unique<int[]>(mol.getNumAtoms());
  detail::rankWithFunctor(ftor, breakTies, order.get(), true, includeChirality,
                          &atomsInPlay, &bondsInPlay);

  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res[order[i]] = atoms[order[i]].index;
  }

  if (clearRings) {
    mol.getRingInfo()->reset();
  }
}  // end of rankFragmentAtoms()

void chiralRankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res) {
  if (!mol.getNumAtoms()) {
    return;
  }

  bool clearRings = false;
  if (!mol.getRingInfo()->isFindFastOrBetter()) {
    MolOps::fastFindRings(mol);
    clearRings = true;
  }
  res.resize(mol.getNumAtoms());

  std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
  detail::initChiralCanonAtoms(mol, atoms);
  ChiralAtomCompareFunctor ftor(&atoms.front(), mol);

  auto order = std::make_unique<int[]>(mol.getNumAtoms());
  detail::rankWithFunctor(ftor, false, order.get());

  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res[order[i]] = atoms[order[i]].index;
  }

  if (clearRings) {
    mol.getRingInfo()->reset();
  }
}
}  // namespace Canon
}  // namespace RDKit
