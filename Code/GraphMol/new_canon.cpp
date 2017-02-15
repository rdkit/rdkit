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
#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include <iostream>
#include <cassert>

namespace RDKit {
namespace Canon {
void CreateSinglePartition(unsigned int nAtoms, int *order, int *count,
                           canon_atom *atoms) {
  for (unsigned int i = 0; i < nAtoms; i++) {
    atoms[i].index = 0;
    order[i] = i;
    count[i] = 0;
  }
  count[0] = nAtoms;
}

void ActivatePartitions(unsigned int nAtoms, int *order, int *count,
                        int &activeset, int *next, int *changed) {
  unsigned int i, j;
  activeset = -1;
  for (i = 0; i < nAtoms; i++) next[i] = -2;

  i = 0;
  do {
    j = order[i];
    if (count[j] > 1) {
      next[j] = activeset;
      activeset = j;
      i += count[j];
    } else
      i++;
  } while (i < nAtoms);

  for (i = 0; i < nAtoms; i++) {
    j = order[i];
    int flag = 1;
    //#define SKIP_NODE_CHANGED_OPTIMIZATION 0
    //#ifndef SKIP_NODE_CHANGED_OPTIMIZATION
    //        if(count[j]){
    //          std::cout << "j " << j << std::endl;
    //          flag=(next[j]!=-2);
    //        }
    //#endif
    changed[j] = flag;
  }
}

void compareRingAtomsConcerningNumNeighbors(Canon::canon_atom *atoms,
                                            unsigned int nAtoms,
                                            const ROMol &mol) {
  RingInfo *ringInfo = mol.getRingInfo();
  if (!ringInfo->isInitialized()) {
    ringInfo->initialize();
  }
  for (unsigned idx = 0; idx < nAtoms; ++idx) {
    const Canon::canon_atom &a = atoms[idx];
    if (ringInfo->numAtomRings(a.atom->getIdx()) < 1) {
      continue;
    }
    std::deque<int> neighbors;
    neighbors.push_back(idx);
    unsigned currentRNIdx = 0;
    atoms[idx].neighborNum.reserve(1000);
    atoms[idx].revistedNeighbors.assign(1000, 0);
    char *visited = (char *)malloc(nAtoms * sizeof(char));
    memset(visited, 0, nAtoms * sizeof(char));
    unsigned count = 1;
    std::vector<int> nextLevelNbrs;
    char *lastLevelNbrs = (char *)malloc(nAtoms * sizeof(char));
    memset(lastLevelNbrs, 0, nAtoms * sizeof(char));
    char *currentLevelNbrs = (char *)malloc(nAtoms * sizeof(char));
    memset(currentLevelNbrs, 0, nAtoms * sizeof(char));
    int *revisitedNeighbors = (int *)malloc(nAtoms * sizeof(int));
    memset(revisitedNeighbors, 0, nAtoms * sizeof(int));
    while (!neighbors.empty()) {
      unsigned int numLevelNbrs = 0;
      nextLevelNbrs.resize(0);
      while (!neighbors.empty()) {
        int nidx = neighbors.front();
        neighbors.pop_front();
        const Canon::canon_atom &atom = atoms[nidx];
        if (ringInfo->numAtomRings(atom.atom->getIdx()) < 1) {
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
      memset(lastLevelNbrs, 0, nAtoms * sizeof(char));
      for (unsigned i = 0; i < nAtoms; ++i) {
        if (currentLevelNbrs[i]) {
          lastLevelNbrs[i] = 1;
        }
      }
      memset(currentLevelNbrs, 0, nAtoms * sizeof(char));
      std::vector<int> tmp;
      tmp.reserve(30);
      for (unsigned i = 0; i < nAtoms; ++i) {
        if (revisitedNeighbors[i] > 0) {
          tmp.push_back(revisitedNeighbors[i]);
        }
      }
      std::sort(tmp.begin(), tmp.end());
      tmp.push_back(-1);
      for (unsigned i = 0; i < tmp.size(); ++i) {
        if (currentRNIdx >= atoms[idx].revistedNeighbors.size()) {
          atoms[idx].revistedNeighbors.resize(
              atoms[idx].revistedNeighbors.size() + 1000);
        }
        atoms[idx].revistedNeighbors[currentRNIdx] = tmp[i];
        currentRNIdx++;
      }
      memset(revisitedNeighbors, 0, nAtoms * sizeof(int));

      atoms[idx].neighborNum.push_back(numLevelNbrs);
      atoms[idx].neighborNum.push_back(-1);

      neighbors.insert(neighbors.end(), nextLevelNbrs.begin(),
                       nextLevelNbrs.end());
      count++;
    }
    atoms[idx].revistedNeighbors.resize(currentRNIdx);

    free(visited);
    free(currentLevelNbrs);
    free(lastLevelNbrs);
    free(revisitedNeighbors);
  }
}

template <typename T>
void rankWithFunctor(T &ftor, bool breakTies, int *order,
                     bool useSpecial = false, bool useChirality = false,
                     const boost::dynamic_bitset<> *atomsInPlay = NULL,
                     const boost::dynamic_bitset<> *bondsInPlay = NULL) {
  const ROMol &mol = *ftor.dp_mol;
  canon_atom *atoms = ftor.dp_atoms;
  unsigned int nAts = mol.getNumAtoms();
  int *count = (int *)malloc(nAts * sizeof(int));
  int activeset;
  int *next = (int *)malloc(nAts * sizeof(int));
  int *changed = (int *)malloc(nAts * sizeof(int));
  char *touched = (char *)malloc(nAts * sizeof(char));
  memset(touched, 0, nAts * sizeof(char));
  memset(changed, 1, nAts * sizeof(int));
  CreateSinglePartition(nAts, order, count, atoms);
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
  ActivatePartitions(nAts, order, count, activeset, next, changed);
#ifdef VERBOSE_CANON
  std::cerr << "1a--------" << std::endl;
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    std::cerr << order[i] + 1 << " "
              << " index: " << atoms[order[i]].index
              << " count: " << count[order[i]] << std::endl;
  }
#endif
  RefinePartitions(mol, atoms, ftor, true, order, count, activeset, next,
                   changed, touched);
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
    ActivatePartitions(nAts, order, count, activeset, next, changed);
    RefinePartitions(mol, atoms, scftor, true, order, count, activeset, next,
                     changed, touched);
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
  if (!ringInfo->isInitialized()) {
    ringInfo->initialize();
  }
  for (unsigned i = 0; i < nAts; ++i) {
    if (ringInfo->numAtomRings(order[i])) {
      if (count[order[i]] > 2) {
        symRingAtoms += count[order[i]];
      }
      ringAtoms++;
      if (ringInfo->numAtomRings(order[i]) > 1 && count[order[i]] > 1) {
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
    ActivatePartitions(nAts, order, count, activeset, next, changed);
    RefinePartitions(mol, atoms, sftor, true, order, count, activeset, next,
                     changed, touched);
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
    BreakTies(mol, atoms, ftor, true, order, count, activeset, next, changed,
              touched);
#ifdef VERBOSE_CANON
    std::cerr << "3--------" << std::endl;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      std::cerr << order[i] + 1 << " "
                << " index: " << atoms[order[i]].index
                << " count: " << count[order[i]] << std::endl;
    }
#endif
  }
  free(count);
  free(next);
  free(touched);
  free(changed);
}

namespace {
bool hasRingNbr(const ROMol &mol, const Atom *at) {
  ROMol::ADJ_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomNeighbors(at);
  while (beg != end) {
    const ATOM_SPTR nbr = mol[*beg];
    ++beg;
    if ((nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
         nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
        nbr->hasProp(common_properties::_ringStereoAtoms)) {
      return true;
    }
  }
  return false;
}

void getNbrs(const ROMol &mol, const Atom *at, int *ids) {
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(at);
  unsigned int idx = 0;

  while (beg != end) {
    const BOND_SPTR bond = (mol)[*beg];
    ++beg;
    unsigned int nbrIdx = bond->getOtherAtomIdx(at->getIdx());
    ids[idx] = nbrIdx;
    ++idx;
  }
}

void getBonds(const ROMol &mol, const Atom *at, std::vector<bondholder> &nbrs,
              bool includeChirality) {
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(at);
  while (beg != end) {
    const BOND_SPTR bond = (mol)[*beg];
    ++beg;
    Bond::BondStereo stereo = Bond::STEREONONE;
    if (includeChirality) {
      stereo = bond->getStereo();
      if (stereo == Bond::STEREOANY) {
        stereo = Bond::STEREONONE;
      }
    }
    unsigned int idx = bond->getOtherAtomIdx(at->getIdx());
    Bond::BondType bt =
        bond->getIsAromatic() ? Bond::AROMATIC : bond->getBondType();
    bondholder bh(bondholder(bt, stereo, idx, 0));
    nbrs.push_back(bh);
  }
  std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
}

void getChiralBonds(const ROMol &mol, const Atom *at,
                    std::vector<bondholder> &nbrs) {
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = mol.getAtomBonds(at);
  while (beg != end) {
    const BOND_SPTR bond = (mol)[*beg];
    ++beg;
    unsigned int nbrIdx = bond->getOtherAtomIdx(at->getIdx());
    const Atom *nbr = mol.getAtomWithIdx(nbrIdx);
    unsigned int degreeNbr = nbr->getDegree();
    unsigned int nReps = 1;
    unsigned int stereo = 0;
    // FIX: Since the user can set the stereoatoms, the use of STEREOCIS and
    // STEREOTRANS here isn't actually canonical. In order for that to work we
    // would need to be able to canonicalize the CIS/TRANS assignment, which is
    // not currently being done
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
      // a special case for chiral phophorous compounds
      // (this was leading to incorrect assignment of
      // R/S labels ):
      nReps = 1;
      // general justification of this is:
      // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
      // Multiple-Bond Unsaturation and Aromaticity". It contains several
      // conventions of which convention (b) is the one applying here:
      // "(b) Contibutions by d orbitals to bonds of quadriligant atoms are
      // neglected."
      // FIX: this applies to more than just P
    } else {
      nReps =
          static_cast<unsigned int>(floor(2. * bond->getBondTypeAsDouble()));
    }
    unsigned int symclass =
        nbr->getAtomicNum() * ATNUM_CLASS_OFFSET + nbrIdx + 1;
    bondholder bh(bondholder(Bond::SINGLE, stereo, nbrIdx, symclass));
    std::vector<bondholder>::iterator iPos =
        std::lower_bound(nbrs.begin(), nbrs.end(), bh);
    nbrs.insert(iPos, nReps, bh);
  }
  std::reverse(nbrs.begin(), nbrs.end());

  if (!at->needsUpdatePropertyCache()) {
    for (unsigned int ii = 0; ii < at->getTotalNumHs(); ++ii) {
      nbrs.push_back(bondholder(Bond::SINGLE, Bond::STEREONONE,
                                ATNUM_CLASS_OFFSET, ATNUM_CLASS_OFFSET));
      nbrs.push_back(bondholder(Bond::SINGLE, Bond::STEREONONE,
                                ATNUM_CLASS_OFFSET, ATNUM_CLASS_OFFSET));
    }
  }
}

void basicInitCanonAtom(const ROMol &mol, Canon::canon_atom &atom,
                        const int &idx) {
  atom.atom = mol.getAtomWithIdx(idx);
  atom.index = idx;
  atom.p_symbol = NULL;
  atom.degree = atom.atom->getDegree();
  atom.nbrIds = (int *)malloc(atom.degree * sizeof(int));
  getNbrs(mol, atom.atom, atom.nbrIds);
}

void advancedInitCanonAtom(const ROMol &mol, Canon::canon_atom &atom,
                           const int &idx) {
  RDUNUSED_PARAM(idx);
  atom.totalNumHs = atom.atom->getTotalNumHs();
  atom.isRingStereoAtom =
      (atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
       atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
      atom.atom->hasProp(common_properties::_ringStereoAtoms);
  atom.hasRingNbr = hasRingNbr(mol, atom.atom);
}
}  // end anonymous namespace

void initCanonAtoms(const ROMol &mol, std::vector<Canon::canon_atom> &atoms,
                    bool includeChirality) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    basicInitCanonAtom(mol, atoms[i], i);
    advancedInitCanonAtom(mol, atoms[i], i);
    atoms[i].bonds.reserve(atoms[i].degree);
    getBonds(mol, atoms[i].atom, atoms[i].bonds, includeChirality);
  }
}

void initFragmentCanonAtoms(const ROMol &mol,
                            std::vector<Canon::canon_atom> &atoms,
                            bool includeChirality,
                            const std::vector<std::string> *atomSymbols) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    basicInitCanonAtom(mol, atoms[i], i);
    atoms[i].degree = 0;
    if (atomSymbols) {
      atoms[i].p_symbol = &(*atomSymbols)[i];
    } else {
      atoms[i].p_symbol = 0;
    }
    advancedInitCanonAtom(mol, atoms[i], i);
    atoms[i].bonds.reserve(atoms[i].degree);
    getBonds(mol, atoms[i].atom, atoms[i].bonds, includeChirality);
  }
}

void initChiralCanonAtoms(const ROMol &mol,
                          std::vector<Canon::canon_atom> &atoms) {
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    basicInitCanonAtom(mol, atoms[i], i);
    getChiralBonds(mol, atoms[i].atom, atoms[i].bonds);
  }
}

void freeCanonAtoms(std::vector<Canon::canon_atom> &atoms) {
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    if (atoms[i].nbrIds) {
      free(atoms[i].nbrIds);
      atoms[i].nbrIds = NULL;
    }
  }
}

void updateAtomNeighborIndex(canon_atom *atoms, std::vector<bondholder> &nbrs) {
  for (unsigned j = 0; j < nbrs.size(); ++j) {
    unsigned nbrIdx = nbrs[j].nbrIdx;
    unsigned newSymClass = atoms[nbrIdx].index;
    nbrs.at(j).nbrSymClass = newSymClass;
  }
  std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
}

void updateAtomNeighborNumSwaps(
    canon_atom *atoms, std::vector<bondholder> &nbrs, unsigned int atomIdx,
    std::vector<std::pair<unsigned int, unsigned int> > &result) {
  for (unsigned j = 0; j < nbrs.size(); ++j) {
    unsigned nbrIdx = nbrs[j].nbrIdx;

    if (atoms[nbrIdx].atom->getChiralTag() != 0) {
      std::vector<int> ref, probe;
      for (unsigned i = 0; i < atoms[nbrIdx].degree; ++i) {
        ref.push_back(atoms[nbrIdx].nbrIds[i]);
      }

      probe.push_back(atomIdx);
      for (unsigned i = 0; i < atoms[nbrIdx].bonds.size(); ++i) {
        if (atoms[nbrIdx].bonds.at(i).nbrIdx != atomIdx) {
          probe.push_back(atoms[nbrIdx].bonds.at(i).nbrIdx);
        }
      }

      int nSwaps = static_cast<int>(countSwapsToInterconvert(ref, probe));
      if (atoms[nbrIdx].atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
        if (nSwaps % 2) {
          result.push_back(std::make_pair(nbrs[j].nbrSymClass, 2));
        } else {
          result.push_back(std::make_pair(nbrs[j].nbrSymClass, 1));
        }
      } else if (atoms[nbrIdx].atom->getChiralTag() ==
                 Atom::CHI_TETRAHEDRAL_CCW) {
        if (nSwaps % 2) {
          result.push_back(std::make_pair(nbrs[j].nbrSymClass, 1));
        } else {
          result.push_back(std::make_pair(nbrs[j].nbrSymClass, 2));
        }
      }
    } else {
      result.push_back(std::make_pair(nbrs[j].nbrSymClass, 0));
    }
  }
  sort(result.begin(), result.end());
}

void rankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res,
                  bool breakTies, bool includeChirality, bool includeIsotopes) {
  std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
  initCanonAtoms(mol, atoms, includeChirality);
  AtomCompareFunctor ftor(&atoms.front(), mol);
  ftor.df_useIsotopes = includeIsotopes;
  ftor.df_useChirality = includeChirality;
  ftor.df_useChiralityRings = includeChirality;

  int *order = (int *)malloc(mol.getNumAtoms() * sizeof(int));
  rankWithFunctor(ftor, breakTies, order, true, includeChirality);

  res.resize(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res[order[i]] = atoms[order[i]].index;
  }
  free(order);
  freeCanonAtoms(atoms);
}  // end of rankMolAtoms()

void rankFragmentAtoms(const ROMol &mol, std::vector<unsigned int> &res,
                       const boost::dynamic_bitset<> &atomsInPlay,
                       const boost::dynamic_bitset<> &bondsInPlay,
                       const std::vector<std::string> *atomSymbols,
                       bool breakTies, bool includeChirality,
                       bool includeIsotopes) {
  PRECONDITION(atomsInPlay.size() == mol.getNumAtoms(), "bad atomsInPlay size");
  PRECONDITION(bondsInPlay.size() == mol.getNumBonds(), "bad bondsInPlay size");
  PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
               "bad atomSymbols size");

  std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
  initFragmentCanonAtoms(mol, atoms, includeChirality, atomSymbols);
  for (ROMol::ConstBondIterator bI = mol.beginBonds(); bI != mol.endBonds();
       ++bI) {
    if (!bondsInPlay[(*bI)->getIdx()]) continue;
    atoms[(*bI)->getBeginAtomIdx()].degree++;
    atoms[(*bI)->getEndAtomIdx()].degree++;
  }
  AtomCompareFunctor ftor(&atoms.front(), mol, &atomsInPlay, &bondsInPlay);
  ftor.df_useIsotopes = includeIsotopes;
  ftor.df_useChirality = includeChirality;

  int *order = (int *)malloc(mol.getNumAtoms() * sizeof(int));

  rankWithFunctor(ftor, breakTies, order, true, includeChirality, &atomsInPlay,
                  &bondsInPlay);

  res.resize(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res[order[i]] = atoms[order[i]].index;
  }
  free(order);
  freeCanonAtoms(atoms);
}  // end of rankFragmentAtoms()

void chiralRankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res) {
  std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
  initChiralCanonAtoms(mol, atoms);
  ChiralAtomCompareFunctor ftor(&atoms.front(), mol);

  int *order = (int *)malloc(mol.getNumAtoms() * sizeof(int));
  rankWithFunctor(ftor, false, order);

  res.resize(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    res[order[i]] = atoms[order[i]].index;
  }
  free(order);
  freeCanonAtoms(atoms);
}  // end of rankMolAtoms()
}  // end of Canon namespace
}  // end of RDKit namespace
