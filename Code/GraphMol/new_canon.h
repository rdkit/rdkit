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

#include <RDGeneral/hanoiSort.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RingInfo.h>
#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/dynamic_bitset.hpp>
#include <cstring>
#include <iostream>
#include <cassert>
#include <cstring>
#include <vector>

//#define VERBOSE_CANON 1

namespace RDKit {
namespace Canon {

struct bondholder {
  Bond::BondType bondType;
  unsigned int bondStereo;
  unsigned int nbrSymClass;
  unsigned int nbrIdx;
  bondholder()
      : bondType(Bond::UNSPECIFIED),
        bondStereo(static_cast<unsigned int>(Bond::STEREONONE)),
        nbrSymClass(0),
        nbrIdx(0){};
  bondholder(Bond::BondType bt, Bond::BondStereo bs, unsigned int ni,
             unsigned int nsc)
      : bondType(bt),
        bondStereo(static_cast<unsigned int>(bs)),
        nbrSymClass(nsc),
        nbrIdx(ni){};
  bondholder(Bond::BondType bt, unsigned int bs, unsigned int ni,
             unsigned int nsc)
      : bondType(bt), bondStereo(bs), nbrSymClass(nsc), nbrIdx(ni){};
  bool operator<(const bondholder &o) const {
    if (bondType != o.bondType) return bondType < o.bondType;
    if (bondStereo != o.bondStereo) return bondStereo < o.bondStereo;
    return nbrSymClass < o.nbrSymClass;
  }
  static bool greater(const bondholder &lhs, const bondholder &rhs) {
    if (lhs.bondType != rhs.bondType) return lhs.bondType > rhs.bondType;
    if (lhs.bondStereo != rhs.bondStereo)
      return lhs.bondStereo > rhs.bondStereo;
    return lhs.nbrSymClass > rhs.nbrSymClass;
  }

  static int compare(const bondholder &x, const bondholder &y,
                     unsigned int div = 1) {
    if (x.bondType < y.bondType)
      return -1;
    else if (x.bondType > y.bondType)
      return 1;
    if (x.bondStereo < y.bondStereo)
      return -1;
    else if (x.bondStereo > y.bondStereo)
      return 1;
    return x.nbrSymClass / div - y.nbrSymClass / div;
  }
};

struct canon_atom {
  const Atom *atom;
  int index;
  unsigned int degree;
  unsigned int totalNumHs;
  bool hasRingNbr;
  bool isRingStereoAtom;
  int *nbrIds;
  const std::string *p_symbol;  // if provided, this is used to order atoms
  std::vector<int> neighborNum;
  std::vector<int> revistedNeighbors;
  std::vector<bondholder> bonds;

  canon_atom()
      : atom(NULL),
        index(-1),
        degree(0),
        totalNumHs(0),
        hasRingNbr(false),
        isRingStereoAtom(false),
        nbrIds(NULL),
        p_symbol(NULL){};
};

void updateAtomNeighborIndex(canon_atom *atoms, std::vector<bondholder> &nbrs);

void updateAtomNeighborNumSwaps(
    canon_atom *atoms, std::vector<bondholder> &nbrs, unsigned int atomIdx,
    std::vector<std::pair<unsigned int, unsigned int> > &result);

/*
 * Different types of atom compare functions:
 *
 * - SpecialChiralityAtomCompareFunctor: Allows canonizing molecules exhibiting
 *dependent chirality
 * - SpecialSymmetryAtomCompareFunctor: Very specialized, allows canonizing
 *highly symmetrical graphs/molecules
 * - AtomCompareFunctor: Basic atom compare function which also allows to
 *include neighbors within the ranking
 */

class SpecialChiralityAtomCompareFunctor {
 public:
  Canon::canon_atom *dp_atoms;
  const ROMol *dp_mol;
  const boost::dynamic_bitset<> *dp_atomsInPlay, *dp_bondsInPlay;

  SpecialChiralityAtomCompareFunctor()
      : dp_atoms(NULL),
        dp_mol(NULL),
        dp_atomsInPlay(NULL),
        dp_bondsInPlay(NULL){};
  SpecialChiralityAtomCompareFunctor(
      Canon::canon_atom *atoms, const ROMol &m,
      const boost::dynamic_bitset<> *atomsInPlay = NULL,
      const boost::dynamic_bitset<> *bondsInPlay = NULL)
      : dp_atoms(atoms),
        dp_mol(&m),
        dp_atomsInPlay(atomsInPlay),
        dp_bondsInPlay(bondsInPlay){};
  int operator()(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(i != j, "bad call");
    if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
      return 0;
    }

    if ((dp_atomsInPlay && (*dp_atomsInPlay)[i]) || !dp_atomsInPlay) {
      updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
    }
    if ((dp_atomsInPlay && (*dp_atomsInPlay)[j]) || !dp_atomsInPlay) {
      updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
    }
    for (unsigned int ii = 0;
         ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size(); ++ii) {
      int cmp =
          bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
      if (cmp) return cmp;
    }

    std::vector<std::pair<unsigned int, unsigned int> > swapsi;
    std::vector<std::pair<unsigned int, unsigned int> > swapsj;
    if ((dp_atomsInPlay && (*dp_atomsInPlay)[i]) || !dp_atomsInPlay) {
      updateAtomNeighborNumSwaps(dp_atoms, dp_atoms[i].bonds, i, swapsi);
    }
    if ((dp_atomsInPlay && (*dp_atomsInPlay)[j]) || !dp_atomsInPlay) {
      updateAtomNeighborNumSwaps(dp_atoms, dp_atoms[j].bonds, j, swapsj);
    }
    for (unsigned int ii = 0; ii < swapsi.size() && ii < swapsj.size(); ++ii) {
      int cmp = swapsi[ii].second - swapsj[ii].second;
      if (cmp) return cmp;
    }
    return 0;
  }
};

class SpecialSymmetryAtomCompareFunctor {
 public:
  Canon::canon_atom *dp_atoms;
  const ROMol *dp_mol;
  const boost::dynamic_bitset<> *dp_atomsInPlay, *dp_bondsInPlay;

  SpecialSymmetryAtomCompareFunctor()
      : dp_atoms(NULL),
        dp_mol(NULL),
        dp_atomsInPlay(NULL),
        dp_bondsInPlay(NULL){};
  SpecialSymmetryAtomCompareFunctor(
      Canon::canon_atom *atoms, const ROMol &m,
      const boost::dynamic_bitset<> *atomsInPlay = NULL,
      const boost::dynamic_bitset<> *bondsInPlay = NULL)
      : dp_atoms(atoms),
        dp_mol(&m),
        dp_atomsInPlay(atomsInPlay),
        dp_bondsInPlay(bondsInPlay){};
  int operator()(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(i != j, "bad call");
    if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
      return 0;
    }

    if (dp_atoms[i].neighborNum < dp_atoms[j].neighborNum) {
      return -1;
    } else if (dp_atoms[i].neighborNum > dp_atoms[j].neighborNum) {
      return 1;
    }

    if (dp_atoms[i].revistedNeighbors < dp_atoms[j].revistedNeighbors) {
      return -1;
    } else if (dp_atoms[i].revistedNeighbors > dp_atoms[j].revistedNeighbors) {
      return 1;
    }

    if ((dp_atomsInPlay && (*dp_atomsInPlay)[i]) || !dp_atomsInPlay) {
      updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
    }
    if ((dp_atomsInPlay && (*dp_atomsInPlay)[j]) || !dp_atomsInPlay) {
      updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
    }
    for (unsigned int ii = 0;
         ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size(); ++ii) {
      int cmp =
          bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
      if (cmp) return cmp;
    }

    if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
      return -1;
    } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
      return 1;
    }
    return 0;
  }
};

class AtomCompareFunctor {
  unsigned int getAtomRingNbrCode(unsigned int i) const {
    if (!dp_atoms[i].hasRingNbr) return 0;

    int *nbrs = dp_atoms[i].nbrIds;
    unsigned int code = 0;
    for (unsigned j = 0; j < dp_atoms[i].degree; ++j) {
      if (dp_atoms[nbrs[j]].isRingStereoAtom) {
        code += dp_atoms[nbrs[j]].index * 10000 + 1;  // j;
      }
    }
    return code;
  }

  int basecomp(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    unsigned int ivi, ivj;

    // always start with the current class:
    ivi = dp_atoms[i].index;
    ivj = dp_atoms[j].index;
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // use the atom-mapping numbers if they were assigned
    /* boost::any_cast FILED here:
            std::string molAtomMapNumber_i="";
            std::string molAtomMapNumber_j="";
    */
    int molAtomMapNumber_i = 0;
    int molAtomMapNumber_j = 0;
    dp_atoms[i].atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                       molAtomMapNumber_i);
    dp_atoms[j].atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                       molAtomMapNumber_j);
    if (molAtomMapNumber_i < molAtomMapNumber_j)
      return -1;
    else if (molAtomMapNumber_i > molAtomMapNumber_j)
      return 1;

    // start by comparing degree
    ivi = dp_atoms[i].degree;
    ivj = dp_atoms[j].degree;
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    if (dp_atoms[i].p_symbol && dp_atoms[j].p_symbol) {
      if (*(dp_atoms[i].p_symbol) < *(dp_atoms[j].p_symbol))
        return -1;
      else if (*(dp_atoms[i].p_symbol) > *(dp_atoms[j].p_symbol))
        return 1;
      else
        return 0;
    }
    // move onto atomic number
    ivi = dp_atoms[i].atom->getAtomicNum();
    ivj = dp_atoms[j].atom->getAtomicNum();
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // isotopes if we're using them
    if (df_useIsotopes) {
      ivi = dp_atoms[i].atom->getIsotope();
      ivj = dp_atoms[j].atom->getIsotope();
      if (ivi < ivj)
        return -1;
      else if (ivi > ivj)
        return 1;
    }

    // nHs
    ivi = dp_atoms[i].totalNumHs;
    ivj = dp_atoms[j].totalNumHs;
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // charge
    ivi = dp_atoms[i].atom->getFormalCharge();
    ivj = dp_atoms[j].atom->getFormalCharge();
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // chirality if we're using it
    if (df_useChirality) {
      // first atom stereochem:
      ivi = 0;
      ivj = 0;
      std::string cipCode;
      if (dp_atoms[i].atom->getPropIfPresent(common_properties::_CIPCode,
                                             cipCode)) {
        ivi = cipCode == "R" ? 2 : 1;
      }
      if (dp_atoms[j].atom->getPropIfPresent(common_properties::_CIPCode,
                                             cipCode)) {
        ivj = cipCode == "R" ? 2 : 1;
      }
      if (ivi < ivj)
        return -1;
      else if (ivi > ivj)
        return 1;

      // can't actually use values here, because they are arbitrary
      ivi = dp_atoms[i].atom->getChiralTag() != 0;
      ivj = dp_atoms[j].atom->getChiralTag() != 0;
      if (ivi < ivj)
        return -1;
      else if (ivi > ivj)
        return 1;
    }
    if (df_useChiralityRings) {
      // ring stereochemistry
      ivi = getAtomRingNbrCode(i);
      ivj = getAtomRingNbrCode(j);
      if (ivi < ivj)
        return -1;
      else if (ivi > ivj)
        return 1;
      // bond stereo is taken care of in the neighborhood comparison
    }
    return 0;
  }

 public:
  Canon::canon_atom *dp_atoms;
  const ROMol *dp_mol;
  const boost::dynamic_bitset<> *dp_atomsInPlay, *dp_bondsInPlay;
  bool df_useNbrs;
  bool df_useIsotopes;
  bool df_useChirality;
  bool df_useChiralityRings;

  AtomCompareFunctor()
      : dp_atoms(NULL),
        dp_mol(NULL),
        dp_atomsInPlay(NULL),
        dp_bondsInPlay(NULL),
        df_useNbrs(false),
        df_useIsotopes(true),
        df_useChirality(true),
        df_useChiralityRings(true){};
  AtomCompareFunctor(Canon::canon_atom *atoms, const ROMol &m,
                     const boost::dynamic_bitset<> *atomsInPlay = NULL,
                     const boost::dynamic_bitset<> *bondsInPlay = NULL)
      : dp_atoms(atoms),
        dp_mol(&m),
        dp_atomsInPlay(atomsInPlay),
        dp_bondsInPlay(bondsInPlay),
        df_useNbrs(false),
        df_useIsotopes(true),
        df_useChirality(true),
        df_useChiralityRings(true){};
  int operator()(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(i != j, "bad call");
    if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
      return 0;
    }

    int v = basecomp(i, j);
    if (v) {
      return v;
    }

    if (df_useNbrs) {
      if ((dp_atomsInPlay && (*dp_atomsInPlay)[i]) || !dp_atomsInPlay) {
        updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
      }
      if ((dp_atomsInPlay && (*dp_atomsInPlay)[j]) || !dp_atomsInPlay) {
        updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
      }

      for (unsigned int ii = 0;
           ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size();
           ++ii) {
        int cmp =
            bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
        if (cmp) return cmp;
      }

      if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
        return -1;
      } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
        return 1;
      }
    }
    return 0;
  }
};

/*
 * A compare function to discriminate chiral atoms, similar to the CIP rules.
 * This functionality is currently not used.
 */

const unsigned int ATNUM_CLASS_OFFSET = 10000;
class ChiralAtomCompareFunctor {
  void getAtomNeighborhood(std::vector<bondholder> &nbrs) const {
    for (unsigned j = 0; j < nbrs.size(); ++j) {
      unsigned int nbrIdx = nbrs[j].nbrIdx;
      if (nbrIdx == ATNUM_CLASS_OFFSET) {
        // Ignore the Hs
        continue;
      }
      const Atom *nbr = dp_atoms[nbrIdx].atom;
      nbrs[j].nbrSymClass =
          nbr->getAtomicNum() * ATNUM_CLASS_OFFSET + dp_atoms[nbrIdx].index + 1;
    }
    std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
    // FIX: don't want to be doing this long-term
  }

  int basecomp(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    unsigned int ivi, ivj;

    // always start with the current class:
    ivi = dp_atoms[i].index;
    ivj = dp_atoms[j].index;
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // move onto atomic number
    ivi = dp_atoms[i].atom->getAtomicNum();
    ivj = dp_atoms[j].atom->getAtomicNum();
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // isotopes:
    ivi = dp_atoms[i].atom->getIsotope();
    ivj = dp_atoms[j].atom->getIsotope();
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // atom stereochem:
    ivi = 0;
    ivj = 0;
    std::string cipCode;
    if (dp_atoms[i].atom->getPropIfPresent(common_properties::_CIPCode,
                                           cipCode)) {
      ivi = cipCode == "R" ? 2 : 1;
    }
    if (dp_atoms[j].atom->getPropIfPresent(common_properties::_CIPCode,
                                           cipCode)) {
      ivj = cipCode == "R" ? 2 : 1;
    }
    if (ivi < ivj)
      return -1;
    else if (ivi > ivj)
      return 1;

    // bond stereo is taken care of in the neighborhood comparison
    return 0;
  }

 public:
  Canon::canon_atom *dp_atoms;
  const ROMol *dp_mol;
  bool df_useNbrs;
  ChiralAtomCompareFunctor()
      : dp_atoms(NULL), dp_mol(NULL), df_useNbrs(false){};
  ChiralAtomCompareFunctor(Canon::canon_atom *atoms, const ROMol &m)
      : dp_atoms(atoms), dp_mol(&m), df_useNbrs(false){};
  int operator()(int i, int j) const {
    PRECONDITION(dp_atoms, "no atoms");
    PRECONDITION(dp_mol, "no molecule");
    PRECONDITION(i != j, "bad call");
    int v = basecomp(i, j);
    if (v) return v;

    if (df_useNbrs) {
      getAtomNeighborhood(dp_atoms[i].bonds);
      getAtomNeighborhood(dp_atoms[j].bonds);

      // we do two passes through the neighbor lists. The first just uses the
      // atomic numbers (by passing the optional 10000 to bondholder::compare),
      // the second takes the already-computed index into account
      for (unsigned int ii = 0;
           ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size();
           ++ii) {
        int cmp = bondholder::compare(
            dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii], ATNUM_CLASS_OFFSET);
        if (cmp) return cmp;
      }
      for (unsigned int ii = 0;
           ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size();
           ++ii) {
        int cmp =
            bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
        if (cmp) return cmp;
      }
      if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
        return -1;
      } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
        return 1;
      }
    }
    return 0;
  }
};

/*
 * Basic canonicalization function to organize the partitions which will be
 * sorted next.
 * */

template <typename CompareFunc>
void RefinePartitions(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
                      int mode, int *order, int *count, int &activeset,
                      int *next, int *changed, char *touchedPartitions) {
  unsigned int nAtoms = mol.getNumAtoms();
  int partition;
  int symclass = 0;
  int *start;
  int offset;
  int index;
  int len;
  int i;
  // std::vector<char> touchedPartitions(mol.getNumAtoms(),0);

  // std::cerr<<"&&&&&&&&&&&&&&&& RP"<<std::endl;
  while (activeset != -1) {
    // std::cerr<<"ITER: "<<activeset<<" next: "<<next[activeset]<<std::endl;
    // std::cerr<<" next: ";
    // for(unsigned int ii=0;ii<nAtoms;++ii){
    //   std::cerr<<ii<<":"<<next[ii]<<" ";
    // }
    // std::cerr<<std::endl;
    // for(unsigned int ii=0;ii<nAtoms;++ii){
    //   std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" index:
    //   "<<atoms[order[ii]].index<<std::endl;
    // }

    partition = activeset;
    activeset = next[partition];
    next[partition] = -2;

    len = count[partition];
    offset = atoms[partition].index;
    start = order + offset;
    // std::cerr<<"\n\n**************************************************************"<<std::endl;
    // std::cerr<<"  sort - class:"<<atoms[partition].index<<" len: "<<len<<":";
    // for(unsigned int ii=0;ii<len;++ii){
    //   std::cerr<<" "<<order[offset+ii]+1;
    // }
    // std::cerr<<std::endl;
    // for(unsigned int ii=0;ii<nAtoms;++ii){
    //   std::cerr<<order[ii]+1<<" count: "<<count[order[ii]]<<" index:
    //   "<<atoms[order[ii]].index<<std::endl;
    // }
    hanoisort(start, len, count, changed, compar);
    // std::cerr<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<std::endl;
    // std::cerr<<"  result:";
    // for(unsigned int ii=0;ii<nAtoms;++ii){
    //    std::cerr<<order[ii]+1<<" count: "<<count[order[ii]]<<" index:
    //    "<<atoms[order[ii]].index<<std::endl;
    //  }
    for (int k = 0; k < len; ++k) {
      changed[start[k]] = 0;
    }

    index = start[0];
    // std::cerr<<"  len:"<<len<<" index:"<<index<<"
    // count:"<<count[index]<<std::endl;
    for (i = count[index]; i < len; i++) {
      index = start[i];
      if (count[index]) symclass = offset + i;
      atoms[index].index = symclass;
      // std::cerr<<" "<<index+1<<"("<<symclass<<")";
      // if(mode && (activeset<0 || count[index]>count[activeset]) ){
      //  activeset=index;
      //}
      for (unsigned j = 0; j < atoms[index].degree; ++j) {
        changed[atoms[index].nbrIds[j]] = 1;
      }
    }
    // std::cerr<<std::endl;

    if (mode) {
      index = start[0];
      for (i = count[index]; i < len; i++) {
        index = start[i];
        for (unsigned j = 0; j < atoms[index].degree; ++j) {
          unsigned int nbor = atoms[index].nbrIds[j];
          touchedPartitions[atoms[nbor].index] = 1;
        }
      }
      for (unsigned int ii = 0; ii < nAtoms; ++ii) {
        if (touchedPartitions[ii]) {
          partition = order[ii];
          if ((count[partition] > 1) && (next[partition] == -2)) {
            next[partition] = activeset;
            activeset = partition;
          }
          touchedPartitions[ii] = 0;
        }
      }
    }
  }
}  // end of RefinePartitions()

template <typename CompareFunc>
void BreakTies(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
               int mode, int *order, int *count, int &activeset, int *next,
               int *changed, char *touchedPartitions) {
  unsigned int nAtoms = mol.getNumAtoms();
  int partition;
  int offset;
  int index;
  int len;
  int oldPart = 0;

  for (unsigned int i = 0; i < nAtoms; i++) {
    partition = order[i];
    oldPart = atoms[partition].index;
    while (count[partition] > 1) {
      len = count[partition];
      offset = atoms[partition].index + len - 1;
      index = order[offset];
      atoms[index].index = offset;
      count[partition] = len - 1;
      count[index] = 1;

      // test for ions, water molecules with no
      if (atoms[index].degree < 1) {
        continue;
      }
      for (unsigned j = 0; j < atoms[index].degree; ++j) {
        unsigned int nbor = atoms[index].nbrIds[j];
        touchedPartitions[atoms[nbor].index] = 1;
        changed[nbor] = 1;
      }

      for (unsigned int ii = 0; ii < nAtoms; ++ii) {
        if (touchedPartitions[ii]) {
          int npart = order[ii];
          if ((count[npart] > 1) && (next[npart] == -2)) {
            next[npart] = activeset;
            activeset = npart;
          }
          touchedPartitions[ii] = 0;
        }
      }
      RefinePartitions(mol, atoms, compar, mode, order, count, activeset, next,
                       changed, touchedPartitions);
    }
    // not sure if this works each time
    if (atoms[partition].index != oldPart) {
      i -= 1;
    }
  }
}  // end of BreakTies()

void CreateSinglePartition(unsigned int nAtoms, int *order, int *count,
                           canon_atom *atoms);

void ActivatePartitions(unsigned int nAtoms, int *order, int *count,
                        int &activeset, int *next, int *changed);

void rankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res,
                  bool breakTies = true, bool includeChirality = true,
                  bool includeIsotopes = true);

void rankFragmentAtoms(const ROMol &mol, std::vector<unsigned int> &res,
                       const boost::dynamic_bitset<> &atomsInPlay,
                       const boost::dynamic_bitset<> &bondsInPlay,
                       const std::vector<std::string> *atomSymbols = NULL,
                       bool breakTies = true, bool includeChirality = true,
                       bool includeIsotopes = true);

void chiralRankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res);

void initCanonAtoms(const ROMol &mol, std::vector<Canon::canon_atom> &atoms,
                    bool includeChirality = true);

}  // end of Canon namespace
}  // end of RDKit namespace
