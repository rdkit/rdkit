//
//  Copyright (C) 2002-2016 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SmilesWrite.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/Canon.h>
#include <GraphMol/new_canon.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <sstream>
#include <map>
#include <list>

//#define VERBOSE_CANON 1

namespace RDKit {

namespace SmilesWrite {
const int atomicSmiles[] = {5, 6, 7, 8, 9, 15, 16, 17, 35, 53, -1};
bool inOrganicSubset(int atomicNumber) {
  unsigned int idx = 0;
  while (atomicSmiles[idx] < atomicNumber && atomicSmiles[idx] != -1) {
    ++idx;
  }
  if (atomicSmiles[idx] == atomicNumber) {
    return true;
  }
  return false;
}

std::string GetAtomSmiles(const Atom *atom, bool doKekule, const Bond *bondIn,
                          bool allHsExplicit) {
  RDUNUSED_PARAM(bondIn);
  PRECONDITION(atom, "bad atom");
  INT_VECT atomicSmilesVect(
      atomicSmiles,
      atomicSmiles + (sizeof(atomicSmiles) - 1) / sizeof(atomicSmiles[0]));
  std::string res;
  int fc = atom->getFormalCharge();
  int num = atom->getAtomicNum();
  int isotope = atom->getIsotope();

  bool needsBracket = false;
  std::string symb;
  if (!atom->getPropIfPresent(common_properties::smilesSymbol, symb)) {
    symb = PeriodicTable::getTable()->getElementSymbol(num);
  }

  // check for atomic stereochemistry
  std::string atString = "";
  if (atom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        !atom->hasProp(common_properties::_brokenChirality)) {
      switch (atom->getChiralTag()) {
        case Atom::CHI_TETRAHEDRAL_CW:
          atString = "@@";
          break;
        case Atom::CHI_TETRAHEDRAL_CCW:
          atString = "@";
          break;
        default:
          break;
      }
    }
  }

  if (!allHsExplicit && inOrganicSubset(num)) {
    // it's a member of the organic subset

    // -----
    // figure out if we need to put a bracket around the atom,
    // the conditions for this are:
    //   - formal charge specified
    //   - the atom has a nonstandard valence
    //   - chirality present and writing isomeric smiles
    //   - non-default isotope and writing isomeric smiles
    //   - atom-map information present
    const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    int totalValence = atom->getTotalValence();
    bool nonStandard = false;

    if (atom->getNumRadicalElectrons()) {
      nonStandard = true;
    } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
               atom->getNumExplicitHs()) {
      // another type of "nonstandard" valence is an aromatic N or P with
      // explicit Hs indicated:
      nonStandard = true;
    } else {
      nonStandard =
          (totalValence != defaultVs.front() && atom->getTotalNumHs());
    }

    if (fc || nonStandard ||
        atom->hasProp(common_properties::molAtomMapNumber)) {
      needsBracket = true;
    } else if (atom->getOwningMol().hasProp(common_properties::_doIsoSmiles) &&
               (isotope || atString != "")) {
      needsBracket = true;
    }
  } else {
    needsBracket = true;
  }
  if (needsBracket) res += "[";

  if (isotope &&
      atom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    res += boost::lexical_cast<std::string>(isotope);
  }
  // this was originally only done for the organic subset,
  // applying it to other atom-types is a fix for Issue 3152751:
  if (!doKekule && atom->getIsAromatic() && symb[0] >= 'A' && symb[0] <= 'Z') {
    symb[0] -= ('A' - 'a');
  }
  res += symb;

  res += atString;

  if (needsBracket) {
    unsigned int totNumHs = atom->getTotalNumHs();
    if (totNumHs > 0) {
      res += "H";
      if (totNumHs > 1) res += boost::lexical_cast<std::string>(totNumHs);
    }
    if (fc > 0) {
      res += "+";
      if (fc > 1) res += boost::lexical_cast<std::string>(fc);
    } else if (fc < 0) {
      if (fc < -1)
        res += boost::lexical_cast<std::string>(fc);
      else
        res += "-";
    }

    int mapNum;
    if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
      res += ":";
      res += boost::lexical_cast<std::string>(mapNum);
    }
    res += "]";
  }

  // If the atom has this property, the contained string will
  // be inserted directly in the SMILES:
  std::string label;
  if (atom->getPropIfPresent(common_properties::_supplementalSmilesLabel,
                             label)) {
    res += label;
  }

  return res;
}

std::string GetBondSmiles(const Bond *bond, int atomToLeftIdx, bool doKekule,
                          bool allBondsExplicit) {
  PRECONDITION(bond, "bad bond");
  if (atomToLeftIdx < 0) atomToLeftIdx = bond->getBeginAtomIdx();

  std::string res = "";
  bool aromatic = false;
  if (!doKekule && (bond->getBondType() == Bond::SINGLE ||
                    bond->getBondType() == Bond::DOUBLE ||
                    bond->getBondType() == Bond::AROMATIC)) {
    Atom *a1, *a2;
    a1 = bond->getOwningMol().getAtomWithIdx(atomToLeftIdx);
    a2 = bond->getOwningMol().getAtomWithIdx(
        bond->getOtherAtomIdx(atomToLeftIdx));
    if ((a1->getIsAromatic() && a2->getIsAromatic()) &&
        (a1->getAtomicNum() || a2->getAtomicNum()))
      aromatic = true;
  }

  Bond::BondDir dir = bond->getBondDir();

  if (bond->hasProp(common_properties::_TraversalRingClosureBond)) {
    // std::cerr<<"FLIP: "<<bond->getIdx()<<"
    // "<<bond->getBeginAtomIdx()<<"-"<<bond->getEndAtomIdx()<<std::endl;
    // if(dir==Bond::ENDDOWNRIGHT) dir=Bond::ENDUPRIGHT;
    // else if(dir==Bond::ENDUPRIGHT) dir=Bond::ENDDOWNRIGHT;
    bond->clearProp(common_properties::_TraversalRingClosureBond);
  }

  switch (bond->getBondType()) {
    case Bond::SINGLE:
      if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        switch (dir) {
          case Bond::ENDDOWNRIGHT:
            if (bond->getOwningMol().hasProp(common_properties::_doIsoSmiles))
              res = "\\";
            break;
          case Bond::ENDUPRIGHT:
            if (bond->getOwningMol().hasProp(common_properties::_doIsoSmiles))
              res = "/";
            break;
          default:
            break;
        }
      } else {
        // if the bond is marked as aromatic and the two atoms
        //  are aromatic, we need no marker (this arises in kekulized
        //  molecules).
        // FIX: we should be able to dump kekulized smiles
        //   currently this is possible by removing all
        //   isAromatic flags, but there should maybe be another way
        if (allBondsExplicit)
          res = "-";
        else if (aromatic && !bond->getIsAromatic())
          res = "-";
      }
      break;
    case Bond::DOUBLE:
      // see note above
      if (!aromatic || !bond->getIsAromatic()) res = "=";
      break;
    case Bond::TRIPLE:
      res = "#";
      break;
    case Bond::AROMATIC:
      if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        switch (dir) {
          case Bond::ENDDOWNRIGHT:
            if (bond->getOwningMol().hasProp(common_properties::_doIsoSmiles))
              res = "\\";
            break;
          case Bond::ENDUPRIGHT:
            if (bond->getOwningMol().hasProp(common_properties::_doIsoSmiles))
              res = "/";
            break;
          default:
            break;
        }
      } else if (allBondsExplicit || !aromatic) {
        res = ":";
      }
      break;
    case Bond::DATIVE:
      if (atomToLeftIdx >= 0 &&
          bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx))
        res = "->";
      else
        res = "<-";
      break;
    default:
      res = "~";
  }
  return res;
}

std::string FragmentSmilesConstruct(
    ROMol &mol, int atomIdx, std::vector<Canon::AtomColors> &colors,
    const UINT_VECT &ranks, bool doKekule, bool canonical,
    bool doIsomericSmiles, bool allBondsExplicit, bool allHsExplicit,
    std::vector<unsigned int> &atomOrdering,
    const boost::dynamic_bitset<> *bondsInPlay = 0,
    const std::vector<std::string> *atomSymbols = 0,
    const std::vector<std::string> *bondSymbols = 0) {
  PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
               "bad bondsInPlay");
  PRECONDITION(!atomSymbols || atomSymbols->size() >= mol.getNumAtoms(),
               "bad atomSymbols");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bad bondSymbols");

  Canon::MolStack molStack;
  // try to prevent excessive reallocation
  molStack.reserve(mol.getNumAtoms() + mol.getNumBonds());
  std::stringstream res;

  std::map<int, int> ringClosureMap;
  int ringIdx, closureVal;
  if (!canonical) mol.setProp(common_properties::_StereochemDone, 1);
  std::list<unsigned int> ringClosuresToErase;

  Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack,
                              bondsInPlay, bondSymbols, doIsomericSmiles);
  Bond *bond = 0;
  BOOST_FOREACH (Canon::MolStackElem mSE, molStack) {
    switch (mSE.type) {
      case Canon::MOL_STACK_ATOM:
        if (!ringClosuresToErase.empty()) {
          BOOST_FOREACH (unsigned int rclosure, ringClosuresToErase) {
            ringClosureMap.erase(rclosure);
          }
          ringClosuresToErase.clear();
        }
        // std::cout<<"\t\tAtom: "<<mSE.obj.atom->getIdx()<<std::endl;
        if (!atomSymbols) {
          res << GetAtomSmiles(mSE.obj.atom, doKekule, bond, allHsExplicit);
        } else {
          res << (*atomSymbols)[mSE.obj.atom->getIdx()];
        }
        atomOrdering.push_back(mSE.obj.atom->getIdx());

        break;
      case Canon::MOL_STACK_BOND:
        bond = mSE.obj.bond;
        // std::cout<<"\t\tBond: "<<bond->getIdx()<<std::endl;
        if (!bondSymbols) {
          res << GetBondSmiles(bond, mSE.number, doKekule, allBondsExplicit);
        } else {
          res << (*bondSymbols)[bond->getIdx()];
        }
        break;
      case Canon::MOL_STACK_RING:
        ringIdx = mSE.number;
        // std::cout<<"\t\tRing: "<<ringIdx;
        if (ringClosureMap.count(ringIdx)) {
          // the index is already in the map ->
          //   we're closing a ring, so grab
          //   the index and then delete the value:
          closureVal = ringClosureMap[ringIdx];
          // ringClosureMap.erase(ringIdx);
          ringClosuresToErase.push_back(ringIdx);
        } else {
          // we're opening a new ring, find the index for it:
          closureVal = 1;
          bool done = false;
          // EFF: there's got to be a more efficient way to do this
          while (!done) {
            std::map<int, int>::iterator mapIt;
            for (mapIt = ringClosureMap.begin(); mapIt != ringClosureMap.end();
                 mapIt++) {
              if (mapIt->second == closureVal) break;
            }
            if (mapIt == ringClosureMap.end()) {
              done = true;
            } else {
              closureVal += 1;
            }
          }
          ringClosureMap[ringIdx] = closureVal;
        }
        if (closureVal >= 10) {
          res << "%";
        }
        // std::cerr << " > " << closureVal <<std::endl;
        res << closureVal;
        break;
      case Canon::MOL_STACK_BRANCH_OPEN:
        res << "(";
        break;
      case Canon::MOL_STACK_BRANCH_CLOSE:
        res << ")";
        break;
      default:
        break;
    }
  }
  return res.str();
}

}  // end of namespace SmilesWrite

static bool SortBasedOnFirstElement(
    const std::pair<std::string, std::vector<unsigned int> > &a,
    const std::pair<std::string, std::vector<unsigned int> > &b) {
  return a.first < b.first;
}

std::string MolToSmiles(const ROMol &mol, bool doIsomericSmiles, bool doKekule,
                        int rootedAtAtom, bool canonical, bool allBondsExplicit,
                        bool allHsExplicit) {
  if (!mol.getNumAtoms()) return "";
  PRECONDITION(rootedAtAtom < 0 ||
                   static_cast<unsigned int>(rootedAtAtom) < mol.getNumAtoms(),
               "rootedAtomAtom must be less than the number of atoms");

  std::vector<std::vector<int> > fragsMolAtomMapping;
  std::vector<ROMOL_SPTR> mols =
      MolOps::getMolFrags(mol, false, NULL, &fragsMolAtomMapping, false);
  std::vector<std::string> vfragsmi;

  //    for(unsigned i=0; i<fragsMolAtomMapping.size(); i++){
  //      std::cout << i << ": ";
  //      for(unsigned j=0; j<fragsMolAtomMapping[i].size(); j++){
  //        std::cout << j <<"("<<fragsMolAtomMapping[i][j]<<") ";
  //      }
  //      std::cout << std::endl;
  //    }

  std::vector<std::vector<RDKit::UINT> > allAtomOrdering;
  for (unsigned i = 0; i < mols.size(); i++) {
    ROMol *tmol = mols[i].get();

    // update property cache
    for (ROMol::AtomIterator atomIt = tmol->beginAtoms();
         atomIt != tmol->endAtoms(); ++atomIt) {
      (*atomIt)->updatePropertyCache(false);
    }

    // clean up the chirality on any atom that is marked as chiral,
    // but that should not be:
    if (doIsomericSmiles) {
      tmol->setProp(common_properties::_doIsoSmiles, 1);
      if (!mol.hasProp(common_properties::_StereochemDone)) {
        MolOps::assignStereochemistry(*tmol, true);
      }
    }
#if 0
      std::cout << "----------------------------" << std::endl;
      std::cout << "MolToSmiles:"<< std::endl;
      tmol->debugMol(std::cout);
      std::cout << "----------------------------" << std::endl;
#endif

    std::string res;
    unsigned int nAtoms = tmol->getNumAtoms();
    UINT_VECT ranks(nAtoms);
    std::vector<unsigned int> atomOrdering;

    if (canonical) {
      if (tmol->hasProp("_canonicalRankingNumbers")) {
        for (unsigned int i = 0; i < tmol->getNumAtoms(); ++i) {
          unsigned int rankNum = 0;
          tmol->getAtomWithIdx(i)->getPropIfPresent("_canonicalRankingNumber",
                                                    rankNum);
          ranks[i] = rankNum;
        }
      } else {
        Canon::rankMolAtoms(*tmol, ranks, true, doIsomericSmiles,
                            doIsomericSmiles);
      }
    } else {
      for (unsigned int i = 0; i < tmol->getNumAtoms(); ++i) ranks[i] = i;
    }
#ifdef VERBOSE_CANON
    for (unsigned int tmpI = 0; tmpI < ranks.size(); tmpI++) {
      std::cout << tmpI << " " << ranks[tmpI] << " "
                << *(tmol.getAtomWithIdx(tmpI)) << std::endl;
    }
#endif

    std::vector<Canon::AtomColors> colors(nAtoms, Canon::WHITE_NODE);
    std::vector<Canon::AtomColors>::iterator colorIt;
    colorIt = colors.begin();
    // loop to deal with the possibility that there might be disconnected
    // fragments
    while (colorIt != colors.end()) {
      int nextAtomIdx = -1;
      std::string subSmi;

      // find the next atom for a traverse
      if (rootedAtAtom >= 0) {
        nextAtomIdx = rootedAtAtom;
        rootedAtAtom = -1;
      } else {
        unsigned int nextRank = nAtoms + 1;
        for (unsigned int i = 0; i < nAtoms; i++) {
          if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
            nextRank = ranks[i];
            nextAtomIdx = i;
          }
        }
      }
      CHECK_INVARIANT(nextAtomIdx >= 0, "no start atom found");

      subSmi = SmilesWrite::FragmentSmilesConstruct(
          *tmol, nextAtomIdx, colors, ranks, doKekule, canonical,
          doIsomericSmiles, allBondsExplicit, allHsExplicit, atomOrdering);

      res += subSmi;
      colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
      if (colorIt != colors.end()) {
        res += ".";
      }
    }
    vfragsmi.push_back(res);

    for (std::vector<RDKit::UINT>::iterator vit = atomOrdering.begin();
         vit != atomOrdering.end(); ++vit) {
      *vit = fragsMolAtomMapping[i][*vit];  // Lookup the Id in the original
                                            // molecule
    }
    allAtomOrdering.push_back(atomOrdering);
  }

  std::string result;
  std::vector<unsigned int> flattenedAtomOrdering;
  if (canonical) {
    // Sort the vfragsmi, but also sort the atom order vectors into the same
    // order
    typedef std::pair<std::string, std::vector<unsigned int> > PairStrAndVec;
    std::vector<PairStrAndVec> tmp(vfragsmi.size());
    for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti)
      tmp[ti] = PairStrAndVec(vfragsmi[ti], allAtomOrdering[ti]);

    std::sort(tmp.begin(), tmp.end(), SortBasedOnFirstElement);

    for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti) {
      result += tmp[ti].first;
      if (ti < vfragsmi.size() - 1) result += ".";
      flattenedAtomOrdering.insert(flattenedAtomOrdering.end(),
                                   tmp[ti].second.begin(),
                                   tmp[ti].second.end());
    }
  } else {  // Not canonical
    for (unsigned int i = 0; i < allAtomOrdering.size(); ++i)
      flattenedAtomOrdering.insert(flattenedAtomOrdering.end(),
                                   allAtomOrdering[i].begin(),
                                   allAtomOrdering[i].end());
    for (unsigned i = 0; i < vfragsmi.size(); ++i) {
      result += vfragsmi[i];
      if (i < vfragsmi.size() - 1) {
        result += ".";
      }
    }
  }
  mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
              true);
  return result;
}  // end of MolToSmiles()

std::string MolFragmentToSmiles(const ROMol &mol,
                                const std::vector<int> &atomsToUse,
                                const std::vector<int> *bondsToUse,
                                const std::vector<std::string> *atomSymbols,
                                const std::vector<std::string> *bondSymbols,
                                bool doIsomericSmiles, bool doKekule,
                                int rootedAtAtom, bool canonical,
                                bool allBondsExplicit, bool allHsExplicit) {
  PRECONDITION(atomsToUse.size(), "no atoms provided");
  PRECONDITION(rootedAtAtom < 0 ||
                   static_cast<unsigned int>(rootedAtAtom) < mol.getNumAtoms(),
               "rootedAtomAtom must be less than the number of atoms");
  PRECONDITION(rootedAtAtom < 0 ||
                   std::find(atomsToUse.begin(), atomsToUse.end(),
                             rootedAtAtom) != atomsToUse.end(),
               "rootedAtomAtom not found in atomsToUse");
  PRECONDITION(!atomSymbols || atomSymbols->size() >= mol.getNumAtoms(),
               "bad atomSymbols vector");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bad bondSymbols vector");
  if (!mol.getNumAtoms()) return "";

  ROMol tmol(mol, true);
  if (doIsomericSmiles) {
    tmol.setProp(common_properties::_doIsoSmiles, 1);
  }
  std::string res;

  boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms(), 0);
  BOOST_FOREACH (int aidx, atomsToUse) { atomsInPlay.set(aidx); }
  // figure out which bonds are actually in play:
  boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds(), 0);
  if (bondsToUse) {
    BOOST_FOREACH (int bidx, *bondsToUse) { bondsInPlay.set(bidx); }
  } else {
    BOOST_FOREACH (int aidx, atomsToUse) {
      ROMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = mol.getAtomBonds(mol.getAtomWithIdx(aidx));
      while (beg != end) {
        const BOND_SPTR bond = mol[*beg];
        if (atomsInPlay[bond->getOtherAtomIdx(aidx)])
          bondsInPlay.set(bond->getIdx());
        ++beg;
      }
    }
  }

  // copy over the rings that only involve atoms/bonds in this fragment:
  if (mol.getRingInfo()->isInitialized()) {
    tmol.getRingInfo()->reset();
    tmol.getRingInfo()->initialize();
    for (unsigned int ridx = 0; ridx < mol.getRingInfo()->numRings(); ++ridx) {
      const INT_VECT &aring = mol.getRingInfo()->atomRings()[ridx];
      const INT_VECT &bring = mol.getRingInfo()->bondRings()[ridx];
      bool keepIt = true;
      BOOST_FOREACH (int aidx, aring) {
        if (!atomsInPlay[aidx]) {
          keepIt = false;
          break;
        }
      }
      if (keepIt) {
        BOOST_FOREACH (int bidx, bring) {
          if (!bondsInPlay[bidx]) {
            keepIt = false;
            break;
          }
        }
      }
      if (keepIt) {
        tmol.getRingInfo()->addRing(aring, bring);
      }
    }
  }
  if (tmol.needsUpdatePropertyCache()) {
    for (ROMol::AtomIterator atIt = tmol.beginAtoms(); atIt != tmol.endAtoms();
         atIt++) {
      (*atIt)->updatePropertyCache(false);
    }
  }

  UINT_VECT ranks(tmol.getNumAtoms());

  std::vector<unsigned int> atomOrdering;

  // clean up the chirality on any atom that is marked as chiral,
  // but that should not be:
  if (doIsomericSmiles) {
    if (!mol.hasProp(common_properties::_StereochemDone)) {
      MolOps::assignStereochemistry(tmol, true);
    } else {
      tmol.setProp(common_properties::_StereochemDone, 1);
      // we need the CIP codes:
      BOOST_FOREACH (int aidx, atomsToUse) {
        const Atom *oAt = mol.getAtomWithIdx(aidx);
        std::string cipCode;
        if (oAt->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
          tmol.getAtomWithIdx(aidx)->setProp(common_properties::_CIPCode,
                                             cipCode);
        }
      }
    }
  }
  if (canonical) {
    Canon::rankFragmentAtoms(tmol, ranks, atomsInPlay, bondsInPlay, atomSymbols,
                             true, doIsomericSmiles, doIsomericSmiles);
    // MolOps::rankAtomsInFragment(tmol,ranks,atomsInPlay,bondsInPlay,atomSymbols,bondSymbols);
  } else {
    for (unsigned int i = 0; i < tmol.getNumAtoms(); ++i) ranks[i] = i;
  }
#ifdef VERBOSE_CANON
  for (unsigned int tmpI = 0; tmpI < ranks.size(); tmpI++) {
    std::cout << tmpI << " " << ranks[tmpI] << " "
              << *(tmol.getAtomWithIdx(tmpI)) << std::endl;
  }
#endif

  std::vector<Canon::AtomColors> colors(tmol.getNumAtoms(), Canon::BLACK_NODE);
  BOOST_FOREACH (int aidx, atomsToUse) { colors[aidx] = Canon::WHITE_NODE; }
  std::vector<Canon::AtomColors>::iterator colorIt;
  colorIt = colors.begin();
  // loop to deal with the possibility that there might be disconnected
  // fragments
  while (colorIt != colors.end()) {
    int nextAtomIdx = -1;
    std::string subSmi;

    // find the next atom for a traverse
    if (rootedAtAtom >= 0) {
      nextAtomIdx = rootedAtAtom;
      rootedAtAtom = -1;
    } else {
      unsigned int nextRank = rdcast<unsigned int>(tmol.getNumAtoms()) + 1;
      BOOST_FOREACH (int i, atomsToUse) {
        if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
          nextRank = ranks[i];
          nextAtomIdx = i;
        }
      }
    }
    CHECK_INVARIANT(nextAtomIdx >= 0, "no start atom found");

    subSmi = SmilesWrite::FragmentSmilesConstruct(
        tmol, nextAtomIdx, colors, ranks, doKekule, canonical, doIsomericSmiles,
        allBondsExplicit, allHsExplicit, atomOrdering, &bondsInPlay,
        atomSymbols, bondSymbols);

    res += subSmi;
    colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    if (colorIt != colors.end()) {
      res += ".";
    }
  }
  mol.setProp(common_properties::_smilesAtomOutputOrder, atomOrdering, true);
  return res;
}  // end of MolFragmentToSmiles()
}
