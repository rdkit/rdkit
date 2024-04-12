//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SmilesWrite.h"
#include "SmilesParseOps.h"
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/Canon.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/utils.h>
#include <RDGeneral/BoostEndInclude.h>
#include <boost/format.hpp>

#include <sstream>
#include <map>
#include <list>

// #define VERBOSE_CANON 1

namespace RDKit {

namespace SmilesWrite {
const int atomicSmiles[] = {0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53, -1};
bool inOrganicSubset(int atomicNumber) {
  unsigned int idx = 0;
  while (atomicSmiles[idx] < atomicNumber && atomicSmiles[idx] != -1) {
    ++idx;
  }
  return atomicSmiles[idx] == atomicNumber;
}

namespace {
std::string getAtomChiralityInfo(const Atom *atom) {
  auto allowNontet = Chirality::getAllowNontetrahedralChirality();
  std::string atString;
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
  if (atString.empty() && allowNontet) {
    switch (atom->getChiralTag()) {
      case Atom::CHI_SQUAREPLANAR:
        atString = "@SP";
        break;
      case Atom::CHI_TRIGONALBIPYRAMIDAL:
        atString = "@TB";
        break;
      case Atom::CHI_OCTAHEDRAL:
        atString = "@OH";
        break;
      default:
        break;
    }
    if (!atString.empty()) {
      // we added info about non-tetrahedral stereo, so check whether or not
      // we need to also add permutation info
      int permutation = 0;
      if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER &&
          atom->getPropIfPresent(common_properties::_chiralPermutation,
                                 permutation) &&
          !SmilesParseOps::checkChiralPermutation(atom->getChiralTag(),
                                                  permutation)) {
        throw ValueErrorException("bad chirality spec");
      } else if (permutation) {
        atString += std::to_string(permutation);
      }
    }
  }
  return atString;
}
}  // namespace

std::string GetAtomSmiles(const Atom *atom, bool doKekule, const Bond *,
                          bool allHsExplicit, bool isomericSmiles) {
  PRECONDITION(atom, "bad atom");
  std::string res;
  int fc = atom->getFormalCharge();
  int num = atom->getAtomicNum();
  int isotope = atom->getIsotope();

  bool needsBracket = false;
  std::string symb;
  bool hasCustomSymbol =
      atom->getPropIfPresent(common_properties::smilesSymbol, symb);
  if (!hasCustomSymbol) {
    symb = PeriodicTable::getTable()->getElementSymbol(num);
  }

  // check for atomic stereochemistry
  std::string atString;
  if (isomericSmiles ||
      (atom->hasOwningMol() &&
       atom->getOwningMol().hasProp(common_properties::_doIsoSmiles))) {
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        !atom->hasProp(common_properties::_brokenChirality)) {
      atString = getAtomChiralityInfo(atom);
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

    if (hasCustomSymbol || atom->getNumRadicalElectrons()) {
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
    } else if ((isomericSmiles || (atom->hasOwningMol() &&
                                   atom->getOwningMol().hasProp(
                                       common_properties::_doIsoSmiles))) &&
               (isotope || atString != "")) {
      needsBracket = true;
    }
  } else {
    needsBracket = true;
  }
  if (needsBracket) {
    res += "[";
  }

  if (isotope && (isomericSmiles || (atom->hasOwningMol() &&
                                     atom->getOwningMol().hasProp(
                                         common_properties::_doIsoSmiles)))) {
    res += std::to_string(isotope);
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
      if (totNumHs > 1) {
        res += std::to_string(totNumHs);
      }
    }
    if (fc > 0) {
      res += "+";
      if (fc > 1) {
        res += std::to_string(fc);
      }
    } else if (fc < 0) {
      if (fc < -1) {
        res += std::to_string(fc);
      } else {
        res += "-";
      }
    }

    int mapNum;
    if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
      res += ":";
      res += std::to_string(mapNum);
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
  if (atomToLeftIdx < 0) {
    atomToLeftIdx = bond->getBeginAtomIdx();
  }

  std::string res = "";
  bool aromatic = false;
  if (!doKekule && (bond->getBondType() == Bond::SINGLE ||
                    bond->getBondType() == Bond::DOUBLE ||
                    bond->getBondType() == Bond::AROMATIC)) {
    if (bond->hasOwningMol()) {
      auto a1 = bond->getOwningMol().getAtomWithIdx(atomToLeftIdx);
      auto a2 = bond->getOwningMol().getAtomWithIdx(
          bond->getOtherAtomIdx(atomToLeftIdx));
      if ((a1->getIsAromatic() && a2->getIsAromatic()) &&
          (a1->getAtomicNum() || a2->getAtomicNum())) {
        aromatic = true;
      }
    } else {
      aromatic = false;
    }
  }

  Bond::BondDir dir = bond->getBondDir();

  bond->clearProp(common_properties::_TraversalRingClosureBond);

  switch (bond->getBondType()) {
    case Bond::SINGLE:
      if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        switch (dir) {
          case Bond::ENDDOWNRIGHT:
            if (allBondsExplicit || (bond->hasOwningMol() &&
                                     bond->getOwningMol().hasProp(
                                         common_properties::_doIsoSmiles))) {
              res = "\\";
            }
            break;
          case Bond::ENDUPRIGHT:
            if (allBondsExplicit || (bond->hasOwningMol() &&
                                     bond->getOwningMol().hasProp(
                                         common_properties::_doIsoSmiles))) {
              res = "/";
            }
            break;
          default:
            if (allBondsExplicit) {
              res = "-";
            }
            break;
        }
      } else {
        // if the bond is marked as aromatic and the two atoms
        //  are aromatic, we need no marker (this arises in kekulized
        //  molecules).
        // FIX: we should be able to dump kekulized smiles
        //   currently this is possible by removing all
        //   isAromatic flags, but there should maybe be another way
        if (allBondsExplicit) {
          res = "-";
        } else if (aromatic && !bond->getIsAromatic()) {
          res = "-";
        }
      }
      break;
    case Bond::DOUBLE:
      // see note above
      if (!aromatic || !bond->getIsAromatic() || allBondsExplicit) {
        res = "=";
      }
      break;
    case Bond::TRIPLE:
      res = "#";
      break;
    case Bond::QUADRUPLE:
      res = "$";
      break;
    case Bond::AROMATIC:
      if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        switch (dir) {
          case Bond::ENDDOWNRIGHT:
            if (allBondsExplicit || (bond->hasOwningMol() &&
                                     bond->getOwningMol().hasProp(
                                         common_properties::_doIsoSmiles))) {
              res = "\\";
            }
            break;
          case Bond::ENDUPRIGHT:
            if (allBondsExplicit || (bond->hasOwningMol() &&
                                     bond->getOwningMol().hasProp(
                                         common_properties::_doIsoSmiles))) {
              res = "/";
            }
            break;
          default:
            if (allBondsExplicit || !aromatic) {
              res = ":";
            }
            break;
        }
      } else if (allBondsExplicit || !aromatic) {
        res = ":";
      }
      break;
    case Bond::DATIVE:
      if (atomToLeftIdx >= 0 &&
          bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx)) {
        res = "->";
      } else {
        res = "<-";
      }
      break;
    default:
      res = "~";
  }
  return res;
}

std::string FragmentSmilesConstruct(
    ROMol &mol, int atomIdx, std::vector<Canon::AtomColors> &colors,
    const UINT_VECT &ranks, const SmilesWriteParams &params,
    std::vector<unsigned int> &atomOrdering,
    std::vector<unsigned int> &bondOrdering,
    const boost::dynamic_bitset<> *atomsInPlay = nullptr,
    const boost::dynamic_bitset<> *bondsInPlay = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr) {
  PRECONDITION(!atomsInPlay || atomsInPlay->size() >= mol.getNumAtoms(),
               "bad atomsInPlay");
  PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
               "bad bondsInPlay");
  PRECONDITION(!atomSymbols || atomSymbols->size() >= mol.getNumAtoms(),
               "bad atomSymbols");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bad bondSymbols");
  if (params.doKekule) {
    if (atomsInPlay && bondsInPlay) {
      MolOps::details::KekulizeFragment(static_cast<RWMol &>(mol), *atomsInPlay,
                                        *bondsInPlay);
    } else {
      MolOps::Kekulize(static_cast<RWMol &>(mol));
    }
  }

  Canon::MolStack molStack;
  // try to prevent excessive reallocation
  molStack.reserve(mol.getNumAtoms() + mol.getNumBonds());
  std::stringstream res;

  std::map<int, int> ringClosureMap;
  int ringIdx, closureVal;
  if (!params.canonical) {
    mol.setProp(common_properties::_StereochemDone, 1);
  }
  std::list<unsigned int> ringClosuresToErase;

  if (params.canonical && params.doIsomericSmiles) {
    Canon::canonicalizeEnhancedStereo(mol, &ranks);
  }
  Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack,
                              bondsInPlay, bondSymbols, params.doIsomericSmiles,
                              params.doRandom);
  Bond *bond = nullptr;
  for (auto &mSE : molStack) {
    switch (mSE.type) {
      case Canon::MOL_STACK_ATOM:
        for (auto rclosure : ringClosuresToErase) {
          ringClosureMap.erase(rclosure);
        }
        ringClosuresToErase.clear();
        // std::cout << "\t\tAtom: " << mSE.obj.atom->getIdx() << std::endl;
        if (!atomSymbols) {
          res << GetAtomSmiles(mSE.obj.atom, params.doKekule, bond,
                               params.allHsExplicit, params.doIsomericSmiles);
        } else {
          res << (*atomSymbols)[mSE.obj.atom->getIdx()];
        }
        atomOrdering.push_back(mSE.obj.atom->getIdx());
        break;
      case Canon::MOL_STACK_BOND:
        bond = mSE.obj.bond;
        // std::cout << "\t\tBond: " << bond->getIdx() << std::endl;
        if (!bondSymbols) {
          res << GetBondSmiles(bond, mSE.number, params.doKekule,
                               params.allBondsExplicit);
        } else {
          res << (*bondSymbols)[bond->getIdx()];
        }
        bondOrdering.push_back(bond->getIdx());
        break;
      case Canon::MOL_STACK_RING:
        ringIdx = mSE.number;
        // std::cout << "\t\tRing: " << ringIdx << std::endl;
        if (ringClosureMap.count(ringIdx)) {
          // the index is already in the map ->
          //   we're closing a ring, so grab
          //   the index and then delete the value:
          closureVal = ringClosureMap[ringIdx];
          ringClosuresToErase.push_back(ringIdx);
        } else {
          // we're opening a new ring, find the index for it:
          closureVal = 1;
          bool done = false;
          // EFF: there's got to be a more efficient way to do this
          while (!done) {
            std::map<int, int>::iterator mapIt;
            for (mapIt = ringClosureMap.begin(); mapIt != ringClosureMap.end();
                 ++mapIt) {
              if (mapIt->second == closureVal) {
                break;
              }
            }
            if (mapIt == ringClosureMap.end()) {
              done = true;
            } else {
              closureVal += 1;
            }
          }
          ringClosureMap[ringIdx] = closureVal;
        }
        if (closureVal < 10) {
          res << (char)(closureVal + '0');
        } else if (closureVal < 100) {
          res << '%' << closureVal;
        } else {  // use extension to OpenSMILES
          res << "%(" << closureVal << ')';
        }
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
    const std::pair<std::string, std::vector<unsigned int>> &a,
    const std::pair<std::string, std::vector<unsigned int>> &b) {
  return a.first < b.first;
}

namespace SmilesWrite {
namespace detail {
std::string MolToSmiles(const ROMol &mol, const SmilesWriteParams &params,
                        bool doingCXSmiles) {
  if (!mol.getNumAtoms()) {
    return "";
  }
  PRECONDITION(
      params.rootedAtAtom < 0 ||
          static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
      "rootedAtomAtom must be less than the number of atoms");
  int rootedAtAtom = params.rootedAtAtom;
  std::vector<std::vector<int>> fragsMolAtomMapping;
  auto mols =
      MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
  // we got the mapping between fragments and atoms; repeat that for bonds
  std::vector<std::vector<int>> fragsMolBondMapping;
  boost::dynamic_bitset<> atsPresent(mol.getNumAtoms());
  std::vector<int> bondsInFrag;
  bondsInFrag.reserve(mol.getNumBonds());
  for (const auto &atsInFrag : fragsMolAtomMapping) {
    atsPresent.reset();
    bondsInFrag.clear();
    for (auto aidx : atsInFrag) {
      atsPresent.set(aidx);
    }
    for (const auto bnd : mol.bonds()) {
      if (atsPresent[bnd->getBeginAtomIdx()] &&
          atsPresent[bnd->getEndAtomIdx()]) {
        bondsInFrag.push_back(bnd->getIdx());
      }
    }
    fragsMolBondMapping.push_back(bondsInFrag);
  }

  std::vector<std::string> vfragsmi(mols.size());

  //    for(unsigned i=0; i<fragsMolAtomMapping.size(); i++){
  //      std::cout << i << ": ";
  //      for(unsigned j=0; j<fragsMolAtomMapping[i].size(); j++){
  //        std::cout << j <<"("<<fragsMolAtomMapping[i][j]<<") ";
  //      }
  //      std::cout << std::endl;
  //    }

  std::vector<std::vector<RDKit::UINT>> allAtomOrdering;
  std::vector<std::vector<RDKit::UINT>> allBondOrdering;
  for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    ROMol *tmol = mols[fragIdx].get();

    // update property cache
    for (auto atom : tmol->atoms()) {
      atom->updatePropertyCache(false);
    }

    // clean up the chirality on any atom that is marked as chiral,
    // but that should not be:
    if (params.doIsomericSmiles) {
      tmol->setProp(common_properties::_doIsoSmiles, 1);
      if (!tmol->hasProp(common_properties::_StereochemDone)) {
        MolOps::assignStereochemistry(*tmol, true);
      }
    }
    if (!doingCXSmiles) {
      // remove any stereo groups that may be present. Otherwise they will be
      // used in the canonicalization
      std::vector<StereoGroup> noStereoGroups;
      tmol->setStereoGroups(noStereoGroups);
      // remove any wiggle bonds or unspecified double bond stereochemistry
      for (auto bond : tmol->bonds()) {
        if (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
            bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
          bond->setBondDir(Bond::BondDir::NONE);
        }
        if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
          bond->setStereo(Bond::BondStereo::STEREONONE);
        }
      }

      // if other CXSMILES features are added to the canonicalization code
      // in the future, they should be removed here.
    }
#if 0
      std::cout << "----------------------------" << std::endl;
      std::cout << "MolToSmiles:"<< std::endl;
      tmol->debugMol(std::cout);
      std::cout << "----------------------------" << std::endl;
#endif

    if (params.doRandom && rootedAtAtom == -1) {
      // need to find a random atom id between 0 and mol.getNumAtoms()
      // exclusively
      rootedAtAtom = getRandomGenerator()() % tmol->getNumAtoms();
    }

    std::string res;
    unsigned int nAtoms = tmol->getNumAtoms();
    std::vector<unsigned int> ranks(nAtoms);
    std::vector<unsigned int> atomOrdering;
    std::vector<unsigned int> bondOrdering;

    if (params.canonical) {
      if (tmol->hasProp("_canonicalRankingNumbers")) {
        for (const auto atom : tmol->atoms()) {
          unsigned int rankNum = 0;
          atom->getPropIfPresent("_canonicalRankingNumber", rankNum);
          ranks[atom->getIdx()] = rankNum;
        }
      } else {
        bool breakTies = true;
        Canon::rankMolAtoms(*tmol, ranks, breakTies, params.doIsomericSmiles,
                            params.doIsomericSmiles);
      }
    } else {
      std::iota(ranks.begin(), ranks.end(), 0);
    }
#ifdef VERBOSE_CANON
    for (unsigned int tmpI = 0; tmpI < ranks.size(); tmpI++) {
      std::cout << tmpI << " " << ranks[tmpI] << " "
                << *(tmol->getAtomWithIdx(tmpI)) << std::endl;
    }
#endif

    std::vector<Canon::AtomColors> colors(nAtoms, Canon::WHITE_NODE);
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
        *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);

    res += subSmi;
    vfragsmi[fragIdx] = res;

    for (unsigned int &vit : atomOrdering) {
      vit = fragsMolAtomMapping[fragIdx][vit];  // Lookup the Id in the original
                                                // molecule
    }
    allAtomOrdering.push_back(atomOrdering);
    for (unsigned int &vit : bondOrdering) {
      vit = fragsMolBondMapping[fragIdx][vit];  // Lookup the Id in the original
                                                // molecule
    }
    allBondOrdering.push_back(bondOrdering);
  }

  std::string result;
  std::vector<unsigned int> flattenedAtomOrdering;
  flattenedAtomOrdering.reserve(mol.getNumAtoms());
  std::vector<unsigned int> flattenedBondOrdering;
  flattenedBondOrdering.reserve(mol.getNumBonds());
  if (params.canonical) {
    // Sort the vfragsmi, but also sort the atom and bond order vectors into the
    // same order
    typedef std::tuple<std::string, std::vector<unsigned int>,
                       std::vector<unsigned int>>
        tplType;
    std::vector<tplType> tmp(vfragsmi.size());
    for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti) {
      tmp[ti] = std::make_tuple(vfragsmi[ti], allAtomOrdering[ti],
                                allBondOrdering[ti]);
    }
    std::sort(tmp.begin(), tmp.end());

    for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti) {
      result += std::get<0>(tmp[ti]);
      if (ti < vfragsmi.size() - 1) {
        result += ".";
      }
      flattenedAtomOrdering.insert(flattenedAtomOrdering.end(),
                                   std::get<1>(tmp[ti]).begin(),
                                   std::get<1>(tmp[ti]).end());
      flattenedBondOrdering.insert(flattenedBondOrdering.end(),
                                   std::get<2>(tmp[ti]).begin(),
                                   std::get<2>(tmp[ti]).end());
    }
  } else {  // Not canonical
    for (auto &i : allAtomOrdering) {
      flattenedAtomOrdering.insert(flattenedAtomOrdering.end(), i.begin(),
                                   i.end());
    }
    for (auto &i : allBondOrdering) {
      flattenedBondOrdering.insert(flattenedBondOrdering.end(), i.begin(),
                                   i.end());
    }
    for (unsigned i = 0; i < vfragsmi.size(); ++i) {
      result += vfragsmi[i];
      if (i < vfragsmi.size() - 1) {
        result += ".";
      }
    }
  }
  mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
              true);
  mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
              true);
  return result;
}
}  // namespace detail
}  // namespace SmilesWrite
std::string MolToSmiles(const ROMol &mol, const SmilesWriteParams &params) {
  bool doingCXSmiles = false;
  return SmilesWrite::detail::MolToSmiles(mol, params, doingCXSmiles);
}

std::string MolToCXSmiles(const ROMol &romol, const SmilesWriteParams &params,
                          std::uint32_t flags,
                          RestoreBondDirOption restoreBondDirs) {
  RWMol trwmol(romol);

  bool doingCXSmiles = true;

  auto res = SmilesWrite::detail::MolToSmiles(trwmol, params, doingCXSmiles);
  if (res.empty()) {
    return res;
  }
  if (restoreBondDirs == RestoreBondDirOptionTrue) {
    RDKit::Chirality::reapplyMolBlockWedging(trwmol);
  } else if (restoreBondDirs == RestoreBondDirOptionClear) {
    for (auto bond : trwmol.bonds()) {
      if (!canHaveDirection(*bond)) {
        continue;
      }
      if (bond->getBondDir() != Bond::BondDir::NONE) {
        bond->setBondDir(Bond::BondDir::NONE);
      }
      unsigned int cfg;
      if (bond->getPropIfPresent<unsigned int>(
              common_properties::_MolFileBondCfg, cfg)) {
        bond->clearProp(common_properties::_MolFileBondCfg);
      }
    }
  }

  if (!params.doIsomericSmiles) {
    flags &= ~(SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO |
               SmilesWrite::CXSmilesFields::CX_BOND_CFG);
  }

  auto cxext = SmilesWrite::getCXExtensions(trwmol, flags);
  if (!cxext.empty()) {
    res += " " + cxext;
  }
  return res;
}

std::vector<std::string> MolToRandomSmilesVect(
    const ROMol &mol, unsigned int numSmiles, unsigned int randomSeed,
    bool doIsomericSmiles, bool doKekule, bool allBondsExplicit,
    bool allHsExplicit) {
  if (randomSeed > 0) {
    getRandomGenerator(rdcast<int>(randomSeed));
  }
  std::vector<std::string> res;
  res.reserve(numSmiles);
  for (unsigned int i = 0; i < numSmiles; ++i) {
    bool canonical = false;
    int rootedAtAtom = -1;
    bool doRandom = true;
    res.push_back(MolToSmiles(mol, doIsomericSmiles, doKekule, rootedAtAtom,
                              canonical, allBondsExplicit, allHsExplicit,
                              doRandom));
  }
  return res;
};
std::string MolFragmentToSmiles(const ROMol &mol,
                                const SmilesWriteParams &params,
                                const std::vector<int> &atomsToUse,
                                const std::vector<int> *bondsToUse,
                                const std::vector<std::string> *atomSymbols,
                                const std::vector<std::string> *bondSymbols) {
  PRECONDITION(atomsToUse.size(), "no atoms provided");
  PRECONDITION(
      params.rootedAtAtom < 0 ||
          static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
      "rootedAtomAtom must be less than the number of atoms");
  PRECONDITION(params.rootedAtAtom < 0 ||
                   std::find(atomsToUse.begin(), atomsToUse.end(),
                             params.rootedAtAtom) != atomsToUse.end(),
               "rootedAtAtom not found in atomsToUse");
  PRECONDITION(!atomSymbols || atomSymbols->size() >= mol.getNumAtoms(),
               "bad atomSymbols vector");
  PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
               "bad bondSymbols vector");
  if (!mol.getNumAtoms()) {
    return "";
  }
  int rootedAtAtom = params.rootedAtAtom;

  ROMol tmol(mol, true);
  if (params.doIsomericSmiles) {
    tmol.setProp(common_properties::_doIsoSmiles, 1);
  }
  std::string res;

  boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms(), 0);
  for (auto aidx : atomsToUse) {
    atomsInPlay.set(aidx);
  }
  // figure out which bonds are actually in play:
  boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds(), 0);
  if (bondsToUse) {
    for (auto bidx : *bondsToUse) {
      bondsInPlay.set(bidx);
    }
  } else {
    for (auto aidx : atomsToUse) {
      for (const auto &bndi : boost::make_iterator_range(
               mol.getAtomBonds(mol.getAtomWithIdx(aidx)))) {
        const Bond *bond = mol[bndi];
        if (atomsInPlay[bond->getOtherAtomIdx(aidx)]) {
          bondsInPlay.set(bond->getIdx());
        }
      }
    }
  }

  // copy over the rings that only involve atoms/bonds in this fragment:
  if (mol.getRingInfo()->isInitialized()) {
    tmol.getRingInfo()->reset();
    tmol.getRingInfo()->initialize();
    for (unsigned int ridx = 0; ridx < mol.getRingInfo()->numRings(); ++ridx) {
      const INT_VECT &aring = mol.getRingInfo()->atomRings()[ridx];
      bool keepIt = true;
      for (auto aidx : aring) {
        if (!atomsInPlay[aidx]) {
          keepIt = false;
          break;
        }
      }
      if (keepIt) {
        const INT_VECT &bring = mol.getRingInfo()->bondRings()[ridx];
        for (auto bidx : bring) {
          if (!bondsInPlay[bidx]) {
            keepIt = false;
            break;
          }
        }
        if (keepIt) {
          tmol.getRingInfo()->addRing(aring, bring);
        }
      }
    }
  }
  if (tmol.needsUpdatePropertyCache()) {
    for (auto atom : tmol.atoms()) {
      atom->updatePropertyCache(false);
    }
  }

  UINT_VECT ranks(tmol.getNumAtoms());

  std::vector<unsigned int> atomOrdering;
  std::vector<unsigned int> bondOrdering;

  // clean up the chirality on any atom that is marked as chiral,
  // but that should not be:
  if (params.doIsomericSmiles) {
    if (!mol.hasProp(common_properties::_StereochemDone)) {
      MolOps::assignStereochemistry(tmol, true);
    } else {
      tmol.setProp(common_properties::_StereochemDone, 1);
      // we need the CIP codes:
      for (auto aidx : atomsToUse) {
        const Atom *oAt = mol.getAtomWithIdx(aidx);
        std::string cipCode;
        if (oAt->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
          tmol.getAtomWithIdx(aidx)->setProp(common_properties::_CIPCode,
                                             cipCode);
        }
      }
    }
  }
  if (params.canonical) {
    bool breakTies = true;
    Canon::rankFragmentAtoms(tmol, ranks, atomsInPlay, bondsInPlay, atomSymbols,
                             breakTies, params.doIsomericSmiles,
                             params.doIsomericSmiles);
    // std::cerr << "RANKS: ";
    // std::copy(ranks.begin(), ranks.end(),
    //           std::ostream_iterator<int>(std::cerr, " "));
    // std::cerr << std::endl;
    // MolOps::rankAtomsInFragment(tmol,ranks,atomsInPlay,bondsInPlay,atomSymbols,bondSymbols);
  } else {
    for (unsigned int i = 0; i < tmol.getNumAtoms(); ++i) {
      ranks[i] = i;
    }
  }
#ifdef VERBOSE_CANON
  for (unsigned int tmpI = 0; tmpI < ranks.size(); tmpI++) {
    std::cout << tmpI << " " << ranks[tmpI] << " "
              << *(tmol.getAtomWithIdx(tmpI)) << std::endl;
  }
#endif

  std::vector<Canon::AtomColors> colors(tmol.getNumAtoms(), Canon::BLACK_NODE);
  for (auto aidx : atomsToUse) {
    colors[aidx] = Canon::WHITE_NODE;
  }
  std::vector<Canon::AtomColors>::iterator colorIt;
  colorIt = colors.begin();
  // loop to deal with the possibility that there might be disconnected
  // fragments
  while (colorIt != colors.end()) {
    int nextAtomIdx = -1;

    // find the next atom for a traverse
    if (rootedAtAtom >= 0) {
      nextAtomIdx = rootedAtAtom;
      rootedAtAtom = -1;
    } else {
      unsigned int nextRank = rdcast<unsigned int>(tmol.getNumAtoms()) + 1;
      for (auto i : atomsToUse) {
        if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
          nextRank = ranks[i];
          nextAtomIdx = i;
        }
      }
    }
    CHECK_INVARIANT(nextAtomIdx >= 0, "no start atom found");
    auto subSmi = SmilesWrite::FragmentSmilesConstruct(
        tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering,
        &atomsInPlay, &bondsInPlay, atomSymbols, bondSymbols);

    res += subSmi;
    colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    if (colorIt != colors.end()) {
      res += ".";
    }
  }

  mol.setProp(common_properties::_smilesAtomOutputOrder, atomOrdering, true);
  mol.setProp(common_properties::_smilesBondOutputOrder, bondOrdering, true);

  return res;
}  // end of MolFragmentToSmiles()

std::string MolFragmentToCXSmiles(const ROMol &mol,
                                  const SmilesWriteParams &params,
                                  const std::vector<int> &atomsToUse,
                                  const std::vector<int> *bondsToUse,
                                  const std::vector<std::string> *atomSymbols,
                                  const std::vector<std::string> *bondSymbols) {
  auto res = MolFragmentToSmiles(mol, params, atomsToUse, bondsToUse,
                                 atomSymbols, bondSymbols);
  auto cxext = SmilesWrite::getCXExtensions(mol);
  if (!cxext.empty()) {
    res += " " + cxext;
  }
  return res;
}
}  // namespace RDKit
