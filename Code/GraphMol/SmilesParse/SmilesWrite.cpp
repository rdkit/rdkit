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
#include "SmilesParse.h"
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
#include <boost/algorithm/string.hpp>

#include <RDGeneral/utils.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <boost/format.hpp>

#include <sstream>
#include <map>
#include <list>

// #define VERBOSE_CANON 1

namespace RDKit {

void updateSmilesWriteParamsFromJSON(SmilesWriteParams &params,
                                     const char *details_json) {
  if (details_json && strlen(details_json)) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    params.doIsomericSmiles =
        pt.get("doIsomericSmiles", params.doIsomericSmiles);
    params.doKekule = pt.get("doKekule", params.doKekule);
    params.rootedAtAtom = pt.get("rootedAtAtom", params.rootedAtAtom);
    params.canonical = pt.get("canonical", params.canonical);
    params.allBondsExplicit =
        pt.get("allBondsExplicit", params.allBondsExplicit);
    params.allHsExplicit = pt.get("allHsExplicit", params.allHsExplicit);
    params.doRandom = pt.get("doRandom", params.doRandom);
  }
}

void updateSmilesWriteParamsFromJSON(SmilesWriteParams &params,
                                     const std::string &details_json) {
  updateSmilesWriteParamsFromJSON(params, details_json.c_str());
}

void updateCXSmilesFieldsFromJSON(SmilesWrite::CXSmilesFields &cxSmilesFields,
                                  RestoreBondDirOption &restoreBondDirs,
                                  const char *details_json) {
  static const auto cxSmilesFieldsKeyValuePairs = CXSMILESFIELDS_ITEMS_MAP;
  static const auto restoreBondDirOptionKeyValuePairs =
      RESTOREBONDDIROPTION_ITEMS_MAP;
  if (details_json && strlen(details_json)) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    auto cxSmilesFieldsFromJson =
        static_cast<std::underlying_type<SmilesWrite::CXSmilesFields>::type>(
            SmilesWrite::CXSmilesFields::CX_NONE);
    for (const auto &keyValuePair : cxSmilesFieldsKeyValuePairs) {
      cxSmilesFieldsFromJson |= (pt.get(keyValuePair.first, false)
                                     ? keyValuePair.second
                                     : SmilesWrite::CXSmilesFields::CX_NONE);
    }
    if (cxSmilesFieldsFromJson) {
      cxSmilesFields =
          static_cast<SmilesWrite::CXSmilesFields>(cxSmilesFieldsFromJson);
    }
    std::string restoreBondDirOption;
    restoreBondDirOption = pt.get("restoreBondDirOption", restoreBondDirOption);
    auto it = restoreBondDirOptionKeyValuePairs.find(restoreBondDirOption);
    if (it != restoreBondDirOptionKeyValuePairs.end()) {
      restoreBondDirs = it->second;
    }
  }
}

void updateCXSmilesFieldsFromJSON(SmilesWrite::CXSmilesFields &cxSmilesFields,
                                  RestoreBondDirOption &restoreBondDirs,
                                  const std::string &details_json) {
  updateCXSmilesFieldsFromJSON(cxSmilesFields, restoreBondDirs,
                               details_json.c_str());
}

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

std::string GetAtomSmiles(const Atom *atom, const SmilesWriteParams &params) {
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
  if (params.doIsomericSmiles) {
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        !atom->hasProp(common_properties::_brokenChirality)) {
      atString = getAtomChiralityInfo(atom);
    }
  }
  if (!params.allHsExplicit && inOrganicSubset(num)) {
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
    } else if (params.doIsomericSmiles && (isotope || atString != "")) {
      needsBracket = true;
    }
  } else {
    needsBracket = true;
  }
  if (needsBracket) {
    res += "[";
  }

  if (isotope && params.doIsomericSmiles) {
    res += std::to_string(isotope);
  }
  // this was originally only done for the organic subset,
  // applying it to other atom-types is a fix for Issue 3152751:
  // Only accept for atom->getAtomicNum() in [5, 6, 7, 8, 14, 15, 16, 33, 34, 52]
  if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' && symb[0] <= 'Z') {
    switch (atom->getAtomicNum()) {
      case 5:
      case 6:
      case 7:
      case 8:
      case 14:
      case 15:
      case 16:
      case 33:
      case 34:
      case 52:
        symb[0] -= ('A' - 'a');
    }
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

std::string GetBondSmiles(const Bond *bond, const SmilesWriteParams &params,
                          int atomToLeftIdx) {
  PRECONDITION(bond, "bad bond");
  if (atomToLeftIdx < 0) {
    atomToLeftIdx = bond->getBeginAtomIdx();
  }

  std::string res = "";
  bool aromatic = false;
  if (!params.doKekule && (bond->getBondType() == Bond::SINGLE ||
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
            if (params.allBondsExplicit || params.doIsomericSmiles) {
              res = "\\";
            }
            break;
          case Bond::ENDUPRIGHT:
            if (params.allBondsExplicit || params.doIsomericSmiles) {
              res = "/";
            }
            break;
          default:
            if (params.allBondsExplicit) {
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
        if (params.allBondsExplicit) {
          res = "-";
        } else if (aromatic && !bond->getIsAromatic()) {
          res = "-";
        }
      }
      break;
    case Bond::DOUBLE:
      // see note above
      if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
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
            if (params.allBondsExplicit || params.doIsomericSmiles) {
              res = "\\";
            }
            break;
          case Bond::ENDUPRIGHT:
            if (params.allBondsExplicit || params.doIsomericSmiles) {
              res = "/";
            }
            break;
          default:
            if (params.allBondsExplicit || !aromatic) {
              res = ":";
            }
            break;
        }
      } else if (params.allBondsExplicit || !aromatic) {
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
          res << GetAtomSmiles(mSE.obj.atom, params);
        } else {
          res << (*atomSymbols)[mSE.obj.atom->getIdx()];
        }
        atomOrdering.push_back(mSE.obj.atom->getIdx());
        break;
      case Canon::MOL_STACK_BOND:
        bond = mSE.obj.bond;
        // std::cout << "\t\tBond: " << bond->getIdx() << std::endl;
        if (!bondSymbols) {
          res << GetBondSmiles(bond, params, mSE.number);
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
        MolOps::assignStereochemistry(*tmol, params.cleanStereo);
      }
    }
    if (!doingCXSmiles) {
      // remove any stereo groups that may be present. Otherwise they will be
      // used in the canonicalization
      std::vector<StereoGroup> noStereoGroups;
      tmol->setStereoGroups(noStereoGroups);
      // remove any wiggle bonds, unspecified double bond stereochemistry, or
      // dative bonds (if we aren't doing dative bonds in the standard SMILES)
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

    if (doingCXSmiles || !params.includeDativeBonds) {
      // do not output dative bonds in the SMILES if we are doing CXSmiles (we
      // output coordinate bonds there) or if the flag is set to ignore them
      for (auto bond : tmol->bonds()) {
        if (bond->getBondType() == Bond::DATIVE) {
          // we are intentionally only handling DATIVE here. The other weird
          // RDKit dative alternatives really shouldn't ever show up.
          bond->setBondType(Bond::SINGLE);
          // update the explicit valence of the begin atom since the implicit
          // valence will no longer be properly perceived
          bond->getBeginAtom()->calcExplicitValence(false);
        }
      }
    }

#if 0
    std::cout << "----------------------------" << std::endl;
    std::cout << "MolToSmiles:" << std::endl;
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
        if (params.useStereoToBreakTies) {
          // get the ranking WITHOUT stereochemistry, and save them
          // this is used in RigourousEnhancedStereo to make sure all
          // smiles generatesd have the same basic smiles structure
          // except tot stereochemistry
          Canon::rankMolAtoms(*tmol, ranks, false, false, true, true, false);
          for (auto atom : tmol->atoms()) {
            atom->setProp("_canonicalRankingNumber", ranks[atom->getIdx()]);
          }
        }
        bool breakTies = true;
        Canon::rankMolAtoms(*tmol, ranks, breakTies, params.doIsomericSmiles,
                            params.doIsomericSmiles, true,
                            params.useStereoToBreakTies);
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
      vit = fragsMolAtomMapping[fragIdx][vit];  // Lookup the Id in the
                                                // original molecule
    }
    allAtomOrdering.push_back(atomOrdering);
    for (unsigned int &vit : bondOrdering) {
      vit = fragsMolBondMapping[fragIdx][vit];  // Lookup the Id in the
                                                // original molecule
    }
    allBondOrdering.push_back(bondOrdering);
  }

  std::string result;
  std::vector<unsigned int> flattenedAtomOrdering;
  flattenedAtomOrdering.reserve(mol.getNumAtoms());
  std::vector<unsigned int> flattenedBondOrdering;
  flattenedBondOrdering.reserve(mol.getNumBonds());
  if (params.canonical) {
    // Sort the vfragsmi, but also sort the atom and bond order vectors into
    // the same order
    typedef std::tuple<std::string, std::vector<unsigned int>,
                       std::vector<unsigned int>>
        tplType;
    std::vector<tplType> tmp(vfragsmi.size());
    for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti) {
      tmp[ti] = std::make_tuple(vfragsmi[ti], allAtomOrdering[ti],
                                allBondOrdering[ti]);
    }

    if (!params.doNotSortFragments) {
      std::sort(tmp.begin(), tmp.end());
    }

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
std::string MolToCXSmiles_internal(const ROMol &romol,
                                   const SmilesWriteParams &params,
                                   std::uint32_t flags,
                                   RestoreBondDirOption restoreBondDirs) {
  std::unique_ptr<RWMol> trwmol(new RWMol(romol));

  bool doingCXSmiles = true;

  auto res = SmilesWrite::detail::MolToSmiles(*trwmol, params, doingCXSmiles);
  if (res.empty()) {
    return res;
  }

  std::vector<unsigned int> atomOrder =
      trwmol->getProp<std::vector<unsigned int>>(
          common_properties::_smilesAtomOutputOrder);
  romol.setProp(common_properties::_smilesAtomOutputOrder, atomOrder, true);

  std::vector<unsigned int> bondOrder =
      trwmol->getProp<std::vector<unsigned int>>(
          common_properties::_smilesBondOutputOrder);
  romol.setProp(common_properties::_smilesBondOutputOrder, bondOrder, true);

  if (restoreBondDirs == RestoreBondDirOptionTrue) {
    RDKit::Chirality::reapplyMolBlockWedging(*trwmol);
  } else if (restoreBondDirs == RestoreBondDirOptionClear) {
    for (auto bond : trwmol->bonds()) {
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

  if (params.cleanStereo) {
    if (trwmol->needsUpdatePropertyCache()) {
      trwmol->updatePropertyCache();
    }
    MolOps::assignStereochemistry(*trwmol, true);
    Chirality::cleanupStereoGroups(*trwmol);
  }

  auto cxext = SmilesWrite::getCXExtensions(*trwmol, flags);
  if (!cxext.empty()) {
    res += " " + cxext;
  }
  return res;
}

bool doesAtomChiralityVary(std::vector<std::unique_ptr<RWMol>> &allMols,
                           unsigned int atomIndex) {
  PRECONDITION(allMols.size() != 0, "bad allMols size");

  auto firstChiralVal = allMols[0]->getAtomWithIdx(atomIndex)->getChiralTag();

  for (auto &oneMol : allMols) {
    auto atom = oneMol->getAtomWithIdx(atomIndex);
    if (firstChiralVal != atom->getChiralTag()) {
      return true;
    }
  }

  return false;
}

bool doesBondStereoVary(std::vector<std::unique_ptr<RWMol>> &allMols,
                        unsigned int bondIndex) {
  PRECONDITION(allMols.size() != 0, "bad allMols size");

  auto firstStereoVal = allMols[0]->getBondWithIdx(bondIndex)->getStereo();

  for (auto &oneMol : allMols) {
    auto bond = oneMol->getBondWithIdx(bondIndex);
    if (firstStereoVal != bond->getStereo()) {
      return true;
    }
  }

  return false;
}

bool doTwoAtomsVaryTheSame(std::vector<std::unique_ptr<RWMol>> &allMols,
                           unsigned int atomIndex1, unsigned int atomIndex2) {
  PRECONDITION(allMols.size() != 0, "bad allMols size");

  auto firstChiralVal1 = allMols[0]->getAtomWithIdx(atomIndex1)->getChiralTag();
  auto firstChiralVal2 = allMols[0]->getAtomWithIdx(atomIndex2)->getChiralTag();

  for (auto &oneMol : allMols) {
    auto atom1 = oneMol->getAtomWithIdx(atomIndex1);
    auto atom2 = oneMol->getAtomWithIdx(atomIndex2);
    if ((firstChiralVal1 == atom1->getChiralTag()) !=
        (firstChiralVal2 == atom2->getChiralTag())) {
      return false;
    }
  }

  return true;
}

bool doAtomAndBondVaryTheSame(std::vector<std::unique_ptr<RWMol>> &allMols,
                              unsigned int atomIndex1,
                              unsigned int bondIndex2) {
  PRECONDITION(allMols.size() != 0, "bad allMols size");

  auto firstChiralVal1 = allMols[0]->getAtomWithIdx(atomIndex1)->getChiralTag();
  auto firstStereoVal2 = allMols[0]->getBondWithIdx(bondIndex2)->getStereo();

  for (auto &oneMol : allMols) {
    auto atom1 = oneMol->getAtomWithIdx(atomIndex1);
    auto bond2 = oneMol->getBondWithIdx(bondIndex2);
    if ((firstChiralVal1 == atom1->getChiralTag()) !=
        (firstStereoVal2 == bond2->getStereo())) {
      return false;
    }
  }

  return true;
}

bool doTwoBondsVaryTheSame(std::vector<std::unique_ptr<RWMol>> &allMols,
                           unsigned int bondIndex1, unsigned int bondIndex2) {
  PRECONDITION(allMols.size() != 0, "bad allMols size");

  auto firstStereoVal1 = allMols[0]->getBondWithIdx(bondIndex1)->getStereo();
  auto firstStereoVal2 = allMols[0]->getBondWithIdx(bondIndex2)->getStereo();

  for (auto &oneMol : allMols) {
    auto bond1 = oneMol->getBondWithIdx(bondIndex1);
    auto bond2 = oneMol->getBondWithIdx(bondIndex2);
    if ((firstStereoVal1 == bond1->getStereo()) !=
        (firstStereoVal2 == bond2->getStereo())) {
      return false;
    }
  }

  return true;
}

std::unique_ptr<RWMol> canonicalizeStereoGroups_internal(
    const std::unique_ptr<RWMol> &mol, SmilesWriteParams params,
    std::uint32_t flags, RestoreBondDirOption restoreBondDirs,
    const std::vector<RDKit::StereoGroup> *groupsToProcess,
    const std::vector<RDKit::StereoGroup> *groupsToKeep,
    bool saveNewAbsoluteGroups) {
  // this expanded a mol with stereo groups to a vector of mols that have no
  // stereo groups
  //

  std::set<std::string> allSmiles;
  auto stereoGroupType = (*groupsToProcess)[0].getGroupType();
  SmilesWriteParams wps;
  wps.rigorousEnhancedStereo = false;  // avoid innfinite loop
  wps.canonical = true;
  wps.cleanStereo = false;
  wps.useStereoToBreakTies = true;
  wps.doNotSortFragments = true;

  params.rigorousEnhancedStereo = false;  // avoid infinite loop
  params.canonical = true;
  params.cleanStereo = false;
  params.useStereoToBreakTies = true;

  SmilesParserParams rps;
  rps.sanitize = false;
  rps.removeHs = false;

  std::string firstSmiles;
  auto protoMol = std::unique_ptr<RWMol>(new RWMol(*mol));
  auto protoGroups = std::vector<RDKit::StereoGroup>();
  if (groupsToKeep != nullptr) {
    for (auto grp : *groupsToKeep) {
      std::vector<Atom *> atomsToAdd;
      std::vector<Bond *> bondsToAdd;
      for (auto atomPtr : grp.getAtoms()) {
        atomsToAdd.push_back(protoMol->getAtomWithIdx(atomPtr->getIdx()));
      }
      for (auto bondPtr : grp.getBonds()) {
        bondsToAdd.push_back(protoMol->getBondWithIdx(bondPtr->getIdx()));
      }

      protoGroups.push_back(StereoGroup(grp.getGroupType(), atomsToAdd,
                                        bondsToAdd, grp.getReadId()));
    }
  }

  protoMol->setStereoGroups(protoGroups);
  std::unique_ptr<RWMol> bestNewMol;
  auto newMolCount = std::pow(2, groupsToProcess->size());
  for (unsigned int molIndex = 0; molIndex < newMolCount; ++molIndex) {
    auto newMol = std::unique_ptr<RWMol>(new RWMol(*(protoMol.get())));

    for (unsigned int grpIndex = 0; grpIndex < groupsToProcess->size();
         ++grpIndex) {
      if (molIndex & (1 << grpIndex)) {
        for (auto atomPtr : (*groupsToProcess).at(grpIndex).getAtoms()) {
          if (atomPtr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
            newMol->getAtomWithIdx(atomPtr->getIdx())
                ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
          } else if (atomPtr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
            newMol->getAtomWithIdx(atomPtr->getIdx())
                ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
          }
        }
        // do any  atropisomer bonds in this stereo group

        for (auto bond : groupsToProcess->at(grpIndex).getBonds()) {
          if (bond->getStereo() == Bond::STEREOATROPCW) {
            newMol->getBondWithIdx(bond->getIdx())
                ->setStereo(Bond::Bond::STEREOATROPCCW);
          } else if (bond->getStereo() == Bond::STEREOATROPCCW) {
            newMol->getBondWithIdx(bond->getIdx())
                ->setStereo(Bond::Bond::STEREOATROPCW);
          }
        }
      }
    }

    std::string newSmi;
    if (groupsToKeep != nullptr && groupsToKeep->size() > 0) {
      newMol = SmilesWrite::detail::canonicalizeStereoGroups_internal(
          newMol, params, flags, restoreBondDirs, &newMol->getStereoGroups(),
          nullptr, false);
    }

    newSmi = MolToCXSmiles_internal(
        *newMol.get(), wps, SmilesWrite::CX_ALL_BUT_COORDS,
        RestoreBondDirOption::RestoreBondDirOptionClear);

    auto insertResult = allSmiles.insert(newSmi);
    if (insertResult.second && newSmi == *allSmiles.begin()) {
      bestNewMol = std::move(newMol);
    }

    boost::erase_all(newSmi, "@");
    boost::erase_all(newSmi, "/");
    boost::erase_all(newSmi, "\\");
    boost::erase_all(newSmi, "[");
    boost::erase_all(newSmi, "H3]");
    boost::erase_all(newSmi, "H2]");
    boost::erase_all(newSmi, "H]");

    auto pos = newSmi.find(" |");
    if (pos) {
      newSmi = newSmi.substr(0, pos);
    }
    if (firstSmiles.size() == 0) {
      firstSmiles = newSmi;

    } else {
      if (firstSmiles != newSmi) {
        throw RigorousEnhancedStereoException("smiles are not the same");
      }
    }
  }

  // get the mol to return, including all info from the original mol

  params.doNotSortFragments = true;
  auto bestFullSmiles =
      // MolToCXSmiles(*bestNewMol, params, flags, restoreBondDirs);
      MolToCXSmiles(*bestNewMol, params, flags, RestoreBondDirOptionClear);
  params.doNotSortFragments = false;
  std::unique_ptr<RWMol> molToReturn(SmilesToMol(bestFullSmiles, rps));
  // first get all atoms and bonds in stereo groups - these will not be
  // checked in this pass.   In the second pass, there should be none also
  // while we are at it, copy the stereo groups to the new stereoGroups vector

  std::vector<unsigned int> atomIndices;
  std::vector<unsigned int> bondIndices;
  std::vector<StereoGroup> newGroups;

  std::vector<unsigned int> atomIndicesInStereoGroups;
  std::vector<unsigned int> bondIndicesInStereoGroups;

  for (auto grp : molToReturn->getStereoGroups()) {
    auto newGroup = StereoGroup(grp);
    newGroups.push_back(newGroup);

    for (auto atomPtr : grp.getAtoms()) {
      atomIndicesInStereoGroups.push_back(atomPtr->getIdx());
    }
    for (auto bondPtr : grp.getBonds()) {
      bondIndicesInStereoGroups.push_back(bondPtr->getIdx());
    }
  }

  // now get all stereo centers and atrop bonds that are not in the stereo
  // groups

  // for (auto atom : bestMol->atoms()) {
  for (auto atom : molToReturn->atoms()) {
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        std::find(atomIndicesInStereoGroups.begin(),
                  atomIndicesInStereoGroups.end(),
                  atom->getIdx()) == atomIndicesInStereoGroups.end()) {
      atomIndices.push_back(atom->getIdx());
    }
  }

  for (auto bond : molToReturn->bonds()) {
    if (bond->getStereo() != Bond::STEREONONE &&
        std::find(bondIndicesInStereoGroups.begin(),
                  bondIndicesInStereoGroups.end(),
                  bond->getIdx()) == bondIndicesInStereoGroups.end()) {
      bondIndices.push_back(bond->getIdx());
    }
  }

  atomIndicesInStereoGroups.clear();  // not needed past here
  bondIndicesInStereoGroups.clear();  // not needed past here

  std::vector<bool> atomsDone(atomIndices.size(), false);
  std::vector<bool> bondsDone(bondIndices.size(), false);

  std::vector<Atom *> absGroupAtoms;
  std::vector<Bond *> absGroupBonds;

  // if there is only one smiles, then there is no variation and
  // add stereo groups are actual abs (and there is only one group)

  if (allSmiles.size() == 1) {
    if (!saveNewAbsoluteGroups) {
      return molToReturn;
    }

    for (unsigned int atomIndex : atomIndices) {
      absGroupAtoms.push_back(molToReturn->getAtomWithIdx(atomIndex));
    }
    for (unsigned int bondIndex : bondIndices) {
      absGroupBonds.push_back(molToReturn->getBondWithIdx(bondIndex));
    }

  } else {
    // Now make the new stereo-enhanced canonical smiles

    std::vector<std::unique_ptr<RWMol>> allMols;
    for (auto smi : allSmiles) {
      allMols.push_back(std::unique_ptr<RWMol>(SmilesToMol(smi, rps)));
    }

    unsigned int groupCount = 0;

    // set up the one ABS group

    for (unsigned int index1 = 0; index1 < atomIndices.size(); ++index1) {
      if (atomsDone[index1]) {
        continue;
      }

      unsigned int atomIndex1 = atomIndices[index1];

      if (!doesAtomChiralityVary(allMols, atomIndex1)) {
        if (saveNewAbsoluteGroups) {
          absGroupAtoms.push_back(molToReturn->getAtomWithIdx(atomIndex1));
        }
        atomsDone[index1] = true;
        continue;
      }

      std::vector<Atom *> atomsToAdd;
      std::vector<Bond *> bondsToAdd;
      atomsToAdd.push_back(molToReturn->getAtomWithIdx(atomIndex1));
      atomsDone[index1] = true;

      // now look through all other possible atoms and bonds to see if they
      // vary the same way as the first one in the group

      for (unsigned int index2 = index1 + 1; index2 < atomIndices.size();
           ++index2) {
        if (atomsDone[index2]) {
          continue;
        }
        unsigned int atomIndex2 = atomIndices[index2];

        if (doTwoAtomsVaryTheSame(allMols, atomIndex1, atomIndex2)) {
          atomsToAdd.push_back(molToReturn->getAtomWithIdx(atomIndex2));
          atomsDone[index2] = true;
        }
      }

      for (unsigned int index2 = 0; index2 < bondIndices.size(); ++index2) {
        if (bondsDone[index2]) {
          continue;
        }

        auto bondIndex2 = bondIndices[index2];

        if (doAtomAndBondVaryTheSame(allMols, atomIndex1, bondIndex2)) {
          bondsToAdd.push_back(molToReturn->getBondWithIdx(bondIndex2));
          bondsDone[index2] = true;
        }
      }

      std::sort(atomsToAdd.begin(), atomsToAdd.end(),
                [](Atom *a, Atom *b) { return a->getIdx() < b->getIdx(); });
      std::sort(bondsToAdd.begin(), bondsToAdd.end(),
                [](Bond *a, Bond *b) { return a->getIdx() < b->getIdx(); });
      newGroups.emplace_back(stereoGroupType, atomsToAdd, bondsToAdd,
                             ++groupCount);
    }

    // now any groups that only involve bonds

    for (unsigned int index1 = 0; index1 < bondIndices.size(); ++index1) {
      if (bondsDone[index1]) {
        continue;
      }

      if (!doesBondStereoVary(allMols, index1)) {
        if (saveNewAbsoluteGroups) {
          absGroupBonds.push_back(
              molToReturn->getBondWithIdx(bondIndices[index1]));
        }
        bondsDone[index1] = true;
        continue;
      }

      std::vector<Atom *> atomsToAdd;  // nothing added to this one here
      std::vector<Bond *> bondsToAdd;
      bondsToAdd.push_back(molToReturn->getBondWithIdx(bondIndices[index1]));
      bondsDone[index1] = true;

      // now look through all other possible bonds to see if they vary
      // the same way as the first one in the group

      for (unsigned int index2 = index1 + 1; index2 < bondIndices.size();
           ++index2) {
        if (bondsDone[index2]) {
          continue;
        }
        unsigned int bondIndex2 = bondIndices[index2];

        if (doTwoBondsVaryTheSame(allMols, bondIndices[index1], bondIndex2)) {
          bondsToAdd.push_back(molToReturn->getBondWithIdx(bondIndex2));
          bondsDone[index2] = true;
        }
      }

      std::sort(atomsToAdd.begin(), atomsToAdd.end(),
                [](Atom *a, Atom *b) { return a->getIdx() < b->getIdx(); });
      std::sort(bondsToAdd.begin(), bondsToAdd.end(),
                [](Bond *a, Bond *b) { return a->getIdx() < b->getIdx(); });
      newGroups.emplace_back(stereoGroupType, atomsToAdd, bondsToAdd,
                             ++groupCount);
    }
  }
  // if the abs group is not empty, add it

  if (saveNewAbsoluteGroups &&
      (absGroupAtoms.size() != 0 || absGroupBonds.size() != 0)) {
    std::sort(absGroupAtoms.begin(), absGroupAtoms.end(),
              [](Atom *a, Atom *b) { return a->getIdx() < b->getIdx(); });
    std::sort(absGroupBonds.begin(), absGroupBonds.end(),
              [](Bond *a, Bond *b) { return a->getIdx() < b->getIdx(); });

    newGroups.emplace_back(StereoGroupType::STEREO_ABSOLUTE, absGroupAtoms,
                           absGroupBonds, 0);
  }

  molToReturn->setStereoGroups(newGroups);

  return molToReturn;
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
  std::unique_ptr<RWMol> trwmol(new RWMol(romol));

  if (params.canonical && params.rigorousEnhancedStereo) {
    return canonicalizeStereoGroups(trwmol, params, flags, restoreBondDirs);
  }

  return SmilesWrite::detail::MolToCXSmiles_internal(romol, params, flags,
                                                     restoreBondDirs);
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

namespace {

enum class ChiralItemType {
  ATOM,
  BOND,
};
class ChiralItem {
 public:
  ChiralItemType chiralAtomType;
  unsigned int id;
};

}  // namespace

std::string canonicalizeStereoGroups(const std::unique_ptr<RWMol> &mol,
                                     const SmilesWriteParams &params,
                                     std::uint32_t flags,
                                     RestoreBondDirOption restoreBondDirs) {
  // this expanded a mol with stereo groups to a vector of mols that have no
  // stereo grouos

  if (mol->getStereoGroups().empty()) {
    return SmilesWrite::detail::MolToCXSmiles_internal(*mol, params, flags,
                                                       restoreBondDirs);
  }

  // if there is only one group and it is absolute,
  // simply return

  if (mol.get()->getStereoGroups().size() == 1 &&
      mol.get()->getStereoGroups()[0].getGroupType() ==
          StereoGroupType::STEREO_ABSOLUTE) {
    return SmilesWrite::detail::MolToCXSmiles_internal(*mol, params, flags,
                                                       restoreBondDirs);
  }

  std::vector<RDKit::StereoGroup> orGroups;  // empty vector of new groups
  std::vector<RDKit::StereoGroup> andGroups;
  for (auto &stg : mol->getStereoGroups()) {
    if (stg.getGroupType() == StereoGroupType::STEREO_OR) {
      orGroups.push_back(stg);
    } else if (stg.getGroupType() == StereoGroupType::STEREO_AND) {
      andGroups.push_back(stg);
    }
  }

  try {
    std::unique_ptr<RWMol> canonMol;
    if (orGroups.size() == 0) {
      canonMol = SmilesWrite::detail::canonicalizeStereoGroups_internal(
          mol, params, flags, restoreBondDirs, &andGroups, nullptr, true);
    } else {
      canonMol = SmilesWrite::detail::canonicalizeStereoGroups_internal(
          mol, params, flags, restoreBondDirs, &orGroups, &andGroups, true);
    }

    return SmilesWrite::detail::MolToCXSmiles_internal(*canonMol, params, flags,
                                                       restoreBondDirs);
  } catch (const RigorousEnhancedStereoException &e) {
    SmilesWriteParams newParams(params);
    newParams.rigorousEnhancedStereo = false;
    return SmilesWrite::detail::MolToCXSmiles_internal(*mol, newParams, flags,
                                                       restoreBondDirs);
  }
}

}  // namespace RDKit
