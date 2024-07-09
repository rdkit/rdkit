//
//  Copyright (C) 2006-2021 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>

namespace RDKit {
namespace details {
void updateSubMolConfs(const ROMol &mol, RWMol &res,
                       boost::dynamic_bitset<> &removedAtoms) {
  // update conformer information:
  res.clearConformers();
  for (auto citer = mol.beginConformers(); citer != mol.endConformers();
       ++citer) {
    auto *newConf = new Conformer(res.getNumAtoms());
    newConf->setId((*citer)->getId());
    newConf->set3D((*citer)->is3D());
    int aIdx = 0;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!removedAtoms[i]) {
        newConf->setAtomPos(aIdx, (*citer)->getAtomPos(i));
        ++aIdx;
      }
    }
    res.addConformer(newConf, false);
  }
}
void copyStereoGroups(const std::map<const Atom *, Atom *> &molAtomMap,
                      const ROMol &mol, RWMol &newMol) {
  // Copy over any stereo groups that lie in the new molecule
  if (!mol.getStereoGroups().empty()) {
    std::vector<StereoGroup> newStereoGroups;
    for (auto &stereoGroup : mol.getStereoGroups()) {
      std::vector<Atom *> newStereoAtoms;
      std::vector<Bond *> newStereoBonds;
      for (const auto stereoGroupAtom : stereoGroup.getAtoms()) {
        if (auto found = molAtomMap.find(stereoGroupAtom);
            found != molAtomMap.end()) {
          newStereoAtoms.push_back(found->second);
        }
      }
      for (const auto stereoGroupBond : stereoGroup.getBonds()) {
        auto foundFirst = molAtomMap.find(stereoGroupBond->getBeginAtom());
        auto foundSecond = molAtomMap.find(stereoGroupBond->getEndAtom());

        if (foundFirst != molAtomMap.end() && foundSecond != molAtomMap.end()) {
          newStereoBonds.push_back(newMol.getBondBetweenAtoms(
              foundFirst->second->getIdx(), foundSecond->second->getIdx()));
        }
      }
      if (!newStereoAtoms.empty()) {
        newStereoGroups.emplace_back(stereoGroup.getGroupType(), newStereoAtoms,
                                     newStereoBonds, stereoGroup.getReadId());
      }
    }
    newMol.setStereoGroups(std::move(newStereoGroups));
  }
}
}  // namespace details

namespace {
struct SideChainMapping {
  int molIndex;
  int coreIndex;
  bool useMatch;

  SideChainMapping(int molIndex)
      : molIndex(molIndex), coreIndex(-1), useMatch(false) {}

  SideChainMapping(int molIndex, int coreIndex, bool useMatch)
      : molIndex(molIndex), coreIndex(coreIndex), useMatch(useMatch) {}
};
}  // namespace

ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query, bool onlyFrags,
                        bool useChirality) {
  auto *res = new RWMol(mol, false);
  std::vector<MatchVectType> fgpMatches;
  // do the substructure matching and get the atoms that match the query
  const bool uniquify = true;
  const bool recursionPossible = true;
  SubstructMatch(*res, query, fgpMatches, uniquify, recursionPossible,
                 useChirality);

  // if didn't find any matches nothing to be done here
  // simply return a copy of the molecule
  if (fgpMatches.empty()) {
    return res;
  }

  // all matches on the molecule - list of list of atom ids
  VECT_INT_VECT matches;
  matches.reserve(fgpMatches.size());
  for (const auto& mati : fgpMatches) {
    INT_VECT match;  // each match onto the molecule - list of atoms ids
    match.reserve(mati.size());
    for (const auto&  mi : mati) {
      match.push_back(mi.second);
    }
    matches.push_back(std::move(match));
  }

  // now loop over the list of matches and check if we can delete any of them
  INT_VECT delList;
  if (onlyFrags) {
    VECT_INT_VECT frags;
    MolOps::getMolFrags(*res, frags);
    for (auto& fi : frags) {
      std::sort(fi.begin(), fi.end());
      for (auto& mxi : matches) {
        std::sort(mxi.begin(), mxi.end());
        if (fi == mxi) {
          INT_VECT tmp;
          Union(mxi, delList, tmp);
          delList = tmp;
          break;
        }  // end of if we found a matching fragment
      }    // end of loop over matches
    }      // end of loop over fragments
  } else {
    // in this case we want to delete any matches we find
    // simply loop over the matches and collect the atoms that need to
    // be removed
    for (const auto& mxi : matches) {
      INT_VECT tmp;
      Union(mxi, delList, tmp);
      delList = tmp;
    }
  }

  if (delList.empty()) {
    return res;
  }

  // now loop over the union list and delete the atoms
  res->beginBatchEdit();
  boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
  for (auto idx : delList) {
    removedAtoms.set(idx);
    res->removeAtom(idx);
  }
  res->commitBatchEdit();

  details::updateSubMolConfs(mol, *res, removedAtoms);

  res->clearComputedProps(true);
  // update our properties, but allow unhappiness:
  res->updatePropertyCache(false);

  return res;
}

std::vector<ROMOL_SPTR> replaceSubstructs(
    const ROMol &mol, const ROMol &query, const ROMol &replacement,
    bool replaceAll, unsigned int replacementConnectionPoint,
    bool useChirality) {
  PRECONDITION(replacementConnectionPoint < replacement.getNumAtoms(),
               "bad replacementConnectionPoint");
  std::vector<ROMOL_SPTR> res;
  std::vector<MatchVectType> fgpMatches;

  boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
  // do the substructure matching and get the atoms that match the query
  const bool uniquify = true;
  const bool recursionPossible = true;
  SubstructMatch(mol, query, fgpMatches, uniquify, recursionPossible,
                 useChirality);

  // if we didn't find any matches, there's nothing to be done here
  // simply return a list with a copy of the starting molecule
  if (fgpMatches.size() == 0) {
    res.push_back(ROMOL_SPTR(new ROMol(mol, false)));
    res[0]->clearComputedProps(false);
    return res;
  }

  INT_VECT delList;

  // now loop over the list of matches and replace them:
  for (const auto &fgpMatche : fgpMatches) {
    INT_VECT match;  // each match onto the molecule - list of atoms ids
    for (const auto &mi : fgpMatche) {
      match.push_back(mi.second);
    }

    INT_VECT sortMatch = match;
    std::sort(sortMatch.begin(), sortMatch.end());

    if (!replaceAll || !res.size()) {
      res.push_back(ROMOL_SPTR(new ROMol(mol, false)));
    }
    RWMol *newMol = static_cast<RWMol *>(res.rbegin()->get());

    // we need a tab to the orig number of atoms because the
    // new molecule will start numbered above this:
    int numOrigAtoms = newMol->getNumAtoms();

    // Add the atoms and bonds from the replacement:
    newMol->insertMol(replacement);

    // loop over the central atom's (the first atom in match) bonds
    // and duplicate any that connect to the remainder of the molecule:
    Atom *origAtom = newMol->getAtomWithIdx(match[0]);
    ROMol::ADJ_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) = newMol->getAtomNeighbors(origAtom);
    while (nbrIdx != endNbrs) {
      // we don't want to duplicate any "intra-match" bonds:
      if (!std::binary_search(sortMatch.begin(), sortMatch.end(),
                              int(*nbrIdx))) {
        Bond *oBond = newMol->getBondBetweenAtoms(match[0], *nbrIdx);
        CHECK_INVARIANT(oBond, "required bond not found");
        newMol->addBond(numOrigAtoms + replacementConnectionPoint, *nbrIdx,
                        oBond->getBondType());
      }
      nbrIdx++;
    }

    if (replaceAll) {
      // we'll accumulate  a list of atoms to be removed:
      INT_VECT tmp;
      Union(sortMatch, delList, tmp);
      delList = tmp;
    } else {
      // just delete the atoms now:
      newMol->beginBatchEdit();
      for (auto idx : sortMatch) {
        removedAtoms.set(idx);
        newMol->removeAtom(idx);
      }
      newMol->commitBatchEdit();
    }
  }

  if (delList.size()) {
    if (replaceAll) {
      // remove the atoms from the delList:
      auto *newMol = static_cast<RWMol *>(res[0].get());
      newMol->beginBatchEdit();
      for (auto idx : delList) {
        removedAtoms.set(idx);
        newMol->removeAtom(idx);
      }
      newMol->commitBatchEdit();
    }

    // clear conformers and computed props and do basic updates
    // on the resulting molecules, but allow unhappiness:
    for (auto &re : res) {
      details::updateSubMolConfs(mol, *(RWMol *)re.get(), removedAtoms);
      re->clearComputedProps(true);
      re->updatePropertyCache(false);
    }
  }
  return res;
}

ROMol *replaceSidechains(const ROMol &mol, const ROMol &coreQuery,
                         bool useChirality) {
  MatchVectType matchV;

  // do the substructure matching and get the atoms that match the query
  const bool recursionPossible = true;
  bool matchFound =
      SubstructMatch(mol, coreQuery, matchV, recursionPossible, useChirality);

  // if we didn't find any matches, there's nothing to be done here
  // simply return null to indicate the problem
  if (!matchFound || matchV.size() == 0) {
    return nullptr;
  }

  boost::dynamic_bitset<> matchingIndices(mol.getNumAtoms());
  for (auto mvit : matchV) {
    matchingIndices[mvit.second] = 1;
  }

  auto *newMol = new RWMol(mol);
  boost::dynamic_bitset<> keepSet(newMol->getNumAtoms());
  std::vector<unsigned int> dummyIndices;
  for (auto mvit : matchV) {
    keepSet.set(mvit.second);
    // if the atom in the molecule has higher degree than the atom in the
    // core, we have an attachment point:
    if (newMol->getAtomWithIdx(mvit.second)->getDegree() >
        coreQuery.getAtomWithIdx(mvit.first)->getDegree()) {
      ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) =
          newMol->getAtomNeighbors(newMol->getAtomWithIdx(mvit.second));
      while (nbrIdx != endNbrs) {
        if (!matchingIndices[*nbrIdx]) {
          // this neighbor isn't in the match, convert it to a dummy atom and
          // save it
          keepSet.set(*nbrIdx);
          dummyIndices.push_back(*nbrIdx);
          Atom *at = newMol->getAtomWithIdx(*nbrIdx);
          Bond *b = newMol->getBondBetweenAtoms(mvit.second, *nbrIdx);
          if (b) {
            b->setIsAromatic(false);
            b->setBondType(Bond::SINGLE);
          }
          at->setAtomicNum(0);
          at->setNumExplicitHs(0);
          at->setIsAromatic(false);
          at->setIsotope(dummyIndices.size());
        }
        ++nbrIdx;
      }
    }
  }
  boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
  newMol->beginBatchEdit();
  for (const auto at : newMol->atoms()) {
    if (!keepSet.test(at->getIdx())) {
      newMol->removeAtom(at);
      removedAtoms.set(at->getIdx());
    }
  }
  // Remove bonds between newly added dummies (if any)
  if (dummyIndices.size() > 1) {
    for (size_t i = 0; i < dummyIndices.size() - 1; ++i) {
      for (size_t j = i + 1; j < dummyIndices.size(); ++j) {
        const auto b =
            newMol->getBondBetweenAtoms(dummyIndices.at(i), dummyIndices.at(j));
        if (b) {
          newMol->removeBond(dummyIndices.at(i), dummyIndices.at(j));
        }
      }
    }
  }
  newMol->commitBatchEdit();

  details::updateSubMolConfs(mol, *newMol, removedAtoms);

  // clear computed props and do basic updates on the
  // the resulting molecule, but allow unhappiness:
  newMol->clearComputedProps(true);
  newMol->updatePropertyCache(false);
  return static_cast<ROMol *>(newMol);
}

ROMol *replaceCore(const ROMol &mol, const ROMol &coreQuery,
                   bool replaceDummies, bool labelByIndex,
                   bool requireDummyMatch, bool useChirality) {
  MatchVectType matchV;

  // do the substructure matching and get the atoms that match the query
  const bool recursionPossible = true;
  bool matchFound =
      SubstructMatch(mol, coreQuery, matchV, recursionPossible, useChirality);

  // if we didn't find any matches, there's nothing to be done here
  // simply return null to indicate the problem
  if (!matchFound || matchV.size() == 0) {
    return nullptr;
  }

  return replaceCore(mol, coreQuery, matchV, replaceDummies, labelByIndex,
                     requireDummyMatch);
}

namespace {
const std::string replaceCoreDummyBond = "_replaceCoreDummyBond";

int findNbrBond(RWMol &mol, Bond *bond, Atom *bondAtom, const INT_VECT &bring,
                const boost::dynamic_bitset<> &removedAtoms) {
  int res = -1;
  for (const auto nbrBond : mol.atomBonds(bondAtom)) {
    if (nbrBond != bond &&
        (nbrBond->hasProp(replaceCoreDummyBond) ||
         (!removedAtoms[nbrBond->getOtherAtomIdx(bondAtom->getIdx())] &&
          std::find(bring.begin(), bring.end(), nbrBond->getIdx()) !=
              bring.end()))) {
      res = nbrBond->getOtherAtomIdx(bondAtom->getIdx());
      break;
    }
  }
  return res;
}

void setSubMolBrokenRingStereo(RWMol &mol,
                               const boost::dynamic_bitset<> &removedAtoms) {
  PRECONDITION(mol.getRingInfo() && mol.getRingInfo()->isInitialized(),
               "bad ringinfo");

  for (const auto &bring : mol.getRingInfo()->bondRings()) {
    // check whether or not this bond ring is affected by the removal
    if (std::find_if(bring.begin(), bring.end(),
                     [&mol, &removedAtoms](auto idx) {
                       const auto bond = mol.getBondWithIdx(idx);
                       return removedAtoms[bond->getBeginAtomIdx()] ||
                              removedAtoms[bond->getEndAtomIdx()];
                     }) != bring.end()) {
      for (auto bidx : bring) {
        const auto bond = mol.getBondWithIdx(bidx);
        // is this a bond we can reasonably set cis/trans for and where neither
        // atom is being removed?
        if ((bond->getIsAromatic() ||
             bond->getBondType() == Bond::BondType::DOUBLE) &&
            bond->getStereo() == Bond::BondStereo::STEREONONE &&
            mol.getRingInfo()->minBondRingSize(bond->getIdx()) <
                Chirality::minRingSizeForDoubleBondStereo &&
            !removedAtoms[bond->getBeginAtomIdx()] &&
            !removedAtoms[bond->getEndAtomIdx()]) {
          // find the two neighboring bonds which are in the ring or to the
          // newly added dummy atoms. Make sure they aren't affected by the atom
          // removal
          int beginAtomNbrIdx =
              findNbrBond(mol, bond, bond->getBeginAtom(), bring, removedAtoms);
          if (beginAtomNbrIdx >= 0) {
            int endAtomNbrIdx =
                findNbrBond(mol, bond, bond->getEndAtom(), bring, removedAtoms);
            if (endAtomNbrIdx >= 0) {
              // we can set stereo on the bond
              bond->setStereoAtoms(beginAtomNbrIdx, endAtomNbrIdx);
              bond->setStereo(Bond::BondStereo::STEREOCIS);
            }
          }
        }
      }
    }
  }
}
}  // namespace

ROMol *replaceCore(const ROMol &mol, const ROMol &core,
                   const MatchVectType &matchV, bool replaceDummies,
                   bool labelByIndex, bool requireDummyMatch) {
  unsigned int origNumAtoms = mol.getNumAtoms();
  std::vector<std::pair<int, SideChainMapping>> matches;
  matches.reserve(origNumAtoms);

  std::vector<int> matchingIndices(origNumAtoms, -1);
  std::vector<int> allIndices(origNumAtoms, -1);
  boost::dynamic_bitset<> molAtomsMapped(origNumAtoms);
  boost::dynamic_bitset<> multipleMappedMolAtoms(origNumAtoms);
  for (const auto &mvit : matchV) {
    if (mvit.first < 0 || mvit.first >= rdcast<int>(core.getNumAtoms())) {
      throw ValueErrorException(
          "Supplied MatchVect indices out of bounds of the core molecule");
    }
    if (mvit.second < 0 || mvit.second >= rdcast<int>(mol.getNumAtoms())) {
      throw ValueErrorException(
          "Supplied MatchVect indices out of bounds of the target molecule");
    }
    bool useMatch = false;
    if (replaceDummies || core.getAtomWithIdx(mvit.first)->getAtomicNum() > 0) {
      matchingIndices[mvit.second] = mvit.first;
      useMatch = true;
    }
    allIndices[mvit.second] = mvit.first;
    SideChainMapping mapping(mvit.second, mvit.first, useMatch);
    matches.emplace_back(mvit.second, mapping);
    if (molAtomsMapped[mvit.second]) {
      multipleMappedMolAtoms.set(mvit.second);
    }
    molAtomsMapped.set(mvit.second);
  }

  boost::dynamic_bitset<> multipleOwnedBonds(mol.getNumBonds());
  if (multipleMappedMolAtoms.any()) {
    for (const auto &match : matches) {
      const auto &mappingInfo = match.second;
      if (multipleMappedMolAtoms[mappingInfo.molIndex]) {
        auto coreAtom = core.getAtomWithIdx(mappingInfo.coreIndex);
        CHECK_INVARIANT(
            coreAtom->getDegree() == 1,
            "Multiple core atoms match a mol atom, but one of the core "
            "atoms has degree > 1 ");
        auto coreNeighborIdx =
            core[*core.getAtomNeighbors(coreAtom).first]->getIdx();
        auto molNeighborIdx =
            std::find_if(matchV.cbegin(), matchV.cend(),
                         [coreNeighborIdx](std::pair<int, int> p) {
                           return p.first == static_cast<int>(coreNeighborIdx);
                         })
                ->second;
        if (molNeighborIdx > -1) {
          auto connectingBond =
              mol.getBondBetweenAtoms(mappingInfo.molIndex, molNeighborIdx);
          CHECK_INVARIANT(connectingBond,
                          "expected bond in molecule not found");
          multipleOwnedBonds.set(connectingBond->getIdx());
        }
      }
    }
  }

  auto *newMol = new RWMol(mol);
  std::vector<Atom *> keepList;
  std::map<Atom *, int> dummyAtomMap;
  std::map<const Atom *, Atom *> molAtomMap;

  // go through the matches in query order, not target molecule
  //  order
  for (unsigned int i = 0; i < origNumAtoms; ++i) {
    if (!molAtomsMapped[i]) {
      SideChainMapping mapping(i);
      matches.emplace_back(i, mapping);
    }
  }

  std::sort(matches.begin(), matches.end(),
            [](const std::pair<int, SideChainMapping> &p1,
               const std::pair<int, SideChainMapping> &p2) {
              if (p1.second.coreIndex == p2.second.coreIndex) {
                return p1.first < p2.first;
              }
              return p1.second.coreIndex < p2.second.coreIndex;
            });
  std::vector<std::pair<int, Atom *>> dummies;

  std::list<Bond *> allNewBonds;
  for (const auto &match : matches) {
    const auto &mappingInfo = match.second;

    if (!mappingInfo.useMatch) {
      Atom *sidechainAtom = newMol->getAtomWithIdx(mappingInfo.molIndex);
      // we're keeping the sidechain atoms:
      keepList.push_back(sidechainAtom);
      molAtomMap[mol.getAtomWithIdx(mappingInfo.molIndex)] = sidechainAtom;

      // loop over our neighbors and see if any are in the match:
      std::list<unsigned int> nbrList;
      ROMol::ADJ_ITER nbrIter, endNbrs;
      boost::tie(nbrIter, endNbrs) = newMol->getAtomNeighbors(sidechainAtom);
      while (nbrIter != endNbrs && (*nbrIter) < origNumAtoms) {
        // we need to add bonds and atoms to the molecule while looping
        // over neighbors. This invalidates iterators, so collect a list
        // of our neighbors now:
        nbrList.push_back(*nbrIter);
        ++nbrIter;
      }
      unsigned int whichNbr = 0;
      std::list<Bond *> newBonds;
      for (unsigned int nbrIdx : nbrList) {
        Bond *connectingBond =
            newMol->getBondBetweenAtoms(mappingInfo.molIndex, nbrIdx);
        bool bondToCore = matchingIndices[nbrIdx] > -1;
        auto coreBond =
            bondToCore && allIndices[nbrIdx] > -1 && mappingInfo.coreIndex > -1
                ? core.getBondBetweenAtoms(mappingInfo.coreIndex,
                                           allIndices[nbrIdx])
                : nullptr;
        if (bondToCore && multipleMappedMolAtoms[mappingInfo.molIndex] &&
            mappingInfo.coreIndex > -1) {
          // The core has multiple atoms that map onto this mol atom - check we
          // have matched correct core bond.
          // Otherwise we can use this bond only if nobody else owns it.
          if (coreBond == nullptr &&
              multipleOwnedBonds[connectingBond->getIdx()]) {
            bondToCore = false;
          }
        }
        if (bondToCore) {
          // we've matched an atom in the core.
          if (requireDummyMatch &&
              core.getAtomWithIdx(matchingIndices[nbrIdx])->getAtomicNum() !=
                  0) {
            delete newMol;
            return nullptr;
          }
          auto *newAt = new Atom(0);

          // if we were not in the matching list, still keep
          //  the original indices (replaceDummies=False)
          int mapping = mappingInfo.coreIndex;
          // If we don't have a core bond, the label belongs to the neighbor
          if (coreBond == nullptr) {
            mapping = allIndices[nbrIdx];
          }

          // we want to order the dummies in the same orders as
          //  the mappings, if not labelling by Index they are in arbitrary
          //  order
          //  right now so save and sort later.
          if (mapping != -1) {
            dummies.emplace_back(mapping, newAt);
          } else {
            dummies.emplace_back(matchingIndices[nbrIdx], newAt);
          }

          newMol->addAtom(newAt, false, true);
          dummyAtomMap[newAt] = nbrIdx;
          keepList.push_back(newAt);
          Bond *bnd = connectingBond->copy();
          // If the connecting bond has stereo settings those cannot be
          // preserved
          if (bnd->getStereo() > Bond::STEREOANY) {
            bnd->setStereo(Bond::STEREOANY);
          }
          if (bnd->getBeginAtomIdx() ==
              static_cast<size_t>(mappingInfo.molIndex)) {
            bnd->setEndAtomIdx(newAt->getIdx());
          } else {
            bnd->setBeginAtomIdx(newAt->getIdx());
          }
          newBonds.push_back(bnd);
          allNewBonds.push_back(bnd);

          // Check to see if we are breaking a stereo bond definition, by
          // removing one of the stereo atoms If so, set to the new atom
          for (const auto bond : newMol->atomBonds(sidechainAtom)) {
            if (bond->getIdx() == connectingBond->getIdx()) {
              continue;
            }

            if (bond->getStereo() > Bond::STEREOANY) {
              auto &stereoAtoms = bond->getStereoAtoms();
              for (int &stereoAtom : stereoAtoms) {
                if (stereoAtom == static_cast<int>(nbrIdx)) {
                  stereoAtom = static_cast<int>(newAt->getIdx());
                }
              }
            }
          }

          // we may be changing the bond ordering at the atom.
          // e.g. replacing the N in C[C@](Cl)(N)F gives an atom ordering of
          // C[C?](Cl)(F)[X]
          // so we need the SMILES C[C@@](Cl)(F)[X] to maintain the appropriate
          // chirality
          // check for these cases and adjust our chirality flags as
          // appropriate.
          //
          if (sidechainAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
              sidechainAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
            bool switchIt = false;
            switch (newMol->getAtomDegree(sidechainAtom)) {
              case 4:
                //     start:         ordering:        swap?
                //   N[C@](F)(Cl)C -> F[C@@](Cl)(C)X    yes
                //   F[C@](N)(Cl)C -> F[C@](Cl)(C)X     no
                //   F[C@](Cl)(N)C -> F[C@@](Cl)(C)X    yes
                //   F[C@](Cl)(C)N -> F[C@](Cl)(C)X     no
                if (!(whichNbr % 2)) {
                  switchIt = true;
                }
                break;
              case 3:
                // things are different in the degree three case because of the
                // implicit H:
                //     start:         ordering:     swap?
                //   N[C@H](F)C -> [C@H](F)(C)X     yes
                //   [C@H](N)(F)C -> [C@H](F)(C)X   no
                //   F[C@H](N)C -> F[C@@H](C)X      yes
                //   F[C@H](C)N -> F[C@H](C)X       no
                if (whichNbr == 1) {
                  switchIt = true;
                }
                break;
            }
            if (switchIt) {
              sidechainAtom->invertChirality();
            }
            // the neighbor count needs to be decremented, in case an additional
            // bond to the core is present at this chiral center (unlikely!),
            // because this bond will now be at the end of the neighbor list.
            --whichNbr;
          }
        }
        ++whichNbr;
      }
      // add the bonds now, after we've finished the loop over neighbors:
      for (auto &newBond : newBonds) {
        newBond->setProp(replaceCoreDummyBond, 1);
        newMol->addBond(newBond, true);
      }
    }
  }

  if (!labelByIndex) {
    // sort the mapping indices, but label from 1..N
    std::stable_sort(dummies.begin(), dummies.end());
    for (size_t nDummy = 0; nDummy < dummies.size(); ++nDummy) {
      dummies[nDummy].second->setIsotope(nDummy + 1);
    }
  } else {
    // don't sort, just label by the index
    for (auto &dummy : dummies) {
      dummy.second->setIsotope(dummy.first);
    }
  }

  std::vector<Atom *> delList;
  boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
  bool removedRingAtom = false;
  newMol->beginBatchEdit();
  for (const auto at : newMol->atoms()) {
    if (std::find(keepList.begin(), keepList.end(), at) == keepList.end()) {
      newMol->removeAtom(at);
      removedAtoms.set(at->getIdx());
      if (!removedRingAtom && mol.getRingInfo() &&
          mol.getRingInfo()->isInitialized() &&
          mol.getRingInfo()->numAtomRings(at->getIdx())) {
        removedRingAtom = true;
      }
    }
  }
  if (removedRingAtom) {
    setSubMolBrokenRingStereo(*newMol, removedAtoms);
    for (auto bond : newMol->bonds()) {
      bond->clearProp(replaceCoreDummyBond);
    }
  }
  newMol->commitBatchEdit();

  details::updateSubMolConfs(mol, *newMol, removedAtoms);

  // Update any terminal dummy atom coordinates after removing atoms not in the
  // keeplist and calling updateSubMolConfs
  for (auto &newBond : allNewBonds) {
    auto beginAtom = newBond->getBeginAtom();
    auto endAtom = newBond->getEndAtom();
    CHECK_INVARIANT(beginAtom->getDegree() == 1 || endAtom->getDegree() == 1,
                    "neither atom has degree one");
    if (newMol->getNumConformers()) {
      if (endAtom->getAtomicNum() == 0 && endAtom->getDegree() == 1) {
        MolOps::setTerminalAtomCoords(*newMol, endAtom->getIdx(),
                                      beginAtom->getIdx());
      } else {
        MolOps::setTerminalAtomCoords(*newMol, beginAtom->getIdx(),
                                      endAtom->getIdx());
      }
    }
  }

  // make a guess at the position of the dummy atoms showing the attachment
  // point:
  for (auto citer = mol.beginConformers(); citer != mol.endConformers();
       ++citer) {
    Conformer &newConf = newMol->getConformer((*citer)->getId());
    for (auto &iter : dummyAtomMap) {
      newConf.setAtomPos(iter.first->getIdx(),
                         (*citer)->getAtomPos(iter.second));
    }
  }

  // copy over stereo groups
  details::copyStereoGroups(molAtomMap, mol, *newMol);

  // clear computed props and do basic updates on
  // the resulting molecule, but allow unhappiness:
  newMol->clearComputedProps(true);
  newMol->updatePropertyCache(false);

  return static_cast<ROMol *>(newMol);
}

ROMol *MurckoDecompose(const ROMol &mol) {
  auto *res = new RWMol(mol);
  unsigned int nAtoms = res->getNumAtoms();
  if (!nAtoms) {
    return res;
  }

  // start by getting the shortest paths matrix:
  MolOps::getDistanceMat(mol, false, false, true);
  boost::shared_array<int> pathMat;
  mol.getProp(common_properties::DistanceMatrix_Paths, pathMat);

  boost::dynamic_bitset<> keepAtoms(nAtoms);
  const RingInfo *ringInfo = res->getRingInfo();
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (ringInfo->numAtomRings(i)) {
      keepAtoms[i] = 1;
    }
  }
  const VECT_INT_VECT &rings = ringInfo->atomRings();
  // std::cerr<<"  rings: "<<rings.size()<<std::endl;
  // now find the shortest paths between each ring system and mark the atoms
  // along each as being keepers:
  for (auto ringsItI = rings.begin(); ringsItI != rings.end(); ++ringsItI) {
    for (auto ringsItJ = ringsItI + 1; ringsItJ != rings.end(); ++ringsItJ) {
      int atomI = (*ringsItI)[0];
      int atomJ = (*ringsItJ)[0];
      // std::cerr<<atomI<<" -> "<<atomJ<<": ";
      while (atomI != atomJ) {
        keepAtoms[atomI] = 1;
        atomI = pathMat[atomJ * nAtoms + atomI];
        // test for the disconnected case:
        if (atomI < 0) {
          break;
        }
        // std::cerr<<atomI<<" ";
      }
      // std::cerr<<std::endl;
    }
  }

  boost::dynamic_bitset<> removedAtoms(nAtoms);
  res->beginBatchEdit();
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (!keepAtoms[i]) {
      Atom *atom = res->getAtomWithIdx(i);
      bool removeIt = true;

      // check if the atom has a neighboring keeper:
      for (auto nbr : res->atomNeighbors(atom)) {
        if (keepAtoms[nbr->getIdx()]) {
          if (res->getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
                  ->getBondType() == Bond::DOUBLE) {
            removeIt = false;
            break;
          } else if (nbr->getIsAromatic() && nbr->getAtomicNum() != 6) {
            // fix aromatic heteroatoms:
            nbr->setNumExplicitHs(1);
          } else if (nbr->getIsAromatic() && nbr->getAtomicNum() == 6 &&
                     nbr->getFormalCharge() == 1) {
            // fix aromatic carbocations
            nbr->setNumExplicitHs(1);
          } else if (nbr->getNoImplicit() ||
                     nbr->getChiralTag() != Atom::CHI_UNSPECIFIED) {
            nbr->setNoImplicit(false);
            nbr->setNumExplicitHs(0);
            nbr->setChiralTag(Atom::CHI_UNSPECIFIED);
          }
        }
      }

      if (removeIt) {
        res->removeAtom(atom);
        removedAtoms.set(atom->getIdx());
      }
    }
  }
  res->commitBatchEdit();

  details::updateSubMolConfs(mol, *res, removedAtoms);
  res->clearComputedProps();

  return (ROMol *)res;
}

ROMol *combineMols(const ROMol &mol1, const ROMol &mol2,
                   RDGeom::Point3D offset) {
  auto *res = new RWMol(mol1);
  int nAtoms1 = res->getNumAtoms();
  res->insertMol(mol2);

  // copy over coordinates
  if (mol1.getNumConformers() && mol2.getNumConformers()) {
    if (mol1.getNumConformers() != mol2.getNumConformers()) {
      BOOST_LOG(rdWarningLog)
          << "combineMols: molecules have unequal numbers of conformers"
          << std::endl;
    }
    for (auto conf1It = res->beginConformers(); conf1It != res->endConformers();
         ++conf1It) {
      Conformer *conf1 = (*conf1It).get();
      try {
        const Conformer *conf2 = &mol2.getConformer(conf1->getId());
        for (unsigned int i = 0; i < mol2.getNumAtoms(); ++i) {
          conf1->setAtomPos(i + nAtoms1, conf2->getAtomPos(i) + offset);
        }
      } catch (ConformerException &) {
        BOOST_LOG(rdWarningLog) << "combineMols: conformer id "
                                << conf1->getId() << " not found in mol2";
      }
    }
  }
  res->clearComputedProps(true);
  return (ROMol *)res;
}

void addRecursiveQueries(
    ROMol &mol, const std::map<std::string, ROMOL_SPTR> &queries,
    const std::string &propName,
    std::vector<std::pair<unsigned int, std::string>> *reactantLabels) {
  std::string delim = ",";
  boost::char_separator<char> sep(delim.c_str());
  if (reactantLabels != nullptr) {
    (*reactantLabels).resize(0);
  }

  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    Atom *at = mol[*atBegin];
    ++atBegin;
    if (!at->hasProp(propName)) {
      continue;
    }
    std::string pval;
    at->getProp(propName, pval);
    std::string maybeSmarts =
        pval;  // keep unmodified in case we are a smarts string
    boost::algorithm::to_lower(pval);
    if (reactantLabels != nullptr) {
      std::pair<unsigned int, std::string> label(at->getIdx(), pval);
      (*reactantLabels).push_back(label);
    }

    QueryAtom::QUERYATOM_QUERY *qToAdd = nullptr;
    bool notFound = false;
    if (pval.find(delim) != std::string::npos) {
      boost::tokenizer<boost::char_separator<char>> tokens(pval, sep);
      boost::tokenizer<boost::char_separator<char>>::iterator token;
      qToAdd = new ATOM_OR_QUERY();
      qToAdd->setDescription("AtomOr");

      for (token = tokens.begin(); token != tokens.end(); ++token) {
        auto iter = queries.find(*token);
        if (iter == queries.end()) {
          delete qToAdd;
          notFound = true;
          break;
        }
        auto *tqp = new RecursiveStructureQuery(new ROMol(*(iter->second)));
        std::shared_ptr<RecursiveStructureQuery> nq(tqp);
        qToAdd->addChild(nq);
      }
    } else {
      auto iter = queries.find(pval);
      if (iter == queries.end()) {
        notFound = true;
      } else {
        qToAdd = new RecursiveStructureQuery(new ROMol(*(iter->second)));
      }
    }

    if (notFound) {
      // See if we are actually a smarts expression already
      RWMol *m = nullptr;
      try {
        m = SmartsToMol(maybeSmarts);
        if (!m) {
          throw KeyErrorException(pval);
        }
        qToAdd = new RecursiveStructureQuery(m);
      } catch (...) {
        throw KeyErrorException(pval);
      }
    }

    if (!at->hasQuery()) {
      QueryAtom qAt(*at);
      unsigned int idx = at->getIdx();
      static_cast<RWMol &>(mol).replaceAtom(idx, &qAt);
      at = mol.getAtomWithIdx(idx);
    }
    at->expandQuery(qToAdd, Queries::COMPOSITE_AND);
  }
}

void parseQueryDefFile(std::istream *inStream,
                       std::map<std::string, ROMOL_SPTR> &queryDefs,
                       bool standardize, const std::string &delimiter,
                       const std::string &comment, unsigned int nameColumn,
                       unsigned int smartsColumn) {
  PRECONDITION(inStream, "no stream");
  queryDefs.clear();

  boost::char_separator<char> sep(delimiter.c_str());
  unsigned int line = 0;
  std::string tempStr;
  while (!inStream->eof() && !inStream->fail()) {
    line++;
    tempStr = getLine(inStream);
    if (tempStr == "" || tempStr.find(comment) == 0) {
      continue;
    }
    boost::tokenizer<boost::char_separator<char>> tokens(tempStr, sep);
    unsigned int tpos;
    boost::tokenizer<boost::char_separator<char>>::iterator token;
    std::string qname = "";
    std::string qval = "";
    for (token = tokens.begin(), tpos = 0; token != tokens.end();
         ++token, ++tpos) {
      if (tpos == nameColumn) {
        qname = *token;
      } else if (tpos == smartsColumn) {
        qval = *token;
      }
    }
    boost::trim_if(qname, boost::is_any_of(" \t"));
    boost::trim_if(qval, boost::is_any_of(" \t"));
    if (qname == "" || qval == "") {
      continue;
    }
    RWMol *m = nullptr;
    try {
      m = SmartsToMol(qval);
    } catch (...) {
      m = nullptr;
    }
    if (!m) {
      BOOST_LOG(rdWarningLog) << "cannot convert SMARTS " << qval
                              << " to molecule at line " << line << std::endl;
      continue;
    }
    ROMOL_SPTR msptr(m);
    if (standardize) {
      boost::algorithm::to_lower(qname);
    }
    queryDefs[qname] = msptr;
  }
}

void parseQueryDefFile(const std::string &filename,
                       std::map<std::string, ROMOL_SPTR> &queryDefs,
                       bool standardize, const std::string &delimiter,
                       const std::string &comment, unsigned int nameColumn,
                       unsigned int smartsColumn) {
  std::ifstream inStream(filename.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << filename;
    throw BadFileException(errout.str());
  }
  parseQueryDefFile(&inStream, queryDefs, standardize, delimiter, comment,
                    nameColumn, smartsColumn);
}

void parseQueryDefText(const std::string &queryDefText,
                       std::map<std::string, ROMOL_SPTR> &queryDefs,
                       bool standardize, const std::string &delimiter,
                       const std::string &comment, unsigned int nameColumn,
                       unsigned int smartsColumn) {
  std::stringstream inStream(queryDefText);
  parseQueryDefFile(&inStream, queryDefs, standardize, delimiter, comment,
                    nameColumn, smartsColumn);
}
}  // end of namespace RDKit
