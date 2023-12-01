//
//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RGroupCore.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"
#include "GraphMol/ChemTransforms/ChemTransforms.h"
#include "GraphMol/Substruct/SubstructUtils.h"
#include "GraphMol/TautomerQuery/TautomerQuery.h"

namespace RDKit {
namespace {
// From answer 12 in
// https://stackoverflow.com/questions/5279051/how-can-i-create-cartesian-product-of-vector-of-vectors
// by anumi
// modified to exclude duplicates
static std::vector<std::vector<int>> cartesianProduct(
    const std::vector<std::vector<int>> &v, bool allowDuplicates) {
  std::vector<std::vector<int>> s = {{}};
  for (const auto &u : v) {
    std::vector<std::vector<int>> r;
    for (const auto &x : s) {
      for (const auto y : u) {
        // check for duplicates
        if (allowDuplicates || std::find(x.begin(), x.end(), y) == x.end()) {
          r.push_back(x);
          r.back().push_back(y);
        }
      }
    }
    if (r.empty()) {
      // unable to extend
      return r;
    }
    s = std::move(r);
  }
  return s;
}
}  // namespace

// move this to constructor if the create new core path can be removed from
// RGroupDecomposition::add
void RCore::init() {
  findIndicesWithRLabel();
  countUserRGroups();
  buildMatchingMol();
}

void RCore::findIndicesWithRLabel() {
  // Find all the core atoms that have user
  // label and set their indices to 1 in core_atoms_with_user_labels
  core_atoms_with_user_labels.resize(core->getNumAtoms());
  for (const auto atom : core->atoms()) {
    int label;
    if (atom->getPropIfPresent(RLABEL, label) && label > 0) {
      core_atoms_with_user_labels.set(atom->getIdx());
    }
  }
}

RWMOL_SPTR RCore::extractCoreFromMolMatch(
    const ROMol &mol, const MatchVectType &match,
    const RGroupDecompositionParameters &params) const {
  auto extractedCore = boost::make_shared<RWMol>(mol);
  std::set<int> atomIndicesToKeep;
  std::vector<Bond *> newBonds;
  std::map<Atom *, int> dummyAtomMap;
  std::map<const Atom *, Atom *> molAtomMap;
  for (auto &pair : match) {
    const auto queryAtom = core->getAtomWithIdx(pair.first);
    auto const targetAtom = extractedCore->getAtomWithIdx(pair.second);
    if (int rLabel; queryAtom->getPropIfPresent(RLABEL, rLabel)) {
      targetAtom->setProp(RLABEL, rLabel);
    }
    if (int rLabelType; queryAtom->getPropIfPresent(RLABEL_TYPE, rLabelType)) {
      targetAtom->setProp(RLABEL_TYPE, rLabelType);
    }

    if (queryAtom->getAtomicNum() == 0 && queryAtom->hasProp(RLABEL) &&
        queryAtom->getDegree() == 1) {
      continue;
    } else {
      atomIndicesToKeep.insert(pair.second);
      molAtomMap[mol.getAtomWithIdx(pair.second)] = targetAtom;
      int neighborNumber = -1;
#ifdef VERBOSE
      std::cerr << "Atom Chirality In " << targetAtom->getChiralTag()
                << std::endl;
#endif
      bool isChiral = targetAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
                      targetAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW;
      if (isChiral && !params.substructmatchParams.useChirality) {
        // if we're not doing chiral matching don't copy chirality to the
        // extracted core.
        targetAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
        isChiral = false;
      }
      // collect neighbors in vector, so we can add atoms while looping
      std::vector<Atom *> targetNeighborAtoms;
      for (auto targetNeighborAtom : extractedCore->atomNeighbors(targetAtom)) {
        targetNeighborAtoms.push_back(targetNeighborAtom);
      }
      // explicit hydrogens to keep in extracted core to preserve chiralty
      std::vector<int> hydrogensToAdd;
      for (const auto targetNeighborAtom : targetNeighborAtoms) {
        ++neighborNumber;
        auto targetNeighborIndex = targetNeighborAtom->getIdx();
        auto queryNeighborMapping = std::find_if(
            match.begin(), match.end(),
            [this, targetNeighborIndex, pair](const auto &p) {
              return p.second == static_cast<int>(targetNeighborIndex) &&
                     core->getBondBetweenAtoms(pair.first, p.first);
            });
        if (queryNeighborMapping == match.end()) {
          if (targetNeighborAtom->getAtomicNum() == 1) {
            // Hydrogen needed to define chirality is present in target but
            // not mapped to core.  Copy it to the extracted core
            hydrogensToAdd.push_back(static_cast<int>(targetNeighborIndex));
            molAtomMap[mol.getAtomWithIdx(targetNeighborIndex)] =
                targetNeighborAtom;
          } else if (isChiral) {
            // There is a heavy sidechain in the decomp that is connected to
            // the core by an unknown bond (onlyMatchAtRGroups = False and
            // allowMultipleRGroupsOnUnlabelled = False).
            // As there is no explicit target bond chriality is not preserved.
            // In some cases we could handle chirality, but that has not been
            // implemented (if there is only one free bond on the query that
            // would be easy- or we could arbitrarily assign bonds).
            isChiral = false;
            targetAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
          }
          continue;
        }
        const auto queryNeighbor =
            core->getAtomWithIdx((*queryNeighborMapping).first);
        if (queryNeighbor->getAtomicNum() == 0 &&
            queryNeighbor->hasProp(RLABEL) && queryNeighbor->getDegree() == 1) {
          auto newDummy = new Atom(*queryNeighbor);
          dummyAtomMap[newDummy] = static_cast<int>(targetNeighborIndex);
          newDummy->clearComputedProps();
          const auto newDummyIdx =
              extractedCore->addAtom(newDummy, false, true);
          atomIndicesToKeep.insert(static_cast<int>(newDummyIdx));
          auto connectingBond =
              extractedCore
                  ->getBondBetweenAtoms(pair.second, targetNeighborIndex)
                  ->copy();
          if (connectingBond->getStereo() > Bond::BondStereo::STEREOANY) {
            // stereo double bonds
            connectingBond->setStereo(Bond::BondStereo::STEREOANY);
          }

          if (connectingBond->getBeginAtomIdx() == targetNeighborIndex) {
            connectingBond->setBeginAtomIdx(newDummyIdx);
          } else {
            connectingBond->setEndAtomIdx(newDummyIdx);
          }
          newBonds.push_back(connectingBond);

          // Check to see if we are breaking a stereo bond definition, by removing one of the stereo atoms
          // If so, set to the new atom
          for (auto bond: extractedCore->atomBonds(targetAtom)) {
            if (bond->getIdx() == connectingBond->getIdx()) {
              continue;
            }

            if (bond->getStereo() > Bond::STEREOANY) {
              auto &stereoAtoms = bond->getStereoAtoms();
              for (int& stereoAtom : stereoAtoms) {
                if (stereoAtom == static_cast<int>(targetNeighborIndex)) {
                  stereoAtom = static_cast<int>(newDummyIdx);
                }
              }
            }
          }

          // Chirality parity stuff see RDKit::replaceCore in
          // Code/GraphMol/ChemTransforms/ChemTransforms.cpp
          if (isChiral) {
            bool switchIt = false;
            switch (extractedCore->getAtomDegree(targetAtom)) {
              case 4:
                if (!(neighborNumber % 2)) {
                  switchIt = true;
                }
                break;
              case 3:
                if (neighborNumber == 1) {
                  switchIt = true;
                }
                break;
            }
            // because this neighbor will be moved at the end of the neighbor
            // list, we need to decrement the neighbor number in case additional
            // dummy atoms are to be connected to the targetAtom
            --neighborNumber;
            if (switchIt) {
              targetAtom->invertChirality();
            }
          }
        }
      }
      for (const int index : hydrogensToAdd) {
        atomIndicesToKeep.insert(index);
      }
#ifdef VERBOSE
      std::cerr << "Atom Chirality Out " << targetAtom->getChiralTag()
                << std::endl;
#endif
    }
  }
  for (const auto newBond : newBonds) {
    extractedCore->addBond(newBond, true);
  }

  // Now delete atom's that are not in the core.
  extractedCore->beginBatchEdit();
  boost::dynamic_bitset<> removedAtoms(mol.getNumAtoms());
  for (const auto atom : extractedCore->atoms()) {
    if (atomIndicesToKeep.find(atom->getIdx()) == atomIndicesToKeep.end()) {
      extractedCore->removeAtom(atom);
      removedAtoms.set(atom->getIdx());
    }
  }
  extractedCore->commitBatchEdit();

  // Copy molecule coordinates to extracted core
  details::updateSubMolConfs(mol, *extractedCore, removedAtoms);

  for (auto citer = mol.beginConformers(); citer != mol.endConformers();
       ++citer) {
    Conformer &newConf =
        extractedCore->getConformer(static_cast<int>((*citer)->getId()));
    for (const auto &[fst, snd] : dummyAtomMap) {
      newConf.setAtomPos(fst->getIdx(), (*citer)->getAtomPos(snd));
    }
  }

  // If the molecule has no coordinates and the core does, copy those over
  if (!mol.getNumConformers() && core->getNumConformers()) {
    ROMol molCopy(mol);
    for (auto citer = core->beginConformers(); citer != core->endConformers();
         ++citer) {
      const auto newConf = new Conformer(mol.getNumAtoms());
      newConf->setId((*citer)->getId());
      newConf->set3D((*citer)->is3D());
      for (const auto &[fst, snd] : match) {
        newConf->setAtomPos(snd, (*citer)->getAtomPos(fst));
      }
      molCopy.addConformer(newConf);
    }
    details::updateSubMolConfs(molCopy, *extractedCore, removedAtoms);
    molCopy.clearConformers();

    for (const auto atom : extractedCore->atoms()) {
      if (isUserRLabel(*atom)) {
        int rLabel = atom->getProp<int>(RLABEL);
        for (const auto coreAtom : core->atoms()) {
          if (int l; coreAtom->getPropIfPresent(RLABEL, l) && l == rLabel) {
            int i = 0;
            for (auto citer = core->beginConformers();
                 citer != core->endConformers(); ++citer, ++i) {
              extractedCore->getConformer(i).setAtomPos(
                  atom->getIdx(), (*citer)->getAtomPos(coreAtom->getIdx()));
            }
            break;
          }
        }
      }
    }
  }

  // Copy over any stereo groups that lie in the extracted core
  details::copyStereoGroups(molAtomMap, mol, *extractedCore);

  extractedCore->clearComputedProps(true);
  extractedCore->updatePropertyCache(false);

#ifdef VERBOSE
  std::cerr << "Extracted core smiles " << MolToSmiles(*extractedCore)
            << std::endl;
  std::cerr << "Extracted core smarts " << MolToSmarts(*extractedCore)
            << std::endl;
#endif

  try {
    unsigned int failed;
    MolOps::sanitizeMol(*extractedCore, failed,
                        MolOps::SANITIZE_SYMMRINGS | MolOps::SANITIZE_CLEANUP);
  } catch (const MolSanitizeException &) {
  }

  return extractedCore;
}

// Return a copy of core where dummy atoms are replaced by
// the respective matching atom in mol, while other atoms have
// their aromatic flag and formal charge copied from
// the respective matching atom in mol
ROMOL_SPTR RCore::replaceCoreAtomsWithMolMatches(
    const ROMol &mol, const MatchVectType &match) const {
  auto coreReplacedAtoms = boost::make_shared<RWMol>(*core);
  for (const auto &p : match) {
    auto atom = coreReplacedAtoms->getAtomWithIdx(p.first);
    if (isAtomWithMultipleNeighborsOrNotDummyRGroupAttachment(*atom)) {
      auto molAtom = mol.getAtomWithIdx(p.second);
      replaceCoreAtom(*coreReplacedAtoms, *atom, *molAtom);
    }
  }

  std::map<int, int> matchLookup(match.cbegin(), match.cend());
  for (auto bond : coreReplacedAtoms->bonds()) {
    if (bond->hasQuery()) {
      const auto molBond =
          mol.getBondBetweenAtoms(matchLookup[bond->getBeginAtomIdx()],
                                  matchLookup[bond->getEndAtomIdx()]);
      if (molBond == nullptr) {
        // this can happen if we have a user-defined R group that is not
        // matched in the query
        CHECK_INVARIANT(bond->getBeginAtom()->getAtomicNum() == 0 ||
                            bond->getEndAtom()->getAtomicNum() == 0,
                        "Failed to find core bond in molecule");
      } else {
        Bond newBond(molBond->getBondType());
        newBond.setIsAromatic(molBond->getIsAromatic());
        coreReplacedAtoms->replaceBond(bond->getIdx(), &newBond, true);
      }
    }
  }

#ifdef VERBOSE
  std::cerr << "Original core smarts  " << MolToSmarts(*core) << std::endl;
  std::cerr << "Dummy replaced core smarts  " << MolToSmarts(*coreReplacedAtoms)
            << std::endl;
#endif

  if (mol.getNumConformers() > 0) {
    // if the input structure has coordinates copy them to the core
    if (!coreReplacedAtoms->getNumConformers()) {
      coreReplacedAtoms->addConformer(
          new Conformer(coreReplacedAtoms->getNumAtoms()));
    }
    auto &replacedConformer = coreReplacedAtoms->getConformer();
    const auto &molConformer = mol.getConformer();

    for (const auto &p : match) {
      auto molPoint = molConformer.getAtomPos(p.second);
      replacedConformer.setAtomPos(p.first, molPoint);
    }
  } else {
    // otherwise, delete all core coordinates from the replaced core
    coreReplacedAtoms->clearConformers();
  }

  return coreReplacedAtoms;
}

void RCore::replaceCoreAtom(RWMol &mol, Atom &atom, const Atom &other) const {
  auto atomicNumber = other.getAtomicNum();
  auto targetAtom = &atom;
  bool wasDummy = (atom.getAtomicNum() == 0);
  if (wasDummy || atom.hasQuery()) {
    if (atom.hasQuery()) {
      Atom newAtom(atomicNumber);
      auto atomIdx = atom.getIdx();
      mol.replaceAtom(atomIdx, &newAtom, false, true);
      targetAtom = mol.getAtomWithIdx(atomIdx);
    } else {
      atom.setAtomicNum(atomicNumber);
    }
  }
  targetAtom->setIsAromatic(other.getIsAromatic());
  targetAtom->setFormalCharge(other.getFormalCharge());
  if (wasDummy) {
    targetAtom->setNoImplicit(true);
    unsigned int numHs = 0;
    const auto &otherMol = other.getOwningMol();
    for (const auto &nbri :
         boost::make_iterator_range(otherMol.getAtomNeighbors(&other))) {
      const auto nbrAtom = otherMol[nbri];
      if (nbrAtom->getAtomicNum() == 1) {
        ++numHs;
      }
    }
    targetAtom->setNumExplicitHs(numHs + other.getTotalNumHs());
    targetAtom->updatePropertyCache(false);
  }
}

// matching the core to the target molecule is a two step process
// First match to a reduced representation (the core minus terminal user
// R-groups). Next, match the R-groups. We do this as the core may not be a
// substructure match for the molecule if a single molecule atom matches 2
// terminal user defined RGroup attachments (see
// https://github.com/rdkit/rdkit/pull/4002) buildMatchingMol() creates the
// reduced representation from the core and matchTerminalUserRGroups() adds in
// the terminal R Groups to the match

// Builds a matching molecule which is the core with terminal user R groups
// removed Also creates the data structures used in matching the R groups
void RCore::buildMatchingMol() {
  matchingMol = boost::make_shared<RWMol>(*core);
  terminalRGroupDummyAtoms.clear();
  terminalRGroupAtomToNeighbor.clear();
  RWMol::ATOM_PTR_VECT atomsToRemove;
  for (auto atom : matchingMol->atoms()) {
    // keep track of the original core index in the matching molecule atom
    atom->setProp<int>(RLABEL_CORE_INDEX, atom->getIdx());
    // TODO for unlabelled core attachments if the heavy neighbor is not dummy
    // then keep one attachment
    if (atom->getAtomicNum() == 0 && atom->getDegree() == 1 &&
        isDummyRGroupAttachment(*atom)) {
      // remove terminal user R groups and save core atom index and mapping to
      // heavy neighbor
      atomsToRemove.push_back(atom);
      terminalRGroupDummyAtoms.insert(atom->getIdx());
      const int neighborIdx = *matchingMol->getAtomNeighbors(atom).first;
      terminalRGroupAtomToNeighbor.emplace(atom->getIdx(), neighborIdx);
    }
  }

  for (auto atom : atomsToRemove) {
    matchingMol->removeAtom(atom);
  }
}

// Given a matching molecule substructure match to a target molecule, return
// core matches with terminal user R groups matched
std::vector<MatchVectType> RCore::matchTerminalUserRGroups(
    const RWMol &target, MatchVectType match,
    const SubstructMatchParameters &sssParams) const {
  // Transform match indexed by matching molecule atoms to a map
  // indexed by core atoms
  std::transform(match.begin(), match.end(), match.begin(),
                 [this](const std::pair<int, int> &mapping) {
                   auto queryIdx =
                       this->matchingIndexToCoreIndex(mapping.first);
                   std::pair<int, int> newMapping(queryIdx, mapping.second);
                   return newMapping;
                 });
  std::map<int, int> matchMap(match.cbegin(), match.cend());

  std::vector<MatchVectType> allMappings;
  if (terminalRGroupDummyAtoms.empty()) {
    allMappings.push_back(match);
    return allMappings;
  }

  // build a set of target atoms currently mapped
  std::set<int> mappedTargetIdx;
  std::transform(
      match.cbegin(), match.cend(),
      std::inserter(mappedTargetIdx, mappedTargetIdx.begin()),
      [](const std::pair<int, int> &mapping) { return mapping.second; });

  // Dummy atoms/r group attachments that cannot be mapped to target atoms
  std::vector<int> missingDummies;
  // A map of terminal dummies/R group attachment points to a list of possible
  // target atoms that the R group can map to
  std::map<int, std::vector<int>> availableMappingsForDummyMap;

  std::set<int> symmetricHydrogens;

  // keep a count of target atoms and the number of connections to the core
  std::map<int, int> targetAtomBondsToCoreCounts;
  for (const auto coreAtom : core->atoms()) {
    if (coreAtom->getAtomicNum() == 1) {
      continue;
    }
    if (terminalRGroupDummyAtoms.find(coreAtom->getIdx()) !=
        terminalRGroupDummyAtoms.end()) {
      continue;
    }
    std::vector<int> dummyIndexes;
    for (auto neighbor : core->atomNeighbors(coreAtom)) {
      if (terminalRGroupDummyAtoms.find(neighbor->getIdx()) !=
          terminalRGroupDummyAtoms.end()) {
        dummyIndexes.push_back(neighbor->getIdx());
      }
    }
    if (dummyIndexes.empty()) {
      continue;
    }

    // Sort dummies based on user RLABEL (ascending) then unlabeled
    // (descending as they are negative)
    std::sort(dummyIndexes.begin(), dummyIndexes.end(), [this](int a, int b) {
      auto dummy = core->getAtomWithIdx(a);
      auto otherDummy = core->getAtomWithIdx(b);
      auto l1 = dummy->getProp<int>(RLABEL);
      auto l2 = otherDummy->getProp<int>(RLABEL);
      return l1 > 0 && l2 > 0 ? l1 < l2 : l1 > l2;
    });

    // Find what target atoms will bond to these dummies
    std::vector<std::vector<int>> neighborDummyLists;
    for (auto dummyIndex : dummyIndexes) {
      const int neighborIdx =
          terminalRGroupAtomToNeighbor.find(dummyIndex)->second;
      const auto coreBond = core->getBondBetweenAtoms(dummyIndex, neighborIdx);
      // find the atom in the target mapped to the neighbor in the core
      const int targetIdx = matchMap[neighborIdx];
      const auto targetAtom = target.getAtomWithIdx(targetIdx);
      ROMol::ADJ_ITER nbrIter, endNbrs;
      std::vector<int> available;
      // now look for neighbors of that target atom that are not mapped to a
      // core atom- the dummy atom can potentially be mapped to each of those
      boost::tie(nbrIter, endNbrs) = target.getAtomNeighbors(targetAtom);
      while (nbrIter != endNbrs) {
        if (mappedTargetIdx.find(*nbrIter) == mappedTargetIdx.end()) {
          const auto targetBond =
              target.getBondBetweenAtoms(targetIdx, *nbrIter);
          // check for bond compatibility
          if (bondCompat(coreBond, targetBond, sssParams)) {
            available.push_back(*nbrIter);
            targetAtomBondsToCoreCounts[*nbrIter]++;
          }
        }
        ++nbrIter;
      }

      if (available.size() > 1) {
        bool allHydrogens = std::all_of(
            available.begin(), available.end(), [&target](const int idx) {
              return target.getAtomWithIdx(idx)->getAtomicNum() == 1;
            });
        if (allHydrogens) {
          // If all neighbors are hydrogens we don't need to iterate through
          // them- just assign the first free hydrogen. Could extend to cover
          // symmetric groups in general
          auto hydrogen = std::find_if(available.begin(), available.end(),
                                       [&symmetricHydrogens](const int idx) {
                                         return symmetricHydrogens.find(idx) ==
                                                symmetricHydrogens.end();
                                       });
          int singleHydrogen =
              hydrogen == available.end() ? *available.begin() : *hydrogen;
          std::vector<int> hydrogenVec{singleHydrogen};
          available = hydrogenVec;
          symmetricHydrogens.insert(singleHydrogen);
        }
      }
      neighborDummyLists.push_back(available);
    }

    // Now search through all the dummies and see if we need to exclude any as
    // it may not be possible to assign target atoms to all R groups.  If that
    // is the case exclude R Groups with higher labels

    // use negative numbers to indicate a group is excluded
    std::vector<std::vector<int>> neighborListsWithUnmapped(
        neighborDummyLists.size());
    int start = 0;
    std::transform(neighborDummyLists.begin(), neighborDummyLists.end(),
                   neighborListsWithUnmapped.begin(),
                   [&start](std::vector<int> v) {
                     v.push_back(--start);
                     return v;
                   });
    // could optimize this product as depth search to return as soon as we
    // have a permutation of all positive numbers
    auto cp = cartesianProduct(neighborListsWithUnmapped, false);
    // now sort
    bool foundAll = false;
    std::sort(cp.begin(), cp.end(),
              [&foundAll](const std::vector<int> a, std::vector<int> b) {
                auto isNegative = [](int v) -> bool { return v < 0; };
                // firstly to minimize number of unmapped dummies
                const int numUnmappedA =
                    std::count_if(a.begin(), a.end(), isNegative);
                const int numUnmappedB =
                    std::count_if(b.begin(), b.end(), isNegative);
                if (numUnmappedA < numUnmappedB) {
                  return true;
                }
                if (numUnmappedB == 0 || numUnmappedA == 0) {
                  foundAll = true;
                }
                if (!foundAll && numUnmappedA == numUnmappedB && numUnmappedA) {
                  // in a tie prefer the permutation that excludes the r group
                  // with the lowest label
                  // don't need to sort these if we've found an all mapped
                  // combination
                  for (size_t i = 0; i < a.size(); ++i) {
                    const int v1 = a[i];
                    const int v2 = b[i];
                    if ((v1 > 0 && v2 > 0) || (v1 < 0 && v2 < 0)) {
                      continue;
                    }
                    if (v1 > 0 && v2 < 0) {
                      return true;
                    }
                    if (v1 < 0 && v2 > 0) {
                      break;
                    }
                  }
                }
                return false;
              });
    auto best = cp[0];
    for (size_t i = 0; i < dummyIndexes.size(); i++) {
      auto dummyIdx = dummyIndexes[i];
      if (best[i] < 0) {
        // We can't map this dummy to a target atom.  That is OK if it is an
        // unlabeled core attachment atom or a user R label connected to a
        // query or wildcard atom
        const auto dummy = core->getAtomWithIdx(dummyIdx);
        if (dummy->hasProp(UNLABELED_CORE_ATTACHMENT)) {
          missingDummies.push_back(dummyIdx);
        } else if (isUserRLabel(*dummy) &&
                   (coreAtom->getAtomicNum() == 0 || coreAtom->hasQuery())) {
          // https://github.com/rdkit/rdkit/issues/4505
          missingDummies.push_back(dummyIdx);
        } else {
          return allMappings;
        }
      } else {
        availableMappingsForDummyMap[dummyIdx] = neighborDummyLists[i];
      }
    }
  }
  if (availableMappingsForDummyMap.empty()) {
    allMappings.push_back(match);
    return allMappings;
  }

  std::vector<int> dummiesWithMapping;
  std::vector<std::vector<int>> availableMappingsForDummy;
  for (const auto &mapping : availableMappingsForDummyMap) {
    dummiesWithMapping.push_back(mapping.first);
    availableMappingsForDummy.push_back(mapping.second);
  }

  // enumerate over all available atoms using a cartesian product.
  // only allow duplicates if a target atom can be bonded to more than one query
  // atom
  const bool allowDuplicates =
      std::find_if(targetAtomBondsToCoreCounts.begin(),
                   targetAtomBondsToCoreCounts.end(),
                   [](const std::pair<int, int> &p) { return p.second > 1; }) !=
      targetAtomBondsToCoreCounts.end();
  const auto allAvailableMappings =
      cartesianProduct(availableMappingsForDummy, allowDuplicates);
  if (allAvailableMappings.empty()) {
    allMappings.push_back(match);
    return allMappings;
  }
  // the size of the final mapping
  size_t size = allAvailableMappings[0].size() + match.size();
  // these indices are needed for the whole molecule match check functor

  std::unique_ptr<RWMol> checkCore;
  std::map<size_t, size_t> coreToCheck;
  const std::string indexProp("__core_index__");
  bool hasMissing = !missingDummies.empty();
  if (hasMissing) {
    // if there are dummies that we can't map these need to be removed from the
    // query before atom-by-atom matching.  Create a copy of the query for that
    // and use properties to map atoms back to the core
    for (auto atom : core->atoms()) {
      atom->setProp(indexProp, atom->getIdx());
    }
    checkCore = std::make_unique<RWMol>(*core);
    std::sort(missingDummies.begin(), missingDummies.end(),
              std::greater<int>());
    for (int index : missingDummies) {
      auto [nbrIdx, endNbrs] =
          checkCore->getAtomNeighbors(checkCore->getAtomWithIdx(index));
      auto neighborAtom = checkCore->getAtomWithIdx(*nbrIdx);
      checkCore->removeAtom(index);
      neighborAtom->updatePropertyCache(false);
    }
    size_t index = 0U;
    for (const auto atom : checkCore->atoms()) {
      auto coreIndex = atom->getProp<int>(indexProp);
      coreToCheck[coreIndex] = index++;
    }
    for (auto atom : core->atoms()) {
      atom->clearProp(indexProp);
    }
  }

  auto queryIndices = new std::uint32_t[size];
  auto targetIndices = new std::uint32_t[size];
  for (size_t position = 0; position < match.size(); position++) {
    const auto &pair = match[position];
    auto queryIndex = hasMissing ? coreToCheck[pair.first] : pair.first;
    queryIndices[position] = queryIndex;
    targetIndices[position] = pair.second;
  }

  auto queryMatchingMol = hasMissing ? checkCore.get() : core.get();
  MolMatchFinalCheckFunctor molMatchFunctor(*queryMatchingMol, target,
                                            sssParams);
  boost::dynamic_bitset<> targetBondsPresent(target.getNumBonds());

  // Filter all available mappings removing those that violate chirality or have
  // duplicate bonds
  for (const auto &dummyMapping : allAvailableMappings) {
    CHECK_INVARIANT(match.size() + dummyMapping.size() == size,
                    "Size error in dummy mapping");
    auto duplicateBonds = false;
    targetBondsPresent.reset();
    for (size_t i = 0; i < dummyMapping.size(); i++) {
      size_t position = match.size() + i;
      auto queryIndex = hasMissing ? coreToCheck[dummiesWithMapping[i]]
                                   : dummiesWithMapping[i];
      queryIndices[position] = queryIndex;
      targetIndices[position] = dummyMapping[i];
      if (allowDuplicates) {
        const int neighborIdx =
            terminalRGroupAtomToNeighbor.find(dummiesWithMapping[i])->second;
        const int targetNeighborIdx = matchMap[neighborIdx];
        const auto targetBond =
            target.getBondBetweenAtoms(dummyMapping[i], targetNeighborIdx);
        CHECK_INVARIANT(targetBond != nullptr,
                        "Matching target bond not found");
        const auto targetBondIdx = targetBond->getIdx();
        // check for duplicates
        if (targetBondsPresent[targetBondIdx]) {
          duplicateBonds = true;
          break;
        }
        targetBondsPresent[targetBondIdx] = 1;
      }
    }
    // use MolMatchFinalCheckFunctor to check this match works with chirality
    if (!duplicateBonds && molMatchFunctor(queryIndices, targetIndices)) {
      MatchVectType matchWithDummy(match);
      for (size_t i = 0; i < dummyMapping.size(); i++) {
        matchWithDummy.emplace_back(dummiesWithMapping[i], dummyMapping[i]);
      }
      allMappings.push_back(matchWithDummy);
    }
  }

  delete[] queryIndices;
  delete[] targetIndices;
  return allMappings;
}

// This function checks the bond environment is valid when onlyMatchAtRGroups
// is set.  When this function is called attachmentIdx is the index of an atom
// in the target that is mapped to an R group.  Validates that all core bonds to
// the attachment point are present - in certain circumstances there may be two
// core bonds to a target attachment point and when onlyMatchAtRGroups is set
// both bonds should be present.
bool RCore::checkAllBondsToAttachmentPointPresent(
    const ROMol &mol, const int attachmentIdx,
    const MatchVectType &mapping) const {
  const auto atom = mol.getAtomWithIdx(attachmentIdx);
  std::set<int> coreNeighbors;
  for (const auto &nbri :
       boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
    const auto &nbr = mol[nbri];
    // could a neighbor to an r group attachment match another r group
    // attachment?  I don't think so.
    if (nbr->getAtomicNum() >= 1) {
      const auto match = std::find_if(
          mapping.cbegin(), mapping.cend(), [nbri](std::pair<int, int> p) {
            return p.second == static_cast<int>(nbri);
          });
      if (match != mapping.end()) {
        auto coreAtom = core->getAtomWithIdx(match->first);
        // don't need to match a non terminal user R group
        // if (!(coreAtom->getDegree() > 1 && isUserRLabel(*coreAtom))) {
        if (!(coreAtom->getAtomicNum() == 0 && isUserRLabel(*coreAtom))) {
          coreNeighbors.insert(match->first);
        }
      }
    }
  }

  CHECK_INVARIANT(
      coreNeighbors.size() >= 1,
      "Unable to find target atom(s) matching core for attachment point");
  if (coreNeighbors.size() == 1) {
    // currently this routine is only called when we know the attachment to
    // one core atom exists.
    return true;
  }

  // at this point we know the target atom is connected to two or more core
  // atoms.  Now check we have attachment points for each.
  //  There should be a map entry for a different attachment point for each of
  //  the core atoms.
  std::vector<int> coreDummyIndices;
  std::for_each(mapping.cbegin(), mapping.cend(),
                [&coreDummyIndices, attachmentIdx](std::pair<int, int> p) {
                  if (p.second == attachmentIdx) {
                    coreDummyIndices.push_back(p.first);
                  }
                });
  if (coreDummyIndices.size() != coreNeighbors.size()) {
    return false;
  }
  // check that there is a bond to an attachment point for every neighbor
  for (const auto neighborIdx : coreNeighbors) {
    bool foundBond = false;
    for (const auto dummyIdx : coreDummyIndices) {
      if (core->getBondBetweenAtoms(neighborIdx, dummyIdx) != nullptr) {
        foundBond = true;
        break;
      }
    }
    if (!foundBond) {
      return false;
    }
  }
  return true;
}

// Convert a matching molecule index to a core index
int RCore::matchingIndexToCoreIndex(int matchingIndex) const {
  auto atom = matchingMol->getAtomWithIdx(matchingIndex);
  CHECK_INVARIANT(atom->hasProp(RLABEL_CORE_INDEX),
                  "Matched atom missing core index");
  return atom->getProp<int>(RLABEL_CORE_INDEX);
}

// Create tautomer query for the matching mol on demand and cache for performance
// If the tautomer query cannot be created (because we can't kekulize the query)
// then nullptr will be returned and we revert to non-tautomer match
std::shared_ptr<TautomerQuery> RCore::getMatchingTautomerQuery() {
  if (!checkedForTautomerQuery) {
    try {
      // Enumerate tautomers from a sanitized copy of the matching molecule
      RWMol copy(*matchingMol);
      // If the core has had rgroup labels removed when creating the matching mol
      // then we need to update properties.  Should a full sanitization be done?
      // MolOps::sanitizeMol(*copy);
      copy.updatePropertyCache(false);
      std::shared_ptr<TautomerQuery> tautomerQuery(
          TautomerQuery::fromMol(copy));
      matchingTautomerQuery = tautomerQuery;
    } catch (const MolSanitizeException &) {
      matchingTautomerQuery = nullptr;
    }
    checkedForTautomerQuery = true;
  }
  return matchingTautomerQuery;
}

}  // namespace RDKit
