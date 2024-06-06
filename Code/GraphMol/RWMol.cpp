//
//  Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/tokenizer.hpp>

// our stuff
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "RWMol.h"
#include "Atom.h"
#include "Bond.h"
#include "BondIterators.h"
#include "RingInfo.h"
#include "SubstanceGroup.h"

namespace RDKit {

namespace {
void insertStereoGroups(RWMol &mol, const ROMol &other,
                        unsigned int origNumAtoms, unsigned int origNumBonds) {
  std::vector<RDKit::Atom *> abs_atoms;
  std::vector<RDKit::Bond *> abs_bonds;
  std::vector<RDKit::StereoGroup> new_groups;
  new_groups.reserve(mol.getStereoGroups().size());
  for (const auto &sg : mol.getStereoGroups()) {
    // The sdf specification forbids more than one ABS stereo group, but we
    // don't enforce that in our code. But if we see more than one ABS groups
    // here, just merge the atoms and bonds in them into one group. Other stereo
    // groups are just forwarded.
    if (sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE) {
      auto &atoms = sg.getAtoms();
      auto &bonds = sg.getBonds();
      abs_atoms.insert(abs_atoms.end(), atoms.begin(), atoms.end());
      abs_bonds.insert(abs_bonds.end(), bonds.begin(), bonds.end());
    } else {
      new_groups.emplace_back(sg);
    }
  }

  for (const auto &sg : other.getStereoGroups()) {
    // update the stereo group's atom and bond indices
    std::vector<RDKit::Atom *> new_atoms;
    std::vector<RDKit::Bond *> new_bonds;
    for (auto atom : sg.getAtoms()) {
      auto idx = atom->getIdx() + origNumAtoms;
      new_atoms.push_back(mol.getAtomWithIdx(idx));
    }
    for (auto bond : sg.getBonds()) {
      auto idx = bond->getIdx() + origNumBonds;
      new_bonds.push_back(mol.getBondWithIdx(idx));
    }

    // Collect all ABS atoms and bonds so they are added as a single group
    if (sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE) {
      abs_atoms.insert(abs_atoms.end(), new_atoms.begin(), new_atoms.end());
      abs_bonds.insert(abs_bonds.end(), new_bonds.begin(), new_bonds.end());
    } else {
      RDKit::StereoGroup new_group(sg.getGroupType(), new_atoms, new_bonds,
                                   sg.getReadId());
      // default write ID to 0 to avoid id clashes. We can use
      // assignStereoGroupIds() later on to assign new IDs
      new_group.setWriteId(0);
      new_groups.push_back(new_group);
    }
  }
  new_groups.emplace_back(RDKit::StereoGroupType::STEREO_ABSOLUTE, abs_atoms,
                          abs_bonds);
  mol.setStereoGroups(new_groups);
}

void insertSubstanceGroups(RWMol &mol, const RWMol &other,
                           unsigned int origNumAtoms,
                           unsigned int origNumBonds) {
  for (auto sgroup : getSubstanceGroups(other)) {
    sgroup.setOwningMol(&mol);

    // update the sgroup's atom and bond indices
    auto atom_indices = sgroup.getAtoms();
    std::transform(atom_indices.begin(), atom_indices.end(),
                   atom_indices.begin(),
                   [&origNumAtoms](unsigned int old_index) {
                     return origNumAtoms + old_index;
                   });
    sgroup.setAtoms(atom_indices);

    auto bond_indices = sgroup.getBonds();
    std::transform(bond_indices.begin(), bond_indices.end(),
                   bond_indices.begin(),
                   [&origNumBonds](unsigned int old_index) {
                     return origNumBonds + old_index;
                   });
    sgroup.setBonds(bond_indices);

    // patoms
    auto patom_indices = sgroup.getParentAtoms();
    std::transform(patom_indices.begin(), patom_indices.end(),
                   patom_indices.begin(),
                   [&origNumAtoms](unsigned int old_index) {
                     return origNumAtoms + old_index;
                   });
    sgroup.setParentAtoms(patom_indices);

    // cstates (these are references, can be updated in place)
    for (auto &cstate : sgroup.getCStates()) {
      cstate.bondIdx = origNumBonds + cstate.bondIdx;
    }

    // attachment points (can also be updated in place)
    for (auto &sap : sgroup.getAttachPoints()) {
      sap.aIdx = origNumAtoms + sap.aIdx;
      if (sap.lvIdx != -1) {
        sap.lvIdx = static_cast<int>(origNumAtoms + sap.lvIdx);
      }
    }

    addSubstanceGroup(mol, sgroup);
  }
}
}  // namespace

RWMol &RWMol::operator=(const RWMol &other) {
  if (this != &other) {
    this->clear();
    numBonds = 0;
    initFromOther(other, false, -1);
  }
  return *this;
}

void RWMol::insertMol(const ROMol &other) {
  auto origNumAtoms = getNumAtoms();
  auto origNumBonds = getNumBonds();
  for (const auto oatom : other.atoms()) {
    Atom *newAt = oatom->copy();
    const bool updateLabel = false;
    const bool takeOwnership = true;
    addAtom(newAt, updateLabel, takeOwnership);
    // take care of atom-numbering-dependent properties:
    INT_VECT nAtoms;
    if (newAt->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
      for (auto &val : nAtoms) {
        if (val < 0) {
          val = -1 * (-val + origNumAtoms);
        } else {
          val += origNumAtoms;
        }
      }
      newAt->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    }
  }

  for (const auto obond : other.bonds()) {
    Bond *bond_p = obond->copy();
    unsigned int idx1, idx2;
    idx1 = bond_p->getBeginAtomIdx() + origNumAtoms;
    idx2 = bond_p->getEndAtomIdx() + origNumAtoms;
    bond_p->setOwningMol(this);
    bond_p->setBeginAtomIdx(idx1);
    bond_p->setEndAtomIdx(idx2);
    for (auto &v : bond_p->getStereoAtoms()) {
      v += origNumAtoms;
    }
    const bool takeOwnership = true;
    addBond(bond_p, takeOwnership);
  }

  // add atom to any conformers as well, if we have any
  if (other.getNumConformers() && !getNumConformers()) {
    for (const auto &oconf : other.d_confs) {
      auto *nconf = new Conformer(getNumAtoms());
      nconf->set3D(oconf->is3D());
      nconf->setId(oconf->getId());
      for (unsigned int i = 0; i < oconf->getNumAtoms(); ++i) {
        nconf->setAtomPos(i + origNumAtoms, oconf->getAtomPos(i));
      }
      const bool assignId = false;
      addConformer(nconf, assignId);
    }
  } else if (getNumConformers()) {
    if (other.getNumConformers() == getNumConformers()) {
      ConformerIterator cfi;
      ConstConformerIterator ocfi;
      for (cfi = beginConformers(), ocfi = other.beginConformers();
           cfi != endConformers(); ++cfi, ++ocfi) {
        for (unsigned int i = 0; i < (*ocfi)->getNumAtoms(); ++i) {
          (*cfi)->setAtomPos(i + origNumAtoms, (*ocfi)->getAtomPos(i));
        }
      }
    }
  }

  // add stereo groups
  insertStereoGroups(*this, other, origNumAtoms, origNumBonds);
  // add substance groups
  insertSubstanceGroups(*this, other, origNumAtoms, origNumBonds);
}

unsigned int RWMol::addAtom(bool updateLabel) {
  auto *atom_p = new Atom();
  atom_p->setOwningMol(this);
  auto which = boost::add_vertex(d_graph);
  d_graph[which] = atom_p;
  atom_p->setIdx(which);
  if (updateLabel) {
    clearAtomBookmark(ci_RIGHTMOST_ATOM);
    setAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
  }

  // add atom to any conformers as well, if we have any
  for (auto &conf : d_confs) {
    conf->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
  }
  return rdcast<unsigned int>(which);
}

void RWMol::replaceAtom(unsigned int idx, Atom *atom_pin, bool,
                        bool preserveProps) {
  PRECONDITION(atom_pin, "bad atom passed to replaceAtom");
  URANGE_CHECK(idx, getNumAtoms());
  auto atom_p = atom_pin->copy();
  atom_p->setOwningMol(this);
  atom_p->setIdx(idx);
  auto vd = boost::vertex(idx, d_graph);
  if (preserveProps) {
    const bool replaceExistingData = false;
    atom_p->updateProps(*d_graph[vd], replaceExistingData);
  }

  const auto orig_p = d_graph[vd];
  delete orig_p;
  d_graph[vd] = atom_p;

  // handle bookmarks
  for (auto &ab : d_atomBookmarks) {
    for (auto &elem : ab.second) {
      if (elem == orig_p) {
        elem = atom_p;
      }
    }
  }

  // handle stereo group
  for (auto &group : d_stereo_groups) {
    auto groupId = group.getReadId();
    auto atoms = group.getAtoms();
    auto bonds = group.getBonds();
    auto aiter = std::find(atoms.begin(), atoms.end(), orig_p);
    while (aiter != atoms.end()) {
      *aiter = atom_p;
      ++aiter;
      aiter = std::find(aiter, atoms.end(), orig_p);
    }
    group = StereoGroup(group.getGroupType(), std::move(atoms),
                        std::move(bonds), groupId);
  }
};

void RWMol::replaceBond(unsigned int idx, Bond *bond_pin, bool preserveProps,
                        bool keepSGroups) {
  PRECONDITION(bond_pin, "bad bond passed to replaceBond");
  URANGE_CHECK(idx, getNumBonds());
  auto bIter = getEdges();
  for (unsigned int i = 0; i < idx; i++) {
    ++bIter.first;
  }
  const auto *obond = d_graph[*(bIter.first)];
  auto *bond_p = bond_pin->copy();
  bond_p->setOwningMol(this);
  bond_p->setIdx(idx);
  bond_p->setBeginAtomIdx(obond->getBeginAtomIdx());
  bond_p->setEndAtomIdx(obond->getEndAtomIdx());

  // Update explicit Hs, if set, on both ends. This was github #7128
  auto orderDifference =
      bond_p->getBondTypeAsDouble() - obond->getBondTypeAsDouble();
  if (orderDifference > 0) {
    for (auto atom : {bond_p->getBeginAtom(), bond_p->getEndAtom()}) {
      if (auto explicit_hs = atom->getNumExplicitHs(); explicit_hs > 0) {
        auto new_hs = static_cast<int>(explicit_hs - orderDifference);
        atom->setNumExplicitHs(std::max(new_hs, 0));
      }
    }
  }

  if (preserveProps) {
    const bool replaceExistingData = false;
    bond_p->updateProps(*d_graph[*(bIter.first)], replaceExistingData);
  }

  const auto orig_p = d_graph[*(bIter.first)];
  delete orig_p;
  d_graph[*(bIter.first)] = bond_p;

  if (!keepSGroups) {
    removeSubstanceGroupsReferencingBond(*this, idx);
  }

  // handle bookmarks
  for (auto &ab : d_bondBookmarks) {
    for (auto &elem : ab.second) {
      if (elem == orig_p) {
        elem = bond_p;
      }
    }
  }
};

Atom *RWMol::getActiveAtom() {
  if (hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
    return getAtomWithBookmark(ci_RIGHTMOST_ATOM);
  } else {
    return getLastAtom();
  }
};

void RWMol::setActiveAtom(Atom *at) {
  PRECONDITION(at, "NULL atom provided");
  clearAtomBookmark(ci_RIGHTMOST_ATOM);
  setAtomBookmark(at, ci_RIGHTMOST_ATOM);
};
void RWMol::setActiveAtom(unsigned int idx) {
  setActiveAtom(getAtomWithIdx(idx));
};

void RWMol::removeAtom(unsigned int idx) { removeAtom(getAtomWithIdx(idx)); }

void RWMol::removeAtom(Atom *atom, bool clearProps) {
  PRECONDITION(atom, "NULL atom provided");
  PRECONDITION(static_cast<RWMol *>(&atom->getOwningMol()) == this,
               "atom not owned by this molecule");
  unsigned int idx = atom->getIdx();

  // remove bonds attached to the atom
  //  In batch mode this will schedule bond removal
  std::vector<std::pair<unsigned int, unsigned int>> nbrs;
  ADJ_ITER b1, b2;
  boost::tie(b1, b2) = getAtomNeighbors(atom);
  while (b1 != b2) {
    nbrs.emplace_back(atom->getIdx(), rdcast<unsigned int>(*b1));
    ++b1;
  }
  for (auto &nbr : nbrs) {
    removeBond(nbr.first, nbr.second);
  }

  if (dp_delAtoms) {
    // we're in a batch edit
    // if atoms have been added since we started, resize dp_delAtoms
    if (dp_delAtoms->size() < getNumAtoms()) {
      dp_delAtoms->resize(getNumAtoms());
    }
    dp_delAtoms->set(idx);
    return;
  }

  // remove any bookmarks which point to this atom:
  ATOM_BOOKMARK_MAP *marks = getAtomBookmarks();
  auto markI = marks->begin();
  while (markI != marks->end()) {
    const ATOM_PTR_LIST &atoms = markI->second;
    // we need to copy the iterator then increment it, because the
    // deletion we're going to do in clearAtomBookmark will invalidate
    // it.
    auto tmpI = markI;
    ++markI;
    if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
      clearAtomBookmark(tmpI->first, atom);
    }
  }

  // loop over all atoms with higher indices and update their indices
  for (unsigned int i = idx + 1; i < getNumAtoms(); i++) {
    Atom *higher_index_atom = getAtomWithIdx(i);
    higher_index_atom->setIdx(i - 1);
  }

  // do the same with the coordinates in the conformations
  for (auto conf : d_confs) {
    RDGeom::POINT3D_VECT &positions = conf->getPositions();
    auto pi = positions.begin();
    for (unsigned int i = 0; i < getNumAtoms() - 1; i++) {
      ++pi;
      if (i >= idx) {
        positions[i] = positions[i + 1];
      }
    }
    positions.erase(pi);
  }
  // now deal with bonds:
  //   their end indices may need to be decremented and their
  //   indices will need to be handled and if they have an
  //   ENDPTS prop that includes idx, it will need updating.
  unsigned int nBonds = 0;
  EDGE_ITER beg, end;
  boost::tie(beg, end) = getEdges();
  std::string sprop;
  while (beg != end) {
    Bond *bond = d_graph[*beg++];
    if (bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
                               sprop)) {
      // This would ideally use ParseV3000Array but I'm buggered if I can get
      // the linker to find it.
      //      std::vector<unsigned int> oats =
      //          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(sprop);
      if ('(' == sprop.front() && ')' == sprop.back()) {
        sprop = sprop.substr(1, sprop.length() - 2);

        // This is doing what ParseV3000Array would do.
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
        unsigned int num_ats = std::stod(*tokens.begin());
        std::vector<unsigned int> oats;
        auto beg = tokens.begin();
        ++beg;
        std::transform(beg, tokens.end(), std::back_inserter(oats),
                       [](const std::string &a) { return std::stod(a); });

        auto idx_pos = std::find(oats.begin(), oats.end(), idx + 1);
        if (idx_pos != oats.end()) {
          oats.erase(idx_pos);
          --num_ats;
        }
        if (!num_ats) {
          bond->clearProp(RDKit::common_properties::_MolFileBondEndPts);
          bond->clearProp(common_properties::_MolFileBondAttach);
        } else {
          sprop = "(" + std::to_string(num_ats) + " ";
          for (auto &i : oats) {
            if (i > idx + 1) {
              --i;
            }
            sprop += std::to_string(i) + " ";
          }
          sprop[sprop.length() - 1] = ')';
          bond->setProp(RDKit::common_properties::_MolFileBondEndPts, sprop);
        }
      }
    }
    unsigned int tmpIdx = bond->getBeginAtomIdx();
    if (tmpIdx > idx) {
      bond->setBeginAtomIdx(tmpIdx - 1);
    }
    tmpIdx = bond->getEndAtomIdx();
    if (tmpIdx > idx) {
      bond->setEndAtomIdx(tmpIdx - 1);
    }
    bond->setIdx(nBonds++);
    for (auto bsi = bond->getStereoAtoms().begin();
         bsi != bond->getStereoAtoms().end(); ++bsi) {
      if ((*bsi) == rdcast<int>(idx)) {
        bond->getStereoAtoms().clear();
        break;
      } else if ((*bsi) > rdcast<int>(idx)) {
        --(*bsi);
      }
    }
  }

  removeSubstanceGroupsReferencingAtom(*this, idx);

  // Remove this atom from any stereo group
  removeAtomFromGroups(atom, d_stereo_groups);

  // clear computed properties and reset our ring info structure
  // they are pretty likely to be wrong now:
  if (clearProps) {
    clearComputedProps(true);
  }

  atom->setOwningMol(nullptr);

  // remove all connections to the atom:
  MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
  boost::clear_vertex(vd, d_graph);
  // finally remove the vertex itself
  boost::remove_vertex(vd, d_graph);
  delete atom;
}

void RWMol::removeAtom(Atom *atom) { removeAtom(atom, true); }

unsigned int RWMol::addBond(unsigned int atomIdx1, unsigned int atomIdx2,
                            Bond::BondType bondType) {
  URANGE_CHECK(atomIdx1, getNumAtoms());
  URANGE_CHECK(atomIdx2, getNumAtoms());
  PRECONDITION(atomIdx1 != atomIdx2, "attempt to add self-bond");
  PRECONDITION(!(boost::edge(atomIdx1, atomIdx2, d_graph).second),
               "bond already exists");

  auto *b = new Bond(bondType);
  b->setOwningMol(this);
  if (bondType == Bond::AROMATIC) {
    b->setIsAromatic(1);
    //
    // assume that aromatic bonds connect aromatic atoms
    //   This is relevant for file formats like MOL, where there
    //   is no such thing as an aromatic atom, but bonds can be
    //   marked aromatic.
    //
    getAtomWithIdx(atomIdx1)->setIsAromatic(1);
    getAtomWithIdx(atomIdx2)->setIsAromatic(1);
  }
  auto [which, ok] = boost::add_edge(atomIdx1, atomIdx2, d_graph);
  d_graph[which] = b;
  // unsigned int res = rdcast<unsigned int>(boost::num_edges(d_graph));
  ++numBonds;
  b->setIdx(numBonds - 1);
  b->setBeginAtomIdx(atomIdx1);
  b->setEndAtomIdx(atomIdx2);

  // if both atoms have a degree>1, reset our ring info structure,
  // because there's a non-trivial chance that it's now wrong.
  if (dp_ringInfo && dp_ringInfo->isInitialized() &&
      boost::out_degree(atomIdx1, d_graph) > 1 &&
      boost::out_degree(atomIdx2, d_graph) > 1) {
    dp_ringInfo->reset();
  }

  return numBonds;  // res;
}

unsigned int RWMol::addBond(Atom *atom1, Atom *atom2, Bond::BondType bondType) {
  PRECONDITION(atom1 && atom2, "NULL atom passed in");
  return addBond(atom1->getIdx(), atom2->getIdx(), bondType);
}

void RWMol::removeBond(unsigned int aid1, unsigned int aid2) {
  URANGE_CHECK(aid1, getNumAtoms());
  URANGE_CHECK(aid2, getNumAtoms());
  auto *bnd = getBondBetweenAtoms(aid1, aid2);
  if (!bnd) {
    return;
  }
  auto idx = bnd->getIdx();
  if (dp_delBonds) {
    // we're in a batch edit
    // if bonds have been added since we started, resize dp_delBonds
    if (dp_delBonds->size() < getNumBonds()) {
      dp_delBonds->resize(getNumBonds());
    }
    dp_delBonds->set(idx);
    return;
  }

  // remove any bookmarks which point to this bond:
  BOND_BOOKMARK_MAP *marks = getBondBookmarks();
  auto markI = marks->begin();
  while (markI != marks->end()) {
    BOND_PTR_LIST &bonds = markI->second;
    // we need to copy the iterator then increment it, because the
    // deletion we're going to do in clearBondBookmark will invalidate
    // it.
    auto tmpI = markI;
    ++markI;
    if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
      clearBondBookmark(tmpI->first, bnd);
    }
  }

  // loop over neighboring double bonds and remove their stereo atom
  //  information. This is definitely now invalid (was github issue 8)
  auto beginAtm = bnd->getBeginAtom();
  auto endAtm = bnd->getEndAtom();
  std::vector<std::vector<Atom *>> bond_atoms = {{beginAtm, endAtm},
                                                 {endAtm, beginAtm}};
  for (const auto &atoms : bond_atoms) {
    for (auto obnd : this->atomBonds(atoms[0])) {
      if (obnd == bnd) {
        continue;
      }
      if (std::find(obnd->getStereoAtoms().begin(),
                    obnd->getStereoAtoms().end(),
                    atoms[1]->getIdx()) != obnd->getStereoAtoms().end()) {
        // github #6900 if we remove stereo atoms we need to remove
        //  the CIS and or TRANS since this requires stereo atoms
        if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
            obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
          obnd->setStereo(Bond::BondStereo::STEREONONE);
        }

        obnd->getStereoAtoms().clear();
      }
    }
  }
  // reset our ring info structure, because it is pretty likely
  // to be wrong now:
  dp_ringInfo->reset();

  removeSubstanceGroupsReferencingBond(*this, idx);

  // loop over all bonds with higher indices and update their indices
  for (auto bond : bonds()) {
    if (bond->getIdx() > idx) {
      bond->setIdx(bond->getIdx() - 1);
    }
  }
  bnd->setOwningMol(nullptr);

  auto vd1 = boost::vertex(bnd->getBeginAtomIdx(), d_graph);
  auto vd2 = boost::vertex(bnd->getEndAtomIdx(), d_graph);
  boost::remove_edge(vd1, vd2, d_graph);
  delete bnd;
  --numBonds;
}

Bond *RWMol::createPartialBond(unsigned int atomIdx1, Bond::BondType bondType) {
  URANGE_CHECK(atomIdx1, getNumAtoms());

  auto *b = new Bond(bondType);
  b->setOwningMol(this);
  b->setBeginAtomIdx(atomIdx1);

  return b;
}

unsigned int RWMol::finishPartialBond(unsigned int atomIdx2, int bondBookmark,
                                      Bond::BondType bondType) {
  PRECONDITION(hasBondBookmark(bondBookmark), "no such partial bond");
  URANGE_CHECK(atomIdx2, getNumAtoms());

  auto *bsp = getBondWithBookmark(bondBookmark);
  if (bondType == Bond::UNSPECIFIED) {
    bondType = bsp->getBondType();
  }

  return addBond(bsp->getBeginAtomIdx(), atomIdx2, bondType);
}

void RWMol::beginBatchEdit() {
  if (dp_delAtoms || dp_delBonds) {
    throw ValueErrorException("Attempt to re-enter batchEdit mode");
  }
  dp_delAtoms.reset(new boost::dynamic_bitset<>(getNumAtoms()));
  dp_delBonds.reset(new boost::dynamic_bitset<>(getNumBonds()));
}

void RWMol::commitBatchEdit() {
  if (!(dp_delBonds || dp_delAtoms)) {
    return;
  } else if (dp_delBonds->none() && dp_delAtoms->none()) {
    // no need to reset ring info & calculated properties,
    // since nothing gets removed
    dp_delBonds.reset();
    dp_delAtoms.reset();
    return;
  }

  batchRemoveBonds();
  batchRemoveAtoms();

  // remove ring info
  dp_ringInfo->reset();

  // fix properties
  clearComputedProps(true);
  dp_delBonds.reset();
  dp_delAtoms.reset();
}

void RWMol::batchRemoveBonds() {
  if (!dp_delBonds || dp_delBonds->none()) {
    return;
  }

  auto &delBonds = *dp_delBonds;
  unsigned int min_idx = getNumBonds();
  for (unsigned int i = rdcast<unsigned int>(delBonds.size()); i > 0; --i) {
    if (!delBonds[i - 1]) continue;
    unsigned int idx = rdcast<unsigned int>(i - 1);
    Bond *bnd = getBondWithIdx(idx);
    if (!bnd) {
      continue;
    }

    min_idx = idx;
    // remove any bookmarks which point to this bond:
    BOND_BOOKMARK_MAP *marks = getBondBookmarks();
    auto markI = marks->begin();
    while (markI != marks->end()) {
      BOND_PTR_LIST &bonds = markI->second;
      // we need to copy the iterator then increment it, because the
      // deletion we're going to do in clearBondBookmark will invalidate
      // it.
      auto tmpI = markI;
      ++markI;
      if (std::find(bonds.begin(), bonds.end(), bnd) != bonds.end()) {
        clearBondBookmark(tmpI->first, bnd);
      }
    }

    // loop over neighboring double bonds and remove their stereo atom
    //  information. This is definitely now invalid (was github issue 8)
    auto beginAtm = bnd->getBeginAtom();
    auto endAtm = bnd->getEndAtom();
    std::vector<std::vector<Atom *>> bond_atoms = {{beginAtm, endAtm},
                                                   {endAtm, beginAtm}};
    for (const auto &atoms : bond_atoms) {
      for (auto obnd : atomBonds(atoms[0])) {
        if (obnd == bnd) {
          continue;
        }
        if (std::find(obnd->getStereoAtoms().begin(),
                      obnd->getStereoAtoms().end(),
                      atoms[1]->getIdx()) != obnd->getStereoAtoms().end()) {
          // github #6900 if we remove stereo atoms we need to remove
          //  the CIS and or TRANS since this requires stereo atoms
          if (obnd->getStereo() == Bond::BondStereo::STEREOCIS ||
              obnd->getStereo() == Bond::BondStereo::STEREOTRANS) {
            obnd->setStereo(Bond::BondStereo::STEREONONE);
          }
          obnd->getStereoAtoms().clear();
        }
      }
    }

    removeSubstanceGroupsReferencingBond(*this, idx);

    bnd->setOwningMol(nullptr);

    auto vd1 = boost::vertex(bnd->getBeginAtomIdx(), d_graph);
    auto vd2 = boost::vertex(bnd->getEndAtomIdx(), d_graph);
    boost::remove_edge(vd1, vd2, d_graph);
    delete bnd;
    --numBonds;
  }

  // loop over all bonds with higher indices than the minimum modified and
  // update their indices
  auto [firstB, lastB] = this->getEdges();
  unsigned int next_idx = min_idx;
  while (firstB != lastB) {
    Bond *bond = (*this)[*firstB];
    if (bond->getIdx() > min_idx) {
      bond->setIdx(next_idx++);
    }
    ++firstB;
  }
}

void RWMol::batchRemoveAtoms() {
  if (!dp_delAtoms || dp_delAtoms->none()) {
    return;
  }
  std::vector<Atom *> oldIndices(getNumAtoms());
  for (auto *atom : atoms()) {
    oldIndices[atom->getIdx()] = atom;
  }

  auto &delAtoms = *dp_delAtoms;
  for (unsigned int i = rdcast<unsigned int>(delAtoms.size()); i > 0; --i) {
    if (!delAtoms[i - 1]) continue;
    unsigned int idx = i - 1;
    Atom *atom = getAtomWithIdx(idx);
    if (!atom) continue;

    // remove any bookmarks which point to this atom:
    ATOM_BOOKMARK_MAP *marks = getAtomBookmarks();
    auto markI = marks->begin();
    while (markI != marks->end()) {
      const ATOM_PTR_LIST &atoms = markI->second;
      // we need to copy the iterator then increment it, because the
      // deletion we're going to do in clearAtomBookmark will invalidate
      // it.
      auto tmpI = markI;
      ++markI;
      if (std::find(atoms.begin(), atoms.end(), atom) != atoms.end()) {
        clearAtomBookmark(tmpI->first, atom);
      }
    }

    // now deal with bonds:
    //   their end indices may need to be decremented and their
    //   indices will need to be handled and if they have an
    //   ENDPTS prop that includes idx, it will need updating.
    EDGE_ITER beg, end;
    boost::tie(beg, end) = getEdges();
    std::string sprop;
    while (beg != end) {
      Bond *bond = d_graph[*beg++];
      if (bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
                                 sprop)) {
        // This would ideally use ParseV3000Array but I'm buggered if I can get
        // the linker to find it.
        //      std::vector<unsigned int> oats =
        //          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(sprop);
        if ('(' == sprop.front() && ')' == sprop.back()) {
          sprop = sprop.substr(1, sprop.length() - 2);

          // This is doing what ParseV3000Array would do.
          boost::char_separator<char> sep(" ");
          boost::tokenizer<boost::char_separator<char>> tokens(sprop, sep);
          unsigned int num_ats = std::stod(*tokens.begin());
          std::vector<unsigned int> oats;
          auto beg = tokens.begin();
          ++beg;
          std::transform(beg, tokens.end(), std::back_inserter(oats),
                         [](const std::string &a) { return std::stod(a); });

          auto idx_pos = std::find(oats.begin(), oats.end(), idx + 1);
          if (idx_pos != oats.end()) {
            oats.erase(idx_pos);
            --num_ats;
          }
          if (!num_ats) {
            bond->clearProp(RDKit::common_properties::_MolFileBondEndPts);
            bond->clearProp(common_properties::_MolFileBondAttach);
          } else {
            sprop = "(" + std::to_string(num_ats) + " ";
            for (auto &i : oats) {
              if (i > idx + 1) {
                --i;
              }
              sprop += std::to_string(i) + " ";
            }
            sprop[sprop.length() - 1] = ')';
            bond->setProp(RDKit::common_properties::_MolFileBondEndPts, sprop);
          }
        }
      }
    }

    removeSubstanceGroupsReferencingAtom(*this, idx);

    // Remove this atom from any stereo group
    removeAtomFromGroups(atom, d_stereo_groups);
    atom->setOwningMol(nullptr);

    // remove all connections to the atom:
    MolGraph::vertex_descriptor vd = boost::vertex(idx, d_graph);
    boost::clear_vertex(vd, d_graph);
    // finally remove the vertex itself
    boost::remove_vertex(vd, d_graph);
    delete atom;
    oldIndices[idx] = nullptr;
  }

  // reassign atom indices
  for (unsigned int i = 0; i < getNumAtoms(); i++) {
    Atom *atm = getAtomWithIdx(i);
    atm->setIdx(i);
  }

  // reassign atom indices in bonds
  for (auto bond : bonds()) {
    auto bgnidx = bond->getBeginAtomIdx();
    auto endidx = bond->getEndAtomIdx();
    Atom *bgn = oldIndices[bgnidx];
    Atom *end = oldIndices[endidx];
    CHECK_INVARIANT(bgn, "Atom mapping failed");
    CHECK_INVARIANT(end, "Atom mapping failed");
    bond->setBeginAtomIdx(bgn->getIdx());
    bond->setEndAtomIdx(end->getIdx());
    INT_VECT stereoAtoms;
    INT_VECT &oldStereoAtoms = bond->getStereoAtoms();
    if (oldStereoAtoms.size()) {
      for (auto &idx : oldStereoAtoms) {
        if (oldIndices[idx]) {
          stereoAtoms.push_back(oldIndices[idx]->getIdx());
        }
      }
      bond->getStereoAtoms().swap(stereoAtoms);
    }
  }

  // do the same with the coordinates in the conformations
  for (auto conf : d_confs) {
    RDGeom::POINT3D_VECT &positions = conf->getPositions();
    RDGeom::POINT3D_VECT newPositions;
    newPositions.reserve(getNumAtoms());

    for (RDGeom::POINT3D_VECT::size_type i = 0; i < positions.size(); ++i) {
      if (oldIndices[i] != nullptr) {
        newPositions.push_back(positions[i]);
      }
    }
    CHECK_INVARIANT(newPositions.size() == getNumAtoms(), "Lost coordinates!");
    positions.swap(newPositions);
  }
}

}  // namespace RDKit
