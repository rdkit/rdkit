//
//  Copyright (C) 2001-2017 Greg Landrum and Rational Discovery LLC
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>

#include <vector>
#include <algorithm>

#include <RDGeneral/BoostStartInclude.h>

#include <boost/graph/connected_components.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <RDGeneral/BoostEndInclude.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <GraphMol/ROMol.h>

const int ci_LOCAL_INF = static_cast<int>(1e8);

namespace RDKit {
namespace MolOps {
namespace {
void nitrogenCleanup(RWMol &mol, Atom *atom) {
  // conversions here:
  // - neutral 5 coordinate Ns with double bonds to Os to the
  //   zwitterionic form.  e.g.:
  //   CN(=O)=O -> C[N+](=O)[O-]
  //   and:
  //   C1=CC=CN(=O)=C1 -> C1=CC=C[N+]([O-])=C1
  // - neutral 5 coordinate Ns with triple bonds to Ns to the
  //   zwitterionic form.  e.g.:
  //   C-N=N#N -> C-N=[N+]=[N-]

  PRECONDITION(atom, "bad atom");
  bool aromHolder;

  // we only want to do neutrals so that things like this don't get
  // munged:
  //  O=[n+]1occcc1
  // this was sf.net issue 1811276
  if (atom->getFormalCharge()) {
    return;
  }

  // we need to play this little aromaticity game because the
  // explicit valence code modifies its results for aromatic
  // atoms.
  aromHolder = atom->getIsAromatic();
  atom->setIsAromatic(0);
  // NOTE that we are calling calcExplicitValence() here, we do
  // this because we cannot be sure that it has already been
  // called on the atom (cleanUp() gets called pretty early in
  // the sanitization process):
  if (atom->calcExplicitValence(false) == 5) {
    unsigned int aid = atom->getIdx();
    for (const auto nbr : mol.atomNeighbors(atom)) {
      if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
          (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
           Bond::DOUBLE)) {
        // here's the double bonded oxygen
        auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        b->setBondType(Bond::SINGLE);
        atom->setFormalCharge(1);
        nbr->setFormalCharge(-1);
        break;
      } else if ((nbr->getAtomicNum() == 7) && (nbr->getFormalCharge() == 0) &&
                 (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
                  Bond::TRIPLE)) {
        // here's the triple bonded nitrogen
        auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        b->setBondType(Bond::DOUBLE);
        atom->setFormalCharge(1);
        nbr->setFormalCharge(-1);
        break;
      }
    }  // end of loop over the first neigh
  }    // if this atom is 5 coordinate nitrogen
  // force a recalculation of the explicit valence here
  atom->setIsAromatic(aromHolder);
  atom->calcExplicitValence(false);
}

void phosphorusCleanup(RWMol &mol, Atom *atom) {
  // conversions here:
  // - neutral 5 coordinate Ps with one double bonds to an Os
  //   and one to a C or N to the zwitterionic form.  e.g.:
  //   C=P(=O)X -> C=[P+]([O-])X
  PRECONDITION(atom, "bad atom");

  // we only want to do neutrals
  if (atom->getFormalCharge()) {
    return;
  }

  // NOTE that we are calling calcExplicitValence() here, we do
  // this because we cannot be sure that it has already been
  // called on the atom (cleanUp() gets called pretty early in
  // the sanitization process):
  if (atom->calcExplicitValence(false) == 5 && atom->getDegree() == 3) {
    unsigned int aid = atom->getIdx();
    Bond *dbl_to_O = nullptr;
    Atom *O_atom = nullptr;
    bool hasDoubleToCorN = false;
    for (const auto nbr : mol.atomNeighbors(atom)) {
      if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
          (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
           Bond::DOUBLE)) {
        // here's the double bonded oxygen
        dbl_to_O = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        O_atom = nbr;
      } else if ((nbr->getAtomicNum() == 6 || nbr->getAtomicNum() == 7) &&
                 (nbr->getDegree() >= 2) &&
                 (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
                  Bond::DOUBLE)) {
        hasDoubleToCorN = true;
      }
    }  // end of loop over the first neigh
    if (hasDoubleToCorN && dbl_to_O != nullptr) {
      TEST_ASSERT(O_atom != nullptr);
      O_atom->setFormalCharge(-1);
      dbl_to_O->setBondType(Bond::SINGLE);
      atom->setFormalCharge(1);
    }
  }
  // force a recalculation of the explicit valence here
  atom->calcExplicitValence(false);
}

void halogenCleanup(RWMol &mol, Atom *atom) {
  PRECONDITION(atom, "bad atom");
  // Conversions done:
  //    X(=O)(=O)(=O)O -> [X+3]([O-])([O-])([O-])O
  //    X(=O)(=O)O -> [X+2]([O-])([O-])O
  //    X(=O)O -> [X+]([O-])O
  int ev = atom->calcExplicitValence(false);
  if (atom->getFormalCharge() == 0 && (ev == 7 || ev == 5 || ev == 3)) {
    bool neighborsAllO = true;
    for (const auto nbr : mol.atomNeighbors(atom)) {
      if (nbr->getAtomicNum() != 8) {
        neighborsAllO = false;
        break;
      }
    }
    if (neighborsAllO) {
      int formalCharge = 0;
      for (auto bond : mol.atomBonds(atom)) {
        if (bond->getBondType() == Bond::DOUBLE) {
          bond->setBondType(Bond::SINGLE);
          auto otherAtom = bond->getOtherAtom(atom);
          formalCharge++;
          otherAtom->setFormalCharge(-1);
          otherAtom->calcExplicitValence(false);
        }
      }
      atom->setFormalCharge(formalCharge);
      atom->calcExplicitValence(false);
    }
  }
}
}  // namespace

void cleanUp(RWMol &mol) {
  for (auto atom : mol.atoms()) {
    switch (atom->getAtomicNum()) {
      case 7:
        nitrogenCleanup(mol, atom);
        break;
      case 15:
        phosphorusCleanup(mol, atom);
        break;
      case 17:
      case 35:
      case 53:
        halogenCleanup(mol, atom);
        break;
    }
  }
}

void adjustHs(RWMol &mol) {
  //
  //  Go through and adjust the number of implicit and explicit Hs
  //  on each atom in the molecule.
  //
  //  Atoms that do not *need* explicit Hs
  //
  //  Assumptions: this is called after the molecule has been
  //  sanitized, aromaticity has been perceived, and the implicit
  //  valence of everything has been calculated.
  //
  for (auto atom : mol.atoms()) {
    int origImplicitV = atom->getImplicitValence();
    atom->calcExplicitValence(false);
    int origExplicitV = atom->getNumExplicitHs();

    int newImplicitV = atom->calcImplicitValence(false);
    //
    //  Case 1: The disappearing Hydrogen
    //    Smiles:  O=C1NC=CC2=C1C=CC=C2
    //
    //    after perception is done, the N atom has two aromatic
    //    bonds to it and a single implicit H.  When the Smiles is
    //    written, we get: n1ccc2ccccc2c1=O.  Here the nitrogen has
    //    no implicit Hs (because there are two aromatic bonds to
    //    it, giving it a valence of 3).  Also: this SMILES is bogus
    //    (un-kekulizable).  The correct SMILES would be:
    //    [nH]1ccc2ccccc2c1=O.  So we need to loop through the atoms
    //    and find those that have lost implicit H; we'll add those
    //    back as explicit Hs.
    //
    //    <phew> that takes way longer to comment than it does to
    //    write:
    if (newImplicitV < origImplicitV) {
      atom->setNumExplicitHs(origExplicitV + (origImplicitV - newImplicitV));
      atom->calcExplicitValence(false);
    }
  }
}

void assignRadicals(RWMol &mol) {
  for (auto atom : mol.atoms()) {
    // we only put automatically assign radicals to things that
    // don't have them already and don't have implicit Hs:
    if (!atom->getNoImplicit() || atom->getNumRadicalElectrons() ||
        !atom->getAtomicNum()) {
      continue;
    }
    const auto &valens =
        PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
    int chg = atom->getFormalCharge();
    int nOuter =
        PeriodicTable::getTable()->getNouterElecs(atom->getAtomicNum());
    if (valens.size() != 1 || valens[0] != -1) {
      double accum = 0.0;
      RWMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = mol.getAtomBonds(atom);
      while (beg != end) {
        accum += mol[*beg]->getValenceContrib(atom);
        ++beg;
      }
      accum += atom->getNumExplicitHs();
      int totalValence = static_cast<int>(accum + 0.1);
      int baseCount = 8;
      if (atom->getAtomicNum() == 1 || atom->getAtomicNum() == 2) {
        baseCount = 2;
      }

      // applies to later (more electronegative) elements:
      int numRadicals = baseCount - nOuter - totalValence + chg;
      if (numRadicals < 0) {
        numRadicals = 0;
        // can the atom be "hypervalent"?  (was github #447)
        const INT_VECT &valens =
            PeriodicTable::getTable()->getValenceList(atom->getAtomicNum());
        if (valens.size() > 1) {
          for (auto val : valens) {
            if (val - totalValence + chg >= 0) {
              numRadicals = val - totalValence + chg;
              break;
            }
          }
        }
      }
      // applies to earlier elements:
      int numRadicals2 = nOuter - totalValence - chg;
      if (numRadicals2 >= 0) {
        numRadicals = std::min(numRadicals, numRadicals2);
      }
      atom->setNumRadicalElectrons(numRadicals);
    } else {
      //  if this is an atom where we have no preferred valence info at all,
      //  e.g. for transition metals, then we shouldn't be guessing. This was
      //  #3330
      auto nValence = nOuter - chg;
      if (nValence < 0) {
        // this was github #5462
        nValence = 0;
        BOOST_LOG(rdWarningLog)
            << "Unusual charge on atom " << atom->getIdx()
            << " number of radical electrons set to zero" << std::endl;
      }
      atom->setNumRadicalElectrons(nValence % 2);
    }
  }
}

void sanitizeMol(RWMol &mol) {
  unsigned int failedOp = 0;
  sanitizeMol(mol, failedOp, SANITIZE_ALL);
}
void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
                 unsigned int sanitizeOps) {
  // clear out any cached properties
  mol.clearComputedProps();

  operationThatFailed = SANITIZE_CLEANUP;
  if (sanitizeOps & operationThatFailed) {
    // clean up things like nitro groups
    cleanUp(mol);
  }

  // update computed properties on atoms and bonds:
  operationThatFailed = SANITIZE_PROPERTIES;
  if (sanitizeOps & operationThatFailed) {
    mol.updatePropertyCache(true);
  } else {
    mol.updatePropertyCache(false);
  }

  operationThatFailed = SANITIZE_SYMMRINGS;
  if (sanitizeOps & operationThatFailed) {
    VECT_INT_VECT arings;
    MolOps::symmetrizeSSSR(mol, arings);
  }

  // kekulizations
  operationThatFailed = SANITIZE_KEKULIZE;
  if (sanitizeOps & operationThatFailed) {
    Kekulize(mol);
  }

  // look for radicals:
  // We do this now because we need to know
  // that the N in [N]1C=CC=C1 has a radical
  // before we move into setAromaticity().
  // It's important that this happen post-Kekulization
  // because there's no way of telling what to do
  // with the same molecule if it's in the form
  // [n]1cccc1
  operationThatFailed = SANITIZE_FINDRADICALS;
  if (sanitizeOps & operationThatFailed) {
    assignRadicals(mol);
  }

  // then do aromaticity perception
  operationThatFailed = SANITIZE_SETAROMATICITY;
  if (sanitizeOps & operationThatFailed) {
    setAromaticity(mol);
  }

  // set conjugation
  operationThatFailed = SANITIZE_SETCONJUGATION;
  if (sanitizeOps & operationThatFailed) {
    setConjugation(mol);
  }

  // set hybridization
  operationThatFailed = SANITIZE_SETHYBRIDIZATION;
  if (sanitizeOps & operationThatFailed) {
    setHybridization(mol);
  }

  // remove bogus chirality specs:
  operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
  if (sanitizeOps & operationThatFailed) {
    cleanupChirality(mol);
  }

  // adjust Hydrogen counts:
  operationThatFailed = SANITIZE_ADJUSTHS;
  if (sanitizeOps & operationThatFailed) {
    adjustHs(mol);
  }

  // now that everything has been cleaned up, go through and check/update the
  // computed valences on atoms and bonds one more time
  operationThatFailed = SANITIZE_PROPERTIES;
  if (sanitizeOps & operationThatFailed) {
    mol.updatePropertyCache(true);
  }
  operationThatFailed = 0;
}

std::vector<std::unique_ptr<MolSanitizeException>> detectChemistryProblems(
    const ROMol &imol, unsigned int sanitizeOps) {
  RWMol mol(imol);
  std::vector<std::unique_ptr<MolSanitizeException>> res;

  // clear out any cached properties
  mol.clearComputedProps();

  int operation;
  operation = SANITIZE_CLEANUP;
  if (sanitizeOps & operation) {
    // clean up things like nitro groups
    cleanUp(mol);
  }

  // update computed properties on atoms and bonds:
  operation = SANITIZE_PROPERTIES;
  if (sanitizeOps & operation) {
    for (auto &atom : mol.atoms()) {
      try {
        bool strict = true;
        atom->updatePropertyCache(strict);
      } catch (const MolSanitizeException &e) {
        res.emplace_back(e.copy());
      }
    }
  } else {
    mol.updatePropertyCache(false);
  }

  // kekulizations
  operation = SANITIZE_KEKULIZE;
  if (sanitizeOps & operation) {
    try {
      Kekulize(mol);
    } catch (const MolSanitizeException &e) {
      res.emplace_back(e.copy());
    }
  }
  return res;
}

std::vector<ROMOL_SPTR> getMolFrags(const ROMol &mol, bool sanitizeFrags,
                                    INT_VECT *frags,
                                    VECT_INT_VECT *fragsMolAtomMapping,
                                    bool copyConformers) {
  bool ownIt = false;
  INT_VECT *mapping;
  if (frags) {
    mapping = frags;
  } else {
    mapping = new INT_VECT;
    ownIt = true;
  }
  unsigned int nFrags = getMolFrags(mol, *mapping);
  std::vector<RWMOL_SPTR> res;
  if (nFrags == 1) {
    auto *tmp = new RWMol(mol);
    RWMOL_SPTR sptr(tmp);
    res.push_back(sptr);
    if (fragsMolAtomMapping) {
      INT_VECT comp;
      for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
        comp.push_back(idx);
      }
      (*fragsMolAtomMapping).push_back(comp);
    }
  } else {
    std::vector<int> ids(mol.getNumAtoms(), -1);
    boost::dynamic_bitset<> copiedAtoms(mol.getNumAtoms(), 0);
    boost::dynamic_bitset<> copiedBonds(mol.getNumBonds(), 0);
    res.reserve(nFrags);
    for (unsigned int frag = 0; frag < nFrags; ++frag) {
      auto *tmp = new RWMol();
      RWMOL_SPTR sptr(tmp);
      res.push_back(sptr);
    }

    // copy atoms
    INT_INT_VECT_MAP comMap;
    for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
      const Atom *oAtm = mol.getAtomWithIdx(idx);
      ids[idx] = res[(*mapping)[idx]]->addAtom(oAtm->copy(), false, true);
      copiedAtoms[idx] = 1;
      if (fragsMolAtomMapping) {
        if (comMap.find((*mapping)[idx]) == comMap.end()) {
          INT_VECT comp;
          comMap[(*mapping)[idx]] = comp;
        }
        comMap[(*mapping)[idx]].push_back(idx);
      }
      // loop over neighbors and add bonds in the fragment to all atoms
      // that are already in the same fragment
      for (const auto nbr : mol.atomNeighbors(oAtm)) {
        if (copiedAtoms[nbr->getIdx()]) {
          copiedBonds[mol.getBondBetweenAtoms(idx, nbr->getIdx())->getIdx()] =
              1;
        }
      }
    }
    // update ring stereochemistry information
    for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
      const Atom *oAtm = mol.getAtomWithIdx(idx);
      INT_VECT ringStereoAtomsMol;
      if (oAtm->getPropIfPresent(common_properties::_ringStereoAtoms,
                                 ringStereoAtomsMol)) {
        INT_VECT ringStereoAtomsCopied;
        for (int rnbr : ringStereoAtomsMol) {
          int ori_ridx = abs(rnbr) - 1;
          int ridx = ids[ori_ridx] + 1;
          if (rnbr < 0) {
            ridx *= (-1);
          }
          ringStereoAtomsCopied.push_back(ridx);
        }
        res[(*mapping)[idx]]->getAtomWithIdx(ids[idx])->setProp(
            common_properties::_ringStereoAtoms, ringStereoAtomsCopied);
      }
    }

    // copy bonds and bond stereochemistry information
    ROMol::EDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getEdges();
    while (beg != end) {
      const Bond *bond = (mol)[*beg];
      ++beg;
      if (!copiedBonds[bond->getIdx()]) {
        continue;
      }
      Bond *nBond = bond->copy();
      RWMol *tmp = res[(*mapping)[nBond->getBeginAtomIdx()]].get();
      nBond->setOwningMol(tmp);
      nBond->setBeginAtomIdx(ids[nBond->getBeginAtomIdx()]);
      nBond->setEndAtomIdx(ids[nBond->getEndAtomIdx()]);
      nBond->getStereoAtoms().clear();
      INT_VECT stereoAtoms = bond->getStereoAtoms();
      for (int stereoAtom : stereoAtoms) {
        nBond->getStereoAtoms().push_back(ids[stereoAtom]);
      }
      tmp->addBond(nBond, true);
    }

    // copy RingInfo
    if (mol.getRingInfo()->isInitialized()) {
      for (const auto &i : mol.getRingInfo()->atomRings()) {
        INT_VECT aids;
        auto tmp = res[(*mapping)[i[0]]].get();
        if (!tmp->getRingInfo()->isInitialized()) {
          tmp->getRingInfo()->initialize();
        }
        for (int j : i) {
          aids.push_back(ids[j]);
        }
        INT_VECT bids;
        INT_VECT_CI lastRai = aids.begin();
        for (INT_VECT_CI rai = aids.begin() + 1; rai != aids.end(); ++rai) {
          const Bond *bnd = tmp->getBondBetweenAtoms(*rai, *lastRai);
          if (!bnd) {
            throw ValueErrorException("expected bond not found");
          }
          bids.push_back(bnd->getIdx());
          lastRai = rai;
        }
        const Bond *bnd = tmp->getBondBetweenAtoms(*lastRai, *(aids.begin()));
        if (!bnd) {
          throw ValueErrorException("expected bond not found");
        }
        bids.push_back(bnd->getIdx());
        tmp->getRingInfo()->addRing(aids, bids);
      }
    }

    if (copyConformers) {
      // copy conformers
      for (auto cit = mol.beginConformers(); cit != mol.endConformers();
           ++cit) {
        for (auto &re : res) {
          ROMol *newM = re.get();
          auto *conf = new Conformer(newM->getNumAtoms());
          conf->setId((*cit)->getId());
          conf->set3D((*cit)->is3D());
          newM->addConformer(conf);
        }
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
          if (ids[i] < 0) {
            continue;
          }
          res[(*mapping)[i]]
              ->getConformer((*cit)->getId())
              .setAtomPos(ids[i], (*cit)->getAtomPos(i));
        }
      }
    }

    if (fragsMolAtomMapping) {
      for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end();
           mci++) {
        (*fragsMolAtomMapping).push_back((*mci).second);
      }
    }
  }

  if (sanitizeFrags) {
    for (auto &re : res) {
      sanitizeMol(*re);
    }
  }

  if (ownIt) {
    delete mapping;
  }
  return std::vector<ROMOL_SPTR>(res.begin(), res.end());
}

unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
  unsigned int natms = mol.getNumAtoms();
  mapping.resize(natms);
  return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
               : 0;
};

unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
  frags.clear();
  INT_VECT mapping;
  getMolFrags(mol, mapping);

  INT_INT_VECT_MAP comMap;
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    int mi = mapping[i];
    if (comMap.find(mi) == comMap.end()) {
      INT_VECT comp;
      comMap[mi] = comp;
    }
    comMap[mi].push_back(i);
  }

  for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    frags.push_back((*mci).second);
  }
  return rdcast<unsigned int>(frags.size());
}

template <typename T>
std::map<T, boost::shared_ptr<ROMol>> getMolFragsWithQuery(
    const ROMol &mol, T (*query)(const ROMol &, const Atom *),
    bool sanitizeFrags, const std::vector<T> *whiteList, bool negateList) {
  PRECONDITION(query, "no query");

  std::vector<T> assignments(mol.getNumAtoms());
  std::vector<int> ids(mol.getNumAtoms(), -1);
  std::map<T, boost::shared_ptr<ROMol>> res;
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    T where = query(mol, mol.getAtomWithIdx(i));
    if (whiteList) {
      bool found = std::find(whiteList->begin(), whiteList->end(), where) !=
                   whiteList->end();
      if (!found && !negateList) {
        continue;
      } else if (found && negateList) {
        continue;
      }
    }
    assignments[i] = where;
    if (res.find(where) == res.end()) {
      res[where] = boost::shared_ptr<ROMol>(new ROMol());
    }
    auto *frag = static_cast<RWMol *>(res[where].get());
    ids[i] = frag->addAtom(mol.getAtomWithIdx(i)->copy(), false, true);
    // loop over neighbors and add bonds in the fragment to all atoms
    // that are already in the same fragment
    ROMol::ADJ_ITER nbrIdx, endNbrs;
    boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(i));
    while (nbrIdx != endNbrs) {
      if (*nbrIdx < i && assignments[*nbrIdx] == where) {
        Bond *nBond = mol.getBondBetweenAtoms(i, *nbrIdx)->copy();
        nBond->setOwningMol(static_cast<ROMol *>(frag));
        nBond->setBeginAtomIdx(ids[nBond->getBeginAtomIdx()]);
        nBond->setEndAtomIdx(ids[nBond->getEndAtomIdx()]);
        frag->addBond(nBond, true);
      }
      ++nbrIdx;
    }
  }
  // update conformers
  for (auto cit = mol.beginConformers(); cit != mol.endConformers(); ++cit) {
    for (auto iter = res.begin(); iter != res.end(); ++iter) {
      ROMol *newM = iter->second.get();
      auto *conf = new Conformer(newM->getNumAtoms());
      conf->setId((*cit)->getId());
      conf->set3D((*cit)->is3D());
      newM->addConformer(conf);
    }
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (ids[i] < 0) {
        continue;
      }
      res[assignments[i]]
          ->getConformer((*cit)->getId())
          .setAtomPos(ids[i], (*cit)->getAtomPos(i));
    }
  }
  if (sanitizeFrags) {
    for (auto iter = res.begin(); iter != res.end(); ++iter) {
      sanitizeMol(*static_cast<RWMol *>(iter->second.get()));
    }
  }
  return res;
}
template RDKIT_GRAPHMOL_EXPORT std::map<std::string, boost::shared_ptr<ROMol>>
getMolFragsWithQuery(const ROMol &mol,
                     std::string (*query)(const ROMol &, const Atom *),
                     bool sanitizeFrags, const std::vector<std::string> *,
                     bool);
template RDKIT_GRAPHMOL_EXPORT std::map<int, boost::shared_ptr<ROMol>>
getMolFragsWithQuery(const ROMol &mol,
                     int (*query)(const ROMol &, const Atom *),
                     bool sanitizeFrags, const std::vector<int> *, bool);
template RDKIT_GRAPHMOL_EXPORT std::map<unsigned int, boost::shared_ptr<ROMol>>
getMolFragsWithQuery(const ROMol &mol,
                     unsigned int (*query)(const ROMol &, const Atom *),
                     bool sanitizeFrags, const std::vector<unsigned int> *,
                     bool);

#if 0
    void findSpanningTree(const ROMol &mol,INT_VECT &mst){
      //
      //  The BGL provides Prim's and Kruskal's algorithms for finding
      //  the MST of a graph.  Prim's is O(n2) (n=# of atoms) while
      //  Kruskal's is O(e log e) (e=# of bonds).  For molecules, where
      //  e << n2, Kruskal's should be a win.
      //
      const MolGraph *mgraph = &mol.getTopology();
      MolGraph *molGraph = const_cast<MolGraph *> (mgraph);

      std::vector<MolGraph::edge_descriptor> treeEdges;
      treeEdges.reserve(boost::num_vertices(*molGraph));

      boost::property_map < MolGraph, edge_wght_t >::type w = boost::get(edge_wght_t(), *molGraph);
      boost::property_map < MolGraph, edge_bond_t>::type bps = boost::get(edge_bond_t(), *molGraph);
      boost::graph_traits < MolGraph >::edge_iterator e, e_end;
      Bond* bnd;
      for (boost::tie(e, e_end) = boost::edges(*molGraph); e != e_end; ++e) {
        bnd = bps[*e];

        if(!bnd->getIsAromatic()){
          w[*e] = (bnd->getBondTypeAsDouble());
        } else {
          w[*e] = 3.0/2.0;
        }
      }

      // FIX: this is a hack due to problems with MSVC++
#if 1
      typedef boost::graph_traits<MolGraph>::vertices_size_type size_type;
      typedef boost::graph_traits<MolGraph>::vertex_descriptor vertex_t;
      typedef boost::property_map<MolGraph,boost::vertex_index_t>::type index_map_t;
      boost::graph_traits<MolGraph>::vertices_size_type
        n = boost::num_vertices(*molGraph);
      std::vector<size_type> rank_map(n);
      std::vector<vertex_t> pred_map(n);

      boost::detail::kruskal_mst_impl
        (*molGraph, std::back_inserter(treeEdges),
         boost::make_iterator_property_map(rank_map.begin(),
                                           boost::get(boost::vertex_index, *molGraph),
                                           rank_map[0]),
         boost::make_iterator_property_map(pred_map.begin(),
                                           boost::get(boost::vertex_index, *molGraph),
                                           pred_map[0]),
         w);

#else
      boost::kruskal_minimum_spanning_tree(*molGraph,std::back_inserter(treeEdges),
                                           w, *molGraph);
      //boost::weight_map(static_cast<boost::property_map<MolGraph,edge_wght_t>::const_type>(boost::get(edge_wght_t(),*molGraph))));
#endif
      mst.resize(0);
      for(std::vector<MolGraph::edge_descriptor>::iterator edgeIt=treeEdges.begin();
          edgeIt!=treeEdges.end();edgeIt++){
        mst.push_back(mol[*edgeIt]->getIdx());
      }
    }
#endif
int getFormalCharge(const ROMol &mol) {
  int accum = 0;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    accum += (*atomIt)->getFormalCharge();
  }
  return accum;
};

unsigned getNumAtomsWithDistinctProperty(const ROMol &mol, std::string prop) {
  unsigned numPropAtoms = 0;
  for (const auto atom : mol.atoms()) {
    if (atom->hasProp(prop)) {
      ++numPropAtoms;
    }
  }
  return numPropAtoms;
}
};  // end of namespace MolOps
};  // end of namespace RDKit
