//
//  Copyright (C) 2001-2023 Greg Landrum and other RDKit contributors
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
#include <GraphMol/Chirality.h>
#include <GraphMol/RDKitQueries.h>

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
#include <boost/tokenizer.hpp>
#include <Geometry/point.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/FileParsers/MolSGroupParsing.h>

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

bool isMetal(const Atom *atom) {
  static const std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
  return q->Match(atom);
}

bool isHypervalentNonMetal(Atom *atom) {
  if (isMetal(atom)) {
    return false;
  }
  auto ev = atom->calcExplicitValence(false);
  // Check the explicit valence of the non-metal against the allowed
  // valences of the atom, adjusted by its formal charge.  This means that
  // N+ is treated the same as C, O+ the same as N.  This allows for,
  // for example, c1cccc[n+]1-[Fe] to be acceptable and not turned into
  // c1cccc[n+]1->[Fe].  After all, c1cccc[n+]1-C is ok.  Although this is
  // a poor example because c1ccccn1->[Fe] appears to be the normal
  // way that pyridine complexes with transition metals.  Heme b in
  // CHEBI:26355 is an example of when this is required.
  int effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
  if (effAtomicNum <= 0) {
    return false;
  }
  // atom is a non-metal, so if its atomic number is > 2, it should
  // obey the octet rule (2*ev <= 8).  If it doesn't, don't set the bond to
  // dative, so the molecule is later flagged as having bad valence.
  // [CH4]-[Na] being a case in point, line 1800 of
  // FileParsers/file_parsers_catch.cpp.
  const auto &otherValens =
      PeriodicTable::getTable()->getValenceList(effAtomicNum);
  if (otherValens.back() > 0 && ev > otherValens.back() && ev <= 4) {
    return true;
  }

  return false;
}

int numDativeBonds(const Atom *atom) {
  int numDatives = 0;
  auto &mol = atom->getOwningMol();
  for (auto bond : mol.atomBonds(atom)) {
    if (bond->getBondType() == Bond::BondType::DATIVE ||
        bond->getBondType() == Bond::BondType::DATIVEONE ||
        bond->getBondType() == Bond::BondType::DATIVEL ||
        bond->getBondType() == Bond::BondType::DATIVER) {
      ++numDatives;
    }
  }
  return numDatives;
}

// Returns true if the atom shouldn't do dative bonds.
bool noDative(const Atom *a) {
  static const std::set<int> noD{1, 2, 9, 10};
  return (noD.find(a->getAtomicNum()) != noD.end());
};

void metalBondCleanup(RWMol &mol, Atom *atom,
                      const std::vector<unsigned int> &ranks) {
  PRECONDITION(atom, "bad atom in metalBondCleanup");
  // The IUPAC recommendation for ligand->metal coordination bonds is that they
  // be single.  This upsets the RDKit valence model, as seen in CHEBI:26355,
  // heme b.  If the valence of a non-metal atom is above the maximum in the
  // RDKit model, and there are single bonds from it to metal
  // change those bonds to atom->metal dative.
  // If the atom is bonded to more than 1 metal atom, choose the one
  // with the fewer dative bonds incident on it, with the canonical
  // rank of the atoms as a tie-breaker.
  if (isHypervalentNonMetal(atom) && !noDative(atom)) {
    std::vector<Atom *> metals;
    // see if there are any metals bonded to it by a single bond
    for (auto bond : mol.atomBonds(atom)) {
      if (bond->getBondType() == Bond::BondType::SINGLE &&
          isMetal(bond->getOtherAtom(atom))) {
        metals.push_back(bond->getOtherAtom(atom));
      }
    }
    if (!metals.empty()) {
      std::sort(metals.begin(), metals.end(),
                [&](const Atom *a1, const Atom *a2) -> bool {
                  int nda1 = numDativeBonds(a1);
                  int nda2 = numDativeBonds(a2);
                  if (nda1 == nda2) {
                    return ranks[a1->getIdx()] > ranks[a2->getIdx()];
                  } else {
                    return nda1 < nda2;
                  }
                });
      auto bond =
          mol.getBondBetweenAtoms(atom->getIdx(), metals.front()->getIdx());
      if (bond) {
        bond->setBondType(RDKit::Bond::BondType::DATIVE);
        bond->setBeginAtom(atom);
        bond->setEndAtom(metals.front());
      }
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

void cleanUpOrganometallics(RWMol &mol) {
  // At present all this does is look for single bonds between
  // non-metals and metals where the non-metal exceeds one of
  // its normal valence states, and replaces that bond with
  // a dative one from the non-metal to the metal.
  bool needsFixing = false;
  for (const auto atom : mol.atoms()) {
    if (isHypervalentNonMetal(atom) && !noDative(atom)) {
      // see if there are any metals bonded to it by a single bond
      for (auto bond : mol.atomBonds(atom)) {
        if (bond->getBondType() == Bond::BondType::SINGLE &&
            isMetal(bond->getOtherAtom(atom))) {
          needsFixing = true;
          break;
        }
      }
    }
    if (needsFixing) {
      break;
    }
  }
  if (!needsFixing) {
    return;
  }

  mol.updatePropertyCache(false);
  // First see if anything needs doing
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  RDKit::Canon::rankMolAtoms(mol, ranks);
  std::vector<std::pair<int, int>> atom_ranks;
  for (size_t i = 0; i < ranks.size(); ++i) {
    atom_ranks.push_back(std::make_pair(i, ranks[i]));
  }
  std::sort(atom_ranks.begin(), atom_ranks.end(),
            [](const std::pair<int, int> &p1, std::pair<int, int> &p2) -> bool {
              return p1.second < p2.second;
            });
  for (auto ar : atom_ranks) {
    auto atom = mol.getAtomWithIdx(ar.first);
    metalBondCleanup(mol, atom, ranks);
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
    // we only automatically assign radicals to atoms that
    // don't have implicit Hs:
    if (!atom->getNoImplicit() || !atom->getAtomicNum()) {
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
      // #7122: if there's a bond to the metal center, then don't assign
      // radicals:
      if (atom->getDegree() > 0) {
        atom->setNumRadicalElectrons(0);
      } else {
        auto nValence = nOuter - chg;
        //  if this is an atom where we have no preferred valence info at all,
        //  e.g. for transition metals, then we shouldn't be guessing. This was
        //  #3330
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
}

MolOps::Hybridizations::Hybridizations(const ROMol &mol) {
  d_hybridizations.clear();
  // see if the mol already has computed hybridizations:

  if (mol.getNumAtoms() == 0) {
    return;
  }

  if ((*mol.atoms().begin())->getHybridization() !=
      Atom::HybridizationType::UNSPECIFIED) {
    for (auto atom : mol.atoms()) {
      d_hybridizations.push_back((int)atom->getHybridization());
    }
    return;
  }

  // compute them in a copy of the mol, so as not to change the mol passed in

  RWMol molCopy(mol);
  unsigned int operationThatFailed;
  unsigned int santitizeOps =
      MolOps::SANITIZE_SETCONJUGATION | MolOps::SANITIZE_SETHYBRIDIZATION;
  MolOps::sanitizeMol(molCopy, operationThatFailed, santitizeOps);
  for (auto atom : molCopy.atoms()) {
    // determine hybridization and remove chiral atoms that are not sp3
    d_hybridizations.push_back((int)atom->getHybridization());
  }
  return;
}

void cleanupAtropisomers(RWMol &mol) {
  auto hybs = MolOps::Hybridizations(mol);

  MolOps::cleanupAtropisomers(mol, hybs);
}

void cleanupAtropisomers(RWMol &mol, MolOps::Hybridizations &hybs) {
  // make sure that ring info is available
  // (defensive, current calls have it available)
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }
  const RingInfo *ri = mol.getRingInfo();
  for (auto bond : mol.bonds()) {
    switch (bond->getStereo()) {
      case Bond::BondStereo::STEREOATROPCW:
      case Bond::BondStereo::STEREOATROPCCW:
        if (hybs[bond->getBeginAtomIdx()] != Atom::SP2 ||
            hybs[bond->getEndAtomIdx()] != Atom::SP2 ||
            // do not clear bonds that part of a macrocycle
            // because they can be linking actual atropisomeric portions
            (ri->numBondRings(bond->getIdx()) > 0 &&
             ri->minBondRingSize(bond->getIdx()) < 8)) {
          bond->setStereo(Bond::BondStereo::STEREONONE);
        }
        break;

      default:
        break;
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

  // fix things like non-metal to metal bonds that should be dative.
  operationThatFailed = SANITIZE_CLEANUP_ORGANOMETALLICS;
  if (sanitizeOps & operationThatFailed) {
    cleanUpOrganometallics(mol);
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

  operationThatFailed = SANITIZE_CLEANUPATROPISOMERS;
  if (sanitizeOps & operationThatFailed) {
    cleanupAtropisomers(mol);
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
  std::unique_ptr<INT_VECT> mappingStorage;
  if (!frags) {
    mappingStorage.reset(new INT_VECT);
    frags = mappingStorage.get();
  }
  int nFrags = getMolFrags(mol, *frags);
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
    res.reserve(nFrags);
    for (int i = 0; i < nFrags; ++i) {
      boost::dynamic_bitset<> atomsInFrag(mol.getNumAtoms());
      INT_VECT comp;
      for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
        if ((*frags)[idx] == i) {
          comp.push_back(idx);
          atomsInFrag.set(idx);
        }
      }
      auto fragmentHasChallengingFeatures =
          [&](const INT_VECT &comp,
              const boost::dynamic_bitset<> &atomsInFrag) -> bool {
        for (auto idx : comp) {
          // check for atoms with stereochem:
          const auto atom = mol.getAtomWithIdx(idx);
          if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
              atom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
            return true;
          }
          for (auto bnd : mol.atomBonds(atom)) {
            if (atomsInFrag[bnd->getOtherAtomIdx(idx)]) {
              if (bnd->getStereo() != Bond::BondStereo::STEREONONE &&
                  bnd->getStereo() != Bond::BondStereo::STEREOANY) {
                return true;
              }
            }
          }
        }
        for (auto sgroup : getSubstanceGroups(mol)) {
          for (auto aid : sgroup.getAtoms()) {
            if (atomsInFrag[aid]) {
              return true;
            }
          }
          for (auto aid : sgroup.getParentAtoms()) {
            if (atomsInFrag[aid]) {
              return true;
            }
          }
        }
        // doesn't seem like this should be necessary, but in case
        // we ever need stereogroups where the atoms aren't marked
        // with stereo...
        for (auto stereoGroup : mol.getStereoGroups()) {
          for (auto atom : stereoGroup.getAtoms()) {
            if (atomsInFrag[atom->getIdx()]) {
              return true;
            }
          }
        }
        return false;
      };
      if (comp.size() == 1 ||
          (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
        // special case for a small, simple fragments when a bunch of fragments
        // are present. The check on the number of fragments is purely
        // empirical. This is mainly intended to catch situations like proteins
        // where you have a bunch of single-atom fragments (waters); the
        // standard approach below ends up being horribly inefficient there
        RWMOL_SPTR frag(new RWMol());
        res.push_back(frag);
        std::map<unsigned int, unsigned int> atomIdxMap;
        for (auto aid : comp) {
          atomIdxMap[aid] =
              frag->addAtom(mol.getAtomWithIdx(aid)->copy(), false, true);
        }
        for (auto bond : mol.bonds()) {
          if (atomsInFrag[bond->getBeginAtomIdx()] &&
              atomsInFrag[bond->getEndAtomIdx()]) {
            auto bondCopy = bond->copy();
            bondCopy->setBeginAtomIdx(atomIdxMap[bond->getBeginAtomIdx()]);
            bondCopy->setEndAtomIdx(atomIdxMap[bond->getEndAtomIdx()]);
            frag->addBond(bondCopy, true);
          }
        }
        if (copyConformers) {
          for (auto cit = mol.beginConformers(); cit != mol.endConformers();
               ++cit) {
            auto *conf = new Conformer(frag->getNumAtoms());
            conf->setId((*cit)->getId());
            conf->set3D((*cit)->is3D());
            unsigned int cidx = 0;
            for (auto ai : comp) {
              conf->setAtomPos(cidx++, (*cit)->getAtomPos(ai));
            }
            frag->addConformer(conf);
          }
        }
      } else {
        RWMOL_SPTR frag(new RWMol(mol));
        res.push_back(frag);
        frag->beginBatchEdit();
        for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
          if (!atomsInFrag[idx]) {
            frag->removeAtom(idx);
          }
        }
        frag->commitBatchEdit();
      }
      if (fragsMolAtomMapping) {
        (*fragsMolAtomMapping).push_back(comp);
      }
    }
  }
  if (!copyConformers) {
    for (auto &frag : res) {
      frag->clearConformers();
    }
  }

  if (sanitizeFrags) {
    for (auto &frag : res) {
      sanitizeMol(*frag);
    }
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

ROMol *hapticBondsToDative(const ROMol &mol) {
  auto *res = new RWMol(mol);
  hapticBondsToDative(*res);
  return static_cast<ROMol *>(res);
}

void hapticBondsToDative(RWMol &mol) {
  std::vector<unsigned int> dummiesToGo;
  std::vector<std::pair<unsigned int, unsigned int>> bondsToAdd;
  mol.beginBatchEdit();
  for (const auto &bond : mol.bonds()) {
    if (bond->getBondType() == Bond::BondType::DATIVE) {
      auto oats = details::hapticBondEndpoints(bond);
      if (oats.empty()) {
        continue;
      }
      Atom *dummy = nullptr;
      Atom *metal = nullptr;
      if (bond->getBeginAtom()->getAtomicNum() == 0) {
        dummy = bond->getBeginAtom();
        metal = bond->getEndAtom();
      } else if (bond->getEndAtom()->getAtomicNum() == 0) {
        metal = bond->getBeginAtom();
        dummy = bond->getEndAtom();
      }
      if (dummy == nullptr) {
        continue;
      }
      for (auto oat : oats) {
        auto atom = mol.getAtomWithIdx(oat);
        if (atom) {
          mol.addBond(atom, metal, Bond::DATIVE);
        }
      }
      mol.removeAtom(dummy);
    }
  }
  mol.commitBatchEdit();
}

ROMol *dativeBondsToHaptic(const ROMol &mol) {
  auto *res = new RWMol(mol);
  dativeBondsToHaptic(*res);
  return static_cast<ROMol *>(res);
}

namespace {
// return sets of contiguous atoms of more than 1 atom that are in
// allAts.
std::vector<std::vector<unsigned int>> contiguousAtoms(
    const ROMol &mol, const std::vector<unsigned int> &allAts) {
  std::vector<std::vector<unsigned int>> contigAts;
  std::vector<char> doneAts(mol.getNumAtoms(), 0);
  std::vector<char> inAllAts(mol.getNumAtoms(), 0);
  for (auto a : allAts) {
    inAllAts[a] = 1;
  }
  for (size_t i = 0; i < allAts.size(); ++i) {
    if (doneAts[allAts[i]]) {
      continue;
    }
    contigAts.push_back(std::vector<unsigned int>());
    std::list<const Atom *> toDo{mol.getAtomWithIdx(allAts[i])};
    while (!toDo.empty()) {
      auto nextAt = toDo.front();
      toDo.pop_front();
      if (!doneAts[nextAt->getIdx()]) {
        doneAts[nextAt->getIdx()] = 1;
        contigAts.back().push_back(nextAt->getIdx());
      }
      for (const auto &nbri :
           boost::make_iterator_range(mol.getAtomNeighbors(nextAt))) {
        if (inAllAts[nbri] && !doneAts[nbri]) {
          toDo.push_back(mol.getAtomWithIdx(nbri));
        }
      }
    }
    if (contigAts.back().size() < 2) {
      contigAts.pop_back();
    }
  }
  return contigAts;
}

// add to the molecule a dummy atom centred on the
// atoms passed in, with a dative bond from it to the metal atom.
void addHapticBond(RWMol &mol, unsigned int metalIdx,
                   std::vector<unsigned int> hapticAtoms) {
  // So there is a * in the V3000 file as the symbol for the atom.
  auto dummyAt = new QueryAtom(0);
  dummyAt->setQuery(makeAtomNullQuery());

  bool updateLabel = true;
  bool takeOwnwership = true;
  unsigned int dummyIdx = mol.addAtom(dummyAt, updateLabel, takeOwnwership);
  for (auto i = 0u; i < mol.getNumConformers(); ++i) {
    auto &conf = mol.getConformer(i);
    RDGeom::Point3D dummyPos;
    for (auto ha : hapticAtoms) {
      auto haPos = conf.getAtomPos(ha);
      dummyPos += haPos;
    }
    dummyPos /= hapticAtoms.size();
    conf.setAtomPos(dummyIdx, dummyPos);
  }
  unsigned int numbonds = mol.addBond(dummyIdx, metalIdx, Bond::DATIVE);
  auto bond = mol.getBondWithIdx(numbonds - 1);

  // Get the atom numbers for the end points.  First number is the
  // count, the rest count from 1.
  std::ostringstream oss;
  oss << "(" << hapticAtoms.size() << " ";
  for (auto ha : hapticAtoms) {
    oss << ha + 1 << " ";
  }
  std::string endpts{oss.str()};
  if (endpts.back() == ' ') {
    endpts = endpts.substr(0, endpts.length() - 1);
  }
  endpts += ")";
  bond->setProp(common_properties::_MolFileBondEndPts, endpts);
  bond->setProp<std::string>(common_properties::_MolFileBondAttach, "ALL");
}
}  // namespace

void dativeBondsToHaptic(RWMol &mol) {
  // First collect all the atoms that have a dative bond to them.
  // Assume that the ones of interest will have a metal as their
  // end atoms.
  std::map<unsigned int, std::vector<unsigned int>> dativeAtoms;
  for (const auto &b : mol.bonds()) {
    if (b->getBondType() == Bond::DATIVE) {
      auto ins = dativeAtoms.find(b->getEndAtomIdx());
      if (ins == dativeAtoms.end()) {
        dativeAtoms.insert(
            std::make_pair(b->getEndAtomIdx(),
                           std::vector<unsigned int>{b->getBeginAtomIdx()}));
      } else {
        ins->second.push_back(b->getBeginAtomIdx());
      }
    }
  }

  mol.beginBatchEdit();
  for (auto &dativeSet : dativeAtoms) {
    // Find the sets of contiguous atoms in the dativeAtoms lists.  Each one
    // will be the EndPts of a haptic bond going to the metal atom that is
    // dativeSet.first.
    auto contigAtoms = contiguousAtoms(mol, dativeSet.second);
    for (const auto &ca : contigAtoms) {
      addHapticBond(mol, dativeSet.first, ca);
      for (auto cat : ca) {
        mol.removeBond(dativeSet.first, cat);
      }
    }
  }
  mol.commitBatchEdit();
}

namespace details {
std::vector<int> hapticBondEndpoints(const Bond *bond) {
  // This would ideally use ParseV3000Array but I'm buggered if I can get
  // the linker to find it.  The issue, I think, is that it's in the
  // FileParsers library which is built after GraphMol so not available
  // to link in.  It can't be built first because it needs GraphMol.
  //      std::vector<unsigned int> oats =
  //          RDKit::SGroupParsing::ParseV3000Array<unsigned int>(endpts);
  // Returns the atom indices i.e. subtracts 1 from the numbers in the prop.
  std::vector<int> oats;
  std::string endpts;
  if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, endpts)) {
    if ('(' == endpts.front() && ')' == endpts.back()) {
      endpts = endpts.substr(1, endpts.length() - 2);
      boost::char_separator<char> sep(" ");
      boost::tokenizer<boost::char_separator<char>> tokens(endpts, sep);
      auto beg = tokens.begin();
      ++beg;
      std::transform(beg, tokens.end(), std::back_inserter(oats),
                     [](const std::string &a) { return std::stod(a) - 1; });
    }
  }
  return oats;
}
}  // end of namespace details

namespace details {
unsigned int addExplicitAttachmentPoint(RWMol &mol, unsigned int atomIdx,
                                        unsigned int val, bool addAsQuery,
                                        bool addCoords) {
  Atom *newAtom = nullptr;
  if (addAsQuery) {
    newAtom = new QueryAtom(0);
    newAtom->setQuery(RDKit::makeAtomNullQuery());
  } else {
    newAtom = new Atom(0);
  }
  newAtom->setProp(common_properties::_fromAttachPoint, val);
  bool updateLabel = false;
  bool takeOwnership = true;
  auto idx = mol.addAtom(newAtom, updateLabel, takeOwnership);
  mol.addBond(atomIdx, idx, Bond::SINGLE);
  mol.getAtomWithIdx(idx)->updatePropertyCache(false);
  if (addCoords) {
    setTerminalAtomCoords(mol, idx, atomIdx);
  }
  return idx;
}

bool isAttachmentPoint(const Atom *atom, bool markedOnly) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(atom->hasOwningMol(), "atom not associated with a molecule");
  if (atom->getAtomicNum() != 0 || atom->getDegree() != 1) {
    return false;
  }
  if (markedOnly && !atom->hasProp(common_properties::_fromAttachPoint)) {
    return false;
  }
  // we know that the atom is degree 1
  const auto bond = *atom->getOwningMol().atomBonds(atom).begin();
  if ((bond->getBondType() != Bond::BondType::SINGLE &&
       bond->getBondType() != Bond::BondType::UNSPECIFIED) ||
      bond->getBondDir() != Bond::BondDir::NONE) {
    return false;
  }

  if (atom->hasQuery()) {
    // a * from SMARTS
    if (!atom->getQuery()->getNegation() &&
        atom->getQuery()->getDescription() == "AtomNull") {
      return true;
    }
    // a * from CXSMILES
    if (atom->getQuery()->getNegation() &&
        atom->getQuery()->getDescription() == "AtomAtomicNum" &&
        static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() == 1) {
      return true;
    }
    return false;
  }

  return true;
}

}  // namespace details

void expandAttachmentPoints(RWMol &mol, bool addAsQueries, bool addCoords) {
  for (auto atom : mol.atoms()) {
    int value;
    if (atom->getPropIfPresent(common_properties::molAttachPoint, value)) {
      std::vector<int> tgtVals;
      if (value == 1 || value == -1) {
        tgtVals.push_back(1);
      }
      if (value == 2 || value == -1) {
        tgtVals.push_back(2);
      }
      if (tgtVals.empty()) {
        BOOST_LOG(rdWarningLog)
            << "Invalid value for molAttachPoint: " << value << " on atom "
            << atom->getIdx() << ". Not expanding this atttachment point."
            << std::endl;
        continue;
      }
      for (auto tval : tgtVals) {
        atom->clearProp(common_properties::molAttachPoint);
        details::addExplicitAttachmentPoint(mol, atom->getIdx(), tval,
                                            addAsQueries, addCoords);
      }
    }
  }
}

void collapseAttachmentPoints(RWMol &mol, bool markedOnly) {
  bool removedAny = false;
  std::vector<int> attachLabels(mol.getNumAtoms(), 0);

  for (auto atom : mol.atoms()) {
    if (details::isAttachmentPoint(atom, markedOnly)) {
      int value = 0;
      atom->getPropIfPresent(common_properties::_fromAttachPoint, value);
      if (markedOnly && (value < 0 || value > 2)) {
        BOOST_LOG(rdWarningLog)
            << "Invalid value for _fromAttachPoint: " << value << " on atom "
            << atom->getIdx() << ". Not collapsing this atom" << std::endl;
        continue;
      }
      if (!markedOnly && !value) {
        value = 1;
      }
      auto bond = *mol.atomBonds(atom).begin();
      if ((bond->getBondType() != Bond::BondType::SINGLE &&
           bond->getBondType() != Bond::BondType::UNSPECIFIED) ||
          bond->getBondDir() != Bond::BondDir::NONE) {
        continue;
      }
      auto oAtomIdx = bond->getOtherAtom(atom)->getIdx();
      if (attachLabels[oAtomIdx]) {
        if (attachLabels[oAtomIdx] != -1) {
          value = -1;
        } else {
          BOOST_LOG(rdWarningLog)
              << "More than two attachment points on atom " << oAtomIdx
              << ". Attachment point " << atom->getIdx()
              << " will not be collapsed." << std::endl;
          continue;
        }
      }
      if (!removedAny) {
        mol.beginBatchEdit();
        removedAny = true;
      }
      attachLabels[oAtomIdx] = value;
      mol.removeAtom(atom);
    }
  }
  // set the attachment point labels
  for (auto atom : mol.atoms()) {
    if (attachLabels[atom->getIdx()]) {
      atom->setProp(common_properties::molAttachPoint,
                    attachLabels[atom->getIdx()]);
    }
  }
  if (removedAny) {
    mol.commitBatchEdit();
  }
}
}  // end of namespace MolOps
}  // end of namespace RDKit
