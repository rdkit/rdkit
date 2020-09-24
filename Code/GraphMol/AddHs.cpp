//
//  Copyright (C) 2003-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDKitBase.h"
#include <list>
#include "QueryAtom.h"
#include "QueryOps.h"
#include "MonomerInfo.h"
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {

// Local utility functionality:
namespace {
Atom *getAtomNeighborNot(ROMol *mol, const Atom *atom, const Atom *other) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  PRECONDITION(atom->getDegree() > 1, "bad degree");
  PRECONDITION(other, "bad atom");
  Atom *res = nullptr;

  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(atom);
  while (nbrIdx != endNbrs) {
    if (*nbrIdx != other->getIdx()) {
      res = mol->getAtomWithIdx(*nbrIdx);
      break;
    }
    ++nbrIdx;
  }

  POSTCONDITION(res, "no neighbor found");
  return res;
}

void setHydrogenCoords(ROMol *mol, unsigned int hydIdx, unsigned int heavyIdx) {
  // we will loop over all the coordinates
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(heavyIdx != hydIdx, "degenerate atoms");
  Atom *hydAtom = mol->getAtomWithIdx(hydIdx);
  PRECONDITION(mol->getAtomDegree(hydAtom) == 1, "bad atom degree");
  const Bond *bond = mol->getBondBetweenAtoms(heavyIdx, hydIdx);
  PRECONDITION(bond, "no bond between atoms");

  const Atom *heavyAtom = mol->getAtomWithIdx(heavyIdx);
  double bondLength =
      PeriodicTable::getTable()->getRb0(1) +
      PeriodicTable::getTable()->getRb0(heavyAtom->getAtomicNum());

  RDGeom::Point3D dirVect(0, 0, 0);

  RDGeom::Point3D perpVect, rotnAxis, nbrPerp;
  RDGeom::Point3D nbr1Vect, nbr2Vect, nbr3Vect;
  RDGeom::Transform3D tform;
  RDGeom::Point3D heavyPos, hydPos;

  const Atom *nbr1 = nullptr, *nbr2 = nullptr, *nbr3 = nullptr;
  const Bond *nbrBond;
  ROMol::ADJ_ITER nbrIdx, endNbrs;

  switch (heavyAtom->getDegree()) {
    case 1:
      // --------------------------------------------------------------------------
      //   No other atoms present:
      // --------------------------------------------------------------------------
      // loop over the conformations and set the coordinates
      for (auto cfi = mol->beginConformers(); cfi != mol->endConformers();
           cfi++) {
        if ((*cfi)->is3D()) {
          dirVect.z = 1;
        } else {
          dirVect.x = 1;
        }
        heavyPos = (*cfi)->getAtomPos(heavyIdx);
        hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
        (*cfi)->setAtomPos(hydIdx, hydPos);
      }
      break;

    case 2:
      // --------------------------------------------------------------------------
      //  One other neighbor:
      // --------------------------------------------------------------------------
      nbr1 = getAtomNeighborNot(mol, heavyAtom, hydAtom);
      for (auto cfi = mol->beginConformers(); cfi != mol->endConformers();
           ++cfi) {
        heavyPos = (*cfi)->getAtomPos(heavyIdx);
        RDGeom::Point3D nbr1Pos = (*cfi)->getAtomPos(nbr1->getIdx());
        // get a normalized vector pointing away from the neighbor:
        nbr1Vect = nbr1Pos - heavyPos;
        if (fabs(nbr1Vect.lengthSq()) < 1e-4) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(hydIdx, heavyPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr1Vect *= -1;

        // ok, nbr1Vect points away from the other atom, figure out where
        // this H goes:
        switch (heavyAtom->getHybridization()) {
          case Atom::SP3:
            // get a perpendicular to nbr1Vect:
            if ((*cfi)->is3D()) {
              perpVect = nbr1Vect.getPerpendicular();
            } else {
              perpVect.z = 1.0;
            }
            // and move off it:
            tform.SetRotation((180 - 109.471) * M_PI / 180., perpVect);
            dirVect = tform * nbr1Vect;
            hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP2:
            // default 3D position is to just take an arbitrary perpendicular
            // for 2D we take the normal to the xy plane
            if ((*cfi)->is3D()) {
              perpVect = nbr1Vect.getPerpendicular();
            } else {
              perpVect.z = 1.0;
            }
            if (nbr1->getDegree() > 1) {
              // can we use the neighboring atom to establish a perpendicular?
              nbrBond = mol->getBondBetweenAtoms(heavyIdx, nbr1->getIdx());
              if (nbrBond->getIsAromatic() ||
                  nbrBond->getBondType() == Bond::DOUBLE ||
                  nbrBond->getIsConjugated()) {
                nbr2 = getAtomNeighborNot(mol, nbr1, heavyAtom);
                nbr2Vect =
                    nbr1Pos.directionVector((*cfi)->getAtomPos(nbr2->getIdx()));
                perpVect = nbr2Vect.crossProduct(nbr1Vect);
              }
            }
            perpVect.normalize();
            // rotate the nbr1Vect 60 degrees about perpVect and we're done:
            tform.SetRotation(60. * M_PI / 180., perpVect);
            dirVect = tform * nbr1Vect;
            hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP:
            // just lay the H along the vector:
            dirVect = nbr1Vect;
            hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the vector:
            dirVect = nbr1Vect;
            hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(hydIdx, hydPos);
        }
      }
      break;
    case 3:
      // --------------------------------------------------------------------------
      // Two other neighbors:
      // --------------------------------------------------------------------------
      boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(heavyAtom);
      while (nbrIdx != endNbrs) {
        if (*nbrIdx != hydIdx) {
          if (!nbr1) {
            nbr1 = mol->getAtomWithIdx(*nbrIdx);
          } else {
            nbr2 = mol->getAtomWithIdx(*nbrIdx);
          }
        }
        ++nbrIdx;
      }
      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      for (auto cfi = mol->beginConformers(); cfi != mol->endConformers();
           ++cfi) {
        // start along the average of the two vectors:
        heavyPos = (*cfi)->getAtomPos(heavyIdx);
        nbr1Vect = heavyPos - (*cfi)->getAtomPos(nbr1->getIdx());
        nbr2Vect = heavyPos - (*cfi)->getAtomPos(nbr2->getIdx());
        if (fabs(nbr1Vect.lengthSq()) < 1e-4 ||
            fabs(nbr2Vect.lengthSq()) < 1e-4) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(hydIdx, heavyPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr2Vect.normalize();
        dirVect = nbr1Vect + nbr2Vect;

        dirVect.normalize();
        if ((*cfi)->is3D()) {
          switch (heavyAtom->getHybridization()) {
            case Atom::SP3:
              // get the perpendicular to the neighbors:
              nbrPerp = nbr1Vect.crossProduct(nbr2Vect);
              // and the perpendicular to that:
              rotnAxis = nbrPerp.crossProduct(dirVect);
              // and then rotate about that:
              rotnAxis.normalize();
              tform.SetRotation((109.471 / 2) * M_PI / 180., rotnAxis);
              dirVect = tform * dirVect;
              hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(hydIdx, hydPos);
              break;
            case Atom::SP2:
              // don't need to do anything here, the H atom goes right on the
              // direction vector
              hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(hydIdx, hydPos);
              break;
            default:
              // FIX: handle other hybridizations
              // for now, just lay the H along the neighbor vector;
              hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(hydIdx, hydPos);
              break;
          }
        } else {
          // don't need to do anything here, the H atom goes right on the
          // direction vector
          hydPos = heavyPos + dirVect;
          (*cfi)->setAtomPos(hydIdx, hydPos);
        }
      }
      break;
    case 4:
      // --------------------------------------------------------------------------
      // Three other neighbors:
      // --------------------------------------------------------------------------
      boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(heavyAtom);

      if (heavyAtom->hasProp(common_properties::_CIPCode)) {
        // if the central atom is chiral, we'll order the neighbors
        // by CIP rank:
        std::vector<std::pair<unsigned int, int>> nbrs;
        while (nbrIdx != endNbrs) {
          if (*nbrIdx != hydIdx) {
            const Atom *tAtom = mol->getAtomWithIdx(*nbrIdx);
            unsigned int cip = 0;
            tAtom->getPropIfPresent<unsigned int>(common_properties::_CIPRank,
                                                  cip);
            nbrs.emplace_back(cip, rdcast<int>(*nbrIdx));
          }
          ++nbrIdx;
        }
        std::sort(nbrs.begin(), nbrs.end());
        nbr1 = mol->getAtomWithIdx(nbrs[0].second);
        nbr2 = mol->getAtomWithIdx(nbrs[1].second);
        nbr3 = mol->getAtomWithIdx(nbrs[2].second);
      } else {
        // central atom isn't chiral, so the neighbor ordering isn't important:
        while (nbrIdx != endNbrs) {
          if (*nbrIdx != hydIdx) {
            if (!nbr1) {
              nbr1 = mol->getAtomWithIdx(*nbrIdx);
            } else if (!nbr2) {
              nbr2 = mol->getAtomWithIdx(*nbrIdx);
            } else {
              nbr3 = mol->getAtomWithIdx(*nbrIdx);
            }
          }
          ++nbrIdx;
        }
      }
      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      TEST_ASSERT(nbr3);
      for (auto cfi = mol->beginConformers(); cfi != mol->endConformers();
           ++cfi) {
        // use the average of the three vectors:
        heavyPos = (*cfi)->getAtomPos(heavyIdx);
        nbr1Vect = heavyPos - (*cfi)->getAtomPos(nbr1->getIdx());
        nbr2Vect = heavyPos - (*cfi)->getAtomPos(nbr2->getIdx());
        nbr3Vect = heavyPos - (*cfi)->getAtomPos(nbr3->getIdx());
        if (fabs(nbr1Vect.lengthSq()) < 1e-4 ||
            fabs(nbr2Vect.lengthSq()) < 1e-4 ||
            fabs(nbr3Vect.lengthSq()) < 1e-4) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(hydIdx, heavyPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr2Vect.normalize();
        nbr3Vect.normalize();

        // if three neighboring atoms are more or less planar, this
        // is going to be in a quasi-random (but almost definitely bad)
        // direction...
        // correct for this (issue 2951221):
        if ((*cfi)->is3D()) {
          if (fabs(nbr3Vect.dotProduct(nbr1Vect.crossProduct(nbr2Vect))) <
              0.1) {
            // compute the normal:
            dirVect = nbr1Vect.crossProduct(nbr2Vect);
            std::string cipCode;
            if (heavyAtom->getPropIfPresent(common_properties::_CIPCode,
                                            cipCode)) {
              // the heavy atom is a chiral center, make sure
              // that we went go the right direction to preserve
              // its chirality. We use the chiral volume for this:
              RDGeom::Point3D v1 = dirVect - nbr3Vect;
              RDGeom::Point3D v2 = nbr1Vect - nbr3Vect;
              RDGeom::Point3D v3 = nbr2Vect - nbr3Vect;
              double vol = v1.dotProduct(v2.crossProduct(v3));
              // FIX: this is almost certainly wrong and should use the chiral
              // tag
              if ((cipCode == "S" && vol < 0) || (cipCode == "R" && vol > 0)) {
                dirVect *= -1;
              }
            }
          } else {
            dirVect = nbr1Vect + nbr2Vect + nbr3Vect;
          }
        } else {
          // we're in flatland
          // this was github #908
          // We're in a 2D conformation, put the H between the two neighbors
          // that have the widest angle between them:
          double minDot = nbr1Vect.dotProduct(nbr2Vect);
          dirVect = nbr1Vect + nbr2Vect;
          if (nbr2Vect.dotProduct(nbr3Vect) < minDot) {
            minDot = nbr2Vect.dotProduct(nbr3Vect);
            dirVect = nbr2Vect + nbr3Vect;
          }
          if (nbr1Vect.dotProduct(nbr3Vect) < minDot) {
            minDot = nbr1Vect.dotProduct(nbr3Vect);
            dirVect = nbr1Vect + nbr3Vect;
          }
          dirVect *= -1;
        }
        dirVect.normalize();
        hydPos = heavyPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
        (*cfi)->setAtomPos(hydIdx, hydPos);
      }
      break;
    default:
      // --------------------------------------------------------------------------
      // FIX: figure out what to do here
      // --------------------------------------------------------------------------
      hydPos = heavyPos + dirVect * bondLength;
      for (auto cfi = mol->beginConformers(); cfi != mol->endConformers();
           ++cfi) {
        (*cfi)->setAtomPos(hydIdx, hydPos);
      }
      break;
  }
}

void AssignHsResidueInfo(RWMol &mol) {
  int max_serial = 0;
  unsigned int stopIdx = mol.getNumAtoms();
  for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    auto *info =
        (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
    if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE &&
        info->getSerialNumber() > max_serial) {
      max_serial = info->getSerialNumber();
    }
  }

  AtomPDBResidueInfo *current_info = nullptr;
  int current_h_id = 0;
  for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    Atom *newAt = mol.getAtomWithIdx(aidx);
    auto *info = (AtomPDBResidueInfo *)(newAt->getMonomerInfo());
    if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
      ROMol::ADJ_ITER begin, end;
      boost::tie(begin, end) = mol.getAtomNeighbors(newAt);
      while (begin != end) {
        if (mol.getAtomWithIdx(*begin)->getAtomicNum() == 1) {
          // Make all Hs unique - increment id even for existing
          ++current_h_id;
          // skip if hydrogen already has PDB info
          auto *h_info = (AtomPDBResidueInfo *)mol.getAtomWithIdx(*begin)
                             ->getMonomerInfo();
          if (h_info &&
              h_info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
            continue;
          }
          // the hydrogens have unique names on residue basis (H1, H2, ...)
          if (!current_info ||
              current_info->getResidueNumber() != info->getResidueNumber() ||
              current_info->getChainId() != info->getChainId()) {
            current_h_id = 1;
            current_info = info;
          }
          std::string h_label = std::to_string(current_h_id);
          if (h_label.length() > 3) {
            h_label = h_label.substr(h_label.length() - 3, 3);
          }
          while (h_label.length() < 3) {
            h_label = h_label + " ";
          }
          h_label = "H" + h_label;
          // wrap around id to '3H12'
          h_label = h_label.substr(3, 1) + h_label.substr(0, 3);
          AtomPDBResidueInfo *newInfo = new AtomPDBResidueInfo(
              h_label, max_serial, "", info->getResidueName(),
              info->getResidueNumber(), info->getChainId(), "", 1.0, 0.0,
              info->getIsHeteroAtom());
          mol.getAtomWithIdx(*begin)->setMonomerInfo(newInfo);

          ++max_serial;
        }
        ++begin;
      }
    }
  }
}

std::string isoHsToString(const std::vector<unsigned int> &isoHs) {
  std::stringstream ss;
  std::copy(isoHs.begin(), isoHs.end(),
            std::ostream_iterator<unsigned int>(ss, " "));
  std::string res(std::move(ss.str()));
  boost::trim(res);
  return res;
}

std::map<unsigned int, std::vector<unsigned int>> getIsoMap(const ROMol &mol) {
  std::map<unsigned int, std::vector<unsigned int>> isoMap;
  for (auto atom : mol.atoms()) {
    if (atom->hasProp(common_properties::_isotopicHs)) {
      atom->clearProp(common_properties::_isotopicHs);
    }
  }
  for (auto bond : mol.bonds()) {
    auto ba = bond->getBeginAtom();
    auto ea = bond->getEndAtom();
    int ha = -1;
    unsigned int iso;
    if (ba->getAtomicNum() == 1 && ba->getIsotope() &&
        ea->getAtomicNum() != 1) {
      ha = ea->getIdx();
      iso = ba->getIsotope();
    } else if (ea->getAtomicNum() == 1 && ea->getIsotope() &&
               ba->getAtomicNum() != 1) {
      ha = ba->getIdx();
      iso = ea->getIsotope();
    }
    if (ha == -1) {
      continue;
    }
    auto &v = isoMap[ha];
    v.push_back(iso);
  }
  return isoMap;
}

}  // end of unnamed namespace

namespace MolOps {

void addHs(RWMol &mol, bool explicitOnly, bool addCoords,
           const UINT_VECT *onlyOnAtoms, bool addResidueInfo) {
  // when we hit each atom, clear its computed properties
  // NOTE: it is essential that we not clear the ring info in the
  // molecule's computed properties.  We don't want to have to
  // regenerate that.  This caused Issue210 and Issue212:
  mol.clearComputedProps(false);

  // precompute the number of hydrogens we are going to add so that we can
  // pre-allocate the necessary space on the conformations of the molecule
  // for their coordinates
  unsigned int numAddHyds = 0;
  for (auto at : mol.atoms()) {
    if (!onlyOnAtoms || std::find(onlyOnAtoms->begin(), onlyOnAtoms->end(),
                                  at->getIdx()) != onlyOnAtoms->end()) {
      numAddHyds += at->getNumExplicitHs();
      if (!explicitOnly) {
        numAddHyds += at->getNumImplicitHs();
      }
    }
  }
  unsigned int nSize = mol.getNumAtoms() + numAddHyds;

  // loop over the conformations of the molecule and allocate new space
  // for the H locations (need to do this even if we aren't adding coords so
  // that the conformers have the correct number of atoms).
  for (auto cfi = mol.beginConformers(); cfi != mol.endConformers(); ++cfi) {
    (*cfi)->reserve(nSize);
  }

  unsigned int stopIdx = mol.getNumAtoms();
  for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    if (onlyOnAtoms && std::find(onlyOnAtoms->begin(), onlyOnAtoms->end(),
                                 aidx) == onlyOnAtoms->end()) {
      continue;
    }

    Atom *newAt = mol.getAtomWithIdx(aidx);

    std::vector<unsigned int> isoHs;
    std::string isotopicHsProp;
    if (newAt->getPropIfPresent(common_properties::_isotopicHs,
                                isotopicHsProp)) {
      newAt->clearProp(common_properties::_isotopicHs);
      // be lenient on input, even if we write only space-separated
      // strings of indices
      boost::trim_if(isotopicHsProp, boost::is_any_of(" \t\r\n,()[]{}"));
      boost::tokenizer<> tokens(isotopicHsProp);
      std::transform(tokens.begin(), tokens.end(), std::back_inserter(isoHs),
                     [](const std::string &t) {
                       return boost::lexical_cast<unsigned int>(t);
                     });
    }
    std::vector<unsigned int>::const_iterator isoH = isoHs.begin();
    unsigned int newIdx;
    newAt->clearComputedProps();
    // always convert explicit Hs
    unsigned int onumexpl = newAt->getNumExplicitHs();
    for (unsigned int i = 0; i < onumexpl; i++) {
      newIdx = mol.addAtom(new Atom(1), false, true);
      mol.addBond(aidx, newIdx, Bond::SINGLE);
      auto hAtom = mol.getAtomWithIdx(newIdx);
      hAtom->updatePropertyCache();
      if (addCoords) {
        setHydrogenCoords(&mol, newIdx, aidx);
      }
      if (isoH != isoHs.end()) {
        hAtom->setIsotope(*isoH);
        ++isoH;
      }
    }
    // clear the local property
    newAt->setNumExplicitHs(0);

    if (!explicitOnly) {
      // take care of implicits
      for (unsigned int i = 0; i < mol.getAtomWithIdx(aidx)->getNumImplicitHs();
           i++) {
        newIdx = mol.addAtom(new Atom(1), false, true);
        mol.addBond(aidx, newIdx, Bond::SINGLE);
        // set the isImplicit label so that we can strip these back
        // off later if need be.
        auto hAtom = mol.getAtomWithIdx(newIdx);
        hAtom->setProp(common_properties::isImplicit, 1);
        hAtom->updatePropertyCache();
        if (addCoords) {
          setHydrogenCoords(&mol, newIdx, aidx);
        }
        if (isoH != isoHs.end()) {
          hAtom->setIsotope(*isoH);
          ++isoH;
        }
      }
      // be very clear about implicits not being allowed in this
      // representation
      newAt->setProp(common_properties::origNoImplicit, newAt->getNoImplicit(),
                     true);
      newAt->setNoImplicit(true);
    }
    // update the atom's derived properties (valence count, etc.)
    // no sense in being strict here (was github #2782)
    newAt->updatePropertyCache(false);
    if (isoH != isoHs.end()) {
      BOOST_LOG(rdWarningLog) << "extra H isotope information found on atom "
                              << newAt->getIdx() << std::endl;
    }
  }
  // take care of AtomPDBResidueInfo for Hs if root atom has it
  if (addResidueInfo) {
    AssignHsResidueInfo(mol);
  }
}

ROMol *addHs(const ROMol &mol, bool explicitOnly, bool addCoords,
             const UINT_VECT *onlyOnAtoms, bool addResidueInfo) {
  auto *res = new RWMol(mol);
  addHs(*res, explicitOnly, addCoords, onlyOnAtoms, addResidueInfo);
  return static_cast<ROMol *>(res);
};

namespace {
// returns whether or not an adjustment was made, in case we want that info
bool adjustStereoAtomsIfRequired(RWMol &mol, const Atom *atom,
                                 const Atom *heavyAtom) {
  PRECONDITION(atom != nullptr, "bad atom");
  PRECONDITION(heavyAtom != nullptr, "bad heavy atom");
  // nothing we can do if the degree is only 2 (and we should have covered
  // that earlier anyway)
  if (heavyAtom->getDegree() == 2) {
    return false;
  }
  const auto &cbnd =
      mol.getBondBetweenAtoms(atom->getIdx(), heavyAtom->getIdx());
  if (!cbnd) {
    return false;
  }
  for (const auto &nbri :
       boost::make_iterator_range(mol.getAtomBonds(heavyAtom))) {
    Bond *bnd = mol[nbri];
    if (bnd->getBondType() == Bond::DOUBLE &&
        bnd->getStereo() > Bond::STEREOANY) {
      auto sAtomIt = std::find(bnd->getStereoAtoms().begin(),
                               bnd->getStereoAtoms().end(), atom->getIdx());
      if (sAtomIt != bnd->getStereoAtoms().end()) {
        // sAtomIt points to the position of this atom's index in the list.
        // find the index of another atom attached to the heavy atom and
        // use it to update sAtomIt
        unsigned int dblNbrIdx = bnd->getOtherAtomIdx(heavyAtom->getIdx());
        for (const auto &nbri :
             boost::make_iterator_range(mol.getAtomNeighbors(heavyAtom))) {
          const auto &nbr = mol[nbri];
          if (nbr->getIdx() == dblNbrIdx || nbr->getIdx() == atom->getIdx()) {
            continue;
          }
          *sAtomIt = nbr->getIdx();
          bool madeAdjustment = true;
          switch (bnd->getStereo()) {
            case Bond::STEREOCIS:
              bnd->setStereo(Bond::STEREOTRANS);
              break;
            case Bond::STEREOTRANS:
              bnd->setStereo(Bond::STEREOCIS);
              break;
            default:
              // I think we shouldn't need to do anything with E and Z...
              madeAdjustment = false;
              break;
          }
          return madeAdjustment;
        }
      }
    }
  }
  return false;
}

void molRemoveH(RWMol &mol, unsigned int idx, bool updateExplicitCount) {
  auto atom = mol.getAtomWithIdx(idx);
  PRECONDITION(atom->getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");
  for (const auto &nbri : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    const Bond *bond = mol[nbri];
    Atom *heavyAtom = bond->getOtherAtom(atom);
    int heavyAtomNum = heavyAtom->getAtomicNum();

    // we'll update the neighbor's explicit H count if we were told to
    // *or* if the neighbor is chiral, in which case the H is needed
    // in order to complete the coordination
    // *or* if the neighbor has the noImplicit flag set:
    if (updateExplicitCount || heavyAtom->getNoImplicit() ||
        heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
    } else {
      // this is a special case related to Issue 228 and the
      // "disappearing Hydrogen" problem discussed in MolOps::adjustHs
      //
      // If we remove a hydrogen from an aromatic N or P, or if
      // the heavy atom it is connected to is not in its default
      // valence state, we need to be *sure* to increment the
      // explicit count, even if the H itself isn't marked as explicit
      const INT_VECT &defaultVs =
          PeriodicTable::getTable()->getValenceList(heavyAtomNum);
      if (((heavyAtomNum == 7 || heavyAtomNum == 15) &&
           heavyAtom->getIsAromatic()) ||
          (std::find(defaultVs.begin() + 1, defaultVs.end(),
                     heavyAtom->getTotalValence()) != defaultVs.end())) {
        heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
      }
    }

    // One other consequence of removing the H from the graph is
    // that we may change the ordering of the bonds about a
    // chiral center.  This may change the chiral label at that
    // atom.  We deal with that by explicitly checking here:
    if (heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      INT_LIST neighborIndices;
      for (const auto &nbri :
           boost::make_iterator_range(mol.getAtomBonds(heavyAtom))) {
        Bond *nbnd = mol[nbri];
        if (nbnd->getIdx() != bond->getIdx()) {
          neighborIndices.push_back(nbnd->getIdx());
        }
      }
      neighborIndices.push_back(bond->getIdx());

      int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
      // std::cerr << "H: "<<atom->getIdx()<<" hvy:
      // "<<heavyAtom->getIdx()<<" swaps: " << nSwaps<<std::endl;
      if (nSwaps % 2) {
        heavyAtom->invertChirality();
      }
    }

    // if it's a wavy bond, then we need to
    // mark the beginning atom with the _UnknownStereo tag.
    // so that we know later that something was affecting its
    // stereochem
    if (bond->getBondDir() == Bond::UNKNOWN &&
        bond->getBeginAtomIdx() == heavyAtom->getIdx()) {
      heavyAtom->setProp(common_properties::_UnknownStereo, 1);
    } else if (bond->getBondDir() == Bond::ENDDOWNRIGHT ||
               bond->getBondDir() == Bond::ENDUPRIGHT) {
      // if the direction is set on this bond and the atom it's connected to
      // has no other single bonds with directions set, then we need to set
      // direction on one of the other neighbors in order to avoid double
      // bond stereochemistry possibly being lost. This was github #754
      bool foundADir = false;
      Bond *oBond = nullptr;
      for (const auto &nbri :
           boost::make_iterator_range(mol.getAtomBonds(heavyAtom))) {
        Bond *nbnd = mol[nbri];
        if (nbnd->getIdx() != bond->getIdx() &&
            nbnd->getBondType() == Bond::SINGLE) {
          if (nbnd->getBondDir() == Bond::NONE) {
            oBond = nbnd;
          } else {
            foundADir = true;
          }
        }
      }
      if (!foundADir && oBond != nullptr) {
        bool flipIt = (oBond->getBeginAtom() == heavyAtom) &&
                      (bond->getBeginAtom() == heavyAtom);
        if (flipIt) {
          oBond->setBondDir(bond->getBondDir() == Bond::ENDDOWNRIGHT
                                ? Bond::ENDUPRIGHT
                                : Bond::ENDDOWNRIGHT);
        } else {
          oBond->setBondDir(bond->getBondDir());
        }
      }
      // if this atom is one of the stereoatoms for a double bond we need
      // to switch the stereo atom on this end to be the other neighbor
      // This was part of github #1810
      adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
    } else {
      // if this atom is one of the stereoatoms for a double bond we need
      // to switch the stereo atom on this end to be the other neighbor
      // This was part of github #1810
      adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
    }
  }
  mol.removeAtom(atom);
}

}  // end of anonymous namespace

void removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize) {
  if (ps.removeAndTrackIsotopes) {
    // if there are any non-isotopic Hs remove them first
    // to make sure chirality is preserved
    bool needRemoveHs = false;
    for (auto atom : mol.atoms()) {
      if (atom->getAtomicNum() == 1 && atom->getIsotope() == 0) {
        needRemoveHs = true;
        break;
      }
    }
    if (needRemoveHs) {
      RemoveHsParameters psCopy(ps);
      psCopy.removeAndTrackIsotopes = false;
      psCopy.removeIsotopes = false;
      removeHs(mol, psCopy, false);
    }
  }
  for (auto atom : mol.atoms()) {
    atom->updatePropertyCache(false);
  }
  if (ps.removeAndTrackIsotopes) {
    for (const auto &pair : getIsoMap(mol)) {
      mol.getAtomWithIdx(pair.first)
          ->setProp(common_properties::_isotopicHs, isoHsToString(pair.second));
    }
  }
  boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
  for (auto atom : mol.atoms()) {
    if (atom->getAtomicNum() != 1) {
      continue;
    }
    if (!ps.removeWithQuery && atom->hasQuery()) {
      continue;
    }
    if (!ps.removeDegreeZero && !atom->getDegree()) {
      if (ps.showWarnings) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: not removing hydrogen atom without neighbors"
            << std::endl;
      }
      continue;
    }
    if (!ps.removeHigherDegrees && atom->getDegree() > 1) {
      continue;
    }
    if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes &&
        atom->getIsotope()) {
      continue;
    }
    if (!ps.removeNonimplicit &&
        !atom->hasProp(common_properties::isImplicit)) {
      continue;
    }
    if (!ps.removeMapped && atom->getAtomMapNum()) {
      continue;
    }
    if (!ps.removeInSGroups) {
      bool skipIt = false;
      for (const auto &sg : getSubstanceGroups(mol)) {
        if (sg.includesAtom(atom->getIdx())) {
          skipIt = true;
          break;
        }
      }
      if (skipIt) {
        continue;
      }
    }
    if (!ps.removeHydrides && atom->getFormalCharge() == -1) {
      continue;
    }
    bool removeIt = true;
    if (atom->getDegree() &&
        (!ps.removeDummyNeighbors || !ps.removeDefiningBondStereo ||
         !ps.removeOnlyHNeighbors)) {
      bool onlyHNeighbors = true;
      ROMol::ADJ_ITER begin, end;
      boost::tie(begin, end) = mol.getAtomNeighbors(atom);
      while (begin != end && removeIt) {
        auto nbr = mol.getAtomWithIdx(*begin);
        // is it a dummy?
        if (!ps.removeDummyNeighbors && nbr->getAtomicNum() < 1) {
          removeIt = false;
          if (ps.showWarnings) {
            BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
                                       "with dummy atom neighbors"
                                    << std::endl;
          }
        }
        if (!ps.removeOnlyHNeighbors && nbr->getAtomicNum() != 1) {
          onlyHNeighbors = false;
        }
        if (!ps.removeWithWedgedBond) {
          const auto bnd =
              mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
          if (bnd->getBondDir() == Bond::BEGINDASH ||
              bnd->getBondDir() == Bond::BEGINWEDGE) {
            removeIt = false;
            if (ps.showWarnings) {
              BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
                                         "with wedged bond"
                                      << std::endl;
            }
          }
        }
        // Check to see if the neighbor has a double bond and we're the only
        // neighbor at this end.  This was part of github #1810
        if (!ps.removeDefiningBondStereo && nbr->getDegree() == 2) {
          for (const auto &nbri :
               boost::make_iterator_range(mol.getAtomBonds(nbr))) {
            const Bond *bnd = mol[nbri];
            if (bnd->getBondType() == Bond::DOUBLE &&
                (bnd->getStereo() > Bond::STEREOANY ||
                 mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
                         ->getBondDir() > Bond::NONE)) {
              removeIt = false;
              break;
            }
          }
        }
        ++begin;
      }
      if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
        removeIt = false;
      }
    }
    if (removeIt) {
      atomsToRemove.set(atom->getIdx());
    }
  }  // end of the loop over atoms
  // now that we know which atoms need to be removed, go ahead and remove them
  for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    if (atomsToRemove[idx]) {
      molRemoveH(mol, idx, ps.updateExplicitCount);
    }
  }
  //
  //  If we didn't only remove implicit Hs, which are guaranteed to
  //  be the highest numbered atoms, we may have altered atom indices.
  //  This can screw up derived properties (such as ring members), so
  //  do some checks:
  //
  if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    sanitizeMol(mol);
  }
};
ROMol *removeHs(const ROMol &mol, const RemoveHsParameters &ps, bool sanitize) {
  auto *res = new RWMol(mol);
  try {
    removeHs(*res, ps, sanitize);
  } catch (const MolSanitizeException &) {
    delete res;
    throw;
  }
  return static_cast<ROMol *>(res);
}
void removeHs(RWMol &mol, bool implicitOnly, bool updateExplicitCount,
              bool sanitize) {
  RemoveHsParameters ps;
  ps.removeNonimplicit = !implicitOnly;
  ps.updateExplicitCount = updateExplicitCount;
  removeHs(mol, ps, sanitize);
};
ROMol *removeHs(const ROMol &mol, bool implicitOnly, bool updateExplicitCount,
                bool sanitize) {
  auto *res = new RWMol(mol);
  try {
    removeHs(*res, implicitOnly, updateExplicitCount, sanitize);
  } catch (const MolSanitizeException &) {
    delete res;
    throw;
  }
  return static_cast<ROMol *>(res);
}

void removeAllHs(RWMol &mol, bool sanitize) {
  RemoveHsParameters ps;
  ps.removeDegreeZero = true;
  ps.removeHigherDegrees = true;
  ps.removeOnlyHNeighbors = true;
  ps.removeIsotopes = true;
  ps.removeDummyNeighbors = true;
  ps.removeDefiningBondStereo = true;
  ps.removeWithWedgedBond = true;
  ps.removeWithQuery = true;
  ps.removeNonimplicit = true;
  ps.removeInSGroups = true;
  ps.showWarnings = false;
  ps.removeHydrides = true;
  removeHs(mol, ps, sanitize);
};
ROMol *removeAllHs(const ROMol &mol, bool sanitize) {
  auto *res = new RWMol(mol);
  try {
    removeAllHs(*res, sanitize);
  } catch (const MolSanitizeException &) {
    delete res;
    throw;
  }
  return static_cast<ROMol *>(res);
}

namespace {
bool isQueryH(const Atom *atom) {
  PRECONDITION(atom, "bogus atom");
  if (atom->getAtomicNum() == 1) {
    // the simple case: the atom is flagged as being an H and
    // has no query
    if (!atom->hasQuery() ||
        (!atom->getQuery()->getNegation() &&
         atom->getQuery()->getDescription() == "AtomAtomicNum")) {
      return true;
    }
  }

  if (atom->getDegree() != 1) {
    // only degree 1
    return false;
  }

  if (atom->hasQuery() && atom->getQuery()->getNegation()) {
    // we will not merge negated queries
    return false;
  }

  bool hasHQuery = false, hasOr = false;
  if (atom->hasQuery()) {
    if (atom->getQuery()->getDescription() == "AtomOr") {
      hasOr = true;
    }
    std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
        atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
    // the logic gets too complicated if there's an OR in the children, so
    // just punt on those (with a warning)
    while (!(hasHQuery && hasOr) && childStack.size()) {
      QueryAtom::QUERYATOM_QUERY::CHILD_TYPE query = childStack.front();
      childStack.pop_front();
      if (query->getDescription() == "AtomOr") {
        hasOr = true;
      } else if (query->getDescription() == "AtomAtomicNum") {
        if (static_cast<ATOM_EQUALS_QUERY *>(query.get())->getVal() == 1 &&
            !query->getNegation()) {
          hasHQuery = true;
        }
      } else {
        QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI child1;
        for (child1 = query->beginChildren(); child1 != query->endChildren();
             ++child1) {
          childStack.push_back(*child1);
        }
      }
    }
    // std::cerr<<"   !!!1 "<<atom->getIdx()<<" "<<hasHQuery<<"
    // "<<hasOr<<std::endl;
    if (hasHQuery && hasOr) {
      BOOST_LOG(rdWarningLog) << "WARNING: merging explicit H queries involved "
                                 "in ORs is not supported. This query will not "
                                 "be merged"
                              << std::endl;
      return false;
    }
  }
  return hasHQuery;
}
}  // namespace

//
//  This routine removes explicit hydrogens (and bonds to them) from
//  the molecular graph and adds them as queries to the heavy atoms
//  to which they are bound.  If the heavy atoms (or atom queries)
//  already have hydrogen-count queries, they will be updated.
//
//  NOTE:
//   - Hydrogens which aren't connected to a heavy atom will not be
//     removed.  This prevents molecules like "[H][H]" from having
//     all atoms removed.
//
//   - By default all hydrogens are removed, however if
//     merge_unmapped_only is true, any hydrogen participating
//     in an atom map will be retained
void mergeQueryHs(RWMol &mol, bool mergeUnmappedOnly) {
  std::vector<unsigned int> atomsToRemove;

  boost::dynamic_bitset<> hatoms(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    hatoms[i] = isQueryH(mol.getAtomWithIdx(i));
  }
  unsigned int currIdx = 0, stopIdx = mol.getNumAtoms();
  while (currIdx < stopIdx) {
    Atom *atom = mol.getAtomWithIdx(currIdx);
    if (!hatoms[currIdx]) {
      unsigned int numHsToRemove = 0;
      ROMol::ADJ_ITER begin, end;
      boost::tie(begin, end) = mol.getAtomNeighbors(atom);

      while (begin != end) {
        if (hatoms[*begin]) {
          Atom &bgn = *mol.getAtomWithIdx(*begin);
          if (!mergeUnmappedOnly ||
              !bgn.hasProp(common_properties::molAtomMapNumber)) {
            atomsToRemove.push_back(rdcast<unsigned int>(*begin));
            ++numHsToRemove;
          }
        }
        ++begin;
      }
      if (numHsToRemove) {
        //
        //  We have H neighbors:
        //   Add the appropriate queries to compensate for their removal.
        //
        //  Examples:
        //    C[H] -> [C;!H0]
        //    C([H])[H] -> [C;!H0;!H1]
        //
        //  It would be more efficient to do this using range queries like:
        //    C([H])[H] -> [C;H{2-}]
        //  but that would produce non-standard SMARTS without the user
        //  having started with a non-standard SMARTS.
        //
        if (!atom->hasQuery()) {
          // it wasn't a query atom, we need to replace it so that we can add
          // a query:
          ATOM_EQUALS_QUERY *tmp = makeAtomNumQuery(atom->getAtomicNum());
          auto *newAt = new QueryAtom;
          newAt->setQuery(tmp);
          newAt->updateProps(*atom);
          mol.replaceAtom(atom->getIdx(), newAt);
          delete newAt;
          atom = mol.getAtomWithIdx(currIdx);
        }
        for (unsigned int i = 0; i < numHsToRemove; ++i) {
          ATOM_EQUALS_QUERY *tmp = makeAtomHCountQuery(i);
          tmp->setNegation(true);
          atom->expandQuery(tmp);
        }
      }  // end of numHsToRemove test

      // recurse if needed (was github isusue 544)
      if (atom->hasQuery()) {
        // std::cerr<<"  q: "<<atom->getQuery()->getDescription()<<std::endl;
        if (atom->getQuery()->getDescription() == "RecursiveStructure") {
          auto *rqm = static_cast<RWMol *>(const_cast<ROMol *>(
              static_cast<RecursiveStructureQuery *>(atom->getQuery())
                  ->getQueryMol()));
          mergeQueryHs(*rqm, mergeUnmappedOnly);
        }

        // FIX: shouldn't be repeating this code here
        std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
            atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
        while (childStack.size()) {
          QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
          childStack.pop_front();
          // std::cerr<<"      child: "<<qry->getDescription()<<std::endl;
          if (qry->getDescription() == "RecursiveStructure") {
            // std::cerr<<"    recurse"<<std::endl;
            auto *rqm = static_cast<RWMol *>(const_cast<ROMol *>(
                static_cast<RecursiveStructureQuery *>(qry.get())
                    ->getQueryMol()));
            mergeQueryHs(*rqm, mergeUnmappedOnly);
            // std::cerr<<"    back"<<std::endl;
          } else if (qry->beginChildren() != qry->endChildren()) {
            childStack.insert(childStack.end(), qry->beginChildren(),
                              qry->endChildren());
          }
        }
      }  // end of recursion loop
    }
    ++currIdx;
  }
  std::sort(atomsToRemove.begin(), atomsToRemove.end());
  for (std::vector<unsigned int>::const_reverse_iterator aiter =
           atomsToRemove.rbegin();
       aiter != atomsToRemove.rend(); ++aiter) {
    Atom *atom = mol.getAtomWithIdx(*aiter);
    mol.removeAtom(atom);
  }
};
ROMol *mergeQueryHs(const ROMol &mol, bool mergeUnmappedOnly) {
  auto *res = new RWMol(mol);
  mergeQueryHs(*res, mergeUnmappedOnly);
  return static_cast<ROMol *>(res);
};

};  // end of namespace MolOps
};  // namespace RDKit
