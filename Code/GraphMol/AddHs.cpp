//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
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
#include "Chirality.h"
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <boost/algorithm/string/classification.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/iterator_range.hpp>

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

bool may_need_extra_H(const ROMol &mol, const Atom *atom) {
  unsigned single_bonds = 0;
  unsigned aromatic_bonds = 0;
  for (auto bond : mol.atomBonds(atom)) {
    if (bond->getBondType() == Bond::SINGLE) {
      ++single_bonds;
    } else if (bond->getBondType() == Bond::AROMATIC) {
      ++aromatic_bonds;
    } else {
      return false;
    }
  }
  return single_bonds == 1 && aromatic_bonds == 2 &&
         atom->getTotalValence() == 3;
}

}  // end of unnamed namespace

namespace MolOps {

namespace {
RDGeom::Point3D pickBisector(const RDGeom::Point3D &nbr1Vect,
                             const RDGeom::Point3D &nbr2Vect,
                             const RDGeom::Point3D &nbr3Vect) {
  auto dirVect = nbr2Vect + nbr3Vect;
  if (dirVect.lengthSq() < 1e-4) {
    // nbr2Vect and nbr3Vect are anti-parallel (was #3854)
    dirVect = nbr2Vect;
    std::swap(dirVect.x, dirVect.y);
    dirVect.x *= -1;
  }
  if (dirVect.dotProduct(nbr1Vect) < 0) {
    dirVect *= -1;
  }

  return dirVect;
}
}  // namespace

void setTerminalAtomCoords(ROMol &mol, unsigned int idx,
                           unsigned int otherIdx) {
  // we will loop over all the coordinates
  PRECONDITION(otherIdx != idx, "degenerate atoms");
  Atom *atom = mol.getAtomWithIdx(idx);
  PRECONDITION(mol.getAtomDegree(atom) == 1, "bad atom degree");
  const Bond *bond = mol.getBondBetweenAtoms(otherIdx, idx);
  PRECONDITION(bond, "no bond between atoms");

  const Atom *otherAtom = mol.getAtomWithIdx(otherIdx);
  double bondLength =
      PeriodicTable::getTable()->getRb0(1) +
      PeriodicTable::getTable()->getRb0(otherAtom->getAtomicNum());

  RDGeom::Point3D dirVect(0, 0, 0);

  RDGeom::Point3D perpVect, rotnAxis, nbrPerp;
  RDGeom::Point3D nbr1Vect, nbr2Vect, nbr3Vect;
  RDGeom::Transform3D tform;
  RDGeom::Point3D otherPos, atomPos;

  const Atom *nbr1 = nullptr, *nbr2 = nullptr, *nbr3 = nullptr;
  const Bond *nbrBond;
  ROMol::ADJ_ITER nbrIdx, endNbrs;

  switch (otherAtom->getDegree()) {
    case 1:
      // --------------------------------------------------------------------------
      //   No other atoms present:
      // --------------------------------------------------------------------------
      // loop over the conformations and set the coordinates
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           cfi++) {
        if ((*cfi)->is3D()) {
          dirVect.z = 1;
        } else {
          dirVect.x = 1;
        }
        otherPos = (*cfi)->getAtomPos(otherIdx);
        atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
        (*cfi)->setAtomPos(idx, atomPos);
      }
      break;

    case 2:
      // --------------------------------------------------------------------------
      //  One other neighbor:
      // --------------------------------------------------------------------------
      nbr1 = getAtomNeighborNot(&mol, otherAtom, atom);
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        otherPos = (*cfi)->getAtomPos(otherIdx);
        RDGeom::Point3D nbr1Pos = (*cfi)->getAtomPos(nbr1->getIdx());
        // get a normalized vector pointing away from the neighbor:
        nbr1Vect = nbr1Pos - otherPos;
        if (nbr1Vect.lengthSq() < 1e-4) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(idx, otherPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr1Vect *= -1;

        // ok, nbr1Vect points away from the other atom, figure out where
        // this H goes:
        switch (otherAtom->getHybridization()) {
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
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
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
              nbrBond = mol.getBondBetweenAtoms(otherIdx, nbr1->getIdx());
              if (nbrBond->getIsAromatic() ||
                  nbrBond->getBondType() == Bond::DOUBLE ||
                  nbrBond->getIsConjugated()) {
                nbr2 = getAtomNeighborNot(&mol, nbr1, otherAtom);
                nbr2Vect =
                    nbr1Pos.directionVector((*cfi)->getAtomPos(nbr2->getIdx()));
                perpVect = nbr2Vect.crossProduct(nbr1Vect);
              }
            }
            perpVect.normalize();
            // rotate the nbr1Vect 60 degrees about perpVect and we're done:
            tform.SetRotation(60. * M_PI / 180., perpVect);
            dirVect = tform * nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
            break;
          case Atom::SP:
            // just lay the H along the vector:
            dirVect = nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the vector:
            dirVect = nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
        }
      }
      break;
    case 3:
      // --------------------------------------------------------------------------
      // Two other neighbors:
      // --------------------------------------------------------------------------
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(otherAtom);
      while (nbrIdx != endNbrs) {
        if (*nbrIdx != idx) {
          if (!nbr1) {
            nbr1 = mol.getAtomWithIdx(*nbrIdx);
          } else {
            nbr2 = mol.getAtomWithIdx(*nbrIdx);
          }
        }
        ++nbrIdx;
      }
      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        // start along the average of the two vectors:
        otherPos = (*cfi)->getAtomPos(otherIdx);
        nbr1Vect = otherPos - (*cfi)->getAtomPos(nbr1->getIdx());
        nbr2Vect = otherPos - (*cfi)->getAtomPos(nbr2->getIdx());
        if (nbr1Vect.lengthSq() < 1e-4 || nbr2Vect.lengthSq() < 1e-4) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(idx, otherPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr2Vect.normalize();
        dirVect = nbr1Vect + nbr2Vect;

        dirVect.normalize();
        if ((*cfi)->is3D()) {
          switch (otherAtom->getHybridization()) {
            case Atom::SP3:
              // get the perpendicular to the neighbors:
              nbrPerp = nbr1Vect.crossProduct(nbr2Vect);
              // and the perpendicular to that:
              rotnAxis = nbrPerp.crossProduct(dirVect);
              // and then rotate about that:
              rotnAxis.normalize();
              tform.SetRotation((109.471 / 2) * M_PI / 180., rotnAxis);
              dirVect = tform * dirVect;
              atomPos =
                  otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(idx, atomPos);
              break;
            case Atom::SP2:
              // don't need to do anything here, the H atom goes right on the
              // direction vector
              atomPos =
                  otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(idx, atomPos);
              break;
            default:
              // FIX: handle other hybridizations
              // for now, just lay the H along the neighbor vector;
              atomPos =
                  otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(idx, atomPos);
              break;
          }
        } else {
          // don't need to do anything here, the H atom goes right on the
          // direction vector
          atomPos = otherPos + dirVect;
          (*cfi)->setAtomPos(idx, atomPos);
        }
      }
      break;
    case 4:
      // --------------------------------------------------------------------------
      // Three other neighbors:
      // --------------------------------------------------------------------------
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(otherAtom);

      // We're using chiral tag for checking chirality, so we just take the
      // initial order
      while (nbrIdx != endNbrs) {
        if (*nbrIdx != idx) {
          if (!nbr1) {
            nbr1 = mol.getAtomWithIdx(*nbrIdx);
          } else if (!nbr2) {
            nbr2 = mol.getAtomWithIdx(*nbrIdx);
          } else {
            nbr3 = mol.getAtomWithIdx(*nbrIdx);
          }
        }
        ++nbrIdx;
      }

      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      TEST_ASSERT(nbr3);

      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        otherPos = (*cfi)->getAtomPos(otherIdx);
        nbr1Vect = otherPos - (*cfi)->getAtomPos(nbr1->getIdx());
        nbr2Vect = otherPos - (*cfi)->getAtomPos(nbr2->getIdx());
        nbr3Vect = otherPos - (*cfi)->getAtomPos(nbr3->getIdx());
        if (nbr1Vect.lengthSq() < 1e-4 || nbr2Vect.lengthSq() < 1e-4 ||
            nbr3Vect.lengthSq() < 1e-4) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(idx, otherPos);
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
            if (otherAtom->getPropIfPresent(common_properties::_CIPCode,
                                            cipCode)) {
              // the heavy atom is a chiral center, make sure
              // that we went go the right direction to preserve
              // its chirality. We use the chiral volume for this:
              RDGeom::Point3D v1 = dirVect - nbr3Vect;
              RDGeom::Point3D v2 = nbr1Vect - nbr3Vect;
              RDGeom::Point3D v3 = nbr2Vect - nbr3Vect;
              double vol = v1.dotProduct(v2.crossProduct(v3));

              if ((otherAtom->getChiralTag() ==
                       Atom::ChiralType::CHI_TETRAHEDRAL_CCW &&
                   vol < 0) ||
                  (otherAtom->getChiralTag() ==
                       Atom::ChiralType::CHI_TETRAHEDRAL_CW &&
                   vol > 0)) {
                dirVect *= -1;
              }
            }
          } else {
            dirVect = nbr1Vect + nbr2Vect + nbr3Vect;
          }
        } else {
          // we're in flatland

          // github #3879 and #908: find the two neighbors with the largest
          // outer angle between them and then place the H to bisect that angle
          // This is recommendation ST-1.1.4 from the 2006 IUPAC "Graphical
          // representation of stereochemical configuration" guideline
          auto angle12 = nbr1Vect.angleTo(nbr2Vect);
          auto angle13 = nbr1Vect.angleTo(nbr3Vect);
          auto angle23 = nbr2Vect.angleTo(nbr3Vect);
          auto accum1 = angle12 + angle13;
          auto accum2 = angle12 + angle23;
          auto accum3 = angle13 + angle23;
          if (accum1 <= accum2 && accum1 <= accum3) {
            dirVect = pickBisector(nbr1Vect, nbr2Vect, nbr3Vect);
          } else if (accum2 <= accum1 && accum2 <= accum3) {
            dirVect = pickBisector(nbr2Vect, nbr1Vect, nbr3Vect);
          } else {
            dirVect = pickBisector(nbr3Vect, nbr1Vect, nbr2Vect);
          }
        }

        dirVect.normalize();
        atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
        (*cfi)->setAtomPos(idx, atomPos);
      }
      break;
    default:
      // --------------------------------------------------------------------------
      // FIX: figure out what to do here
      // --------------------------------------------------------------------------
      atomPos = otherPos + dirVect * bondLength;
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        (*cfi)->setAtomPos(idx, atomPos);
      }
      break;
  }
}

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
    if (newAt->getPropIfPresent(common_properties::_isotopicHs, isoHs)) {
      newAt->clearProp(common_properties::_isotopicHs);
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
        setTerminalAtomCoords(mol, newIdx, aidx);
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
          setTerminalAtomCoords(mol, newIdx, aidx);
        }
        if (isoH != isoHs.end()) {
          hAtom->setIsotope(*isoH);
          ++isoH;
        }
      }
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
  for (const auto bond : mol.atomBonds(atom)) {
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
      if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
            may_need_extra_H(mol, heavyAtom)) &&
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
      for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
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

    // If we are removing a H atom that defines bond stereo (e.g. imines),
    // Then also remove the bond stereo information, as it is no longer valid.
    if (heavyAtom->getDegree() == 2) {
      for (auto &nbnd : mol.atomBonds(heavyAtom)) {
        if (nbnd != bond) {
          if (nbnd->getStereo() > Bond::STEREOANY) {
            nbnd->setStereo(Bond::STEREONONE);
            nbnd->getStereoAtoms().clear();
          }
          break;
        }
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

    // remove the bond from any SGroups that might include it.
    for (auto &sg : getSubstanceGroups(mol)) {
      sg.removeBondWithIdx(bond->getIdx());
    }
  }

  // Finally, remove the atom from any SGroups that might include it, so that
  // the SGroups don't get removed in removeAtom(). Since we allow removing
  // SGroup SAP lvidx H atoms, we need to check for those and update them.
  for (auto &sg : getSubstanceGroups(mol)) {
    sg.removeAtomWithIdx(idx);
    sg.removeParentAtomWithIdx(idx);

    for (auto &sap : sg.getAttachPoints()) {
      if (sap.lvIdx == static_cast<int>(idx)) {
        sap.lvIdx = -1;
      }
    }
  }
  // computed properties will be cleared after all hydrogens are removed
  bool clearProps = false;
  mol.removeAtom(atom, clearProps);
}

bool shouldRemoveH(const RWMol &mol, const Atom *atom,
                   const RemoveHsParameters &ps) {
  if (atom->getAtomicNum() != 1) {
    return false;
  }
  if (!ps.removeWithQuery && atom->hasQuery()) {
    return false;
  }
  if (!ps.removeDegreeZero && !atom->getDegree()) {
    if (ps.showWarnings) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: not removing hydrogen atom without neighbors"
          << std::endl;
    }
    return false;
  }
  if (!ps.removeHigherDegrees && atom->getDegree() > 1) {
    return false;
  }
  if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes && atom->getIsotope()) {
    return false;
  }
  if (!ps.removeNonimplicit && !atom->hasProp(common_properties::isImplicit)) {
    return false;
  }
  if (!ps.removeMapped && atom->getAtomMapNum()) {
    return false;
  }

  if (ps.removeInSGroups) {
    // If removing H in SGroups, do not remove H atoms in special
    // roles in the SGroup
    for (const auto &sg : getSubstanceGroups(mol)) {
      // The H atom is one of the "caps" of the SGroup. Technically,
      // it's not part of the group, but it defines its boundaries.
      for (const auto &bond_idx : sg.getBonds()) {
        if (sg.getBondType(bond_idx) == SubstanceGroup::BondType::XBOND) {
          auto bond = mol.getBondWithIdx(bond_idx);
          if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
            return false;
          }
        }
      }

      for (const auto &sap : sg.getAttachPoints()) {
        // The H atoms is an attach point. This would be weird, but is possible.
        // (if it is a 'leaving atom' we don't care, though)
        if (sap.aIdx == atom->getIdx()) {
          return false;
        }
      }

      for (const auto &cs : sg.getCStates()) {
        // The bond to the H atom defines a CState
        auto bond = mol.getBondWithIdx(cs.bondIdx);
        if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
          return false;
        }
      }
    }
  } else {
    for (const auto &sg : getSubstanceGroups(mol)) {
      if (sg.includesAtom(atom->getIdx())) {
        return false;
      }
    }
  }

  if (!ps.removeHydrides && atom->getFormalCharge() == -1) {
    return false;
  }
  bool removeIt = true;
  if (atom->getDegree() &&
      (!ps.removeDummyNeighbors || !ps.removeDefiningBondStereo ||
       !ps.removeOnlyHNeighbors || !ps.removeNontetrahedralNeighbors ||
       !ps.removeWithWedgedBond)) {
    bool onlyHNeighbors = true;
    for (const auto nbr : mol.atomNeighbors(atom)) {
      // is it a dummy?
      if (!ps.removeDummyNeighbors && nbr->getAtomicNum() < 1) {
        if (ps.showWarnings) {
          BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
                                     "with dummy atom neighbors"
                                  << std::endl;
        }
        return false;
      }
      // does it have non-tetrahedral stereo:
      if (!ps.removeNontetrahedralNeighbors &&
          Chirality::hasNonTetrahedralStereo(nbr)) {
        if (ps.showWarnings) {
          BOOST_LOG(rdWarningLog)
              << "WARNING: not removing hydrogen atom "
                 "with neighbor that has non-tetrahedral stereochemistry"
              << std::endl;
        }
        return false;
      }
      if (!ps.removeOnlyHNeighbors && nbr->getAtomicNum() != 1) {
        onlyHNeighbors = false;
      }
      if (!ps.removeWithWedgedBond) {
        const auto bnd = mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
        if (bnd->getBondDir() == Bond::BEGINDASH ||
            bnd->getBondDir() == Bond::BEGINWEDGE) {
          if (ps.showWarnings) {
            BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
                                       "with wedged bond"
                                    << std::endl;
          }
          return false;
        }
      }
      // Check to see if the neighbor has a double bond and we're the only
      // neighbor at this end.  This was part of github #1810
      if (!ps.removeDefiningBondStereo && nbr->getDegree() == 2) {
        for (const auto bnd : mol.atomBonds(nbr)) {
          if (bnd->getBondType() == Bond::DOUBLE &&
              (bnd->getStereo() > Bond::STEREOANY ||
               mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
                       ->getBondDir() > Bond::NONE)) {
            return false;
          }
        }
      }
    }
    if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
      return false;
    }
  }
  return removeIt;
}

// Do not remove H atoms that are part of SGroups that only contain H atoms.
void filter_sgroup_emptying_hydrogens(const ROMol &mol,
                                      boost::dynamic_bitset<> &atomsToRemove) {
  for (const auto &sg : getSubstanceGroups(mol)) {
    const auto &atoms = sg.getAtoms();
    const auto &patoms = sg.getParentAtoms();

    // If the SGroup already didn't have atoms, we don't care about it
    if (atoms.empty() && patoms.empty()) {
      continue;
    }

    auto would_remove_atom = [&atomsToRemove](const auto idx) {
      return atomsToRemove[idx];
    };

    auto no_atoms = atoms.empty() ||
                    std::all_of(atoms.begin(), atoms.end(), would_remove_atom);
    if (no_atoms) {
      auto no_patoms =
          patoms.empty() ||
          std::all_of(patoms.begin(), patoms.end(), would_remove_atom);
      if (no_patoms) {
        for (auto atom : atoms) {
          atomsToRemove.set(atom, false);
        }
        for (auto patom : patoms) {
          atomsToRemove.set(patom, false);
        }
      }
    }
  }
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
          ->setProp(common_properties::_isotopicHs, pair.second);
    }
  }
  boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};

  for (auto atom : mol.atoms()) {
    if (shouldRemoveH(mol, atom, ps)) {
      atomsToRemove.set(atom->getIdx());
    }
  }  // end of the loop over atoms

  // Once we know which H atoms would be removed, filter out those that
  // would cause any SGroups to become empty
  if (ps.removeInSGroups) {
    filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
  }

  // now that we know which atoms need to be removed, go ahead and remove them
  // NOTE: there's too much complexity around stereochemistry here
  // to be able to safely use batch editing.
  for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    if (atomsToRemove[idx]) {
      molRemoveH(mol, idx, ps.updateExplicitCount);
    }
  }
  mol.clearComputedProps(true);
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
  ps.removeNontetrahedralNeighbors = true;
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
enum class HydrogenType {
  NotAHydrogen,
  UnMergableQueryHydrogen,
  QueryHydrogen
};

HydrogenType isQueryH(const Atom *atom) {
  PRECONDITION(atom, "bogus atom");
  if (atom->getAtomicNum() == 1) {
    // the simple case: the atom is flagged as being an H and
    // has no query
    if (!atom->hasQuery() ||
        (!atom->getQuery()->getNegation() &&
         atom->getQuery()->getDescription() == "AtomAtomicNum")) {
      return HydrogenType::QueryHydrogen;
    }
  }

  if (!(atom->getDegree() <= 1)) {
    // bonded and unbonded H atoms will continue rest will be returned
    return HydrogenType::NotAHydrogen;
  }

  if (atom->hasQuery() && atom->getQuery()->getNegation()) {
    // we will not merge negated queries
    return HydrogenType::NotAHydrogen;
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
      return HydrogenType::UnMergableQueryHydrogen;
    }
  }
  return hasHQuery ? HydrogenType::QueryHydrogen : HydrogenType::NotAHydrogen;
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
void mergeQueryHs(RWMol &mol, bool mergeUnmappedOnly, bool mergeIsotopes) {
  std::vector<unsigned int> atomsToRemove;

  boost::dynamic_bitset<> hatoms(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    hatoms[i] = isQueryH(mol.getAtomWithIdx(i)) == HydrogenType::QueryHydrogen;
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
          bool checkUnmapped =
              !mergeUnmappedOnly ||
              !bgn.hasProp(common_properties::molAtomMapNumber);
          bool checkIsotope = mergeIsotopes || bgn.getIsotope() == 0;
          if (checkUnmapped && checkIsotope) {
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
        if (atom->getQuery()->getDescription() == "RecursiveStructure") {
          auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
          CHECK_INVARIANT(rsq, "could not convert recursive structure query");
          RWMol *rqm = new RWMol(*rsq->getQueryMol());
          mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
          rsq->setQueryMol(rqm);
        }

        // FIX: shouldn't be repeating this code here
        std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
            atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
        while (childStack.size()) {
          QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
          childStack.pop_front();
          if (qry->getDescription() == "RecursiveStructure") {
            auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
            CHECK_INVARIANT(rsq, "could not convert recursive structure query");
            RWMol *rqm = new RWMol(*rsq->getQueryMol());
            mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
            rsq->setQueryMol(rqm);
          } else if (qry->beginChildren() != qry->endChildren()) {
            childStack.insert(childStack.end(), qry->beginChildren(),
                              qry->endChildren());
          }
        }
      }  // end of recursion loop
    }
    ++currIdx;
  }
  mol.beginBatchEdit();
  for (auto aidx : atomsToRemove) {
    mol.removeAtom(aidx);
  }
  mol.commitBatchEdit();
};
ROMol *mergeQueryHs(const ROMol &mol, bool mergeUnmappedOnly,
                    bool mergeIsotopes) {
  auto *res = new RWMol(mol);
  mergeQueryHs(*res, mergeUnmappedOnly, mergeIsotopes);
  return static_cast<ROMol *>(res);
};

bool needsHs(const ROMol &mol) {
  for (const auto atom : mol.atoms()) {
    unsigned int nHNbrs = 0;
    for (const auto nbri :
         boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
      const auto nbr = mol[nbri];
      if (nbr->getAtomicNum() == 1) {
        ++nHNbrs;
      }
    }
    bool noNeighbors = false;
    if (atom->getTotalNumHs(noNeighbors) > nHNbrs) {
      return true;
    }
  }
  return false;
}

std::pair<bool, bool> hasQueryHs(const ROMol &mol) {
  bool queryHs = false;
  // We don't care about announcing ORs or other items during isQueryH
  RDLog::LogStateSetter blocker;

  for (const auto atom : mol.atoms()) {
    switch (isQueryH(atom)) {
      case HydrogenType::UnMergableQueryHydrogen:
        return std::make_pair(true, true);
      case HydrogenType::QueryHydrogen:
        queryHs = true;
        break;
      default:  // HydrogenType::NotAHydrogen:
        break;
    }
    if (atom->hasQuery()) {
      if (atom->getQuery()->getDescription() == "RecursiveStructure") {
        auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
        CHECK_INVARIANT(rsq, "could not convert recursive structure query");
        auto res = hasQueryHs(*rsq->getQueryMol());
        if (res.second) {  // unmergableH implies queryH
          return res;
        }
        queryHs |= res.first;
      }

      // FIX: shouldn't be repeating this code here -- yet again!
      std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
          atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
      while (!childStack.empty()) {
        QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
        childStack.pop_front();
        if (qry->getDescription() == "RecursiveStructure") {
          auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
          CHECK_INVARIANT(rsq, "could not convert recursive structure query");
          auto res = hasQueryHs(*rsq->getQueryMol());
          if (res.second) {
            return res;
          }
          queryHs |= res.first;
        } else {
          childStack.insert(childStack.end(), qry->beginChildren(),
                            qry->endChildren());
        }
      }
    }
  }  // end of recursion loop

  return std::make_pair(queryHs, false);
}

}  // namespace MolOps
}  // namespace RDKit
