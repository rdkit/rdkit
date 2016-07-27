// $Id$
//
//  Copyright (C) 2003-2013 Greg Landrum and Rational Discovery LLC
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
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <boost/foreach.hpp>

namespace RDKit {

// Local utility functionality:
namespace {
Atom *getAtomNeighborNot(ROMol *mol, const Atom *atom, const Atom *other) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  PRECONDITION(atom->getDegree() > 1, "bad degree");
  PRECONDITION(other, "bad atom");
  Atom *res = 0;

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

  const Atom *nbr1 = 0, *nbr2 = 0, *nbr3 = 0;
  const Bond *nbrBond;
  ROMol::ADJ_ITER nbrIdx, endNbrs;

  switch (heavyAtom->getDegree()) {
    case 1:
      // --------------------------------------------------------------------------
      //   No other atoms present:
      // --------------------------------------------------------------------------
      dirVect.z = 1;
      // loop over the conformations and set the coordinates
      for (ROMol::ConformerIterator cfi = mol->beginConformers();
           cfi != mol->endConformers(); cfi++) {
        heavyPos = (*cfi)->getAtomPos(heavyIdx);
        hydPos = heavyPos + dirVect * bondLength;
        (*cfi)->setAtomPos(hydIdx, hydPos);
      }
      break;

    case 2:
      // --------------------------------------------------------------------------
      //  One other neighbor:
      // --------------------------------------------------------------------------
      nbr1 = getAtomNeighborNot(mol, heavyAtom, hydAtom);
      for (ROMol::ConformerIterator cfi = mol->beginConformers();
           cfi != mol->endConformers(); ++cfi) {
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
            perpVect = nbr1Vect.getPerpendicular();
            // and move off it:
            tform.SetRotation((180 - 109.471) * M_PI / 180., perpVect);
            dirVect = tform * nbr1Vect;
            hydPos = heavyPos + dirVect * bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP2:
            // default position is to just take an arbitrary perpendicular:
            perpVect = nbr1Vect.getPerpendicular();

            if (nbr1->getDegree() > 1) {
              // can we use the neighboring atom to establish a perpendicular?
              nbrBond = mol->getBondBetweenAtoms(heavyIdx, nbr1->getIdx());
              if (nbrBond->getIsAromatic() ||
                  nbrBond->getBondType() == Bond::DOUBLE) {
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
            hydPos = heavyPos + dirVect * bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP:
            // just lay the H along the vector:
            dirVect = nbr1Vect;
            hydPos = heavyPos + dirVect * bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the vector:
            dirVect = nbr1Vect;
            hydPos = heavyPos + dirVect * bondLength;
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
          if (!nbr1)
            nbr1 = mol->getAtomWithIdx(*nbrIdx);
          else
            nbr2 = mol->getAtomWithIdx(*nbrIdx);
        }
        ++nbrIdx;
      }
      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      for (ROMol::ConformerIterator cfi = mol->beginConformers();
           cfi != mol->endConformers(); ++cfi) {
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
            hydPos = heavyPos + dirVect * bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          case Atom::SP2:
            // don't need to do anything here, the H atom goes right on the
            // direction vector
            hydPos = heavyPos + dirVect * bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the neighbor vector;
            hydPos = heavyPos + dirVect * bondLength;
            (*cfi)->setAtomPos(hydIdx, hydPos);
            break;
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
        std::vector<std::pair<unsigned int, int> > nbrs;
        while (nbrIdx != endNbrs) {
          if (*nbrIdx != hydIdx) {
            const Atom *tAtom = mol->getAtomWithIdx(*nbrIdx);
            unsigned int cip = 0;
            tAtom->getPropIfPresent<unsigned int>(common_properties::_CIPRank,
                                                  cip);
            nbrs.push_back(std::make_pair(cip, rdcast<int>(*nbrIdx)));
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
      for (ROMol::ConformerIterator cfi = mol->beginConformers();
           cfi != mol->endConformers(); ++cfi) {
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
        if (fabs(nbr3Vect.dotProduct(nbr1Vect.crossProduct(nbr2Vect))) < 0.1) {
          if ((*cfi)->is3D()) {
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
        } else {
          dirVect = nbr1Vect + nbr2Vect + nbr3Vect;
        }
        dirVect.normalize();
        hydPos = heavyPos + dirVect * bondLength;
        (*cfi)->setAtomPos(hydIdx, hydPos);
      }
      break;
    default:
      // --------------------------------------------------------------------------
      // FIX: figure out what to do here
      // --------------------------------------------------------------------------
      hydPos = heavyPos + dirVect * bondLength;
      for (ROMol::ConformerIterator cfi = mol->beginConformers();
           cfi != mol->endConformers(); ++cfi) {
        (*cfi)->setAtomPos(hydIdx, hydPos);
      }
      break;
  }
}
}  // end of unnamed namespace

namespace MolOps {
void addHs(RWMol &mol, bool explicitOnly, bool addCoords,
           const UINT_VECT *onlyOnAtoms) {
  // when we hit each atom, clear its computed properties
  // NOTE: it is essential that we not clear the ring info in the
  // molecule's computed properties.  We don't want to have to
  // regenerate that.  This caused Issue210 and Issue212:
  mol.clearComputedProps(false);

  // precompute the number of hydrogens we are going to add so that we can
  // pre-allocate the necessary space on the conformations of the molecule
  // for their coordinates
  unsigned int numAddHyds = 0;
  for (ROMol::AtomIterator at = mol.beginAtoms(); at != mol.endAtoms(); ++at) {
    if (!onlyOnAtoms ||
        std::find(onlyOnAtoms->begin(), onlyOnAtoms->end(), (*at)->getIdx()) !=
            onlyOnAtoms->end()) {
      numAddHyds += (*at)->getNumExplicitHs();
      if (!explicitOnly) {
        numAddHyds += (*at)->getNumImplicitHs();
      }
    }
  }
  unsigned int nSize = mol.getNumAtoms() + numAddHyds;

  // loop over the conformations of the molecule and allocate new space
  // for the H locations (need to do this even if we aren't adding coords so
  // that the conformers have the correct number of atoms).
  for (ROMol::ConformerIterator cfi = mol.beginConformers();
       cfi != mol.endConformers(); ++cfi) {
    (*cfi)->reserve(nSize);
  }

  unsigned int stopIdx = mol.getNumAtoms();
  for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    if (onlyOnAtoms &&
        std::find(onlyOnAtoms->begin(), onlyOnAtoms->end(), aidx) ==
            onlyOnAtoms->end()) {
      continue;
    }

    Atom *newAt = mol.getAtomWithIdx(aidx);
    unsigned int newIdx;
    newAt->clearComputedProps();
    // always convert explicit Hs
    unsigned int onumexpl = newAt->getNumExplicitHs();
    for (unsigned int i = 0; i < onumexpl; i++) {
      newIdx = mol.addAtom(new Atom(1), false, true);
      mol.addBond(aidx, newIdx, Bond::SINGLE);
      mol.getAtomWithIdx(newIdx)->updatePropertyCache();
      if (addCoords) setHydrogenCoords(&mol, newIdx, aidx);
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
        mol.getAtomWithIdx(newIdx)->setProp(common_properties::isImplicit, 1);
        mol.getAtomWithIdx(newIdx)->updatePropertyCache();
        if (addCoords) setHydrogenCoords(&mol, newIdx, aidx);
      }
      // be very clear about implicits not being allowed in this representation
      newAt->setProp(common_properties::origNoImplicit, newAt->getNoImplicit(),
                     true);
      newAt->setNoImplicit(true);
    }
    // update the atom's derived properties (valence count, etc.)
    newAt->updatePropertyCache();
  }
}
ROMol *addHs(const ROMol &mol, bool explicitOnly, bool addCoords,
             const UINT_VECT *onlyOnAtoms) {
  RWMol *res = new RWMol(mol);
  addHs(*res, explicitOnly, addCoords, onlyOnAtoms);
  return static_cast<ROMol *>(res);
};

//
//  This routine removes hydrogens (and bonds to them) from the molecular graph.
//  Other Atom and bond indices may be affected by the removal.
//
//  NOTES:
//   - Hydrogens which aren't connected to a heavy atom will not be
//     removed.  This prevents molecules like "[H][H]" from having
//     all atoms removed.
//   - Labelled hydrogen (e.g. atoms with atomic number=1, but isotope > 1),
//     will not be removed.
//   - two coordinate Hs, like the central H in C[H-]C, will not be removed
//   - Hs connected to dummy atoms will not be removed
//
void removeHs(RWMol &mol, bool implicitOnly, bool updateExplicitCount,
              bool sanitize) {
  unsigned int currIdx = 0, origIdx = 0;
  std::map<unsigned int, unsigned int> idxMap;
  for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
       ++atIt) {
    if ((*atIt)->getAtomicNum() == 1) continue;
    (*atIt)->updatePropertyCache(false);
  }
  while (currIdx < mol.getNumAtoms()) {
    Atom *atom = mol.getAtomWithIdx(currIdx);
    idxMap[origIdx] = currIdx;
    ++origIdx;
    if (atom->getAtomicNum() == 1) {
      bool removeIt = false;

      if (atom->hasProp(common_properties::isImplicit)) {
        removeIt = true;
      } else if (!implicitOnly && !atom->getIsotope() &&
                 atom->getDegree() == 1) {
        ROMol::ADJ_ITER begin, end;
        boost::tie(begin, end) = mol.getAtomNeighbors(atom);
        if (mol.getAtomWithIdx(*begin)->getAtomicNum() > 1) {
          removeIt = true;
        }
      }

      if (removeIt) {
        ROMol::OEDGE_ITER beg, end;
        boost::tie(beg, end) = mol.getAtomBonds(atom);
        // note the assumption that the H only has one neighbor... I
        // feel no need to handle the case of hypervalent hydrogen!
        // :-)
        const BOND_SPTR bond = mol[*beg];
        Atom *heavyAtom = bond->getOtherAtom(atom);
        int heavyAtomNum = heavyAtom->getAtomicNum();
        const INT_VECT &defaultVs =
            PeriodicTable::getTable()->getValenceList(heavyAtomNum);

        // we'll update the atom's explicit H count if we were told to
        // *or* if the atom is chiral, in which case the H is needed
        // in order to complete the coordination
        // *or* if the atom has the noImplicit flag set:
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
          boost::tie(beg, end) = mol.getAtomBonds(heavyAtom);
          while (beg != end) {
            if (mol[*beg]->getIdx() != bond->getIdx()) {
              neighborIndices.push_back(mol[*beg]->getIdx());
            }
            ++beg;
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
        }

        // if the direction is set on this bond and the atom it's connected to
        // has no other single bonds with directions set, then we need to set
        // direction on one of the other neighbors in order to avoid double bond
        // stereochemistry possibly being lost. This was github #754
        if (bond->getBondDir() == Bond::ENDDOWNRIGHT ||
            bond->getBondDir() == Bond::ENDUPRIGHT) {
          bool foundADir = false;
          Bond *oBond = NULL;
          boost::tie(beg, end) = mol.getAtomBonds(heavyAtom);
          while (beg != end) {
            if (mol[*beg]->getIdx() != bond->getIdx() &&
                mol[*beg]->getBondType() == Bond::SINGLE) {
              if (mol[*beg]->getBondDir() == Bond::NONE) {
                oBond = mol[*beg].get();
              } else {
                foundADir = true;
              }
            }
            ++beg;
          }
          if (!foundADir && oBond != NULL) {
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
        }

        mol.removeAtom(atom);
      } else {
        // only increment the atom idx if we don't remove the atom
        currIdx++;
      }
    } else {
      // only increment the atom idx if we don't remove the atom
      currIdx++;
      bool origNoImplicit;
      if (atom->getPropIfPresent(common_properties::origNoImplicit,
                                 origNoImplicit)) {
        // we'll get in here if we haven't already processed the atom's implicit
        //  hydrogens. (this is protection for the case that removeHs() is
        //  called
        //  multiple times on a single molecule without intervening addHs()
        //  calls)
        atom->setNoImplicit(origNoImplicit);
        atom->clearProp(common_properties::origNoImplicit);
      }
    }
  }
  //
  //  If we didn't only remove implicit Hs, which are guaranteed to
  //  be the highest numbered atoms, we may have altered atom indices.
  //  This can screw up derived properties (such as ring members), so
  //  do some checks:
  //
  if (!implicitOnly) {
    if (sanitize) {
      sanitizeMol(mol);
    }
  }
};

ROMol *removeHs(const ROMol &mol, bool implicitOnly, bool updateExplicitCount,
                bool sanitize) {
  RWMol *res = new RWMol(mol);
  try {
    removeHs(*res, implicitOnly, updateExplicitCount, sanitize);
  } catch (MolSanitizeException &se) {
    delete res;
    throw se;
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
    // the logic gets too complicated if there's an OR in the children, so just
    // punt on those (with a warning)
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
}

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
        //   If we have no H query already:
        //        - add a generic H query
        //      else:
        //        - do nothing
        //
        //  Examples:
        //    C[H] -> [C;!H0]
        //    [C;H1][H] -> [C;H1]
        //    [C;H2][H] -> [C;H2]
        //
        // FIX: this is going to behave oddly in the case of a contradictory
        //  SMARTS like: [C;H0][H], where it will give the equivalent of:
        //  [C;H0]  I think this is actually correct, but I can be persuaded
        //  otherwise.
        //
        //  First we'll search for an H query:
        bool hasHQuery = false;
        if (!atom->hasQuery()) {
          // it wasn't a query atom, we need to replace it so that we can add a
          // query:
          ATOM_EQUALS_QUERY *tmp = makeAtomNumQuery(atom->getAtomicNum());
          QueryAtom *newAt = new QueryAtom;
          newAt->setQuery(tmp);
          mol.replaceAtom(atom->getIdx(), newAt);
          delete newAt;
          atom = mol.getAtomWithIdx(currIdx);
        }
        if (!hasHQuery) {
          for (unsigned int i = 0; i < numHsToRemove; ++i) {
            ATOM_EQUALS_QUERY *tmp = makeAtomHCountQuery(i);
            tmp->setNegation(true);
            atom->expandQuery(tmp);
          }
        }
      }  // end of numHsToRemove test

      // recurse if needed (was github isusue 544)
      if (atom->hasQuery()) {
        // std::cerr<<"  q: "<<atom->getQuery()->getDescription()<<std::endl;
        if (atom->getQuery()->getDescription() == "RecursiveStructure") {
          RWMol *rqm = static_cast<RWMol *>(const_cast<ROMol *>(
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
            RWMol *rqm = static_cast<RWMol *>(const_cast<ROMol *>(
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
  RWMol *res = new RWMol(mol);
  mergeQueryHs(*res, mergeUnmappedOnly);
  return static_cast<ROMol *>(res);
};

};  // end of namespace MolOps
};  // end of namespace RDKit
