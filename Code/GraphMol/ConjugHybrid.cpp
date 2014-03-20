// $Id$
//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ROMol.h"
#include "RWMol.h"
#include "Atom.h"
#include "Bond.h"
#include "MolOps.h"
#include "PeriodicTable.h"
#include "AtomIterators.h"
#include "BondIterators.h"


namespace RDKit {
  // local utility namespace:
  namespace {
    bool isAtomConjugCand(const Atom *at){
      // the second check here is for Issue211, where the c-P bonds in
      // Pc1ccccc1 were being marked as conjugated.  This caused the P atom
      // itself to be SP2 hybridized.  This is wrong.  For now we'll do a quick
      // hack and forbid this check from adding conjugation to anything out of
      // the first row of the periodic table.  (Conjugation in aromatic rings
      // has already been attended to, so this is safe.)
      int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
      if (((at->getAtomicNum() <= 10) || (nouter != 5 && nouter != 6)) &&
          (MolOps::countAtomElec(at) > 0)
          ) {
        return true;
      }
      return false;
    }

    void markConjAtomBonds(Atom *at) {
      if(!isAtomConjugCand(at)) return;
      ROMol &mol = at->getOwningMol();
      Atom* at2;

      int atx = at->getIdx();
      // make sure that have either 2 or 3 subtitutions on this atom
      int sbo = at->getDegree() + at->getTotalNumHs();
      if ( (sbo < 2) || (sbo > 3) ) {
	return;
      }

      ROMol::OEDGE_ITER bnd1, end1, bnd2, end2;
      boost::tie(bnd1,end1) = mol.getAtomBonds(at);
      while (bnd1 != end1) {
	if (mol[*bnd1]->getValenceContrib(at) < 1.5) {
	  bnd1++;
	  continue;
	}
	boost::tie(bnd2,end2) = mol.getAtomBonds(at);
	while (bnd2 != end2) {
	  if (bnd1 == bnd2) {
	    bnd2++;
	    continue;
	  }
	  at2 = mol.getAtomWithIdx(mol[*bnd2]->getOtherAtomIdx(atx));
	  sbo = at2->getDegree() + at2->getTotalNumHs();
	  if (sbo > 3) {
	    bnd2++;
	    continue;
	  }
	  if (isAtomConjugCand(at2)) {
	    mol[*bnd1]->setIsConjugated(true);
	    mol[*bnd2]->setIsConjugated(true);
	  }
	  bnd2++;
	}
	bnd1++;
      }
    }

    int numBondsPlusLonePairs(Atom *at) {
      PRECONDITION(at,"bad atom");
      int deg = at->getTotalDegree();
      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
      while(beg!=end){
        BOND_SPTR bond=at->getOwningMol()[*beg];
        if(bond->getBondType()==Bond::ZERO) --deg;
        ++beg;
      }

      if (at->getAtomicNum()<=1){
        return deg;
      }
      int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
      int totalValence = at->getExplicitValence() + at->getImplicitValence();
      int chg = at->getFormalCharge();

      int numFreeElectrons=nouter - (totalValence+chg);
      if(totalValence + nouter - chg < 8){
        // we're below an octet, so we need to think
        // about radicals:
        int numRadicals = at->getNumRadicalElectrons();
        int numLonePairs = (numFreeElectrons - numRadicals)/2;
        return deg + numLonePairs + numRadicals;
      } else {
        int numLonePairs= numFreeElectrons/2;
        return deg + numLonePairs;
      }
    }
  } //end of utility namespace


  namespace MolOps {
    bool atomHasConjugatedBond(const Atom *at){
      PRECONDITION(at,"bad atom");

      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
      while(beg!=end){
	if(at->getOwningMol()[*beg]->getIsConjugated()) return true;
	beg++;
      }
      return false;
    }

    void setConjugation(ROMol &mol) {
    
      // start with all bonds being marked unconjugated
      // except for aromatic bonds
      ROMol::BondIterator bi;
      for (bi = mol.beginBonds(); bi != mol.endBonds(); bi++) {
	if ((*bi)->getIsAromatic()) {
	  (*bi)->setIsConjugated(true);
	}
	else {
	  (*bi)->setIsConjugated(false);
	}
      }

      ROMol::AtomIterator ai;
      // loop over each atom and check if the bonds connecting to it can 
      // be conjugated
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
	markConjAtomBonds(*ai);
      }
    
    }

  
    void setHybridization(ROMol &mol) {
      ROMol::AtomIterator ai;
      int norbs;
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        if((*ai)->getAtomicNum()==0){
	  (*ai)->setHybridization(Atom::UNSPECIFIED);
        } else {
          norbs = numBondsPlusLonePairs(*ai);
          switch(norbs) {
          case 0:
            // This occurs for things like Na+
            (*ai)->setHybridization(Atom::S);
            break;
          case 1:
            (*ai)->setHybridization(Atom::S);
            break;
          case 2:
            (*ai)->setHybridization(Atom::SP);
            break;
          case 3:
            (*ai)->setHybridization(Atom::SP2);
            break;
          case 4:
            // potentially SP3, but we'll set it down to SP2
            // if we have a conjugated bond (like the second O
            // in O=CO)
            // we'll also avoid setting the hybridization down to
            // SP2 in the case of an atom with degree higher than 3
            // (e.g. things like CP1(C)=CC=CN=C1C, where the P
            //   has norbs = 4, and a conjugated bond, but clearly should
            //   not be SP2)
            // This is Issue276
            if((*ai)->getDegree()>3 || !MolOps::atomHasConjugatedBond(*ai)){
              (*ai)->setHybridization(Atom::SP3);
            } else {
              (*ai)->setHybridization(Atom::SP2);
            }
            break;
          case 5:
            (*ai)->setHybridization(Atom::SP3D);
            break;
          case 6:
            (*ai)->setHybridization(Atom::SP3D2);
            break;
          default :
            (*ai)->setHybridization(Atom::UNSPECIFIED);
          }
        }
      }
    }
  } // end of namespace MolOps
} // end of namespace RDKit
      
