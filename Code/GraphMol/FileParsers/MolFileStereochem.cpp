// $Id$
//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
//
#include <list>
#include <RDGeneral/RDLog.h>
#include "MolFileStereochem.h"
#include <Geometry/point.h>

namespace RDKit {
  typedef std::list<double> DOUBLE_LIST;

  // ----------------------------------- -----------------------------------
  // This algorithm is identical to that used in the CombiCode Mol file
  //  parser (also developed by RD).
  //
  //
  // SUMMARY:
  //   Derive a chiral code for an atom that has a wedged (or dashed) bond
  //   drawn to it.
  //
  // RETURNS:
  //   The chiral type
  //
  // CAVEATS:
  //   This is careful to ensure that the central atom has 4 neighbors and
  //   only single bonds to it, but that's about it.
  // 
  // NOTE: this isn't careful at all about checking to make sure that
  // things actually *should* be chiral. e.g. if the file has a
  // 3-coordinate N with a wedged bond, it will make some erroneous
  // assumptions about the chirality.
  //   
  // ----------------------------------- -----------------------------------

  Atom::ChiralType FindAtomStereochemistry(const RWMol &mol,const Bond *bond, 
                                           const Conformer *conf){
    PRECONDITION(bond,"no bond");
    PRECONDITION(conf,"no conformer");
    Bond::BondDir bondDir=bond->getBondDir();
    PRECONDITION(bondDir==Bond::BEGINWEDGE || bondDir==Bond::BEGINDASH,
                 "bad bond direction");

    // NOTE that according to the CT file spec, wedging assigns chirality
    // to the atom at the point of the wedge, (atom 1 in the bond).
    const Atom *atom = bond->getBeginAtom();
    PRECONDITION(atom,"no atom");
    const Atom *bondAtom = bond->getEndAtom();

    Atom::ChiralType res=Atom::CHI_UNSPECIFIED;

    INT_LIST neighborBondIndices;
    RDGeom::Point3D centerLoc, tmpPt;
    centerLoc=conf->getAtomPos(atom->getIdx());
    tmpPt=conf->getAtomPos(bondAtom->getIdx());
    centerLoc.z=0.0;
    tmpPt.z = 0.0;

    RDGeom::Point3D refVect=centerLoc.directionVector(tmpPt);

    //----------------------------------------------------------
    //
    //  start by ensuring that all the bonds to neighboring atoms
    //  are single bonds and collecting a list of neighbor indices:
    //
    //----------------------------------------------------------
    bool hSeen=false;

    neighborBondIndices.push_back(bond->getIdx());
    if(bondAtom->getAtomicNum()==1) hSeen=true;

    bool allSingle=true;
    bool hasTruePrecedingAtom=false;
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = mol.getAtomBonds(atom);
    ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap = mol.getBondPMap();
    while(beg!=end){
      Bond *nbrBond=pMap[*beg];
      if(nbrBond->getBondType() != Bond::SINGLE){
        allSingle=false;
        break;
      }
      if(nbrBond != bond){
        if(nbrBond->getOtherAtom(atom)->getAtomicNum()==1) hSeen=true;
        neighborBondIndices.push_back(nbrBond->getIdx());
      }

      unsigned int oIdx=nbrBond->getOtherAtomIdx(atom->getIdx());
      if(oIdx<atom->getIdx() && nbrBond->getBeginAtomIdx()==oIdx){
        hasTruePrecedingAtom=true;
      }
      beg++;
    }
    int nNbrs = neighborBondIndices.size();

    //----------------------------------------------------------
    //
    //  Return now if there aren't at least 3 non-H bonds to the atom.
    //  (we can implicitly add a single H to 3 coordinate atoms, but
    //  we're horked otherwise).
    //
    //----------------------------------------------------------
    if(nNbrs<3 || (hSeen && nNbrs<4) ){
      return Atom::CHI_UNSPECIFIED;
    }

    //----------------------------------------------------------
    //
    //  Continue if there are all single bonds:
    //
    //----------------------------------------------------------
    if(allSingle){
      //------------------------------------------------------------
      //
      //  Here we need to figure out the rotation direction between
      //  the neighbor bonds and the wedged bond:
      //
      //------------------------------------------------------------
      bool isCCW;
      double angle0,angle1,angle2;
      RDGeom::Point3D atomVect;
      INT_LIST::const_iterator bondIter=neighborBondIndices.begin();
      bondIter++;
      int oaid = mol.getBondWithIdx(*bondIter)->getOtherAtom(atom)->getIdx();
      tmpPt = conf->getAtomPos(oaid);
      atomVect = centerLoc.directionVector(tmpPt);
      angle0 = refVect.signedAngleTo(atomVect);
      if(angle0<0) angle0 += 2.*M_PI;

      bondIter++;
      oaid = mol.getBondWithIdx(*bondIter)->getOtherAtom(atom)->getIdx();
      tmpPt = conf->getAtomPos(oaid);
      atomVect = centerLoc.directionVector(tmpPt);
      angle1 = refVect.signedAngleTo(atomVect);
      if(angle1<0) angle1 += 2.*M_PI;

      // We proceed differently for 3 and 4 coordinate atoms:
      double firstAngle,secondAngle;
      if(nNbrs==4){
        bool flipIt=false;
        // grab the angle to the last neighbor:
        bondIter++;
        oaid = mol.getBondWithIdx(*bondIter)->getOtherAtom(atom)->getIdx();
        tmpPt = conf->getAtomPos(oaid);
        atomVect = centerLoc.directionVector(tmpPt);
        angle2 = refVect.signedAngleTo(atomVect);
        if(angle2<0) angle2 += 2.*M_PI;

        // find the lowest and second-lowest angle and keep track of
        // whether or not we have to do a non-cyclic permutation to
        // get there:
        if(angle0<angle1){
          if(angle1<angle2){
            // order is angle0 -> angle1 -> angle2
            firstAngle = angle0;
            secondAngle = angle1;
          } else if(angle0 < angle2){
            // order is angle0 -> angle2 -> angle1
            firstAngle = angle0;
            secondAngle = angle2;
            flipIt=true;
          } else {
            // order is angle2 -> angle0 -> angle1
            firstAngle = angle2;
            secondAngle = angle0;
          }
        } else if(angle0 < angle2) {
          // order is angle1 -> angle0 -> angle2
          firstAngle=angle1;
          secondAngle=angle0;
          flipIt=true;
        } else {
          if(angle1 < angle2){
            // order is angle1 -> angle2 -> angle0
            firstAngle=angle1;
            secondAngle=angle2;
          } else {
            // order is angle2 -> angle1 -> angle0
            firstAngle=angle2;
            secondAngle=angle1;
            flipIt=true;
          }
        }
        if(secondAngle - firstAngle < M_PI){
          isCCW = true;
        } else {
          isCCW = false;
        }
        if(flipIt){
          isCCW = !isCCW;
        }
      } else {
        // it's three coordinate.  Things are a bit different here
        // because we have to at least kind of figure out where the
        // hydrogen might be:
        if(angle0<angle1){
          firstAngle = angle0;
          secondAngle = angle1;
          isCCW = true;
        } else {
          firstAngle = angle1;
          secondAngle = angle0;
          isCCW = false;
        }
        if(secondAngle-firstAngle >= (M_PI-1e-4)){
          // it's a situation like one of these:
          //      
          //      0        1 0 2
          //      *         \*/
          //  1 - C - 2      C
          //
          // In each of these cases, the implicit H is between atoms 1
          // and 2, so we need to flip the rotation direction (go
          // around the back).
          isCCW = !isCCW;
        }
      }
      // reverse the rotation direction if the reference is wedged down:
      if(bondDir==Bond::BEGINDASH){
        isCCW = !isCCW;
      }

      // ----------------
      // 
      // We now have the rotation direction using mol-file order.
      // We need to convert that into the appropriate label for the
      // central atom
      // 
      // ----------------
      int nSwaps = atom->getPerturbationOrder(neighborBondIndices);
      if(nNbrs==3 && !hasTruePrecedingAtom)++nSwaps;
      if(nSwaps%2) isCCW = !isCCW;
      if(isCCW) res = Atom::CHI_TETRAHEDRAL_CCW;
      else res = Atom::CHI_TETRAHEDRAL_CW;
    }

    return res;
  }
      
  void ComputeBondStereoChemistry(ROMol &mol, Bond *dblBond,
                                  const Conformer *conf) {
    // ok this function got a lot easier than before 
    //   - now we do not have to set the directions on the single
    //     bonds surrounding a double bond
    //   - we simply have to figure the if the high ranking neighbor
    //     of the double bond are in cis or trans position and set 
    //     the BondStereo on the Bond
    // we want to deal only with double bonds:
    PRECONDITION(dblBond, "");
    PRECONDITION(dblBond->getBondType() == Bond::DOUBLE, "");
    PRECONDITION(conf,"no conformer");

    const INT_VECT &ctNbrs = dblBond->getStereoAtoms();
    if (ctNbrs.size() >= 2) {
      int cnbr1 = ctNbrs[0];
      int cnbr2 = ctNbrs[1];
      RDGeom::Point3D begNbrLoc, begAtmLoc, endAtmLoc, endNbrLoc;
      begNbrLoc = conf->getAtomPos(cnbr1);
      begAtmLoc = conf->getAtomPos(dblBond->getBeginAtomIdx());
      endAtmLoc = conf->getAtomPos(dblBond->getEndAtomIdx());
      endNbrLoc = conf->getAtomPos(cnbr2);
      double ang = RDGeom::computeDihedralAngle(begNbrLoc, begAtmLoc,
                                                endAtmLoc, endNbrLoc);
      if (ang < RDKit::PI/2) {
        dblBond->setStereo(Bond::STEREOZ);
      } else {
        dblBond->setStereo(Bond::STEREOE);
      }
    }
  }


  void WedgeMolBonds(ROMol &mol, const Conformer *conf){
    PRECONDITION(conf,"no conformer");
    INT_MAP_INT wedgeBonds=pickBondsToWedge(mol);
    for(ROMol::BondIterator bondIt=mol.beginBonds();
        bondIt!=mol.endBonds();
        bondIt++){
      Bond *bond=*bondIt;
      if(bond->getBondType()==Bond::SINGLE) {
        Bond::BondDir dir=DetermineBondWedgeState(bond,wedgeBonds,conf);
        if(dir==Bond::BEGINWEDGE || dir==Bond::BEGINDASH){
          bond->setBondDir(dir);
        }
      }
    }

  }
  
  
  INT_MAP_INT pickBondsToWedge(const ROMol &mol) {
    // picks a bond for each atom that we will wedge when we write the mol file
    // here is what we are going to do
    // - at each chiral center look for a bond that is begins at the atom and
    //   is not yet picked to be wedged for a different chiral center
    // - if we do not find a bond that begins at the chiral center - we will take
    //   the first bond that is not yet picked by any other chiral centers
    ROMol::ConstAtomIterator cai;
    INT_MAP_INT res;
    Atom::ChiralType type;
    RDKit::ROMol::OBOND_ITER_PAIR atomBonds;
    RDKit::ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap = mol.getBondPMap();
    int bid;
    const Bond *bond;
    for (cai = mol.beginAtoms(); cai != mol.endAtoms(); cai++) {
      type = (*cai)->getChiralTag();
      if ((type == Atom::CHI_TETRAHEDRAL_CW) || (type == Atom::CHI_TETRAHEDRAL_CCW)) {
        int pick = -1;
        atomBonds = mol.getAtomBonds(*cai);
        while (atomBonds.first != atomBonds.second ){
          bond = pMap[*atomBonds.first];
          bid = bond->getIdx();
          if (res.find(bid) == res.end()) {
            if(bond->getBeginAtomIdx() == (*cai)->getIdx() &&
               mol.getRingInfo()->numBondRings(bid)==0 ){
              // it's a non-ring bond starting at this atom, that's
              // enough to declare ourselves finished here.
              pick = bid;
              break;
            } else if (pick == -1) {
              // be sure we get at least something:
              pick = bid;
            }
          }
          atomBonds.first++;
        }
        CHECK_INVARIANT(pick >= 0, "");
        res[pick] = (*cai)->getIdx();
      }
    }
    return res;
  }

  //
  // Determine bond wedge state
  ///
  Bond::BondDir DetermineBondWedgeState(const Bond *bond,
                                        const INT_MAP_INT &wedgeBonds, 
                                        const Conformer *conf){
    PRECONDITION(bond,"no bond");
    PRECONDITION(bond->getBondType()==Bond::SINGLE,
                 "bad bond order for wedging");
    const ROMol *mol=&(bond->getOwningMol());
    PRECONDITION(mol,"no mol");

    Bond::BondDir res=Bond::NONE;
    if(!conf){
      return res;
    }
    
    int bid = bond->getIdx();
    INT_MAP_INT_CI wbi = wedgeBonds.find(bid);
    if (wbi == wedgeBonds.end()) {
      return res;
    }

    unsigned int waid = wbi->second;
    
    Atom *atom, *bondAtom; // = bond->getBeginAtom();
    if (bond->getBeginAtom()->getIdx() == waid) {
      atom = bond->getBeginAtom();
      bondAtom = bond->getEndAtom();
    } else {
      atom = bond->getEndAtom();
      bondAtom = bond->getBeginAtom();
    }
      
    Atom::ChiralType chiralType=atom->getChiralTag();
    CHECK_INVARIANT( chiralType==Atom::CHI_TETRAHEDRAL_CW ||
                     chiralType==Atom::CHI_TETRAHEDRAL_CCW, "");
 
    // if we got this far, we really need to think about it:
    INT_LIST neighborBondIndices;
    DOUBLE_LIST neighborBondAngles;
    RDGeom::Point3D centerLoc, tmpPt;
    centerLoc=conf->getAtomPos(atom->getIdx());
    tmpPt=conf->getAtomPos(bondAtom->getIdx());
    centerLoc.z=0.0;
    tmpPt.z = 0.0;
    RDGeom::Point3D refVect=centerLoc.directionVector(tmpPt);

    neighborBondIndices.push_back(bond->getIdx());
    neighborBondAngles.push_back(0.0);

    bool hasTruePrecedingAtom=false;
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = mol->getAtomBonds(atom);
    ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap = mol->getBondPMap();
    while(beg!=end){
      Bond *nbrBond=pMap[*beg];
      Atom *otherAtom = nbrBond->getOtherAtom(atom);
      if(nbrBond != bond){
        tmpPt = conf->getAtomPos(otherAtom->getIdx());
        tmpPt.z = 0.0;
        RDGeom::Point3D tmpVect=centerLoc.directionVector(tmpPt);
        double angle=refVect.signedAngleTo(tmpVect);
        if(angle<0.0) angle += 2.*M_PI;
        INT_LIST::iterator nbrIt=neighborBondIndices.begin();
        DOUBLE_LIST::iterator angleIt=neighborBondAngles.begin();
        // find the location of this neighbor in our angle-sorted list
        // of neighbors:
        while(angleIt!=neighborBondAngles.end() && angle>(*angleIt)){
          angleIt++;
          nbrIt++;
        }
        neighborBondAngles.insert(angleIt,angle);
        neighborBondIndices.insert(nbrIt,nbrBond->getIdx());
      }
      unsigned int oIdx=otherAtom->getIdx();
      if(oIdx<atom->getIdx() && nbrBond->getBeginAtomIdx()==oIdx){
        hasTruePrecedingAtom=true;
      }
      beg++;
    }
    
    // at this point, neighborBondIndices contains a list of bond
    // indices from the central atom.  They are arranged starting
    // at the reference bond in CCW order (based on the current
    // depiction).
    int nSwaps = atom->getPerturbationOrder(neighborBondIndices);

    // in the case of three-coordinated atoms we may have to worry about
    // the location of the implicit hydrogen - Issue 209
    // Check if we have one of these situation
    //      
    //      0        1 0 2
    //      *         \*/
    //  1 - C - 2      C
    //     
    // here the hydrogen will be between 1 and 2 and we need to add an additional swap
    if (neighborBondAngles.size() == 3) {
      // three coordinated
      DOUBLE_LIST::iterator angleIt = neighborBondAngles.begin(); 
      angleIt++; // the first is the 0 (or reference bond - we will ignoire that
      double angle1 = (*angleIt);
      angleIt++;
      double angle2 = (*angleIt);
      if(angle2 - angle1 >= (M_PI-1e-4)){
        // we have the above situation
        nSwaps++;
      }
    }

#ifdef VERBOSE_STEREOCHEM
    BOOST_LOG(rdDebugLog) << "--------- " << nSwaps << std::endl;
    std::copy(neighborBondIndices.begin(),neighborBondIndices.end(),
              std::ostream_iterator<int>(BOOST_LOG(rdDebugLog)," "));
    BOOST_LOG(rdDebugLog) << std::endl;
    std::copy(neighborBondAngles.begin(),neighborBondAngles.end(),
              std::ostream_iterator<double>(BOOST_LOG(rdDebugLog)," "));
    BOOST_LOG(rdDebugLog) << std::endl;
#endif
    if(neighborBondIndices.size()==3 && !hasTruePrecedingAtom) ++nSwaps;
    if(chiralType==Atom::CHI_TETRAHEDRAL_CCW){
      if(nSwaps%2==1) {// ^ reverse) {
        res=Bond::BEGINDASH;
      } else {
        res=Bond::BEGINWEDGE;
      }
    } else {
      if (nSwaps%2==1) { // ^ reverse) {
        res=Bond::BEGINWEDGE;
      } else {
        res=Bond::BEGINDASH;
      }
    }

    return res;
  }



  // handles stereochem markers set by the Mol file parser and
  // converts them to the RD standard:
  void DetectAtomStereoChemistry(RWMol &mol, const Conformer *conf){
    PRECONDITION(conf,"no conformer");      

    // make sure we've calculated the implicit valence on each atom:
    for(RWMol::AtomIterator atomIt=mol.beginAtoms();
        atomIt!=mol.endAtoms();
        ++atomIt) {
      (*atomIt)->calcImplicitValence(false);
    }

    for(RWMol::BondIterator bondIt=mol.beginBonds();
        bondIt != mol.endBonds();
        ++bondIt){
      Bond *bond=*bondIt;
      if(bond->getBondDir() != Bond::UNKNOWN){
        Bond::BondDir dir=bond->getBondDir();
        // the bond is marked as chiral:
        if(dir==Bond::BEGINWEDGE || dir==Bond::BEGINDASH){
          Atom *atom = bond->getBeginAtom();
          if(atom->getImplicitValence()==-1){
            atom->calcExplicitValence();
            atom->calcImplicitValence(); 
          }
          Atom::ChiralType code=FindAtomStereochemistry(mol,bond, conf);
          atom->setChiralTag(code);
          // within the RD representation, if a three-coordinate atom
          // is chiral and has an implicit H, that H needs to be made explicit:
          if(atom->getDegree()==3 && !atom->getNumExplicitHs() &&
             atom->getNumImplicitHs()==1){
            atom->setNumExplicitHs(1);
            // recalculated number of implicit Hs:
            atom->updatePropertyCache();
          }
        }
        }
    }
    
    // now cleanup any bad specification we might have seen:
    MolOps::assignAtomChiralCodes(mol,true);
  }
  
  void DetectBondStereoChemistry(ROMol &mol, const Conformer *conf) {
    PRECONDITION(conf,"no conformer");      

    // mark all the non-ring double bonds that can be cis/trans:
    MolOps::findPotentialStereoBonds(mol);
    for (RWMol::BondIterator bondIt = mol.beginBonds();
         bondIt != mol.endBonds(); ++bondIt) {
      if ((*bondIt)->getBondType() == Bond::DOUBLE) {
        ComputeBondStereoChemistry(mol, *bondIt, conf);
      }
    }
    MolOps::assignBondStereoCodes(mol);
  }
}
