// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RWMol.h>
#include <math.h>
#include <GraphMol/MolOps.h>
#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include "EmbeddedFrag.h"
#include "DepictUtils.h"
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/Bond.h>
#include "RDDepictor.h"
#include <list>
#include <algorithm>

const double NEIGH_RADIUS = 2.5;

namespace RDDepict {

  EmbeddedFrag::EmbeddedFrag(unsigned int aid, const RDKit::ROMol *mol) {
    PRECONDITION(mol,"");
    PRECONDITION(aid < mol->getNumAtoms(), "");
    
    EmbeddedAtom eatm;
    eatm.aid = aid;
    RDGeom::Point2D org(0.0, 0.0);
    RDGeom::Point2D normal(1.0, 0.0);
    eatm.loc = org;
    eatm.normal = normal;
    eatm.angle = -1.0;
    eatm.ccw = true;
    eatm.neighs.clear();
    d_eatoms.clear();
    d_attachPts.clear();
    d_eatoms[aid] = eatm;
    d_done = false;
    dp_mol = mol;
    this->updateNewNeighs(aid); 
  }
    
  EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol, const RDKit::VECT_INT_VECT &fusedRings) {
    PRECONDITION(mol,"");
    dp_mol = mol;
    d_eatoms.clear();
    d_attachPts.clear();
    this->embedFusedRings(fusedRings);
    d_done = false;
  }

  
  EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol, const RDGeom::INT_POINT2D_MAP &coordMap) {
    // constructor of a case where the user specifies the coordinates for a portion of the 
    // atoms in the molecule - we will use these coordinates blindly without testing for any 
    // kind of correctness - user is GOD :)
   
    // we are not going to do much here simply add the atoms we have coordinates for to this fragment;
    // as a result this fragment may not be as ready to add new neighbors etc. for the following reason. 
    // - the user may have specified coords for only a part of the atoms in a fused ring systems
    // - once we use these coordinates we need to set up the atoms properly so that new 
    //   neighbors can be added to them
    PRECONDITION(mol,"");
    dp_mol = mol;
    d_eatoms.clear();
    d_attachPts.clear();
    RDGeom::INT_POINT2D_MAP_CI cri; 
    unsigned int na = mol->getNumAtoms();
    for (cri = coordMap.begin(); cri != coordMap.end(); cri++) {
      unsigned int aid = cri->first;
      CHECK_INVARIANT(aid < na, "");
      EmbeddedAtom eatom(aid, cri->second);
      eatom.neighs.clear();
      d_eatoms[aid] = eatom;
      d_done = false;
    }
    this->setupNewNeighs();
    
    // now for points that new atoms will be added to later on we need to do some setup
    //RDKit::INT_DEQUE_CI dai;
    RDKit::INT_LIST_CI dai;
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    for (dai = d_attachPts.begin(); dai != d_attachPts.end(); dai++) {
      boost::tie(nbrIdx,endNbrs) = dp_mol->getAtomNeighbors(dp_mol->getAtomWithIdx(*dai));
      // find the neighbors that are already embedded for each of these atoms
      RDKit::INT_VECT doneNbrs;
      const RDKit::INT_VECT &enbrs = d_eatoms[*dai].neighs;
      while (nbrIdx != endNbrs) {
        if (std::find(enbrs.begin(), enbrs.end(), static_cast<int>(*nbrIdx)) == enbrs.end()) {
          // we found a neighbor that is part of this embedded system
          doneNbrs.push_back(*nbrIdx);
        }
        nbrIdx++;
      }
      if (doneNbrs.size() == 0) {
        d_eatoms[*dai].normal = RDGeom::Point2D(1.,0.);
        d_eatoms[*dai].angle = -1.;
      } else if (doneNbrs.size() == 1) {
        int nbid = doneNbrs.front();
        d_eatoms[*dai].nbr1 = nbid; 
        d_eatoms[*dai].normal = computeNormal(d_eatoms[*dai].loc, 
                                              d_eatoms[nbid].loc);
      } else if (doneNbrs.size() == 2){
        int nb1 = doneNbrs[0]; 
        int nb2 = doneNbrs[1];
        d_eatoms[*dai].nbr1 = nb1; 
        d_eatoms[*dai].nbr2 = nb2;
        d_eatoms[*dai].angle = computeAngle(d_eatoms[*dai].loc, 
                                            d_eatoms[nb1].loc, d_eatoms[nb2].loc);
      } else if (doneNbrs.size() >= 3) {
        // this is a pain - delegate it to a utility function
        this->computeNbrsAndAng((*dai), doneNbrs);
      }
    }
  }
  
  int _anglComp(const DOUBLE_INT_PAIR &arg1, const DOUBLE_INT_PAIR &arg2) {
    return (arg1.first < arg2.first);
  }

  void EmbeddedFrag::computeNbrsAndAng(unsigned int aid, const RDKit::INT_VECT &doneNbrs) {
    //                                     const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol,"");
    PRECONDITION(aid<dp_mol->getNumAtoms(),"");
    
    PRECONDITION(doneNbrs.size() >= 3, "");
    // we will find all the inter nbr angles, pick the one with the largest angle 
    // make those neighbors the nbr1 and nbr2 of aid
    std::list<DOUBLE_INT_PAIR> anglePairs;
    RDKit::INT_VECT_CI nbi1, nbi2, nbi3;
    double ang;
    for (nbi1 = doneNbrs.begin(); nbi1 != doneNbrs.end(); nbi1++) {
      nbi3 = nbi1;
      for (nbi2 = nbi3++; nbi2 != doneNbrs.end(); nbi2++) {
        ang = computeAngle(d_eatoms[aid].loc, 
                           d_eatoms[*nbi1].loc, d_eatoms[*nbi2].loc);
        INT_PAIR nbrPair = std::make_pair((*nbi1), (*nbi2));
        anglePairs.push_back(std::make_pair(ang, nbrPair));
      }
    }     
    anglePairs.sort(_anglComp);
    std::list<DOUBLE_INT_PAIR>::reverse_iterator apcri;

    // more pain, more pain
    // we unfortunately cannot right away pick the largest angle - it is possible that
    // we pick an angle that is in a fused ring - see if I can explain this with a diagram
    //        _     _
    //       / B   C \                                this space
    //      /   \ /   \                               intentionally left blank
    //     |     A     |
    //     |     |     |
    //      \    D    /  
    //       \_/   \_/  
    //
    //  Let's say we are sitting on A with nbrs B, C, D - it is possible that we find
    //  ang(BAD) to be largest, but a new neighbor in this case will be added inside the ring
    //  We want to find ang(BAC) instead - which we will this do by checking that both our neighbors 
    //  are not involved in more than one ring. Bridged systems - don't even go there
    DOUBLE_INT_PAIR winner = anglePairs.back();
    for (apcri = anglePairs.rbegin(); apcri != anglePairs.rend(); apcri++) {
      INT_PAIR nbrPair = apcri->second;
      if ((dp_mol->getRingInfo()->numAtomRings(nbrPair.first) <= 1) &&
          (dp_mol->getRingInfo()->numAtomRings(nbrPair.second) <= 1)) {
        winner = (*apcri);
        break;
      }
    }

    INT_PAIR winPair = winner.second;
    int wnb1 = winPair.first;
    int wnb2 = winPair.second;
    
    // now find the smallest angle that contains one of these nbrs
    std::list<DOUBLE_INT_PAIR>::const_iterator apci;
    int nb2=-1, nb1=-1;
    for (apci = anglePairs.begin(); apci != anglePairs.end(); apci++) {
      INT_PAIR nbrPair = apci->second;
      if (wnb1 == nbrPair.first) {
        nb2 = wnb1;
        nb1 = nbrPair.second;
        break;
      } else if (wnb1 == nbrPair.second) {
        nb2 = wnb1;
        nb1 = nbrPair.first; 
        break;
      } else if (wnb2 == nbrPair.first) {
        nb2 = wnb2;
        nb1 = nbrPair.second;
        break;
      } else if (wnb2 == nbrPair.second) {
        nb2 = wnb2;
        nb1 = nbrPair.first; 
        break;
      }
    }

    // now find the rotation between nb1 and nb2
    double wAng = winner.first;
    d_eatoms[aid].rotDir = rotationDir(d_eatoms[aid].loc, d_eatoms[nb1].loc,
                                       d_eatoms[nb2].loc, wAng);
    d_eatoms[aid].nbr1 = nb1;
    d_eatoms[aid].nbr2 = nb2;
    d_eatoms[aid].angle = 2*M_PI - wAng;
  }

  // constructor to embed a cis/trans system
  EmbeddedFrag::EmbeddedFrag(const RDKit::Bond* dblBond) {
    // Earlier embedding a cis/trans system meant to assign coordinates to the
    // atoms on the double bond as well as the neighboring atoms connected by the 
    // single bond for which the cis/trans code has been specified. 
    // this causes some ugliness in cases where these neighboring atoms are either
    // part of a different cis/trans system or a ring system. The function "merge"
    // used to deal with this ugliness. 
    // Now we will just embed the atoms on the double bonds and mark at these atoms
    // the direction in which the in comming single bonds should go.
    // Makes the merge function easier and address issue 171 simultaneously.
    PRECONDITION(dblBond,"");
    PRECONDITION(dblBond->getBondType()==RDKit::Bond::DOUBLE, "");
    RDKit::Bond::BondStereo stype = dblBond->getStereo();
    PRECONDITION(stype > RDKit::Bond::STEREOANY, "");
    const RDKit::INT_VECT &nbrAtms = dblBond->getStereoAtoms();
    PRECONDITION(nbrAtms.size() == 2, "");
    dp_mol = &(dblBond->getOwningMol());

    int begAtm = dblBond->getBeginAtomIdx();
    int endAtm = dblBond->getEndAtomIdx();
    
    // the begin atom goes at the origin and the normal goes along -ve y-axis 
    // to be rotate clock to add the cis/trans single bond
    EmbeddedAtom beatm;
    beatm.aid = begAtm;
    beatm.loc = RDGeom::Point2D(0.0, 0.0);
    beatm.nbr1 = endAtm;
    
    beatm.normal = RDGeom::Point2D(0.0, -1.0);
    beatm.ccw = false;
    beatm.CisTransNbr = nbrAtms[0];
    d_eatoms[begAtm] = beatm;
    
    // the end atom goes on the x-axis 
    EmbeddedAtom eeatm;
    eeatm.aid = endAtm;
    eeatm.loc = RDGeom::Point2D(BOND_LEN, 0.0);
    eeatm.nbr1 = begAtm;
    eeatm.CisTransNbr = nbrAtms[1];
    if (stype == RDKit::Bond::STEREOZ) {
      eeatm.normal = RDGeom::Point2D(0.0, -1.0);
      eeatm.ccw = true;
    } else {
      eeatm.normal = RDGeom::Point2D(0.0, 1.0);
      eeatm.ccw = false;
    }
    d_eatoms[endAtm] = eeatm;

    d_done = false;
  }

  int EmbeddedFrag::findNumNeigh(const RDGeom::Point2D &pt, double radius) {
    // find the number of atoms in the current embedded system that are within 
    // 'radius' of the specified point
    INT_EATOM_MAP_CI efi;
    int res = 0;
    for (efi = d_eatoms.begin(); efi != d_eatoms.end(); efi++) {
      RDGeom::Point2D rloc = efi->second.loc;
      if ((rloc - pt).length() < radius) {
        res++;
      }
    }
    return res;
  }

  void EmbeddedFrag::updateNewNeighs(unsigned int aid) {//, const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol, "");
    
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    d_eatoms[aid].neighs.clear();
    boost::tie(nbrIdx,endNbrs) = dp_mol->getAtomNeighbors(dp_mol->getAtomWithIdx(aid));
    while (nbrIdx != endNbrs) {
      if (d_eatoms.find(*nbrIdx) == d_eatoms.end()) {
        d_eatoms[aid].neighs.push_back(*nbrIdx);
      }
      nbrIdx++;
    }

    int deg = dp_mol->getAtomWithIdx(aid)->getDegree();
    // order the neigbors by their CIPranks, if the number is between > 0 but less than 3
    if ((d_eatoms[aid].neighs.size() > 0) && 
        ((deg < 4) || (d_eatoms[aid].neighs.size() < 3) ) ) {
      d_eatoms[aid].neighs = rankAtomsByRank(*dp_mol, d_eatoms[aid].neighs);
    }
    else if ((deg >= 4) && (d_eatoms[aid].neighs.size() >= 3)){
      // now if we have more more than 2 neighbors change the order so that atoms with
      // the highest rank fall on opposite sides of each other
      d_eatoms[aid].neighs = setNbrOrder(aid, d_eatoms[aid].neighs, *dp_mol);
    }

    if (d_eatoms[aid].neighs.size() > 0) {
      if (std::find(d_attachPts.begin(), d_attachPts.end(),
                    static_cast<int>(aid)) == d_attachPts.end()) {
        d_attachPts.push_back(aid);
      }
    }
  }

  void EmbeddedFrag::setupNewNeighs() { //const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol, "");

    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    INT_EATOM_MAP_I eci;
    d_attachPts.clear();
    for (eci = d_eatoms.begin(); eci != d_eatoms.end(); eci++) {
      unsigned int aid = eci->first;
      this->updateNewNeighs(aid); 
    }
    // arrange the d_attachPts so that they are traversed in the order of CIPRanks
    d_attachPts = rankAtomsByRank(*dp_mol, d_attachPts);
    
  }
      
  int EmbeddedFrag::findNeighbor(unsigned int aid) { //, const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol, "");

    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    const RDKit::Atom *atm = dp_mol->getAtomWithIdx(aid);
    boost::tie(nbrIdx,endNbrs) = dp_mol->getAtomNeighbors(atm);
    while (nbrIdx != endNbrs) {
      if (d_eatoms.find(*nbrIdx) != d_eatoms.end()) {
        return (*nbrIdx);
      }
      nbrIdx++;
    }
    return -1;
  }

  //
  // NOTE: the individual rings in fusedRings must appear in traversal order.
  //    This is what is provided by the current ring-finding code.
  //
  void EmbeddedFrag::embedFusedRings(const RDKit::VECT_INT_VECT &fusedRings) {
    PRECONDITION(dp_mol,"");
    // ok this is what we are going to do here
    // embed each of the individual rings. Then 
    // find the largest ring , leave that at the origin 
    // and fuse each of remaining rings

    // get the union of the atoms in the rings
    RDKit::INT_VECT funion; 
    RDKit::Union(fusedRings, funion);

    // embed each of the rings independenty and find the largest ring
    std::vector<RDGeom::INT_POINT2D_MAP> coords;
    coords.reserve(fusedRings.size());
    RDKit::VECT_INT_VECT_CI ri;
    // FIX for issue 197
    // find the ring with the max substituents
    // If there are multiple pick the largest
    int firstRingId = pickFirstRingToEmbed(*dp_mol, fusedRings);
    
    for (ri = fusedRings.begin(); ri != fusedRings.end(); ri++) {
      coords.push_back(embedRing((*ri)));
    }
    
    this->initFromRingCoords(fusedRings[firstRingId], coords[firstRingId]);
    
    RDKit::INT_VECT doneRings;
    doneRings.push_back(firstRingId);

    // now loop over the remaining rings and attach then one at a time
    // the order is determined by how many atoms a ring has in common with
    // the atoms already embedded
    while (d_eatoms.size() < funion.size()) { // ) {
      int nextId;
      // we will take the ring with maximum number of common atoms with 
      // with atoms already done
      RDKit::INT_VECT commonAtomIds = findNextRingToEmbed(doneRings, fusedRings, nextId);

      RDGeom::Transform2D trans;
      EmbeddedFrag embRing;
      embRing.initFromRingCoords(fusedRings[nextId], coords[nextId]);
      RDKit::INT_VECT pinAtoms;
      // REVIEW: using the average position of the shared atoms and the
      // centroid vector, we can make this a single case.
      if (commonAtomIds.size() == 1) {
        trans.assign(this->computeOneAtomTrans(commonAtomIds[0], embRing));
        embRing.Transform(trans);
        pinAtoms.push_back(commonAtomIds.front());
      }
      else {
        // if the common atoms form a chain they are going to be in order - we try to 
        // do that in findNextRingToEmbed 
        // we will therfore try to use the last and the first atoms in the chain to 
        // fuse the rings - will hopefully fix issue 177
        int aid1 = commonAtomIds.front();
        int aid2 = commonAtomIds.back();
        pinAtoms.push_back(aid1);
        pinAtoms.push_back(aid2);
        trans.assign(this->computeTwoAtomTrans(aid1, aid2, coords[nextId]));
        embRing.Transform(trans);
        reflectIfNecessaryDensity(embRing, aid1, aid2);
      }
      
      this->mergeRing(embRing, commonAtomIds.size(), pinAtoms);
      doneRings.push_back(nextId);
    }
    
  }
  
  RDGeom::Transform2D EmbeddedFrag::computeOneAtomTrans(unsigned int commAid, 
                                                        const EmbeddedFrag &other) {
      
    // find the coordinates for the same atom in the embedded system
    RDGeom::Point2D rcr = d_eatoms[commAid].loc;
    
    // find the coordinate for the same atom in the other system
    //INT_EATOM_MAP_CI eati = other.d_eatoms.find(commAid);
    const EmbeddedAtom &oeatm = other.GetEmbeddedAtom(commAid);
    RDGeom::Point2D ccr = oeatm.loc;
    int onb1 = oeatm.nbr1;
    int onb2 = oeatm.nbr2; 
    CHECK_INVARIANT ((onb1>=0) && (onb2>=0), "");
    RDGeom::Point2D midPt = other.GetEmbeddedAtom(onb1).loc;
    midPt += other.GetEmbeddedAtom(onb2).loc;
    midPt *= 0.5;

    // get the coordinates for the neighboring atoms
    int nb1 = d_eatoms[commAid].nbr1;
    int nb2 = d_eatoms[commAid].nbr2;
    RDGeom::Point2D nbp1 = d_eatoms[nb1].loc;
    RDGeom::Point2D nbp2 = d_eatoms[nb2].loc;

    double ang = d_eatoms[commAid].angle;
    double largestAngle = 2*M_PI - ang;
  
    RDGeom::Point2D bpt = computeBisectPoint(rcr, largestAngle, nbp1, nbp2);

    // now that we have the bisect point compute the transform that will take ccr to coincide with rcr 
    // and the mid point between the neighbors of ccr to fall on the line from rcr to bpt
    RDGeom::Transform2D trans;
    trans.SetTransform(rcr, bpt, ccr, midPt);
    return trans;
  }

  RDGeom::Transform2D EmbeddedFrag::computeTwoAtomTrans(unsigned int aid1, unsigned int aid2, 
                                                        const RDGeom::INT_POINT2D_MAP &nringCor) {
    // this is an easier thing to do than computeOneAtomTrans
    // we know that there are atleast two atoms in common between the new ring and the
    // rings that have already been embedded.
    // 
    // we are going to simply use the first two atoms on the commIds list and
    // use those to compute a transforms
    RDGeom::INT_POINT2D_MAP_CI posi;
    posi = nringCor.find(aid1);
    RDGeom::Point2D loc1 = posi->second;
    posi = nringCor.find(aid2);
    RDGeom::Point2D loc2 = posi->second;

    // get the coordinates for the same atoms in the already embedded ring system
    CHECK_INVARIANT(d_eatoms.find(aid1) != d_eatoms.end(), "");
    CHECK_INVARIANT(d_eatoms.find(aid2) != d_eatoms.end(), "");
    RDGeom::Point2D ref1 = d_eatoms[aid1].loc;
    RDGeom::Point2D ref2 = d_eatoms[aid2].loc;
    RDGeom::Transform2D trans;
    trans.SetTransform(ref1, ref2, loc1, loc2);
    return trans;
  }

  void EmbeddedFrag::Reflect(const RDGeom::Point2D &loc1,
                             const RDGeom::Point2D &loc2) {
    INT_EATOM_MAP_I ei;
    RDGeom::Point2D temp;
    for (ei = d_eatoms.begin(); ei != d_eatoms.end(); ei++) {
      ei->second.Reflect(loc1, loc2);
    }
  }

  void EmbeddedFrag::reflectIfNecessaryCisTrans(EmbeddedFrag &embFrag, unsigned int ctCase, 
                                                unsigned int aid1, unsigned int aid2) {
    // ok this is a cis/trans case - we may have violated the cis/trans specification
    // so lets try to correct it with a reflection
    double dot;
    RDGeom::Point2D p1Loc = d_eatoms[aid1].loc;
    int ringAtm;
    RDGeom::Point2D p1norm, rAtmLoc;
    if (ctCase == 1) {
      // embObj is the cis/trans case - find the normal at aid1 - this shoujd tell us
      // where the ring single bond in the cis/trans system should have gone
      p1norm = embFrag.d_eatoms[aid1].normal;
      ringAtm = embFrag.d_eatoms[aid1].CisTransNbr;
      if(d_eatoms.find(ringAtm)!=d_eatoms.end()){
        rAtmLoc = d_eatoms[ringAtm].loc;
      } else {
        // FIX: this is a work-around arising from issue 3135833
        BOOST_LOG(rdWarningLog)<<"Warning: stereochemistry around double bond may be incorrect in depiction."<<std::endl;
        return;
      }
    } else {
      // this is the cis/trans object
      p1norm = d_eatoms[aid1].normal;
      ringAtm = d_eatoms[aid1].CisTransNbr;
      rAtmLoc = embFrag.d_eatoms[ringAtm].loc;
    }
    rAtmLoc -= p1Loc;
    dot = rAtmLoc.dotProduct(p1norm);
    RDGeom::Point2D p2Loc = d_eatoms[aid2].loc;
    if (dot < 0.0) {
      embFrag.Reflect(p1Loc, p2Loc);
    }
  }

  void EmbeddedFrag::reflectIfNecessaryThirdPt(EmbeddedFrag &embFrag, unsigned int aid1, 
                                               unsigned int aid2, unsigned int aid3) {
    RDGeom::Point2D oth3 = embFrag.GetEmbeddedAtom(aid3).loc;
    RDGeom::Point2D pt3 = d_eatoms[aid3].loc;
    RDGeom::Point2D pt1 = d_eatoms[aid1].loc;
    RDGeom::Point2D pt2 = d_eatoms[aid2].loc;
    
    RDGeom::Point2D normal = pt2;
    normal -= pt1;
    normal.rotate90();
        
    pt3 -= pt1;
    oth3 -= pt1;

    double dot1 = normal.dotProduct(pt3);
    double dot2 = normal.dotProduct(oth3);
    if (dot1*dot2 < 0.0) {
      // the third atom is on either sides of the line between aid1 and aid2 in the 
      // two fragment - let us reflect to correct it
      embFrag.Reflect(pt1, pt2);
    }
  }

  void EmbeddedFrag::reflectIfNecessaryDensity(EmbeddedFrag &embFrag, 
                                               unsigned int aid1, unsigned int aid2) {
    // ok we will do this the new way by measuring a density function
    RDGeom::Point2D pin1 = d_eatoms[aid1].loc;
    RDGeom::Point2D pin2 = d_eatoms[aid2].loc;
    double densityNormal = 0.0;
    double densityReflect = 0.0;
    INT_EATOM_MAP_CI oci, tci;
    int oaid;
    RDGeom::Point2D loc1, rloc1, loc2, t1, rt1;
    const INT_EATOM_MAP &oatoms = embFrag.GetEmbeddedAtoms();
    for (oci = oatoms.begin(); oci != oatoms.end(); oci++) {
      oaid = oci->first;
      if (d_eatoms.find(oaid) == d_eatoms.end()) {
        loc1 = oci->second.loc;
        rloc1 = reflectPoint(loc1, pin1, pin2);
        for (tci = d_eatoms.begin(); tci != d_eatoms.end(); tci++) {
          t1 = tci->second.loc;
          t1 -= loc1;
          double td = t1.length();
          rt1 = tci->second.loc;
          rt1 -= rloc1;
          double rtd = rt1.length();
          if (td > 1.0e-3) {
            densityNormal += (1.0/td);
          } else {
            densityNormal += 1000.0;
          }
          if (rtd > 1.0e-3) {
            densityReflect += (1.0/rtd);
          } else {
            densityReflect += 1000.0;
          }
        }
      }
    }
    if (densityNormal - densityReflect > 1.0e-4) {
      embFrag.Reflect(pin1, pin2);
    }
  }

  void EmbeddedFrag::initFromRingCoords(const RDKit::INT_VECT &ring, 
                                        const RDGeom::INT_POINT2D_MAP &nringMap) {
    double largestAngle =  M_PI*( 1 - (2.0/ring.size()));
    RDKit::INT_VECT_CI ai;
    int prev = ring.back();
    int cnt = 0;
    RDGeom::INT_POINT2D_MAP_CI coord;
    for (ai = ring.begin(); ai != ring.end(); ai++) {
      EmbeddedAtom eatm;
      // this sucks - the following find is because of the constness of nringMap
      // nringMap[*ai] will not work
      coord = nringMap.find(*ai);
      eatm.loc = coord->second;
      eatm.aid = (*ai);
      eatm.angle = largestAngle;
      eatm.nbr1 = prev;
      if (cnt > 0) {
        d_eatoms[prev].nbr2 = (*ai);
      }
      d_eatoms[(*ai)] = eatm;
      prev = (*ai);
      cnt++;
    }
    d_eatoms[prev].nbr2 = ring.front();
  }

  void EmbeddedFrag::mergeRing(const EmbeddedFrag &embRing, unsigned int nCommon,
                               const RDKit::INT_VECT &pinAtoms) {
    const INT_EATOM_MAP &oatoms = embRing.GetEmbeddedAtoms();
    INT_EATOM_MAP_CI ori;
    for (ori = oatoms.begin(); ori != oatoms.end(); ori++) {
      int aid = ori->first;
      if (d_eatoms.find(aid) == d_eatoms.end()) {
        d_eatoms[aid] = ori->second;
      }
      else {
        // update the neighbor only on atoms that were used to compute the transform to merge the
        // and only if the the two are the only common atoms 
        // i.e. we are doing bridged systems we will leave the nbrs untouched
        if (nCommon <= 2) { 
          if (std::find(pinAtoms.begin(), pinAtoms.end(), aid) != pinAtoms.end()) {
            d_eatoms[aid].angle += ori->second.angle;
            if (d_eatoms[aid].nbr1 == ori->second.nbr1) {
              d_eatoms[aid].nbr1 = ori->second.nbr2;
            } else if (d_eatoms[aid].nbr1 == ori->second.nbr2) {
              d_eatoms[aid].nbr1 = ori->second.nbr1;
            } else if (d_eatoms[aid].nbr2 == ori->second.nbr1) {
              d_eatoms[aid].nbr2 = ori->second.nbr2;
            } else if (d_eatoms[aid].nbr2 == ori->second.nbr2) {
              d_eatoms[aid].nbr2 = ori->second.nbr1;
            }
          }
        }
      }
    }
  }

  void EmbeddedFrag::addNonRingAtom(unsigned int aid, unsigned int toAid) {
    //const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol, "");
    // check that aid does not belong the the embedded fragment yet
    PRECONDITION(d_eatoms.find(aid) == d_eatoms.end(), "");
    // and that toAid is already in the embedded system
    PRECONDITION(d_eatoms.find(toAid) != d_eatoms.end(), "");
    if (d_eatoms[toAid].angle > 0.0) {
      addAtomToAtomWithAng(aid, toAid);
    } else {
      addAtomToAtomWithNoAng(aid, toAid); //, mol);
    }
    // remove aid from the neighbor list of toAid
    d_eatoms[toAid].neighs.erase(std::remove(d_eatoms[toAid].neighs.begin(),
                                             d_eatoms[toAid].neighs.end(),
                                             static_cast<int>(aid)));
    this->updateNewNeighs(aid); //, mol);
  }

  void EmbeddedFrag::addAtomToAtomWithAng(unsigned int aid, unsigned int toAid) {
    
    EmbeddedAtom refAtom = d_eatoms[toAid];
    RDGeom::Point2D refLoc = refAtom.loc;
    RDGeom::Point2D origin(0.0, 0.0);
    PRECONDITION(refAtom.angle > 0.0, "");

    // we are adding to either to a ring atom or an atom to which we added atleast one 
    //substituent previously
    
    // determine the angle at which we want to add the new atom based on the number 
    // of remaining substituents
    int nnbr = refAtom.neighs.size();
    double remAngle = 2*M_PI - refAtom.angle;
    double currAngle = remAngle/(1 + nnbr);
    d_eatoms[toAid].angle += currAngle;
    
    RDGeom::Point2D nb1 = d_eatoms[refAtom.nbr1].loc;
    RDGeom::Point2D nb2 = d_eatoms[refAtom.nbr2].loc;
    RDGeom::Point2D rotar;
    if (d_eatoms[toAid].rotDir == 0) {
      d_eatoms[toAid].rotDir = rotationDir(refLoc, nb1, nb2, remAngle);
    }
    
    currAngle *= d_eatoms[toAid].rotDir;
    
    RDGeom::Transform2D rtrans;
    rtrans.SetTransform(refLoc, currAngle);
    RDGeom::Point2D currLoc = nb2;
    rtrans.TransformPoint(currLoc);
    
    // set the neighbors for the current point
    d_eatoms[toAid].nbr2 = aid;
    
    EmbeddedAtom eatm;
    eatm.aid = aid;
    eatm.loc = currLoc;
    eatm.nbr1 = toAid;
    eatm.angle = -1.0;
    // now compute the normal at this atom - which gives the direction in which we want to
    // add the next atom. We will go in the direction that seem to be least explored 
    RDGeom::Point2D tpt = currLoc - refLoc;
    RDGeom::Point2D norm, tp1, tp2;
    norm.x = -tpt.y;
    norm.y = tpt.x;
    tp1 = currLoc + norm;
    tp2 = currLoc - norm;
    
    int nccw = findNumNeigh(tp1, NEIGH_RADIUS); // number of neighbors if we go counter-clockwise
    int ncw = findNumNeigh(tp2, NEIGH_RADIUS); // number of neighbors if we go clockwise
    
    norm.normalize();
    if (nccw < ncw) {
      eatm.normal = norm;
      eatm.ccw = false;
    } else {
      eatm.normal = (-norm);
      eatm.ccw = true;
    }
    
    d_eatoms[aid] = eatm;
  }

  void EmbeddedFrag::addAtomToAtomWithNoAng(unsigned int aid, unsigned int toAid) {
    //const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol,"");
    EmbeddedAtom refAtom = d_eatoms[toAid];
    PRECONDITION(refAtom.angle <= 0.0, "");
    RDGeom::Point2D refLoc = refAtom.loc;
    RDGeom::Point2D origin(0.0, 0.0);

    // -----------------------------------------------------------------------
    // we are adding to a non-ring atom,
    // the direction in which we add the new atom matters here
    RDGeom::Point2D currLoc = refAtom.normal;
    if (refAtom.CisTransNbr >= 0) {
      // ok this atom is part of a cis/trans dbl bond 
      if (static_cast<unsigned int>(refAtom.CisTransNbr) != aid) {
        // but we are note adding the single bond atom to which the cis/trans specification was
        // made, inthis case reverse the normal and the ccw
        refAtom.ccw = !(refAtom.ccw);
        currLoc *= -1.0;  
      }
    }
    
    CHECK_INVARIANT(currLoc.lengthSq() > 1.0e-8, "");
    
    // find out what angle we want to add bond at
    const RDKit::Atom *atm = dp_mol->getAtomWithIdx(toAid);
    int deg = atm->getDegree();
    
    double angle = computeSubAngle(deg, atm->getHybridization()); 
    
    // update the current atom 
    // we already have a nbr1 set on the current atom update the angle etc
    //d_eatoms[toAid].nbr2 = aid;
    bool flipNorm = false;
    if (d_eatoms[toAid].nbr1 >= 0) {
      d_eatoms[toAid].angle = angle;
      d_eatoms[toAid].nbr2 = aid;
    }
    else {
      // ------------------
      // We'll be here for the first atom in a system with no rings, we have nothing
      // else set up, so we will deal with this case carefully.
      //  - if the angle is 120 deg we will add the first atom at 30 deg angle to the x-axis 
      //  - for any other angle we will use the x-axis to add the new atom
      //  - we will set the normal perpendicular to this first bond in teh counter clockwis direction
      //
      //RDGeom::Point2D norm;
      
      RDGeom::Point2D norm = d_eatoms[toAid].normal;
      double ang = angle;
      
      RDGeom::Transform2D rtrans;
      rtrans.SetTransform(origin, ang);
      rtrans.TransformPoint(norm);
      d_eatoms[toAid].normal = norm;
      d_eatoms[toAid].nbr1 = aid;
      flipNorm = true;
      
    }
    
    angle -= M_PI/2;
    if (!refAtom.ccw) {
      // we want to rotate cloackwise
      angle *= -1.0;
    }
    
    RDGeom::Transform2D trans;
    trans.SetTransform(origin, angle);
    trans.TransformPoint(currLoc);
    currLoc *= BOND_LEN;
    currLoc += refLoc;
    
    // now compute the normal at this new point for the next addition
    RDGeom::Point2D tpt = refLoc - currLoc;
    RDGeom::Point2D norm;
    // This is the lazy man's rotation by 90 degrees about the origin:
    norm.x = -tpt.y;
    norm.y = tpt.x;
    if ((refAtom.ccw) ^ (flipNorm)){
      norm *= -1.0;
    }
    norm.normalize();
    EmbeddedAtom eatm;
    eatm.loc = currLoc;
    eatm.normal = norm;
    eatm.nbr1 = toAid;
    
    eatm.angle = -1.0;
    
    eatm.ccw = (!refAtom.ccw)^(flipNorm);
    d_eatoms[aid] = eatm;
  }
  
  RDKit::INT_VECT EmbeddedFrag::findCommonAtoms(const EmbeddedFrag &efrag2) {
    RDKit::INT_VECT res;
    const INT_EATOM_MAP &eatms1 = this->GetEmbeddedAtoms();
    const INT_EATOM_MAP &eatms2 = efrag2.GetEmbeddedAtoms();

    INT_EATOM_MAP_CI eri1, eri2;
    for (eri1 = eatms1.begin(); eri1 != eatms1.end(); eri1++) {
      for (eri2 = eatms2.begin(); eri2 != eatms2.end(); eri2++) {
        if (eri1->first == eri2->first) {
          res.push_back(eri1->first);
        }
      }
    }
    return res;
  }
  
  void EmbeddedFrag::mergeNoCommon(EmbeddedFrag &embObj, unsigned int toAid, unsigned int nbrAid){
    // merge embObj to this fragment when there are no common atoms between the two fragments
    PRECONDITION(dp_mol, "");
    // check that both this fragment and the one we are merging with belong to the same molecule
    PRECONDITION(dp_mol == embObj.getMol(), "Molecule mismatch");
    RDKit::INT_VECT commAtms;
    this->addNonRingAtom(nbrAid, toAid); //, mol);
    embObj.addNonRingAtom(toAid, nbrAid); //, mol);
    commAtms.push_back(toAid);
    commAtms.push_back(nbrAid);
    this->mergeWithCommon(embObj, commAtms); //, mol);
  }

  void EmbeddedFrag::mergeWithCommon(EmbeddedFrag &embObj, RDKit::INT_VECT &commAtms) {
    PRECONDITION(dp_mol, "");
    PRECONDITION(dp_mol == embObj.getMol(), "Molecule mismatch");
    PRECONDITION(commAtms.size() >= 1, "");
    

    // we already have one or more common atoms between this fragment  
    // One atom in common can happen (look at issue 173)
    // - for cases where a cis/trans double bond is being merged with a 
    //   ring system that shares one of atoms on the double bond.
    // - or if 'this' fragment was created by user specified coordinates - where only
    //    part of a fused ring system or cis/trans system was specified
    
    // if we have one atom in common, we have to deal with it carefully - 
    unsigned int ctCase = 0; // book-keeper - if we have to merge a ring with a cis/trans dbl bond 
    // what kind is it
    // 0 - if we are doing a cis/trans merge, 1 - cis/trans and embObj is the
    //  dblBond, 2 - cis/trans merge and 'this' is the dblBond
    
    if (commAtms.size() == 1) {
      // couple of possibilities here
      // 1. we are merging a ring system with a cis/trans dbl bond
      // 2. We are merging with a fused ring system out of which one of the atoms
      //    has already been embedded beacause the user specified its coordinates
      // First deal with the cis/trans case
      int commAid = commAtms.front();
      int otherAtom = -1;
      if (d_eatoms[commAid].CisTransNbr >= 0) {
        ctCase = 2;
        // this fragment is the cis/trans dbl bond
        otherAtom = d_eatoms[commAid].nbr1; // this is the other atom on the double bnd
        // now add this atom to the other fragment
        embObj.addNonRingAtom(otherAtom, commAid); //, mol);
      } else if (embObj.d_eatoms[commAid].CisTransNbr >= 0) {
        ctCase = 1;
        // otherwise embObj is the cis/trans dbl bond
        otherAtom = embObj.d_eatoms[commAid].nbr1;
        this->addNonRingAtom(otherAtom, commAid); //, mol);
      } else {
        otherAtom = d_eatoms[commAid].nbr1;
        if(otherAtom>=0){
          embObj.addNonRingAtom(otherAtom, commAid); //, mol);
        }
      }
      if (otherAtom >= 0) {
        commAtms.push_back(otherAtom);
      }
    }

    RDGeom::Transform2D rtrans;
    if (commAtms.size() == 1) {
      // if we have only one atom in common we will use a one atom transform 
      rtrans.assign(this->computeOneAtomTrans(commAtms.front(), embObj));
    } else { 
      // if we have more than one we will use a two point transform
      RDGeom::Point2D ref1, ref2, oth1, oth2;
      int cid1 = commAtms[0]; 
      int cid2 = commAtms[1];
      ref1 = d_eatoms[cid1].loc;
      ref2 = d_eatoms[cid2].loc;
      oth1 = embObj.GetEmbeddedAtom(cid1).loc;
      oth2 = embObj.GetEmbeddedAtom(cid2).loc;
      // now compute the transform
      rtrans.SetTransform(ref1, ref2, oth1, oth2);
    }

    // transform the second fragment
    embObj.Transform(rtrans);

    // check to see if this transform screws up any cis/trans specifications
    if (commAtms.size() >= 2) {
      if (ctCase > 0) {
        // we have a cis/trans case we may have violated the specification
        // check and correct it with a reflection
        reflectIfNecessaryCisTrans(embObj, ctCase, commAtms[0], commAtms[1]);
      } else if (commAtms.size() == 2) {
        // we have just two atoms in common but we may a simply overcrowed one side
        // check for crowding and reflect
        reflectIfNecessaryDensity(embObj, commAtms[0], commAtms[1]);
      } else {
        // finally if we have more than two atoms in common - we will use the third
        // atom to figure out if we need a reflection12
        reflectIfNecessaryThirdPt(embObj, commAtms[0], commAtms[1], commAtms[2]);
      }
    }

    // finally merge the fragment by copying the non common atoms
    const INT_EATOM_MAP &oatoms = embObj.GetEmbeddedAtoms();
    INT_EATOM_MAP_CI ori;
    // copy the eatoms in embObj to this fragment 
    for (ori = oatoms.begin(); ori != oatoms.end(); ori++) {
      int aid = ori->first;
      if (std::find(commAtms.begin(), commAtms.end(), aid) == commAtms.end()) { 
        d_eatoms[aid] = ori->second;
        // also if any of these atoms have unattached neighbors add them to the queue
        if (ori->second.neighs.size() > 0) {
          if (std::find(d_attachPts.begin(), d_attachPts.end(), aid) == d_attachPts.end()) {
            d_attachPts.push_back(aid);
          }
        }
      } else {
        if (ori->second.CisTransNbr >= 0) {
          d_eatoms[aid].CisTransNbr = ori->second.CisTransNbr;
          d_eatoms[aid].normal = ori->second.normal;
          d_eatoms[aid].ccw = ori->second.ccw;
        }
        if (ori->second.angle > 0.0) {
          d_eatoms[aid].angle = ori->second.angle;
          d_eatoms[aid].nbr1 = ori->second.nbr1;
          d_eatoms[aid].nbr2 = ori->second.nbr2;
        }
      }
    }

    // remember to update the not yet done neighbor of nbrAid
    RDKit::INT_VECT_CI cai;
    for (cai = commAtms.begin(); cai != commAtms.end(); cai++) {
      this->updateNewNeighs((*cai)); //, mol);
    }
  }

  void EmbeddedFrag::mergeFragsWithComm(std::list<EmbeddedFrag> &efrags) { //, const RDKit::ROMol *mol) {
    PRECONDITION(dp_mol,"");
    // first merge any fragments what share atoms in common
    std::list<EmbeddedFrag>::iterator efri, nfri;
    while (1) {
      RDKit::INT_VECT commAtms;
      for (efri = efrags.begin(); efri != efrags.end(); efri++) {
        if (!efri->isDone()) {
          commAtms = this->findCommonAtoms(*efri);
          if (commAtms.size() > 0) {
            nfri = efri;
            break;
          }
        }
      }
      if (commAtms.size() == 0) {
        break;
      }

      this->mergeWithCommon((*nfri), commAtms); //, mol);
      RDKit::INT_VECT_CI cai;
      for (cai = commAtms.begin(); cai != commAtms.end(); cai++) {
        if ((d_eatoms[*cai].neighs.size() == 0) &&
            (std::find(d_attachPts.begin(), d_attachPts.end(), (*cai)) != d_attachPts.end()) ) {
          
          d_attachPts.erase(std::remove(d_attachPts.begin(), d_attachPts.end(), (*cai)));
        }
      }
      efrags.erase(nfri);
    }
  }

  void EmbeddedFrag::expandEfrag(RDKit::INT_LIST &nratms, std::list<EmbeddedFrag> &efrags) {
                                
    PRECONDITION(dp_mol, "");
    
    // first merge any fragments that share atoms in common
    std::list<EmbeddedFrag>::iterator efri, nfri;

    this->mergeFragsWithComm(efrags); //, dp_mol);

    while (d_attachPts.size() > 0) {
      int aid = d_attachPts.front();
      RDKit::INT_VECT nbrs = d_eatoms[aid].neighs;
      CHECK_INVARIANT(nbrs.size() > 0, "");
      RDKit::INT_VECT_I nbri;
      RDKit::INT_LIST_I nratmi;
      for (nbri = nbrs.begin(); nbri != nbrs.end(); nbri++) {
        nratmi = std::find(nratms.begin(), nratms.end(), (*nbri));
        if (nratmi != nratms.end()) {
          // the neighbor we have to add is a non ring atoms
          this->addNonRingAtom((*nbri), aid); //, mol);
          // remove this atom we just added from the nnratms list
          nratms.erase(nratmi);
        }
        else {
          // the neighbor atom must be part of a different embedded fragment - 
          // merge that fragment with this one
          nfri = efrags.end();
          for (efri = efrags.begin(); efri != efrags.end(); efri++) {
            // don't search fragments that are done
            if (!efri->isDone()) {
              const INT_EATOM_MAP &eatoms = efri->GetEmbeddedAtoms();
              if (eatoms.find(*nbri) != eatoms.end()) {
                nfri = efri;
                break;
              }
            }
          }
          if (nfri != efrags.end()) {
            this->mergeNoCommon((*nfri), aid, (*nbri)); //, mol);
            if ((d_eatoms[*nbri].neighs.size() == 0) &&
                (std::find(d_attachPts.begin(), d_attachPts.end(), (*nbri)) != d_attachPts.end()) ) {
              d_attachPts.erase(std::remove(d_attachPts.begin(), d_attachPts.end(), (*nbri)));
            }
            // remove this fragment from the list of embedded fragments
            efrags.erase(nfri);
          }
        }
      }
      
      // ok we are done with this atom forever
      d_attachPts.pop_front();
      d_eatoms[aid].neighs.clear();
      // now that we added new atoms to the this fragments - check if there are new
      // fragment we have common atoms with and merge with them
      this->mergeFragsWithComm(efrags); //, mol);
    }
  }

  void EmbeddedFrag::Transform(const RDGeom::Transform2D &trans) {
    INT_EATOM_MAP_I eri;
    for (eri = d_eatoms.begin(); eri != d_eatoms.end(); eri++) {
      eri->second.Transform(trans);
    }
  }

  void EmbeddedFrag::computeBox() {
    INT_EATOM_MAP_I eri;
    d_px = -1.0e8;
    d_nx = 1.0e8;
    d_py = -1.0e8;
    d_ny = 1.0e8;

    for (eri = d_eatoms.begin(); eri != d_eatoms.end(); eri++) {
      const RDGeom::Point2D &loc = eri->second.loc;
      if (loc.x > d_px) {
        d_px = loc.x;
      } 
      if (loc.x < d_nx) {
        d_nx = loc.x;
      }
      if (loc.y > d_py) {
        d_py = loc.y;
      }
      if (loc.y < d_ny) {
        d_ny = loc.y;
      }
    }
    d_nx *= -1.0;
    d_ny *= -1.0;
  }

  void EmbeddedFrag::canonicalizeOrientation() {
    // fix for issue 198
    // no need to canonicalize if we are dealing with a single atm
    if (d_eatoms.size() <= 1) return;
    
    RDGeom::Point2D cent(0.0, 0.0);
    INT_EATOM_MAP_I eri;
    for (eri = d_eatoms.begin(); eri != d_eatoms.end(); eri++) {
      cent += eri->second.loc;
    }
    cent *= (1.0/d_eatoms.size());
      
    double xx, xy, yy;
    xx = 0.0; xy = 0.0; yy = 0.0; 
      
    // shift the center of the fragment to the origin and compute the covariance matrix
    for (eri = d_eatoms.begin(); eri != d_eatoms.end(); eri++) {
      eri->second.loc -= cent;
      xx += (eri->second.loc.x)*(eri->second.loc.x);
      xy += (eri->second.loc.x)*(eri->second.loc.y);
      yy += (eri->second.loc.y)*(eri->second.loc.y);
    }
      
    RDGeom::Point2D eig1, eig2;
    // the eigen vectors are given by (2*xy, (yy - xx) + d) and (2*xy, (yy - xx) - d) 
    // where d = sqrt((xx - yy)^2 + 4*xy^2)
    double d = (xx - yy)*(xx - yy) + 4*xy*xy;
    d = sqrt(d);
    RDGeom::Transform2D trans;
    eig1.x = 2*xy;
    eig1.y = (yy - xx) + d;
    if(eig1.length()<=1e-4) return;
    double eVal1=(xx+yy+d)/2;
    eig1.normalize();
    
    eig2.x = 2*xy;
    eig2.y = (yy - xx) - d;
    if(eig2.length()<=1e-4) return;
    double eVal2=(xx+yy-d)/2;
    eig2.normalize();

    // make sure eig1 corresponds to the larger eigenvalue:
    if(eVal2>eVal1){
      RDGeom::Point2D tmp=eig1;
      eig1=eig2;
      eig2=tmp;
    }
    // now rotate eig1 onto the X axis:
    trans.setVal(0,0, eig1.x);
    trans.setVal(1,0, -eig1.y);
    trans.setVal(0,1, eig1.y);
    trans.setVal(1,1, eig1.x);
    this->Transform(trans);
  }

  void _recurseAtomOneSide(unsigned int endAid, unsigned int begAid, const RDKit::ROMol *mol,
                           RDKit::INT_VECT &flipAids) {
    PRECONDITION(mol,"");
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    flipAids.push_back(endAid);
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(mol->getAtomWithIdx(endAid));
    while (nbrIdx != endNbrs) {
      if (((*nbrIdx) != begAid) && 
          (std::find(flipAids.begin(), flipAids.end(), static_cast<int>(*nbrIdx)) == flipAids.end()) ) {
        _recurseAtomOneSide(*nbrIdx, begAid, mol, flipAids);
      }
      nbrIdx++;
    }
    return;
  }

  double _crossVal(const RDGeom::Point2D &v1, const RDGeom::Point2D &v2) {
    double res = (v1.x)*(v2.y) - (v2.x)*(v1.y);
    return res;
  }

  int _pairDIICompAscending(const  PAIR_D_I_I &arg1, const PAIR_D_I_I &arg2) {
    return (arg1.first < arg2.first);
  }

  PAIR_I_I _findClosestPair(unsigned int beg1, unsigned int end1, 
                            unsigned int beg2, unsigned int end2, 
                            const RDKit::ROMol &mol, const double *dmat) {
    unsigned int na = mol.getNumAtoms();
    double d1 = dmat[beg1*na + beg2];
    double d2 = dmat[beg1*na + end2];
    double d3 = dmat[end1*na + beg2];
    double d4 = dmat[end1*na + end2];
    LIST_PAIR_DII dAtomsList;
    dAtomsList.push_back(PAIR_D_I_I(d1, PAIR_I_I(beg1, beg2)));
    dAtomsList.push_back(PAIR_D_I_I(d2, PAIR_I_I(beg1, end2)));
    dAtomsList.push_back(PAIR_D_I_I(d3, PAIR_I_I(end1, beg2)));
    dAtomsList.push_back(PAIR_D_I_I(d4, PAIR_I_I(end1, end2)));
    dAtomsList.sort(_pairDIICompAscending);
    return dAtomsList.front().second;
  }
  
  void EmbeddedFrag::computeDistMat(DOUBLE_SMART_PTR &dmat) {
    unsigned ai, aj;
    INT_EATOM_MAP_I efi, efj;
    RDGeom::Point2D pti, ptj;

    INT_EATOM_MAP_I tempi = d_eatoms.begin();
    tempi++;
    double *dmatPtr = dmat.get();
    for (efi = tempi; efi != d_eatoms.end(); efi++) {
      pti = efi->second.loc;
      ai = efi->first;
      for (efj = d_eatoms.begin(); efj != efi; efj++) {
        ptj = efj->second.loc;
        aj = efj->first;
        ptj -= pti;
        if (ai < aj) {
          dmatPtr[(aj*(aj-1)/2) + ai] = ptj.length();
        } else {
          dmatPtr[(ai*(ai-1)/2) + aj] = ptj.length();
        }
      }
    }
  }
  
  double EmbeddedFrag::mimicDistMatAndDensityCostFunc(const DOUBLE_SMART_PTR *dmat, double mimicDmatWt) {

    const double *ddata;
    if (dmat) {
      ddata = dmat->get();
    } else {
      ddata = 0;
    }
    unsigned int na = dp_mol->getNumAtoms(); 
    unsigned int dsize = na*(na-1)/2;
    double *ddata2D = new double[dsize];
    DOUBLE_SMART_PTR dmat2D(ddata2D);
    this->computeDistMat(dmat2D);
    double res1 = 0.0;
    double res2 = 0.0;
    double d, d2, dd;

    for (unsigned int i = 0; i < dsize; ++i) {
      d = ddata2D[i];
      d2 = d*d;
      if (d2 > 1.e-3){
        res1 += 1.0/d2;
      } else {
        res1 += 1000.0;
      }
      if ((ddata) && (ddata[i] >= 0.0)) {
        dd = d - ddata[i];
        res2 += dd*dd;
      }
    }
    
    double wt = mimicDmatWt;
    if (wt > 1.0) {
      wt = 1.0;
    } else if (wt < 0.0) {
      wt = 0.0;
    }
    
    return ((1.0-wt)*res1) + (wt*res2);
  }
  
  
  // Permute the bonds at a degree 4 node
  //
  //      A                    B
  //      |                    |
  //   B--C--D     to       A--C--D
  //      |                    |
  //      E                    E
  //
  // Note that everything attached to B and A are also effected. This is what
  // happnds here
  // 1. Find the line "l" bisecting the angle BCA
  // 2. Find the atoms in the fragment generated by breaking the bond between C and A
  //    that includes A. Lets call is Fa
  // 3. Similarly find the fragment Fb that includes B by breaking the bond CB
  // 4. Reflect Fb and Fa through "l"
  void EmbeddedFrag::permuteBonds(unsigned int aid, unsigned int aid1, unsigned int aid2) {
    PRECONDITION(dp_mol, "");
    //std::cerr<<"permute "<<aid<<" "<<aid1<<" "<<aid2<<std::endl;
    RDGeom::Point2D rl1 = d_eatoms[aid].loc;
    RDGeom::Point2D rl2 = d_eatoms[aid1].loc + d_eatoms[aid2].loc;
    rl2 *= 0.5;
    
    RDKit::INT_VECT fragA, fragB;

    // now find the fragment that contains aid1 but not aid
    _recurseAtomOneSide(aid1, aid, dp_mol, fragA); 

    // now find the fragment that contains aid2 but not aid
    _recurseAtomOneSide(aid2, aid, dp_mol, fragB);

    // now just loop through these atoms and reflect them
    RDKit::INT_VECT_CI fi;
    for (fi = fragA.begin(); fi != fragA.end(); fi++) {
      d_eatoms[*fi].Reflect(rl1, rl2);
    }

    for (fi = fragB.begin(); fi != fragB.end(); fi++) {
      d_eatoms[*fi].Reflect(rl1, rl2);
    }
  }

  void EmbeddedFrag::randomSampleFlipsAndPermutations(unsigned int nBondsPerSample,
                                                      unsigned int nSamples, int seed, 
                                                      const DOUBLE_SMART_PTR *dmat, 
                                                      double mimicDmatWt, bool permuteDeg4Nodes) {
    PRECONDITION(dp_mol, "");

    RDKit::rng_type &generator = RDKit::getRandomGenerator();
    if (seed > 0) {
      generator.seed(seed);
    }

    RDKit::INT_VECT rotBonds = getAllRotatableBonds(*dp_mol);

    unsigned int nb = rotBonds.size(); // number of rotatable bonds that can be flipped

    // if we also want to permute deg 4 nodes, find out how many of these are 
    // around and can be permuted
    unsigned int nt, nd4;
    nd4 = 0;
    RDKit::INT_VECT deg4nodes;
    RDKit::VECT_INT_VECT deg4NbrBids, deg4NbrAids;
    
    if (permuteDeg4Nodes) {
      for (RDKit::ROMol::ConstAtomIterator ai = dp_mol->beginAtoms(); ai != dp_mol->endAtoms(); ai++) {
        unsigned int caid = (*ai)->getIdx();
        if ( ((*ai)->getDegree() == 4) && (!(dp_mol->getRingInfo()->numAtomRings(caid))) ) {
          RDKit::INT_VECT aids, bids;
          getNbrAtomAndBondIds(caid, dp_mol, aids, bids);
          // make sure all the atoms in aids are in this embeddedfrag
          bool allin = true;
          for (RDKit::INT_VECT_CI ivci = aids.begin(); ivci != aids.end(); ivci++) {
            if (d_eatoms.find(*ivci) == d_eatoms.end()) {
              allin = false;
            }
          }
          if (allin) {
            deg4nodes.push_back(caid);
            deg4NbrBids.push_back(bids);
            deg4NbrAids.push_back(aids);
          }
        }
      }
      nd4 = deg4nodes.size();
    }

    nt = nb + nd4;

    unsigned int nPerSample = std::min(nt, nBondsPerSample);

    RDKit::uniform_int dist(0, nt-1);
    RDKit::int_source_type intRandomSrc(generator, dist);

    unsigned int si, fi, bi, ai;
    RDGeom::INT_POINT2D_MAP bestCrdMap;
    double bestDens = this->mimicDistMatAndDensityCostFunc(dmat, mimicDmatWt); 
    INT_EATOM_MAP_I efi;
    for (efi = d_eatoms.begin(); efi != d_eatoms.end(); efi++) {
      bestCrdMap[efi->first] = efi->second.loc;
    }
    for (si = 0; si < nSamples; ++si) {
      // randomly pick nPerSample bonds and flip them 
      for (fi = 0; fi < nPerSample; ++fi) {
        unsigned int ri = intRandomSrc();
        // if ri is less than the number of rotatable bonds (nb), we will flip a rot bond
        if (ri < nb) {
          bi = rotBonds[ri];
          this->flipAboutBond(bi);

        } else { // ri is >= nb we permute the bonds at a deg 4 node
          unsigned int d4i = ri - nb; // so we will permute at the 'di'th degree 4 node
          ai = deg4nodes[d4i];
          // collect the locations for the neighbors
          VECT_C_POINT nbrLocs;
          for (RDKit::INT_VECT_CI aci = deg4NbrAids[d4i].begin();
               aci != deg4NbrAids[d4i].end(); aci++) {
            nbrLocs.push_back(&(d_eatoms[*aci].loc));
          }
          INT_PAIR_VECT bndPairs = findBondsPairsToPermuteDeg4(d_eatoms[ai].loc, deg4NbrBids[d4i], nbrLocs);
          
          double rval = RDKit::getRandomVal();
          unsigned int fbi = 0;
          if (rval > 0.5) {
            fbi = 1;
          }
          unsigned int aid1, aid2;
          aid1 = dp_mol->getBondWithIdx(bndPairs[fbi].first)->getOtherAtomIdx(ai);
          aid2 = dp_mol->getBondWithIdx(bndPairs[fbi].second)->getOtherAtomIdx(ai);
          this->permuteBonds(ai, aid1, aid2);
        }
      }
      
      // compute the density of the stucture and check if it improved
      double density = this->mimicDistMatAndDensityCostFunc(dmat, mimicDmatWt);
      //if (density < bestDens) {
      if (bestDens-density>1e-4){
        bestDens = density;
        for (efi = d_eatoms.begin(); efi != d_eatoms.end(); efi++) {
          bestCrdMap[efi->first] = efi->second.loc;
        }
      }
    }
    // now copy the best coordinates to the fragment
    for (efi = d_eatoms.begin(); efi != d_eatoms.end(); efi++) {
      efi->second.loc = bestCrdMap[efi->first];
    }
  }

  std::vector<PAIR_I_I> EmbeddedFrag::findCollisions(const double *dmat,
                                                     bool includeBonds) {
    // find a pair of atoms that are too close to each other
    INT_EATOM_MAP_I efi, efj, tempi;
    RDGeom::Point2D pti, ptj;
    double d2;
    std::vector<PAIR_I_I> res;
    for (efi = d_eatoms.begin(); efi != d_eatoms.end(); ++efi) {
      efi->second.d_density = 0.0;
    }

    tempi = d_eatoms.begin();
    ++tempi;
    double colThres2 = COLLISION_THRES*COLLISION_THRES;
    // if we a re dealing with non carbon atoms we will increase the collision threshold.
    // This is because only hetero atoms are typically drawn in a depiction.
    double atomTypeFactor1, atomTypeFactor2;
    for (efi = tempi; efi != d_eatoms.end(); efi++) {
      pti = efi->second.loc;
      atomTypeFactor1 = 1.0;
      if (dp_mol->getAtomWithIdx(efi->first)->getAtomicNum() != 6) {
        atomTypeFactor1 = HETEROATOM_COLL_SCALE;
      }
      for (efj = d_eatoms.begin(); efj != efi; efj++) {
        if (efj == efi) {
          continue;
        }
        atomTypeFactor2 = 1.0;
        if (dp_mol->getAtomWithIdx(efj->first)->getAtomicNum() != 6) {
          atomTypeFactor2 = HETEROATOM_COLL_SCALE;
        }
        ptj = efj->second.loc;
        ptj -= pti;
        d2 = ptj.lengthSq();
        if (d2 > 1.0e-3) {
          efi->second.d_density += (1/d2);
          efj->second.d_density += (1/d2);
        } else {
          efi->second.d_density += 1000.0;
          efj->second.d_density += 1000.0;
        }
        d2 /= (atomTypeFactor1*atomTypeFactor2);
        //std::cerr<<" "<<efi->first<<"-"<<efj->first<<": "<<d2<<" "<<colThres2<<std::endl;
        if (d2 < colThres2) {
          PAIR_I_I cAids(efi->first, efj->first);
          res.push_back(cAids);
        }
      }
    }
    if (includeBonds) {
      // now find bond collisions
      RDKit::ROMol::ConstBondIterator bi1, bi2;
      unsigned int bid1, bid2;
      unsigned int beg1, end1, beg2, end2;
      RDGeom::Point2D avg1, avg2, v1, v2, v3;
      double BOND_THRES2 = BOND_THRES*BOND_THRES;
      double valProd;
      for (bi1 = dp_mol->beginBonds(); bi1 != dp_mol->endBonds(); bi1++) {
        bid1 = (*bi1)->getIdx();
        beg1 = (*bi1)->getBeginAtomIdx();
        end1 = (*bi1)->getEndAtomIdx();
        if ((d_eatoms.find(beg1) != d_eatoms.end()) &&
            (d_eatoms.find(end1) != d_eatoms.end())) {
          v1 = d_eatoms[end1].loc - d_eatoms[beg1].loc;
          avg1 = d_eatoms[end1].loc + d_eatoms[beg1].loc;
          avg1 *= 0.5;
          for (bi2 = dp_mol->beginBonds(); bi2 != dp_mol->endBonds(); bi2++) {
            bid2 = (*bi2)->getIdx();
            if (bid2 <= bid1) {
              continue;
            }
            
            beg2 = (*bi2)->getBeginAtomIdx();
            end2 = (*bi2)->getEndAtomIdx();
            if ((d_eatoms.find(beg2) != d_eatoms.end()) && 
                (d_eatoms.find(end2) != d_eatoms.end())) {
              avg2 = d_eatoms[end2].loc + d_eatoms[beg2].loc;
              avg2 *= 0.5;
              avg2 -= avg1;
              if(avg2.lengthSq()<0.5 &&
                 avg2.lengthSq() < BOND_THRES2) {
                v2 = d_eatoms[beg2].loc - d_eatoms[beg1].loc;
                v3 = d_eatoms[end2].loc - d_eatoms[beg1].loc;
                valProd = _crossVal(v1, v2)*_crossVal(v1,v3);
                if (valProd < -1e-6) {
                  // we have a collision, find the closest two atoms
                  PAIR_I_I cAids = _findClosestPair(beg1, end1, beg2, end2,
                                                    *dp_mol, dmat);
                  res.push_back(cAids);
                }
              }
            }
          }
        }
      }
    }
    return res;
  }


  double EmbeddedFrag::totalDensity() {
    INT_EATOM_MAP_I efi;
    double res = 0.0;
    for (efi = d_eatoms.begin(); efi != d_eatoms.end(); efi++) {
      res += efi->second.d_density;
    }
    return res;
  }

  void _recurseDegTwoRingAtoms(unsigned int aid, const RDKit::ROMol *mol, RDKit::INT_VECT &rPath,
                               RDKit::INT_INT_VECT_MAP &nbrMap) {
    PRECONDITION(mol,"");
    // find all atoms along a path that have two ring atoms on them
    // aid is where will start looking and then we will recurse
    RDKit::ROMol::OBOND_ITER_PAIR atomBonds;
    atomBonds = mol->getAtomBonds(mol->getAtomWithIdx(aid));
    int bondId;
    RDKit::INT_VECT nbrs;
    while( atomBonds.first != atomBonds.second ){
      const RDKit::BOND_SPTR bnd = (*mol)[*atomBonds.first];
      bondId = bnd->getIdx();
      if (mol->getRingInfo()->numBondRings(bondId)) {
        nbrs.push_back(bnd->getOtherAtomIdx(aid));
      }
      atomBonds.first++;
    }
    if (nbrs.size() != 2) {
      return;
    } else {
      rPath.push_back(aid);
      nbrMap[aid] = nbrs;
      RDKit::INT_VECT_CI nbi;
      for (nbi = nbrs.begin(); nbi != nbrs.end(); nbi++) {
        if (std::find(rPath.begin(), rPath.end(), (*nbi)) == rPath.end()) {
          _recurseDegTwoRingAtoms((*nbi), mol, rPath, nbrMap);
        }
      }
    }
  }

  int _anyNonRingBonds(unsigned int aid, RDKit::INT_LIST path, const RDKit::ROMol *mol) {
    PRECONDITION(mol,"");
    // check if there are any non-ring bonds on the path starting at aid
    const RDKit::Bond *bond;
    int prev = aid;
    int nOpen=0;
    RDKit::INT_LIST_CI pi;
    for (pi = path.begin(); pi != path.end(); pi++) {
      bond = mol->getBondBetweenAtoms(prev, (*pi));
      if (!mol->getRingInfo()->numBondRings(bond->getIdx())) {
        nOpen++;
      }
      prev = (*pi);
    }
    return nOpen;
  }

  void EmbeddedFrag::flipAboutBond(unsigned int bondId,bool flipEnd) {
    PRECONDITION(dp_mol, "");
    PRECONDITION(bondId < dp_mol->getNumBonds(), "");
    //std::cerr<<"  flip about: "<<bondId<<" "<<flipEnd<<std::endl;
    // reflect all the atoms on one side of a bond using the bond as the mirror
    const RDKit::Bond *bond = dp_mol->getBondWithIdx(bondId);
    
    // we should not be flip things around a ring bond
    CHECK_INVARIANT(!(dp_mol->getRingInfo()->numBondRings(bondId)), "");
      
    int begAid = bond->getBeginAtomIdx();
    int endAid = bond->getEndAtomIdx();

    if(!flipEnd){
      int tmp=begAid;
      begAid=endAid;
      endAid=tmp;
    }
    
    RDGeom::Point2D begLoc = d_eatoms[begAid].loc;
    RDGeom::Point2D endLoc = d_eatoms[endAid].loc;

    // arbitrary choice here - find all atoms on one side of the bond
    // endAtom side - we will do this recursively
    RDKit::INT_VECT endSideAids;
    endSideAids.clear();
    _recurseAtomOneSide(endAid, begAid, dp_mol, endSideAids);

    // noew we have the molecule split into two groups of atoms 
    // atom on the side of endAid and the rest. 
    // we will flip the side that is smaller
    int nats = d_eatoms.size();
    int nEndSide = endSideAids.size();
    bool endSideFlip = true;
    if ((nats - nEndSide) < nEndSide) {
      endSideFlip = false;
    }

    for (INT_EATOM_MAP_I efi = d_eatoms.begin(); efi != d_eatoms.end(); efi++) {
      RDKit::INT_VECT_CI fii = std::find(endSideAids.begin(), endSideAids.end(),
                                         static_cast<int>(efi->first));
      if (endSideFlip ^ (fii == endSideAids.end()) ) {
        efi->second.Reflect(begLoc, endLoc);
      }
    }
  }  

  unsigned int _findDeg1Neighbor(const RDKit::ROMol *mol, unsigned int aid) {
    PRECONDITION(mol, "");
    unsigned int deg = mol->getAtomWithIdx(aid)->getDegree();
    CHECK_INVARIANT(deg == 1, "");
    unsigned int res=0;
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(mol->getAtomWithIdx(aid));
    res = (*nbrIdx);
#if 0      
    while (nbrIdx != endNbrs) {
      if(mol->getAtomWithIdx(*nbrIdx)->getDegree()==1){
        res = (*nbrIdx);
        break;
      }
      ++nbrIdx;
    }
#endif
    return res;
  }

  unsigned int _findClosestNeighbor(const RDKit::ROMol *mol, const double *dmat,
                                    unsigned int aid1, unsigned int aid2) {
    PRECONDITION(mol, "");
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(mol->getAtomWithIdx(aid2));
    unsigned int res=0;
    double d, mdist = 1.e8;
    unsigned int naid = aid1*(mol->getNumAtoms());
    while (nbrIdx != endNbrs) {
      d = dmat[naid + (*nbrIdx)];
      if (d < mdist){
        mdist = d;
        res = (*nbrIdx);
      }
      nbrIdx++;
    }
    return res;
  }
    
  void EmbeddedFrag::openAngles(const double *dmat, unsigned int aid1, unsigned int aid2) {
    // Assuming that either aid1, and/or aid2 are degree 1 atoms, we will open up the angles 
    //
    //     1 2
    //    /   \                                                   this space
    //   /     \                                            intentionally left blank
    //  a-------b
    //
    // If 1 and 2 are too close to each other we open up angle(1ab) if 1 is a degree 1 node and
    // angle(2ba) if 2 is a degree 1 node. Say 1 is a degree 1 node but 2 is not. 
    // Then from the neighbors of 2 we need to choose which one should be b. Also keep in mind 
    // that a need not be a neighbor of b. In this case we will pick b to be the closest neighbor of a
    
    PRECONDITION(dp_mol, "");
    PRECONDITION(dmat, "");
    unsigned int deg1 = dp_mol->getAtomWithIdx(aid1)->getDegree();
    unsigned int deg2 = dp_mol->getAtomWithIdx(aid2)->getDegree();
    if ((deg1 > 1) && (deg2 > 1) ) {
      return;
    }
    unsigned int aidA;
    unsigned int aidB;
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    int type = 0;
    if ((deg1 == 1) && (deg2 == 1)) {
      aidA = _findDeg1Neighbor(dp_mol, aid1);
      aidB = _findDeg1Neighbor(dp_mol, aid2);
      type = 1;
    } else if ((deg1 == 1) && (deg2 > 1)) {
      aidA = _findDeg1Neighbor(dp_mol, aid1);
      aidB = _findClosestNeighbor(dp_mol, dmat, aidA, aid2);
      type = 2;
    } else {
      aidB = _findDeg1Neighbor(dp_mol, aid2);
      aidA = _findClosestNeighbor(dp_mol, dmat, aidB, aid1);
      type = 3;
    }

    //std::cerr<<" openAngles: "<<aid1<<"-"<<aidA<<"-"<<aidB<<"-"<<aid2<<"  type: "<<type<<std::endl;
    //std::cerr<<"             len: "<<(d_eatoms[aid1].loc - d_eatoms[aid2].loc).length()<<std::endl;
    
    RDGeom::Point2D v2 = d_eatoms[aid1].loc - d_eatoms[aidA].loc;
    RDGeom::Point2D v1 = d_eatoms[aidB].loc - d_eatoms[aidA].loc;
    double cross = (v1.x)*(v2.y) - (v1.y)*(v2.x);
    double angle;
    RDGeom::Transform2D trans1, trans2;
    switch (type) {
    case 1:
      angle = ANGLE_OPEN;
      if (cross < 0) {
        angle *= -1.0;
      }
      trans1.SetTransform(d_eatoms[aidA].loc, angle);
      trans2.SetTransform(d_eatoms[aidB].loc, -1.0*angle);
      trans1.TransformPoint(d_eatoms[aid1].loc);
      trans2.TransformPoint(d_eatoms[aid2].loc);
      break;
    case 2:
      angle = 2.0*ANGLE_OPEN;
      if (cross < 0) {
        angle *= -1.0;
      }
      trans1.SetTransform(d_eatoms[aidA].loc, angle);
      trans1.TransformPoint(d_eatoms[aid1].loc);
      break;
    case 3:
      angle = -2.0*ANGLE_OPEN;
      if (cross < 0) {
        angle *= -1.0;
      }
      trans2.SetTransform(d_eatoms[aidB].loc, angle);
      trans2.TransformPoint(d_eatoms[aid2].loc);
      break;
    default:
      break;
    }
    //std::cerr<<"             post len: "<<(d_eatoms[aid1].loc - d_eatoms[aid2].loc).length()<<std::endl;
  }
  
  void EmbeddedFrag::removeCollisionsBondFlip() {
    // try to remove collisions in a structure by flipping rotatable bonds
    // along the shortest path between the colliding atoms.
    unsigned int iter = 0;
    // we will limit the number of times we are going to do this since
    // we may fall into spiral where removing a collision may
    // create a new one
    std::vector<PAIR_I_I> colls;
    double *dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
    colls = this->findCollisions(dmat);
    unsigned int ncols;
    //std::cerr<<"removeCollisionsBondFlip(): "<<colls.size()<<std::endl;
    std::map<int,unsigned int> doneBonds;
    while (iter < MAX_COLL_ITERS && colls.size()) {
      ncols = colls.size();
      //std::cerr<<"iter: "<<iter<<" "<<ncols<<std::endl;
      if (ncols > 0) {
        // we have a collision
        PAIR_I_I cAids = colls[0];
        RDKit::INT_VECT rotBonds = getRotatableBonds(*dp_mol, cAids.first,
                                                     cAids.second);
        RDKit::INT_VECT_CI ri;
        double prevDensity = this->totalDensity();
        //std::cerr<<"   density: "<<prevDensity<<std::endl;
        for (ri = rotBonds.begin(); ri != rotBonds.end(); ri++) {
          if ((doneBonds.find(*ri) == doneBonds.end()) ||
              (doneBonds[*ri] < NUM_BONDS_FLIPS)) { 
            if (doneBonds.find(*ri) == doneBonds.end()) {
              doneBonds[*ri] = 1;
            } else {
              doneBonds[*ri] += 1;
            }

            flipAboutBond((*ri));
            colls = this->findCollisions(dmat);
            double newDensity = this->totalDensity();
            //std::cerr<<"  newcolls: "<<colls.size()<<" "<<newDensity<<std::endl;
            if (colls.size() < ncols) {
              doneBonds[*ri] = NUM_BONDS_FLIPS; // lock this rotatable bond
              break;
            }
            else if (colls.size()==ncols && newDensity<prevDensity) {
              break;
            } else  {
              // we made the wrong move earlier - reject the flip move it back
              flipAboutBond((*ri));
              colls = this->findCollisions(dmat);
              // and try the other end:
              flipAboutBond((*ri),false);
              colls = this->findCollisions(dmat);
              newDensity = this->totalDensity();
              //std::cerr<<"  newcolls2: "<<colls.size()<<" "<<newDensity<<std::endl;
              if (colls.size() < ncols) {
                doneBonds[*ri] = NUM_BONDS_FLIPS; // lock this rotatable bond
                break;
              } else if (colls.size()==ncols && newDensity<prevDensity) {
                break;
              } else {
                flipAboutBond((*ri),false);
                colls = this->findCollisions(dmat);
              }
            }
          }
        }
      }
      iter++;
    }
  }
  
  void EmbeddedFrag::removeCollisionsOpenAngles() {
    double *dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
    std::vector<PAIR_I_I> colls = this->findCollisions(dmat, 0);
    // try opening up angles
    std::vector<PAIR_I_I>::const_iterator cpi;
    for (cpi = colls.begin(); cpi != colls.end(); cpi++) {
      // find out which of the two offending atoms we want to move
      // we will use the one with the smallest degree
      int aid1 = cpi->first;
      int aid2 = cpi->second;
      this->openAngles(dmat, aid1, aid2);
    }
  }

  void EmbeddedFrag::removeCollisionsShortenBonds() {
    double *dmat = RDKit::MolOps::getDistanceMat(*dp_mol);
    // if there are still some collision points left - flipping rotatable bonds
    // and opening angles is not doing it - we will try two last things
    //  - if all the bonds between the colliding atoms are rings bonds,
    //    we most likely have a collision within a bridged system (Issue 199).
    //    In this case we will try to find a path of colliding atoms (in one
    //    of the rings) and shorten all the bond in the path
    //  - on the other hand if we have non-ring bonds as well in the path
    //    between the colliding atoms we will simply shorten each one of
    //    them by a little bit.
    std::vector<PAIR_I_I> colls = this->findCollisions(dmat, 0);
    unsigned int ncols = colls.size();
    unsigned int iter = 0;
    while ( ncols && iter < MAX_COLL_ITERS ) {
      PAIR_I_I cAids = colls.front();
      // find out which of the two offending atoms we want to move
      // we will use the one with the smallest degree
      int aid1 = cAids.first;
      int aid2 = cAids.second;
      int deg1 = dp_mol->getAtomWithIdx(aid1)->getDegree();
      int deg2 = dp_mol->getAtomWithIdx(aid2)->getDegree();
      if (deg2 > deg1) {
        // reverse the order
        int temp = aid1;
        aid1 = aid2;
        aid2 = temp;
        temp=deg1;
        deg1=deg2;
        deg2=temp;
      }
      // now find the path between the two ends
      RDKit::INT_LIST path = RDKit::MolOps::getShortestPath(*dp_mol, aid1, aid2);
      //std::cerr<<" collide! "<<aid1<<" "<<aid2<<" "<<path.size()<<std::endl;
      if(!path.size()){
        // there's no path between the ends, so there's nothing
        // we can really do about this collision.
        colls.erase(colls.begin());
      } else {
       // aid1 is on the front of the path, pop it off:
       CHECK_INVARIANT(path.front()==aid1,"bad path head");
       path.pop_front();

       int nOpen = _anyNonRingBonds(aid1, path, dp_mol);
       //std::cerr<<"     nOpen: "<<nOpen<<std::endl;
       if (nOpen > 0) {
         if(deg1==1){
           RDGeom::Point2D loc = d_eatoms[aid1].loc;
           int aidA=_findDeg1Neighbor(dp_mol,aid1);
           loc -= d_eatoms[aidA].loc;
           loc *= .9;
           //std::cerr<<"  >>> "<<aid1<<" "<<loc.length()<<std::endl;
           if(loc.length()>.75){
             loc += d_eatoms[aidA].loc;
             d_eatoms[aid1].loc=loc;
           }
         }
         if(deg2==1){
           RDGeom::Point2D loc = d_eatoms[aid2].loc;
           int aidA=_findDeg1Neighbor(dp_mol,aid2);
           loc -= d_eatoms[aidA].loc;
           loc *= .9;
           //std::cerr<<"  >>> "<<aid2<<" "<<loc.length()<<std::endl;
           if(loc.length()>.75){
             loc += d_eatoms[aidA].loc;
             d_eatoms[aid2].loc=loc;
           }
         }
       } else {
         // we probably have a bridged system
         // lets hope that aids has only two ring bond on it
         RDKit::INT_VECT rPath;
         RDKit::INT_INT_VECT_MAP nbrMap;
         _recurseDegTwoRingAtoms(aid1, dp_mol, rPath, nbrMap);
         if (rPath.size() == 0) {
           _recurseDegTwoRingAtoms(aid2, dp_mol, rPath, nbrMap);
         }
         // now we will take each of the atoms in rPath and
         // "move them in" a little bit this is what "move them
         //  in" means (what we need is hand drawn picture in the comments)
         // - let r1 and r2 be the ring neighbor of the current atom r0
         // - we will find the vector that bisects angle(r1, r0, r2) 
         // - we will move r0 along this vector
         RDKit::INT_VECT_CI rpi;
         RDGeom::INT_POINT2D_MAP moveMap;
         for (rpi = rPath.begin(); rpi != rPath.end(); rpi++) {
           RDGeom::Point2D move;
           move = d_eatoms[nbrMap[*rpi][0]].loc;
           move += d_eatoms[nbrMap[*rpi][1]].loc;
           move *= 0.5;
           move -= d_eatoms[*rpi].loc;
           move.normalize();
           move *= COLLISION_THRES;
           moveMap[*rpi] = move;
         } 
         for (rpi = rPath.begin(); rpi != rPath.end(); rpi++) {
           d_eatoms[*rpi].loc += moveMap[*rpi];
         }
       }
       colls = this->findCollisions(dmat,0);
      }
      ncols = colls.size();
      ++iter;
    }
  }

}
    
    
