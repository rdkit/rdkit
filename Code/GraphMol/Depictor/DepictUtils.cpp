// $Id$
//
//  Copyright (C) 2003-2010 greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/types.h>
#include <math.h>
#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include "DepictUtils.h"
#include <iostream>
#include <RDGeneral/Invariant.h>
#include <algorithm>

namespace RDDepict {
  double BOND_LEN=1.5;
  double COLLISION_THRES=0.70;
  double BOND_THRES = 0.50;
  double ANGLE_OPEN=0.1222; // that is about 7 deg
  unsigned int MAX_COLL_ITERS=15;
  double HETEROATOM_COLL_SCALE=1.3;
  unsigned int NUM_BONDS_FLIPS=3;
  
  RDGeom::INT_POINT2D_MAP embedRing(const RDKit::INT_VECT &ring) {
    // The process here is very straight forward
    // we take the center of the ring to lies at the origin put the first 
    // point at the orgin anf then sweep
    // anticlock wise so by an angle A = 360/n for the next point
    // the length of the arm (l) we want to sweep is easy to compute given the 
    // bond lenght (b) we want to use for each bond in the ring (for now
    // we will assume that this bond legnth is the same for all bonds in the ring
    //  l = b/sqrt(2*(1 - cos(A)) 
    // the above formular derives from the traingle formula, where side 'c' is given
    // interms of sides 'a' and 'b' as 
    // c = a^2 + b^2 - 2.a.b.cos(A) 
    // where A is the angle between a and b
    
    // compute the sweep angle
    unsigned int na = ring.size();
    double ang = 2*M_PI/na;
    
    // compute the arm length
    double al = BOND_LEN/(sqrt(2*(1 - cos(ang))));

    RDGeom::INT_POINT2D_MAP res;

    unsigned int i, aid;
    double x,y;

    for (i = 0; i < na; i++) {
      x = al*cos(i*ang);
      y = al*sin(i*ang);
      RDGeom::Point2D loc(x,y);
      aid = ring[i];
      res[aid] = loc;
    }
    return res;
  }

  void transformPoints(RDGeom::INT_POINT2D_MAP &nringCor, const RDGeom::Transform2D &trans) {
    RDGeom::INT_POINT2D_MAP_I nrci;
    for (nrci = nringCor.begin(); nrci != nringCor.end(); nrci++) {
      RDGeom::Point2D loc = nrci->second;
      trans.TransformPoint(loc);
      nrci->second = loc;
    }
    
  }

  RDGeom::Point2D computeBisectPoint(const RDGeom::Point2D &rcr, 
                                     double ang, const RDGeom::Point2D &nb1,
                                     const RDGeom::Point2D &nb2) {
    
    RDGeom::Point2D cloc = nb1;
    cloc += nb2;
    cloc *= 0.5;
    if (ang > M_PI) {
      // invert the cloc
      cloc -= rcr;
      cloc *= -1.0;
      cloc += rcr;
    }
    return cloc;
  }
  
  RDGeom::Point2D reflectPoint(const RDGeom::Point2D &point, const RDGeom::Point2D &loc1,
                   const RDGeom::Point2D &loc2) {
    RDGeom::Point2D org(0.0, 0.0);
    RDGeom::Point2D xaxis(1.0, 0.0);
    RDGeom::Point2D cent = (loc1 + loc2);
    cent *= 0.5;
    
    RDGeom::Transform2D trans;
    trans.SetTransform(org, xaxis, cent, loc1);
    
    /// reverse transform
    RDGeom::Transform2D itrans;
    itrans.SetTransform(cent, loc1, org, xaxis);

    RDGeom::INT_POINT2D_MAP_I nci;
    RDGeom::Point2D res;
    res = point;
    trans.TransformPoint(res);
    res.y = -res.y;
    itrans.TransformPoint(res);
    return res;
  }
    
  void reflectPoints(RDGeom::INT_POINT2D_MAP &coordMap, const RDGeom::Point2D &loc1,
                   const RDGeom::Point2D &loc2) {
    RDGeom::INT_POINT2D_MAP_I nci;
    for (nci = coordMap.begin(); nci != coordMap.end(); nci++) {
      nci->second = reflectPoint(nci->second, loc1, loc2);
    }
  }

  RDKit::INT_VECT setNbrOrder(unsigned int aid, const RDKit::INT_VECT &nbrs,
                              const RDKit::ROMol &mol) {
    PRECONDITION(aid<mol.getNumAtoms(), "");
    PR_QUEUE subsAid;
    int ref=-1;
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(aid));
    // find the neighbor of aid that is not in nbrs i.e. atom A from the comments in the header file
    // and the store the pair <degree, aid> in the order of increasing degree
    while (nbrIdx != endNbrs) {
      // We used to use degree here instead we will start using the CIP rank here
      if (std::find(nbrs.begin(), nbrs.end(), static_cast<int>(*nbrIdx)) == nbrs.end()) {
        ref = (*nbrIdx);
      }
      nbrIdx++;
    }
    
    RDKit::INT_VECT thold = nbrs;
    if (ref >= 0) {
      thold.push_back(ref);
    }
    // we should be here unless we have more than 3 atoms to worry about
    CHECK_INVARIANT(thold.size() > 3, "");
    thold = rankAtomsByRank(mol, thold);
        
    // swap the position of the 3rd to last and second to last items in sorted list
    unsigned int ln = thold.size();
    int tint = thold[ln-3];
    thold[ln-3] = thold[ln-2];
    thold[ln-2] = tint;
    

    // go clock wise along the list from this position for the arranged neighbor list
    RDKit::INT_VECT res;
    res.reserve(thold.size());
    RDKit::INT_VECT_I pos=std::find(thold.begin(),thold.end(),ref);
    if(pos!=thold.end()){
      res.insert(res.end(),pos+1,thold.end());
    }
    if(pos!=thold.begin()){
      res.insert(res.end(),thold.begin(),pos);
    }

    POSTCONDITION(res.size() == nbrs.size(), "");
    return res;
  } 

  int pickFirstRingToEmbed(const RDKit::ROMol &mol, const RDKit::VECT_INT_VECT &fusedRings) {
    // ok this is what we will do here 
    // we will pick the ring with the smallest number of substiuents
    int res=-1;
    unsigned int maxSize = 0;
    int subs, minsubs = static_cast<int>(1e8);
    int cnt = 0;
    for (RDKit::VECT_INT_VECT_CI ri = fusedRings.begin();
         ri != fusedRings.end(); ri++) {
      subs = 0;
      for (RDKit::INT_VECT_CI rii = ri->begin(); rii != ri->end(); rii++) {
        int deg = mol.getAtomWithIdx(*rii)->getDegree();
        if (deg > 2) {
          subs++;
        }
      }
      if (subs < minsubs) {
        res = cnt;
        minsubs = subs;
        maxSize = ri->size();
      } else if (subs == minsubs) {
        if (ri->size() > maxSize) {
          res = cnt;
          maxSize = ri->size();
        }
      }
      cnt++;
    }
    return res;

  }

  RDKit::INT_VECT findNextRingToEmbed(const RDKit::INT_VECT &doneRings, 
                                        const RDKit::VECT_INT_VECT &fusedRings, 
                                       int &nextId) {
    // REVIEW: We are changing this after Issue166 
    // Originally the ring that have maximum number of atoms in common with the atoms
    // that have already been embedded will be the ring that will get embedded. But 
    // if we can find a ring with two atoms in common with the embedded atoms, we will 
    // choose that first before systems with more than 2 atoms in common. Cases with two atoms
    // in common are in general flat systems to start with and can be embedded cleanly. 
    // when there are more than 2 atoms in common, these are most likely bridged syste, which are 
    // screwed up anyway, might as well screw them up later
    // if we do not have a system with two rings in common then we will return the ring with max,
    // common atoms
    PRECONDITION(doneRings.size() > 0, "");
    PRECONDITION(fusedRings.size() > 1, "");

    RDKit::INT_VECT commonAtoms, res, doneAtoms, notDone;
    RDKit::INT_VECT_CI dri;
    for (unsigned int i = 0; i < fusedRings.size(); i++) {
      if (std::find(doneRings.begin(), doneRings.end(), static_cast<int>(i)) == doneRings.end()) {
        notDone.push_back(i);
      }
    }
    
    RDKit::Union(fusedRings, doneAtoms, &notDone);
    
    int maxCommonAtoms = 0;
    
    int currRingId = 0;
    for (RDKit::VECT_INT_VECT_CI ri = fusedRings.begin();
         ri != fusedRings.end(); ri++) {
      if (std::find(doneRings.begin(), doneRings.end(), currRingId) != doneRings.end()) {
        currRingId++;
        continue;
      }
      commonAtoms.clear();
      int numCommonAtoms = 0;
      for (RDKit::INT_VECT_CI rii = ri->begin(); rii != ri->end(); rii++) {
        if (std::find(doneAtoms.begin(), doneAtoms.end(), (*rii)) != doneAtoms.end()) {
          commonAtoms.push_back(*rii);
          numCommonAtoms++;
        }
      }
      if (numCommonAtoms == 2) {
        // if we found a ring with two atoms in common get out 
        nextId = currRingId;
        return commonAtoms; // FIX: this causes the rendering to be non-canonical
      }
      if (numCommonAtoms > maxCommonAtoms) {
        maxCommonAtoms = numCommonAtoms;
        nextId = currRingId;
        res = commonAtoms;
      }
      currRingId++;
    }
    // here is an additional constrain we will put on the common atoms 
    // it is quite likely that the common atoms form a chain (it is possible we can 
    // construct some weird cases where this does not hold true - but for now we will
    // assume this is true. However the IDs in the res may not be in the order of going
    // from one end of the chain to the other - here is an example
    // C1CCC(CC12)CCC2 - two rings here with three atoms in common
    // let ring1:(0,1,2,3,4,5) be a ring that is already embedded, then let ring2:(4,3,6,7,8,5) be the ring
    // that we found to be the next ring we should embed. The commonAtoms are (4,3,5) - note that
    // they will be in this order since the rings are always traversed in order. Now we would like these 
    // common atoms to be returned in the order (5,4,3) - then we have a continuous chain, we can 
    // do this by simply looking at the original ring order (4,3,6,7,8,5) and observing that 5 need to come to
    // the front
    
    // find out how many atoms from the end we need to move to the front
    unsigned int cmnLst = 0;
    unsigned int nCmn  = res.size();
    for (unsigned int i = 0; i < nCmn; i++) {
      if (res[i] == fusedRings[nextId][i]) {
        cmnLst++;
      } else {
        break;
      }
    }
    // now do the moving if we have to
    if ((cmnLst > 0) && (cmnLst < res.size()) ) {
      RDKit::INT_VECT tempV = res;
      
      for (unsigned int i = cmnLst; i < nCmn; i++) {
        res[i - cmnLst] = tempV[i];
      }
      unsigned int nMov = nCmn - cmnLst;
      for (unsigned int i = 0; i < cmnLst; i++) {
        res[nMov + i] = tempV[i];
      }
    }
   
    POSTCONDITION(res.size()>0,"");
    return res;
  }

  RDKit::INT_VECT getAllRotatableBonds(const RDKit::ROMol &mol) {
    RDKit::INT_VECT res;
    RDKit::ROMol::ConstBondIterator bondIt;
    for(bondIt=mol.beginBonds();bondIt!=mol.endBonds();bondIt++){
      int bid = (*bondIt)->getIdx();
      if ( ((*bondIt)->getStereo() <= RDKit::Bond::STEREOANY) &&
           (!(mol.getRingInfo()->numBondRings(bid))))  {
          res.push_back(bid);
      }
    }
    return res;
  }

  
  RDKit::INT_VECT getRotatableBonds(const RDKit::ROMol &mol, unsigned int aid1, unsigned int aid2) {
    PRECONDITION(aid1<mol.getNumAtoms(), "");
    PRECONDITION(aid2<mol.getNumAtoms(), "");

    RDKit::INT_LIST path = RDKit::MolOps::getShortestPath(mol, aid1, aid2);
    RDKit::INT_VECT res;
    if (path.size() >= 4) {
      // remove the first atom (aid1) and last atom (aid2)
      CHECK_INVARIANT(static_cast<unsigned int>(path.front())==aid1,"bad first element");
      path.pop_front();
      CHECK_INVARIANT(static_cast<unsigned int>(path.back())==aid2,"bad last element");
      path.pop_back();
      
      RDKit::INT_LIST_CI pi = path.begin();
      int pid = (*pi);
      ++pi;
      while (pi != path.end()) {
        int aid = (*pi);
        const RDKit::Bond *bond = mol.getBondBetweenAtoms(pid, aid);
        int bid = bond->getIdx();
        if ( (bond->getStereo() <= RDKit::Bond::STEREOANY) &&
             (!(mol.getRingInfo()->numBondRings(bid))) ) {
          res.push_back(bid);
        }
        pid = aid;
        ++pi;
      }
    }
    return res;
  }

  void getNbrAtomAndBondIds(unsigned int aid, const RDKit::ROMol *mol,
                            RDKit::INT_VECT &aids, RDKit::INT_VECT &bids) {
    CHECK_INVARIANT(mol, "");
    unsigned int na = mol->getNumAtoms();
    RANGE_CHECK(0, aid, na-1);
    
    RDKit::ROMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = mol->getAtomNeighbors(mol->getAtomWithIdx(aid));

    unsigned int ai, bi;
    while (nbrIdx != endNbrs) {
      ai = (*nbrIdx);
      bi = mol->getBondBetweenAtoms(aid, ai)->getIdx();
      aids.push_back(ai);
      bids.push_back(bi);
      nbrIdx++;
    }
    
  }

  // find pairs of bonds that can be permuted at a non-ring degree 4 
  // node. This function will return only those pairs that cannot be 
  // permuted by flipping a rotatble bond
  // 
  //       D
  //       |
  //       b3
  //       |
  //  A-b1-B-b2-C
  //       |
  //       b4
  //       |
  //       E
  // For example in teh above situation on the pairs (b1, b3) and (b1, b4) will be returned
  // All other permutations can be achieved via a rotatable bond flip.
  INT_PAIR_VECT findBondsPairsToPermuteDeg4(const RDGeom::Point2D &center, const RDKit::INT_VECT &nbrBids, 
                                            const VECT_C_POINT &nbrLocs) {

    INT_PAIR_VECT res;
    
    // make sure there are four of them
    CHECK_INVARIANT(nbrBids.size() == 4, "");
    CHECK_INVARIANT(nbrLocs.size() == 4, "");
    
    VECT_C_POINT::const_iterator npi;
    
    std::vector<RDGeom::Point2D> nbrPts;
    RDKit::INT_VECT_CI aci;
    
    for (npi = nbrLocs.begin(); npi != nbrLocs.end(); npi++) {
      RDGeom::Point2D v = (*(*npi)) - center;
      nbrPts.push_back(v);
    }

    // now find the lay out of the bonds and return the bonds that are 90deg to the 
    // the bond to the first neighbor; i.e. we want to find b3 and b4 in the above picture
    double dp1 = nbrPts[0].dotProduct(nbrPts[1]);
    if (fabs(dp1) < 1.e-3) {
      // the first two vectors are perpendicular to each other. We now have b1 and b3 we need to
      // find b4
      INT_PAIR p1(nbrBids[0], nbrBids[1]);
      res.push_back(p1);

      double dp2 = nbrPts[0].dotProduct(nbrPts[2]);
      if (fabs(dp2) < 1.e-3) {
        // now we found b4 as well return the results
        INT_PAIR p2(nbrBids[0], nbrBids[2]);
        res.push_back(p2);
      } else {
        // bids[0] and bids[2] are opposite to each other and we know bids[1] is
        // perpendicular to bids[0]. So bids[3] is also perpendicular to bids[0]
        INT_PAIR p2(nbrBids[0], nbrBids[3]);
        res.push_back(p2);
      }
      return res;
    } else {
      // bids[0] and bids[1] are oppostie to each other, so bids[2] and bids[3] must
      // be perpendicular to bids[0]
      INT_PAIR p1(nbrBids[0], nbrBids[2]);
      res.push_back(p1);
      INT_PAIR p2(nbrBids[0], nbrBids[3]);
      res.push_back(p2);
      return res;
    }
  }

  // compare the first elements of two pairs of integers/
  
  int _pairCompDescending(const INT_PAIR &arg1,const INT_PAIR &arg2){
    return (arg1.first!=arg2.first ? arg1.first>arg2.first : arg1.second>arg2.second);
  }

  int _pairCompAscending(const INT_PAIR &arg1,const INT_PAIR &arg2){
    return (arg1.first!=arg2.first ? arg1.first<arg2.first : arg1.second<arg2.second);
  }

  template<class T> T rankAtomsByRank(const RDKit::ROMol &mol, const T &commAtms,
                                      bool ascending) {
    int natms = commAtms.size();
    INT_PAIR_VECT rankAid;
    rankAid.reserve(natms);
    T res;
    //res.reserve(natms);
    //RDKit::INT_VECT_CI ci;
    typename T::const_iterator ci;
    int rank;
    for (ci = commAtms.begin(); ci != commAtms.end(); ci++) {
      const RDKit::Atom *at=mol.getAtomWithIdx(*ci);
      if(at->hasProp("_CIPRank")){
        at->getProp("_CIPRank", rank);
      } else {
        rank = 2*mol.getNumAtoms()+(*ci);
      }
      rankAid.push_back(std::make_pair(rank, (*ci)));
    }
    if (ascending) {
      std::stable_sort(rankAid.begin(),rankAid.end(),_pairCompAscending);
    }
    else {
      std::stable_sort(rankAid.begin(),rankAid.end(),_pairCompDescending);
    }
    INT_PAIR_VECT_CI rai;
    for (rai = rankAid.begin(); rai != rankAid.end(); rai++) {
      res.push_back(rai->second);
    }

    return res;
  }

  template RDKit::INT_VECT rankAtomsByRank(const RDKit::ROMol &mol, 
                                           const RDKit::INT_VECT &commAtms,
                                           bool ascending);
  template RDKit::INT_DEQUE rankAtomsByRank(const RDKit::ROMol &mol, 
                                            const RDKit::INT_DEQUE &commAtms,
                                            bool ascending);
  template RDKit::INT_LIST rankAtomsByRank(const RDKit::ROMol &mol, 
                                           const RDKit::INT_LIST &commAtms,
                                           bool ascending);
}
