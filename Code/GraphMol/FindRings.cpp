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
#include <GraphMol/Rings.h>
#include <RDGeneral/RDLog.h>
#include <RDBoost/Exceptions.h>

#include <RDGeneral/utils.h>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <boost/cstdint.hpp>
#include <RDGeneral/hash/hash.hpp>

typedef std::set< boost::uint32_t > RINGINVAR_SET;
typedef RINGINVAR_SET::const_iterator RINGINVAR_SET_CI;
typedef std::vector< boost::uint32_t > RINGINVAR_VECT;

namespace RingUtils {
  using namespace RDKit;

  boost::uint32_t computeRingInvariant(INT_VECT ring,unsigned int nAtoms){
    std::sort(ring.begin(),ring.end());
    boost::uint32_t res=gboost::hash_range(ring.begin(),ring.end());
    return res;
  }

  void convertToBonds(const VECT_INT_VECT &res, VECT_INT_VECT &brings, const ROMol &mol) {
    for (VECT_INT_VECT_CI ri=res.begin(); ri!=res.end(); ++ri) {
      unsigned int rsiz = ri->size();
      INT_VECT bring(rsiz);
      for (unsigned int i = 0; i < (rsiz-1); i++) {
        const Bond *bnd=mol.getBondBetweenAtoms((*ri)[i],(*ri)[i+1]);
        if(!bnd) throw ValueErrorException("expected bond not found");
        bring[i]=bnd->getIdx();
      }
      // bond from last to first atom
      const Bond *bnd=mol.getBondBetweenAtoms((*ri)[rsiz-1],(*ri)[0]);
      if(!bnd) throw ValueErrorException("expected bond not found");

      bring[rsiz-1]=bnd->getIdx();
      brings.push_back(bring);
    }
  }

} // end of namespace RingUtils

namespace FindRings {
  using namespace RDKit;
  int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings,
                       boost::dynamic_bitset<> &activeBonds,
                       INT_VECT *forbidden=0);
  void trimBonds(unsigned int cand, const ROMol &tMol, INT_SET &changed,
                 INT_VECT &atomDegrees,boost::dynamic_bitset<> &activeBonds);
  void storeRingInfo(const ROMol &mol, const INT_VECT &ring) {
    INT_VECT bondIndices;
    INT_VECT_CI lastRai;
    for(INT_VECT_CI rai=ring.begin();rai != ring.end();rai++){
      if(rai!=ring.begin()){
        const Bond *bnd=mol.getBondBetweenAtoms(*rai,*lastRai);
        if(!bnd) throw ValueErrorException("expected bond not found");
        bondIndices.push_back(bnd->getIdx());
      }
      lastRai = rai;
    }
    const Bond *bnd=mol.getBondBetweenAtoms(*lastRai,*(ring.begin()));
    if(!bnd) throw ValueErrorException("expected bond not found");
    bondIndices.push_back(bnd->getIdx());
    mol.getRingInfo()->addRing(ring,bondIndices);
  }

  void storeRingsInfo(const ROMol &mol, const VECT_INT_VECT &rings) {
    for (VECT_INT_VECT_CI ri = rings.begin(); ri != rings.end(); ri++) {
      storeRingInfo(mol,*ri);
    }
  }

  void markUselessD2s(unsigned int root,const ROMol &tMol, boost::dynamic_bitset<> &forb,
                      const INT_VECT &atomDegrees, const boost::dynamic_bitset<> &activeBonds) {
    // recursive function to mark any degree 2 nodes that are already represnted 
    // by root for the purpose of finding smallest rings.
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = tMol.getAtomBonds(tMol.getAtomWithIdx(root));
    while(beg!=end){
      BOND_SPTR bond=tMol[*beg];
      ++beg;
      if(!activeBonds[bond->getIdx()]) continue;
      unsigned int oIdx=bond->getOtherAtomIdx(root);
      if(!forb[oIdx] && atomDegrees[oIdx]==2){
        forb[oIdx]=1;
        markUselessD2s(oIdx,tMol,forb,atomDegrees,activeBonds);
      }
    }
  }


  void pickD2Nodes(const ROMol &tMol, INT_VECT &d2nodes, const INT_VECT &currFrag,
                   const INT_VECT &atomDegrees, const boost::dynamic_bitset<> &activeBonds) {
    d2nodes.resize(0);

    // forb contains all d2 nodes, not just the ones we want to keep
    boost::dynamic_bitset<> forb(tMol.getNumAtoms());
    while (1) {
      int root = -1;
      for (INT_VECT_CI axci = currFrag.begin(); axci != currFrag.end(); ++axci) {
        if ( atomDegrees[*axci]==2 && !forb[*axci] ){
          root = (*axci);
          d2nodes.push_back(*axci);
          forb[*axci]=1;
          break;
        }
      }
      if (root == -1){
        break;
      }
      else {
        markUselessD2s(root, tMol, forb,atomDegrees,activeBonds);
      }  
    }
  }

#if 0
  typedef std::map<double, INT_VECT> DOUBLE_INT_VECT_MAP;
  typedef DOUBLE_INT_VECT_MAP::iterator DOUBLE_INT_VECT_MAP_I;
  typedef DOUBLE_INT_VECT_MAP::const_iterator DOUBLE_INT_VECT_MAP_CI;
#else
  typedef std::map<boost::uint32_t, INT_VECT> RINGINVAR_INT_VECT_MAP;
  typedef RINGINVAR_INT_VECT_MAP::iterator RINGINVAR_INT_VECT_MAP_I;
  typedef RINGINVAR_INT_VECT_MAP::const_iterator RINGINVAR_INT_VECT_MAP_CI;
#endif


  void findSSSRforDupCands(const ROMol &mol, VECT_INT_VECT &res, 
                           RINGINVAR_SET &invars, const INT_INT_VECT_MAP dupMap, 
                           const RINGINVAR_INT_VECT_MAP &dupD2Cands,
                           INT_VECT &atomDegrees, boost::dynamic_bitset<> activeBonds){
    for (RINGINVAR_INT_VECT_MAP_CI dvmi = dupD2Cands.begin();
         dvmi != dupD2Cands.end(); ++dvmi) {
      const INT_VECT &dupCands = dvmi->second;
      if (dupCands.size() > 1) {
        // we have duplicate candidates.
        VECT_INT_VECT nrings;
        unsigned int minSiz = static_cast<unsigned int>(MAX_INT);
        for (INT_VECT_CI dupi = dupCands.begin(); dupi != dupCands.end(); ++dupi) {
          // now break bonds for all the d2 nodes for that give the same rings as 
          // with (*dupi) and recompute smallest ring with (*dupi)
          INT_VECT atomDegreesCopy=atomDegrees;
          boost::dynamic_bitset<>  activeBondsCopy=activeBonds;
          INT_SET changed;
          INT_INT_VECT_MAP_CI dmci = dupMap.find(*dupi);
          for (INT_VECT_CI dni = dmci->second.begin(); dni != dmci->second.end(); ++dni) {
            trimBonds((*dni), mol, changed, atomDegreesCopy, activeBondsCopy);
          }

          // now find the smallest ring/s around (*dupi)
          VECT_INT_VECT srings;
          smallestRingsBfs(mol, (*dupi), srings, activeBondsCopy);
          for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end(); ++sri) {
            if (sri->size() < minSiz) {
              minSiz = sri->size();
            }
            nrings.push_back((*sri));
          }
        }
      
        for (VECT_INT_VECT_CI nri = nrings.begin(); nri != nrings.end(); ++nri) {
          if (nri->size() == minSiz) {
            boost::uint32_t invr = RingUtils::computeRingInvariant(*nri,mol.getNumAtoms());
            if (invars.find(invr) == invars.end()) {
              res.push_back((*nri));
              invars.insert(invr);
            }
          }
        } // end of loop over new rings found
      } // end if (dupCand.size() > 1) 
    } // end of loop over all set of duplicate candidates
  }
  
  struct compRingSize : public std::binary_function<INT_VECT,INT_VECT,bool> {
    bool operator()(const INT_VECT &v1, const INT_VECT &v2) const {
      return v1.size() < v2.size();
    }
  };

  void removeExtraRings(VECT_INT_VECT &res, unsigned int nexpt, const ROMol &mol) {
    // sort on size
    std::sort(res.begin(), res.end(), compRingSize());

    // change the rings from atom IDs to bondIds
    VECT_INT_VECT brings;
    RingUtils::convertToBonds(res, brings, mol);
    std::vector< boost::dynamic_bitset<> > bitBrings;
    bitBrings.reserve(brings.size());
    for(VECT_INT_VECT_CI vivi=brings.begin();vivi!=brings.end();++vivi){
      boost::dynamic_bitset<> lring(mol.getNumBonds());
      for(INT_VECT_CI ivi=vivi->begin();ivi!=vivi->end();++ivi){
        lring.set(*ivi);
      }
      bitBrings.push_back(lring);
    }

    boost::dynamic_bitset<> availRings(res.size());
    availRings.set();
    boost::dynamic_bitset<> keepRings(res.size());

    for(unsigned int i=0;i<res.size();++i){
      if(!availRings[i]) continue;
      keepRings.set(i);
      boost::dynamic_bitset<> munion(mol.getNumBonds());
      munion = bitBrings[i];
      for(unsigned int j=i+1;j<res.size();++j){
        if(!availRings[j]) continue;
        if(bitBrings[j].is_subset_of(munion)){
          availRings.set(j,0);
        } else {
          keepRings.set(j);
          availRings.set(j,0);
          munion |= bitBrings[j];
        }
      }
    }
    // remove the extra rings from res and store them on the molecule in case we wish 
    // symmetrize the SSSRs later
    VECT_INT_VECT extras;
    VECT_INT_VECT temp = res;
    res.resize(0);
    for (unsigned int i = 0; i < temp.size(); i++) {
      if(keepRings[i]){
        res.push_back(temp[i]);
      } else {
        extras.push_back(temp[i]);
      }
    }
  
    mol.setProp("extraRings", extras, true);
  }

  void findRingsD2nodes(const ROMol &tMol, VECT_INT_VECT &res, 
                        RINGINVAR_SET &invars, const INT_VECT &d2nodes,
                        INT_VECT &atomDegrees, boost::dynamic_bitset<> &activeBonds,
                        boost::dynamic_bitset<> &ringBonds,
                        boost::dynamic_bitset<> &ringAtoms
                        ) {
    // place to record any duplicate rings discovered from the current d2 nodes
    RINGINVAR_INT_VECT_MAP dupD2Cands;
    int cand;
    INT_VECT_CI d2i;

    INT_INT_VECT_MAP dupMap;
    // here is an example of molecule where the this scheme of finding other node that 
    // result in duplicates is necessary : C12=CON=C1C(C4)CC3CC2CC4C3
    // It would help to draw this molecule, and number the atoms but here is what happen
    //  - there are 6 d2 node - 1, 6, 7, 9, 11, 13
    //  - both 6 and 7 find the same ring (5,6,12,13,8,7) but we do not find the 7 membered ring
    //    (5,7,8,9,10,0,4)
    //  - similarly 9 and 11 find a duplicate ring (9,10,11,12,13)
    //  - when we move to 13 both the above duplicate rings are found
    //  - so we will keep track for each d2 all the other node that resulted in duplicate rings
    //  - the bonds to these nodes will be broken and we attempt to find a new ring, for e.g. by breaking
    //    bonds to 7 and 13, we will find a 7 membered ring with 6 (this is done in findSSSRforDupCands)
    std::map<int, RINGINVAR_VECT> nodeInvars;
    std::map<int, RINGINVAR_VECT>::const_iterator nici;
    DOUBLE_VECT_CI ici;
    for (d2i = d2nodes.begin(); d2i != d2nodes.end(); ++d2i) {
      cand = (*d2i);
      //std::cerr<<"    smallest rings bfs: "<<cand<<std::endl;
      VECT_INT_VECT srings;
      // we have to find all non duplicate possible smallest rings for each node
      smallestRingsBfs(tMol, cand, srings, activeBonds);
      for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end(); ++sri) {
        const INT_VECT &nring = (*sri);
        boost::uint32_t invr = RingUtils::computeRingInvariant(nring,tMol.getNumAtoms());
        if (invars.find(invr) == invars.end()) {
          res.push_back(nring);
          invars.insert(invr);
          for(unsigned int i=0;i<nring.size()-1;++i){
            unsigned int bIdx=tMol.getBondBetweenAtoms(nring[i],nring[i+1])->getIdx();
            ringBonds.set(bIdx);
            ringAtoms.set(nring[i]);
          }
          ringBonds.set(tMol.getBondBetweenAtoms(nring[0],nring[nring.size()-1])->getIdx());
          ringAtoms.set(nring[nring.size()-1]);
#if 0
          std::cerr<<"    res: "<<invr<<" | ";
          std::copy(nring.begin(),nring.end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;
#endif
        }

        nodeInvars[cand].push_back(invr);
        // check if this ring is duplicate with something else
        for (nici = nodeInvars.begin(); nici != nodeInvars.end(); nici++) {
          if (nici->first != cand) {
            if (std::find(nici->second.begin(), nici->second.end(), invr) != nici->second.end()) {
              // ok we discovered this ring via another node before
              // add that node as duplicate to this node and and vice versa
              dupMap[cand].push_back(nici->first);
              dupMap[nici->first].push_back(cand);
            }
          }
        }
        dupD2Cands[invr].push_back(cand);
      }
    
      // We don't want to trim the bonds connecting cand here - this can disrupt 
      // a second small ring. Here is an example SC(C3C1CC(C3)CC(C2S)(O)C1)2S
      // by trimming the bond connecting to atom #4 , we loose the smallest ring that 
      // contains atom #7. Issue 134
      //MolOps::trimBonds(cand, tMol, changed);
    }
  
    // now deal with any d2 nodes that resulted in duplicate rings before trimming their bonds.
    // it is possible that one of these nodes is involved a different small ring, that is not found 
    // because the first nodes has not be trimmed. Here is an example molecule: 
    // CC1=CC=C(C=C1)S(=O)(=O)O[CH]2[CH]3CO[CH](O3)[CH]4OC(C)(C)O[CH]24
    findSSSRforDupCands(tMol, res, invars, dupMap, dupD2Cands, atomDegrees, activeBonds);
  }

  void findRingsD3Node(const ROMol &tMol, VECT_INT_VECT &res, RINGINVAR_SET &invars, int cand,
                       INT_VECT &atomDegrees, boost::dynamic_bitset<> activeBonds ) {
    // this is brutal - we have no degree 2 nodes - find the first possible degree 3 node
    int nsmall;

    // We've got a degree three node. The goal of what follows is to find the
    // three rings in which it's involved, push those onto our results, and
    // then remove the node from consideration.  This will create a bunch of degree
    // 2 nodes, which we can then chew off the next time around the loop.
  
    // this part is a bit different from the Figueras algorithm 
    // here we try to find all the rings the rings that have a potential for contributing to
    // SSSR - i.e. we try to find 3 rings for this node. 
    // - each bond (that contributes to the degree 3 ) is allowed to participate in exactly
    //    two of these rings.
    // - also any rings that are included in already found rings are ignored
  
  
    // ASSUME: every connection from a degree three node at this point is a
    //         ring bond
    // REVIEW: Is this valid?
  
    // first find all smallest possible rings
    VECT_INT_VECT srings;
    nsmall = smallestRingsBfs(tMol, cand, srings, activeBonds);
  
    for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end(); ++sri) {
      const INT_VECT &nring = (*sri);
      boost::uint32_t invr = RingUtils::computeRingInvariant(nring,tMol.getNumAtoms());
      if (invars.find(invr) == invars.end()) {
        res.push_back(nring);
        invars.insert(invr);
      }
    }
  
    // if already found >3 rings we are done with this degree 3 node
    // if we found less than 3 we have to find other potential ring/s
    if (nsmall < 3) {
      int n1=-1, n2=-1, n3=-1;

      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = tMol.getAtomBonds(tMol.getAtomWithIdx(cand));
      while(beg!=end && !activeBonds[tMol[*beg]->getIdx()]) ++beg;
      CHECK_INVARIANT(beg!=end,"neighbor not found");
      n1 = tMol[*beg]->getOtherAtomIdx(cand);

      ++beg;
      while(beg!=end && !activeBonds[tMol[*beg]->getIdx()]) ++beg;
      CHECK_INVARIANT(beg!=end,"neighbor not found");
      n2 = tMol[*beg]->getOtherAtomIdx(cand);

      ++beg;
      while(beg!=end && !activeBonds[tMol[*beg]->getIdx()]) ++beg;
      CHECK_INVARIANT(beg!=end,"neighbor not found");
      n3 = tMol[*beg]->getOtherAtomIdx(cand);

      if (nsmall == 2) {
        // we found two rings find the third one
        // first find the neighbor that is common to the two ring we found so far
        int f;
      
        if ( (std::find(srings[0].begin(), srings[0].end(), n1) != srings[0].end()) 
             && (std::find(srings[1].begin(), srings[1].end(), n1) != srings[1].end()) ) {
          f = n1;
        }
        else if ( (std::find(srings[0].begin(), srings[0].end(), n2) != srings[0].end()) 
                  && (std::find(srings[1].begin(), srings[1].end(), n2) != srings[1].end()) ) {
          f = n2;
        }
        else if ( (std::find(srings[0].begin(), srings[0].end(), n3) != srings[0].end()) 
                  && (std::find(srings[1].begin(), srings[1].end(), n3) != srings[1].end()) ) {
          f = n3;
        }
      
        // now find the smallest possible ring that does not contain f
        VECT_INT_VECT trings;
        INT_VECT forb;
        forb.push_back(f);
        smallestRingsBfs(tMol, cand, trings, activeBonds,&forb);
        for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
          const INT_VECT &nring = (*sri);
          boost::uint32_t invr = RingUtils::computeRingInvariant(nring,tMol.getNumAtoms());

          if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
          }
        }
      } // doing degree 3 node  - end of 2 smallest rings found for cand
      if (nsmall == 1) {
        // we found 1 ring - we need to find two more that involve the 3rd neighbor
        int f1, f2;
        // Which of our three neighbors are in the small ring?
        //   these are f1 and f2
        if (std::find(srings[0].begin(), srings[0].end(), n1) == srings[0].end()) {
          f1 = n2, f2 = n3;
        }
        else if (std::find(srings[0].begin(), srings[0].end(), n2) == srings[0].end()) {
          f1 = n1; f2 = n3;
        }
        else if (std::find(srings[0].begin(), srings[0].end(), n3) == srings[0].end()) {
          f1 = n1; f2 = n2;
        }
      
        // now find two rings that include cand, one of these rings should include f1
        // and the other should include f2 
      
        // first ring with f1 and no f2
        VECT_INT_VECT trings;
        INT_VECT forb;
        forb.push_back(f2);
        smallestRingsBfs(tMol, cand, trings, activeBonds,&forb);
        for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
          const INT_VECT &nring = (*sri);
          boost::uint32_t invr = RingUtils::computeRingInvariant(nring,tMol.getNumAtoms());
          if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
          }
        }
      
        // next the ring with f2 and no f1
        trings.clear();
        forb.clear();
        forb.push_back(f1);
        smallestRingsBfs(tMol, cand, trings, activeBonds,&forb);
        for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
          const INT_VECT &nring = (*sri);
          boost::uint32_t invr = RingUtils::computeRingInvariant(nring,tMol.getNumAtoms());
          if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
          }
        }
      } // doing node of degree 3 - end of found only 1 smallest ring
    } // end of found less than 3 smallest ring for the degree 3 node
  }
 


  int greatestComFac(long curfac, long nfac) {
    long small;
    long large;
    long rem;

    // Determine which of the numbers is the larger, and which is the smaller
    large = (curfac > nfac) ? curfac : nfac;
    small = (curfac < nfac) ? curfac : nfac;

    // Keep looping until no remainder, as this means it is a factor of both
    while (small != 0){
      // Set the larger var to the smaller, and set the smaller to the remainder of (large / small)
      rem = (large % small);
      large = small;
      small = rem;
    }
  
    // By here nLarge will hold the largest common factor, so just return it
    return large;
  }



  /******************************************************************************
   * SUMMARY:
   *  remove the bond in the molecule that connect to the spcified atom
   *  
   * ARGUMENTS:
   *  cand - the node(atom) of interest
   *  tMol - molecule of interest
   *  changed - list of the atoms that are effected the bond removal
   *             this may be accumulated over multiple calls to trimBonds
   *             it basically forms a list of atom that need to be searched for 
   *             the next round of pruning
   *
   ******************************************************************************/
  void trimBonds(unsigned int cand, const ROMol &tMol, INT_SET &changed,
                 INT_VECT &atomDegrees,boost::dynamic_bitset<> &activeBonds) {

    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = tMol.getAtomBonds(tMol.getAtomWithIdx(cand));
    while(beg!=end){
      BOND_SPTR bond=tMol[*beg];
      ++beg;
      if(!activeBonds[bond->getIdx()]) continue;
      unsigned int oIdx=bond->getOtherAtomIdx(cand);
      if(atomDegrees[oIdx]<=2) changed.insert(oIdx);
      activeBonds[bond->getIdx()]=0;
      atomDegrees[oIdx]-=1;
      atomDegrees[cand]-=1;
    }
  }
    
  /*******************************************************************************
   * SUMMARY:
   *  this again is a modified version of the BFS algorihtm  in Figueras paper to find
   *  the smallest ring with a specified root atom.
   *    JCICS, Vol. 30, No. 5, 1996, 986-991
   *  The follwing are changes from the original algorithm
   *   - find all smallest rings around a node not just one
   *   - once can provided a list of node IDs that should not be include in the discovered rings
   * 
   * ARGUMENTS:
   *  mol - molecule of interest
   *  root - Atom ID of the node of interest
   *  rings - list of rings into which the results are entered
   *  forbidden - list of atoms ID that should be avoided
   * 
   * RETURNS:
   *  number of smallest rings found
   ***********************************************************************************/
  int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings,
                       boost::dynamic_bitset<> &activeBonds,
                       INT_VECT *forbidden) {
    // this function finds the smallest ring with the given root atom. 
    // if multiple smallest rings are found all of them are return
    // if any atoms are specified in the forbidden list, those atoms are avoided.

    // FIX: this should be number of atoms in the fragment (if it's required at all, see below)
    const int WHITE=0,GRAY=1,BLACK=2;
    INT_VECT done(mol.getNumAtoms(),WHITE);

    if (forbidden) {
      for (INT_VECT_CI dci = forbidden->begin(); dci != forbidden->end(); dci++) {
        done[*dci]=BLACK;
      }
    }

    // it would be "nicer" to use a map for this, but that ends up being too slow:
    VECT_INT_VECT atPaths(mol.getNumAtoms());
    INT_VECT rpath(1,root);
    atPaths[root] = rpath;

    std::deque<int> bfsq;
    bfsq.push_back(root);
    int curr=-1;
    unsigned int curSize=256;
    while (bfsq.size() > 0) {
      curr = bfsq.front();
      bfsq.pop_front();

      done[curr]=BLACK;

      INT_VECT &cpath = atPaths[curr];

      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = mol.getAtomBonds(mol.getAtomWithIdx(curr));
      while(beg!=end){
        BOND_SPTR bond=mol[*beg];
        ++beg;
        if(!activeBonds[bond->getIdx()]) continue;
        int nbrIdx=bond->getOtherAtomIdx(curr);
        if ((std::find(cpath.begin(), cpath.end(), nbrIdx) == cpath.end())
            && done[nbrIdx]!=BLACK ){
          // i.e. we are not at a node that is making up the current path
          // and we are not a node that has been completely explored before
          // (it has been a curr node before)

          // FIX: can we avoid this find by coloring atoms gray when they go into
          // the queue and just looking up the colors?
          if (done[nbrIdx]==WHITE) {
            // we have never been to this node before through via any path
            atPaths[nbrIdx] = cpath;
            atPaths[nbrIdx].push_back(nbrIdx);
            done[nbrIdx]=GRAY;
            bfsq.push_back(nbrIdx);
          } // end of found a untouched node
          else {
            // we have been here via a different path 
            // there is a potential for ring closure here
            INT_VECT npath = atPaths[nbrIdx];
            // make sure that the intersections of cpath and npath give exactl one 
            // element and that should be the root element for correct ring closure
            int id=-1;
            unsigned int com = 0;
            for (INT_VECT_CI ci = cpath.begin(); ci != cpath.end(); ++ci) {
              if (std::find(npath.begin(), npath.end(), (*ci)) != npath.end()) {
                com++;
                id = (*ci);
                if(id!=root) break;
              }
            } // end of found stuff in common with neighbor

            if (id == root){ // we found a ring     
              // make the ring
              INT_VECT ring = cpath;

              // remove the root node and attach the other half of the ring from npath 
              // reverse this piece so that the ring is traversed correctly

              // FIX: we're probably assured that root is the first node, so we can
              //  just pop it from the front
              npath.erase(std::remove(npath.begin(), npath.end(), root));

#ifndef WIN32
              ring.insert(ring.end(), npath.rbegin(), npath.rend());
#else // I <heart> MSVC++ v6
              std::reverse(npath.begin(), npath.end());
              ring.insert(ring.end(), npath.begin(), npath.end());
#endif
              if (ring.size() <= curSize) {
                curSize = ring.size();
                rings.push_back(ring) ;
              }
              else {
                // we are done with the smallest rings
                return rings.size();
              }
            } // end of found a ring
          } // end of we have seen this neighbor before
        } // end of nbrIdx not part of current path and not a done atom
      } // end of loop over neighbors of current atom 
    } // moving to the next node
    return rings.size(); // if we are here we should have found everything around the node
  }

  bool _atomSearchBFS(const ROMol &tMol,
                      unsigned int startAtomIdx,
                      unsigned int endAtomIdx,
                      boost::dynamic_bitset<> &ringAtoms,
                      INT_VECT &res,
                      RINGINVAR_SET &invars){
    res.clear();
    std::deque<INT_VECT> bfsq;

    INT_VECT tv;
    tv.push_back(startAtomIdx);
    bfsq.push_back(tv);
    while(!bfsq.empty()){
      tv = bfsq.front();
      bfsq.pop_front();

      unsigned int currAtomIdx=tv.back();
      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs) = tMol.getAtomNeighbors(tMol.getAtomWithIdx(currAtomIdx));
      while(nbrIdx!=endNbrs){
        if(*nbrIdx==endAtomIdx) {
          if(currAtomIdx!=startAtomIdx){
            tv.push_back(*nbrIdx);
            // make sure the ring we just found isn't already in our set
            // of rings (this was an extension of sf.net issue 249)
            boost::uint32_t invr = RingUtils::computeRingInvariant(tv,tMol.getNumAtoms());
            if (invars.find(invr) == invars.end()) {
              // we're done!
              res.resize(tv.size());
              std::copy(tv.begin(),tv.end(),res.begin());
              return true;
            }
          } else {
            // ignore this one
          }
        } else if(ringAtoms[*nbrIdx] && std::find(tv.begin(),tv.end(),*nbrIdx)==tv.end()){
          //} else if(ringAtoms[*nbrIdx]){
          INT_VECT nv(tv);
          nv.push_back(*nbrIdx);
          bfsq.push_back(nv);
        }
        ++nbrIdx;
      }
    }
    return false;
  }
                      
  bool findRingConnectingAtoms(const ROMol &tMol,
                               const Bond *bond,
                               VECT_INT_VECT &res,
                               RINGINVAR_SET &invars,
                               boost::dynamic_bitset<> &ringBonds,
                               boost::dynamic_bitset<> &ringAtoms ){
    PRECONDITION(bond,"bad bond");
    PRECONDITION(!ringBonds[bond->getIdx()],"not a ring bond");
    PRECONDITION(ringAtoms[bond->getBeginAtomIdx()],"not a ring atom");
    PRECONDITION(ringAtoms[bond->getEndAtomIdx()],"not a ring atom");

    INT_VECT nring;
    if(_atomSearchBFS(tMol,bond->getBeginAtomIdx(),bond->getEndAtomIdx(),
                      ringAtoms,nring,invars)){
      boost::uint32_t invr = RingUtils::computeRingInvariant(nring,tMol.getNumAtoms());
      if (invars.find(invr) == invars.end()) {
        res.push_back(nring);
        invars.insert(invr);
#if 0
        std::cerr<<"    local: "<<invr<<" | ";
        std::copy(nring.begin(),nring.end(),std::ostream_iterator<int>(std::cerr," "));
        std::cerr<<std::endl;
#endif
        for(unsigned int i=0;i<nring.size()-1;++i){
          unsigned int bIdx=tMol.getBondBetweenAtoms(nring[i],nring[i+1])->getIdx();
          ringBonds.set(bIdx);
          ringAtoms.set(nring[i]);
        }
        ringBonds.set(tMol.getBondBetweenAtoms(nring[0],nring[nring.size()-1])->getIdx());
        ringAtoms.set(nring[nring.size()-1]);
      }
    } else {
      return false;
    }
    return true;
  }


} // end of FindRings namespace

namespace RDKit {
  namespace MolOps {
    int findSSSR(const ROMol &mol, VECT_INT_VECT *res) {
      if (!res) {
        VECT_INT_VECT rings;
        return findSSSR(mol, rings);
      }
      else {
        return findSSSR(mol,(*res));
      }
    }

    int findSSSR(const ROMol &mol, VECT_INT_VECT &res) {
      res.resize(0);
      // check if SSSR's are already on the molecule
      if(mol.getRingInfo()->isInitialized()){
        res = mol.getRingInfo()->atomRings();
        return res.size();
      } else {
        mol.getRingInfo()->initialize();
      }

      RINGINVAR_SET invars;

      unsigned int nats = mol.getNumAtoms();
      boost::dynamic_bitset<> activeAtoms(nats);
      activeAtoms.set();
      int nbnds = mol.getNumBonds();
      boost::dynamic_bitset<> activeBonds(nbnds);
      activeBonds.set();

      // Zero-order bonds are not candidates for rings
      ROMol::EDGE_ITER firstB,lastB;
      boost::tie(firstB,lastB) = mol.getEdges();
      while(firstB!=lastB){
        BOND_SPTR bond = mol[*firstB];
        if(bond->getBondType()==Bond::ZERO) activeBonds[bond->getIdx()]=0;
        ++firstB;
      }

      
      boost::dynamic_bitset<> ringBonds(nbnds);
      boost::dynamic_bitset<> ringAtoms(nats);

      INT_VECT atomDegrees(nats);
      for(unsigned int i=0;i<nats;++i){
        const Atom *atom=mol.getAtomWithIdx(i);
        atomDegrees[i] = atom->getDegree();
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(atom);
        while(beg!=end){
          BOND_SPTR bond=mol[*beg];
          if(bond->getBondType()==Bond::ZERO) atomDegrees[i]--;
          ++beg;
        }
      }
      
      // find the number of fragments in the molecule - we will loop over them
      VECT_INT_VECT frags;
      INT_VECT curFrag;
      unsigned int nfrags = getMolFrags(mol, frags);
      for (unsigned int fi = 0; fi < nfrags; fi++) { // loop over the fragments in a molecule
        VECT_INT_VECT fragRes;
        curFrag = frags[fi];
    
        // the following is the list of atoms that are useful in the next round of trimming
        // basically atoms that become degree 0 or 1 because of bond removals
        // initialized with atoms of degrees 0 and 1
        INT_SET changed;
        for (INT_VECT_CI aidi = curFrag.begin(); aidi != curFrag.end(); aidi++) {
          unsigned int deg = atomDegrees[*aidi];
          if (deg<2) {
            changed.insert((*aidi));
          }
        }
    
        boost::dynamic_bitset<> doneAts(nats);
        unsigned int nAtomsDone=0;
        while (nAtomsDone < curFrag.size()) {
          //std::cerr<<" ndone: "<<nAtomsDone<<std::endl;
          //std::cerr<<" activeBonds: "<<activeBonds<<std::endl;
          //std::cerr<<"  done: ";
          // trim all bonds that connect to degree 0 and 1 atoms
          while (changed.size() > 0) {
            int cand = *(changed.begin());
            changed.erase(changed.begin());
            if (!doneAts[cand]){
              //std::cerr<<cand<<" ";
              doneAts.set(cand);
              ++nAtomsDone;
              FindRings::trimBonds(cand, mol, changed,atomDegrees,activeBonds);
            }
          }
          //std::cerr<<std::endl;
          //std::cerr<<"activeBonds2: "<<activeBonds<<std::endl;

          // all atoms left in the fragment should atleast have a degree >= 2
          // collect all the degree two nodes;
          INT_VECT d2nodes;
          FindRings::pickD2Nodes(mol, d2nodes, curFrag, atomDegrees, activeBonds);
#if 0
          std::cerr<<"d2nodes: ";
          std::copy(d2nodes.begin(),d2nodes.end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;          
#endif
          if (d2nodes.size() > 0) { // deal with the current degree two nodes
            // place to record any duplicate rings discovered from the current d2 nodes
            FindRings::findRingsD2nodes(mol, fragRes, invars, d2nodes, atomDegrees, activeBonds,
                                        ringBonds,ringAtoms);
#if 0
            std::cerr<<"  d2nodes post: ";
            std::copy(d2nodes.begin(),d2nodes.end(),std::ostream_iterator<int>(std::cerr," "));
            std::cerr<<std::endl;          
            std::cerr<<"  ring bonds: "<<ringBonds<<std::endl;
#endif
            INT_VECT_CI d2i;
            // trim after we have dealt with all the current d2 nodes, 
            for (d2i = d2nodes.begin(); d2i != d2nodes.end(); d2i++) {
              doneAts.set(*d2i);
              ++nAtomsDone;
              FindRings::trimBonds((*d2i), mol, changed, atomDegrees, activeBonds);
            }
          } // end of degree two nodes
          else if ( nAtomsDone < curFrag.size() ) { // now deal with higher degree nodes
            // this is brutal - we have no degree 2 nodes - find the first possible degree 3 node
            int cand = -1;
            for (INT_VECT_CI aidi = curFrag.begin(); aidi != curFrag.end(); aidi++) {
              unsigned int deg = atomDegrees[*aidi];
              if (deg == 3){ 
                cand = (*aidi); 
                break;
              }
            }
        
            // if we did not find a degree 3 node we are done
            // REVIEW:
            if (cand == -1) {
              break;
            }
            FindRings::findRingsD3Node(mol, fragRes, invars, cand, atomDegrees, activeBonds);
            doneAts.set(cand);
            ++nAtomsDone;
            FindRings::trimBonds(cand, mol, changed, atomDegrees, activeBonds); 
          } // done with degree 3 node
        } // done finding rings in this fragement

        // calculate the cyclomatic number for the fragment:
        unsigned int nbnds=0;
        for(ROMol::ConstBondIterator bndIt=mol.beginBonds();
            bndIt!=mol.endBonds();++bndIt){
          if(std::find(curFrag.begin(),curFrag.end(),(*bndIt)->getBeginAtomIdx())!=curFrag.end() &&
             std::find(curFrag.begin(),curFrag.end(),(*bndIt)->getEndAtomIdx())!=curFrag.end() &&
             (*bndIt)->getBondType()!=Bond::ZERO ) {
            ++nbnds;
          }
        }

#if 0
        std::cerr<<"\n\nFOUND:\n";
        for(VECT_INT_VECT::const_iterator iter=fragRes.begin();
            iter!=fragRes.end();++iter){
          std::copy(iter->begin(),iter->end(),std::ostream_iterator<int>(std::cerr," "));
          std::cerr<<std::endl;          
        }
#endif
        int nexpt = (nbnds - curFrag.size()+1);
        int ssiz = fragRes.size();

        // first check that we got at least the number of expected rings
        if(ssiz<nexpt){
          // Issue 3514824: in certain highly fused ring systems, the algorithm
          // above would miss rings.
          // for this fix to apply we have to have at least one non-ring bond
          // that terminates in ring atoms. Find those bonds:
          std::vector <const Bond *>  possibleBonds;
          for(unsigned int i=0;i<nbnds;++i){
            if(!ringBonds[i]){
              const Bond *bnd=mol.getBondWithIdx(i);
              if(ringAtoms[bnd->getBeginAtomIdx()] &&
                 ringAtoms[bnd->getEndAtomIdx()]){
                possibleBonds.push_back(bnd);
                break;
              }
            }
          }
          boost::dynamic_bitset<> deadBonds(mol.getNumBonds());
          while(possibleBonds.size()){
            bool ringFound=FindRings::findRingConnectingAtoms(mol,possibleBonds[0],
                                                              fragRes,invars,ringBonds,ringAtoms);
            if(!ringFound) deadBonds.set(possibleBonds[0]->getIdx(),1);
            possibleBonds.clear();
            // check if we need to repeat the process:
            for(unsigned int i=0;i<nbnds;++i){
              if(!ringBonds[i]){
                const Bond *bnd=mol.getBondWithIdx(i);
                if(!deadBonds[bnd->getIdx()] &&
                   ringAtoms[bnd->getBeginAtomIdx()] &&
                   ringAtoms[bnd->getEndAtomIdx()]
                   ){
                  possibleBonds.push_back(bnd);
                  break;
                }
              }
            }
          }
          ssiz = fragRes.size();
          if(ssiz<nexpt){
            throw ValueErrorException("could not find number of expected rings.");
          }
        }
        // if we have more than expected we need to do some cleanup
        // otherwise do som clean up work
        if (ssiz > nexpt) {
          FindRings::removeExtraRings(fragRes, nexpt, mol);
        }

        res.reserve(res.size()+fragRes.size());
        for(VECT_INT_VECT::const_iterator iter=fragRes.begin();
            iter!=fragRes.end();++iter){
          res.push_back(*iter);
        }
      } // done with all fragments
  
      FindRings::storeRingsInfo(mol,res);

      // update the ring memberships of atoms and bonds in the molecule:
      // store the SSSR rings on the the molecule as a property
      // we will ignore any existing SSSRs ont eh molecule - simply overwrite
      return res.size();
    }
  
    int symmetrizeSSSR(ROMol &mol) {
      VECT_INT_VECT tmp;
      return symmetrizeSSSR(mol,tmp);
    };

    int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res) {
      res.clear();res.resize(0);
      unsigned int nsssr;
      VECT_INT_VECT sssrs;

      // FIX: need to set flag here the symmetrization has been done in order to avoid
      //    repeating this work
      if(!mol.getRingInfo()->isInitialized()){
        nsssr = findSSSR(mol, sssrs);
      } else {
        sssrs = mol.getRingInfo()->atomRings();
        nsssr = sssrs.size();
      }

      VECT_INT_VECT_CI srci;
      INT_VECT copr;
      for (srci = sssrs.begin(); srci != sssrs.end(); srci++) {
        copr = (*srci);
        res.push_back(copr);
      }
 
      // now check if there are any extra rings on the molecule 
      if (!mol.hasProp("extraRings")) {
        // no extra rings nothign to be done
        return res.size();
      }
      const VECT_INT_VECT &extras=mol.getProp<VECT_INT_VECT>("extraRings");


      // convert the rings to bond ids
      VECT_INT_VECT bsrs, bextra;
      RingUtils::convertToBonds(sssrs, bsrs, mol);
      RingUtils::convertToBonds(extras, bextra, mol);
      INT_VECT munion, nunion, symids;
      Union(bsrs, munion);
      INT_VECT sr, exr;
      INT_VECT_CI eri;
      unsigned int eid, srid, ssiz;
      unsigned int next = bextra.size();
      // now the trick is the following
      // we will replace each ring of size ssiz from the SSSR with 
      // one of the same size rings in the extras. Compute the union of of the new set
      // if all the union elements of the new set if same as munion we found a  symmetric ring
      for (srid = 0; srid < nsssr; srid++) {
        sr = bsrs[srid];
        ssiz = sr.size();
        INT_VECT exrid;
        exrid.push_back(srid);
        Union(bsrs, nunion, &exrid);
        for (eid = 0; eid < next; eid++) {
          // if we already added this ring continue
          // FIX: if the ring has already been added,it probably shouldn't be
          // in the list at all?  Is this perhaps the most efficient way?
          if (std::find(symids.begin(), symids.end(), static_cast<int>(eid)) != symids.end()){
            continue;
          }
          exr = bextra[eid];
          if (ssiz == exr.size()) {
            INT_VECT eunion;
            Union(nunion, exr, eunion);
            // now check if the eunion is same as the original union from the SSSRs
            if (eunion.size() == munion.size()) {
              //we found a symmetric ring
              symids.push_back(eid);
            }
          }
        }
      }

      // add the symmertic rings
      for (eri = symids.begin(); eri != symids.end(); eri++) {
        exr = extras[*eri];
        res.push_back(exr);
        FindRings::storeRingInfo(mol, exr);
      }
      if (mol.hasProp("extraRings")) {
        mol.clearProp("extraRings");
      }
      return res.size();
    }

    namespace {
      void _DFS(const ROMol &mol,const Atom *atom,INT_VECT &atomColors,std::vector<const Atom *> &traversalOrder,
                VECT_INT_VECT &res,const Atom *fromAtom=0){
        //std::cerr<<"  dfs: "<<atom->getIdx()<<" from "<<(fromAtom?fromAtom->getIdx():-1)<<std::endl;
        PRECONDITION(atom,"bad atom");
        PRECONDITION(atomColors[atom->getIdx()]==0,"bad color");
        atomColors[atom->getIdx()]=1;
        traversalOrder.push_back(atom);


        ROMol::ADJ_ITER nbrIter,endNbrs;
        boost::tie(nbrIter,endNbrs) = mol.getAtomNeighbors(atom);
        while(nbrIter!=endNbrs){
          const Atom *nbr=mol[*nbrIter].get();
          unsigned int nbrIdx=nbr->getIdx();
          //std::cerr<<"   "<<atom->getIdx()<<"       consider: "<<nbrIdx<<"  "<<atomColors[nbrIdx]<<std::endl;
          if(atomColors[nbrIdx]==0){
            if(nbr->getDegree()<2){
              atomColors[nbr->getIdx()]=2;
            } else {
              _DFS(mol,nbr,atomColors,traversalOrder,res,atom);
            }
          } else if(atomColors[nbrIdx]==1){
            if(fromAtom && nbrIdx!=fromAtom->getIdx()){
              INT_VECT cycle;
              std::vector<const Atom *>::reverse_iterator lastElem=std::find(traversalOrder.rbegin(),traversalOrder.rend(),atom);
              for(std::vector<const Atom *>::reverse_iterator rIt=lastElem;//traversalOrder.rbegin();
                  rIt!=traversalOrder.rend() && (*rIt)->getIdx()!=nbrIdx;
                  ++rIt){
                cycle.push_back((*rIt)->getIdx());
              }
              cycle.push_back(nbrIdx);
              res.push_back(cycle);
              //std::cerr<<"    cycle from "<<atom->getIdx()<<" :";
              //std::copy(cycle.begin(),cycle.end(),std::ostream_iterator<int>(std::cerr," "));
              //std::cerr<<std::endl;
            }
          }
          ++nbrIter;
        }
        atomColors[atom->getIdx()]=2;
        traversalOrder.pop_back();
        //std::cerr<<"  done "<<atom->getIdx()<<std::endl;
      }
    } // end of anonymous namespace
    void fastFindRings(const ROMol &mol){
      //std::cerr<<"ffr"<<std::endl;
      VECT_INT_VECT res;
      res.resize(0);
      // check if SSSR's are already on the molecule
      if(mol.getRingInfo()->isInitialized()){
        return;
      } else {
        mol.getRingInfo()->initialize();
      }

      unsigned int nats = mol.getNumAtoms();

      INT_VECT atomColors(nats,0);

      for(unsigned int i=0;i<nats;++i){
        if(atomColors[i]) continue;
        if(mol.getAtomWithIdx(i)->getDegree()<2){
          atomColors[i]=2;
          continue;
        }
        std::vector<const Atom *> traversalOrder;
        _DFS(mol,mol.getAtomWithIdx(i),atomColors,traversalOrder,res);
      }
  
      FindRings::storeRingsInfo(mol,res);
    }

    
  }// end of MolOps namespace  

} // end of RDKit namespace
