// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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

typedef std::set<double> DOUBLE_SET;
typedef DOUBLE_SET::const_iterator DOUBLE_SET_CI;

namespace RingUtils {
  using namespace RDKit;

  void convertToBonds(const VECT_INT_VECT &res, VECT_INT_VECT &brings, const ROMol &mol) {
    for (VECT_INT_VECT_CI ri=res.begin(); ri!=res.end(); ++ri) {
      unsigned int rsiz = ri->size();
      INT_VECT bring(rsiz);
      for (unsigned int i = 0; i < (rsiz-1); i++) {
        bring[i]=mol.getBondBetweenAtoms((*ri)[i],(*ri)[i+1])->getIdx();
      }
      // bond from last to first atom
      bring[rsiz-1]=mol.getBondBetweenAtoms((*ri)[rsiz-1], (*ri)[0])->getIdx();
      brings.push_back(bring);
    }
  }

} // end of namespace RingUtils

namespace FindRings {
  using namespace RDKit;
  int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings,
                       INT_VECT *forbidden=0);
  void trimBonds(int cand, RWMol &tMol, INT_SET &changed);
  void storeRingInfo(const ROMol &mol, const INT_VECT &ring) {
    INT_VECT bondIndices;
    INT_VECT_CI lastRai;
    for(INT_VECT_CI rai=ring.begin();rai != ring.end();rai++){
      if(rai!=ring.begin()){
        bondIndices.push_back(mol.getBondBetweenAtoms(*rai,*lastRai)->getIdx());
      }
      lastRai = rai;
    }
    bondIndices.push_back(mol.getBondBetweenAtoms(*lastRai,*(ring.begin()))->getIdx());
    mol.getRingInfo()->addRing(ring,bondIndices);
  }

  void storeRingsInfo(const ROMol &mol, const VECT_INT_VECT &rings) {
    for (VECT_INT_VECT_CI ri = rings.begin(); ri != rings.end(); ri++) {
      storeRingInfo(mol,*ri);
    }
  }

  void markUselessD2s(int root,const ROMol &tMol, boost::dynamic_bitset<> &forb) {
    // recursive function to mark any degree 2 nodes that are already represnted 
    // by root for the purpose of finding smallest rings.
    ROMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = tMol.getAtomNeighbors(tMol.getAtomWithIdx(root));
    while(nbrIdx != endNbrs) {
      if (!forb[*nbrIdx]){
        const Atom *at = tMol.getAtomWithIdx(*nbrIdx);
        if (at->getDegree() == 2) {
          forb[*nbrIdx]=1;
          markUselessD2s(*nbrIdx, tMol, forb);
        }  
      }
      ++nbrIdx;
    }
  }


  void pickD2Nodes(const ROMol &tMol, INT_VECT &d2nodes, const INT_VECT &currFrag) { 
    d2nodes.resize(0);

    // forb contains all d2 nodes, not just the ones we want to keep
    boost::dynamic_bitset<> forb(tMol.getNumAtoms());
    while (1) {
      int root = -1;
      for (INT_VECT_CI axci = currFrag.begin(); axci != currFrag.end(); ++axci) {
        const Atom *at = tMol.getAtomWithIdx(*axci);
        if ( (at->getDegree() == 2 ) && !forb[*axci] ){
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
        markUselessD2s(root, tMol, forb);
      }  
    }
  }

  typedef std::map<double, INT_VECT> DOUBLE_INT_VECT_MAP;
  typedef DOUBLE_INT_VECT_MAP::iterator DOUBLE_INT_VECT_MAP_I;
  typedef DOUBLE_INT_VECT_MAP::const_iterator DOUBLE_INT_VECT_MAP_CI;

  void findSSSRforDupCands(const RWMol &mol, VECT_INT_VECT &res, 
                           DOUBLE_SET &invars, const INT_INT_VECT_MAP dupMap, 
                           const DOUBLE_INT_VECT_MAP &dupD2Cands) {
    for (DOUBLE_INT_VECT_MAP_CI dvmi = dupD2Cands.begin();
         dvmi != dupD2Cands.end(); ++dvmi) {
      const INT_VECT &dupCands = dvmi->second;
      if (dupCands.size() > 1) {
        // we have duplicate candidates.
        VECT_INT_VECT nrings;
        unsigned int minSiz = static_cast<unsigned int>(MAX_INT);
        for (INT_VECT_CI dupi = dupCands.begin(); dupi != dupCands.end(); ++dupi) {
          // now break bonds for all the d2 nodes for that give the same rings as 
          // with (*dupi) and recompute smallest ring with (*dupi)
          // copy the molecule so that we can break teh bonds
          RWMol tMol(mol,true);
          INT_SET changed;
          INT_INT_VECT_MAP_CI dmci = dupMap.find(*dupi);
          for (INT_VECT_CI dni = dmci->second.begin(); dni != dmci->second.end(); ++dni) {
            trimBonds((*dni), tMol, changed);
          }

          // now find the smallest ring/s around (*dupi)
          VECT_INT_VECT srings;
          smallestRingsBfs(tMol, (*dupi), srings);
          for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end(); ++sri) {
            if (sri->size() < minSiz) {
              minSiz = sri->size();
            }
            nrings.push_back((*sri));
          }
        }
      
        for (VECT_INT_VECT_CI nri = nrings.begin(); nri != nrings.end(); ++nri) {
          if (nri->size() == minSiz) {
            double invr = computeIntVectPrimesProduct((*nri));
            if (invars.find(invr) == invars.end()) {
              res.push_back((*nri));
              invars.insert(invr);
            }
          }
        } // end of loop over new rings found
      } // end if (dupCand.size() > 1) 
    } // end of loop over all set of duplicate candidates
  }

  bool compRingSize(const INT_VECT &ring1, const INT_VECT &ring2) {
    return (ring1.size() < ring2.size());
  }

  void removeExtraRings(VECT_INT_VECT &res, unsigned int nexpt, const ROMol &mol) {
    // convert each ring in res from a list of atom ids to list of bonds id

    // sort on size
    std::sort(res.begin(), res.end(), compRingSize);

    // change the rings from atom IDs to bondIds
    VECT_INT_VECT brings;
    RingUtils::convertToBonds(res, brings, mol);

    unsigned int tot = res.size();

    // the algorithm here is quite straightforward
    // - take the union of bonds from all the rings
    // - since we know how many SSSRs to expect, take the union of 
    //   subsets of expected size.
    // - if the union of bonds from the subset of rings give the entire union we
    //   have the SSSR set
  
    // find the overall union
    boost::dynamic_bitset<> munion(mol.getNumBonds());
    for (VECT_INT_VECT_I ri = brings.begin(); ri != brings.end(); ++ri) {
      for (INT_VECT_CI mi = ri->begin(); mi != ri->end(); ++mi) {
        munion[*mi]=1;
      }
    }

    INT_VECT comb(nexpt);
    for (unsigned int i = 0; i < nexpt; i++) {
      comb[i]=i;
    }

    bool found = false;
    int pos;
    while (!found) {
      boost::dynamic_bitset<> cunion(mol.getNumBonds());
      found = true;
      for (unsigned int i = 0; i < nexpt; ++i) {
        INT_VECT bring = brings[comb[i]];
        for (unsigned int j = 0; j < bring.size(); ++j) {
          cunion[bring[j]]=1;
        }
      }
      if (cunion.count() < munion.count()) {
        pos = nextCombination(comb, tot);
        CHECK_INVARIANT(pos >= 0,""); // we couldn't have run through all the combinations without removing any rings
        found = false;
      }
    }

    // remove the extra rings from res and store them on the molecule in case we wish 
    // symmetrize the SSSRs later
    VECT_INT_VECT extras;
    VECT_INT_VECT temp = res;
    res.resize(0);
    for (unsigned int i = 0; i < temp.size(); i++) {
      if (std::find(comb.begin(), comb.end(), static_cast<int>(i)) != comb.end()) {
        res.push_back(temp[i]);
      } else {
        extras.push_back(temp[i]);
      }
    }
    // store the extra rings on teh molecule for later use like
    // symmetrizing the SSSRs
    mol.setProp("extraRings", extras, true);
  }

  void findRingsD2nodes(RWMol &tMol, VECT_INT_VECT &res, 
                        DOUBLE_SET &invars, const INT_VECT &d2nodes) {
    // place to record any duplicate rings discovered from the current d2 nodes
    DOUBLE_INT_VECT_MAP dupD2Cands;
    int cand, nsmall;
    double invr;
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
    std::map<int, DOUBLE_VECT> nodeInvars;
    std::map<int, DOUBLE_VECT>::const_iterator nici;
    DOUBLE_VECT_CI ici;
    for (d2i = d2nodes.begin(); d2i != d2nodes.end(); d2i++) {
      cand = (*d2i);
    
      VECT_INT_VECT srings;
    
      VECT_INT_VECT_CI sri;
    
      // we have to find all non duplicate possible smallest rings for each node
      //rsiz = MolOps::findSmallestRing(cand, tMol, ring, invars);
      nsmall = smallestRingsBfs(tMol, cand, srings);
      for (sri = srings.begin(); sri != srings.end(); sri++) {
        const INT_VECT &nring = (*sri);
        invr = computeIntVectPrimesProduct(nring);
        if (invars.find(invr) == invars.end()) {
          res.push_back(nring);
          invars.insert(invr);
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
    findSSSRforDupCands(tMol, res, invars, dupMap, dupD2Cands);
  }

  void findRingsD3Node(RWMol &tMol, VECT_INT_VECT &res, DOUBLE_SET &invars, int cand) {
  
    // this is brutal - we have no degree 2 nodes - find the first possible degree 3 node
    int nsmall;

    double invr; 
  
    // We've got a degree three node. The goal of what follows is to find the
    // three rings in which it's involved, push those onto our results, and
    // then remove the node from consideration.  This will create a bunch of degree
    // 2 nodes, which we can then chew off the next time around the loop.
  
    // this part is a bit different fromt he Figueras algorithm 
    // here we try to find all the rings the rings that have a potential for contributing to
    // SSSR - i.e. we try to find 3 rings for this node. 
    // - each bond (that contributres to the degree 3 ) is allowed to participate in exactly
    //    two of these rings.
    // - also any rings that are inclusive in alsready found rings are ingnored
  
  
    // ASSUME: every connection from a degree three node at this point is a
    //         ring bond
    // REVIEW: Is this valid?
  
    // first find all smallest possible rings
    VECT_INT_VECT srings;
    nsmall = smallestRingsBfs(tMol, cand, srings);
  
    for (VECT_INT_VECT_CI sri = srings.begin(); sri != srings.end(); ++sri) {
      const INT_VECT &nring = (*sri);
      invr = computeIntVectPrimesProduct(nring); 
    
      if (invars.find(invr) == invars.end()) {
        res.push_back(nring);
        invars.insert(invr);
      }
    }
  
    // if already found >3 rings we are done with this degree 3 node
    // if we found less than 3 we have to find other potential ring/s
    if (nsmall < 3) {
      int n1, n2, n3;
      ROMol::ADJ_ITER nbrIdx,endNbrs;
      boost::tie(nbrIdx,endNbrs) = tMol.getAtomNeighbors(tMol.getAtomWithIdx(cand));
      n1 = (*nbrIdx); nbrIdx++;
      n2 = (*nbrIdx); nbrIdx++;
      n3 = (*nbrIdx);
    
      if (nsmall == 2) {
        // we found two rings find the third one
        // first find the neighbor that is common to the two ring we found so far
        int f;
      
        if ( (std::find(srings[0].begin(), srings[0].end(), n1) != srings[0].end()) 
             & (std::find(srings[1].begin(), srings[1].end(), n1) != srings[1].end()) ) {
          f = n1;
        }
        else if ( (std::find(srings[0].begin(), srings[0].end(), n2) != srings[0].end()) 
                  & (std::find(srings[1].begin(), srings[1].end(), n2) != srings[1].end()) ) {
          f = n2;
        }
        else if ( (std::find(srings[0].begin(), srings[0].end(), n3) != srings[0].end()) 
                  & (std::find(srings[1].begin(), srings[1].end(), n3) != srings[1].end()) ) {
          f = n3;
        }
      
        // now find the smallest possible ring that does not contain f
        VECT_INT_VECT trings;
        INT_VECT forb;
        forb.push_back(f);
        smallestRingsBfs(tMol, cand, trings, &forb);
        for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
          const INT_VECT &nring = (*sri);
          invr = computeIntVectPrimesProduct(nring); 
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
        int nrngs = smallestRingsBfs(tMol, cand, trings, &forb);
        for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
          const INT_VECT &nring = (*sri);
          invr = computeIntVectPrimesProduct(nring); 
          if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
          }
        }
      
        // next the ring with f2 and no f1
        trings.clear();
        forb.clear();
        forb.push_back(f1);
        nrngs = smallestRingsBfs(tMol, cand, trings, &forb);
        for (VECT_INT_VECT_CI sri = trings.begin(); sri != trings.end(); ++sri) {
          const INT_VECT &nring = (*sri);
          invr = computeIntVectPrimesProduct(nring); 
          if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
          }
        }
      } // doing node of degree 3 - end of found only 1 smallest ring
    } // end of found less than 3 smallest ring for the degree 3 node
    //doneAts.push_back(cand);
    //MolOps::trimBonds(cand, tMol, changed);
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
  void trimBonds(int cand, RWMol &tMol, INT_SET &changed) {
    // basically loop over the bonds for cand and mark the neighbors if any of them after 
    // bond removal become degree 1 or 0
    ROMol::ADJ_ITER nbrIdx,endNbrs;
    boost::tie(nbrIdx,endNbrs) = tMol.getAtomNeighbors(tMol.getAtomWithIdx(cand));
    INT_VECT neighs;
    while (nbrIdx != endNbrs) {
      neighs.push_back(*nbrIdx);
      ++nbrIdx;
    }

    for (INT_VECT_CI nci = neighs.begin(); nci != neighs.end(); ++nci) {
      Atom *nat = tMol.getAtomWithIdx(*nci);
    
      if (nat->getDegree() <= 2) {
        changed.insert(*nci);
      }
      tMol.removeBond(cand, (*nci)); 
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
  int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings, INT_VECT *forbidden) {
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
      ROMol::ADJ_ITER nbr,endNbrs;
      boost::tie(nbr,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(curr));
      while(nbr != endNbrs) {
        int nbrIdx=(*nbr);
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
              // FIX: we can set the reserve on ring here to minimize reallocs

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
        ++nbr;
      } // end of loop over neighbors of current atom 
    } // moving to the next node
    return rings.size(); // if we are here we should have founf everything around the node
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

      DOUBLE_SET invars;
      //DOUBLE_VECT invars;

      // make a copy of the molecule that we can chop around
      // mind you we will never remove atoms to avoid numbering issues
      // only bonds are removed
      RWMol tMol(mol,true);

      int nats = tMol.getNumAtoms();
      int nbnds = tMol.getNumBonds();

      // find the number of fragments in the molecule - we will loop over them
      VECT_INT_VECT frags;
      INT_VECT curFrag;
      int fi, nfrags = getMolFrags(tMol, frags);
      for (fi = 0; fi < nfrags; fi++) { // loop over the fragments in a molecule
    
        curFrag = frags[fi];
    
        // the following is the list of atoms that are useful in the next round of trimming
        // basically atoms that become degree 0 or 1 because of bond removals
        // initialized with atoms of degrees 0 and 1
        INT_VECT doneAts; // atoms that we already dealt with int he fragment
        INT_SET changed;
        INT_VECT_CI aidi;
        int deg, cand;
        for (aidi = curFrag.begin(); aidi != curFrag.end(); aidi++) {
          deg = tMol.getAtomWithIdx((*aidi))->getDegree();
          if ((deg == 0) || (deg == 1)) {
            changed.insert((*aidi));
          }
        }
    
        while (doneAts.size() < curFrag.size()) {
          //trim all bonds that connect to degree 0 and 1 bonds
          while (changed.size() > 0) {
            cand = *(changed.begin());
            changed.erase(changed.begin());
            if (std::find(doneAts.begin(), doneAts.end(), cand) == doneAts.end()) {
              doneAts.push_back(cand);
              FindRings::trimBonds(cand, tMol, changed);
            }
          }

          // all atoms left in the fragment should atleast have a degree >= 2
          // collect all the degree two nodes;
          INT_VECT d2nodes;
    
          // pick all the d2nodes from the current fragment
          FindRings::pickD2Nodes(tMol, d2nodes, curFrag);
          
          if (d2nodes.size() > 0) { // deal with the current degree two nodes
            // place to record any duplicate rings discovered from the current d2 nodes
            FindRings::findRingsD2nodes(tMol, res, invars, d2nodes);

            INT_VECT_CI d2i;
            // trim after we have dealt with all the current d2 nodes, 
            for (d2i = d2nodes.begin(); d2i != d2nodes.end(); d2i++) {
              doneAts.push_back((*d2i));
              FindRings::trimBonds((*d2i), tMol, changed);
            }
          } // end of degree two nodes
    
          else if ( doneAts.size() < curFrag.size() ) { // now deal with higher degree nodes
        
            //INT_VECT ring;
            // this is brutal - we have no degree 2 nodes - find the first possible degree 3 node
            cand = -1;
            for (aidi = curFrag.begin(); aidi != curFrag.end(); aidi++) {
              deg = tMol.getAtomWithIdx((*aidi))->getDegree();
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
            FindRings::findRingsD3Node(tMol, res, invars, cand);
            doneAts.push_back(cand);
            FindRings::trimBonds(cand, tMol, changed); 
          } // done with degree 3 node
        } // done finding rings in this fragement
      } // done with all fragments
  
      // calculate the Frere-Jacque number
      int nexpt = (nbnds - nats + nfrags);
      int ssiz = res.size();
      // first check that we got more than or equal to the number of expected rings
      if(ssiz<nexpt){
        throw ValueErrorException("could not find number of expected rings.");
      }
  
      // if we have more than expected we need to do some cleanup
      // otherwise do som celan up work
      if (ssiz > nexpt) {
        FindRings::removeExtraRings(res, nexpt, mol);
      }

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
      VECT_INT_VECT extras;
      if (!mol.hasProp("extraRings")) {
        // no extra rings nothign to be done
        return res.size();
      }
      else {
        mol.getProp("extraRings", extras);
      }

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

  }// end of MolOps namespace  

} // end of RDKit namespace
