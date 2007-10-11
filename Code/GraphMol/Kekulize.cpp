// $Id$
//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Canon.h>
#include <GraphMol/Rings.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/RDLog.h>
#include <boost/dynamic_bitset.hpp>


// end of namespace Kekulize
namespace RDKit {
  // Local utility namespace 
  namespace {
    void backTrack(RWMol &mol, 
                   INT_INT_DEQ_MAP &options,
                   int lastOpt, 
                   INT_VECT &done, 
                   INT_DEQUE &aqueue, 
                   boost::dynamic_bitset<> &dBndCands, 
                   boost::dynamic_bitset<> &dBndAdds) {
      // so we made a wrong turn at the lastOpt
      //remove on done list that comes after the lastOpt including itself

      INT_VECT_CI ei = std::find(done.begin(), done.end(), lastOpt);
      INT_VECT tdone;
      for (INT_VECT_CI ci = done.begin(); ci != ei; ci++) {
        tdone.push_back(*ci);
      }

      INT_VECT_CRI eri = std::find(done.rbegin(), done.rend(), lastOpt);
      ++eri;
      // and push them back onto the stack
      for (INT_VECT_CRI ri = done.rbegin(); ri != eri; ++ri) {
        aqueue.push_front(*ri);
      }
    
      // remove any double bonds that were add since we passed through lastOpt
      int aid1, aid2;
      Bond *bnd;
      int bi, nbnds = mol.getNumBonds();
      for (bi = 0; bi < nbnds; bi++) {
        if (dBndAdds[bi]) {
          bnd = mol.getBondWithIdx(bi);
          aid1 = bnd->getBeginAtomIdx();
          aid2 = bnd->getEndAtomIdx();
          // if one of these atoms has been dealt with before lastOpt
          // we don't have to chnage the double bond addition
          if ( (std::find(tdone.begin(), tdone.end(), aid1) == tdone.end()) &&
               (std::find(tdone.begin(), tdone.end(), aid2) == tdone.end()) ) {
            // otherwise strip the double bond and set it back to single
            // and add the atoms to candidate for double bonds
            dBndAdds[bi]=0;
            bnd->setBondType(Bond::SINGLE);
            dBndCands[aid1]=1;
            dBndCands[aid2]=1;
          }
        }
      }
      done = tdone;
    }


    void markDbondCands(RWMol &mol, const INT_VECT &allAtms,
                        boost::dynamic_bitset<> &dBndCands, INT_VECT &done){
      // ok this function does more than mark atoms that are candidates for 
      // double bonds during kekulization
      // - check that an non aromatic atom does not have any aromatic bonds
      // - marks all aromatic bonds to single bonds
      // - marks atoms that can take a double bond 
      boost::dynamic_bitset<> makeSingle(mol.getNumBonds());
      RWMol::BondIterator bi;

      Atom *at;
      int tbo, sbo, dv, chrg;
      INT_VECT_CI adx;
      RWMol::OEDGE_ITER beg,end;
      RWMol::GRAPH_MOL_BOND_PMAP::type pMap = mol.getBondPMap();
      for (adx = allAtms.begin(); adx != allAtms.end(); adx++) {
        // if this atom is not aromatic, don't have to kekulize it
        at = mol.getAtomWithIdx(*adx);
        if (!at->getIsAromatic()) {
          done.push_back(*adx);
          // make sure all the bonds on this atom are also non aromatic
          // i.e. can't have aromatic bond onto a non-aromatic atom
          boost::tie(beg,end) = mol.getAtomBonds(at);
        
          while (beg != end) {
            // ok we can't have an aromatic atom
            if (pMap[*beg]->getIsAromatic()) {
              std::ostringstream errout;
              errout << "Aromatic bonds on non aromatic atom " << at->getIdx();
              std::string msg = errout.str();
              BOOST_LOG(rdErrorLog) << msg << std::endl;
              throw MolSanitizeException(msg);
            }
            beg++;
          }
          continue;
        }
      
        sbo = 0;
        dv = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
        chrg = at->getFormalCharge();
        dv += chrg;
        tbo = at->getExplicitValence() + at->getImplicitValence();
        const UINT_VECT &valList =
          PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
        unsigned int vi = 1;
      
        while ((tbo > dv) && (vi < valList.size())) {
          dv = valList[vi] + chrg;
          vi++;
        }
        
        boost::tie(beg,end) = mol.getAtomBonds(at);
        while (beg != end) {
          const Bond *bond=pMap[*beg];
          if (bond->getIsAromatic()) {
            sbo++;
            // mark this bond to be marked single later 
            // we don't want to do right now because it can srew-up the 
            // valence calculation to determine the number of hydrogens below
            makeSingle[bond->getIdx()]=1;
          }
          else {
            sbo += (int)bond->getValenceContrib(at);
          }
          beg++;
        }
        sbo +=  at->getTotalNumHs();
      
        if (dv == (sbo + 1)) {
          dBndCands[*adx]=1;
        } 
      }// loop over all atoms in the fused system

      // now turn all the aromatic bond in this fused system to single
      for (bi = mol.beginBonds(); bi != mol.endBonds(); bi++) {
        if (makeSingle[(*bi)->getIdx()] ) {
          (*bi)->setBondType(Bond::SINGLE);
        }
      }
    }


    void kekulizeFused(RWMol &mol,
                       const VECT_INT_VECT &arings, unsigned int maxBackTracks) { 
      
      // get all teh atoms in the ring system
      INT_VECT allAtms;
      Union(arings, allAtms);

      // get all the atoms that are candidates to recieve a double bond
      // also mark atoms in teh fused system that are not aromatic to begin with
      // as done. Mark all the bonds that are part of the aromatic system
      // to be single bonds
      INT_VECT done;
      int nats = mol.getNumAtoms();
      int nbnds = mol.getNumBonds();
      boost::dynamic_bitset<> dBndCands(nats);
      boost::dynamic_bitset<> dBndAdds(nbnds);

      markDbondCands(mol, allAtms, dBndCands, done);

      //std::cerr << "candidates: ";
      //for(unsigned int i=0;i<nats;++i) std::cerr << dBndCands[i];
      //std::cerr << std::endl;

      INT_DEQUE astack;
      INT_INT_DEQ_MAP options;
      int lastOpt=-1;

      // ok the algorithm goes something like this
      // - start with an atom that has been marked aromatic before
      // - check if it can have a double bond
      // - add its neighbors to the stack
      // - check if one of its neighbors can also have a double bond
      // - if yes add a double bond.
      // - if multiple neighbors can have double bonds - add them to a
      //   options stack we may have to retrace out path if we chose the
      //   wrong neighbor to add the double bond
      // - if double bond added update the candidates for double bond
      // - move to the next atom on the stack and repeat the process
      // - if an atom that can have multiple a double bond has no
      //   neighbors that can take double bond - we made a mistake
      //   earlier by picking a wrong candidate for double bond
      // - in this case back track to where we made the mistake 
    
      int curr;
      INT_VECT_CI ai;
      INT_DEQUE btmoves;
      unsigned int numBT = 0; // number of back tracks so far
      while ( (done.size() < allAtms.size()) || (astack.size() > 0) ) {
        // pick a curr atom to work with
        if (astack.size() > 0) {
          curr = astack.front();
          astack.pop_front();
        }
        else {
          for (ai = allAtms.begin(); ai != allAtms.end(); ai++) {
            if (std::find(done.begin(), done.end(), (*ai)) == done.end()) {
              curr = (*ai);
              break;
            }
          }
        }
        done.push_back(curr);

        // loop over the neighbors if we can add double bonds or
        // simply push them onto the stack
        RWMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(curr));
        INT_DEQUE opts;
        bool cCand = false;
        if (dBndCands[curr]) {
          cCand = true;
        }
        int ncnd;
        // if we are here because of backtracking
        if (options.find(curr) != options.end()) {
          opts = options[curr];
          CHECK_INVARIANT(opts.size() > 0, "");
        }
        else {
          while (nbrIdx != endNbrs) {
            // ignore if the neighbor has already been dealt with before
            if (std::find(done.begin(), done.end(), (*nbrIdx)) != done.end()) {
              ++nbrIdx;
              continue;
            }
            // ignore if the neighbor is not part of the fused system
            if (std::find(allAtms.begin(),allAtms.end(),(*nbrIdx)) == allAtms.end()) {
              ++nbrIdx;
              continue;
            }
          
            // if the neighbor is not on the stack add it
            if (std::find(astack.begin(), astack.end(), (*nbrIdx)) == astack.end()) {
              astack.push_back(*nbrIdx);
            }
          
            // check if the neighbor is also a candidate for a double bond
            if (cCand && dBndCands[*nbrIdx] ){
              opts.push_back(*nbrIdx);
            } // end of curr atoms can have a double bond
            ++nbrIdx;
          } // end of looping over neighbors
        }
        // now add a double bond from current to one of the neighbors if we can
        if (cCand) {
          if (opts.size() > 0) {
            ncnd = opts.front();
            opts.pop_front();
            Bond *bnd = mol.getBondBetweenAtoms(curr, ncnd);
            bnd->setBondType(Bond::DOUBLE);
        
            // remove current and the neighbor from the dBndCands list
            dBndCands[curr]=0;
            dBndCands[ncnd]=0;

            // add them to the list of bonds to which have been made double 
            dBndAdds[bnd->getIdx()]=1;

            // if this is an atom we previously visted and picked we  
            // simply tried a different option now, overwrite the options
            // stored for this atoms 
            if (options.find(curr) != options.end() ) {
              if(opts.size() == 0){
                options.erase(curr);
                btmoves.pop_back();
                if (btmoves.size() > 0) {
                  lastOpt = btmoves.back();
                }
                else {
                  lastOpt = -1;
                }
              }
              else {
                options[curr] = opts;
              }
            }
            else {
              // this is new atoms we are trying and have other
              // neighbors as options to add double bond store this to
              // the options stack, we may have made a mistake in
              // which one we chose and have to return here
              if (opts.size() > 0) {
                lastOpt = curr;
                btmoves.push_back(lastOpt);
                options[curr] = opts;
              }
            }
                      
          } // end of adding a double bond
          else {
            // we have an atom that should be getting a double bond
            // but none of the negihbors can take one. Most likely
            // because of a wrong choice earlier so back track
            if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
              //std::cerr << "PRE BACKTRACK" << std::endl;
              //mol.debugMol(std::cerr);
              backTrack(mol, options, lastOpt, done, astack, dBndCands, dBndAdds);
              //std::cerr << "POST BACKTRACK" << std::endl;
              //mol.debugMol(std::cerr);
              numBT++;
            }
            else {
              // we exhausted all option (or crossed the allowed
              // number of backTracks) and we still need to backtrack
              // can't kekulize this thing
              std::ostringstream errout;
              errout << "Can't kekulize mol " << std::endl; 
              std::string msg = errout.str();
              BOOST_LOG(rdErrorLog) << msg<< std::endl;
              throw MolSanitizeException(msg);
            }
          } // end of else try to backtrack
        } // end of curr atom atom being a cand for double bond
      } // end of while we are not done with all atoms
    }

  }// end of utility namespace

  namespace MolOps {
    void Kekulize(RWMol &mol, bool markAtomsBonds, unsigned int maxBackTracks) {
      ROMol::AtomIterator ai;
      INT_VECT valences;
      // REVIEW
      int numAtoms=mol.getNumAtoms();
      valences.reserve(numAtoms);
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        valences.push_back((*ai)->getImplicitValence(true));
      }

      // before everything do implicit valence calculation and store them
      // we will repeat after kekulizaition and compare for the sake of error
      // checking
    
      // A bit on the state of the molecule at this point
      // - aromatic and non aromatic atoms and bonds may be mixed up

      //  FIX: what does this mean?
      // - no assumption of implicit hydrogen calculation is made.

      // - for all aromatic bonds it is assumed that that both the following
      //   are true:
      //       - getIsAromatic returns true
      //       - getBondType return aromatic
      // - all aromatic atoms return true for "getIsAromatic"

      // first find the all the simple rings in the molecule
      VECT_INT_VECT arings;
    
      int nrings = MolOps::symmetrizeSSSR(mol, arings);

      // convert the rings to bonds ids
      VECT_INT_VECT brings;
      RingUtils::convertToBonds(arings, brings, mol);

      // make a the neighbor map for the rings i.e. a ring is a
      // neighbor a another candidate ring if shares atleast one bond
      // useful to figure out fused systems
      INT_INT_VECT_MAP neighMap;
      RingUtils::makeRingNeighborMap(brings, neighMap);

      INT_VECT doneRs;
      int curr = 0;
      bool allDone=false;
      int rix;
      int cnrs = arings.size();
      boost::dynamic_bitset<> fusDone(cnrs);
      INT_VECT fused;
      while (curr < cnrs) {
        fused.resize(0);
        RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
        VECT_INT_VECT frings;
        INT_VECT_CI ci;
        // FIX: why are we copying into frings instead of just passing
        //  fused to kekulizeFused (where it's a const & anyway)?
        for (ci = fused.begin(); ci != fused.end(); ci++) {
          // REVIEW
          frings.push_back(arings[*ci]);
        }
        kekulizeFused(mol, frings, maxBackTracks);
        for (rix = 0; rix < cnrs; rix++) {
          if (!fusDone[rix]) {
            curr = rix;
            break;
          }
        }
        if (rix == cnrs) {
          break;
        }
      }
    

      RWMol::BondIterator bi;
      if (markAtomsBonds) {
        // if we want the atoms and bonds to be marked non-aromatic do
        // that here.
        for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
          (*ai)->setIsAromatic(false);
        }
        for (bi = mol.beginBonds(); bi != mol.endBonds(); bi++) {
          (*bi)->setIsAromatic(false);
        }
      }
    
      for (bi = mol.beginBonds(); bi != mol.endBonds(); bi++) {
        // by now the bondtype should have already changed from aromatic
        if (markAtomsBonds && (*bi)->getBondType() == Bond::AROMATIC) {
          std::ostringstream errout;
          errout << "Kekulization somehow did not convert bond " << (*bi)->getIdx();
          std::string msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg<< std::endl;
          throw MolSanitizeException(msg);
        }
      }

      // ok some error checking here force a implicit valence
      // calculation that should do some error checking by itself. In
      // addition compare them to what it was before kekulizing
      int i = 0;
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        int val = (*ai)->getImplicitValence(true);
        if (val != valences[i]) {
          std::ostringstream errout;
          errout << "Kekulization somehow screwed up valence on " << (*ai)->getIdx();
          std::string msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg<< std::endl;
          throw MolSanitizeException(msg);
        }
        i++;
      }
    }
  } // end of namespace MolOps
} // end of namespace RDKit
