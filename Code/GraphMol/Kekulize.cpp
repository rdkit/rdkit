// $Id$
//
//  Copyright (C) 2001-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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

    // Determine whether or not a molecule is to the left of Carbon
    bool isEarlyAtom(int atomicNum){
      // FIX: this is duplicated from Atom.cpp. It should be defined once!
      return (4 - PeriodicTable::getTable()->getNouterElecs(atomicNum)) > 0;
    }

    void backTrack(RWMol &mol, 
                   INT_INT_DEQ_MAP &options,
                   int lastOpt, 
                   INT_VECT &done, 
                   INT_DEQUE &aqueue, 
                   boost::dynamic_bitset<> &dBndCands, 
                   boost::dynamic_bitset<> &dBndAdds) {
      // so we made a wrong turn at the lastOpt
      //remove on done list that comes after the lastOpt including itself

      INT_VECT_I ei = std::find(done.begin(), done.end(), lastOpt);
      INT_VECT tdone;
      tdone.insert(tdone.end(),done.begin(),ei);

      INT_VECT_CRI eri = std::find(done.rbegin(), done.rend(), lastOpt);
      ++eri;
      // and push them back onto the stack
      for (INT_VECT_CRI ri = done.rbegin(); ri != eri; ++ri) {
        aqueue.push_front(*ri);
      }
    
      // remove any double bonds that were add since we passed through lastOpt
      Bond *bnd;
      unsigned int nbnds = mol.getNumBonds();
      for (unsigned int bi = 0; bi < nbnds; bi++) {
        if (dBndAdds[bi]) {
          bnd = mol.getBondWithIdx(bi);
          int aid1 = bnd->getBeginAtomIdx();
          int aid2 = bnd->getEndAtomIdx();
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
                        boost::dynamic_bitset<> &dBndCands,
                        INT_VECT &questions,
                        INT_VECT &done){
      // ok this function does more than mark atoms that are candidates for 
      // double bonds during kekulization
      // - check that an non aromatic atom does not have any aromatic bonds
      // - marks all aromatic bonds to single bonds
      // - marks atoms that can take a double bond 
      
      bool hasAromaticOrDummyAtom=false;
      for (INT_VECT_CI adx = allAtms.begin(); adx != allAtms.end();
           ++adx) {
        if(mol.getAtomWithIdx(*adx)->getIsAromatic() ||
           !mol.getAtomWithIdx(*adx)->getAtomicNum()){
          hasAromaticOrDummyAtom=true;
          break;
        }
      }
      // if there's not at least one atom in the ring that's
      // marked as being aromatic or a dummy,
      // there's no point in continuing:
      if(!hasAromaticOrDummyAtom) return;

      std::vector<Bond *> makeSingle;

      for (INT_VECT_CI adx = allAtms.begin(); adx != allAtms.end();
           ++adx) {
        Atom *at = mol.getAtomWithIdx(*adx);

        if (!at->getIsAromatic() && at->getAtomicNum() ) {
          done.push_back(*adx);
          // make sure all the bonds on this atom are also non aromatic
          // i.e. can't have aromatic bond onto a non-aromatic atom
          RWMol::OEDGE_ITER beg,end;
          boost::tie(beg,end) = mol.getAtomBonds(at);
          while (beg != end) {
            // ok we can't have an aromatic atom
            if (mol[*beg]->getIsAromatic()) {
              std::ostringstream errout;
              errout << "Aromatic bonds on non aromatic atom " << at->getIdx();
              std::string msg = errout.str();
              BOOST_LOG(rdErrorLog) << msg << std::endl;
              throw MolSanitizeException(msg);
            }
            ++beg;
          }
          continue;
        }

        // count the number of neighbors connected with single,
        // double, or aromatic bonds. Along the way, mark
        // bonds that we will later mark as being single:
        int sbo = 0;
        RWMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(at);
        while (beg != end) {
          Bond *bond=mol[*beg].get();
          if (bond->getIsAromatic() &&
              (bond->getBondType()==Bond::SINGLE ||
               bond->getBondType()==Bond::DOUBLE ||
               bond->getBondType()==Bond::AROMATIC) ) {
            ++sbo;
            // mark this bond to be marked single later 
            // we don't want to do right now because it can screw-up the 
            // valence calculation to determine the number of hydrogens below
            makeSingle.push_back(bond);
          }
          else {
            sbo += (int)bond->getValenceContrib(at);
          }
          ++beg;
        }


        if(!at->getAtomicNum()){
          // dummies always start as candidates to have a double bond:
          dBndCands[*adx]=1;
          // but they don't have to have one, so mark them as questionable:
          questions.push_back(*adx);
        } else {
          // for non dummies, it's a bit more work to figure out if they
          // can take a double bond:

          sbo +=  at->getTotalNumHs();
          int dv = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
          int chrg = at->getFormalCharge();
          if(isEarlyAtom(at->getAtomicNum())) chrg*=-1; // fix for GitHub #65
          dv += chrg;
          int tbo = at->getTotalValence();
          int nRadicals=at->getNumRadicalElectrons();
          int totalDegree=at->getDegree()+at->getImplicitValence();
        
          const INT_VECT &valList =
            PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
          unsigned int vi = 1;
      
          while ( tbo>dv && vi<valList.size() && valList[vi]>0 ) {
            dv = valList[vi] + chrg;
            ++vi;
          }

          //std::cerr<<"  kek: "<<at->getIdx()<<" tbo:"<<tbo<<" sbo:"<<sbo<<" dv:"<<dv<<" totalDegree:"<<totalDegree<<" nRadicals:"<<nRadicals<<std::endl;
          if(totalDegree+nRadicals>=dv){
            // if our degree + nRadicals exceeds the default valence, 
            // there's no way we can take a double bond, just continue.
            continue;
          }

          // we're a candidate if our total current bond order + nRadicals + 1
          // matches the valence state
          // (including nRadicals here was SF.net issue 3349243)
          if( dv==(sbo+1+nRadicals) ){
            dBndCands[*adx]=1;
          } else if( !nRadicals && at->getNoImplicit() && dv==(sbo+2) ){
            // special case: there is currently no radical on the atom, but if
            // if we allow one then this is a candidate:
            dBndCands[*adx]=1;
          }
        }
      }// loop over all atoms in the fused system

      // now turn all the aromatic bond in this fused system to single
      for (std::vector<Bond *>::iterator bi=makeSingle.begin();
           bi!=makeSingle.end();++bi){
        (*bi)->setBondType(Bond::SINGLE);
      }
    }

    bool kekulizeWorker(RWMol &mol,const INT_VECT &allAtms,
                        boost::dynamic_bitset<> dBndCands,
                        boost::dynamic_bitset<> dBndAdds,
                        INT_VECT done,
                        unsigned int maxBackTracks){
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
      INT_DEQUE btmoves;
      unsigned int numBT = 0; // number of back tracks so far
      while ( (done.size() < allAtms.size()) || (astack.size() > 0) ) {
        // pick a curr atom to work with
        if (astack.size() > 0) {
          curr = astack.front();
          astack.pop_front();
        }
        else {
          for (INT_VECT_CI ai = allAtms.begin();
               ai != allAtms.end(); ++ai){
            if (std::find(done.begin(), done.end(),
                          (*ai)) == done.end()) {
              curr = (*ai);
              break;
            }
          }
        }
        done.push_back(curr);

        // loop over the neighbors if we can add double bonds or
        // simply push them onto the stack
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
          RWMol::ADJ_ITER nbrIdx,endNbrs;
          boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(mol.getAtomWithIdx(curr));
          while (nbrIdx != endNbrs) {
            // ignore if the neighbor has already been dealt with before
            if (std::find(done.begin(), done.end(), 
                          static_cast<int>(*nbrIdx)) != done.end()) {
              ++nbrIdx;
              continue;
            }
            // ignore if the neighbor is not part of the fused system
            if (std::find(allAtms.begin(),allAtms.end(),
                          static_cast<int>(*nbrIdx)) == allAtms.end()) {
              ++nbrIdx;
              continue;
            }
          
            // if the neighbor is not on the stack add it
            if (std::find(astack.begin(), astack.end(),
                          static_cast<int>(*nbrIdx)) == astack.end()) {
              astack.push_back(*nbrIdx);
            }
          
            // check if the neighbor is also a candidate for a double bond
            // the refinement that we'll make to the candidate check we've already
            // done is to make sure that the bond is either flagged as aromatic
            // or involves a dummy atom. This was Issue 3525076.
            // This fix is not really 100% of the way there: a situation like
            // that for Issue 3525076 but involving a dummy atom in the cage
            // could lead to the same failure. The full fix would require
            // a fairly detailed analysis of all bonds in the molecule to determine
            // which of them is eligible to be converted.
            if (cCand && dBndCands[*nbrIdx] &&
                (mol.getBondBetweenAtoms(curr,*nbrIdx)->getIsAromatic() ||
                 mol.getAtomWithIdx(curr)->getAtomicNum()==0 ||
                 mol.getAtomWithIdx(*nbrIdx)->getAtomicNum()==0)
                 ){
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
            // but none of the neighbors can take one. Most likely
            // because of a wrong choice earlier so back track
            if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
              //std::cerr << "PRE BACKTRACK" << std::endl;
              //mol.debugMol(std::cerr);
              backTrack(mol, options, lastOpt, done, astack,
                        dBndCands, dBndAdds);
              //std::cerr << "POST BACKTRACK" << std::endl;
              //mol.debugMol(std::cerr);
              numBT++;
            }
            else {
              return false;
            }
          } // end of else try to backtrack
        } // end of curr atom atom being a cand for double bond
      } // end of while we are not done with all atoms
      return true;
    }

    class QuestionEnumerator{
    public:
      QuestionEnumerator(const INT_VECT &questions) : d_questions(questions) , d_pos(1) {};
      INT_VECT next() {
        INT_VECT res;
        if(d_pos >= (0x1u<<d_questions.size())){
          return res;
        }
        for(unsigned int i=0;i<d_questions.size();++i){
          if(d_pos&(0x1u<<i)){
            res.push_back(d_questions[i]);
          }
        }
        ++d_pos;
        return res;
      };
    private:
      INT_VECT d_questions;
      unsigned int d_pos;
    };

    bool permuteDummiesAndKekulize(RWMol &mol,const INT_VECT &allAtms,
                                   boost::dynamic_bitset<> dBndCands,
                                   INT_VECT &questions,
                                   unsigned int maxBackTracks){
      boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
      for(INT_VECT_CI ai=allAtms.begin();ai!=allAtms.end();++ai){
        atomsInPlay[*ai]=1;
      }
      bool kekulized=false;
      QuestionEnumerator qEnum(questions);
      while(!kekulized && questions.size()){
        boost::dynamic_bitset<> dBndAdds(mol.getNumBonds());
        INT_VECT done;
#if 1
        // reset the state: all aromatic bonds are remarked to single:
        for(RWMol::BondIterator bi=mol.beginBonds();bi!=mol.endBonds();++bi){
          if((*bi)->getIsAromatic() && (*bi)->getBondType()!=Bond::SINGLE &&
             atomsInPlay[(*bi)->getBeginAtomIdx()] &&
             atomsInPlay[(*bi)->getEndAtomIdx()] ){
            (*bi)->setBondType(Bond::SINGLE);
          }
        }
#endif
        // pick a new permutation of the questionable atoms:
        const INT_VECT &switchOff=qEnum.next();
        if(!switchOff.size()) break;
        boost::dynamic_bitset<> tCands=dBndCands;
        for(INT_VECT_CI it=switchOff.begin();it!=switchOff.end();++it){
          tCands[*it]=0;
        }
#if 0
        std::cerr<<"permute: ";
        for (boost::dynamic_bitset<>::size_type i = 0; i < tCands.size(); ++i){
          std::cerr << tCands[i];
        }
        std::cerr<<std::endl;
#endif
        // try kekulizing again:
        kekulized=kekulizeWorker(mol,allAtms,tCands,dBndAdds,done,maxBackTracks);
      }
      return kekulized;
    }

    void kekulizeFused(RWMol &mol,
                       const VECT_INT_VECT &arings,
                       unsigned int maxBackTracks) { 
      // get all the atoms in the ring system
      INT_VECT allAtms;
      Union(arings, allAtms);

      // get all the atoms that are candidates to receive a double bond
      // also mark atoms in the fused system that are not aromatic to begin with
      // as done. Mark all the bonds that are part of the aromatic system
      // to be single bonds
      INT_VECT done;
      INT_VECT questions;
      unsigned int nats = mol.getNumAtoms();
      unsigned int nbnds = mol.getNumBonds();
      boost::dynamic_bitset<> dBndCands(nats);
      boost::dynamic_bitset<> dBndAdds(nbnds);

      markDbondCands(mol, allAtms, dBndCands, questions, done);
#if 0
      std::cerr << "candidates: ";
      for(int i=0;i<nats;++i) std::cerr << dBndCands[i];
      std::cerr << std::endl;
#endif

      bool kekulized;
      kekulized=kekulizeWorker(mol,allAtms,dBndCands,dBndAdds,done,maxBackTracks);
      if(!kekulized && questions.size()){
        // we failed, but there are some dummy atoms we can try permuting.
        kekulized=permuteDummiesAndKekulize(mol,allAtms,dBndCands,questions,maxBackTracks);
      }
      if(!kekulized){
        // we exhausted all option (or crossed the allowed
        // number of backTracks) and we still need to backtrack
        // can't kekulize this thing
        std::ostringstream errout;
        errout << "Can't kekulize mol " << std::endl; 
        std::string msg = errout.str();
        BOOST_LOG(rdErrorLog) << msg<< std::endl;
        throw MolSanitizeException(msg);
      }
    }
  }// end of utility namespace

  namespace MolOps {
    void Kekulize(RWMol &mol, bool markAtomsBonds,
                  unsigned int maxBackTracks) {

      // before everything do implicit valence calculation and store them
      // we will repeat after kekulization and compare for the sake of error
      // checking
      INT_VECT valences;
      // REVIEW
      int numAtoms=mol.getNumAtoms();
      valences.reserve(numAtoms);
      for (ROMol::AtomIterator ai = mol.beginAtoms();
           ai != mol.endAtoms(); ++ai) {
        (*ai)->calcImplicitValence(false);
        valences.push_back((*ai)->getTotalValence());
      }
    
      // A bit on the state of the molecule at this point
      // - aromatic and non aromatic atoms and bonds may be mixed up

      // - for all aromatic bonds it is assumed that that both the following
      //   are true:
      //       - getIsAromatic returns true
      //       - getBondType return aromatic
      // - all aromatic atoms return true for "getIsAromatic"

      // first find the all the simple rings in the molecule
      VECT_INT_VECT arings;
      if(mol.getRingInfo()->isInitialized()){
        arings = mol.getRingInfo()->atomRings();
      } else {
        MolOps::findSSSR(mol, arings);
      }

      VECT_INT_VECT brings;
      brings = mol.getRingInfo()->bondRings();
      //RingUtils::convertToBonds(arings, brings, mol);

      // make a the neighbor map for the rings i.e. a ring is a
      // neighbor to another candidate ring if it shares at least
      // one bond
      // useful to figure out fused systems
      INT_INT_VECT_MAP neighMap;
      RingUtils::makeRingNeighborMap(brings, neighMap);

      int curr = 0;
      int cnrs = arings.size();
      boost::dynamic_bitset<> fusDone(cnrs);
      while (curr < cnrs) {
        INT_VECT fused;
        RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
        VECT_INT_VECT frings;
        for (INT_VECT_CI ci = fused.begin();
             ci != fused.end();++ci) {
          frings.push_back(arings[*ci]);
        }
        kekulizeFused(mol, frings, maxBackTracks);
        int rix;
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

      if (markAtomsBonds) {
        // if we want the atoms and bonds to be marked non-aromatic do
        // that here.
        for (ROMol::BondIterator bi = mol.beginBonds();
             bi != mol.endBonds(); ++bi) {
          (*bi)->setIsAromatic(false);
        }
        for (ROMol::AtomIterator ai = mol.beginAtoms();
             ai != mol.endAtoms(); ++ai) {
          if((*ai)->getIsAromatic()){
            if(!mol.getRingInfo()->numAtomRings((*ai)->getIdx())){
              std::ostringstream errout;
              errout << "non-ring atom " << (*ai)->getIdx()<<" marked aromatic";
              std::string msg = errout.str();
              BOOST_LOG(rdErrorLog) << msg<< std::endl;
              throw MolSanitizeException(msg);
            }
            (*ai)->setIsAromatic(false);
            // make sure "explicit" Hs on things like pyrroles don't hang around
            // this was Github Issue 141
            if(( (*ai)->getAtomicNum()==7 || (*ai)->getAtomicNum()==15 ) &&
               (*ai)->getFormalCharge()==0 &&
               (*ai)->getNumExplicitHs()==1){
              (*ai)->setNoImplicit(false);
              (*ai)->setNumExplicitHs(0);
              (*ai)->updatePropertyCache(false);
            }
          }
        }
      }

      // ok some error checking here force a implicit valence
      // calculation that should do some error checking by itself. In
      // addition compare them to what they were before kekulizing
      int i = 0;
      for (ROMol::AtomIterator ai = mol.beginAtoms();
           ai != mol.endAtoms(); ++ai) {
        int val = (*ai)->getTotalValence();
        if (val != valences[i]) {
          std::ostringstream errout;
          errout << "Kekulization somehow screwed up valence on " << (*ai)->getIdx() <<": "<<val<<"!="<<valences[i]<<std::endl;
          std::string msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg<< std::endl;
          throw MolSanitizeException(msg);
        }
        i++;
      }
    }
  } // end of namespace MolOps
} // end of namespace RDKit
