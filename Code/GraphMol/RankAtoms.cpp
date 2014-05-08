// $Id$
//
//  Copyright (C) 2002-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RankAtoms.h>
#include <boost/cstdint.hpp>
#include <boost/foreach.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/hash/hash.hpp>

#include <boost/lambda/lambda.hpp>
#include <list>
#include <algorithm>
using namespace boost::lambda;

//#define VERBOSE_CANON 1
//#define VERYVERBOSE_CANON 1

namespace RankAtoms{
  using namespace RDKit;
  // --------------------------------------------------
  //
  // grabs the corresponding primes for the rank vector ranks
  //
  // --------------------------------------------------
  void getPrimes(const INT_VECT &ranks,INT_VECT &res){
    PRECONDITION(res.size()==0,"");
    res.reserve(ranks.size());
    for(INT_VECT_CI ivCIt=ranks.begin();ivCIt!=ranks.end();++ivCIt){
      res.push_back(firstThousandPrimes[(*ivCIt)%NUM_PRIMES_AVAIL]);
    }
  }

  // --------------------------------------------------
  //
  // blows out any indices in indicesInPlay which correspond to unique ranks
  //
  // --------------------------------------------------
  void updateInPlayIndices(const INT_VECT &ranks,INT_LIST &indicesInPlay){
    INT_LIST::iterator ivIt=indicesInPlay.begin();    
    while(ivIt!=indicesInPlay.end()){
      // find the first instance of this rank:
      INT_VECT::const_iterator pos=std::find(ranks.begin(),ranks.end(),ranks[*ivIt]);
      ++pos;
      // now check to see if there is at least one more:
      if( std::find(pos,ranks.end(),ranks[*ivIt])==ranks.end()){
        INT_LIST::iterator tmpIt = ivIt;
        ++ivIt;
        indicesInPlay.erase(tmpIt);
      } else {
        ++ivIt;
      }
    }
  }
  
  // --------------------------------------------------
  //
  // for each index in indicesInPlay, generate the products of the adjacent
  // elements
  //
  //  The products are weighted by the order of the bond connecting the atoms.
  //
  // --------------------------------------------------
  void calcAdjacentProducts(unsigned int nAtoms,
                            const INT_VECT &valVect,
                            double const *adjMat,
                            const INT_LIST &indicesInPlay,
                            DOUBLE_VECT &res,
                            bool useSelf=true,
                            double tol=1e-6){
    PRECONDITION(valVect.size() >= nAtoms,"");
    PRECONDITION(res.size() == 0,"");
    PRECONDITION(adjMat,"");
    for(INT_LIST::const_iterator idxIt=indicesInPlay.begin();
        idxIt != indicesInPlay.end();
        ++idxIt){
      double accum;
      if(useSelf)
        accum=valVect[*idxIt];
      else
        accum=1.0;
      const unsigned int iTab = (*idxIt)*nAtoms;
      for(unsigned int j=0;j<nAtoms;++j){
        double elem=adjMat[iTab+j];
        if(elem>tol){
          if(elem<2.-tol){
            accum *= valVect[j];
          } else {
            accum *= pow(static_cast<double>(valVect[j]),
                         static_cast<int>(elem));
          }
        }
      }
      res.push_back(accum);
    }
  }

  template <typename T>
  void debugVect(const std::vector<T> arg){
    typename std::vector<T>::const_iterator viIt;
    for(viIt=arg.begin();viIt!=arg.end();++viIt){
      BOOST_LOG(rdDebugLog)<< *viIt << " ";
    }
    BOOST_LOG(rdDebugLog)<< std::endl;
  }
  
  // --------------------------------------------------
  //
  //  This is one round of the process from Step III in the Daylight
  //  paper
  //
  // --------------------------------------------------
  unsigned int iterateRanks(unsigned int nAtoms,INT_VECT &primeVect,
                            DOUBLE_VECT &atomicVect,
                            INT_LIST &indicesInPlay,
                            double *adjMat,
                            INT_VECT &ranks,
                            VECT_INT_VECT *rankHistory,unsigned int stagnantTol){
    PRECONDITION(!rankHistory||rankHistory->size()>=nAtoms,"bad rankHistory size");
    bool done = false;
    unsigned int numClasses = countClasses(ranks);
    unsigned int lastNumClasses = 0;
    unsigned int nCycles = 0;
    unsigned int nStagnant=0;

    //
    // loop until either we finish or no improvement is seen
    //
#ifdef VERBOSE_CANON
    for(unsigned int i=0;i<nAtoms;i++) {
      BOOST_LOG(rdDebugLog)<< "\t\t>:" << i << " " << ranks[i] << std::endl;
    }
    BOOST_LOG(rdDebugLog)<< "\t\t-*-*-*-*-" << std::endl;
#endif
    while(!done && nCycles < nAtoms){
      // determine which atomic indices are in play (which have duplicate ranks)
      if(rankHistory){
        for(INT_LIST_CI idx=indicesInPlay.begin();idx!=indicesInPlay.end();++idx){
          (*rankHistory)[*idx].push_back(ranks[*idx]);
        }
      }
      updateInPlayIndices(ranks,indicesInPlay);
      if(indicesInPlay.empty()) break;
        
#ifdef VERYVERBOSE_CANON
      BOOST_LOG(rdDebugLog)<< "IN PLAY:" << std::endl;
      BOOST_LOG(rdDebugLog)<< "\t\t->";
      for(INT_LIST::const_iterator tmpI=indicesInPlay.begin();tmpI != indicesInPlay.end();tmpI++){
        BOOST_LOG(rdDebugLog)<< " " << *tmpI;
      }
      BOOST_LOG(rdDebugLog)<< std::endl;
      BOOST_LOG(rdDebugLog)<< "\t\t---------" << std::endl;
#endif
      //-------------------------
      // Step (2):
      //    Get the products of adjacent primes
      //-------------------------
      primeVect.resize(0);
      getPrimes(ranks,primeVect);
      atomicVect.resize(0);
      calcAdjacentProducts(nAtoms,primeVect,adjMat,indicesInPlay,atomicVect,false);
#ifdef VERYVERBOSE_CANON
      BOOST_LOG(rdDebugLog)<< "primes: ";
      debugVect(primeVect);
      BOOST_LOG(rdDebugLog)<< "products: ";
      debugVect(atomicVect);
#endif


      //-------------------------
      // Steps (3) and (4)
      //   sort the products and count classes
      //-------------------------
      sortAndRankVect(nAtoms,atomicVect,indicesInPlay,ranks);
      lastNumClasses = numClasses;
      numClasses = countClasses(ranks);
      if(numClasses == lastNumClasses) nStagnant++;
#ifdef VERYVERBOSE_CANON
      int tmpOff=0;
      for(unsigned int i=0;i<nAtoms;i++){
        //for(INT_LIST::const_iterator tmpI=indicesInPlay.begin();tmpI != indicesInPlay.end();tmpI++){
        BOOST_LOG(rdDebugLog)<< "\t\ti:" << i << "\t" << ranks[i] << "\t" << primeVect[i];
        if(std::find(indicesInPlay.begin(),indicesInPlay.end(),i)!=indicesInPlay.end()){
          BOOST_LOG(rdDebugLog)<< "\t" << atomicVect[tmpOff];
          tmpOff++;
        }
        BOOST_LOG(rdDebugLog)<< std::endl;    }
      BOOST_LOG(rdDebugLog)<< "\t\t---------" << std::endl;
#endif
        
      // terminal condition, we'll allow a single round of stagnancy
      if(numClasses == nAtoms || nStagnant > stagnantTol) done = 1;
      nCycles++;
    }
#ifdef VERBOSE_CANON
    BOOST_LOG(rdDebugLog)<< ">>>>>> done inner iteration. static: "<< nStagnant << " ";
    BOOST_LOG(rdDebugLog)<< nCycles << " " << nAtoms << " " << numClasses << std::endl;
#ifdef VERYVERBOSE_CANON
    for(unsigned int i=0;i<nAtoms;i++) {
      BOOST_LOG(rdDebugLog)<< "\t" << i << " " << ranks[i] << std::endl;
    }
#endif    
    if(nCycles == nAtoms){
      BOOST_LOG(rdWarningLog) << "WARNING: ranking bottomed out" << std::endl;
    }
#endif
    return numClasses;
  }
  unsigned int iterateRanks2(unsigned int nAtoms,INT_VECT &primeVect,
                            DOUBLE_VECT &atomicVect,
                            INT_LIST &indicesInPlay,
                            double *adjMat,
                             INT_VECT &ranks,VECT_DOUBLE_VECT &nRanks,
                            VECT_INT_VECT *rankHistory,unsigned int stagnantTol){
    PRECONDITION(!rankHistory||rankHistory->size()>=nAtoms,"bad rankHistory size");
    bool done = false;
    unsigned int numClasses = countClasses(ranks);
    unsigned int lastNumClasses = 0;
    unsigned int nCycles = 0;
    unsigned int nStagnant=0;

    //
    // loop until either we finish or no improvement is seen
    //
#ifdef VERBOSE_CANON
    for(unsigned int i=0;i<nAtoms;i++) {
      BOOST_LOG(rdDebugLog)<< "\t\t>:" << i << " " << ranks[i] << std::endl;
    }
    BOOST_LOG(rdDebugLog)<< "\t\t-*-*-*-*-" << std::endl;
#endif
    while(!done && nCycles < nAtoms){
      // determine which atomic indices are in play (which have duplicate ranks)
      if(rankHistory){
        BOOST_FOREACH(int idx,indicesInPlay){
          (*rankHistory)[idx].push_back(ranks[idx]);
        }
      }
      updateInPlayIndices(ranks,indicesInPlay);
      if(indicesInPlay.empty()) break;
        
#ifdef VERYVERBOSE_CANON
      BOOST_LOG(rdDebugLog)<< "IN PLAY:" << std::endl;
      BOOST_LOG(rdDebugLog)<< "\t\t->";
      for(INT_LIST::const_iterator tmpI=indicesInPlay.begin();tmpI != indicesInPlay.end();tmpI++){
        BOOST_LOG(rdDebugLog)<< " " << *tmpI;
      }
      BOOST_LOG(rdDebugLog)<< std::endl;
      BOOST_LOG(rdDebugLog)<< "\t\t---------" << std::endl;
#endif
      //-------------------------
      // Step (2):
      //    Get the products of adjacent primes
      //-------------------------
      primeVect.resize(0);
      getPrimes(ranks,primeVect);
      atomicVect.resize(0);
      calcAdjacentProducts(nAtoms,primeVect,adjMat,indicesInPlay,atomicVect,false);
#ifdef VERYVERBOSE_CANON
      BOOST_LOG(rdDebugLog)<< "primes: ";
      debugVect(primeVect);
      BOOST_LOG(rdDebugLog)<< "products: ";
      debugVect(atomicVect);
#endif

      unsigned int p=0;
      BOOST_FOREACH(int idx,indicesInPlay){
        nRanks[idx].push_back(atomicVect[p++]);
      }
#ifdef VERYVERBOSE_CANON
      for(int idx=0;idx<nAtoms;++idx){
        std::cerr<<"  nranks["<<idx<<"]: ";
        std::copy(nRanks[idx].begin(),nRanks[idx].end(),std::ostream_iterator<double>(std::cerr," "));
        std::cerr<<"\n";
      }
#endif

      
      //-------------------------
      // Steps (3) and (4)
      //   sort the products and count classes
      //-------------------------
      rankVect(nRanks,ranks);
      //sortAndRankVect2(nRanks,indicesInPlay,ranks);
      lastNumClasses = numClasses;
      numClasses = countClasses(ranks);
      if(numClasses == lastNumClasses) nStagnant++;
#ifdef VERYVERBOSE_CANON
      int tmpOff=0;
      for(unsigned int i=0;i<nAtoms;i++){
        //for(INT_LIST::const_iterator tmpI=indicesInPlay.begin();tmpI != indicesInPlay.end();tmpI++){
        BOOST_LOG(rdDebugLog)<< "\t\ti:" << i << "\t" << ranks[i] << "\t" << primeVect[i];
        if(std::find(indicesInPlay.begin(),indicesInPlay.end(),i)!=indicesInPlay.end()){
          BOOST_LOG(rdDebugLog)<< "\t" << atomicVect[tmpOff];
          tmpOff++;
        }
        BOOST_LOG(rdDebugLog)<< std::endl;    }
      BOOST_LOG(rdDebugLog)<< "\t\t---------" << std::endl;
#endif
        
      // terminal condition, we'll allow a single round of stagnancy
      if(numClasses == nAtoms || nStagnant > stagnantTol) done = 1;
      nCycles++;
    }
#ifdef VERBOSE_CANON
    BOOST_LOG(rdDebugLog)<< ">>>>>> done inner iteration. static: "<< nStagnant << " ";
    BOOST_LOG(rdDebugLog)<< nCycles << " " << nAtoms << " " << numClasses << std::endl;
#ifdef VERYVERBOSE_CANON
    for(unsigned int i=0;i<nAtoms;i++) {
      BOOST_LOG(rdDebugLog)<< "\t" << i << " " << ranks[i] << std::endl;
    }
#endif    
    if(nCycles == nAtoms){
      BOOST_LOG(rdWarningLog) << "WARNING: ranking bottomed out" << std::endl;
    }
#endif
    return numClasses;
  }

  // --------------------------------------------------
  //
  // Calculates invariants for the atoms of a molecule
  //
  // NOTE: if the atom has not had chirality info pre-calculated, it doesn't
  // much matter what value includeChirality has!
  // --------------------------------------------------
  void buildAtomInvariants(const ROMol &mol,INVAR_VECT &res,
                           bool includeChirality,
                           bool includeIsotopes){
    PRECONDITION(res.size()>=mol.getNumAtoms(),"res vect too small");
    unsigned int atsSoFar=0;
    std::vector<boost::uint64_t> tres(mol.getNumAtoms());
    for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++){
      Atom const *atom = *atIt;
      int nHs = atom->getTotalNumHs() % 8;
      int chg = abs(atom->getFormalCharge()) % 8;
      int chgSign = atom->getFormalCharge() > 0;
      int num =    atom->getAtomicNum() % 128;
      int nConns = atom->getDegree() % 8;
      int deltaMass=0;
      if(includeIsotopes && atom->getIsotope()){
        deltaMass = static_cast<int>(atom->getIsotope() -
                                     PeriodicTable::getTable()->getMostCommonIsotope(atom->getAtomicNum()));
        deltaMass += 128;
        if(deltaMass < 0) deltaMass = 0;
        else deltaMass = deltaMass % 256;
      }

      
      // figure out the minimum-sized ring we're involved in
      int inRing = 0;
      if(atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())){
        RingInfo *ringInfo=atom->getOwningMol().getRingInfo();
        inRing=3;
        while(inRing<256){
          if(ringInfo->isAtomInRingOfSize(atom->getIdx(),inRing)){
            break;
          } else {
            inRing++;
          }
        }
      }
      inRing = inRing % 16;

      boost::uint64_t invariant = 0;
      invariant = (invariant << 3) | nConns;
      // we used to include the number of explicitHs, but that
      // didn't make much sense. TotalValence is another possible
      // discriminator here, but the information is essentially
      // redundant with nCons, num, and nHs.
      // invariant = (invariant << 4) | totalVal; 
      invariant = (invariant << 7) | num;  
      invariant = (invariant << 8) | deltaMass;
      invariant = (invariant << 3) | nHs;
      invariant = (invariant << 4) | inRing;
      invariant = (invariant << 3) | chg;
      invariant = (invariant << 1) | chgSign;
      if(includeChirality ){
        int isR=0;
        if( atom->hasProp("_CIPCode")){
          std::string cipCode;
          atom->getProp("_CIPCode",cipCode);
          if(cipCode=="R"){
            isR=1;
          } else {
            isR=2;
          }
        }
        invariant = (invariant << 2) | isR;
      }

      // now deal with cis/trans - this is meant to address issue 174
      // loop over the bonds on this atom and check if we have a double bond with
      // a chiral code marking 
      if (includeChirality) {
        ROMol::OBOND_ITER_PAIR atomBonds = atom->getOwningMol().getAtomBonds(atom);
        int isT=0;
        while (atomBonds.first != atomBonds.second){
          BOND_SPTR tBond = atom->getOwningMol()[*(atomBonds.first)];
          if( (tBond->getBondType() == Bond::DOUBLE) &&
              (tBond->getStereo()>Bond::STEREOANY )) {
            if (tBond->getStereo()==Bond::STEREOE) {
              isT = 1;
            } else if(tBond->getStereo()==Bond::STEREOZ) {
              isT=2;
            }
            break;
          }
          atomBonds.first++;
        }
        invariant = (invariant << 2) | isT;
      }
      tres[atsSoFar++] = invariant;
    }
    if(includeChirality){
      // ring stereochemistry
      boost::dynamic_bitset<> adjusted(mol.getNumAtoms());
      for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++){
        Atom const *atom = *atIt;
        tres[atom->getIdx()] = tres[atom->getIdx()]<<2;
      }
      for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++){
        Atom const *atom = *atIt;
        if((atom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CW ||
            atom->getChiralTag()==Atom::CHI_TETRAHEDRAL_CCW) &&
           atom->hasProp("_ringStereoAtoms")){
          //atom->hasProp("_CIPRank") &&
          //!atom->hasProp("_CIPCode")){
          ROMol::ADJ_ITER beg,end;
          boost::tie(beg,end) = mol.getAtomNeighbors(atom);
          unsigned int nCount=0;
          while(beg!=end){
            unsigned int nbrIdx=mol[*beg]->getIdx();
            if(!adjusted[nbrIdx]){
              tres[nbrIdx] |= nCount%4;
              adjusted.set(nbrIdx);
            }
            ++nCount;
            ++beg;
          }
        }
      }
    }
    for(unsigned int i=0;i<mol.getNumAtoms();++i) res[i]=tres[i];
  }

  void buildFragmentAtomInvariants(const ROMol &mol,INVAR_VECT &res,
                                   bool includeChirality,
                                   const boost::dynamic_bitset<> &atomsToUse,
                                   const boost::dynamic_bitset<> &bondsToUse,
                                   const std::vector<std::string> *atomSymbols
                                   ){
    PRECONDITION(res.size()>=mol.getNumAtoms(),"res vect too small");

    std::vector<int> degrees(mol.getNumAtoms(),0);
    for(unsigned int i=0;i<bondsToUse.size();++i){
      if(bondsToUse[i]){
        const Bond *bnd=mol.getBondWithIdx(i);
        degrees[bnd->getBeginAtomIdx()]++;
        degrees[bnd->getEndAtomIdx()]++;
      }
    }
    
    for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();++atIt){
      Atom const *atom = *atIt;
      int aIdx=atom->getIdx();
      if(!atomsToUse[aIdx]){
        res[aIdx] = 0;
        continue;
      }
      boost::uint64_t invariant = 0;

      int nConns = degrees[aIdx]% 8;
      invariant = (invariant << 3) | nConns;

      if(!atomSymbols){
        int chg = abs(atom->getFormalCharge()) % 8;
        int chgSign = atom->getFormalCharge() > 0;
        int num =    atom->getAtomicNum() % 128;
        int deltaMass=0;
        if(atom->getIsotope()){
          deltaMass = static_cast<int>(atom->getIsotope() -
                                       PeriodicTable::getTable()->getMostCommonIsotope(atom->getAtomicNum()));
          deltaMass += 128;
          if(deltaMass < 0) deltaMass = 0;
          else deltaMass = deltaMass % 256;
        }
        invariant = (invariant << 7) | num;  
        invariant = (invariant << 8) | deltaMass;
        invariant = (invariant << 3) | chg;
        invariant = (invariant << 1) | chgSign;
        invariant = (invariant << 1) | atom->getIsAromatic();
      } else {
        const std::string &symb=(*atomSymbols)[aIdx];
        boost::uint32_t hsh=gboost::hash_range(symb.begin(),symb.end());
        invariant = (invariant << 20) | (hsh%(1<<20));
      }

      // figure out the minimum-sized ring we're involved in
      int inRing = mol.getRingInfo()->minAtomRingSize(aIdx);
      inRing = inRing % 16;
      invariant = (invariant << 4) | inRing;

      if(includeChirality ){
        int isR=0;
        if( atom->hasProp("_CIPCode")){
          std::string cipCode;
          atom->getProp("_CIPCode",cipCode);
          if(cipCode=="R"){
            isR=1;
          } else {
            isR=2;
          }
        }
        invariant = (invariant << 2) | isR;
      }

      // now deal with cis/trans - this is meant to address issue 174
      // loop over the bonds on this atom and check if we have a double bond with
      // a chiral code marking 
      if (includeChirality) {
        ROMol::OBOND_ITER_PAIR atomBonds = mol.getAtomBonds(atom);
        int isT=0;
        while (atomBonds.first != atomBonds.second){
          BOND_SPTR tBond = mol[*(atomBonds.first)];
          atomBonds.first++;
          if(!bondsToUse[tBond->getIdx()]) continue;
          if( (tBond->getBondType() == Bond::DOUBLE) &&
              (tBond->getStereo()>Bond::STEREOANY )) {
            if (tBond->getStereo()==Bond::STEREOE) {
              isT = 1;
            } else if(tBond->getStereo()==Bond::STEREOZ) {
              isT=2;
            }
            break;
          }
        }
        invariant = (invariant << 2) | isT;
      }
      res[aIdx] = invariant;
    }
  }

  
}// end of RankAtoms namespace

namespace RDKit{
  namespace MolOps {
    // --------------------------------------------------
    //
    //  Daylight canonicalization, loosely based up on algorithm described in
    //     JCICS 29, 97-101, (1989)
    //  When appropriate, specific references are made to the algorithm
    //  description in that paper.  Steps refer to Table III of the paper
    //
    // --------------------------------------------------
    void rankAtoms(const ROMol &mol,INT_VECT &ranks,
                   bool breakTies,
                   bool includeChirality,
                   bool includeIsotopes,
                   VECT_INT_VECT *rankHistory){

      unsigned int i;
      unsigned int nAtoms = mol.getNumAtoms();
      PRECONDITION(ranks.size()>=nAtoms,"");
      PRECONDITION(!rankHistory||rankHistory->size()>=nAtoms,"bad rankHistory size");
      unsigned int stagnantTol=1;

      if(!mol.getRingInfo()->isInitialized()){
        MolOps::findSSSR(mol);
      }
    
      if(nAtoms > 1){
        double *adjMat = MolOps::getAdjacencyMatrix(mol, true);

        // ----------------------
        // generate atomic invariants, Step (1)
        // ----------------------
        INVAR_VECT invariants;
        invariants.resize(nAtoms);
        RankAtoms::buildAtomInvariants(mol,invariants,includeChirality,includeIsotopes);
      
#ifdef VERBOSE_CANON
        BOOST_LOG(rdDebugLog)<< "invariants:" << std::endl;
        for(i=0;i<nAtoms;i++){
          BOOST_LOG(rdDebugLog)<< i << " " << (long)invariants[i]<< std::endl;
        }

#endif

        DOUBLE_VECT atomicVect;
        atomicVect.resize(nAtoms);

        // ----------------------
        // iteration 1: Steps (3) and (4)
        // ----------------------

        // Unlike the original paper, we're going to keep track of the
        // ranks at each iteration and use those vectors to rank
        // atoms. This seems to lead to more stable evolution of the
        // ranks by avoiding ranks oscillating back and forth across
        // iterations. 
        VECT_DOUBLE_VECT nRanks(nAtoms);
        for(i=0;i<nAtoms;i++) nRanks[i].push_back(invariants[i]);

        // start by ranking the atoms using the invariants
        ranks.resize(nAtoms);
        RankAtoms::rankVect(nRanks,ranks);

        if(rankHistory){
          for(i=0;i<nAtoms;i++){
            (*rankHistory)[i].push_back(ranks[i]);
          }
        }

        // how many classes are present?
        unsigned int numClasses = RankAtoms::countClasses(ranks);

        if(numClasses != nAtoms){
          INT_VECT primeVect;
          primeVect.reserve(nAtoms);

          DOUBLE_VECT atomicVect;
          atomicVect.reserve(nAtoms);
      
          // indicesInPlay is used to track the atoms with non-unique ranks
          //  (we'll be modifying these in each step)
          INT_LIST indicesInPlay;
          for(i=0;i<nAtoms;i++) indicesInPlay.push_back(i);

          // if we aren't breaking ties here, allow the rank iteration to
          // go the full number of atoms:
          if(!breakTies) stagnantTol=nAtoms;

          bool done=indicesInPlay.empty();
          while(!done){
            //
            // do one round of iterations
            //
            numClasses = RankAtoms::iterateRanks2(nAtoms,primeVect,atomicVect,
                                                  indicesInPlay,adjMat,ranks,nRanks,
                                                  rankHistory,stagnantTol);

#ifdef VERBOSE_CANON
            BOOST_LOG(rdDebugLog)<< "************************ done outer iteration" << std::endl;
#endif  
#ifdef VERBOSE_CANON
            unsigned int tmpI;
            BOOST_LOG(rdDebugLog)<< "RANKS:" << std::endl;
            for(tmpI=0;tmpI<ranks.size();tmpI++){
              BOOST_LOG(rdDebugLog)<< "\t\t" << tmpI << " " << ranks[tmpI] << std::endl;
            }
            BOOST_LOG(rdDebugLog)<< std::endl;
#endif    
            //
            // This is the tiebreaker stage of things
            //
            if( breakTies && !indicesInPlay.empty() && numClasses<nAtoms){
              INT_VECT newRanks = ranks;

              // Add one to all ranks and multiply by two
              BOOST_FOREACH(int &nr,newRanks) {
                nr=(nr+1)*2;
              }
#ifdef VERBOSE_CANON
              BOOST_LOG(rdDebugLog)<< "postmult:" << std::endl;
              for(tmpI=0;tmpI<newRanks.size();tmpI++){
                BOOST_LOG(rdDebugLog)<< "\t\t" << newRanks[tmpI] << std::endl;
              }
              BOOST_LOG(rdDebugLog)<< std::endl;
#endif    

              //
              // find lowest duplicate rank with lowest invariant:
              //
              int lowestIdx=indicesInPlay.front();
              double lowestInvariant = invariants[lowestIdx];
              int lowestRank=newRanks[lowestIdx];
              BOOST_FOREACH(int ilidx,indicesInPlay){
                if(newRanks[ilidx]<=lowestRank){
                  if(newRanks[ilidx]<lowestRank ||
                     invariants[ilidx] <= lowestInvariant){
                    lowestRank = newRanks[ilidx];
                    lowestIdx = ilidx;
                    lowestInvariant = invariants[ilidx];
                  }
                }
              }

              //
              // subtract one from the lowest index, rerank and proceed
              //
              newRanks[lowestIdx] -= 1;
              RankAtoms::rankVect(newRanks,ranks);
              BOOST_FOREACH(int ilidx,indicesInPlay){
                nRanks[ilidx].push_back(ranks[ilidx]);
              }
#ifdef VERBOSE_CANON
              BOOST_LOG(rdDebugLog)<< "RE-RANKED ON:" << lowestIdx << std::endl;
              for(tmpI=0;tmpI<newRanks.size();tmpI++){
                BOOST_LOG(rdDebugLog)<< "\t\t" << newRanks[tmpI] << " " << ranks[tmpI] << std::endl;
              }
              BOOST_LOG(rdDebugLog)<< std::endl;
#endif    
            } else {
              done = true;
            }
          }
        }
      }
    } // end of function rankAtoms

    void rankAtomsInFragment(const ROMol &mol,INT_VECT &ranks,
                             const boost::dynamic_bitset<> &atomsToUse,
                             const boost::dynamic_bitset<> &bondsToUse,
                             const std::vector<std::string> *atomSymbols,
                             const std::vector<std::string> *bondSymbols,
                             bool breakTies,
                             VECT_INT_VECT *rankHistory){
      unsigned int nAtoms = mol.getNumAtoms();
      unsigned int nActiveAtoms = atomsToUse.count();
      PRECONDITION(ranks.size()>=nAtoms,"");
      PRECONDITION(!atomSymbols||atomSymbols->size()>=nAtoms,"bad atomSymbols");
      PRECONDITION(!rankHistory||rankHistory->size()>=nAtoms,"bad rankHistory size");
      PRECONDITION(mol.getRingInfo()->isInitialized(),"no ring information present");
      PRECONDITION(!rankHistory,"rankHistory not currently supported.");
      unsigned int stagnantTol=1;

      if(nActiveAtoms > 1){

        // ----------------------
        // generate atomic invariants, Step (1)
        // ----------------------
        INVAR_VECT invariants;
        invariants.resize(nAtoms);
        RankAtoms::buildFragmentAtomInvariants(mol,invariants,true,
                                               atomsToUse,bondsToUse,
                                               atomSymbols);
        INVAR_VECT tinvariants;
        tinvariants.resize(nActiveAtoms);
        unsigned int activeIdx=0;
        for(unsigned int aidx=0;aidx<nAtoms;++aidx){
          if(atomsToUse[aidx]){
            tinvariants[activeIdx++]=invariants[aidx];
          }
        }
      
#ifdef VERBOSE_CANON
        BOOST_LOG(rdDebugLog)<< "invariants:" << std::endl;
        for(unsigned int i=0;i<nActiveAtoms;i++){
          BOOST_LOG(rdDebugLog)<< i << " " << (long)tinvariants[i]<< std::endl;
        }

#endif
        // ----------------------
        // iteration 1: Steps (3) and (4)
        // ----------------------

        // start by ranking the atoms using the invariants
        VECT_DOUBLE_VECT nRanks(nActiveAtoms);
        for(unsigned int i=0;i<nActiveAtoms;i++) nRanks[i].push_back(tinvariants[i]);
        INT_VECT tranks(nActiveAtoms,0);
        RankAtoms::rankVect(nRanks,tranks);
#if 0
        if(rankHistory){
          for(unsigned int i=0;i<nAtoms;i++){
            (*rankHistory)[i].push_back(ranks[i]);
          }
        }
#endif
        // how many classes are present?
        unsigned int numClasses = RankAtoms::countClasses(tranks);
        if(numClasses != nActiveAtoms){
          double *tadjMat = new double[nActiveAtoms*nActiveAtoms];
          memset(static_cast<void *>(tadjMat),0,nActiveAtoms*nActiveAtoms*sizeof(double));
          if(!bondSymbols){
            double *adjMat = MolOps::getAdjacencyMatrix(mol,true,0,true,0,&bondsToUse);
            activeIdx=0;
            for(unsigned int aidx=0;aidx<nAtoms;++aidx){
              if(atomsToUse[aidx]){
                unsigned int activeIdx2=activeIdx+1;
                for(unsigned int aidx2=aidx+1;aidx2<nAtoms;++aidx2){
                  if(atomsToUse[aidx2]){
                    tadjMat[activeIdx*nActiveAtoms+activeIdx2]=adjMat[aidx*nAtoms+aidx2];
                    tadjMat[activeIdx2*nActiveAtoms+activeIdx]=adjMat[aidx2*nAtoms+aidx];
                    ++activeIdx2;
                  }
                }
                ++activeIdx;
              }
            }
          } else {
            // rank the bond symbols we have:
            std::vector<boost::uint32_t> tbranks(bondsToUse.size(),
                                                 0);
            for(unsigned int bidx=0;bidx<bondsToUse.size();++bidx){
              if(!bondsToUse[bidx]) continue;
              const std::string &symb=(*bondSymbols)[bidx];
              boost::uint32_t hsh=gboost::hash_range(symb.begin(),symb.end());
              tbranks[bidx]=hsh;
            }
            INT_VECT branks(bondsToUse.size(),1000000);
#ifdef VERBOSE_CANON
            std::cerr<<"  tbranks:";
            std::copy(tbranks.begin(),tbranks.end(),std::ostream_iterator<boost::uint32_t>(std::cerr," "));
            std::cerr<<std::endl;
#endif            
            RankAtoms::rankVect(tbranks,branks);
#ifdef VERBOSE_CANON
            std::cerr<<"  branks:";
            std::copy(branks.begin(),branks.end(),std::ostream_iterator<int>(std::cerr," "));
            std::cerr<<std::endl;
#endif            
            for(unsigned int bidx=0;bidx<bondsToUse.size();++bidx){
              if(!bondsToUse[bidx]) continue;
              const Bond *bond=mol.getBondWithIdx(bidx);
              unsigned int aidx1=bond->getBeginAtomIdx();
              unsigned int aidx2=bond->getEndAtomIdx();
              unsigned int tidx1=0;
              for(unsigned int iidx=0;iidx<aidx1;++iidx){
                if(atomsToUse[iidx]) ++tidx1;
              }
              unsigned int tidx2=0;
              for(unsigned int iidx=0;iidx<aidx2;++iidx){
                if(atomsToUse[iidx]) ++tidx2;
              }
              //const std::string &symb=(*bondSymbols)[bidx];
              //boost::uint32_t hsh=gboost::hash_range(symb.begin(),symb.end());
              //std::cerr<<" ::: "<<bidx<<"->"<<branks[bidx]<<std::endl;
              tadjMat[tidx1*nActiveAtoms+tidx2]=branks[bidx];
              tadjMat[tidx2*nActiveAtoms+tidx1]=branks[bidx];
            }
          }
          INT_VECT primeVect;
          primeVect.reserve(nActiveAtoms);

          DOUBLE_VECT atomicVect;
          atomicVect.reserve(nActiveAtoms);

#ifdef VERBOSE_CANON
          for(unsigned int aidx1=0;aidx1<nActiveAtoms;++aidx1){
            std::cerr<<aidx1<<" : ";
            for(unsigned int aidx2=aidx1+1;aidx2<nActiveAtoms;++aidx2){
              std::cerr<< tadjMat[aidx1*nActiveAtoms+aidx2]<<" ";
            }
            std::cerr<<std::endl;
          }
#endif
      
          // indicesInPlay is used to track the atoms with non-unique ranks
          //  (we'll be modifying these in each step)
          INT_LIST indicesInPlay;
          for(unsigned int i=0;i<nActiveAtoms;i++) indicesInPlay.push_back(i);

          // if we aren't breaking ties here, allow the rank iteration to
          // go the full number of atoms:
          if(!breakTies) stagnantTol=nActiveAtoms;

          bool done=indicesInPlay.empty();
          while(!done){
            //
            // do one round of iterations
            //
            numClasses = RankAtoms::iterateRanks2(nActiveAtoms,primeVect,atomicVect,
                                                  indicesInPlay,tadjMat,tranks,nRanks,
                                                  rankHistory,stagnantTol);

#ifdef VERBOSE_CANON
            BOOST_LOG(rdDebugLog)<< "************************ done outer iteration" << std::endl;
#endif  
#ifdef VERBOSE_CANON
            BOOST_LOG(rdDebugLog)<< "RANKS:" << std::endl;
            for(unsigned int tmpI=0;tmpI<tranks.size();tmpI++){
              BOOST_LOG(rdDebugLog)<< "\t\t" << tmpI << " " << tranks[tmpI] << std::endl;
            }
            BOOST_LOG(rdDebugLog)<< std::endl;
#endif    
            //
            // This is the tiebreaker stage of things
            //
            if( breakTies && !indicesInPlay.empty() && numClasses<nActiveAtoms){
              INT_VECT newRanks = tranks;

              // Add one to all ranks and multiply by two
              std::for_each(newRanks.begin(),newRanks.end(),_1=(_1+1)*2);
#ifdef VERBOSE_CANON
              BOOST_LOG(rdDebugLog)<< "postmult:" << std::endl;
              for(unsigned tmpI=0;tmpI<newRanks.size();tmpI++){
                BOOST_LOG(rdDebugLog)<< "\t\t" << newRanks[tmpI] << std::endl;
              }
              BOOST_LOG(rdDebugLog)<< std::endl;
#endif    

              //
              // find lowest duplicate rank with lowest invariant:
              //
              int lowestIdx=indicesInPlay.front();
              double lowestInvariant = tinvariants[lowestIdx];
              int lowestRank=newRanks[lowestIdx];
              for(INT_LIST_I ilIt=indicesInPlay.begin();
                  ilIt!=indicesInPlay.end();
                  ++ilIt){
                if(newRanks[*ilIt]<=lowestRank){
                  if(newRanks[*ilIt]<lowestRank ||
                     tinvariants[*ilIt] <= lowestInvariant){
                    lowestRank = newRanks[*ilIt];
                    lowestIdx = *ilIt;
                    lowestInvariant = tinvariants[*ilIt];
                  }
                }
              }

              //
              // subtract one from the lowest index, rerank and proceed
              //
              newRanks[lowestIdx] -= 1;
              RankAtoms::rankVect(newRanks,tranks);
              BOOST_FOREACH(int ilidx,indicesInPlay){
                nRanks[ilidx].push_back(ranks[ilidx]);
              }
#ifdef VERBOSE_CANON
              BOOST_LOG(rdDebugLog)<< "RE-RANKED ON:" << lowestIdx << std::endl;
              for(unsigned int tmpI=0;tmpI<newRanks.size();tmpI++){
                BOOST_LOG(rdDebugLog)<< "\t\t" << newRanks[tmpI] << " " << tranks[tmpI] << std::endl;
              }
              BOOST_LOG(rdDebugLog)<< std::endl;
#endif    
            } else {
              done = true;
            }
          }
          delete [] tadjMat;
        }
        unsigned int tidx=0;
        for(unsigned int aidx=0;aidx<nAtoms;++aidx){
          ranks[aidx]=0;
          if(atomsToUse[aidx]){
            ranks[aidx]=tranks[tidx++];
          }
        }
      }
    } // end of function rankAtomsInFragment

    
  } // end of namespace MolOps  
} // End Of RDKit namespace
