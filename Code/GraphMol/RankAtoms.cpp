// $Id: RankAtoms.cpp 4961 2006-02-18 00:14:47Z glandrum $
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RankAtoms.h>

#include <list>

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
    for(INT_VECT_CI ivCIt=ranks.begin();ivCIt!=ranks.end();ivCIt++){
      res.push_back(firstThousandPrimes[*ivCIt]);
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
      pos++;
      // now check to see if there is at least one more:
      if( std::find(pos,ranks.end(),ranks[*ivIt])==ranks.end()){
	INT_LIST::iterator tmpIt = ivIt;
	ivIt++;
	indicesInPlay.erase(tmpIt);
      } else {
	ivIt++;
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
	idxIt++){
      double accum;
      if(useSelf)
	accum=valVect[*idxIt];
      else
	accum=1.0;
      const unsigned int iTab = (*idxIt)*nAtoms;
      for(unsigned int j=0;j<nAtoms;j++){
	double elem=adjMat[iTab+j];
	if(elem>tol){
	  if(elem<2.-tol){
	    accum *= valVect[j];
	  } else {
	    accum *= pow(valVect[j],static_cast<int>(elem));
	  }
	}
      }
      res.push_back(accum);
    }
  }

  template <typename T>
  void debugVect(const std::vector<T> arg){
    typename std::vector<T>::const_iterator viIt;
    for(viIt=arg.begin();viIt!=arg.end();viIt++){
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
			    VECT_INT_VECT *rankHistory=0,unsigned int stagnantTol=1){
    PRECONDITION(!rankHistory||rankHistory->size()>=nAtoms,"bad rankHistory size");
    bool done = false;
    unsigned int numClasses = countClasses(ranks);
    unsigned int lastNumClasses = 0;
    unsigned int nCycles = 0;
    unsigned int nStagnant=0;
    if(stagnantTol<0){
      stagnantTol = nAtoms;
    }
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
	for(INT_LIST_CI idx=indicesInPlay.begin();idx!=indicesInPlay.end();idx++){
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
      calcAdjacentProducts(nAtoms,primeVect,adjMat,indicesInPlay,atomicVect);
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

}// end of RankAtoms namespace

namespace RDKit{

  namespace MolOps {
    // --------------------------------------------------
    //
    // Calculates chiral invariants for the atoms of a molecule
    //  These are based on Labute's proposal in:
    //  "An Efficient Algorithm for the Determination of Topological
    //   RS Chirality" Journal of the CCG (1996)
    //
    // --------------------------------------------------
    void buildChiralAtomInvariants(const ROMol &mol,INVAR_VECT &res){
      PRECONDITION(res.size()>=mol.getNumAtoms(),"res vect too small");
      unsigned int atsSoFar=0;
      for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++){
	Atom const *atom = *atIt;
	unsigned long invariant = 0;
	int num = atom->getAtomicNum() % 128;
	int nMult = (atom->getExplicitValence() -
		     atom->getNumExplicitHs() -
		     atom->getDegree()) % 8;
	// get an int with the deviation in the mass from the default:
	int mass = static_cast<int>(atom->getMass() -
				    PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
	mass += 8;
	if(mass < 0) mass = 0;
	else mass = mass % 16;

	int nHs = atom->getTotalNumHs() % 8;
	int deg = atom->getDegree() % 8;
	// we need to deviate a little bit from the Labute definition because
	// of the possibility of hydrogens:
	invariant = num; // 7 bits here
	invariant = (invariant << 4) | mass;
	invariant = (invariant << 3) | nMult;
	invariant = (invariant << 3) | deg;
	invariant = (invariant << 3) | nHs;
	//BOOST_LOG(rdDebugLog)<< "\t\tAt: " << atom->getIdx() << ": " << invariant << std::endl;
	res[atsSoFar++] = invariant;
      }
    }


    // --------------------------------------------------
    //
    // Calculates invariants for the atoms of a molecule
    //
    // NOTE: if the atom has not had chirality info pre-calculated, it doesn't
    // much matter what value includeChirality has!
    // --------------------------------------------------
    void buildAtomInvariants(const ROMol &mol,INVAR_VECT &res,
				     bool includeChirality){
      PRECONDITION(res.size()>=mol.getNumAtoms(),"res vect too small");
      unsigned int atsSoFar=0;

      for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++){
	Atom const *atom = *atIt;
	// FIX: this is just for debugging purposes:
	unsigned long invariant = atom->getIdx();
	int nHs = atom->getTotalNumHs() % 8;
	int chg = abs(atom->getFormalCharge()) % 8;
	int chgSign = atom->getFormalCharge() > 0;
	int num =    atom->getAtomicNum() % 128;
	int explicitVal = atom->getExplicitValence() % 16;
	int nConns = atom->getDegree() % 8;
	int deltaMass = static_cast<int>(atom->getMass() -
					 PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
	deltaMass += 8;
	if(deltaMass < 0) deltaMass = 0;
	else deltaMass = deltaMass % 16;

      
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

	invariant = 0;
	invariant = (invariant << 3) | nConns;
	invariant = (invariant << 4) | explicitVal; 
	invariant = (invariant << 7) | num;  
	invariant = (invariant << 4) | deltaMass;
	invariant = (invariant << 3) | nHs;
	invariant = (invariant << 4) | inRing;
	invariant = (invariant << 3) | chg;
	invariant = (invariant << 1) | chgSign;
	if(includeChirality ){
	  int isR=0;
	  // FIX this doesn't distinguish between S atoms and those w/o stereochem
	  if( atom->hasProp("_CIPCode")){
	    std::string cipCode;
	    atom->getProp("_CIPCode",cipCode);
	    if(cipCode=="R"){
	      isR=1;
	    }
	  }
	  invariant = (invariant << 1) | isR;
	}

	// now deal with cis/trans - this is meant to address issue 174
	// loop over the bonds on this atom and check if we have a double bond with
	// a chiral code marking 
	if (includeChirality) {
	  ROMol::OBOND_ITER_PAIR atomBonds = atom->getOwningMol().getAtomBonds(atom);
	  ROMol::GRAPH_MOL_BOND_PMAP::type pMap = atom->getOwningMol().getBondPMap();
	  int isT=0;
	  while (atomBonds.first != atomBonds.second){
	    Bond *tBond = pMap[*(atomBonds.first)];
	    if( (tBond->getBondType() == Bond::DOUBLE) &&
		(tBond->getStereo()>Bond::STEREOANY )) {
	      if (tBond->getStereo()==Bond::STEREOE) {
		isT = 1;
	      }
	      break;
	    }
	    atomBonds.first++;
	  }
	  invariant = (invariant << 1) | isT;
	}
	res[atsSoFar++] = invariant;
      }
    }
    // --------------------------------------------------
    //
    //  Daylight canonicalization, based up on algorithm described in
    //     JCICS 29, 97-101, (1989)
    //  When appropriate, specific references are made to the algorithm
    //  description in that paper.  Steps refer to Table III of the paper
    //
    // NOTE: if the atom has not had chirality info pre-calculated, it doesn't
    // much matter what value includeChirality has!
    // --------------------------------------------------
    void rankAtoms(const ROMol &mol,INT_VECT &ranks,
			   AtomInvariantType invariantType,
			   bool breakTies,
			   VECT_INT_VECT *rankHistory,
			   bool includeChirality){
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
	switch(invariantType){
	case ComprehensiveInvariants:
	  MolOps::buildAtomInvariants(mol,invariants,includeChirality);
	  break;
	case ChiralSearchInvariants:
	  MolOps::buildChiralAtomInvariants(mol,invariants);
	  break;
	}
      
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

	// start by ranking the atoms using the invariants
	ranks.resize(nAtoms);
	RankAtoms::rankVect(invariants,ranks);

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
	    numClasses = RankAtoms::iterateRanks(nAtoms,primeVect,atomicVect,
						 indicesInPlay,adjMat,ranks,
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
	      for(INT_VECT_I ivIt=newRanks.begin();ivIt!=newRanks.end();ivIt++)
		*ivIt = (*ivIt+1)*2;
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
	      for(INT_LIST_I ilIt=indicesInPlay.begin();
		  ilIt!=indicesInPlay.end();
		  ilIt++){
		if(newRanks[*ilIt]<=lowestRank){
		  if(newRanks[*ilIt]<lowestRank ||
		     invariants[*ilIt] <= lowestInvariant){
		    lowestRank = newRanks[*ilIt];
		    lowestIdx = *ilIt;
		    lowestInvariant = invariants[*ilIt];
		  }
		}
	      }

	      //
	      // subtract one from the lowest index, rerank and proceed
	      //
	      newRanks[lowestIdx] -= 1;
	      RankAtoms::rankVect(newRanks,ranks);
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
    }
  } // end of namespace MolOps  
} // End Of RDKit namespace
