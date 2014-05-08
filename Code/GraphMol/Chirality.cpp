// $Id$
//
//  Copyright (C) 2004-2014 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RankAtoms.h>
#include <RDGeneral/types.h>
#include <sstream>
#include <algorithm>
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include <boost/dynamic_bitset.hpp>
#include <Geometry/point.h>

// #define VERBOSE_CANON 1


namespace RDKit{
  namespace Chirality {
    typedef std::pair<int,int> INT_PAIR;
    typedef std::vector<INT_PAIR> INT_PAIR_VECT;
    typedef std::vector<INT_PAIR>::iterator INT_PAIR_VECT_I;
    typedef std::vector<INT_PAIR>::const_iterator INT_PAIR_VECT_CI;

    typedef INT_VECT CIP_ENTRY;
    typedef std::vector< CIP_ENTRY > CIP_ENTRY_VECT;

    template <typename T>
    void debugVect(const std::vector<T> arg){
      typename std::vector<T>::const_iterator viIt;
      std::stringstream outS;
      for(viIt=arg.begin();viIt!=arg.end();viIt++){
        outS << *viIt << " ";
      }
      BOOST_LOG(rdDebugLog) << outS.str() << std::endl;
    }
  
    // --------------------------------------------------
    //
    // Calculates chiral invariants for the atoms of a molecule
    //  These are based on Labute's proposal in:
    //  "An Efficient Algorithm for the Determination of Topological
    //   RS Chirality" Journal of the CCG (1996)
    //
    // --------------------------------------------------
    void buildCIPInvariants(const ROMol &mol,DOUBLE_VECT &res){
      PRECONDITION(res.size()>=mol.getNumAtoms(),"res vect too small");
      int atsSoFar=0;
      //
      // NOTE:
      // If you make modifications to this, keep in mind that it is
      // essential that the initial comparison of ranks behave properly.
      // So, though it seems like it would makes sense to include
      // information about the number of Hs (or charge, etc) in the CIP
      // invariants, this will result in bad rankings.  For example, in
      // this molecule: OC[C@H](C)O, including the number of Hs would
      // cause the methyl group (atom 3) to be ranked higher than the CH2
      // connected to O (atom 1).  This is totally wrong.
      //
      // We also don't include any pre-existing stereochemistry information.
      // Though R and S assignments do factor in to the priorities of atoms,
      // we're starting here from scratch and we'll let the R and S stuff
      // be taken into account during the iterations.
      //
      for(ROMol::ConstAtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();++atIt){
        Atom const *atom = *atIt;
        unsigned long invariant = 0;
        int num = atom->getAtomicNum() % 128;
        // get an int with the deviation in the mass from the default:
        int mass = 0;
        if(atom->getIsotope()) {
          mass = atom->getIsotope()-PeriodicTable::getTable()->getMostCommonIsotope(atom->getAtomicNum());
          if(mass>=0) mass += 1;
        }          
        mass += 8;
        if(mass < 0) mass = 0;
        else mass = mass % 16;

#if 0
        // NOTE: the inclusion of hybridization in the invariant (as
        // suggested in the original paper), leads to the situation
        // that 
        //   C[C@@](O)(C=C)C(C)CC 
        // and
        //   C[C@@](O)(C=C)C(C)CO 
        // are assigned S chirality even though the rest of the world
        // seems to agree that they ought to be R (atom 3, sp2, is ranked
        // higher than atom 5, sp3, no matter what their environments)
        int hyb=0;
        switch(atom->getHybridization()) {
        case Atom::SP: hyb=6;break;
        case Atom::SP2: hyb=5;break;
        case Atom::SP3: hyb=1;break;
        case Atom::SP3D: hyb=3;break;
        case Atom::SP3D2: hyb=2;break;
        default: break;
        }
#endif

        invariant = num; // 7 bits here
        invariant = (invariant << 4) | mass;
      
        res[atsSoFar++] = invariant;
      }
    }

    void iterateCIPRanks(const ROMol &mol, DOUBLE_VECT &invars, INT_VECT &ranks,bool seedWithInvars){
      PRECONDITION(invars.size()==mol.getNumAtoms(),"bad invars size");
      PRECONDITION(ranks.size()>=mol.getNumAtoms(),"bad ranks size");

      unsigned int numAtoms = mol.getNumAtoms();
      CIP_ENTRY_VECT cipEntries(numAtoms);
      INT_LIST allIndices;
      for(unsigned int i=0;i<numAtoms;++i){
        allIndices.push_back(i);
      }
#ifdef VERBOSE_CANON
      BOOST_LOG(rdDebugLog) << "invariants:" << std::endl;
      for(unsigned int i=0;i<numAtoms;i++){
        BOOST_LOG(rdDebugLog) << i << ": " << invars[i] << std::endl;
      }
#endif  

      // rank those:
      RankAtoms::rankVect(invars,ranks);
#ifdef VERBOSE_CANON
      BOOST_LOG(rdDebugLog) << "initial ranks:" << std::endl;
      for(unsigned int i=0;i<numAtoms;++i){
        BOOST_LOG(rdDebugLog) << i << ": " << ranks[i] << std::endl;
      }
#endif  
      // Start each atom's rank vector with its atomic number:
      //  Note: in general one should avoid the temptation to
      //  use invariants here, those lead to incorrect answers
      for(unsigned int i=0;i<numAtoms;i++){
        if(!seedWithInvars){
          cipEntries[i].push_back(mol[i]->getAtomicNum());
          cipEntries[i].push_back(static_cast<int>(ranks[i]));
        } else {
          cipEntries[i].push_back(static_cast<int>(invars[i]));
        }
      }

      // Loop until either:
      //   1) all classes are uniquified
      //   2) the number of ranks doesn't change from one iteration to
      //      the next
      //   3) we've gone through maxIts times
      //      maxIts is calculated by dividing the number of atoms
      //      by 2. That's a pessimal version of the
      //      maximum number of steps required for two atoms to 
      //      "feel" each other (each influences one additional 
      //      neighbor shell per iteration). 
      unsigned int maxIts=numAtoms/2+1;
      unsigned int numIts=0;
      int lastNumRanks=-1;
      unsigned int numRanks=*std::max_element(ranks.begin(),ranks.end())+1;
      while( numRanks<numAtoms && numIts<maxIts && (lastNumRanks<0 ||
                                                    static_cast<unsigned int>(lastNumRanks)<numRanks) ){
        unsigned int longestEntry=0;
        // ----------------------------------------------------
        //
        // for each atom, get a sorted list of its neighbors' ranks:
        //
        for(INT_LIST_I it=allIndices.begin();
            it!=allIndices.end();
            ++it){
          CIP_ENTRY localEntry;
          localEntry.reserve(16);

          // start by pushing on our neighbors' ranks:
          ROMol::OEDGE_ITER beg,end;
          boost::tie(beg,end) = mol.getAtomBonds(mol[*it].get());
          while(beg!=end){
            const Bond *bond=mol[*beg].get();
            ++beg;
            unsigned int nbrIdx=bond->getOtherAtomIdx(*it);
            const Atom *nbr=mol[nbrIdx].get();
            
            int rank=ranks[nbrIdx]+1;
            // put the neighbor in 2N times where N is the bond order as a double.
            // this is to treat aromatic linkages on fair footing. i.e. at least in the
            // first iteration --c(:c):c and --C(=C)-C should look the same.
            // this was part of issue 3009911

            unsigned int count;
            if(bond->getBondType()==Bond::DOUBLE &&
               nbr->getAtomicNum()==15 &&
               (nbr->getDegree()==4 ||
                nbr->getDegree()==3) ) {
              // a special case for chiral phophorous compounds
              // (this was leading to incorrect assignment of
              // R/S labels ):
              count=1;

              // general justification of this is:
              // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
              // Multiple-Bond Unsaturation and Aromaticity". It contains several
              // conventions of which convention (b) is the one applying here:
              // "(b) Contibutions by d orbitals to bonds of quadriligant atoms are
              // neglected."
              // FIX: this applies to more than just P
            } else {
              count=static_cast<unsigned int>(floor(2.*bond->getBondTypeAsDouble()+.1));
            }
            CIP_ENTRY::iterator ePos=std::lower_bound(localEntry.begin(),localEntry.end(),rank);
            localEntry.insert(ePos,count,rank);
            ++nbr;
          }
          // add a zero for each coordinated H:
          // (as long as we're not a query atom)
          if(!mol[*it]->hasQuery()){
            localEntry.insert(localEntry.begin(),
                              mol[*it]->getTotalNumHs(),
                              0);
          }

          // we now have a sorted list of our neighbors' ranks,
          // copy it on in reversed order:
          cipEntries[*it].insert(cipEntries[*it].end(),
                                 localEntry.rbegin(),
                                 localEntry.rend());
          if(cipEntries[*it].size() > longestEntry){
            longestEntry = cipEntries[*it].size();
          }
        }
        // ----------------------------------------------------
        //
        // pad the entries so that we compare rounds to themselves:
        // 
        for(INT_LIST_I it=allIndices.begin();it!=allIndices.end();
            ++it){
          unsigned int sz=cipEntries[*it].size();
          if(sz<longestEntry){
            cipEntries[*it].insert(cipEntries[*it].end(),
                                   longestEntry-sz,
                                   -1);
          }
        }
        // ----------------------------------------------------
        //
        // sort the new ranks and update the list of active indices:
        // 
        lastNumRanks=numRanks;

        RankAtoms::rankVect(cipEntries,ranks);
        numRanks = *std::max_element(ranks.begin(),ranks.end())+1;

        // now truncate each vector and stick the rank at the end
        for(unsigned int i=0;i<numAtoms;++i){
          cipEntries[i][numIts+1]=ranks[i];
          cipEntries[i].erase(cipEntries[i].begin()+numIts+2,cipEntries[i].end());
        }
        
        ++numIts;
#ifdef VERBOSE_CANON
        BOOST_LOG(rdDebugLog) << "strings and ranks:" << std::endl;
        for(unsigned int i=0;i<numAtoms;i++){
          BOOST_LOG(rdDebugLog) << i << ": " << ranks[i] << " > ";
          debugVect(cipEntries[i]);
        }
#endif
      }
    }
    // Figure out the CIP ranks for the atoms of a molecule
    void assignAtomCIPRanks(const ROMol &mol, INT_VECT &ranks){
      PRECONDITION((!ranks.size() || ranks.size()>=mol.getNumAtoms()),
                   "bad ranks size");
      if(!ranks.size()) ranks.resize(mol.getNumAtoms());
      unsigned int numAtoms = mol.getNumAtoms();
      // get the initial invariants:
      DOUBLE_VECT invars(numAtoms,0);
      buildCIPInvariants(mol,invars);
      iterateCIPRanks(mol,invars,ranks,false);

      // copy the ranks onto the atoms:
      for(unsigned int i=0;i<numAtoms;++i){
        mol[i]->setProp("_CIPRank",ranks[i],1);
      }
    }
   

    // construct a vector with <atomIdx,direction> pairs for 
    // neighbors of a given atom.  This list will only be
    // non-empty if at least one of the bonds has its direction
    // set.
    void findAtomNeighborDirHelper(const ROMol &mol,const Atom *atom,
                                   const Bond *refBond,
                                   INT_VECT &ranks,
                                   INT_PAIR_VECT &neighbors,
                                   bool& hasExplicitUnknownStereo){
      PRECONDITION(atom,"bad atom");
      PRECONDITION(refBond,"bad bond");

      bool seenDir=false;
      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = mol.getAtomBonds(atom);
      while(beg!=end){
        const BOND_SPTR bond = mol[*beg];
        // check whether this bond is explictly set to have unknown stereo
        if (!hasExplicitUnknownStereo) {
          if (bond->hasProp("_UnknownStereo") &&
              bond->getProp<int>("_UnknownStereo")){
              hasExplicitUnknownStereo = true;
          }
        }

        Bond::BondDir dir=bond->getBondDir();
        if( bond->getIdx() != refBond->getIdx()){
          if(dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT){
            seenDir=true;
            // If we're considering the bond "backwards", (i.e. from end
            // to beginning, reverse the effective direction:
            if(atom != bond->getBeginAtom()){
              if(dir==Bond::ENDDOWNRIGHT) dir=Bond::ENDUPRIGHT;
              else dir=Bond::ENDDOWNRIGHT;
            }
          }
          Atom *nbrAtom=bond->getOtherAtom(atom);
          neighbors.push_back( std::make_pair(nbrAtom->getIdx(), dir) );
        }
        ++beg;
      }
      if(!seenDir) {
        neighbors.clear();
      } else {
        if( neighbors.size() == 2 &&
            ranks[neighbors[0].first] == ranks[neighbors[1].first] ){
          // the two substituents are identical, no stereochemistry here:
          neighbors.clear();
        } else {
          // it's possible that direction was set only one of the bonds, set the other
          // bond's direction to be reversed:
          if( neighbors[0].second != Bond::ENDDOWNRIGHT &&
              neighbors[0].second != Bond::ENDUPRIGHT ){
            CHECK_INVARIANT(neighbors.size()>1,"too few neighbors");
            neighbors[0].second =
              neighbors[1].second==Bond::ENDDOWNRIGHT ? Bond::ENDUPRIGHT : Bond::ENDDOWNRIGHT;
          } else if ( neighbors.size()>1 &&
                      neighbors[1].second != Bond::ENDDOWNRIGHT &&
                      neighbors[1].second != Bond::ENDUPRIGHT ){
            neighbors[1].second =
              neighbors[0].second==Bond::ENDDOWNRIGHT ? Bond::ENDUPRIGHT : Bond::ENDDOWNRIGHT;
          }
        }
      }
    }

    // find the neighbors for an atoms that are not connected by single bond that is not refBond
    // if checkDir is true only neighbor atoms with bonds marked with a direction will be returned
    void findAtomNeighborsHelper(const ROMol &mol,const Atom *atom,const Bond *refBond,
                                 INT_VECT &neighbors, bool checkDir=false) {
      PRECONDITION(atom,"bad atom");
      PRECONDITION(refBond,"bad bond");
      ROMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = mol.getAtomBonds(atom);
      while (beg != end) {
        const BOND_SPTR bond=mol[*beg];
        Bond::BondDir dir = bond->getBondDir();
        if (bond->getBondType()==Bond::SINGLE && bond->getIdx() != refBond->getIdx()) {
          if (checkDir) {
            if ((dir != Bond::ENDDOWNRIGHT) && (dir != Bond::ENDUPRIGHT)) {
              ++beg;
              continue;
            }
          }
          Atom *nbrAtom = bond->getOtherAtom(atom);
          neighbors.push_back(nbrAtom->getIdx());
        }
        ++beg;
      }
    }

    bool atomIsCandidateForRingStereochem(const ROMol &mol,const Atom *atom){
      PRECONDITION(atom,"bad atom");
      bool res=false;
      if(atom->hasProp("_ringStereochemCand")){
        res=atom->getProp<bool>("_ringStereochemCand");
      } else {
        const RingInfo *ringInfo=mol.getRingInfo();
        if(ringInfo->isInitialized() &&
           ringInfo->numAtomRings(atom->getIdx())){
          ROMol::OEDGE_ITER beg,end;
          boost::tie(beg,end) = mol.getAtomBonds(atom);
          std::vector<const Atom *> nonRingNbrs;
          std::vector<const Atom *> ringNbrs;
          while(beg!=end){
            const BOND_SPTR bond=mol[*beg];
            if(!ringInfo->numBondRings(bond->getIdx())){
              nonRingNbrs.push_back(bond->getOtherAtom(atom));
            } else {
              ringNbrs.push_back(bond->getOtherAtom(atom));
            }
            ++beg;
          }

          int rank1=0,rank2=0;
          switch(nonRingNbrs.size()){
          case 0:
            // don't do spiro:
            res=false;
            break;
          case 1:
            if(ringNbrs.size()==2) res=true;
            break;
          case 2:
            if( nonRingNbrs[0]->hasProp("_CIPRank") &&
                nonRingNbrs[1]->hasProp("_CIPRank") ){
              nonRingNbrs[0]->getProp("_CIPRank",rank1);
              nonRingNbrs[1]->getProp("_CIPRank",rank2);
              if(rank1==rank2){
                res=false;
              } else {
                res=true;
              }
            }
            break;
          default:
            res=false;
          }
        }
        atom->setProp("_ringStereochemCand",res,1);
      }
      return res;
    }

    // returns true if the atom is allowed to have stereochemistry specified
    bool checkChiralAtomSpecialCases(ROMol &mol,const Atom *atom){
      PRECONDITION(atom,"bad atom");

      if(!mol.getRingInfo()->isInitialized()){
        VECT_INT_VECT sssrs;
        MolOps::symmetrizeSSSR(mol, sssrs);
      }

      const RingInfo *ringInfo=mol.getRingInfo();
      if(ringInfo->numAtomRings(atom->getIdx()) &&
         atomIsCandidateForRingStereochem(mol,atom) ){
        // the atom is in a ring, so the "chirality" specification may actually
        // be handling ring stereochemistry.

        // check for another chiral tagged
        // atom without stereochem in this atom's rings:
        INT_VECT ringStereoAtoms(0);
        if(atom->hasProp("_ringStereoAtoms")){
          ringStereoAtoms=atom->getProp<INT_VECT>("_ringStereoAtoms");
        }
        const VECT_INT_VECT atomRings=ringInfo->atomRings();
        for(VECT_INT_VECT::const_iterator ringIt=atomRings.begin();
            ringIt!=atomRings.end();++ringIt){
          if(std::find(ringIt->begin(),ringIt->end(),
                       static_cast<int>(atom->getIdx()))!=ringIt->end()){
            for(INT_VECT::const_iterator idxIt=ringIt->begin();
                idxIt!=ringIt->end();++idxIt){
              int same=1;
              if(*idxIt!=static_cast<int>(atom->getIdx()) &&
                 mol.getAtomWithIdx(*idxIt)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
                 !mol.getAtomWithIdx(*idxIt)->hasProp("_CIPCode") &&
                 atomIsCandidateForRingStereochem(mol,mol.getAtomWithIdx(*idxIt)) ){
                // we get to keep the stereochem specification on this atom:
                if(mol.getAtomWithIdx(*idxIt)->getChiralTag()!=atom->getChiralTag()){
                  same=-1;
                }
                ringStereoAtoms.push_back(same*(*idxIt+1));
                INT_VECT oAtoms(0);
                if(mol.getAtomWithIdx(*idxIt)->hasProp("_ringStereoAtoms")){
                  oAtoms=mol.getAtomWithIdx(*idxIt)->getProp<INT_VECT>("_ringStereoAtoms");
                }
                oAtoms.push_back(same*(atom->getIdx()+1));
                mol.getAtomWithIdx(*idxIt)->setProp("_ringStereoAtoms",oAtoms,true);
              }
            }
          }
        }
        if(ringStereoAtoms.size()){
          atom->setProp("_ringStereoAtoms",ringStereoAtoms,true);
          return true;
        }
      }
      return false;
    }

    std::pair<bool,bool> isAtomPotentialChiralCenter(const Atom *atom,const ROMol &mol,const INT_VECT &ranks,
                                     Chirality::INT_PAIR_VECT &nbrs){
      // loop over all neighbors and form a decorated list of their
      // ranks:
      bool legalCenter=true;
      bool hasDupes=false;

      if(atom->getTotalDegree()>4){
        // we only know tetrahedral chirality
        legalCenter=false;
      } else {
        boost::dynamic_bitset<> codesSeen(mol.getNumAtoms());
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(atom);
        while(beg!=end){
          unsigned int otherIdx=mol[*beg]->getOtherAtom(atom)->getIdx();
          CHECK_INVARIANT(ranks[otherIdx]<static_cast<int>(mol.getNumAtoms()),
                          "CIP rank higher than the number of atoms.");
          // watch for neighbors with duplicate ranks, which would mean
          // that we cannot be chiral:
          if(codesSeen[ranks[otherIdx]]){
            // we've already seen this code, it's a dupe
            hasDupes = true;
            break;
          }
          codesSeen[ranks[otherIdx]]=1;
          nbrs.push_back(std::make_pair(ranks[otherIdx],
                                        mol[*beg]->getIdx()));
          ++beg;
        }

        // figure out if this is a legal chiral center or not:
        if(!hasDupes){
          if(nbrs.size()<3){
            // less than three neighbors is never stereogenic
            legalCenter=false;
          } else if(nbrs.size()==3){
            // three-coordinate with a single H we'll accept automatically:
            if(atom->getTotalNumHs()!=1){
              // otherwise we default to not being a legal center
              legalCenter=false;
              // but there are a few special cases we'll accept
              // sulfur or selenium with either a positive charge or a double bond:
              if((atom->getAtomicNum()==16||atom->getAtomicNum()==34) &&
                 (atom->getExplicitValence()==4 ||
                  (atom->getExplicitValence()==3 && atom->getFormalCharge()==1))) {
                legalCenter=true;
              }
            }
          }
        }
      }
      return std::make_pair(legalCenter,hasDupes);
    }
    
    // returns a pair:
    //   1) are there unassigned stereoatoms
    //   2) did we assign any?
    std::pair<bool,bool> assignAtomChiralCodes(ROMol &mol,INT_VECT &ranks,
                                               bool flagPossibleStereoCenters){
      PRECONDITION( (!ranks.size() || ranks.size()==mol.getNumAtoms()),
                    "bad rank vector size");
      bool atomChanged=false;
      unsigned int unassignedAtoms=0;

      // ------------------
      // now loop over each atom and, if it's marked as chiral,
      //  figure out the appropriate CIP label:
      for(ROMol::AtomIterator atIt=mol.beginAtoms();
          atIt!=mol.endAtoms();++atIt){
        Atom *atom=*atIt;
        Atom::ChiralType tag=atom->getChiralTag();

        // only worry about this atom if it has a marked chirality
        // we understand:
        if(flagPossibleStereoCenters || (tag != Atom::CHI_UNSPECIFIED &&
                                         tag != Atom::CHI_OTHER) ){
          if(atom->hasProp("_CIPCode")){
            continue;
          }

          if(!ranks.size()){
            //  if we need to, get the "CIP" ranking of each atom:
            assignAtomCIPRanks(mol,ranks);
          }
          Chirality::INT_PAIR_VECT nbrs;
          bool legalCenter,hasDupes;
          boost::tie(legalCenter,hasDupes)=isAtomPotentialChiralCenter(atom,mol,ranks,nbrs);
          if(legalCenter){
            ++unassignedAtoms;
          }
          if(legalCenter && !hasDupes && flagPossibleStereoCenters){
            atom->setProp("_ChiralityPossible",1);
          }

          if( legalCenter && !hasDupes &&
              tag != Atom::CHI_UNSPECIFIED &&
              tag != Atom::CHI_OTHER ) {
            // stereochem is possible and we have no duplicate neighbors, assign
            // a CIP code:
            atomChanged=true;
            --unassignedAtoms;

            // sort the list of neighbors by their CIP ranks:
            std::sort(nbrs.begin(),nbrs.end(),RankAtoms::pairLess<int,int>());

            // collect the list of neighbor indices:
            std::list<int> nbrIndices;
            for(Chirality::INT_PAIR_VECT_CI nbrIt=nbrs.begin();
                nbrIt!=nbrs.end(); ++nbrIt){
              nbrIndices.push_back((*nbrIt).second);
            }
            // ask the atom how many swaps we have to make:
            int nSwaps = atom->getPerturbationOrder(nbrIndices);

            // if the atom has 3 neighbors and a hydrogen, add a swap:
            if(nbrIndices.size()==3 && atom->getTotalNumHs()==1){
              ++nSwaps;
            }
          
            // if that number is odd, we'll change our chirality:
            if(nSwaps%2){
              if(tag == Atom::CHI_TETRAHEDRAL_CCW) tag=Atom::CHI_TETRAHEDRAL_CW;
              else tag=Atom::CHI_TETRAHEDRAL_CCW;
            }
            // now assign the CIP code:
            std::string cipCode;
            if(tag==Atom::CHI_TETRAHEDRAL_CCW) cipCode="S";
            else cipCode="R";
            atom->setProp("_CIPCode",cipCode,true);
          }
        }
      }
      return std::make_pair((unassignedAtoms>0),atomChanged);
    }

    // returns a pair:
    //   1) are there unassigned stereo bonds?
    //   2) did we assign any?
    std::pair<bool,bool> assignBondStereoCodes(ROMol &mol,INT_VECT &ranks){
      PRECONDITION( (!ranks.size() || ranks.size()==mol.getNumAtoms()),
                    "bad rank vector size");
      bool assignedABond=false;
      unsigned int unassignedBonds=0;

      // find the double bonds:
      for(ROMol::BondIterator bondIt=mol.beginBonds();
          bondIt!=mol.endBonds();
          ++bondIt){
        if( (*bondIt)->getBondType()==Bond::DOUBLE ){
          Bond *dblBond=*bondIt;
          if(dblBond->getStereo()!=Bond::STEREONONE){
            continue;
          }
          if(!ranks.size()){
            assignAtomCIPRanks(mol,ranks);
          }
          dblBond->getStereoAtoms().clear();

          // at the moment we are ignoring stereochem on ring bonds with less than 8 members.
          if(!mol.getRingInfo()->numBondRings(dblBond->getIdx()) ||
             mol.getRingInfo()->minBondRingSize(dblBond->getIdx())>7 ){
            const Atom *begAtom=dblBond->getBeginAtom();
            const Atom *endAtom=dblBond->getEndAtom();
            // we're only going to handle 2 or three coordinate atoms:
            if( (begAtom->getDegree()==2 || begAtom->getDegree()==3) &&
                (endAtom->getDegree()==2 || endAtom->getDegree()==3) ){
              ++unassignedBonds;
                
              // look around each atom and see if it has at least one bond with
              // direction marked:
          
              // the pairs here are: atomrank,bonddir
              Chirality::INT_PAIR_VECT begAtomNeighbors,endAtomNeighbors;
              bool hasExplicitUnknownStereo = false;
              if((dblBond->getBeginAtom()->hasProp("_UnknownStereo") &&
                  dblBond->getBeginAtom()->getProp<int>("_UnknownStereo")) ||
                 (dblBond->getEndAtom()->hasProp("_UnknownStereo") &&
                  dblBond->getEndAtom()->getProp<int>("_UnknownStereo")) ){
                hasExplicitUnknownStereo=true;
              }
              Chirality::findAtomNeighborDirHelper(mol,begAtom,dblBond,
                                                   ranks,begAtomNeighbors,
                                                   hasExplicitUnknownStereo);
              Chirality::findAtomNeighborDirHelper(mol,endAtom,dblBond,
                                                   ranks,endAtomNeighbors,
                                                   hasExplicitUnknownStereo);

              if(begAtomNeighbors.size() && endAtomNeighbors.size()){
                // Each atom has at least one neighboring bond with marked
                // directionality.  Find the highest-ranked directionality
                // on each side:

                int begDir,endDir, endNbrAid, begNbrAid;
                if(begAtomNeighbors.size()==1 ||
                   ranks[begAtomNeighbors[0].first] >
                   ranks[begAtomNeighbors[1].first] ){
                  begDir = begAtomNeighbors[0].second;
                  begNbrAid = begAtomNeighbors[0].first;
                } else {
                  begDir = begAtomNeighbors[1].second;
                  begNbrAid = begAtomNeighbors[1].first;
                }
                if(endAtomNeighbors.size()==1 ||
                   ranks[endAtomNeighbors[0].first] >
                   ranks[endAtomNeighbors[1].first]){
                  endDir = endAtomNeighbors[0].second;
                  endNbrAid = endAtomNeighbors[0].first;
                } else {
                  endDir = endAtomNeighbors[1].second;
                  endNbrAid = endAtomNeighbors[1].first;
                }
                dblBond->getStereoAtoms().push_back(begNbrAid);
                dblBond->getStereoAtoms().push_back(endNbrAid);
                if (hasExplicitUnknownStereo) {
                  dblBond->setStereo(Bond::STEREOANY);
                  assignedABond=true;
                }
                else if( begDir == endDir ){
                  // In findAtomNeighborDirHelper, we've set up the
                  // bond directions here so that they correspond to
                  // having both single bonds START at the double bond.
                  // This means that if the single bonds point in the same
                  // direction, the bond is cis, "Z"
                  dblBond->setStereo(Bond::STEREOZ);
                  assignedABond=true;
                } else {
                  dblBond->setStereo(Bond::STEREOE);
                  assignedABond=true;
                }
                --unassignedBonds;
              }
            }
          }
        }
      }
      return std::make_pair(unassignedBonds>0,assignedABond);
    }

    // reassign atom ranks by supplementing the current ranks
    // with information about known chirality
    void rerankAtoms(const ROMol &mol, INT_VECT &ranks) {
      PRECONDITION(ranks.size()==mol.getNumAtoms(),"bad rank vector size");
      unsigned int factor=100;
      while(factor<mol.getNumAtoms()) factor*=10;

#ifdef VERBOSE_CANON
      BOOST_LOG(rdDebugLog) << "rerank PRE: " << std::endl;
      for(int i=0;i<mol.getNumAtoms();i++){
        BOOST_LOG(rdDebugLog) << "  "<< i << ": " << ranks[i] << std::endl;
      }
#endif

      DOUBLE_VECT invars(mol.getNumAtoms());
      // and now supplement them:
      for(unsigned int i=0;i<mol.getNumAtoms();++i){
        invars[i] = ranks[i]*factor;
        const Atom *atom=mol.getAtomWithIdx(i);
        // Priority order: R > S > nothing
        if(atom->hasProp("_CIPCode")){
          std::string cipCode;
          atom->getProp("_CIPCode",cipCode);
          if(cipCode=="S"){
            invars[i]+=10;
          } else if(cipCode=="R"){
            invars[i]+=20;
          }
        }
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = mol.getAtomBonds(atom);
        while(beg!=end){
          const BOND_SPTR oBond=mol[*beg];
          if(oBond->getBondType()==Bond::DOUBLE){
            if(oBond->getStereo()==Bond::STEREOE){
              invars[i]+=1;
            } else if(oBond->getStereo()==Bond::STEREOZ){
              invars[i]+=2;
            }
          }
          ++beg;
        }
      }
      iterateCIPRanks(mol,invars,ranks,true);
      // copy the ranks onto the atoms:
      for(unsigned int i=0;i<mol.getNumAtoms();i++){
        mol.getAtomWithIdx(i)->setProp("_CIPRank",ranks[i],1);
      }

#ifdef VERBOSE_CANON
      BOOST_LOG(rdDebugLog) << "   post: " << std::endl;
      for(int i=0;i<mol.getNumAtoms();i++){
        BOOST_LOG(rdDebugLog) << "  "<< i << ": " << ranks[i] << std::endl;
      }
#endif

    }
  } // end of chirality namespace

  namespace MolOps {


    /*
        We're going to do this iteratively:
          1) assign atom stereochemistry
          2) assign bond stereochemistry
          3) if there are still unresolved atoms or bonds
             repeat the above steps as necessary
     */
    void assignStereochemistry(ROMol &mol,bool cleanIt,bool force,bool flagPossibleStereoCenters){
      if(!force && mol.hasProp("_StereochemDone")){
        return;
      }

      // later we're going to need ring information, get it now if we don't
      // have it already:
      if(!mol.getRingInfo()->isInitialized()){
        MolOps::symmetrizeSSSR(mol);
      }


#if 0
      std::cerr<<">>>>>>>>>>>>>\n";
      std::cerr<<"assign stereochem\n";
      mol.debugMol(std::cerr);
#endif
      
      if(cleanIt){
        for(ROMol::AtomIterator atIt=mol.beginAtoms();
            atIt!=mol.endAtoms();++atIt){
          if((*atIt)->hasProp("_CIPCode")){
            (*atIt)->clearProp("_CIPCode");
          }
          if((*atIt)->hasProp("_ChiralityPossible")){
            (*atIt)->clearProp("_ChiralityPossible");
          }
        }        
        for(ROMol::BondIterator bondIt=mol.beginBonds();
            bondIt!=mol.endBonds();
            ++bondIt){
          if( (*bondIt)->getBondType()==Bond::DOUBLE &&
	      (*bondIt)->getStereo() != Bond::STEREOANY ){
            (*bondIt)->setStereo(Bond::STEREONONE);
            (*bondIt)->getStereoAtoms().clear();
          }
        }
      }
      INT_VECT atomRanks;
      bool keepGoing=true;
      bool hasStereoAtoms=true,changedStereoAtoms;
      bool hasStereoBonds=true,changedStereoBonds;
      while(keepGoing){
        if(hasStereoAtoms){
          boost::tie(hasStereoAtoms,changedStereoAtoms) = Chirality::assignAtomChiralCodes(mol,atomRanks,
                                                                                           flagPossibleStereoCenters);
        } else {
          changedStereoAtoms=false;
        }
        if(hasStereoBonds){
          boost::tie(hasStereoBonds,changedStereoBonds) = Chirality::assignBondStereoCodes(mol,atomRanks);
        } else {
          changedStereoBonds=false;
        }
        keepGoing=(hasStereoAtoms||hasStereoBonds) && (changedStereoAtoms||changedStereoBonds);

        if(keepGoing){
          // update the atom ranks based on the new information we have:
          Chirality::rerankAtoms(mol,atomRanks);
        }
#if 0
        std::cout<<"*************** done iteration "<<keepGoing<<" ***********"<<std::endl;
        mol.debugMol(std::cout);
        std::cout<<"*************** done iteration "<<keepGoing<<" ***********"<<std::endl;
#endif
      }

      if(cleanIt){
        for(ROMol::AtomIterator atIt=mol.beginAtoms();
            atIt!=mol.endAtoms();++atIt){
          if((*atIt)->hasProp("_ringStereochemCand")) (*atIt)->clearProp("_ringStereochemCand");
          if((*atIt)->hasProp("_ringStereoAtoms")) (*atIt)->clearProp("_ringStereoAtoms");
        }
        for(ROMol::AtomIterator atIt=mol.beginAtoms();
            atIt!=mol.endAtoms();++atIt){
          Atom *atom=*atIt;
          if(atom->getChiralTag()!=Atom::CHI_UNSPECIFIED
             && !atom->hasProp("_CIPCode") &&
             !Chirality::checkChiralAtomSpecialCases(mol,atom) ){
            atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            
            // If the atom has an explicit hydrogen and no charge, that H
            // was probably put there solely because of the chirality.
            // So we'll go ahead and remove it.
            // This was Issue 194
            if(atom->getNumExplicitHs()==1 &&
               atom->getFormalCharge()==0 &&
               !atom->getIsAromatic() ){
              atom->setNumExplicitHs(0);
              atom->setNoImplicit(false);
              atom->calcExplicitValence(false);
              atom->calcImplicitValence(false);
            }
          }
        }        
        for(ROMol::BondIterator bondIt=mol.beginBonds();
            bondIt!=mol.endBonds();++bondIt){
          // wedged bonds to atoms that have no stereochem 
          // should be removed. (github issue 87)
          if(((*bondIt)->getBondDir()==Bond::BEGINWEDGE ||
              (*bondIt)->getBondDir()==Bond::BEGINDASH) &&
             (*bondIt)->getBeginAtom()->getChiralTag()==Atom::CHI_UNSPECIFIED &&
             (*bondIt)->getEndAtom()->getChiralTag()==Atom::CHI_UNSPECIFIED){
            (*bondIt)->setBondDir(Bond::NONE);
          }
        }
      }
      mol.setProp("_StereochemDone",1,true);

#if 0
      std::cerr<<"---\n";
      mol.debugMol(std::cerr);
      std::cerr<<"<<<<<<<<<<<<<<<<\n";
#endif      

    }

    // Find bonds than can be cis/trans in a molecule and mark them as "any"
    // - this function finds any double bonds that can potentially be part 
    //   of a cis/trans system. No attempt is made here to mark them cis or trans
    // 
    // This function is useful in two situations
    //  1) when parsing a mol file; for the bonds marked here, coordinate informations 
    //     on the neighbors can be used to indentify cis or trans states
    //  2) when writing a mol file; bonds that can be cis/trans but not marked as either 
    //     need to be specially marked in the mol file
    //
    //  The CIPranks on the neighboring atoms are check in this function. The _CIPCode property
    //  if set to any on the double bond.
    // 
    // ARGUMENTS:
    //   mol - the molecule of interest
    //   cleanIt - if this option is set to true, any previous marking of _CIPCode 
    //               on the bond is cleared - otherwise it is left untouched
    void findPotentialStereoBonds(ROMol &mol,bool cleanIt) {
      // FIX: The earlier thought was to provide an optional argument to ignore or consider
      //  double bonds in a ring. But I am removing this optional argument and ignoring ring bonds 
      //  completely for now. This is because finding a potential stereo bond in a ring involves
      //  more than just checking the CIPranks for the neighbors - SP 05/04/04

      // make this function callable multiple times
      if ((mol.hasProp("_BondsPotentialStereo")) && (!cleanIt)) {
        return;
      } else {
        INT_VECT ranks;
        ranks.resize(mol.getNumAtoms());
        bool cipDone=false;

        ROMol::BondIterator bondIt;
        for(bondIt=mol.beginBonds(); bondIt!=mol.endBonds(); ++bondIt){
          if( (*bondIt)->getBondType()==Bond::DOUBLE  && 
              !(mol.getRingInfo()->numBondRings((*bondIt)->getIdx())) ) {
            // we are ignoring ring bonds here - read the FIX above
            Bond *dblBond=*bondIt;
            // if the bond is flagged as EITHERDOUBLE, we ignore it:
            if(dblBond->getBondDir()==Bond::EITHERDOUBLE ||
	       dblBond->getStereo()==Bond::STEREOANY ){
              break;
            }
            // proceed only if we either want to clean the stereocode on this bond
            // or if none is set on it yet
            if ( cleanIt || dblBond->getStereo()==Bond::STEREONONE ) {
              dblBond->setStereo(Bond::STEREONONE);
              const Atom *begAtom=dblBond->getBeginAtom(),*endAtom=dblBond->getEndAtom();
              // we're only going to handle 2 or three coordinate atoms:
              if( (begAtom->getDegree()==2 || begAtom->getDegree()==3) &&
                  (endAtom->getDegree()==2 || endAtom->getDegree()==3) ){
                // ------------------
                // get the CIP ranking of each atom if we need it:
                if(!cipDone){
                  Chirality::assignAtomCIPRanks(mol,ranks);
                  cipDone=true;
                }
                // find the neighbors for the begin atom and the endAtom
                INT_VECT begAtomNeighbors,endAtomNeighbors;
                Chirality::findAtomNeighborsHelper(mol,begAtom,dblBond,begAtomNeighbors);
                Chirality::findAtomNeighborsHelper(mol,endAtom,dblBond,endAtomNeighbors);
                if(begAtomNeighbors.size()>0 && endAtomNeighbors.size()>0){
                  if ((begAtomNeighbors.size() == 2) && (endAtomNeighbors.size() == 2) ) {
                    // if both of the atoms have 2 neighbors (other than the one connected
                    // by the double bond) and ....
                    if ( (ranks[begAtomNeighbors[0]] != ranks[begAtomNeighbors[1]]) &&
                         (ranks[endAtomNeighbors[0]] != ranks[endAtomNeighbors[1]]) ) {
                      // the neighbors ranks are different at both the ends,
                      // this bond can be part of a cis/trans system
                      if(ranks[begAtomNeighbors[0]] > ranks[begAtomNeighbors[1]]){
                        dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                      } else {
                        dblBond->getStereoAtoms().push_back(begAtomNeighbors[1]);
                      }
                      if(ranks[endAtomNeighbors[0]] > ranks[endAtomNeighbors[1]]){
                        dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                      } else {
                        dblBond->getStereoAtoms().push_back(endAtomNeighbors[1]);
                      }
                    }
                  } else if (begAtomNeighbors.size() == 2) {
                    // if the begAtom has two neighbors and ....
                    if (ranks[begAtomNeighbors[0]] != ranks[begAtomNeighbors[1]]) {
                      // their ranks are different
                      if(ranks[begAtomNeighbors[0]] > ranks[begAtomNeighbors[1]]){
                        dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                      } else {
                        dblBond->getStereoAtoms().push_back(begAtomNeighbors[1]);
                      }
                      dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                    }
                  } else if (endAtomNeighbors.size() == 2) {
                    // if the endAtom has two neighbors and ...
                    if (ranks[endAtomNeighbors[0]] != ranks[endAtomNeighbors[1]]) {
                      // their ranks are different
                      dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                      if(ranks[endAtomNeighbors[0]] > ranks[endAtomNeighbors[1]]){
                        dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                      } else {
                        dblBond->getStereoAtoms().push_back(endAtomNeighbors[1]);
                      }
                    }
                  } else {
                    // end and beg atoms has only one neighbor each, it doesn't matter what the ranks are:
                    dblBond->getStereoAtoms().push_back(begAtomNeighbors[0]);
                    dblBond->getStereoAtoms().push_back(endAtomNeighbors[0]);
                  } // end of different number of neighbors on beg and end atoms
                } // end of check that beg and end atoms have at least 1 neighbor:
              }// end of 2 and 3 coordinated atoms only
            } // end of we want it or CIP code is not set
          } // end of double bond
        } // end of for loop over all bonds
        mol.setProp("_BondsPotentialStereo", 1, true);
      }
    }

    // removes chirality markers from sp and sp2 hybridized centers:
    void cleanupChirality(RWMol &mol){
      for(ROMol::AtomIterator atomIt=mol.beginAtoms();
          atomIt!=mol.endAtoms();
          ++atomIt){
        if( (*atomIt)->getChiralTag()!=Atom::CHI_UNSPECIFIED &&
            (*atomIt)->getHybridization() < Atom::SP3 ){
          (*atomIt)->setChiralTag(Atom::CHI_UNSPECIFIED);
        }
      }
    }

    void assignChiralTypesFrom3D(ROMol &mol,int confId,bool replaceExistingTags){
      const double ZERO_VOLUME_TOL=0.1;
      if(!mol.getNumConformers()) return;
      const Conformer &conf=mol.getConformer(confId);
      if(!conf.is3D()) return;

      // if the molecule already has stereochemistry
      // perceived, remove the flags that indicate
      // this... what we're about to do will require
      // that we go again.
      if(mol.hasProp("_StereochemDone")){
        mol.clearProp("_StereochemDone");
      }
      
      for(ROMol::AtomIterator atomIt=mol.beginAtoms();atomIt!=mol.endAtoms();++atomIt){
        Atom *atom=*atomIt;
        // if we aren't replacing existing tags and the atom is already tagged, punt:
        if(!replaceExistingTags && atom->getChiralTag()!=Atom::CHI_UNSPECIFIED){
          continue;
        }
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        // additional reasons to skip the atom:
        if(atom->getDegree()<3 || // not enough explicit neighbors
           atom->getTotalDegree()!=4 ||  // not enough total neighbors
           atom->getTotalNumHs(true)>1 // more than two Hs
           ){
          continue;
        }
        const RDGeom::Point3D &p0=conf.getAtomPos(atom->getIdx());
        ROMol::ADJ_ITER nbrIdx,endNbrs;
        boost::tie(nbrIdx,endNbrs) = mol.getAtomNeighbors(atom);
        const RDGeom::Point3D &p1=conf.getAtomPos(*nbrIdx);
        ++nbrIdx;
        const RDGeom::Point3D &p2=conf.getAtomPos(*nbrIdx);
        ++nbrIdx;
        const RDGeom::Point3D &p3=conf.getAtomPos(*nbrIdx);

        RDGeom::Point3D v1=p1-p0;
        RDGeom::Point3D v2=p2-p0;
        RDGeom::Point3D v3=p3-p0;
        
        double chiralVol= v1.dotProduct(v2.crossProduct(v3));
        if(chiralVol<-ZERO_VOLUME_TOL){
          atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        } else if (chiralVol>ZERO_VOLUME_TOL){
          atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
        } else {
          atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        }
      }

    }

    void removeStereochemistry(ROMol &mol){
      if(mol.hasProp("_StereochemDone")){
        mol.clearProp("_StereochemDone");
      }
      for(ROMol::AtomIterator atIt=mol.beginAtoms();
          atIt!=mol.endAtoms();++atIt){
        (*atIt)->setChiralTag(Atom::CHI_UNSPECIFIED);
        if((*atIt)->hasProp("_CIPCode")){
          (*atIt)->clearProp("_CIPCode");
        }
        if((*atIt)->hasProp("_CIPRank")){
          (*atIt)->clearProp("_CIPRank");
        }

      }        
      for(ROMol::BondIterator bondIt=mol.beginBonds();
          bondIt!=mol.endBonds();
          ++bondIt){
        if( (*bondIt)->getBondType()==Bond::DOUBLE ){
          (*bondIt)->setStereo(Bond::STEREONONE);
          (*bondIt)->getStereoAtoms().clear();
        } else if( (*bondIt)->getBondType()==Bond::SINGLE ){
          (*bondIt)->setBondDir(Bond::NONE);
        }
      }
    }
  }  // end of namespace MolOps
}  // end of namespace RDKit
