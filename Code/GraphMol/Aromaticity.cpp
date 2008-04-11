// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Rings.h>
#include <RDGeneral/types.h>
#include <boost/dynamic_bitset.hpp>
#include <set>


/****************************************************************
Here are some molecule that have trouble us aromaticity wise

1. We get 2 aromatic rings for this molecule, but daylight says none. Also 
   if we replace [O+] with the N daylight says we have two aromatic rings
   [O-][N+]1=CC2=CC=C[O+]2[Cu]13[O+]4C=CC=C4C=[N+]3[O-]

2. We the fused ring system with all the rings in it is considered aroamtic by our code.
   This is because we count electrons from burried atoms as well when 
   we are dealing with fused rings

3. Here's a weird fused ring case, Pattern NAN, A overall
   O=C3C2=CC1=CC=COC1=CC2=CC=C3
 ************************************************************/
namespace RingUtils {
  using namespace RDKit;

  void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap,
                      INT_VECT &res,
                      boost::dynamic_bitset<> &done) {
    INT_INT_VECT_MAP::const_iterator pos= neighMap.find(curr);
    PRECONDITION(pos!=neighMap.end(),"bad argument");
    const INT_VECT &neighs = pos->second;
    INT_VECT_CI ni;
    done[curr] = 1;
    res.push_back(curr);
    for (ni = neighs.begin(); ni != neighs.end(); ni++) {
      if (!done[*ni]) {
        pickFusedRings((*ni), neighMap, res, done);
      }
    }
  }

  bool checkFused(const INT_VECT &rids, INT_INT_VECT_MAP &ringNeighs) {
    INT_INT_VECT_MAP_CI nci;
    int nrings = ringNeighs.size();
    boost::dynamic_bitset<> done(nrings);
    int rid;
    INT_VECT fused;

    // mark all rings in the system other than those in rids as done
    for (nci = ringNeighs.begin(); nci != ringNeighs.end(); nci++) {
      rid = (*nci).first;
      if (std::find(rids.begin(), rids.end(), rid) == rids.end()) {
        done[rid]=1;
      }
    }

    // then pick a fused system from the remaining (i.e. rids)
    // If the rings in rids are fused we should get back all of them 
    // in fused
    // if we get a smaller number in fused then rids are not fused
    pickFusedRings(rids.front(), ringNeighs, fused, done);

    CHECK_INVARIANT(fused.size() <= rids.size(), "");
    return (fused.size() == rids.size());
  }

  void makeRingNeighborMap(const VECT_INT_VECT &brings,
                           INT_INT_VECT_MAP &neighMap) {
    int nrings = brings.size();
    int i, j;
    INT_VECT ring1;
    for (i = 0; i < nrings; i++) {
      INT_VECT neighs;
      neighMap[i] = neighs;
    }

    for (i = 0; i < nrings; i++) {
      ring1 = brings[i];
      for (j = i + 1; j < nrings; j++) {
        INT_VECT inter;
        Intersect(ring1, brings[j], inter);
        if (inter.size() > 0) {
          neighMap[i].push_back(j);
          neighMap[j].push_back(i);
        }
      }
    }
  }
  

} // end of namespace RingUtils

// local utility namespace
namespace {
  using namespace RDKit;

  typedef enum {
    VacantElectronDonorType,
    OneElectronDonorType,
    TwoElectronDonorType,
    OneOrTwoElectronDonorType,
    AnyElectronDonorType,
    NoElectronDonorType
  } ElectronDonorType; // used in setting aromaticity
  typedef std::vector<ElectronDonorType> VECT_EDON_TYPE;
  typedef VECT_EDON_TYPE::iterator VECT_EDON_TYPE_I;
  typedef VECT_EDON_TYPE::const_iterator VECT_EDON_TYPE_CI;


  /******************************************************************************
   * SUMMARY:
   *  Apply Huckel rule to the specified ring and mark the bonds and atoms in the 
   *  ring accordingly.
   * 
   * ARGUMENTS:
   *  mol - molecule of interest
   *  ring - list of atoms that form the simple ring
   *  edon - list of electron donor type of the atoms in the molecule
   * 
   * RETURN :
   *  none
   * 
   * ASSUMPTIONS:
   *  the ring has to a simple ring and the electron donor type of each atom
   *  in the ring have been computed
   ******************************************************************************/
  static bool applyHuckel(ROMol &mol, const INT_VECT &ring,
                          const VECT_EDON_TYPE &edon);


  static void applyHuckelToFused(ROMol &mol, // molecule of interets
                                 const VECT_INT_VECT &srings, // list of all ring as atom IDS
                                 const VECT_INT_VECT &brings, // list of all rings as bond ids
                                 const INT_VECT &fused, // list of ring ids in the current fused system
                                 const VECT_EDON_TYPE &edon, // eletron donar state for each atom
                                 INT_INT_VECT_MAP &ringNeighs,
                                 int &narom // number of aromatic ring so far 
                                 );

  void markAtomsBondsArom(ROMol &mol, const VECT_INT_VECT &srings, 
                          const VECT_INT_VECT &brings,
                          const INT_VECT &ringIds,
                          std::set<unsigned int> &doneAtms) {
    INT_VECT aring, bring;
    INT_VECT_CI ri, ai, bi;

    // first mark the atoms in the rings
    for (ri = ringIds.begin(); ri != ringIds.end(); ri++) {
      
      aring = srings[*ri];
      
      // first mark the atoms in the ring
      for (ai = aring.begin(); ai != aring.end(); ai++) {
        mol.getAtomWithIdx(*ai)->setIsAromatic(true);
        doneAtms.insert(*ai);
      }
    }

    // mark the bonds
    // here we want to be more careful. We don't want to mark the fusing bonds 
    // as aromatic - only the outside bonds in a fused system are marked aromatic.
    // - loop through the rings and count the number of times each bond appears in 
    //   all the fused rings. 
    // - bonds that appeard only once are marked aromatic
    INT_MAP_INT bndCntr;
    INT_MAP_INT_I bci;

    for (ri = ringIds.begin(); ri != ringIds.end(); ri++) {
      bring = brings[*ri];
      for (bi = bring.begin(); bi != bring.end(); bi++) {
        if (bndCntr.find(*bi) == bndCntr.end()) {
          bndCntr[*bi] = 1;
        }
        else {
          bndCntr[*bi] += 1;
        }
      }
    }
    int bid;
    
    // now mark bonds that have a count of 1 to be aromatic;
    for (bci = bndCntr.begin(); bci != bndCntr.end(); bci++) {
      if ((*bci).second == 1) {
        bid = (*bci).first;
        mol.getBondWithIdx(bid)->setIsAromatic(true);
        mol.getBondWithIdx(bid)->setBondType(Bond::AROMATIC);
      }
    }
  }

  void getMinMaxAtomElecs(ElectronDonorType dtype, int &atlw, int &atup) {
    switch(dtype) {
    case AnyElectronDonorType :
      atlw = 0;
      atup = 2;
      break;
    case OneOrTwoElectronDonorType :
      atlw = 1;
      atup = 2;
      break;
    case OneElectronDonorType :
      atlw = atup = 1;
      break;
    case TwoElectronDonorType :
      atlw = atup = 2;
      break;
    case NoElectronDonorType :
    case VacantElectronDonorType :
    default :
      atlw = atup = 0;
      break;
    }
  }

  bool incidentNonCyclicMultipleBond(const Atom *at, int &who ) {
    PRECONDITION(at,"bad atom");
    // check if "at" has an non-cyclic multiple bond on it
    // if yes check which atom this bond goes to  
    // and record the atomID in who
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
    ROMol::GRAPH_MOL_BOND_PMAP::type pMap = at->getOwningMol().getBondPMap();
    while(beg!=end){
      if(!at->getOwningMol().getRingInfo()->numBondRings(pMap[*beg]->getIdx())){
        if (pMap[*beg]->getValenceContrib(at) >= 2.0) {
          who = pMap[*beg]->getOtherAtomIdx(at->getIdx());
          return true;
        }
       }
      beg++;
    }
    return false;
  }
   
  bool incidentCyclicMultipleBond(const Atom *at) {
    PRECONDITION(at,"bad atom");
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
    ROMol::GRAPH_MOL_BOND_PMAP::type pMap = at->getOwningMol().getBondPMap();
    while(beg!=end){
      if(at->getOwningMol().getRingInfo()->numBondRings(pMap[*beg]->getIdx())){
        if (pMap[*beg]->getValenceContrib(at) >= 2.0) {
          return true;
         }
       }
       beg++;
    }
    return false;
  }

  bool incidentMultipleBond(const Atom *at) {
    PRECONDITION(at,"bad atom");
    return at->getExplicitValence()!=static_cast<int>(at->getDegree()+at->getNumExplicitHs());
  }

  bool applyHuckel(ROMol &mol, const INT_VECT &ring,
                   const VECT_EDON_TYPE &edon) {
    int atlw, atup, rlw, rup, rie;
    bool aromatic = false;
    rlw = 0;
    rup = 0;
    INT_VECT_CI ivi;
    for (ivi = ring.begin(); ivi != ring.end(); ivi++) {
      ElectronDonorType edonType = edon[*ivi];
      getMinMaxAtomElecs(edonType, atlw, atup);
      rlw += atlw;
      rup += atup;
    }

    if (rup >= 6) {
      for (rie = rlw; rie <=rup; rie++) {
        if ( (rie - 2)%4 == 0) {
          aromatic = true;
          break;
        }
      }
    }
    return aromatic;
  }
    

  
  void applyHuckelToFused(ROMol &mol, // molecule of interets
                          const VECT_INT_VECT &srings, // list of all ring as atom IDS
                          const VECT_INT_VECT &brings, // list of all rings as bond ids
                          const INT_VECT &fused, // list of ring ids in the current fused system
                          const VECT_EDON_TYPE &edon, // eletron donar state for each atom
                          INT_INT_VECT_MAP &ringNeighs, // list of neighbors for eac candidate ring
                          int &narom // number of aromatic ring so far 
                          ) {

    // this function check huckel rule on a fused system it starts
    // with the individual rings in the system and then proceeds to
    // check for larger system i.e. if we have a 3 ring fused system,
    // huckel rule checked first on all teh 1 ring subsystems then 2
    // rung subsystems etc.

    INT_VECT aromRings;
    aromRings.resize(0);
    unsigned int nrings = fused.size();
    INT_VECT curRs;
    INT_VECT_CI mri;
    curRs.push_back(fused.front());
    int pos;
    unsigned int i, curSize = 0;
    INT_VECT comb;
    pos = -1;
    INT_VECT unionAtms;
    Union(srings, unionAtms);
    unsigned int nAtms = unionAtms.size();
    std::set<unsigned int> doneAtoms;
    while (1) {
      if (pos == -1 ) {
        
        curSize += 1;
        // check is we are done with all the atoms in the fused
        // system, if so quit. This is a fix for Issue252 REVIEW: is
        // this check sufficient or should we add an additional
        // contraint on the the number of combinations of rings in a
        // fused system that we will try. The number of combinations
        // can obviously be quite large when the number of rings in
        // the fused system is large
        if ((doneAtoms.size() == nAtms) || (curSize > nrings)) {
          break;
        }
        comb.resize(curSize);
        pos = 0;
        for (i =0; i < curSize; i++) {
          comb[i]=i;
        }
      }
      else {
        pos = nextCombination(comb, nrings);
      }

      if (pos == -1) {
        continue;
      }

      curRs.resize(0);
      for (i = 0; i < comb.size(); i++) {
        curRs.push_back(fused[comb[i]]);
      }

      // check if the picked subsystem is fused
      if (!RingUtils::checkFused(curRs, ringNeighs)) {
        continue;
      }

      // check aromaticity on the current fused system
      INT_VECT exclude;
      for (i = 0; i < srings.size(); i++) {
        if (std::find(curRs.begin(), curRs.end(),
                      static_cast<int>(i)) == curRs.end() ) {
          exclude.push_back(i);
        }
      }
      INT_VECT unon;
      Union(srings, unon, &exclude);
      
      if (applyHuckel(mol, unon, edon)) {
        // mark the atoms and bonds in these rings to be aromatic
        markAtomsBondsArom(mol, srings, brings, curRs, doneAtoms);

        // add the ring IDs to the aromatic rings found so far
        // avoid duplicates
        for (mri = curRs.begin(); mri != curRs.end(); mri++) {
          if (std::find(aromRings.begin(), aromRings.end(),
                        (*mri) ) == aromRings.end()) {
            aromRings.push_back(*mri);
          }
        }
      }// end check huckel rule
    } // end while(1)
    narom += aromRings.size();
  }

  bool isAtomCandForArom(const Atom *at, VECT_EDON_TYPE &edon) {
    PRECONDITION(at,"bad atom");
    // limit aromaticity to:
    //   - the first two rows of the periodic table
    //   - Se and Te
    if (at->getAtomicNum() > 18 &&
        at->getAtomicNum()!=34 &&
        at->getAtomicNum()!=52 ) {
      return false;
    }
    switch (edon[at->getIdx()]) {
    case VacantElectronDonorType:
    case OneElectronDonorType:
    case TwoElectronDonorType:
    case OneOrTwoElectronDonorType:
    case AnyElectronDonorType:
      break;
    default:
      return(false);      
    }

    // We are going to explicitly disallow atoms that have more
    // than one multiple bond. This is to handle the situation:
    //   C1=C=NC=N1 (sf.net bug 1934360) 
    int nUnsaturations=at->getExplicitValence()-at->getDegree();
    if(nUnsaturations>1){
      unsigned int nMult=0;
      const ROMol &mol = at->getOwningMol();
      ROMol::OEDGE_ITER beg,end;
      ROMol::GRAPH_MOL_BOND_PMAP::const_type pMap=mol.getBondPMap();
      boost::tie(beg,end) = mol.getAtomBonds(at);
      while(beg!=end){
        const Bond *bond=pMap[*beg];
        switch(bond->getBondType()){
        case Bond::SINGLE:
        case Bond::AROMATIC:
          break;
        case Bond::DOUBLE:
        case Bond::TRIPLE:
          ++nMult;
          break;
        default:
          // hopefully we had enough sense that we don't even
          // get here with these bonds... If we do land here,
          // just bail... I have no good answer for them.
          break;
        }
        if(nMult>1) break;
        ++beg;
      }
      if(nMult>1) return(false);
    }
    return(true);
  }

  ElectronDonorType getAtomDonorTypeArom(const Atom *at) {
    PRECONDITION(at,"bad atom");
    if(at->getAtomicNum()==0){
      // dummies can be anything:
      return AnyElectronDonorType;
    }

    ElectronDonorType res = NoElectronDonorType;
    int nelec = MolOps::countAtomElec(at);
    int who = -1;
    const ROMol &mol = at->getOwningMol();

    if (nelec < 0) {
      res = NoElectronDonorType;
    }
    else if (nelec == 0) {
      if (incidentNonCyclicMultipleBond(at, who)) {
        // This is border-line no electron to spare but may have an empty p-orbital 
        // Not sure if this should return vacantElectronDonorType 
        // FIX: explicitly doing this as a note for potential problems
        // 
        //res = NoElectronDonorType;
        
        res = VacantElectronDonorType;
      }
      else if (incidentCyclicMultipleBond(at)) {
        // no electron but has one in a in cycle multiple bond
        // FIX: why can't this be two electrons if there is a triple bond
        res = OneElectronDonorType; 
      }
      else {
        // no multiple bonds no electrons
        res = NoElectronDonorType;
      }
    }
    else if (nelec == 1) {
      if (incidentNonCyclicMultipleBond(at, who)) {
        // the only available electron is going to be from the
        // external multiple bond this electron will not be available
        // for aromaticity if this atom is bonded to a more electro
        // negative atom
        const Atom *at2 = mol.getAtomWithIdx(who);
        if (PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
                                                           at->getAtomicNum())) {
          res = VacantElectronDonorType;
        }
        else {
          res = OneElectronDonorType; 
        }
      }
      else {
        // require that the atom have at least one multiple bond
        if(incidentMultipleBond(at)){
          res = OneElectronDonorType;
        }
      }
    }
    else {
      if (incidentNonCyclicMultipleBond(at, who)) {
        // for cases with more than one electron :
        // if there is an incident multiple bond with an element that
        // is more electronegative than the this atom, count one less
        // electron
        const Atom *at2 = mol.getAtomWithIdx(who);
        if (PeriodicTable::getTable()->moreElectroNegative(at2->getAtomicNum(),
                                                           at->getAtomicNum())) {
          nelec--;
        }
      }
      if (nelec%2 == 1) {
        res = OneElectronDonorType;
      }
      else {
        res = TwoElectronDonorType;
      }
    }
    return(res);
  }
}// end of local utility namespace

namespace RDKit {
  namespace MolOps {
    int countAtomElec(const Atom *at) {
      PRECONDITION(at,"bad atom");
      int res;
      int dv; // default valence 
      int nle; // number of lone pair electrons
      int chg; // formal charge
      int tbo; // total bond order count
      int sbo; // number of bond == degree of the atom
  
      tbo = at->getExplicitValence() + at->getImplicitValence();
    
      chg = at->getFormalCharge();
      dv = PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum()); 
      if (dv <=1) {
        // univalent elements can't be either aromatic or conjugated
        return -1;
      }


      // number of lone pair electrons = (outer shell elecs) - (default valence)
      // FIX: this formula is probably wrong for higher order stuff
      nle = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum()) - dv; 

      sbo = at->getDegree() + at->getTotalNumHs();
      // if we are more than 3 coordinated we should not be aromatic
      if (sbo > 3) {
        return -1;
      }

      // num electrons available for donation 
      res = (dv + nle) - sbo - chg;
      return res;
    }


    int setAromaticity(RWMol &mol) {
      // FIX: we will assume for now that if the input molecule came
      // with aromaticity information it is correct and we will not
      // touch it. Loop through the atoms and check if any atom has
      // arom stuff set.  We may want check this more carefully later
      // and start from scratch if necessary
      ROMol::AtomIterator ai; 
      for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ai++) {
        if ((*ai)->getIsAromatic()) {
          // found aromatic info 
          return -1;
        }
      }
    
      // first find the all the simple rings in the molecule
      VECT_INT_VECT srings;
    
      symmetrizeSSSR(mol, srings);
      int narom = 0;
      // loop over all the atoms in the rings that can be candidates
      // for aromaticity
      // Atoms are candidates if 
      //   - it is part of ring
      //   - has one or more electron to donate or has empty p-orbitals
      int natoms = mol.getNumAtoms();
      boost::dynamic_bitset<> acands(natoms);
      VECT_EDON_TYPE edon(natoms);

      int ncands = 0;
      for (VECT_INT_VECT_I vivi = srings.begin();
           vivi != srings.end(); ++vivi) {
        for (INT_VECT_I ivi = (*vivi).begin();
             ivi!=(*vivi).end(); ++ivi) {
          Atom *at = mol.getAtomWithIdx(*ivi);
        
          // now that the atom is part of ring check if it can donate
          // electron or has empty orbitals. Record the donor type
          // information in 'edon' - we will need it when we get to
          // the Huckel rule later
          edon[*ivi] = getAtomDonorTypeArom(at);
          if (isAtomCandForArom(at, edon)) { 
            ncands++;
            acands[*ivi]=1;
          } else {
            acands[*ivi]=0;
          }
        }
      }

      // mark rings in the SSSR list that cannot be candidates for
      // aromaticity 
      // For a ring to be ac andidate for aromaticity: all the atoms
      // in the ring need to be aromatic 
      // Assume here that the SSSR
      // function above will return a continuous list of rings
      VECT_INT_VECT cRings; // holder for the candidate rings
      INT_VECT cring;
      for (VECT_INT_VECT_CI vivi = srings.begin();
           vivi != srings.end(); ++vivi) {
        bool rcand = true;
        for (INT_VECT_CI ivi = vivi->begin();
             ivi != vivi->end(); ++ivi) {
          if (!acands[*ivi]) { // if an atom in the ring is not a candidate
            rcand = false; // mark teh entire ring not a candidate
            break;
          }
        }
        if (rcand) {
          cring = (*vivi);
          cRings.push_back(cring);
        }
      }
      // first convert all rings to bonds ids
      VECT_INT_VECT brings;
      RingUtils::convertToBonds(cRings, brings, mol);
    
      // make the neighbor map for the rings 
      // i.e. a ring is a neighbor a another candidate ring if
      // shares at least one bond
      // useful to figure out fused systems 
      INT_INT_VECT_MAP neighMap;
      RingUtils::makeRingNeighborMap(brings, neighMap);

      // now loop over all the candidate rings and check the
      // huckel rule - of course paying attention to fused systems.
      INT_VECT doneRs;
      int curr = 0;
      int cnrs = cRings.size();
      boost::dynamic_bitset<> fusDone(cnrs);
      INT_VECT fused;
      while (curr < cnrs) {
        fused.resize(0);
        RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
        applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom);

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
    
      mol.setProp("numArom", narom, true);

      return narom;
    }
    

  }; // end of namespace MolOps
};  // end of namespace RDKit
