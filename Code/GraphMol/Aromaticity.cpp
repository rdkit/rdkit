// $Id$
//
//  Copyright (C) 2003-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Rings.h>
#include <RDGeneral/types.h>
#include <boost/dynamic_bitset.hpp>
#include <set>

// introduced for the sake of efficiency
// this is the maximum ring size that will be considered
// as a candiate for fused-ring aromaticity. This is picked to
// be a bit bigger than the outer ring in a porphyrin
// This came up while fixing sf.net issue249
const unsigned int maxFusedAromaticRingSize=24;

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
                      boost::dynamic_bitset<> &done,
                      int depth){
    INT_INT_VECT_MAP::const_iterator pos= neighMap.find(curr);
    PRECONDITION(pos!=neighMap.end(),"bad argument");
    done[curr] = 1;
    res.push_back(curr);

    const INT_VECT &neighs = pos->second;
#if 0
    std::cerr<<"depth: "<<depth<<" ring: "<<curr<<" size: "<<res.size()<<" neighs: "<<neighs.size()<<std::endl;
    std::cerr<<"   ";
    std::copy(neighs.begin(),neighs.end(),std::ostream_iterator<int>(std::cerr," "));
    std::cerr<<"\n";
#endif
    for (INT_VECT_CI ni = neighs.begin(); ni != neighs.end(); ++ni) {
      if (!done[*ni]) {
        pickFusedRings((*ni), neighMap, res, done,depth+1);
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
                           INT_INT_VECT_MAP &neighMap,
                           unsigned int maxSize) {
    int nrings = brings.size();
    int i, j;
    INT_VECT ring1;
    for (i = 0; i < nrings; i++) {
      INT_VECT neighs;
      neighMap[i] = neighs;
    }

    for (i = 0; i < nrings; i++) {
      if(maxSize && brings[i].size()>maxSize) continue;
      ring1 = brings[i];
      for (j = i + 1; j < nrings; j++) {
        if(maxSize && brings[j].size()>maxSize) continue;
        INT_VECT inter;
        Intersect(ring1, brings[j], inter);
        if (inter.size() > 0) {
          neighMap[i].push_back(j);
          neighMap[j].push_back(i);
        }
      }
    }
#if 0
    for (i = 0; i < nrings; i++) {
      std::cerr<<"**************\n    "<<i<<"\n*************\n";
      std::copy(neighMap[i].begin(),neighMap[i].end(),std::ostream_iterator<int>(std::cerr," "));
      std::cerr<<"\n";
    }
#endif
    
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
                                 int &narom, // number of aromatic ring so far
                                 unsigned int maxNumFusedRings
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
    // now mark bonds that have a count of 1 to be aromatic;
    for (bci = bndCntr.begin(); bci != bndCntr.end(); bci++) {
      if ((*bci).second == 1) {
        Bond *bond=mol.getBondWithIdx(bci->first);
        bond->setIsAromatic(true);
        switch(bond->getBondType()){
        case Bond::SINGLE:
        case Bond::DOUBLE:
          bond->setBondType(Bond::AROMATIC);
          break;
        default:
          break;
        }
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
    const ROMol &mol=at->getOwningMol();
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = mol.getAtomBonds(at);
    while(beg!=end){
      if(!mol.getRingInfo()->numBondRings(mol[*beg]->getIdx())){
        if (mol[*beg]->getValenceContrib(at) >= 2.0) {
          who = mol[*beg]->getOtherAtomIdx(at->getIdx());
          return true;
        }
       }
      ++beg;
    }
    return false;
  }
   
  bool incidentCyclicMultipleBond(const Atom *at) {
    PRECONDITION(at,"bad atom");
    const ROMol &mol=at->getOwningMol();
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = mol.getAtomBonds(at);
    while(beg!=end){
      if(mol.getRingInfo()->numBondRings(mol[*beg]->getIdx())){
        if (mol[*beg]->getValenceContrib(at) >= 2.0) {
          return true;
         }
       }
       beg++;
    }
    return false;
  }

  bool incidentMultipleBond(const Atom *at) {
    PRECONDITION(at,"bad atom");
    int deg=at->getDegree()+at->getNumExplicitHs();
    ROMol::OEDGE_ITER beg,end;
    boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
    while(beg!=end){
      BOND_SPTR bond=at->getOwningMol()[*beg];
      if(bond->getBondType()==Bond::ZERO) --deg;
      ++beg;
    }
    return at->getExplicitValence()!=static_cast<int>(deg);
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
    } else if (rup==2) {
      aromatic = true;
    }
    return aromatic;
  }
    

  
  void applyHuckelToFused(ROMol &mol, // molecule of interets
                          const VECT_INT_VECT &srings, // list of all ring as atom IDS
                          const VECT_INT_VECT &brings, // list of all rings as bond ids
                          const INT_VECT &fused, // list of ring ids in the current fused system
                          const VECT_EDON_TYPE &edon, // eletron donor state for each atom
                          INT_INT_VECT_MAP &ringNeighs, // list of neighbors for eac candidate ring
                          int &narom, // number of aromatic ring so far
                          unsigned int maxNumFusedRings
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
        
        curSize++;
        // check is we are done with all the atoms in the fused
        // system, if so quit. This is a fix for Issue252 REVIEW: is
        // this check sufficient or should we add an additional
        // contraint on the the number of combinations of rings in a
        // fused system that we will try. The number of combinations
        // can obviously be quite large when the number of rings in
        // the fused system is large
        if ((doneAtoms.size() == nAtms) || (curSize > std::min(nrings,maxNumFusedRings)) ){
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

  bool isAtomCandForArom(const Atom *at, const ElectronDonorType edon) {
    PRECONDITION(at,"bad atom");
    // limit aromaticity to:
    //   - the first two rows of the periodic table
    //   - Se and Te
    if (at->getAtomicNum() > 18 &&
        at->getAtomicNum()!=34 &&
        at->getAtomicNum()!=52 ) {
      return false;
    }
    switch (edon) {
    case VacantElectronDonorType:
    case OneElectronDonorType:
    case TwoElectronDonorType:
    case OneOrTwoElectronDonorType:
    case AnyElectronDonorType:
      break;
    default:
      return(false);      
    }


    // atoms that aren't in their default valence state also get shut out
    int defVal=PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
    if(defVal>0 &&
       at->getTotalValence()>(PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum()-
                                                                           at->getFormalCharge()))){
      return false;
    }

    // We are going to explicitly disallow atoms that have more
    // than one double or triple bond. This is to handle 
    // the situation:
    //   C1=C=NC=N1 (sf.net bug 1934360) 
    int nUnsaturations=at->getExplicitValence()-at->getDegree();
    if(nUnsaturations>1){
      unsigned int nMult=0;
      const ROMol &mol = at->getOwningMol();
      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = mol.getAtomBonds(at);
      while(beg!=end){
        switch(mol[*beg]->getBondType()){
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
        // This is borderline:  no electron to spare but may have an empty p-orbital 
        // Not sure if this should return vacantElectronDonorType 
        // FIX: explicitly doing this as a note for potential problems
        // 
        res = VacantElectronDonorType;
      }
      else if (incidentCyclicMultipleBond(at)) {
        // no electron but has one in a in cycle multiple bond
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

      // default valence :
      int dv=PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum()); 
      if (dv <=1) {
        // univalent elements can't be either aromatic or conjugated
        return -1;
      }

      // total atom degree:
      int degree=at->getDegree() + at->getTotalNumHs(); 

      ROMol::OEDGE_ITER beg,end;
      boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
      while(beg!=end){
        BOND_SPTR bond=at->getOwningMol()[*beg];
        if(bond->getBondType()==Bond::ZERO) --degree;
        ++beg;
      }

      // if we are more than 3 coordinated we should not be aromatic
      if (degree > 3) {
        return -1;
      }

      // number of lone pair electrons = (outer shell elecs) - (default valence)
      int nlp; 
      nlp = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum()) - dv; 

      // subtract the charge to get the true number of lone pair electrons:
      nlp -= at->getFormalCharge(); 

      int nRadicals=at->getNumRadicalElectrons();
      
      // num electrons available for donation into the pi system:
      int res = (dv-degree) + nlp - nRadicals;

      if(res>1){
        // if we have an incident bond with order higher than 2,
        // (e.g. triple or higher), we only want to return 1 electron
        // we detect this using the total unsaturation, because we
        // know that there aren't multiple unsaturations (detected
        // above in isAtomCandForArom())
        int nUnsaturations=at->getExplicitValence()-at->getDegree();
        if(nUnsaturations>1){
          res=1;
        }
      }

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
      if(mol.getRingInfo()->isInitialized()){
        srings = mol.getRingInfo()->atomRings();
      } else {
        MolOps::symmetrizeSSSR(mol, srings);
      }

      int narom = 0;
      // loop over all the atoms in the rings that can be candidates
      // for aromaticity
      // Atoms are candidates if 
      //   - it is part of ring
      //   - has one or more electron to donate or has empty p-orbitals
      int natoms = mol.getNumAtoms();
      boost::dynamic_bitset<> acands(natoms);
      boost::dynamic_bitset<> aseen(natoms);
      VECT_EDON_TYPE edon(natoms);

      VECT_INT_VECT cRings; // holder for rings that are candidates for aromaticity
      for (VECT_INT_VECT_I vivi = srings.begin();
           vivi != srings.end(); ++vivi) {
        bool allAromatic=true;
        for (INT_VECT_I ivi = (*vivi).begin();
             ivi!=(*vivi).end(); ++ivi) {
          if(aseen[*ivi]){
            if(!acands[*ivi]) allAromatic=false;
            continue;
          }
          aseen[*ivi]=1;
          Atom *at = mol.getAtomWithIdx(*ivi);
        
          // now that the atom is part of ring check if it can donate
          // electron or has empty orbitals. Record the donor type
          // information in 'edon' - we will need it when we get to
          // the Huckel rule later
          edon[*ivi] = getAtomDonorTypeArom(at);
          acands[*ivi]=isAtomCandForArom(at, edon[*ivi]);
          if(!acands[*ivi]) allAromatic=false;
        }
        if(allAromatic){
          cRings.push_back((*vivi));
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
      RingUtils::makeRingNeighborMap(brings, neighMap, maxFusedAromaticRingSize);
      
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
        applyHuckelToFused(mol, cRings, brings, fused, edon, neighMap, narom,6);

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
