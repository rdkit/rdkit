// $Id$
//
//  Copyright (C) 2004-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <cmath>

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <ForceField/ForceField.h>
#include <ForceField/UFF/Params.h>
#include <ForceField/UFF/Contribs.h>

#include "AtomTyper.h"
#include "Builder.h"
namespace RDKit {
  namespace UFF {
    using namespace ForceFields::UFF;

    namespace Tools {
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addBonds(const ROMol &mol,const AtomicParamVect &params,
                    ForceFields::ForceField *field){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        for (ROMol::ConstBondIterator bi=mol.beginBonds();
             bi != mol.endBonds();
             bi++) {
          int idx1=(*bi)->getBeginAtomIdx();
          int idx2=(*bi)->getEndAtomIdx();

          // FIX: recognize amide bonds here.

          if(params[idx1]&&params[idx2]){
            BondStretchContrib *contrib;
            contrib = new BondStretchContrib(field,idx1,idx2,
                                             (*bi)->getBondTypeAsDouble(),
                                             params[idx1],params[idx2]);
            field->contribs().push_back(ForceFields::ContribPtr(contrib));
          }
        }
      }

      void setTwoBitCell(boost::shared_array<boost::uint8_t> &res,
        unsigned int pos, boost::uint8_t value)
      {
        unsigned int twoBitPos = pos / 4;
        unsigned int shift = 2 * (pos % 4);
        boost::uint8_t twoBitMask = 3 << shift;
        res[twoBitPos] = ((res[twoBitPos] & (~twoBitMask)) | (value << shift));
      }
      
      
      boost::uint8_t getTwoBitCell
        (boost::shared_array<boost::uint8_t> &res, unsigned int pos)
      {
        unsigned int twoBitPos = pos / 4;
        unsigned int shift = 2 * (pos % 4);
        boost::uint8_t twoBitMask = 3 << shift;
        
        return ((res[twoBitPos] & twoBitMask) >> shift);
      }
      

      // ------------------------------------------------------------------------
      //
      // the two-bit matrix returned by this contains:
      //   0: if atoms i and j are directly connected
      //   1: if atoms i and j are connected via an atom
      //   3: otherwise
      //
      //  NOTE: the caller is responsible for calling delete []
      //  on the result
      //
      // ------------------------------------------------------------------------
      boost::shared_array<boost::uint8_t> buildNeighborMatrix(const ROMol &mol)
      {
        unsigned int nAtoms = mol.getNumAtoms();
        unsigned nTwoBitCells = (nAtoms * nAtoms - 1) / 4 + 1;
        boost::shared_array<boost::uint8_t> res(new boost::uint8_t[nTwoBitCells]);
        for (unsigned int i = 0; i < nTwoBitCells; ++i) {
          res[i] = 0;
        }
        for (unsigned int i = 0; i < nAtoms; ++i) {
          unsigned int iTab = i * nAtoms;
          for (unsigned int j = i; j < nAtoms; ++j) {
            setTwoBitCell(res, iTab + j, RELATION_1_X);
            setTwoBitCell(res, i + j * nAtoms, RELATION_1_X);
          }
        }
        for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
          const Bond *bondi = mol.getBondWithIdx(i);

          setTwoBitCell(res, bondi->getBeginAtomIdx() * nAtoms + bondi->getEndAtomIdx(), RELATION_1_2);
          setTwoBitCell(res, bondi->getEndAtomIdx() * nAtoms + bondi->getBeginAtomIdx(), RELATION_1_2);

          for (unsigned int j = i + 1; j < mol.getNumBonds(); ++j) {
            const Bond *bondj = mol.getBondWithIdx(j);
            int idx1 = -1;
            int idx3 = -1;
            if (bondi->getBeginAtomIdx() == bondj->getBeginAtomIdx()) {
              idx1 = bondi->getEndAtomIdx();
              idx3 = bondj->getEndAtomIdx();
            }
            else if (bondi->getBeginAtomIdx() == bondj->getEndAtomIdx()) {
              idx1 = bondi->getEndAtomIdx();
              idx3 = bondj->getBeginAtomIdx();
            }
            else if (bondi->getEndAtomIdx() == bondj->getBeginAtomIdx()) {
              idx1 = bondi->getBeginAtomIdx();
              idx3 = bondj->getEndAtomIdx();
            }
            else if (bondi->getEndAtomIdx() == bondj->getEndAtomIdx()) {
              idx1 = bondi->getBeginAtomIdx();
              idx3 = bondj->getBeginAtomIdx();
            }
            if (idx1 > -1) {
              setTwoBitCell(res, idx1 * nAtoms + idx3, RELATION_1_3);
              setTwoBitCell(res, idx3 * nAtoms + idx1, RELATION_1_3);
            }
          }
        }
        return res;
      }
      
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addAngles(const ROMol &mol,const AtomicParamVect &params,
                     ForceFields::ForceField *field){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");
        ROMol::ADJ_ITER nbr1Idx;
        ROMol::ADJ_ITER end1Nbrs;
        ROMol::ADJ_ITER nbr2Idx;
        ROMol::ADJ_ITER end2Nbrs;
        RingInfo *rings=mol.getRingInfo();

        unsigned int nAtoms=mol.getNumAtoms();
        for(unsigned int j=0;j<nAtoms;j++){
          if(!params[j]) continue;
          const Atom *atomJ=mol.getAtomWithIdx(j);
          if(atomJ->getDegree()==1) continue;
          boost::tie(nbr1Idx,end1Nbrs)=mol.getAtomNeighbors(atomJ);
          for (;nbr1Idx!=end1Nbrs;nbr1Idx++) {
            const Atom *atomI=mol[*nbr1Idx].get();
            unsigned int i=atomI->getIdx();
            if(!params[i]) continue;
            boost::tie(nbr2Idx,end2Nbrs)=mol.getAtomNeighbors(atomJ);
            for (;nbr2Idx!=end2Nbrs;nbr2Idx++) {
              if (nbr2Idx<(nbr1Idx+1)) {
                continue;
              }
              const Atom *atomK=mol[*nbr2Idx].get();
              unsigned int k=atomK->getIdx();
              if(!params[k]) continue;
              // skip special cases:
              if( !(atomJ->getHybridization()==Atom::SP3D && atomJ->getDegree()==5) ){
                const Bond *b1 =mol.getBondBetweenAtoms(i,j);
                const Bond *b2 =mol.getBondBetweenAtoms(k,j);
                // FIX: recognize amide bonds here.
                AngleBendContrib *contrib;
                int order=0;
                switch(atomJ->getHybridization()){
                case Atom::SP:
                  order=1;
                  break;
                case Atom::SP2:
                  order=3;
                  // the following is a hack to get decent geometries
                  // with 3- and 4-membered rings incorporating sp2 atoms
                  // if the central atom is in a ring of size 3
                  if (rings->isAtomInRingOfSize(j, 3)) {
                    // if the central atom and one of the bonded atoms, but not the
                    //  other one are inside a ring, then this angle is between a
                    // ring substituent and a ring edge
                    if ((rings->isAtomInRingOfSize(i, 3) && !rings->isAtomInRingOfSize(k, 3))
                      || (!rings->isAtomInRingOfSize(i, 3) && rings->isAtomInRingOfSize(k, 3))) {
                      order = 30;
                    }
                    // if all atoms are inside the ring, then this is one of ring angles
                    else if (rings->isAtomInRingOfSize(i, 3) && rings->isAtomInRingOfSize(k, 3)) {
                      order = 35;
                    }
                  }
                  // if the central atom is in a ring of size 4
                  else if (rings->isAtomInRingOfSize(j, 4)) {
                    // if the central atom and one of the bonded atoms, but not the
                    //  other one are inside a ring, then this angle is between a
                    // ring substituent and a ring edge
                    if ((rings->isAtomInRingOfSize(i, 4) && !rings->isAtomInRingOfSize(k, 4))
                      || (!rings->isAtomInRingOfSize(i, 4) && rings->isAtomInRingOfSize(k, 4))) {
                      order = 40;
                    }
                    // if all atoms are inside the ring, then this is one of ring angles
                    else if (rings->isAtomInRingOfSize(i, 4) && rings->isAtomInRingOfSize(k, 4)) {
                      order = 45;
                    }
                  }
                  // end of the hack
                  break;
                case Atom::SP3D2:
                  order=4;
                  break;
                default:
                  order=0;
                  break;
                } 
                  
                contrib = new AngleBendContrib(field,i,j,k,
                                               b1->getBondTypeAsDouble(),
                                               b2->getBondTypeAsDouble(),
                                               params[i],params[j],params[k],order);
                field->contribs().push_back(ForceFields::ContribPtr(contrib));
              }
            }
          }
        }
      }

      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addTrigonalBipyramidAngles(const Atom *atom,const ROMol &mol, int confId,
                                      const AtomicParamVect &params,
                                      ForceFields::ForceField *field){
        PRECONDITION(atom,"bad atom");
        PRECONDITION(atom->getHybridization()==Atom::SP3D,"bad hybridization");
        PRECONDITION(atom->getDegree()==5,"bad degree");
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        const Bond *ax1=0,*ax2=0;
        const Bond *eq1=0,*eq2=0,*eq3=0;

        const Conformer &conf = mol.getConformer(confId);
        //------------------------------------------------------------
        // identify the axial and equatorial bonds:
        double mostNeg=100.0;
        ROMol::OEDGE_ITER beg1,end1;
        boost::tie(beg1,end1) = mol.getAtomBonds(atom);
        unsigned int aid = atom->getIdx();
        while(beg1!=end1){
          const Bond *bond1=mol[*beg1].get();
          unsigned int oaid = bond1->getOtherAtomIdx(aid);
          RDGeom::Point3D v1=conf.getAtomPos(aid).directionVector(conf.getAtomPos(oaid));
                  
          ROMol::OEDGE_ITER beg2,end2;
          boost::tie(beg2,end2) = mol.getAtomBonds(atom);
          while(beg2 != end2){
            const Bond *bond2=mol[*beg2].get();
            if(bond2->getIdx() > bond1->getIdx()){
              unsigned int oaid2 = bond2->getOtherAtomIdx(aid);
              RDGeom::Point3D v2=conf.getAtomPos(aid).directionVector(conf.getAtomPos(oaid2));
              double dot=v1.dotProduct(v2);
              if(dot<mostNeg){
                mostNeg = dot;
                ax1 = bond1;
                ax2 = bond2;
              }
            }
            ++beg2;
          }
          ++beg1;
        }
        CHECK_INVARIANT(ax1,"axial bond not found");
        CHECK_INVARIANT(ax2,"axial bond not found");
	
        boost::tie(beg1,end1) = mol.getAtomBonds(atom);
        while(beg1!=end1){
          const Bond *bond=mol[*beg1].get();
	  ++beg1;
	  if(bond==ax1 || bond==ax2) continue;
	  if(!eq1) eq1=bond;
	  else if(!eq2) eq2=bond;
	  else if(!eq3) eq3=bond;
	}

        CHECK_INVARIANT(eq1,"equatorial bond not found");
        CHECK_INVARIANT(eq2,"equatorial bond not found");
        CHECK_INVARIANT(eq3,"equatorial bond not found");

        
        //------------------------------------------------------------
        // alright, add the angles:
        AngleBendContrib *contrib;
        int atomIdx=atom->getIdx();
        int i,j;

        // Axial-Axial
        i=ax1->getOtherAtomIdx(atomIdx);
        j=ax2->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax1->getBondTypeAsDouble(),
                                         ax2->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j],2);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }        
        // Equatorial-Equatorial
        i=eq1->getOtherAtomIdx(atomIdx);
        j=eq2->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         eq1->getBondTypeAsDouble(),
                                         eq2->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j],3);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=eq1->getOtherAtomIdx(atomIdx);
        j=eq3->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         eq1->getBondTypeAsDouble(),
                                         eq3->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j],3);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=eq2->getOtherAtomIdx(atomIdx);
        j=eq3->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         eq2->getBondTypeAsDouble(),
                                         eq3->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j],3);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }        

        // Axial-Equatorial
        i=ax1->getOtherAtomIdx(atomIdx);
        j=eq1->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax1->getBondTypeAsDouble(),
                                         eq1->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j]);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=ax1->getOtherAtomIdx(atomIdx);
        j=eq2->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax1->getBondTypeAsDouble(),
                                         eq2->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j]);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=ax1->getOtherAtomIdx(atomIdx);
        j=eq3->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax1->getBondTypeAsDouble(),
                                         eq3->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j]);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=ax2->getOtherAtomIdx(atomIdx);
        j=eq1->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax2->getBondTypeAsDouble(),
                                         eq1->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j]);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=ax2->getOtherAtomIdx(atomIdx);
        j=eq2->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax2->getBondTypeAsDouble(),
                                         eq2->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j]);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
        i=ax2->getOtherAtomIdx(atomIdx);
        j=eq3->getOtherAtomIdx(atomIdx);
        if(params[i]&&params[j]){
          contrib = new AngleBendContrib(field,i,atomIdx,j,
                                         ax2->getBondTypeAsDouble(),
                                         eq3->getBondTypeAsDouble(),
                                         params[i],params[atomIdx],params[j]);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
        }
      }

      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addAngleSpecialCases(const ROMol &mol, int confId, const AtomicParamVect &params,
                                ForceFields::ForceField *field){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        unsigned int nAtoms=mol.getNumAtoms();
        for(unsigned int i=0;i<nAtoms;i++){
          const Atom *atom = mol.getAtomWithIdx(i);
          // trigonal bipyramidal:
          if( (atom->getHybridization()==Atom::SP3D && atom->getDegree()==5) ){
            addTrigonalBipyramidAngles(atom,mol,confId, params,field);
          }
        }
      }
          
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addNonbonded(const ROMol &mol,int confId,const AtomicParamVect &params,
                        ForceFields::ForceField *field,boost::shared_array<boost::uint8_t> neighborMatrix,
                        double vdwThresh,bool ignoreInterfragInteractions){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        INT_VECT fragMapping;
        if(ignoreInterfragInteractions){
          std::vector<ROMOL_SPTR> molFrags=MolOps::getMolFrags(mol,true,&fragMapping);
        }

        unsigned int nAtoms=mol.getNumAtoms();
        const Conformer &conf = mol.getConformer(confId);
        for(unsigned int i=0;i<nAtoms;i++){
          if(!params[i]) continue;
          for(unsigned int j=i+1;j<nAtoms;j++){
            if(!params[j] || (ignoreInterfragInteractions && fragMapping[i]!=fragMapping[j])){
              continue;
            }
            if(getTwoBitCell(neighborMatrix,i*nAtoms+j)>=RELATION_1_4){
              double dist=(conf.getAtomPos(i) - conf.getAtomPos(j)).length();
              if(dist <
                 vdwThresh*Utils::calcNonbondedMinimum(params[i],params[j])){
                vdWContrib *contrib;
                contrib = new vdWContrib(field,i,j,params[i],params[j]);
                field->contribs().push_back(ForceFields::ContribPtr(contrib));
              }
            }
          }
        }
      }

      #if 0
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      bool okToIncludeTorsion(const ROMol &mol,const Bond *bond,
                              int idx1,int idx2,int idx3,int idx4){
        bool res=true;
        RingInfo *rings=mol.getRingInfo();
        // having torsions in small rings makes the solver unstable
        // and tends to yield poor-quality geometries, so filter those out:
        if(rings->isBondInRingOfSize(bond->getIdx(),3)){
          res = false;
        }// else if(rings->isBondInRingOfSize(bond->getIdx(),4)){
         // res = false;
        //}
        return res;
      }
      #endif
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addTorsions(const ROMol &mol,const AtomicParamVect &params,
                       ForceFields::ForceField *field,
                       std::string torsionBondSmarts){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        // find all of the torsion bonds:
        std::vector<MatchVectType> matchVect;
        ROMol *query=SmartsToMol(torsionBondSmarts);
        TEST_ASSERT(query);
        unsigned int nHits=SubstructMatch(mol,*query,matchVect);
        delete query;

        for(unsigned int i=0; i<nHits; i++){
          MatchVectType match=matchVect[i];
          TEST_ASSERT(match.size()==2);
          int idx1=match[0].second;
          int idx2=match[1].second;
          if(!params[idx1]||!params[idx2]) continue;
          const Bond *bond=mol.getBondBetweenAtoms(idx1,idx2);
          std::vector<TorsionAngleContrib *> contribsHere;
          TEST_ASSERT(bond);
          const Atom *atom1=mol.getAtomWithIdx(idx1);
          const Atom *atom2=mol.getAtomWithIdx(idx2);

          if( (atom1->getHybridization()==Atom::SP2||atom1->getHybridization()==Atom::SP3) &&
              (atom2->getHybridization()==Atom::SP2||atom2->getHybridization()==Atom::SP3) ){
            ROMol::OEDGE_ITER beg1,end1;
            boost::tie(beg1,end1) = mol.getAtomBonds(atom1);
            while(beg1!=end1){
              const Bond *tBond1=mol[*beg1].get();
              if(tBond1!=bond){
                int bIdx = tBond1->getOtherAtomIdx(idx1);
                ROMol::OEDGE_ITER beg2,end2;
                boost::tie(beg2,end2) = mol.getAtomBonds(atom2);
                while(beg2 != end2){
                  const Bond *tBond2=mol[*beg2].get();
                  if(tBond2!=bond && tBond2!=tBond1){
                    int eIdx=tBond2->getOtherAtomIdx(idx2);
                    // make sure this isn't a three-membered ring:
                    if(eIdx != bIdx){
                      // we now have a torsion involving atoms (bonds):
                      //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
                      TorsionAngleContrib *contrib;

                      // if either of the end atoms is SP2 hybridized, set a flag
                      // here.  
                      bool hasSP2=false;
                      if(mol.getAtomWithIdx(bIdx)->getHybridization()==Atom::SP2 ||
                         mol.getAtomWithIdx(bIdx)->getHybridization()==Atom::SP2) {
                        hasSP2 = true;
                      }
                      //std::cout << "Torsion: " << bIdx << "-" << idx1 << "-" << idx2 << "-" << eIdx << std::endl;
                      //if(okToIncludeTorsion(mol,bond,bIdx,idx1,idx2,eIdx)){
                        //std::cout << "  INCLUDED" << std::endl;
                        contrib = new TorsionAngleContrib(field,bIdx,idx1,idx2,eIdx,
                                                          bond->getBondTypeAsDouble(),
                                                          atom1->getAtomicNum(),
                                                          atom2->getAtomicNum(),
                                                          atom1->getHybridization(),
                                                          atom2->getHybridization(),
                                                          params[idx1],params[idx2],
                                                          hasSP2);
                        field->contribs().push_back(ForceFields::ContribPtr(contrib));
                        contribsHere.push_back(contrib);
                      //}
                    }
                  }
                  beg2++;
                }
              }
              beg1++;
            }
          }
          // now divide the force constant for each contribution to the torsion energy
          // about this bond by the number of contribs about this bond:
          for(std::vector<TorsionAngleContrib *>::iterator chI=contribsHere.begin();
              chI!=contribsHere.end();++chI){
            (*chI)->scaleForceConstant(contribsHere.size());
          }
        }

      }

      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addInversions(const ROMol &mol,const AtomicParamVect &params,
                       ForceFields::ForceField *field){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        unsigned int idx[4];
        unsigned int n[4];
        const Atom *atom[4];
        ROMol::ADJ_ITER nbrIdx;
        ROMol::ADJ_ITER endNbrs;

        for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
          atom[1] = mol.getAtomWithIdx(idx[1]);
          int at2AtomicNum = atom[1]->getAtomicNum();
          // if the central atom is not carbon, nitrogen, oxygen,
          // phosphorous, arsenic, antimonium or bismuth, skip it
          if (((at2AtomicNum != 6) && (at2AtomicNum != 7) && (at2AtomicNum != 8)
            && (at2AtomicNum != 15) && (at2AtomicNum != 33) && (at2AtomicNum != 51)
            && (at2AtomicNum != 83)) || (atom[1]->getDegree() != 3)) {
            continue;
          }
          // if the central atom is carbon, nitrogen or oxygen
          // but hybridization is not sp2, skip it
          if (((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8))
            && (atom[1]->getHybridization() != Atom::SP2)) {
            continue;
          }
          boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
          unsigned int i = 0;
          bool isBoundToSP2O = false;
          for (; nbrIdx != endNbrs; ++nbrIdx) {
            atom[i] = mol[*nbrIdx].get();
            idx[i] = atom[i]->getIdx();
            // if the central atom is sp2 carbon and is
            // bound to sp2 oxygen, set a flag
            if (!isBoundToSP2O) {
              isBoundToSP2O = ((at2AtomicNum == 6) && (atom[i]->getAtomicNum() == 8)
                && (atom[i]->getHybridization() == Atom::SP2));
            }
            if (!i) {
              ++i;
            }
            ++i;
          }
          for (unsigned int i = 0; i < 3; ++i) {
            n[1] = 1;
            switch (i) {
              case 0:
                n[0] = 0;
                n[2] = 2;
                n[3] = 3;
              break;
              
              case 1:
                n[0] = 0;
                n[2] = 3;
                n[3] = 2;
              break;
              
              case 2:
                n[0] = 2;
                n[2] = 3;
                n[3] = 0;
              break;
            }
            InversionContrib *contrib;
            contrib = new InversionContrib(field, idx[n[0]], idx[n[1]], idx[n[2]],
                                             idx[n[3]], at2AtomicNum, isBoundToSP2O);
            field->contribs().push_back(ForceFields::ContribPtr(contrib));
          }
        }
      }
    } // end of namespace Tools
    
    // ------------------------------------------------------------------------
    //
    //
    //
    // ------------------------------------------------------------------------
    ForceFields::ForceField *constructForceField(ROMol &mol,
                                                 const AtomicParamVect &params,
                                                 double vdwThresh, int confId,
                                                 bool ignoreInterfragInteractions){
      PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        
      ForceFields::ForceField *res=new ForceFields::ForceField();

      // add the atomic positions:
      Conformer &conf = mol.getConformer(confId);
      for(unsigned int i=0;i<mol.getNumAtoms();i++){
        res->positions().push_back(&conf.getAtomPos(i));
      }
      
      Tools::addBonds(mol,params,res);
      Tools::addAngles(mol,params,res);
      Tools::addAngleSpecialCases(mol,confId,params,res);
      boost::shared_array<boost::uint8_t> neighborMat = Tools::buildNeighborMatrix(mol);
      Tools::addNonbonded(mol,confId,params,res,neighborMat,vdwThresh,ignoreInterfragInteractions);
      Tools::addTorsions(mol,params,res);
      Tools::addInversions(mol,params,res);

      return res;
    }
    
    // ------------------------------------------------------------------------
    //
    //
    //
    // ------------------------------------------------------------------------
    ForceFields::ForceField *constructForceField(ROMol &mol,double vdwThresh, int confId,
                                                 bool ignoreInterfragInteractions){
      bool foundAll;
      AtomicParamVect params;
      boost::tie(params,foundAll)=getAtomTypes(mol);
      return constructForceField(mol,params,vdwThresh, confId,ignoreInterfragInteractions);
    }

  }
}
