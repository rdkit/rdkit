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

      // ------------------------------------------------------------------------
      //
      // the matrix returned by this contains:
      //  -1: if atoms i and j are directly connected
      // idx: if atoms i and j are connected via atom idx
      //  -2: otherwise
      //
      //  NOTE: the caller is responsible for calling delete []
      //  on the result
      //
      // ------------------------------------------------------------------------
      boost::shared_array<int> buildNeighborMatrix(const ROMol &mol){
        unsigned int nAtoms = mol.getNumAtoms();
        boost::shared_array<int> res(new int[nAtoms*nAtoms]);
        for(unsigned int i=0;i<nAtoms;i++){
          unsigned int iTab=i*nAtoms;
          for(unsigned int j=i;j<nAtoms;j++){
            res[iTab+j] = -2;
            res[i+j*nAtoms] = -2;
          }
        }
        for(unsigned int i=0;i<mol.getNumBonds();i++){
          const Bond *bondi=mol.getBondWithIdx(i);

          res[bondi->getBeginAtomIdx()*nAtoms+bondi->getEndAtomIdx()] = -1;
          res[bondi->getEndAtomIdx()*nAtoms+bondi->getBeginAtomIdx()] = -1;

          for(unsigned int j=i+1;j<mol.getNumBonds();j++){
            const Bond *bondj=mol.getBondWithIdx(j);
            int idx1=-1,idx2=-1,idx3=-1;
            if(bondi->getBeginAtomIdx()==bondj->getBeginAtomIdx()){
              idx1 = bondi->getEndAtomIdx();
              idx2 = bondi->getBeginAtomIdx();
              idx3 = bondj->getEndAtomIdx();
            } else if(bondi->getBeginAtomIdx()==bondj->getEndAtomIdx()){
              idx1 = bondi->getEndAtomIdx();
              idx2 = bondi->getBeginAtomIdx();
              idx3 = bondj->getBeginAtomIdx();
            } else if(bondi->getEndAtomIdx()==bondj->getBeginAtomIdx()){
              idx1 = bondi->getBeginAtomIdx();
              idx2 = bondi->getEndAtomIdx();
              idx3 = bondj->getEndAtomIdx();
            } else if(bondi->getEndAtomIdx()==bondj->getEndAtomIdx()){
              idx1 = bondi->getBeginAtomIdx();
              idx2 = bondi->getEndAtomIdx();
              idx3 = bondj->getBeginAtomIdx();
            }
            if(idx1>-1){
              res[idx1*nAtoms+idx3] = idx2;
              res[idx3*nAtoms+idx1] = idx2;
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
                     ForceFields::ForceField *field,boost::shared_array<int> neighborMatrix){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        unsigned int nAtoms=mol.getNumAtoms();
        for(unsigned int i=0;i<nAtoms;i++){
          if(!params[i]) continue;
          for(unsigned int j=i+1;j<nAtoms;j++){
            if(!params[j]) continue;
            if(neighborMatrix[i*nAtoms+j]>-1){
              int k = neighborMatrix[i*nAtoms+j];
              if(!params[k]) continue;
              const Atom *atomK = mol.getAtomWithIdx(k);
              // skip special cases:
              if( !(atomK->getHybridization()==Atom::SP3D && atomK->getDegree()==5) ){
                const Bond *b1 =mol.getBondBetweenAtoms(i,k);
                const Bond *b2 =mol.getBondBetweenAtoms(j,k);
                // FIX: recognize amide bonds here.
                AngleBendContrib *contrib;
                int order=0;
                switch(atomK->getHybridization()){
                case Atom::SP:
                  order=2;
                  break;
                case Atom::SP3D2:
                  order=4;
                  break;
                default:
                  order=0;
                  break;
                } 
                  
                contrib = new AngleBendContrib(field,i,k,j,
                                               b1->getBondTypeAsDouble(),
                                               b2->getBondTypeAsDouble(),
                                               params[i],params[k],params[j],order);
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
                        ForceFields::ForceField *field,boost::shared_array<int> neighborMatrix,
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
            if(neighborMatrix[i*nAtoms+j]==-2){
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
                      if(okToIncludeTorsion(mol,bond,bIdx,idx1,idx2,eIdx)){
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
                      }
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

#if 0
      // ------------------------------------------------------------------------
      //
      //
      //
      // ------------------------------------------------------------------------
      void addInversions(const ROMol &mol,const AtomicParamVect &params,
                       ForceFields::ForceField *field){
        PRECONDITION(mol.getNumAtoms()==params.size(),"bad parameters");
        PRECONDITION(field,"bad forcefield");

        unsigned int nAtoms=mol.getNumAtoms();
      }
#endif
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
      boost::shared_array<int> neighborMat = Tools::buildNeighborMatrix(mol);
      Tools::addAngles(mol,params,res,neighborMat);
      Tools::addAngleSpecialCases(mol,confId,params,res);
      Tools::addNonbonded(mol,confId,params,res,neighborMat,vdwThresh,ignoreInterfragInteractions);
      Tools::addTorsions(mol,params,res);
      //Tools::addInversions(mol,params,res);

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
