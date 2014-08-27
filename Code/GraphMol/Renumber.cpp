//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <boost/foreach.hpp>


namespace RDKit {
  namespace MolOps {
    ROMol *renumberAtoms(const ROMol &mol,const std::vector<unsigned int> &newOrder){
      unsigned int nAts=mol.getNumAtoms();
      PRECONDITION(newOrder.size()==nAts,"bad newOrder size");

      std::vector<unsigned int> revOrder(nAts);
      for(unsigned int nIdx=0;nIdx<nAts;++nIdx){
        unsigned int oIdx=newOrder[nIdx];
        if(oIdx>nAts){
          throw ValueErrorException("idx value exceeds numAtoms");
        }
        revOrder[oIdx]=nIdx;
      }

      // ------
      // newOrder[i] : which atom should be in position i of the new mol
      // revOrder[i] : where atom i of the original mol landed in the new mol
      RWMol *res=new RWMol();

      // copy over the atoms:
      for(unsigned int nIdx=0;nIdx<nAts;++nIdx){
        unsigned int oIdx=newOrder[nIdx];
        const Atom *oAtom=mol.getAtomWithIdx(oIdx);
        Atom *nAtom=oAtom->copy();
        res->addAtom(nAtom,false,true);

        // take care of atom-numbering-dependent properties:
        if(nAtom->hasProp("_ringStereoAtoms")){
          // FIX: ought to be able to avoid this copy.
          INT_VECT nAtoms;
          nAtom->getProp<INT_VECT>("_ringStereoAtoms",nAtoms);
          BOOST_FOREACH(int &val,nAtoms){
            if(val<0){
              val=-1*(revOrder[(-val-1)]+1);
            } else {
              val=revOrder[val-1]+1;
            }
          }
          nAtom->setProp("_ringStereoAtoms",nAtoms,true);
        }
      }          
            
      // now the bonds:
      for(ROMol::ConstBondIterator bi=mol.beginBonds();
          bi!=mol.endBonds();++bi){
        const Bond *oBond=(*bi);
        Bond *nBond=oBond->copy();
        nBond->setBeginAtomIdx(revOrder[oBond->getBeginAtomIdx()]);
        nBond->setEndAtomIdx(revOrder[oBond->getEndAtomIdx()]);
        res->addBond(nBond,true);
        // take care of atom-numbering-dependent properties:
        BOOST_FOREACH(int &idx,nBond->getStereoAtoms()){
          idx=revOrder[idx];
        }
      }

      // Conformers:
      for(ROMol::ConstConformerIterator oConf=mol.beginConformers();
          oConf!=mol.endConformers();
          ++oConf){
        Conformer *nConf=new Conformer(nAts);
        for(unsigned int i=0;i<nAts;++i){
          nConf->setAtomPos(i,(*oConf)->getAtomPos(revOrder[i]));
        }
        res->addConformer(nConf);
      }

      // update the ring info:
      const RingInfo *oRings=mol.getRingInfo();
      if(oRings){
        RingInfo *nRings=res->getRingInfo();
        nRings->reset();
        nRings->initialize();
        for(unsigned int i=0;i<oRings->numRings();++i){
          const INT_VECT &oRing=oRings->atomRings()[i];
          INT_VECT nRing(oRing.size());
          for(unsigned int j=0;j<oRing.size();++j){
            nRing[j]=revOrder[oRing[j]];
          }
          nRings->addRing(nRing,oRings->bondRings()[i]);
        }
      }

      return dynamic_cast<ROMol *>(res);
    }
    
  }; // end of namespace MolOps
}; // end of namespace RDKit
