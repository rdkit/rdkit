// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <GraphMol/RDKitBase.h>
#include <ForceField/UFF/Params.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include "AtomTyper.h"

namespace RDKit {
  namespace UFF {
    using namespace ForceFields::UFF;

    namespace Tools {
      // ---------------------------------------------------------------
      void addAtomChargeFlags(const Atom *atom,std::string &atomKey,
                              bool tolerateChargeMismatch){
        PRECONDITION(atom,"bad atom");
        int totalValence=atom->getTotalValence();
          atom->getFormalCharge();

        // FIX: come up with some way of handling metals here
        switch(atom->getAtomicNum()){
        case 12: // Mg
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 13: // Al
          if(totalValence!=3){
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 14: // Si
          if(totalValence!=4){
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 15: // P
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          case 5:
            atomKey += "+5";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+5";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 16: // S
          if(atom->getHybridization() != Atom::SP2){
            switch(totalValence){
            case 2:
              atomKey += "+2";
              break;
            case 4:
              atomKey += "+4";
              break;
            case 6:
              atomKey += "+6";
              break;
            default:
              if(tolerateChargeMismatch) atomKey += "+6";
              BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
            }
          }
          break;
        case 30: // Zn
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 31: // Ga
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 33: // As
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 34: // Se
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 48: // Cd
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 49: // In
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 51: // Sb
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 52: // Te
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 80: // Hg
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 81: // Tl
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 82: // Pb
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 83: // Bi
          switch(totalValence){
          case 3:
            atomKey += "+3";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+3";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        case 84: // Sb
          switch(totalValence){
          case 2:
            atomKey += "+2";
            break;
          default:
            if(tolerateChargeMismatch) atomKey += "+2";
            BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized charge state for atom: "<< atom->getIdx() << std::endl;
          }
          break;
        }
      }
      
      // ---------------------------------------------------------------
      std::string getAtomLabel(const Atom *atom){
        PRECONDITION(atom,"bad atom");
        int atNum = atom->getAtomicNum();
        std::string atomKey=atom->getSymbol();
        if(atomKey.size()==1) atomKey+='_';
        PeriodicTable *table = PeriodicTable::getTable();

        // FIX: handle main group/organometallic cases better:
        if(atNum){
          // do not do hybridization on alkali metals or halogens:
          if( table->getNouterElecs(atNum)!=1 &&
              table->getNouterElecs(atNum)!=7 ){
            switch(atom->getAtomicNum()){
            case 12:
            case 13:
            case 14:
            case 15:
            case 50:
            case 51:
            case 52:
            case 81:
            case 82:
            case 83:
            case 84:
              atomKey += '3';
              if(atom->getHybridization()!=Atom::SP3){
                BOOST_LOG(rdWarningLog) << "UFFTYPER: Warning: hybridization set to SP3 for atom "<< atom->getIdx() << std::endl;
              }
              break;
            case 80:
              atomKey += '1';
              if(atom->getHybridization()!=Atom::SP){
                BOOST_LOG(rdWarningLog) << "UFFTYPER: Warning: hybridization set to SP for atom "<< atom->getIdx() << std::endl;
              }
              break;
            default:
              switch(atom->getHybridization()){
              case Atom::S:
                // don't need to do anything here
                break; 
              case Atom::SP:
                atomKey+='1';
                break;

              case Atom::SP2:
                if((atom->getIsAromatic()||MolOps::atomHasConjugatedBond(atom))
                   && (atNum==6 || atNum==7 || atNum==8 || atNum==16) ){
                  atomKey+='R';
                } else {
                  atomKey+='2';
                }
                break;

              case Atom::SP3:
                atomKey+='3';
                break;

              case Atom::SP3D:
                atomKey+='5';
                break;
          
              case Atom::SP3D2:
                atomKey+='6';
                break;
          
              default:
                BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized hybridization for atom: "<< atom->getIdx() << std::endl;
              }
            }
          }
        }
        // special cases by element type:
        addAtomChargeFlags(atom,atomKey);
        return atomKey;
      }
    } // end of namespace Tools

    // ---------------------------------------------------------------
    std::pair<AtomicParamVect,bool>  getAtomTypes(const ROMol &mol,const std::string &paramData){
      bool foundAll=true;
      ParamCollection *params=ParamCollection::getParams();

      AtomicParamVect paramVect;
      paramVect.resize(mol.getNumAtoms());

      for(unsigned int i=0;i<mol.getNumAtoms();i++){
        const Atom *atom=mol.getAtomWithIdx(i);

        // construct the atom key:
        std::string atomKey = Tools::getAtomLabel(atom);

        // ok, we've got the atom key, now get the parameters:
        const AtomicParams *theParams=(*params)(atomKey);
        if(!theParams){
          foundAll=false;
          BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized atom type: " << atomKey << " (" << i << ")"<< std::endl;
        }

        paramVect[i] = theParams;
      }
      
      return std::make_pair(paramVect,foundAll);
    }

  }
}
