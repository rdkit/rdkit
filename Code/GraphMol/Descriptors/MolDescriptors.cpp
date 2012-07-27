// $Id$
//
//  Copyright (C) 2005-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include "MolDescriptors.h"
#include <map>
#include <list>
#include <algorithm>
#include <sstream>

namespace RDKit{
  namespace Descriptors {

    static std::string _amwVersion="1.0.0";
    double calcAMW(const ROMol &mol,bool onlyHeavy){
      double res=0.0;
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt != mol.endAtoms();++atomIt){
	int atNum=(*atomIt)->getAtomicNum();
	if(atNum!=1 || !onlyHeavy){
	  res += (*atomIt)->getMass();
	}

	// add our implicit Hs if we need to:
	if(!onlyHeavy){
          const PeriodicTable *table=PeriodicTable::getTable();
	  res += (*atomIt)->getTotalNumHs()*table->getAtomicWeight(1);
	}
      }
      return res;
    }

    static std::string _exactmwVersion="1.1.0";
    double calcExactMW(const ROMol &mol,bool onlyHeavy){
      double res=0.0;
      int nHsToCount=0;
      const PeriodicTable *table=PeriodicTable::getTable();
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt != mol.endAtoms();++atomIt){
	int atNum=(*atomIt)->getAtomicNum();
	if(atNum!=1 || !onlyHeavy){
          if(atNum==1) --nHsToCount;
          if(!(*atomIt)->getIsotope()){
            res += table->getMostCommonIsotopeMass(atNum);
          } else {
            res += (*atomIt)->getMass();
          }
	}

	// add our implicit Hs if we need to:
	if(!onlyHeavy){
          nHsToCount += (*atomIt)->getTotalNumHs(true);
	}

      }
      if(!onlyHeavy){
        res += nHsToCount*table->getMostCommonIsotopeMass(1);
      }
      return res;
    }

    static std::string _molFormulaVersion="1.2.0";
    std::string calcMolFormula(const ROMol &mol,bool separateHIsotopes){
      std::ostringstream res;
      std::map<std::string,unsigned int> counts;
      int charge=0;
      unsigned int nHs=0;
      const PeriodicTable *table=PeriodicTable::getTable();
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt != mol.endAtoms();++atomIt){
	int atNum=(*atomIt)->getAtomicNum();
        std::string symb;
        symb=table->getElementSymbol(atNum);
        if(atNum == 1 && separateHIsotopes){
          if((*atomIt)->getIsotope()==2) symb="D";
          else if((*atomIt)->getIsotope()==3) symb="T";
        }
        if(counts.find(symb)!=counts.end()){
          counts[symb]+=1;
        } else {
          counts[symb]=1;
        }
        nHs += (*atomIt)->getTotalNumHs();
        charge += (*atomIt)->getFormalCharge();
      }

      if(nHs){
        if(counts.find(std::string("H"))!=counts.end()){
          counts["H"]+=nHs;
        } else {
          counts["H"]=nHs;
        }
      }
      std::list<std::string> ks;
      for(std::map<std::string,unsigned int>::const_iterator countIt=counts.begin();
          countIt!=counts.end();++countIt){
        ks.push_back(countIt->first);
      }
      ks.sort();
      // put in Hill order:
      if(counts.find(std::string("C"))!=counts.end()){
        ks.remove(std::string("C"));
        if(separateHIsotopes){
          if(counts.find(std::string("T"))!=counts.end()){
            ks.remove(std::string("T"));
            ks.push_front(std::string("T"));
          }
          if(counts.find(std::string("D"))!=counts.end()){
            ks.remove(std::string("D"));
            ks.push_front(std::string("D"));
          }
        }
        if(counts.find(std::string("H"))!=counts.end()){
          ks.remove(std::string("H"));
          ks.push_front(std::string("H"));
        }
        ks.push_front(std::string("C"));
      }

      for(std::list<std::string>::const_iterator kIter=ks.begin();
          kIter!=ks.end();++kIter){
        res << *kIter;
        if(counts[*kIter]>1) res<<counts[*kIter];
      }
      if(charge>0){
        res<<"+";
        if(charge>1) res<<charge;
      } else if(charge<0){
        res<<"-";
        if(charge<-1) res<<-1*charge;
      }
      return res.str();
    }

    
  } // end of namespace Descriptors
} // end of namespace RDKit
