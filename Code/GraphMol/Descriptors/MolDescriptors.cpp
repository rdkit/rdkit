// $Id$
//
//  Copyright (C) 2005-2012 Greg Landrum and Rational Discovery LLC
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
          if(!(*atomIt)->getIsotope()){
            res += table->getMostCommonIsotopeMass(atNum);
          } else {
            res += (*atomIt)->getMass();
          }
	}

	// add our implicit Hs if we need to:
	if(!onlyHeavy){
          nHsToCount += (*atomIt)->getTotalNumHs(false);
	}

      }
      if(!onlyHeavy){
        res += nHsToCount*table->getMostCommonIsotopeMass(1);
      }
      return res;
    }

    namespace {
      bool HillCompare(const std::pair<unsigned int,std::string> &v1,
                       const std::pair<unsigned int,std::string> &v2){

        bool nCompare=(v1.first<v2.first);

        // special cases: Cs, Hs, Ds, and Ts go at the beginning
        if(v1.second=="C"){
          if(v2.second!="C") return true;
          return nCompare;
        } else if(v2.second=="C") return false;

        if(v1.second=="H"){
          if(v2.second!="H") return true;
          return nCompare;
        } else if(v2.second=="H") return false;        

        if(v1.second=="D") return true;
        else if(v2.second=="D") return false;

        if(v1.second=="T") return true;
        else if(v2.second=="T") return false;

        // finally, just compare the symbols and the isotopes:
        if(v1!=v2) return v1<v2;
        else return nCompare;
      }

    }

    static std::string _molFormulaVersion="1.3.0";
    std::string calcMolFormula(const ROMol &mol,bool separateIsotopes,bool abbreviateHIsotopes){
      std::ostringstream res;
      std::map<std::pair<unsigned int,std::string>,unsigned int> counts;
      int charge=0;
      unsigned int nHs=0;
      const PeriodicTable *table=PeriodicTable::getTable();
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt != mol.endAtoms();++atomIt){
	int atNum=(*atomIt)->getAtomicNum();
        std::pair<unsigned int,std::string> key=std::make_pair(0,table->getElementSymbol(atNum));
        if(separateIsotopes){
          unsigned int isotope = (*atomIt)->getIsotope();
          if(abbreviateHIsotopes && atNum==1 && (isotope==2 || isotope==3) ){
            if(isotope==2) key.second="D";
            else key.second="T";
          } else {
            key.first=isotope;
          }
        }
        if(counts.find(key)!=counts.end()){
          counts[key]+=1;
        } else {
          counts[key]=1;
        }
        nHs += (*atomIt)->getTotalNumHs();
        charge += (*atomIt)->getFormalCharge();
      }

      if(nHs){
        std::pair<unsigned int,std::string> key=std::make_pair(0,"H");
        if(counts.find(key)!=counts.end()){
          counts[key]+=nHs;
        } else {
          counts[key]=nHs;
        }
      }
      std::list<std::pair<unsigned int,std::string> > ks;
      for(std::map<std::pair<unsigned int,std::string>,unsigned int>::const_iterator countIt=counts.begin();
          countIt!=counts.end();++countIt){
        ks.push_back(countIt->first);
      }
      ks.sort(HillCompare);
      
      for(std::list<std::pair<unsigned int,std::string> >::const_iterator kIter=ks.begin();
          kIter!=ks.end();++kIter){
        const std::pair<unsigned int,std::string> &key=*kIter;
        if(key.first>0){
          res<<"["<<key.first<<key.second<<"]";
        } else {
          res<<key.second;
        }
        if(counts[key]>1) res<<counts[key];
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
