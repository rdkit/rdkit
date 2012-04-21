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
#include "MolChemicalFeature.h"
#include "MolChemicalFeatureDef.h"
#include "MolChemicalFeatureFactory.h"
#include "FeatureParser.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/ROMol.h>

#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>
#include <sstream>
#include <set>
#include <algorithm>


namespace RDKit {

  FeatSPtrList MolChemicalFeatureFactory::getFeaturesForMol(const ROMol &mol,
                                                            const char* includeOnly) const {
    PRECONDITION(includeOnly,"bad limits");
    std::string limits(includeOnly);
    
#ifdef USE_VFLIB
    AR_MOLGRAPH *molG=getMolGraph(mol);
#endif
    FeatSPtrList res;
    int idx = 1;
    typedef std::vector< std::pair< std::string,std::set<int> > > MatchSetCollection;
    MatchSetCollection matchSets;
    for(MolChemicalFeatureDef::CollectionType::const_iterator featDefIt=beginFeatureDefs();
        featDefIt!=endFeatureDefs();featDefIt++){
      MolChemicalFeatureDef::CollectionType::value_type featDef=*featDefIt;
      if(limits=="" || limits==featDef->getFamily()){
        std::vector< MatchVectType > matches;
#ifdef USE_VFLIB
        unsigned int numMatches=SubstructMatch(molG,*featDef->getPattern(),matches);
#else
        unsigned int numMatches=SubstructMatch(mol,*featDef->getPattern(),matches);
#endif
        for(unsigned int i=0;i<numMatches;i++){
          const MatchVectType &match=matches[i];
          std::set<int> matchSet;
          for(MatchVectType::const_iterator mIt=match.begin();
              mIt!=match.end();
              ++mIt){
            matchSet.insert(mIt->second);
          }
          
          // loop over the matches we've already found and see if this one
          // is unique:
          bool unique=true;
          for(MatchSetCollection::const_iterator vsiCI=matchSets.begin();
              vsiCI!=matchSets.end();
              ++vsiCI){
            if(vsiCI->first==featDef->getFamily() &&
               std::includes(vsiCI->second.begin(),vsiCI->second.end(),
                             matchSet.begin(),matchSet.end())){
              unique=false;
              break;
            }
          }
          if(unique){
            matchSets.push_back(std::make_pair(featDef->getFamily(),matchSet));
    
            // Set up the feature:
            MolChemicalFeature *newFeat=new MolChemicalFeature(&mol,this,featDef.get(),idx++);
            MolChemicalFeature::AtomPtrContainer &atoms=newFeat->d_atoms;
            atoms.resize(match.size());
    
            // set up the atoms:
            for(MatchVectType::const_iterator matchIt=match.begin();
                        matchIt!=match.end();matchIt++){
              int atomIdx=matchIt->second;
              int queryIdx=matchIt->first;
              atoms[queryIdx]=mol.getAtomWithIdx(atomIdx);
            }
    
            // finally, add this to our result:
            res.push_back(FeatSPtrList::value_type(newFeat));
          }
        }
      }
    }
#ifdef USE_VFLIB
#ifndef CACHE_ARMOLGRAPHS
    delete molG;
#endif
#endif
    return res;
  }
  
  MolChemicalFeatureFactory *buildFeatureFactory(const std::string &featureData){
    std::stringstream ss(featureData);
    return buildFeatureFactory(ss);
  }

  MolChemicalFeatureFactory *buildFeatureFactory(std::istream &inStream){
    MolChemicalFeatureFactory *res=0;
    MolChemicalFeatureDef::CollectionType featDefs;

    if(parseFeatureData(inStream,featDefs)==0){
      // everything parsed ok
      res = new MolChemicalFeatureFactory();
      //std::copy(featDefs.begin(),featDefs.end(),res->beginFeatureDefs());
      for(MolChemicalFeatureDef::CollectionType::const_iterator ci=featDefs.begin();
          ci!=featDefs.end();ci++){
        res->addFeatureDef(*ci);
      }
    }

    return res;
  }
}
