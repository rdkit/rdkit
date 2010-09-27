//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __CHEMICALFEATUREFACTORY_H_02122004_1545__
#define __CHEMICALFEATUREFACTORY_H_02122004_1545__

#include "MolChemicalFeatureDef.h"
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace RDKit {
  class MolChemicalFeature;
  typedef boost::shared_ptr<MolChemicalFeature> FeatSPtr;
  typedef std::list< FeatSPtr > FeatSPtrList;
  typedef FeatSPtrList::iterator FeatSPtrList_I;

  //! The class for finding chemical features in molecules
  class MolChemicalFeatureFactory {
  public:

    //! returns the number of feature definitions
    int getNumFeatureDefs() const {return d_featDefs.size();};

    //! returns an iterator referring to the first feature definition
    MolChemicalFeatureDef::CollectionType::iterator
    beginFeatureDefs() { return d_featDefs.begin(); };
    //! returns an iterator referring to the end of the feature definitions
    MolChemicalFeatureDef::CollectionType::iterator
    endFeatureDefs() { return d_featDefs.end(); };
    
    //! returns a const_iterator referring to the first feature definition
    MolChemicalFeatureDef::CollectionType::const_iterator
    beginFeatureDefs() const { return d_featDefs.begin(); };
    //! returns a const_iterator referring to the end of the feature definitions
    MolChemicalFeatureDef::CollectionType::const_iterator
    endFeatureDefs() const { return d_featDefs.end(); };
    
    //! appends a feature definition to the collection of features defs.
    void addFeatureDef(MolChemicalFeatureDef::CollectionType::value_type featDef){
      d_featDefs.push_back(featDef);
    }

    //! returns a list of features on the molecule
    /*!
      \param mol          The molecule of interest
      \param includeOnly  (optional) if this is non-null, only features in this
                          family will be returned
    */
    FeatSPtrList getFeaturesForMol(const ROMol &mol,const char *includeOnly="") const;

  private:
    MolChemicalFeatureDef::CollectionType d_featDefs;
  };

  //! constructs a MolChemicalFeatureFactory from the data in a stream
  MolChemicalFeatureFactory *buildFeatureFactory(std::istream &inStream);
  //! constructs a MolChemicalFeatureFactory from the data in a string
  MolChemicalFeatureFactory *buildFeatureFactory(const std::string &featureData);
  
}// end of namespace RDKit

#endif
