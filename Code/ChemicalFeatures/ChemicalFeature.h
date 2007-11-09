//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __CHEMICALFEATURE_H_11012005_1310__
#define __CHEMICALFEATURE_H_11012005_1310__

#include <Geometry/point.h>
namespace ChemicalFeatures {
  
  //------------------------------------------------------------------
  //! abstract base class for chemical feature 
  class ChemicalFeature {
  public:
    ChemicalFeature() {};
    virtual ~ChemicalFeature() {};
    
    // returns the type of the feature
    virtual const std::string& getType() const = 0;
    
    // returns the family of the feature
    virtual const std::string& getFamily() const = 0;
    
    // returns the position of the feature
    virtual RDGeom::Point3D getPos() const = 0;
  };
}

#endif
    
