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
    
    virtual const std::string& getType() const = 0;
    
    virtual const std::string& getFamily() const = 0;
    
    virtual RDGeom::Point3D getPos() const = 0;
  };
}

#endif
    
