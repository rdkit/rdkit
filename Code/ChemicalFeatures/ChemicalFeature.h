//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
    
    // returns the feature id
    virtual const int getId() const = 0;    
    
    // returns the type of the feature
    virtual const std::string& getType() const = 0;
    
    // returns the family of the feature
    virtual const std::string& getFamily() const = 0;
    
    // returns the position of the feature
    virtual RDGeom::Point3D getPos() const = 0;
  };
}

#endif
    
