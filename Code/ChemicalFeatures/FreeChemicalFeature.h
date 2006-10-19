//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __FREECHEMICALFEATURE_H_13012005_1023__
#define __FREECHEMICALFEATURE_H_13012005_1023__

#include <Geometry/point.h>
#include <ChemicalFeatures/ChemicalFeature.h>

namespace ChemicalFeatures {

  //------------------------------------------------------
  //! Class for chemical features that do orignate from molecules
  //  e.g. pharamcophores, site-maps etc.
  class FreeChemicalFeature : public ChemicalFeature {
    public:
    FreeChemicalFeature(std::string family, std::string type,
                        const RDGeom::Point3D &loc) :
      d_family(family), d_type(type), d_position(loc) {
    }

    FreeChemicalFeature(std::string family, const RDGeom::Point3D &loc) :
      d_family(family), d_type(""), d_position(loc) {
    }

    explicit FreeChemicalFeature(const std::string &pickle) {
      this->initFromString(pickle);
    }
    FreeChemicalFeature() :
      d_family(""), d_type(""), d_position(RDGeom::Point3D(0.0, 0.0, 0.0)) {
    }

    FreeChemicalFeature(const FreeChemicalFeature &other) : 
      d_family(other.getFamily()), d_type(other.getType()), d_position(other.getPos())  {
    }

    ~FreeChemicalFeature() {}

    const std::string& getFamily() const {
      return d_family;
    }

    const std::string& getType() const {
      return d_type;
    }

    RDGeom::Point3D getPos() const {
      return d_position;
    }

    void setFamily(const std::string &family) {
      d_family = family;
    }
    
    void setType(const std::string &type) {
      d_type = type;
    }

    void setPos(const RDGeom::Point3D &loc) {
      //std::cout << loc.x << " " << loc.y << " " << loc.z << "\n";
      d_position = loc;
      //std::cout << d_position.x << " " << d_position.y << " " << d_position.z << "\n";
    }

    std::string toString() const; 
    void initFromString(const std::string &pickle);
    
  private:
    std::string d_family;
    std::string d_type;
    RDGeom::Point3D d_position;
  };
}
  

#endif

