//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __FEATURE_H_30112004_1121__
#define __FEATURE_H_30112004_1121__

#include <vector>
#include <Geometry/point.h>

namespace RDFeatures {
  template <typename FAMILYMARKER, typename TYPEMARKER=FAMILYMARKER, typename LOCTYPE=RDGeom::Point3D>
  class ExplicitFeature {
  public:
    ExplicitFeature() {};
    explicit ExplicitFeature(const FAMILYMARKER &f,const TYPEMARKER &t) :
      d_family(f), d_type(t) {};
    ExplicitFeature(const FAMILYMARKER &f,const TYPEMARKER &t,const LOCTYPE &loc) :
      d_family(f), d_type(t), d_loc(loc){};

    const FAMILYMARKER &getFamily() const { return d_family; };
    void setFamily(const FAMILYMARKER &f)  { d_family=f; };

    const TYPEMARKER &getType() const { return d_type; };
    void setType(const TYPEMARKER &t)  { d_type=t; };

    const LOCTYPE &getLoc() const { return d_loc; };
    void setLoc(const LOCTYPE &loc) { d_loc=loc; };

    const std::vector<LOCTYPE> &getDirs() const { return d_dirs; };
    std::vector<LOCTYPE> &getDirs() { return d_dirs; };
    
  private:
    FAMILYMARKER d_family;
    TYPEMARKER d_type;
    LOCTYPE d_loc;
    std::vector<LOCTYPE> d_dirs;
  };


  template <typename FAMILYMARKER, typename TYPEMARKER=FAMILYMARKER, typename LOCTYPE=RDGeom::Point3D>
  class ImplicitFeature {
  public:
    ImplicitFeature() : d_weightSum(0.0) {};
    explicit ImplicitFeature(const FAMILYMARKER &f,const TYPEMARKER &t) :
      d_weightSum(0.0), d_family(f), d_type(t) {};

    const FAMILYMARKER &getFamily() const { return d_family; };
    void setFamily(const FAMILYMARKER &f)  { d_family=f; };

    const TYPEMARKER &getType() const { return d_type; };
    void setType(const TYPEMARKER &t)  { d_type=t; };

    LOCTYPE getLoc() const {
      PRECONDITION(d_weights.size()==d_locs.size(),"weight/locs mismatch");
      LOCTYPE accum;
      for(unsigned int i=0;i<d_weights.size();i++){
	LOCTYPE tmp=*d_locs[i];
	tmp *= d_weights[i]/d_weightSum;
	accum += tmp;
      }
      return accum;
    };
    void addPoint(const LOCTYPE *p,double weight=1.0){
      d_locs.push_back(p);
      d_weights.push_back(weight);
      d_weightSum += weight;
    }
    void reset() {
      d_locs.clear();
      d_weights.clear();
      d_weightSum=0.0;
    }

    const std::vector<LOCTYPE> &getDirs() const { return d_dirs; };
    std::vector<LOCTYPE> &getDirs() { return d_dirs; };
    

  private:
    double d_weightSum;
    FAMILYMARKER d_family;
    TYPEMARKER d_type;
    std::vector<double>  d_weights;
    std::vector<const LOCTYPE *> d_locs;
    // FIX: add something correct for directions
    std::vector<LOCTYPE> d_dirs;
  };
}
#endif
