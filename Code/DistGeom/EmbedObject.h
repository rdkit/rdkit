//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_EMBED_OBJECT_H__
#define __RD_EMBED_OBJECT_H__

#include <RDGeneral/Invariant.h>
#include <RDGeneral/types.h>
#include <boost/smart_ptr.hpp>

namespace DistGeom {

  /*! \brief Class used to store a quartet of points and chiral volume bounds on them
   * 
   */
  class ChiralSet {
  public:

    unsigned int d_idx1;
    unsigned int d_idx2;
    unsigned int d_idx3;
    unsigned int d_idx4;
    double d_volumeLowerBound;
    double d_volumeUpperBound;

    ChiralSet(unsigned int pid1, unsigned int pid2, unsigned int pid3, unsigned int pid4, 
                 double lowerVolBound, double upperVolBound) {
      CHECK_INVARIANT(lowerVolBound <= upperVolBound, "Inconsistent bounds\n");
      d_idx1 = pid1;
      d_idx2 = pid2;
      d_idx3 = pid3;
      d_idx4 = pid4;
      d_volumeLowerBound = lowerVolBound;
      d_volumeUpperBound = upperVolBound;
    }
    
    inline double getUpperVolumeBound() const {
      return d_volumeUpperBound;
    }
    
    inline double getLowerVolumeBound() const {
      return d_volumeLowerBound;
    }
    
    
  };
 
  typedef boost::shared_ptr<ChiralSet> ChiralSetPtr;
  typedef std::vector<ChiralSetPtr> VECT_CHIRALSET;

  /*! \brief Class that contains a an object that needs to be embedded. 
   *
   *  Basially a pointer to the MetricMatrix matric consisting of the upper and 
   *  lower bounds on the distances, and a list of quartets of points to which chirality 
   *  constraints need to be applied.
   */
  /*
  class EmbedObject {
  public:
    EmbedObject() {};
   
    inline void setMetricMatrix(MetricMatrix *metricMat) {
      PRE_CONDITION(metricMat, "");
      d_metricMat = MetricMatPtr(metricMat);
    }

    inline void setMetricMatrix(MetricMatPtr metricMatPtr) {
      d_metricMat = metricMatPtr;
    }
    
    inline const MetricMatPtr getMetricMatrix() const {
      return d_metricMat;
    }

    inline MetricMatPtr getMetricMatrix() {
      return d_metricMat;
    }

    inline void addChiralSet(ChiralSet *cset) {
      PRE_CONDITION(cset, "");
      d_chiralSets.push_back(ChiralSetPtr(cset));
    }

    inline void addChiralSet(ChiralSetPtr cset) {
      d_chiralSets.push_back(cset);
    }

    const VECT_CHIRALSET &getChiralSets() const {
      return d_chiralSets;
    }

  private:
    MetricMatPtr d_metricMat;
    VECT_CHIRALSET d_chiralSets;
    };*/

}  

#endif
