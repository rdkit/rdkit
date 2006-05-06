// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_METRICMATRIXCAL_H__
#define __RD_METRICMATRIXCAL_H__

#include "MetricFuncs.h"
#include <RDGeneral/Invariant.h>

namespace RDDataManip {
  
  /*! \brief A generic metric matrix calculator (e.g simialrity matrix or
   *         ditance matrix) 
   *
   *  This templated class needs some explanation\n
   *    vectType is a container that can support [] operator \n
   *    entryType is the type of entry that is returned by the [] operator\n
   *  Examples of the container include PySequenceHolder which is wrapper around 
   *  a python sequence objects like lists and tuples. \n
   *  Examples of the entryType include a sequence of double, floats, and Explicit Bit Vectors 
   *
   */
  template <class vectType, class entryType> class MetricMatrixCalc {
  public:
    /*! \brief Default Constructor
     *
     */
    MetricMatrixCalc() {};
    
    /*! \brief Set the metric function
     *
     * Set the pointer to the mertic funvtion to be used by the metric calculator
     *
     * ARGUMENTS:
     *
     *  mFunc - pointer to the metric funtion
     */
    void setMetricFunc(double (*mFunc)(const entryType &, const entryType &, int)) {
      dp_metricFunc = mFunc;
    }

    /*! \brief The calculator function
     *
     * ARGUMENTS:
     *
     *  descrips - vectType container with a entryType for each item\n
     *  nItems - the number of item in the descripts.
     *           In several cases this argument is irrelvant since vectType probably supports
     *           a size() member function, But we would like this interface to take for example 
     *           a double** as the and correctly parse the row and columns. \n
     *  dim - the dimension of the 
     *  distMat - pointer to an array to write the distance matrix to
     *            it is assumed that the right sized array has already be allocated.
     *
     * FIX: we can probably make this function create the correct sized distMat and return
     * it to the caller, but when pushing he result out to a python array not sure how to
     * avoid copy the entire distance matrix in that case
     *
     * RETURNS:
     * 
     *  pointer to a 1D array of doubles. Only the lower triangle elements are
     *  included in the array
     */
    void calcMetricMatrix(const vectType &descripts, int nItems, int dim,
                          double *distMat) {
      int i, j, itab, id;
      CHECK_INVARIANT(distMat, "invalid pointer to a distance matix");
      
      for (i = 1; i < nItems; i++) {
        //itab = nItems*i - (i+1)*(i+2)/2;
        itab = i*(i-1)/2;
        for (j = 0; j < i; j++) {
          id = itab + j;
          distMat[id] = dp_metricFunc(descripts[i], descripts[j], dim);
        }
      }
    };
    
  private:
    // pointer to the metric function
    /*! \brief pointer to the metric function
     *
     * In several cases the last argument 'dim' should be irrelevant, 
     * For example when entryType is a bit vector the size is of the vector
     * or the dimension can be obtained by asking the bit vector itself. However
     * we woul like this interface to support other containers lines double* 
     * in which case the 'dim' value is useful in cumputing the metric.
     */  
    double (*dp_metricFunc)(const entryType &, const entryType &, int);
    
  };
};

#endif
