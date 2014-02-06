//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_CORRMATGENERATOR_H_
#define _RD_CORRMATGENERATOR_H_

#include <RDGeneral/types.h>
#include <DataStructs/BitVects.h>
#include <boost/dynamic_bitset.hpp>

namespace RDInfoTheory {
  //FIX: won't worry about it now, but this class can be templated by the type of 
  // container for the bit list and type of descriptors (fingerprint vs. real valued)
  class BitCorrMatGenerator {
    /*! \brief A class to generate a correlation matrix for a bunch of fingerprints
     *
     *  The correlation matrix is done only for the bit IDs that are set by a call to the 
     *  function setDescriptorIdList
     *  
     *    cr = CorrMatGenerator();
     *    cr.setDescriptorIdList(descList);
     *    for each fingerprint in list of fingerprints {
     *        cr.collectVotes(fingerprint);
     *    }
     *    double *corrMat = cr.getCorrMat()
     *  
     *  The resulting correlation matrix is a one dimension matrix with only the lower triangle elements
     *  of the symmetric matrix
     */
  public:
    BitCorrMatGenerator() {
      this->initGenerator();
    };

    ~BitCorrMatGenerator() {
      delete [] dp_corrMat;
    }

    void initGenerator() {
      dp_corrMat = 0;
      d_descs.resize(0);
      d_nExamples = 0;
    };

    /*! \brief Set the list bits that we are interested in correlating
     *
     *  \param bitIdList is a list of bit ids that need to be correlated e.g. a list top ranked ensemble 
     *  of bits 
     */
    void setBitIdList(const RDKit::INT_VECT &bitIdList) {
      d_descs = bitIdList;
      int i, nd = d_descs.size();
      int nelem = nd*(nd-1)/2;
      delete [] dp_corrMat;

      dp_corrMat = new double[nd*(nd-1)/2];
      for (i = 0; i < nelem; i++) {
        dp_corrMat[i] = 0.0;
      }
    };

    //! \brief get the number of examples we used so far to compute the correlation matrix
    int getNumExamples() const {
      return d_nExamples;
    };

    //! \brief Get the list of bits ID that are used to generate the correlation matrix
    RDKit::INT_VECT getCorrBitList() const {
      return d_descs;
    };

    //! \brief Gets a pointer to the correlation matrix
    double *getCorrMat() {
      return dp_corrMat;
    };
    
    //! \brief For each pair of on bits (bi, bj) in fp increase the correlation count
    //    for the pair by 1
    void collectVotes(const BitVect &fp) {
      unsigned int nd = d_descs.size();
      // use a temporary bit vector to first mask the fingerprint
      ExplicitBitVect ebv(nd);
      int bi;
      for (unsigned int i = 0; i < nd; i++) {
        bi = d_descs[i];
        if (fp[bi]) {
          ebv.setBit(i);
        }
      }
      for (unsigned i = 1; i < nd; i++) {
        unsigned int itab = i*(i-1)/2;
        if (ebv[i]) {
          for (unsigned int j = 0; j < i; j++) {
            if ( ebv[j]) {
              dp_corrMat[itab + j] += 1;
            }
          }
        }
      }
      d_nExamples++;
    };

  private:
    RDKit::INT_VECT d_descs;
    double *dp_corrMat;
    int d_nExamples;
  };

}

#endif
    

