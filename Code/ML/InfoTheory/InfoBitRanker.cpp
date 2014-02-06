// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and  Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "InfoBitRanker.h"
#include "InfoGainFuncs.h"
#include <RDGeneral/Invariant.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <RDGeneral/FileParseException.h>
#include <RDBoost/Exceptions.h>
#include <algorithm>
#include <queue>

namespace RDInfoTheory {
  typedef std::pair<double, int> PAIR_D_I;
  typedef std::vector<PAIR_D_I> VECT_PDI;

  struct gtDIPair {
    bool operator() ( const PAIR_D_I &pd1, const PAIR_D_I &pd2) const {
      return pd1.first > pd2.first;
    }
  };

  typedef std::priority_queue<PAIR_D_I, VECT_PDI, gtDIPair> PR_QUEUE;


  
  void InfoBitRanker::setBiasList(RDKit::INT_VECT &classList) {
    RANGE_CHECK(0, classList.size(), d_classes);
    d_biasList = classList;
    //make sure we don't have any duplicates
    std::sort(d_biasList.begin(), d_biasList.end());
    RDKit::INT_VECT_CI bi = std::unique(d_biasList.begin(), d_biasList.end());
    CHECK_INVARIANT(bi == d_biasList.end(), "There are duplicates in the class bias list");

    // finally make sure all the class ID in d_biasList are within range
    for (bi = d_biasList.begin(); bi != d_biasList.end(); bi++) {
      RANGE_CHECK(0, static_cast<unsigned int>(*bi), d_classes-1);
    }
  }

  void InfoBitRanker::setMaskBits(RDKit::INT_VECT &maskBits) {
    delete dp_maskBits;
    dp_maskBits = new ExplicitBitVect(d_dims);
    for (RDKit::INT_VECT_CI bi = maskBits.begin();
         bi != maskBits.end(); ++bi) {
      dp_maskBits->setBit(*bi);
    }
  }

  bool InfoBitRanker::BiasCheckBit(RDKit::USHORT *resMat) const {
    PRECONDITION(resMat,"bad results pointer");
    if ((d_biasList.size() == 0) || (d_biasList.size() == d_classes)) {
      //we will accept the bit 
      return true;
    }
    RDKit::DOUBLE_VECT fracs;
    fracs.resize(d_classes);

    // compute the fractions of items in each class that hit the bit
    // and record the maximum for the those classes not in the bias list
    double maxCor = 0.0;
    for (unsigned int i = 0; i < d_classes; i++) {
      if (d_clsCount[i] > 0) {
        fracs[i] = ((double)resMat[i])/d_clsCount[i];
      } else {
        fracs[i] = 0.0;
      }
      if (std::find(d_biasList.begin(), d_biasList.end(), i) == d_biasList.end()) {
        // if not in the biasList
        if (fracs[i] > maxCor) {
          // if this is fraction is greater than the previously known maximum
          maxCor = fracs[i];
        }
      }
    }

    bool bitOk = false;
    for (RDKit::INT_VECT_CI bci = d_biasList.begin(); bci !=
           d_biasList.end(); ++bci) {
      if (fracs[*bci] >= maxCor) {
        bitOk = true;
        break;
      }
    }
    return bitOk;
  }

      
  double InfoBitRanker::BiasChiSquareGain(RDKit::USHORT *resMat) const {
    PRECONDITION(resMat,"bad result pointer");
    bool bitOk = this->BiasCheckBit(resMat);
    double info=0.0;
    if (bitOk) {
      info = ChiSquare(resMat, 2, d_classes);
    }
    return info;
  }

  double InfoBitRanker::BiasInfoEntropyGain(RDKit::USHORT *resMat) const {
    PRECONDITION(resMat,"bad result pointer");
    bool bitOk = this->BiasCheckBit(resMat);
    double info=0.0;
    if (bitOk) {
      info = InfoEntropyGain(resMat, 2, d_classes);
    }
    return info;
  }

  void InfoBitRanker::accumulateVotes(const ExplicitBitVect &bv, unsigned int label) {
    RANGE_CHECK(0, label, d_classes-1);
    CHECK_INVARIANT(bv.getNumBits() == d_dims, "Incorrect bit vector size");

    d_nInst += 1;
    d_clsCount[label] += 1;
    for (unsigned int i=0;i<bv.getNumBits();i++){
      if( (*bv.dp_bits)[i] && (!dp_maskBits || dp_maskBits->getBit(i)) ){
        d_counts[label][i] += 1;
      }
    }
  }
  
  void InfoBitRanker::accumulateVotes(const SparseBitVect &bv, unsigned int label) {
    RANGE_CHECK(0, label, d_classes-1);
    CHECK_INVARIANT(bv.getNumBits() == d_dims, "Incorrect bit vector size");

    d_nInst += 1;
    d_clsCount[label] += 1;
    for (IntSet::const_iterator obi = bv.dp_bits->begin();
         obi != bv.dp_bits->end();
        ++obi)  {
      if(!dp_maskBits || dp_maskBits->getBit(*obi)){
        d_counts[label][(*obi)] += 1;
      }
    }
  }
  
  double *InfoBitRanker::getTopN(unsigned int num) {
    // this is a place holder to pass along to infogain function
    // the size of this container should nVals*d_classes, where nVals
    // is the number of values a variable can take.
    // since we are dealing with a binary bit vector nVals = 2
    // in addition the infogain function pretends that this is a 2D matrix
    // with the number of rows equal to nVals and num of columns equal to 
    // d_classes
    if(num>d_dims) throw ValueErrorException("attempt to rank more bits than present in the bit vectors");
    if(dp_maskBits)
      CHECK_INVARIANT(num <= dp_maskBits->getNumOnBits(), "Can't rank more bits than the ensemble size"); 
    RDKit::USHORT *resMat = new RDKit::USHORT[2*d_classes];
    
    PR_QUEUE topN;

    for (unsigned int i = 0; i < d_dims; i++) {
      // we may want to ignore bits that are not turned on in any item of class 
      // "ignoreNoClass"
      /*
      if ((0 <= ignoreNoClass) && (d_classes > ignoreNoClass)) {
        if (d_counts[ignoreNoClass][i] == 0) {
          continue;
        }
        }*/
      
      
      if (dp_maskBits && !dp_maskBits->getBit(i)) {
           continue;
      }

      // fill up dmat
      for (unsigned int j = 0; j < d_classes; j++) {
        // we know that we have only two rows here
        resMat[j] = d_counts[j][i];
        resMat[d_classes + j] = (d_clsCount[j] - d_counts[j][i]);
      }
      double info = 0.0;
      switch (d_type) {
      case ENTROPY:
        info = InfoEntropyGain(resMat, 2, d_classes);
        break;
      case BIASENTROPY:
        info = this->BiasInfoEntropyGain(resMat);
        break;
      case CHISQUARE:
        info = ChiSquare(resMat, 2, d_classes);
        break;
      case BIASCHISQUARE:
        info = BiasChiSquareGain(resMat);
        break;
      default:
        break;
      }

      PAIR_D_I entry(info, i);
      
      if (info >= 0.0) {
        if (topN.size() < num) {
          topN.push(entry);
        }
        else if (info > topN.top().first) {
          topN.pop();
          topN.push(entry);
        }
      }
    }
    
    delete [] resMat;
    
    // now fill up the result matrix for the topN bits
    // the result from this function is a double * of size 
    // num*4. The caller of this function interprets this
    // array as a two dimensional array of size num*(2+d_classes) with each row
    // containing the following entries 
    //   bitId, infogain, 1 additional column for number of hits for each class
    //double *res = new double[num*(2+d_classes)];
    
    d_top = num;
    int ncols = 2+d_classes;
    
    delete [] dp_topBits;
    dp_topBits = new double[num*ncols];
    
    int offset, bid;
    
    RDKit::INT_VECT maskBits;
    if (dp_maskBits && topN.size() < num) {
      dp_maskBits->getOnBits(maskBits);
    }

    for (int i = num - 1; i >= 0; i--) {
      offset = i*ncols;
      if (topN.size() == 0 ) {
        if (dp_maskBits) {
              bid = maskBits[i];
        } else {
              bid = i;
        }
        dp_topBits[offset + 1] = 0.0;
      } else {
        bid = topN.top().second; // bit id
        dp_topBits[offset + 1] = topN.top().first; // value of the infogain
        topN.pop();
      }
      dp_topBits[offset] = (double)bid;
      
      for (unsigned int j = 0; j < d_classes; j++) {
        dp_topBits[offset + 2 + j] = (double)d_counts[j][bid];
      }
    }
    return dp_topBits;
  }

  void InfoBitRanker::writeTopBitsToStream(std::ostream *outStream) const {
    (*outStream) << std::setw(12) << "Bit" << std::setw(12) << "InfoContent";
    for (unsigned int ic = 0; ic < d_classes; ic++) {
      (*outStream) << std::setw(10) << "class" << ic;
    }
    (*outStream) << std::endl;
   
    unsigned int ncols = 2 + d_classes;
    for (unsigned int i = 0; i < d_top; i++) {
      (*outStream) << std::setw(12) << (int)dp_topBits[i*ncols]
                   << std::setw(12) << std::setprecision(5) 
                   << dp_topBits[i*ncols + 1];
      for (unsigned int ic = 0; ic < d_classes; ic++) {
        (*outStream) << std::setw(10) << (int)dp_topBits[i*ncols + 2 + ic];
      }
      (*outStream) << "\n";
    
    }
  }
  
  void InfoBitRanker::writeTopBitsToFile(std::string fileName) const {
    std::ofstream tmpStream(fileName.c_str());
    if ((!tmpStream) || (tmpStream.bad()) ) {
      std::ostringstream errout;
      errout << "Bad output file " << fileName;
      throw RDKit::FileParseException(errout.str());
    }

    std::ostream &outStream = static_cast<std::ostream &>(tmpStream);
    this->writeTopBitsToStream(&outStream);
  }
  
}

    
  
