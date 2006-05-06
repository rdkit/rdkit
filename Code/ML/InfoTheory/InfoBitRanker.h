// $Id$
//
//  Copyright (C) 2003 Rational Discovery LLC
//

#ifndef INFORANKER_H
#define INFORANKER_H

#include <RDGeneral/types.h>
#include <queue>
#include <DataStructs/BitVects.h>
#include <iostream>

/*! \brief Evaluate the info-content for 
 */

typedef std::pair<double, int> PAIR_D_I;
//typedef std::pair<float, int> PAIR_D_I;
typedef std::vector<PAIR_D_I> VECT_PDI;

struct gtDIPair {
  bool operator() ( const PAIR_D_I &pd1, const PAIR_D_I &pd2) const {
    return pd1.first > pd2.first;
  }
};

typedef std::priority_queue<PAIR_D_I, VECT_PDI, gtDIPair> PR_QUEUE;
typedef std::vector<RDKit::USHORT> USHORT_VECT;
typedef std::vector<USHORT_VECT> VECT_USHORT_VECT;


/*! \brief Class used to rank bits based on a specified measure of infomation
 *
 * Basically a premitive mimic of combichem signal
 * The mode of operation 
 *  - create an instance of this calss
 *  - loop over the finger print in the dataset by calling AccumulateVotes function
 *  - call getTopN to get the top n ranked bits
 */
namespace RDInfoTheory {
  class InfoBitRanker {
  public:
    
    /*! \brief the type of measure for information
     * 
     * Currently using only entropy 
   */
    typedef enum {
      ENTROPY=1,
      BIASENTROPY=2,
      CHISQUARE=3,
      BIASCHISQUARE=4} InfoType;
    
    /*! \brief Constructor 
     *
     * ARGUMENTS:
     *
     *   - nBits: the dimension of the bit vectors or the fingerprint length
     *   - nClasses: the number of classes used int eh classification problem (e.g. active,
     *              moderately active, inactive etc.). It is assumed that the classes are 
     *              numbered from 0 to (nClasses - 1) 
     *   - infoType: the type of information metric
     */
    InfoBitRanker(unsigned int nBits, unsigned int nClasses, InfoType infoType=InfoBitRanker::ENTROPY) {
      d_dims = nBits;
      d_classes = nClasses;
      d_type = infoType;
      d_counts.resize(0);
      for (unsigned int i = 0; i < nClasses; i++) {
        USHORT_VECT cCount;
        cCount.resize(d_dims, 0);
        d_counts.push_back(cCount);
      }
      d_clsCount.resize(d_classes, 0);
      d_nInst = 0;
      d_top = 0;
      dp_topBits=0;
      d_biasList.resize(0);
      dp_maskBits=0;
    }
    
    ~InfoBitRanker() {
      delete [] dp_topBits;
      if(dp_maskBits) delete dp_maskBits;
    }

    /*! \brief Accumulate the votes for all the bits turned on in a bit vector
     *  
     *  ARGUMENTS:
     *
     *   - bv : bit vector that supports [] operator
     *   - label : the class label for the bit vector. It is assumed that 0 <= class < nClasses
     */
    void accumulateVotes(const ExplicitBitVect &bv, unsigned int label);
    void accumulateVotes(const SparseBitVect &bv, unsigned int label);
    
    /*! \brief Returns the top n bits ranked by the information metric
     *
     * This is actually the function where most of the work of ranking is happening
     * 
     * ARGUMENTS:
     *  
     *  - num : the number of top ranked bits that are required
     */
    double *getTopN(unsigned int num);//, int ignoreNoClass=-1);
    
    /*! \brief return the number of labelled instances(examples) or fingerprints seen so far
     *
     */
    unsigned int getNumInstances() {
      return d_nInst;
    }
    
    /*! \brief return the number of classes 
     *
     */
    unsigned int getNumClasses() {
      return d_classes;
    }

    /*! \brief Set the classes to which the entropy calculation should be biased
     *
     * This list contains a set of class ids used when in the BIASENTROPY mode of ranking bits. 
     * In this mode, a bit must be correllated higher with one of the biased classes than all the 
     * other classes. For example, in a two class problem with actives and inactives, the fraction of 
     * actives that hit the bit has to be greater than the fraction of inactives that hit the bit
     *
     * ARGUMENTS:
     *   classList - list of class ids that we want a bias towards
     */
    void setBiasList(RDKit::INT_VECT classList);


    /*! \brief Set the bits to be used as a mask
     *
     * If this function is called, only the bits which are present in the
     *   maskBits list will be used.
     *
     * ARGUMENTS:
     *   maskBits - the bits to be considered
     */
    void setMaskBits(RDKit::INT_VECT maskBits);

    /*! \brief Write the top N bits to a stream
     *
     */
    void writeTopBitsToStream(std::ostream *outStream);
    
    /*! \brief Write the top bits to a file
     *
     */
    void writeTopBitsToFile(std::string fileName);

  private:
    /*! \brief check if we want to compute the info content for a bit based on the bias list
     *
     * This what happens here:
     *    - the fraction of items in each class that hit a particular bit are computed
     *    - the maximum of these fractions for classes that are not in the biasList are computed
     *    - If this maximum is less than the fraction for atleast one of classes in the biaslist
     *      the bit is considered good
     * ARGUMENTS:
     *   - resMat : the result matrix, one dimensional matrix of dimension (2*(num of classes))
     *              a 2D structure is assumed with the first row containing number of items of each class 
     *              with the bit set and the second row to entires of each class with the bit turned off
     */
    bool BiasCheckBit(RDKit::USHORT *resMat);

    /*! \brief Compute the biased info entropy gain based on the bias list
     *
     *  This what happens here:
     *    - we call BiasCheckBit to see if the bit qualifies to compute the infocontent
     *    - If this bit is ok then we call InfoEntropyGain otherwise we return 0.0 
     *
     * ARGUMENTS:
     *   - resMat : the result matrix, one dimensional matrix of dimension (2*(num of classes))
     *              a 2D structure is assumed with the first row containing number of items of each class 
     *              with the bit set and the second row to entires of each class with the bit turned off
     */
    double BiasInfoEntropyGain(RDKit::USHORT *resMat);

    /*! \brief Compute the biased chi qsure value based on the bias list
     *
     *  This what happens here:
     *    - we call BiasCheckBit to see if the bit qualifies to compute the infocontent
     *    - If this bit is ok then we call InfoEntropyGain otherwise we return 0.0 
     *
     * ARGUMENTS:
     *   - resMat : the result matrix, one dimensional matrix of dimension (2*(num of classes))
     *              a 2D structure is assumed with the first row containing number of items of each class 
     *              with the bit set and the second row to entires of each class with the bit turned off
     */
    double BiasChiSquareGain(RDKit::USHORT *resMat);

    unsigned int d_dims; // the number of bits in the fingerprints
    unsigned int d_classes; // the number of classes (active, inactive, moderately active etc.)
    InfoType d_type; // the type of information meassure - currently we support only entropy
    VECT_USHORT_VECT d_counts; // place holder of counting the number of hits for each bit for each class
    USHORT_VECT d_clsCount; // counter for the number of instances of each class 
    double *dp_topBits; // storage for the top ranked bits and the corresponding statistics
    unsigned int d_top; // the number of bits that have been ranked
    unsigned int d_nInst; // total number of instances or fingerprints used accumulate votes
    RDKit::INT_VECT d_biasList; // if we want a bias towards certain classes in ranking bits
    ExplicitBitVect *dp_maskBits; // allows only certain bits to be considered
    
  };
  
}
#endif
