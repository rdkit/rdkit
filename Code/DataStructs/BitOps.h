//
//  Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_BITOPS_H__
#define __RD_BITOPS_H__
/*! \file BitOps.h

  \brief Contains general bit-comparison and similarity operations.

  The notation used to document the similarity metrics is:
    - \c V1_n: number of bits in vector 1
    - \c V1_o: number of on bits in vector 1
    - <tt>(V1&V2)_o</tt>: number of on bits in the intersection of vectors 1 and 2
  
 */

#include "BitVects.h"
#include <string>


//! general purpose wrapper for calculating the similarity between two bvs
//! that may be of unequal size (will automatically fold as appropriate)
template <typename T>
double SimilarityWrapper(const T &bv1,const T &bv2,
                         const double (*metric)(const T &,const T &),
                         bool returnDistance=false){
  double res=0.0;
  if(bv1.GetNumBits()>bv2.GetNumBits()){
    T *bv1tmp = FoldFingerprint(bv1,bv1.GetNumBits()/bv2.GetNumBits());
    res = metric(*bv1tmp,bv2);
    delete bv1tmp;
  } else if(bv2.GetNumBits()>bv1.GetNumBits()){
    T *bv2tmp = FoldFingerprint(bv2,bv2.GetNumBits()/bv1.GetNumBits());
    res = metric(bv1,*bv2tmp);
    delete bv2tmp;
  } else {
    res = metric(bv1,bv2);
  }
  if(returnDistance) res = 1.0-res;
  return res;
}


bool AllProbeBitsMatch(const char *probe,const char *ref);
bool AllProbeBitsMatch(const std::string &probe,const std::string &ref);

  
template <typename T1>
bool AllProbeBitsMatch(const T1 &probe,const std::string &pkl);


//! returns the number of on bits in common between two bit vectors
/*!
  \return (bv1&bv2)_o
*/
template <typename T1, typename T2>
int
NumOnBitsInCommon(const T1& bv1,const T2& bv2);

int
NumOnBitsInCommon(const ExplicitBitVect & bv1,const ExplicitBitVect & bv2);

//! returns the Tanimoto similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / [bv1_o + bv2_o - (bv1&bv2)_o]</tt>
*/
template <typename T1, typename T2>
const double
TanimotoSimilarity(const T1& bv1,const T2& bv2);

//! returns the Cosine similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / sqrt(bv1_o + bv2_o)</tt>
*/
template <typename T1, typename T2>
const double
CosineSimilarity(const T1& bv1,
                 const T2& bv2);

//! returns the Kulczynski similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o * [bv1_o + bv2_o] / [2 * bv1_o * bv2_o]</tt>
*/
template <typename T1, typename T2>
const double
KulczynskiSimilarity(const T1& bv1,
                     const T2& bv2);

//! returns the Dice similarity between two bit vects
/*!
  \return <tt>2*(bv1&bv2)_o / [bv1_o + bv2_o]</tt>
*/
template <typename T1, typename T2>
const double
DiceSimilarity(const T1& bv1,
               const T2& bv2);

//! returns the Sokal similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / [2*bv1_o + 2*bv2_o - 3*(bv1&bv2)_o]</tt>
*/
template <typename T1, typename T2>
const double
SokalSimilarity(const T1& bv1,
                const T2& bv2);

//! returns the McConnaughey similarity between two bit vects
/*!
  \return <tt>[(bv1&bv2)_o * (bv1_o + bv2_o) - (bv1_o * bv2_o)] / (bv1_o * bv2_o)</tt>
*/
template <typename T1, typename T2>
const double
McConnaugheySimilarity(const T1& bv1,
                       const T2& bv2);

//! returns the Asymmetric similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / min(bv1_o,bv2_o)</tt>
*/
template <typename T1, typename T2>
const double
AsymmetricSimilarity(const T1& bv1,
                     const T2& bv2);

//! returns the Braun-Blanquet similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / max(bv1_o,bv2_o)</tt>
*/
template <typename T1, typename T2>
const double
BraunBlanquetSimilarity(const T1& bv1,
                        const T2& bv2);

//! returns the Russel similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / bv1_o</tt>

  <b>Note:</b> that this operation is non-commutative:
    RusselSimilarity(bv1,bv2) != RusselSimilarity(bv2,bv1)

*/
template <typename T1, typename T2>
const double
RusselSimilarity(const T1& bv1,
                 const T2& bv2);


//! returns the on bit similarity between two bit vects
/*!
  \return <tt>(bv1&bv2)_o / (bv1|bv2)_o </tt>
*/
template <typename T1, typename T2>
const double
OnBitSimilarity(const T1& bv1,const T2& bv2);

//! returns the number of common bits (on and off) between two bit vects
/*!
  \return <tt>bv1_n - (bv1^bv2)_o</tt>
*/
template <typename T1, typename T2>
const int
NumBitsInCommon(const T1& bv1,const T2& bv2);

//! returns the commong-bit similarity (on and off) between two bit vects
/*!
  \return <tt>[bv1_n - (bv1^bv2)_o] / bv1_n</tt>
*/
template <typename T1, typename T2>
const double
AllBitSimilarity(const T1& bv1,const T2& bv2);

//! returns an IntVect with indices of all on bits in common between two bit vects
template <typename T1, typename T2>
IntVect
OnBitsInCommon(const T1& bv1,const T2& bv2);

//! returns an IntVect with indices of all off bits in common between two bit vects
template <typename T1, typename T2>
IntVect
OffBitsInCommon(const T1& bv1,const T2& bv2);

//! returns the on-bit projected similarities between two bit vects
/*!
  \return two values, as a DoubleVect:
      - <tt>(bv1&bv2)_o / bv1_o</tt> 
      - <tt>(bv1&bv2)_o / bv2_o</tt> 
*/
template <typename T1, typename T2>
DoubleVect
OnBitProjSimilarity(const T1& bv1,const T2& bv2);

//! returns the on-bit projected similarities between two bit vects
/*!
  \return two values, as a DoubleVect:
     - <tt>[bv1_n - (bv1|bv2)_o] / [bv1_n - bv1_o]</tt> 
     - <tt>[bv2_n - (bv1|bv2)_o] / [bv2_n - bv2_o]</tt> 

   <b>Note:</b> <tt>bv1_n = bv2_n</tt>
      
*/
template <typename T1, typename T2>
DoubleVect
OffBitProjSimilarity(const T1& bv1,const T2& bv2);


//! folds a bit vector \c factor times and returns the result
/*!
  \param bv1    the vector to be folded
  \param factor (optional) the number of times to fold it
  
  \return a pointer to the folded fingerprint, which is
     <tt>bv1_n/factor</tt> long.
     
   <b>Note:</b> The caller is responsible for <tt>delete</tt>ing the result.
 */
template <typename T1>
T1 *
FoldFingerprint(const T1& bv1,unsigned int factor=2);

//! returns a text representation of a bit vector (a string of 0s and 1s)
/*!
  \param bv1    the vector to be folded
  
  \return an std::string

 */
template <typename T1>
std::string
BitVectToText(const T1& bv1);



#endif
