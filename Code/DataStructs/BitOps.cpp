// $Id$
//
//  Copyright (C) 2003-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include <RDBoost/Exceptions.h>
#include "BitVects.h"
#include "BitOps.h"
#include <math.h>
#include <string>
#include <iostream>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>

int getBitId(const char *&text,int format,int size,int curr){
  PRECONDITION(text,"no text");
  int res=-1;
  if( (format==0) || 
      ( (format == 1) && (size >= std::numeric_limits<unsigned short>::max()) ) ) {
    int tmp;
    tmp = *(int *)text;
    text += sizeof(tmp);
    res=tmp;
  } else if (format == 1) { // version 16 and on bits sotred as short ints
    unsigned short tmp;
    tmp = *(unsigned short *)text;
    text += sizeof(tmp);
    res=tmp;
  } else if (format == 2) { // run length encoded format
    res = curr + RDKit::pullPackedIntFromString(text);
  }
  return res;
}



bool AllProbeBitsMatch(const std::string &probe,const std::string &ref){
  return AllProbeBitsMatch(probe.c_str(),ref.c_str());
}

bool AllProbeBitsMatch(const char *probe,const char *ref){
  PRECONDITION(probe,"no probe text");
  PRECONDITION(ref,"no probe text");
  int probeFormat=0;
  int refFormat=0;
  int version=0;

  int probeSize = *(int *)probe;
  probe+=sizeof(probeSize);
  if(probeSize<0){
    version = -1*probeSize;
    if (version == 16) {
      probeFormat=1;
    }
    else if (version == 32) {
      probeFormat=2;
    }
    else {
      throw("Unknown version type for the encode bit vect");
    }
    probeSize = *(int *)probe;
    probe+=sizeof(probeSize);
  }

  int refSize = *(int *)ref;
  ref+=sizeof(refSize);
  if(refSize<0){
    version = -1*refSize;
    if (version == 16) {
      refFormat=1;
    }
    else if (version == 32) {
      refFormat=2;
    }
    else {
      throw("Unknown version type for the encode bit vect");
    }
    refSize = *(int *)ref;
    ref+=sizeof(refSize);
  }


  int nProbeOn = *(int *)probe;
  probe+=sizeof(nProbeOn);
  int nRefOn = *(int *)ref;
  ref+=sizeof(nRefOn);

  int currProbeBit=0;
  currProbeBit=getBitId(probe,probeFormat,probeSize,currProbeBit);
  nProbeOn--;

  int currRefBit=0;
  currRefBit=getBitId(ref,refFormat,refSize,currRefBit);
  nRefOn--;

  while(nProbeOn){
    while(currRefBit<currProbeBit && nRefOn>0){
      if(refFormat==2) currRefBit++;
      currRefBit=getBitId(ref,refFormat,refSize,currRefBit);
      nRefOn--;
    }
    if(currRefBit!=currProbeBit) return false;
    if(probeFormat==2) currProbeBit++;
    currProbeBit=getBitId(probe,probeFormat,probeSize,currProbeBit);
    nProbeOn--;
  }
  return true;
}

template <typename T1>
bool AllProbeBitsMatch(const T1 &probe,const std::string &pkl){
  const char *text=pkl.c_str();
  int format=0;
  int nOn=0,size,version=0;
  size = *(int *)text;
  text+=sizeof(size);
  if(size<0){
    version = -1*size;
    if (version == 16) {
      format=1;
    }
    else if (version == 32) {
      format=2;
    }
    else {
      throw("Unknown version type for the encode bit vect");
    }
    size = *(int *)text;
    text+=sizeof(size);
  }
  nOn = *(int *)text;
  text+=sizeof(nOn);

  int currBit=0;
  currBit=getBitId(text,format,size,currBit);
  nOn--;
  std::vector<int> obl;
  probe.GetOnBits(obl);
  

  //for(int i=0;i<probe.GetNumBits();i++){
  //  if(probe.GetBit(i)){
  for(std::vector<int>::const_iterator i=obl.begin();i!=obl.end();i++){
    while(currBit<*i && nOn>0){
        if(format==2) currBit++;
        currBit=getBitId(text,format,size,currBit);
        nOn--;
      }
      if(currBit!=*i) return false;
      //}
  }
  return true;
}
template bool AllProbeBitsMatch(const SparseBitVect& bv1,const std::string &pkl);
template bool AllProbeBitsMatch(const ExplicitBitVect& bv1,const std::string &pkl);

// """ -------------------------------------------------------
//
//  NumOnBitsInCommon(T1,T2)
//  Returns the number of on bits which are set in both T1 and T2.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
int
NumOnBitsInCommon(const T1& bv1,
                const T2& bv2)
{
  return OnBitsInCommon(bv1,bv2).size();
}

int
NumOnBitsInCommon(const ExplicitBitVect& bv1,
                  const ExplicitBitVect& bv2)
{
  //std::cout << "nobic" << std::endl;
  int res = 0;
  unsigned int _sz = bv1.GetNumBits()<bv2.GetNumBits()?bv1.GetNumBits():bv2.GetNumBits();
  for(unsigned int i=0;i<_sz;i++) {
    if((*bv1.dp_bits)[i] && (*bv2.dp_bits)[i]) res+=1;
  }
  return res;
}


// In all these similarity metrics the notation is selected to be
//   consistent with J.W. Raymond and P. Willett, JCAMD _16_ 59-71 (2002)
// """ -------------------------------------------------------
//
//  TanimotoSimilarity(T1,T2)
//   returns the Tanamoto similarity between T1 and T2, a double.
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator&, GetNumBits()
//  and GetOnBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
const double
TanimotoSimilarity(const T1& bv1,
                   const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();
  if((y+z-x)==0.0) return 1.0;
  else return x / (y+z-x);
}


template <typename T1, typename T2>
const double
TverskySimilarity(const T1& bv1,
                  const T2& bv2,
                  double a,
                  double b)
{
  RANGE_CHECK(0,a,1);
  RANGE_CHECK(0,b,1);
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();
  double denom = a*y + b*z + (1-a-b)*x;
  if(denom==0.0) return 1.0;
  else return x / denom;
}

template <typename T1, typename T2>
const double
CosineSimilarity(const T1& bv1,
                   const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();

  if(y*z>0.0){
    return x / sqrt(y*z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
const double
KulczynskiSimilarity(const T1& bv1,
                   const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();

  if(y*z>0.0){
    return x*(y+z)/(2*y*z);
  } else {
    return 0.0;
  }
}


template <typename T1, typename T2>
const double
DiceSimilarity(const T1& bv1,
              const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();

  if(y+z>0.0){
    return 2*x/(y+z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
const double
SokalSimilarity(const T1& bv1,
              const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();

  return x/(2*y+2*z-3*x);
}

template <typename T1, typename T2>
const double
McConnaugheySimilarity(const T1& bv1,
                       const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();

  if(y*z>0.0){
    return (x*(y+z)-(y*z))/(y*z);
  } else {
    return 0.0;
  }
}

template <typename T>
inline T tmin(T v1,T v2) {
  if(v1<v2) return v1;
  return v2;
}

template <typename T>
inline T tmax(T v1,T v2) {
  if(v1>v2) return v1;
  return v2;
}

template <typename T1, typename T2>
const double
AsymmetricSimilarity(const T1& bv1,
                     const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();
  
  if(tmin(y,z)>0){
    return x/tmin(y,z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
const double
BraunBlanquetSimilarity(const T1& bv1,
                     const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  double y = bv1.GetNumOnBits();
  double z = bv2.GetNumOnBits();

  if(tmax(y,z)>0){
    return x/tmax(y,z);
  } else {
    return 0.0;
  }
}

template <typename T1, typename T2>
const double
RusselSimilarity(const T1& bv1,
                     const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  double x = NumOnBitsInCommon(bv1,bv2);
  return x/bv1.GetNumBits();
}




// """ -------------------------------------------------------
//
//  OnBitSimilarity(T1,T2)
//  Returns the percentage of possible on bits in common
//  between T1 and T2 (a double)
//
//  C++ Notes: T1 and T2 must support operator|, operator&
//  and GetOnBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
const double
OnBitSimilarity(const T1& bv1,
                const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");

  double num = NumOnBitsInCommon(bv1,bv2);
  double denom=(bv1|bv2).GetNumOnBits();

  if(denom>0){
    return num/denom;
  } else {
    return 0;
  }
}

// """ -------------------------------------------------------
//
//  NumBitsInCommon(T1,T2)
//  Returns the number of bits in common (on and off)
//  between T1 and T2 (an int)
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator^, GetNumBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
const int
NumBitsInCommon(const T1& bv1,
                const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");

  return bv1.GetNumBits() - (bv1^bv2).GetNumOnBits();
}


// """ -------------------------------------------------------
//
//  AllBitSimilarity(T1,T2)
//  Returns the percentage of bits in common (on and off)
//  between T1 and T2 (a double)
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator^, GetNumBits()
//  and GetNumOnBits().
//
//  Python Notes: T1 and T2 are BitVects.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
const double
AllBitSimilarity(const T1& bv1,
                const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");

  return double(NumBitsInCommon(bv1,bv2))/bv1.GetNumBits();
}




// """ -------------------------------------------------------
//
//  OnBitsInCommon(T1,T2)
//  Returns the on bits which are set in both T1 and T2.
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator&, GetNumBits()
//  and GetOnBits(), the return value is an IntVect.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a tuple of ints.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
IntVect
OnBitsInCommon(const T1& bv1,
                const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  IntVect res;
  (bv1&bv2).GetOnBits(res);
  return res;
}

// """ -------------------------------------------------------
//
//  OffBitsInCommon(T1,T2)
//  Returns the off bits which are set in both T1 and T2.
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator|, operator~,
// GetNumBits() and GetOnBits(), the return value is an IntVect.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a tuple of ints.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
IntVect
OffBitsInCommon(const T1& bv1,
                const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  IntVect res;
  (~(bv1|bv2)).GetOnBits(res);
  return res;
}

// """ -------------------------------------------------------
//
//  OnBitProjSimilarity(T1,T2)
//  Returns the projected similarity between the on bits of
//  T1 and T2.
//
//  The on bit projected similarity of T1 onto T2 is the
//  percentage of T1's on bits which are on in T2.
//
//  This type of measure may be useful for substructure-type
//  searches.
//
//  Two values are returned, the projection of T1 onto T2
//  and the projection of T2 onto T1
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator&, GetNumBits()
//  and GetNumOnBits(), the return value is an DoubleVect with
//  two elements.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a 2-tuple of doubles.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
DoubleVect
OnBitProjSimilarity(const T1& bv1,
                     const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  DoubleVect res(2,0.0);
  double num=NumOnBitsInCommon(bv1,bv2);
  if(num){
    res[0] = num/bv1.GetNumOnBits();
    res[1] = num/bv2.GetNumOnBits();
  }
  return res;
}

// """ -------------------------------------------------------
//
//  OffBitProjSimilarity(T1,T2)
//  Returns the projected similarity between the off bits of
//  T1 and T2.
//
//  The off bit projected similarity of T1 onto T2 is the
//  percentage of T1's off bits which are off in T2.
//
//  This type of measure may be useful for substructure-type
//  searches.
//
//  Two values are returned, the projection of T1 onto T2
//  and the projection of T2 onto T1
//
//  T1 and T2 should be the same length.
//
//  C++ Notes: T1 and T2 must support operator|, GetNumBits()
//  and GetNumOffBits(), the return value is an DoubleVect with
//  two elements.
//
//  Python Notes: T1 and T2 are BitVects, the return value
//  is a 2-tuple of doubles.
//
// """ -------------------------------------------------------
template <typename T1, typename T2>
DoubleVect
OffBitProjSimilarity(const T1& bv1,
                     const T2& bv2)
{
  if(bv1.GetNumBits()!=bv2.GetNumBits())
    throw ValueErrorException("BitVects must be same length");
  DoubleVect res(2,0.0);
  double num=(bv1|bv2).GetNumOffBits();
  if(num){
    res[0] = num/bv1.GetNumOffBits();
    res[1] = num/bv2.GetNumOffBits();
  }
  return res;
}




template <typename T1>
T1 *
FoldFingerprint(const T1 &bv1,unsigned int factor)
{
  if(factor <=0 || factor >= bv1.GetNumBits())
    throw ValueErrorException("invalid fold factor");

  int initSize = bv1.GetNumBits();
  int resSize = initSize/factor;
  T1 *res = new T1(resSize);

  IntVect onBits;
  bv1.GetOnBits(onBits);
  for(IntVectIter iv=onBits.begin();iv!=onBits.end();iv++){
    int pos = (*iv) % resSize;
    res->SetBit(pos);
  }
  return res;
}

template <typename T1>
std::string
BitVectToText(const T1& bv1){
  std::string res(bv1.GetNumBits(),'0');
  for(unsigned int i=0;i<bv1.GetNumBits();i++){
    if(bv1.GetBit(i)) res[i] = '1';
  }
  return res;
}



template const double TanimotoSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double TverskySimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2,double a, double b);
template const double CosineSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double KulczynskiSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double DiceSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double SokalSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double McConnaugheySimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double AsymmetricSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double BraunBlanquetSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double RusselSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double OnBitSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const int NumBitsInCommon(const SparseBitVect& bv1,const SparseBitVect& bv2);
template const double AllBitSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template int NumOnBitsInCommon(const SparseBitVect& bv1,const SparseBitVect& bv2);
template IntVect OnBitsInCommon(const SparseBitVect& bv1,const SparseBitVect& bv2);
template IntVect OffBitsInCommon(const SparseBitVect& bv1,const SparseBitVect& bv2);
template DoubleVect OnBitProjSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);
template DoubleVect OffBitProjSimilarity(const SparseBitVect& bv1,const SparseBitVect& bv2);

template const double TanimotoSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double TverskySimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2,double a, double b);
template const double CosineSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double KulczynskiSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double DiceSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double SokalSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double McConnaugheySimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double AsymmetricSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double BraunBlanquetSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double RusselSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double OnBitSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const int NumBitsInCommon(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template const double AllBitSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template IntVect OnBitsInCommon(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template IntVect OffBitsInCommon(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template DoubleVect OnBitProjSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);
template DoubleVect OffBitProjSimilarity(const ExplicitBitVect& bv1,const ExplicitBitVect& bv2);

template SparseBitVect *FoldFingerprint(const SparseBitVect &,unsigned int);
template ExplicitBitVect *FoldFingerprint(const ExplicitBitVect &,unsigned int);

template std::string BitVectToText(const SparseBitVect &);
template std::string BitVectToText(const ExplicitBitVect &);




