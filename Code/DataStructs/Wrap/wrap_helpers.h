//
//  Copyright (C) 2003-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_WRAPHELPERS_H__
#define __RD_WRAPHELPERS_H__

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <RDBoost/Wrap.h>
#include <DataStructs/base64.h>


namespace python = boost::python;

  template <typename T>
  void InitFromBase64(T& self,const std::string &inD)
  {
    self.InitFromText(inD.c_str(),inD.length(),true);
  };

  template <typename T>
  std::string ToBase64(T& self)
  {
    std::string tmp;
    tmp = self.ToString();
    return Base64Encode(tmp.c_str(),tmp.length());
  };

  template <typename T>
  void SetBitsFromList(T *bv, python::object onBitList) {
    PySequenceHolder<int> bitL(onBitList);
    for (unsigned int i = 0; i < bitL.size(); i++) {
      bv->SetBit(bitL[i]);
    }
  }

  // used to support __getitem__
  template <typename T>
  const int get_VectItem(const T& self,int which)
  {
    if(which<0){
      if(which+self.GetNumBits()<0){
        throw IndexErrorException(which);
      } else {
        which += self.GetNumBits();
      }
    }
    return self.GetBit(static_cast<unsigned int>(which));
  }

  // used to support __setitem__
  template <typename T>
  const int set_VectItem(T& self, int which, const int val)
  {
    if(which<0){
      if(which+self.GetNumBits()<0){
        throw IndexErrorException(which);
      } else {
        which += self.GetNumBits();
      }
    }
    if(val){
      return self.SetBit(static_cast<unsigned int>(which));
    } else {
      return self.UnSetBit(static_cast<unsigned int>(which));
    }
  }

  // used to support GetOnBits()
  template <typename T>
  IntVect GetOnBits(const T& self)
  {
    IntVect res;
    self.GetOnBits(res);
    return res;
  }

#endif
