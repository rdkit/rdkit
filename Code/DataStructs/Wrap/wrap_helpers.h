//
//  Copyright (C) 2003-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
    self.initFromText(inD.c_str(),inD.length(),true);
  };

  template <typename T>
  std::string ToBase64(T& self)
  {
    std::string tmp;
    tmp = self.toString();
    const char *txt=Base64Encode(tmp.c_str(),tmp.length());
    std::string res(txt);
    delete[] txt;
    return res ;
  };

  template <typename T>
  void SetBitsFromList(T *bv, python::object onBitList) {
    PySequenceHolder<int> bitL(onBitList);
    for (unsigned int i = 0; i < bitL.size(); i++) {
      bv->setBit(bitL[i]);
    }
  }

  template <typename T>
  void UnSetBitsFromList(T *bv, python::object offBitList) {
    PySequenceHolder<int> bitL(offBitList);
    for (unsigned int i = 0; i < bitL.size(); i++) {
      bv->unsetBit(bitL[i]);
    }
  }

  // used to support __getitem__
  template <typename T>
  const int get_VectItem(const T& self,int which)
  {
    if(which<0){
      if(which+static_cast<int>(self.getNumBits())<0){
        throw IndexErrorException(which);
      } else {
        which += self.getNumBits();
      }
    }
    return self.getBit(static_cast<unsigned int>(which));
  }

  // used to support __setitem__
  template <typename T>
  const int set_VectItem(T& self, int which, const int val)
  {
    if(which<0){
      if(which+static_cast<int>(self.getNumBits())<0){
        throw IndexErrorException(which);
      } else {
        which += self.getNumBits();
      }
    }
    if(val){
      return self.setBit(static_cast<unsigned int>(which));
    } else {
      return self.unsetBit(static_cast<unsigned int>(which));
    }
  }

  // used to support getOnBits()
  template <typename T>
  IntVect GetOnBits(const T& self)
  {
    IntVect res;
    self.getOnBits(res);
    return res;
  }

  template <typename T>
  python::object BVToBinary(const T &bv){
    std::string res=bv.toString();
    python::object retval = python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length())));
    return retval;
  }


#endif
