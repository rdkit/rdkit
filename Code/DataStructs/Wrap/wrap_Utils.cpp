// $Id$
//
//  Copyright (C) 2003-2012 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitVectUtils.h>
#include <DataStructs/BitOps.h>
namespace python = boost::python;

ExplicitBitVect *createFromBitString(const std::string &bits){
  ExplicitBitVect *res=new ExplicitBitVect(bits.length());
  FromBitString(*res,bits);
  return res;
}

ExplicitBitVect *createFromFPSText(const std::string &fps){
  if(fps.length()%2){
    throw ValueErrorException("input string must have an even number of characters");
  }
  ExplicitBitVect *res=new ExplicitBitVect(fps.length()*4);
  UpdateBitVectFromFPSText(*res,fps);
  return res;
}

ExplicitBitVect *createFromBinaryText(const std::string &fps){
  ExplicitBitVect *res=new ExplicitBitVect(fps.length()*8);
  UpdateBitVectFromBinaryText(*res,fps);
  return res;
}


struct Utils_wrapper {
  static void wrap(){
    python::def("ConvertToExplicit", convertToExplicit, 
                python::return_value_policy<python::manage_new_object>(),
                "Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector");  
    python::def("CreateFromBitString",createFromBitString,
                python::return_value_policy<python::manage_new_object>(),
                "Creates an ExplicitBitVect from a bit string (string of 0s and 1s).");  
    python::def("CreateFromFPSText",createFromFPSText,
                python::return_value_policy<python::manage_new_object>(),
                "Creates an ExplicitBitVect from an FPS string.");  
    python::def("CreateFromBinaryText",createFromBinaryText,
                python::return_value_policy<python::manage_new_object>(),
                "Creates an ExplicitBitVect from a binary string (byte array).");  

    python::def("InitFromDaylightString",
                (void (*)(SparseBitVect &,std::string))FromDaylightString);
    python::def("InitFromDaylightString",
                (void (*)(ExplicitBitVect &,std::string))FromDaylightString,
                "Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.\n\
\n\
   **Arguments**\n\
     - bv: either a _SparseBitVect_ or an _ExplicitBitVect_\n\
     - txt: a string with the Daylight encoding (this is the text that\n\
            the Daylight tools put in the FP field of a TDT)\n\
\n");
  }
};

void wrap_Utils() {
  Utils_wrapper::wrap();
}
  
