// $Id: FreeChemicalFeature.cpp 4943 2006-02-17 01:21:27Z glandrum $
//
// Copyright (c) 2005-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include "FreeChemicalFeature.h"
#include <sstream>


namespace ChemicalFeatures {
  const int ci_FEAT_VERSION=0x0010; //!< version number to use in pickles

  std::string FreeChemicalFeature::toString() const {
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    unsigned int tInt = ci_FEAT_VERSION;
    ss.write((const char *)&(tInt),sizeof(tInt));
    tInt = d_family.size()+1;
    ss.write((const char *)&(tInt),sizeof(tInt));
    ss.write(d_family.c_str(),tInt*sizeof(char));
    tInt = d_type.size()+1;
    ss.write((const char *)&(tInt),sizeof(tInt));
    ss.write(d_type.c_str(),tInt*sizeof(char));
    ss.write((const char *)&d_position.x,sizeof(d_position.x));
    ss.write((const char *)&d_position.y,sizeof(d_position.y));
    ss.write((const char *)&d_position.z,sizeof(d_position.z));
    std::string res(ss.str());
    return res;
  };
  void FreeChemicalFeature::initFromString(const std::string &pickle){
    std::stringstream ss(pickle,
			 std::ios_base::binary|std::ios_base::in|std::ios_base::out);
    int version=0;
    unsigned int tInt;
    ss.read((char *)&tInt,sizeof(tInt));
    switch(tInt){
    case 0x0010:
      version=1;
      break;
    default:
      throw("Unknown version type for FreeChemicalFeature");
    }

    char *tmpChr;
    ss.read((char *)&(tInt),sizeof(tInt));
    tmpChr = new char[tInt];
    ss.read(tmpChr,tInt*sizeof(char));
    d_family = tmpChr;
    delete [] tmpChr;

    ss.read((char *)&(tInt),sizeof(tInt));
    tmpChr = new char[tInt];
    ss.read(tmpChr,tInt*sizeof(char));
    d_type = tmpChr;
    delete [] tmpChr;

    ss.read((char *)&d_position.x,sizeof(d_position.x));
    ss.read((char *)&d_position.y,sizeof(d_position.x));
    ss.read((char *)&d_position.z,sizeof(d_position.x));
  };

  
}
