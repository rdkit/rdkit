// $Id$
//
// Copyright (c) 2005-2008 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FreeChemicalFeature.h"
#include <RDGeneral/StreamOps.h>
#include <sstream>
#include <boost/cstdint.hpp>


namespace ChemicalFeatures {
  using namespace RDKit;
  using boost::int32_t;
  using boost::uint32_t;
  const int ci_FEAT_VERSION=0x0020; //!< version number to use in pickles

  std::string FreeChemicalFeature::toString() const {
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    uint32_t tInt = ci_FEAT_VERSION;
    streamWrite(ss,tInt);
    
    // write the id
    streamWrite(ss,d_id);    
    
    tInt = d_family.size()+1;
    streamWrite(ss,tInt);
    ss.write(d_family.c_str(),tInt*sizeof(char));
    tInt = d_type.size()+1;
    streamWrite(ss,tInt);
    ss.write(d_type.c_str(),tInt*sizeof(char));
    streamWrite(ss,d_position.x);
    streamWrite(ss,d_position.y);
    streamWrite(ss,d_position.z);
    std::string res(ss.str());
    return res;
  };
  
  void FreeChemicalFeature::initFromString(const std::string &pickle){
    std::stringstream ss(pickle,
			 std::ios_base::binary|std::ios_base::in|std::ios_base::out);
    int version=0;
    uint32_t tInt;
    streamRead(ss,tInt);
    switch(tInt){
    case 0x0010:
      version=1;
      break;
    case 0x0020: // ok, I know this is not 2, in hex its 32 but I just want to
                 // keep consistent with the above
      version=2;
      break;
    default:
      throw("Unknown version type for FreeChemicalFeature");
    }

    // read the id, but only if a later/newish version
    if (version == 2) {
        streamRead(ss,d_id);
    }
    
    char *tmpChr;

    streamRead(ss,tInt);
    tmpChr = new char[tInt];
    ss.read(tmpChr,tInt*sizeof(char));
    d_family = tmpChr;
    delete [] tmpChr;

    streamRead(ss,tInt);
    tmpChr = new char[tInt];
    ss.read(tmpChr,tInt*sizeof(char));
    d_type = tmpChr;
    delete [] tmpChr;

    streamRead(ss,d_position.x);
    streamRead(ss,d_position.y);
    streamRead(ss,d_position.z);

  };

  
}
