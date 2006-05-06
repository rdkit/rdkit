// $Id$
//
// Copyright (c) 2003-2006 greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#include "BitVect.h"
#include <sstream>
#include <limits>
#include <RDGeneral/StreamOps.h>
#include "base64.h"
#ifdef WIN32
#include <ios>
#endif

BitVect::~BitVect() {}; // must always implement virtual destructors

void BitVect::InitFromText(const char *data,const unsigned int dataLen,
                   bool isBase64,bool allowOldFormat){
  std::stringstream ss(std::ios_base::binary|std::ios_base::in|std::ios_base::out);
  if(isBase64){
    unsigned int actualLen;
    char *decoded;
    decoded = Base64Decode((const char *)data,&actualLen);
    ss.write(decoded,actualLen);
    free(decoded);
  } else {
    ss.write(data,dataLen);
  }

  int format=0;
  unsigned int nOn=0;
  int size;
  int version=0;
  
  // earlier versions of the code did not have the version number encoded, so
  //  we'll use that to distinguish version 0
  ss.read((char *)&size,sizeof(size));
  if(size<0){
    version = -1*size;
    if (version == 16) {
      format=1;
    }
    else if (version == 32) {
      format=2;
    }
    else {
      throw ValueErrorException("bad version in BitVect pickle");
    }

    ss.read((char *)&size,sizeof(size));
  } else if( !allowOldFormat ) {
    throw ValueErrorException("invalid BitVect pickle");
  }
  ss.read((char *)&nOn,sizeof(nOn));
  _InitForSize(static_cast<int>(size));

  // if the either have older version or or version 16 with ints for on bits
  if( (format==0) || 
      ( (format == 1) && (size >= std::numeric_limits<unsigned short>::max()) ) ) {
    unsigned int tmp;
    for(unsigned int i=0; i<nOn; i++){
      ss.read((char *)&tmp,sizeof(tmp));
      SetBit(tmp);
    }
  } else if (format == 1) { // version 16 and on bits sotred as short ints
    unsigned short tmp;
    for(unsigned int i=0; i<nOn; i++){
      ss.read((char *)&tmp,sizeof(tmp));
      SetBit(tmp);
    }
  } else if (format == 2) { // run length encoded format
    unsigned int curr=0;
    for (unsigned int i=0; i<nOn; i++) {
      curr += RDKit::readPackedIntFromStream(ss);
      SetBit(curr);
      curr++;
    }
  }

}
