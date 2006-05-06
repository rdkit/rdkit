// $Id: PeriodicTable.cpp 5000 2006-02-21 01:42:36Z glandrum $
//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "PeriodicTable.h"
#include <string>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

namespace RDKit {

class PeriodicTable * PeriodicTable::ds_instance = 0;


PeriodicTable::PeriodicTable() {
  // it is assumed that the atomic atomData string constains atoms
  // in sequence and no atoms are missing in between
  byanum.clear();
  byname.clear();

  boost::char_separator<char> eolSep("\n");
  tokenizer tokens(periodicTableAtomData,eolSep);
  for(tokenizer::iterator token=tokens.begin();
      token!=tokens.end();++token){
    if(*token!=" "){
      atomicData adata(*token);
      byanum.push_back(adata);
      std::string enam = adata.Symbol();
      byname[enam] = adata.AtomicNum();
    }
  }
}


PeriodicTable *PeriodicTable::getTable()
{
  if ( ds_instance == 0 ) {
    ds_instance = new PeriodicTable();
  }
  return ds_instance;
}

} // end of namespace
