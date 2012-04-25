// $Id$
//
//  Copyright (C) 2001-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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

  unsigned int lidx=0;
  while(isotopesAtomData[lidx]!="" && isotopesAtomData[lidx]!="EOS"){
   tokenizer lines(isotopesAtomData[lidx++],eolSep);
   boost::char_separator<char> spaceSep(" \t");
   for(tokenizer::iterator line=lines.begin();
       line!=lines.end();++line){
     if(*line!=" "){
       tokenizer tokens(*line,spaceSep);
       tokenizer::iterator token=tokens.begin();
       int anum = atoi(token->c_str());
       atomicData &adata=byanum[anum];
       ++token;if(token==tokens.end()) continue;
       ++token;if(token==tokens.end()) continue;
       unsigned int isotope=static_cast<unsigned int>(atoi(token->c_str()));
       ++token;if(token==tokens.end()) continue;
       double mass=atof(token->c_str());
       ++token;if(token==tokens.end()) continue;
       double abundance=atof(token->c_str());
       adata.d_isotopeInfoMap[isotope]=std::make_pair(mass,abundance);
     }
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
