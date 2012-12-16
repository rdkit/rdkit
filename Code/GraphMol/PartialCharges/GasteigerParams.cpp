// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
#include "GasteigerParams.h"
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <sstream>
#include <locale>

namespace RDKit {
  
  /*! \brief Gasteiger partial charge parameters
   */
  std::string defaultParamData = 
"H       *      7.17    6.24    -0.56 \n \
C       sp3     7.98    9.18    1.88 \n \
C       sp2     8.79    9.32    1.51 \n \
C       sp      10.39   9.45    0.73 \n \
N       sp3     11.54   10.82   1.36 \n \
N       sp2     12.87   11.15   0.85 \n \
N       sp      15.68   11.7    -0.27 \n \
O       sp3     14.18   12.92   1.39 \n \
O       sp2     17.07   13.79   0.47 \n \
F       sp3     14.66   13.85   2.31 \n \
Cl      sp3     11.00   9.69    1.35 \n \
Br      sp3     10.08   8.47    1.16 \n \
I       sp3     9.9     7.96    0.96 \n \
S       sp3     10.14   9.13    1.38 \n \
S       so      10.14   9.13    1.38 \n \
S       so2     12.00   10.81   1.20 \n \
S       sp2     10.88   9.49    1.33 \n \
P       sp3     8.90    8.24    0.96 \n \
X       *       0.00    0.00    0.00 \n \
";

  typedef boost::flyweight<boost::flyweights::key_value<std::string,GasteigerParams>,
                           boost::flyweights::no_tracking > gparam_flyweight;

  GasteigerParams::GasteigerParams(std::string paramData) {
    boost::char_separator<char> eolSep("\n");
    boost::char_separator<char> spaceSep(" \t");
    if(paramData=="") paramData=defaultParamData;
    tokenizer lines(paramData,eolSep);
    d_paramMap.clear();
    std::istringstream istr;
    istr.imbue(std::locale("C"));
    for(tokenizer::iterator lineIter=lines.begin();
	lineIter!=lines.end();++lineIter){
      std::string dataLine = *lineIter;
      tokenizer tokens(dataLine,spaceSep);
      if(tokens.begin()!=tokens.end()){
	tokenizer::iterator tokIter = tokens.begin();

	// read the element and the mode
	std::string elem = *tokIter;
	++tokIter;
	std::string mode = *tokIter;
	++tokIter;

	// read in the parameters
	DOUBLE_VECT params(3);

        istr.clear();
        istr.str(*tokIter);
        istr>>params[0];
        ++tokIter;
        istr.clear();
        istr.str(*tokIter);
        istr>>params[1];
        ++tokIter;
        istr.clear();
        istr.str(*tokIter);
        istr>>params[2];
        ++tokIter;

	std::pair<std::string, std::string> key(elem, mode);
	d_paramMap[key] = params;
      }
    }
  }

  const GasteigerParams *GasteigerParams::getParams(const std::string &paramData) {
    const GasteigerParams *res = &(gparam_flyweight(paramData).get());
    return res;
  }

} 
