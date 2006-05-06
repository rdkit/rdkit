// $Id: Invariant.cpp 4941 2006-02-17 01:17:05Z glandrum $
//
// Copyright (C) 2001-2006 Randal M. Henne and Rational Discovery LLC
// 
//  @@ All Rights Reserved @@
//

#include "Invariant.h"

#include <string>
#include <stdio.h>
#include <iostream>

namespace Invar {

  std::ostream & 
  operator<<( std::ostream & s, Invariant & inv )
  {
    return s << inv.toString().c_str();
  }

  std::string& Invariant::
  toString()
  {
    char 
      line[10];

    sprintf( line, "%d", this->getLine() );
    
    stringRep_d += "\n" + this->getMessage() + "\nViolation occurred on line " 
      + line + " in file " + this->getFile() + "\nFailed Expression: "
      + this->getExpression() + "\n";

    return stringRep_d;
    
  }

};
