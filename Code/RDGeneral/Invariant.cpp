// $Id$
//
// Copyright (C) 2001-2008 Greg Landrum, Randal M. Henne, and Rational Discovery LLC
// 
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Invariant.h"

#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>

namespace Invar {

  std::ostream & 
  operator<<( std::ostream & s, Invariant & inv )
  {
    return s << inv.toString().c_str();
  }

  std::string& Invariant::
  toString()
  {
    std::string line=boost::lexical_cast<std::string>(this->getLine());
    
    stringRep_d += "\n" + this->getMessage() + "\nViolation occurred on line " 
      + line + " in file " + this->getFile() + "\nFailed Expression: "
      + this->getExpression() + "\n";

    return stringRep_d;
    
  }

};
