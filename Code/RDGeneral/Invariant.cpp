// $Id$
//
// Copyright (C) 2001-2013 Greg Landrum, Randal M. Henne, and Rational Discovery LLC
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

#ifdef SHOW_BACKTRACES_WITH_INVARIANT_ERRORS // note: works only with gcc-derived compilers
#include <execinfo.h>
#endif

namespace Invar {

  std::ostream & 
  operator<<( std::ostream & s, const Invariant & inv )
  {
    return s << inv.toString().c_str();
  }

  std::string Invariant::toString() const
  {
    std::string line=boost::lexical_cast<std::string>(this->getLine());
    
    std::string stringRep = this->prefix_d + "\n" + this->getMessage() + "\nViolation occurred on line " 
      + line + " in file " + this->getFile() + "\nFailed Expression: "
      + this->getExpression() + "\n";

#ifdef SHOW_BACKTRACES_WITH_INVARIANT_ERRORS
    void *arr[10];
    size_t sz;
    sz = backtrace(arr,10);
    std::cerr<<" STACK TRACE\n--------------\n"<<std::endl;
    backtrace_symbols_fd(arr,sz,2);
    std::cerr<<"\n--------------\n"<<std::endl;
#endif
    
    return stringRep;
    
  }

};
