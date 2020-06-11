//
// Copyright (C) 2001-2020 Greg Landrum, Randal M. Henne, and Rational Discovery
// LLC
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
#include "versions.h"

#ifdef RDK_USE_BOOST_STACKTRACE
#include <boost/stacktrace.hpp>
#include <sstream>
#endif

namespace Invar {

std::ostream &operator<<(std::ostream &s, const Invariant &inv) {
  return s << inv.toString().c_str();
}

std::string Invariant::toString() const {
  std::string line = std::to_string(this->getLine());

  std::string stringRep =
      this->prefix_d + "\n" + this->what() + "\nViolation occurred on line " +
      line + " in file " + this->getFile() +
      "\nFailed Expression: " + this->getExpression() + "\n";
#ifdef RDK_USE_BOOST_STACKTRACE
  std::stringstream sstr;
  sstr << "----------\n"
       << "Stacktrace:\n"
       << boost::stacktrace::stacktrace() << "----------\n";
  stringRep += sstr.str();
#endif
  return stringRep;
}

std::string Invariant::toUserString() const {
  std::string line = std::to_string(this->getLine());

  std::string filename = this->getFile();

  std::size_t pos = filename.find("Code");  // strip out build directory info
  if (pos != std::string::npos) {
    filename = filename.substr(pos);
  }

  std::string stringRep = this->prefix_d + "\n\t" + this->what() +
                          "\n\tViolation occurred on line " + line +
                          " in file " + filename +
                          "\n\tFailed Expression: " + this->getExpression() +
                          "\n\t" + "RDKIT: " + RDKit::rdkitVersion + "\n\t" +
                          "BOOST: " + RDKit::boostVersion + "\n";

#ifdef SHOW_BACKTRACES_WITH_INVARIANT_ERRORS
  void *arr[10];
  size_t sz;
  sz = backtrace(arr, 10);
  std::cerr << " STACK TRACE\n--------------\n" << std::endl;
  backtrace_symbols_fd(arr, sz, 2);
  std::cerr << "\n--------------\n" << std::endl;
#endif

  return stringRep;
}

};  // namespace Invar
