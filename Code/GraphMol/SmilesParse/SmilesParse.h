//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SMILESPARSE_H_
#define _RD_SMILESPARSE_H_

#include <string>
#include <exception>

namespace RDKit{
  class RWMol;

  RWMol *SmilesToMol(std::string smi,int debugParse=0,bool sanitize=1);
  RWMol *SmartsToMol(std::string sma,int debugParse=0,bool mergeHs=false);

  class SmilesParseException : public std::exception {
  public:
    SmilesParseException(const char *msg) : _msg(msg) {};
    SmilesParseException(const std::string msg) : _msg(msg) {};
    const char *message () const { return _msg.c_str(); };
    ~SmilesParseException () throw () {};
  private:
    std::string _msg;
  };

}

#endif
