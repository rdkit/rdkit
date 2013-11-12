//
// Copyright 2003-2006 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_BADFILEEXCEPTION_H
#define _RD_BADFILEEXCEPTION_H

#include <string>
#include <vector>
#include <stdexcept>

namespace RDKit {
  
  //! used by various file parsing classes to indicate a bad file
  class BadFileException : public std::runtime_error {
  public :
    //! construct with an error message
    explicit BadFileException(const char *msg) : std::runtime_error("BadFileException"), _msg(msg) {};
    //! construct with an error message
    explicit BadFileException(const std::string msg) : std::runtime_error("BadFileException"), _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~BadFileException () throw () {};
    
    private :
      std::string _msg;
  };
}

#endif
