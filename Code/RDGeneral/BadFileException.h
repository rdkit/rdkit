//
// Copyright 2003-2006 Greg Landrum and Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef _RD_BADFILEEXCEPTION_H
#define _RD_BADFILEEXCEPTION_H

#include <string>
#include <vector>
#include <exception>

namespace RDKit {
  
  //! used by various file parsing classes to indicate a bad file
  class BadFileException : public std::exception {
  public :
    //! construct with an error message
    explicit BadFileException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    explicit BadFileException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~BadFileException () throw () {};
    
    private :
      std::string _msg;
  };
}

#endif
