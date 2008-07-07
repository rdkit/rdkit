//
// Copyright 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef _RD_FILEPARSEEXCEPTION_H
#define _RD_FILEPARSEEXCEPTION_H

#include <string>
#include <exception>

namespace RDKit {
  //! used by various file parsing classes to indicate a parse error
  class FileParseException : public std::exception {
  public:
    //! construct with an error message
    explicit FileParseException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    explicit FileParseException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~FileParseException () throw () {};
  private:
    std::string _msg;
  };
}  

#endif

