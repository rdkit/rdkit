//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _DATASTRUCTS_EXCEPTION_H_20050126
#define _DATASTRUCTS_EXCEPTION_H_20050126

class DatastructsException : public std::exception {
 public:
  //! construct with an error message
  DatastructsException(const char *msg) : _msg(msg) {};
  //! construct with an error message
  DatastructsException(const std::string msg) : _msg(msg) {};
  //! get the error message
  const char *message () const { return _msg.c_str(); };
  ~DatastructsException () throw () {};
  private:
  std::string _msg;
};

#endif
